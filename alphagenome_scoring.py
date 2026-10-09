"""
alphagenome_scoring.py
---------------------
Score every variant in variants_all_te.tsv with AlphaGenome's recommended
variant scorers (default: RNA_SEQ only; see --scorers).

- Variants are deduplicated on (Chromosome, Position, Ref, Alt) before
  scoring (the same variant can appear under several elements, e.g. SORT1 and
  SORT1-flip); every output row carries the variant key so it can be joined
  back to the input table.
- Deletions (Alt == "-") are converted to VCF style by left-padding with the
  preceding reference base, fetched from the UCSC REST API.
- Results are written in chunks to <outdir>/scores_chunk_NNNNN.parquet
  (default outdir: data/alphagenome_scores_<scorers>).
  Finished chunks are skipped on rerun, so the job can be restarted after a
  crash or walltime limit without re-scoring anything.

Usage
-----
    export ALPHAGENOME_API_KEY=...
    python alphagenome_scoring.py                       # all variants
    python alphagenome_scoring.py --limit 5             # quick test
    python alphagenome_scoring.py --skip-deletions      # SNVs only
    python alphagenome_scoring.py --scorers DNASE CAGE  # other track types
"""

import argparse
import glob
import os
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"

parser = argparse.ArgumentParser(description="Score variants with AlphaGenome.")
parser.add_argument("--input", default=f"{DATA_DIR}/variants_all_te.tsv")
parser.add_argument("--outdir", default=None,
                    help="Defaults to data/alphagenome_scores_<scorers>.")
parser.add_argument("--sequence-length", default="1MB",
                    choices=["2KB", "16KB", "100KB", "500KB", "1MB"])
parser.add_argument("--chunk-size", type=int, default=100,
                    help="Variants per output file.")
parser.add_argument("--max-workers", type=int, default=5,
                    help="Concurrent AlphaGenome requests.")
parser.add_argument("--limit", type=int, default=None,
                    help="Only score the first N unique variants (for testing).")
parser.add_argument("--skip-deletions", action="store_true")
parser.add_argument("--scorers", nargs="+", default=["RNA_SEQ"],
                    choices=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS),
                    help="Which recommended scorers to run (default: RNA_SEQ).")
args = parser.parse_args()
if args.outdir is None:
    args.outdir = f"{DATA_DIR}/alphagenome_scores_{'_'.join(args.scorers)}"

KEY_COLS = ["Chromosome", "Position", "Ref", "Alt"]


# ── Load and deduplicate variants ───────────────────────────────────────────
df = pd.read_csv(args.input, sep="\t", dtype={"Chromosome": str})
variants_df = df[KEY_COLS].drop_duplicates().reset_index(drop=True)
variants_df["chrom"] = "chr" + variants_df["Chromosome"].str.removeprefix("chr")

is_del = variants_df["Alt"] == "-"
if args.skip_deletions:
    variants_df = variants_df[~is_del].reset_index(drop=True)
    is_del = variants_df["Alt"] == "-"
if args.limit:
    variants_df = variants_df.head(args.limit)
    is_del = is_del.head(args.limit)

print(f"{len(df):,} rows -> {len(variants_df):,} unique variants "
      f"({is_del.sum():,} deletions)")


# ── Resolve deletions to VCF-style alleles ──────────────────────────────────
session = requests.Session()


def fetch_ref(chrom, start0, end0):
    """hg38 sequence for the 0-based half-open interval [start0, end0)."""
    params = {"genome": "hg38", "chrom": chrom, "start": start0, "end": end0}
    for attempt in range(3):
        try:
            resp = session.get("https://api.genome.ucsc.edu/getData/sequence",
                               params=params, timeout=30)
            resp.raise_for_status()
            return resp.json()["dna"].upper()
        except requests.RequestException:
            if attempt == 2:
                raise
            time.sleep(2)


del_positions = list({(r.chrom, r.Position) for r in variants_df[is_del].itertuples()})
if del_positions:
    print(f"Fetching anchor bases for {len(del_positions):,} deletion positions...")
    with ThreadPoolExecutor(max_workers=10) as pool:
        # Two bases: the anchor (Position - 1) and the deleted base itself.
        seqs = pool.map(lambda cp: fetch_ref(cp[0], cp[1] - 2, cp[1]), del_positions)
        ref_by_pos = dict(zip(del_positions, seqs))


def to_variant(row):
    if row.Alt == "-":
        seq = ref_by_pos[(row.chrom, row.Position)]
        if seq[1] != row.Ref.upper():
            print(f"WARNING: Ref mismatch at {row.chrom}:{row.Position} "
                  f"(table {row.Ref}, hg38 {seq[1]})")
        return genome.Variant(chromosome=row.chrom, position=row.Position - 1,
                              reference_bases=seq, alternate_bases=seq[0])
    return genome.Variant(chromosome=row.chrom, position=row.Position,
                          reference_bases=row.Ref, alternate_bases=row.Alt)


variants_df["variant"] = [to_variant(r) for r in variants_df.itertuples()]


# ── Score in chunks ─────────────────────────────────────────────────────────
dna_model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
    f"SEQUENCE_LENGTH_{args.sequence_length}"
]
scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[s] for s in args.scorers]
print(f"Scorers: {', '.join(args.scorers)}")
os.makedirs(args.outdir, exist_ok=True)


def tidy(chunk, scores):
    """Tidy one chunk's scores and tag each row with its input variant key."""
    out = []
    for row, variant_scores in zip(chunk.itertuples(), scores):
        t = variant_scorers.tidy_scores(variant_scores)
        # Variant/Interval objects can't be stored in parquet.
        t["variant_id"] = t["variant_id"].astype(str)
        t["scored_interval"] = t["scored_interval"].astype(str)
        for col in KEY_COLS:
            t[col] = getattr(row, col)
        out.append(t)
    return pd.concat(out, ignore_index=True)


failed = []
n_chunks = (len(variants_df) + args.chunk_size - 1) // args.chunk_size
for i in range(n_chunks):
    out_path = os.path.join(args.outdir, f"scores_chunk_{i:05d}.parquet")
    if os.path.exists(out_path):
        continue
    chunk = variants_df.iloc[i * args.chunk_size:(i + 1) * args.chunk_size]
    print(f"[{time.strftime('%H:%M:%S')}] chunk {i + 1}/{n_chunks}")

    variants = list(chunk["variant"])
    intervals = [v.reference_interval.resize(sequence_length) for v in variants]
    try:
        scores = dna_model.score_variants(
            intervals=intervals, variants=variants, variant_scorers=scorers,
            max_workers=args.max_workers, progress_bar=False,
        )
    except Exception as e:
        # One bad variant fails the whole batch; retry one by one so the rest
        # of the chunk still gets scored.
        print(f"  batch failed ({e}); retrying variants individually")
        keep, scores = [], []
        for j, (iv, v) in enumerate(zip(intervals, variants)):
            try:
                scores.append(dna_model.score_variant(
                    interval=iv, variant=v, variant_scorers=scorers))
                keep.append(j)
            except Exception as e2:
                print(f"  FAILED {v}: {e2}")
                failed.append({**chunk.iloc[j][KEY_COLS].to_dict(),
                               "variant": str(v), "error": str(e2)})
        chunk = chunk.iloc[keep]
        if not scores:
            continue

    result = tidy(chunk, scores)
    # Write to a temp file first so an interrupted write never looks finished.
    result.to_parquet(out_path + ".tmp", index=False)
    os.replace(out_path + ".tmp", out_path)

if failed:
    fail_path = os.path.join(args.outdir, "failed_variants.tsv")
    pd.DataFrame(failed).to_csv(fail_path, sep="\t", index=False)
    print(f"{len(failed)} variants failed; see {fail_path}")

n_done = len(glob.glob(os.path.join(args.outdir, "scores_chunk_*.parquet")))
print(f"Done: {n_done}/{n_chunks} chunks written to {args.outdir}")
