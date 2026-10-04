import argparse
import pandas as pd
import requests

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"

parser = argparse.ArgumentParser(description="Annotate variants with RepeatMasker TE overlap.")
parser.add_argument("--input", default=f"{DATA_DIR}/variants_all_conservation.tsv")
parser.add_argument("--output", default=f"{DATA_DIR}/variants_all_te.tsv")
args = parser.parse_args()
INPUT = args.input
OUTPUT = args.output

# ── Fetch RepeatMasker annotations from UCSC's REST API ─────────────────────
# Variants cluster in ~20 short MPRA element windows, so querying each window
# is far cheaper than downloading the full rmsk table.
API_URL = "https://api.genome.ucsc.edu/getData/track"
GENOME = "hg38"
TRACK = "rmsk"
PAD = 1000

# Transposable element classes. Classes RepeatMasker marks uncertain (e.g.
# "DNA?", "LTR?") are still called TE; rep_class keeps the "?" so they can be
# filtered out downstream. Everything else in rmsk (Simple_repeat,
# Low_complexity, Satellite, small RNAs, Unknown) is "other_repeat".
TE_CLASSES = {"LINE", "SINE", "LTR", "DNA", "RC", "Retroposon"}

session = requests.Session()

def fetch_repeats(chrom, start, end):
    """Returns rmsk records overlapping [start, end) (0-based)."""
    params = {"genome": GENOME, "track": TRACK, "chrom": chrom, "start": start, "end": end}
    for attempt in range(3):
        try:
            resp = session.get(API_URL, params=params, timeout=30)
            resp.raise_for_status()
            return resp.json().get(TRACK, [])
        except requests.RequestException:
            if attempt == 2:
                raise

df = pd.read_csv(INPUT, sep="\t" if INPUT.endswith(".tsv") else ",")
df["chrom"] = df["Chromosome"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")

# Merge positions on each chromosome into windows (gaps >10kb start a new one).
windows = []
for chrom, positions in df.groupby("chrom")["Position"]:
    positions = sorted(positions.unique())
    lo = hi = positions[0]
    for p in positions[1:]:
        if p - hi > 10_000:
            windows.append((chrom, lo, hi))
            lo = p
        hi = p
    windows.append((chrom, lo, hi))
print(f"Fetching RepeatMasker records for {len(windows)} windows...")

repeats = []
for chrom, lo, hi in windows:
    repeats += fetch_repeats(chrom, max(0, lo - 1 - PAD), hi + PAD)

rmsk = pd.DataFrame(repeats).drop_duplicates(subset=["genoName", "genoStart", "genoEnd", "repName"])
print(f"Retrieved {len(rmsk)} repeat records.")

def annotate(chrom, pos):
    """Best-scoring repeat overlapping a 1-based position, or None."""
    if rmsk.empty:
        return None
    hits = rmsk[(rmsk["genoName"] == chrom) & (rmsk["genoStart"] < pos) & (rmsk["genoEnd"] >= pos)]
    return None if hits.empty else hits.loc[hits["swScore"].idxmax()]

unique_pos = df[["chrom", "Position"]].drop_duplicates()
records = {}
for chrom, pos in unique_pos.itertuples(index=False):
    hit = annotate(chrom, pos)
    if hit is None:
        records[(chrom, pos)] = ("non_repeat", None, None, None, None)
    else:
        status = "TE" if hit["repClass"].rstrip("?") in TE_CLASSES else "other_repeat"
        records[(chrom, pos)] = (status, hit["repClass"], hit["repFamily"],
                                 hit["repName"], hit["milliDiv"] / 10)

cols = ["te_status", "rep_class", "rep_family", "rep_name", "rep_divergence"]
ann = pd.DataFrame([records[(c, p)] for c, p in zip(df["chrom"], df["Position"])],
                   columns=cols, index=df.index)
df = pd.concat([df.drop(columns="chrom"), ann], axis=1)

df.to_csv(OUTPUT, sep="\t", index=False)

print("\nVariants per category:")
print(df["te_status"].value_counts().to_string())
print("\nPer element (variant counts):")
print(df.groupby(["Element", "te_status"]).size().unstack(fill_value=0).to_string())
print("\nTE classes/families hit:")
print(df[df["te_status"] != "non_repeat"].groupby(["rep_class", "rep_family", "rep_name"])
      ["Element"].agg(["size", lambda e: ",".join(sorted(set(e)))]).to_string())
print("Saved.")
