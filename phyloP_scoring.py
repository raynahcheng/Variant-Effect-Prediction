import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"
INPUT = f"{DATA_DIR}/GRCh38_significant.csv"
OUTPUT = f"{DATA_DIR}/variants_scored.tsv"

# ── Fetch phyloP447 scores from UCSC's REST API ─────────────────────────────
# Same cactus447way alignment used for species-level conservation checks in
# conservation_level.py, so the two are consistent by construction. The
# bigWig itself is ~10GB and downloads get truncated in this environment, so
# we query it remotely instead.
API_URL = "https://api.genome.ucsc.edu/getData/track"
GENOME = "hg38"
TRACK = "phyloP447wayBW"

session = requests.Session()

def fetch_score(chrom, pos_1based):
    """Returns phyloP447 score at a single base, or None if no data there."""
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    params = {
        "genome": GENOME, "track": TRACK,
        "chrom": chrom, "start": pos_1based - 1, "end": pos_1based,
    }
    for attempt in range(3):
        try:
            resp = session.get(API_URL, params=params, timeout=30)
            resp.raise_for_status()
            items = resp.json().get(TRACK, [])
            return items[0]["value"] if items else None
        except requests.RequestException:
            if attempt == 2:
                raise

df = pd.read_csv(INPUT)

positions = list({(row.Chromosome, row.Position) for row in df.itertuples()})
print(f"Fetching phyloP447 scores for {len(positions)} unique positions "
      f"across {len(df)} variants...")

with ThreadPoolExecutor(max_workers=10) as pool:
    scores = pool.map(lambda cp: fetch_score(*cp), positions)
    score_by_pos = dict(zip(positions, scores))

df["phyloP447"] = [score_by_pos[(row.Chromosome, row.Position)] for row in df.itertuples()]

df.to_csv(OUTPUT, sep="\t", index=False)
print("Saved.")
