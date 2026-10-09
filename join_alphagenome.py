"""
join_alphagenome.py
-------------------
Add an AlphaGenome raw_score column to variants_all_te.tsv and write
data/variants_alphagenome.tsv.

AlphaGenome gives each variant one RNA_SEQ raw_score (log2 fold change, alt vs
ref) per gene in the 1MB window per track. These are collapsed to a single
raw_score per variant: the one with the largest absolute value (sign kept),
i.e. the variant's strongest predicted expression effect. Variants not yet
scored get NaN, so this can be rerun while alphagenome_scoring.py is running.

Usage
-----
    python join_alphagenome.py
"""

import glob
import os

import pandas as pd

DATA_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction/data"
INPUT = f"{DATA_DIR}/variants_all_te.tsv"
SCORES = f"{DATA_DIR}/alphagenome_scores_RNA_SEQ"
OUTPUT = f"{DATA_DIR}/variants_alphagenome.tsv"

KEY_COLS = ["Chromosome", "Position", "Ref", "Alt"]

paths = sorted(glob.glob(os.path.join(SCORES, "scores_chunk_*.parquet")))
print(f"Reading {len(paths)} score chunks from {SCORES}")

best = []
for path in paths:
    s = pd.read_parquet(path, columns=KEY_COLS + ["raw_score"])
    idx = s["raw_score"].abs().groupby([s[c] for c in KEY_COLS]).idxmax()
    best.append(s.loc[idx])
best = pd.concat(best, ignore_index=True)
best["Chromosome"] = best["Chromosome"].astype(str)

variants = pd.read_csv(INPUT, sep="\t", dtype={"Chromosome": str})
out = variants.merge(best, on=KEY_COLS, how="left")
assert len(out) == len(variants)

out.to_csv(OUTPUT, sep="\t", index=False)
print(f"Wrote {OUTPUT}: {out['raw_score'].notna().sum():,}/{len(out):,} rows scored")
