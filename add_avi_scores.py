#!/usr/bin/env python3
import sys
from pathlib import Path

import pandas as pd

KEY = ["Chromosome", "Position", "Ref", "Alt", "Element"]
SCORE_COL = "AVI_SCORE"
PLACEHOLDER = "0.0"


def main(variants_path, atlas_path, out_path=None):
    variants_path, atlas_path = Path(variants_path), Path(atlas_path)
    for p in (variants_path, atlas_path):
        if not p.exists():
            sys.exit(f"ERROR: file not found: {p}")
    out_path = Path(out_path) if out_path else variants_path.with_name(
        variants_path.stem + "_with_AVI.tsv"
    )

    # Read everything as strings so no value is reformatted (e.g. "0.0" -> "0",
    # chromosome "X" mixed with ints, long floats losing digits) when written back.
    variants = pd.read_csv(variants_path, sep="\t", dtype=str, keep_default_na=False)
    atlas = pd.read_csv(atlas_path, dtype=str, keep_default_na=False)

    for name, df in (("variants", variants), ("atlas", atlas)):
        missing = [c for c in KEY if c not in df.columns]
        if missing:
            sys.exit(f"ERROR: {name} file is missing key column(s): {missing}")
    if SCORE_COL not in atlas.columns:
        sys.exit(f"ERROR: atlas file has no {SCORE_COL} column")
    if SCORE_COL in variants.columns:
        sys.exit(f"ERROR: variants file already has a {SCORE_COL} column")

    # Safety: a duplicated key would silently multiply rows in the merge.
    if atlas.duplicated(KEY).any():
        sys.exit("ERROR: duplicate keys in atlas file; matching would be ambiguous")
    if variants.duplicated(KEY).any():
        sys.exit("ERROR: duplicate keys in variants file; matching would be ambiguous")

    n_before = len(variants)
    merged = variants.merge(atlas[KEY + [SCORE_COL]], on=KEY, how="left", validate="1:1")
    assert len(merged) == n_before, "row count changed during merge"

    # Atlas rows that never found a home in the variants file (sanity check).
    unused = len(atlas) - merged[SCORE_COL].notna().sum()

    no_score = merged[SCORE_COL].isna() | (merged[SCORE_COL].str.strip() == "")
    merged[SCORE_COL] = merged[SCORE_COL].where(~no_score, PLACEHOLDER)

    merged.to_csv(out_path, sep="\t", index=False)

    print(f"Variants rows:            {n_before}")
    print(f"Atlas rows:               {len(atlas)}")
    print(f"Rows given an AVI score:  {n_before - no_score.sum()}")
    print(f"Rows given placeholder:   {no_score.sum()} ({PLACEHOLDER})")
    if unused:
        print(f"WARNING: {unused} atlas rows had no matching variant row")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    if len(sys.argv) in (3, 4):
        main(*sys.argv[1:])
    elif len(sys.argv) == 1:
        main(
            r"C:\Users\evilt\PycharmProjects\Capstone\variants_alphagenome_with_CADD.tsv",
            r"C:\Users\evilt\PycharmProjects\Capstone\atlas_successful_lookups.csv",
            r"C:\Users\evilt\PycharmProjects\Capstone\variants_alphagenome_with_CADD_with_AVI.tsv",
        )
    else:
        sys.exit(__doc__)