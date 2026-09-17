import pyBigWig
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from scipy import stats

df = pd.read_csv("/home/chengr/Capstone/variant_effect_prediction/data/GRCh38_significant.csv")

BW_PATH = "/home/chengr/Capstone/variant_effect_prediction/data/cactus241way.phyloP.bw"   # ← update to your actual path
bw_phylop = pyBigWig.open(BW_PATH)

def fetch_score(bw, chrom, pos_1based):
    """Returns score at a single base; pos_1based → 0-based half-open [pos-1, pos)."""
    chrom = f"chr{chrom}"
    vals = bw.values(chrom, pos_1based - 1, pos_1based)
    return vals[0] if vals else None

df["phyloP241"] = df.apply(lambda r: fetch_score(bw_phylop, r.Chromosome, r.Position), axis=1)

df.to_csv("/home/chengr/Capstone/variant_effect_prediction/data/variants_scored.tsv", sep="\t", index=False)

bw_phylop.close()