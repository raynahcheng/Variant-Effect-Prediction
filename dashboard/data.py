"""Data loading for the dashboard.

Variants come from the MPRA saturation-mutagenesis table joined with model
scores (data/variants_alphagenome_with_CADD.tsv). The MPRA measurement
(Value / P-Value) is the ground truth every model is compared against.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "variants_alphagenome_with_CADD.tsv"

# Model name -> score column in the TSV. A model only shows up in the
# dashboard once its column exists, so to add Evo2 just add its column to the
# TSV (rename here if it is not called "evo2_score").
MODEL_COLUMNS = {
    "AlphaGenome": "raw_score",
    "CADD": "CADD_raw_score",
    "Evo2": "evo2_score",
}
# Scores are oriented so higher = larger predicted effect, which the ROC curve
# and Kendall's tau vs |Value| assume. AlphaGenome raw_score is signed
# (direction of effect), so its magnitude is used. Add Evo2's here if needed.
MODEL_TRANSFORMS = {"AlphaGenome": np.abs}

# Every variant is kept. For the ROC curves, variants with a significant MPRA
# effect (the CAGI5 / Kircher et al. 2019 cut-off) are positives and all
# others are negatives.
P_THRESHOLD = 1e-5
LABELS = ["functional", "neutral"]

# Dropdown label -> conservation_level value, most to least conserved.
CLADE_LEVELS = {
    "Mammals": "mammals",
    "Primates": "primates",
    "Apes": "apes",
    "Great Apes": "great_apes",
}
CLADES = ["All", *CLADE_LEVELS]
# Fixed color per species group, shared by every page's plots.
CLADE_COLORS = {"All": "#444444", "Mammals": "#636EFA", "Primates": "#EF553B",
                "Apes": "#00CC96", "Great Apes": "#AB63FA"}

README_DIR = "readmes"


def load_variants(path: Path = DATA_PATH) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": str})
    df.insert(0, "variant_id", df["Chromosome"] + ":" + df["Position"].astype(str)
              + ":" + df["Ref"] + ">" + df["Alt"] + ":" + df["Element"])
    df["chrom"] = "chr" + df["Chromosome"]
    df["abs_value"] = df["Value"].abs()
    df["label"] = np.where(df["P-Value"] < P_THRESHOLD, "functional", "neutral")
    for m, col in MODEL_COLUMNS.items():
        if col in df.columns:
            df[m] = MODEL_TRANSFORMS.get(m, lambda s: s)(df[col])
    return df


VARIANTS = load_variants()

# Models with results in the data; ALL_MODELS also includes pending ones.
ALL_MODELS = list(MODEL_COLUMNS)
MODELS = [m for m in ALL_MODELS if m in VARIANTS.columns]
MODEL_OPTIONS = [
    {"label": m if m in MODELS else f"{m} (no results yet)", "value": m,
     "disabled": m not in MODELS}
    for m in ALL_MODELS
]


def filter_by_clade(df: pd.DataFrame, clade: str) -> pd.DataFrame:
    if clade == "All":
        return df
    return df[df["conservation_level"] == CLADE_LEVELS[clade]]


def roc_curve(y_true: np.ndarray, score: np.ndarray):
    """Exact ROC curve for scores on any scale; NaN scores are dropped."""
    keep = ~np.isnan(score)
    y_true, score = y_true[keep], score[keep]
    order = np.argsort(-score, kind="mergesort")
    y, s = y_true[order], score[order]
    # One point per distinct threshold.
    last = np.r_[np.flatnonzero(np.diff(s)), len(s) - 1]
    tps = np.cumsum(y)[last]
    fps = (last + 1) - tps
    tpr = np.r_[0, tps / max(y.sum(), 1)]
    fpr = np.r_[0, fps / max((1 - y).sum(), 1)]
    return fpr, tpr


def auc(fpr: np.ndarray, tpr: np.ndarray) -> float:
    return float(np.trapezoid(tpr, fpr))


def kendall_tau(df: pd.DataFrame, model: str, target: str = "abs_value") -> tuple[float, int]:
    """Kendall's tau between a model score and MPRA |Value|, plus the n used."""
    sub = df[[target, model]].dropna()
    if len(sub) < 3:
        return np.nan, len(sub)
    return float(kendalltau(sub[target], sub[model])[0]), len(sub)
