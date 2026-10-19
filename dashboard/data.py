"""Mock data for the dashboard skeleton.

Replace the functions below with real loaders (CSV / parquet / DB) later.
"""
import numpy as np
import pandas as pd

MODELS = ["CADD", "AlphaGenome", "Evo2"]

# Dropdown groups from the whiteboard, broadest last.
CLADES = ["All", "Great Apes", "Apes", "Primates", "Mammals"]

# Which species belong to each clade filter.
SPECIES = {
    "Human": ["All", "Great Apes", "Apes", "Primates", "Mammals"],
    "Chimpanzee": ["All", "Great Apes", "Apes", "Primates", "Mammals"],
    "Gorilla": ["All", "Great Apes", "Apes", "Primates", "Mammals"],
    "Orangutan": ["All", "Great Apes", "Apes", "Primates", "Mammals"],
    "Gibbon": ["All", "Apes", "Primates", "Mammals"],
    "Macaque": ["All", "Primates", "Mammals"],
    "Marmoset": ["All", "Primates", "Mammals"],
    "Mouse": ["All", "Mammals"],
    "Dog": ["All", "Mammals"],
    "Cow": ["All", "Mammals"],
}

README_DIR = "readmes"


def species_for(clade: str) -> list[str]:
    return [s for s, groups in SPECIES.items() if clade in groups]



def make_variants(n: int = 600, seed: int = 7) -> pd.DataFrame:
    """Synthetic variant-level dataset with one score column per model."""
    rng = np.random.default_rng(seed)
    species = rng.choice(list(SPECIES), size=n)
    label = rng.choice(["pathogenic", "benign"], size=n, p=[0.35, 0.65])
    base = np.where(label == "pathogenic", 0.7, 0.3)
    df = pd.DataFrame(
        {
            "variant_id": [f"var_{i:05d}" for i in range(n)],
            "species": species,
            "chrom": rng.choice([f"chr{c}" for c in range(1, 23)], size=n),
            "position": rng.integers(1_000_000, 200_000_000, size=n),
            "label": label,
        }
    )
    for i, m in enumerate(MODELS):
        noise = 0.18 + 0.04 * i
        df[m] = np.clip(base + rng.normal(0, noise, size=n), 0, 1).round(4)
    return df


VARIANTS = make_variants()


def filter_by_clade(df: pd.DataFrame, clade: str) -> pd.DataFrame:
    return df[df["species"].isin(species_for(clade))]


def roc_curve(y_true: np.ndarray, score: np.ndarray, points: int = 100):
    thr = np.linspace(1, 0, points)
    pos = y_true == 1
    neg = ~pos
    tpr = np.array([(score[pos] >= t).mean() if pos.any() else 0 for t in thr])
    fpr = np.array([(score[neg] >= t).mean() if neg.any() else 0 for t in thr])
    return fpr, tpr


def auc(fpr: np.ndarray, tpr: np.ndarray) -> float:
    return float(np.trapezoid(tpr, fpr))
