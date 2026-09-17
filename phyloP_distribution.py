import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import stats

# ── user settings ──────────────────────────────────────────────────────────
INPUT_FILE  = "/home/chengr/Capstone/variant_effect_prediction/data/variants_scored.tsv"   # path to your scored variant file
SCORE_COL   = "phyloP241"             # column containing PhyloP scores
SEP         = "\t"                    # "\t" for TSV, "," for CSV
OUTPUT_PNG  = "phyloP241_distribution.png"

# Significance thresholds (Zoonomia paper: q ≤ 0.05 FDR → phyloP ≥ 2.27)
SIG_CONSERVED   =  2.27
SIG_ACCELERATED = -2.0
# ── end settings ───────────────────────────────────────────────────────────

def classify(score):
    if score >= SIG_CONSERVED:
        return "Conserved"
    elif score <= SIG_ACCELERATED:
        return "Accelerated"
    else:
        return "Neutral"


def main():
    # ── load data ──────────────────────────────────────────────────────────
    df = pd.read_csv(INPUT_FILE, sep=SEP)
    assert SCORE_COL in df.columns, (
        f"Column '{SCORE_COL}' not found. Available: {list(df.columns)}"
    )

    # drop rows where score is missing (sites with no alignment coverage)
    n_before = len(df)
    df = df.dropna(subset=[SCORE_COL])
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"Note: dropped {n_dropped} rows with missing PhyloP scores "
              f"(no alignment coverage at those positions).")

    scores = df[SCORE_COL].astype(float)
    df["category"] = scores.apply(classify)
    counts = df["category"].value_counts()
    total  = len(df)

    print(f"\nVariant counts (n={total}):")
    for cat in ["Conserved", "Neutral", "Accelerated"]:
        n = counts.get(cat, 0)
        print(f"  {cat:12s}: {n:5d}  ({n/total*100:.1f}%)")
    print(f"  Median phyloP : {scores.median():.3f}")
    print(f"  Mean phyloP   : {scores.mean():.3f}")

    # ── colours (CVD-safe: no red/green pair) ──────────────────────────────
    pal = {
        "Conserved":   "#2166AC",   # blue
        "Neutral":     "#999999",   # grey
        "Accelerated": "#D6604D",   # orange-red
    }

    # ── figure ─────────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.family": "sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    s_min, s_max = scores.min(), scores.max()
    bins = np.linspace(s_min - 0.3, s_max + 0.3, 65)

    # stacked histogram by category
    for cat in ["Accelerated", "Neutral", "Conserved"]:   # draw order
        subset = df.loc[df["category"] == cat, SCORE_COL].astype(float)
        n_cat  = counts.get(cat, 0)
        label  = f"{cat}  (n={n_cat}, {n_cat/total*100:.1f}%)"
        ax.hist(subset, bins=bins, color=pal[cat], alpha=0.80,
                label=label, zorder=3)

    # KDE overlay on twin axis
    kde_x = np.linspace(s_min - 1.5, s_max + 1.5, 600)
    kde   = stats.gaussian_kde(scores, bw_method=0.25)
    kde_y = kde(kde_x)

    ax2 = ax.twinx()
    ax2.plot(kde_x, kde_y, color="black", lw=1.4, ls="--",
             zorder=5, label="KDE")
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("Density", fontsize=8)
    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.3f"))
    ax2.tick_params(axis="y", labelsize=7)
    ax2.spines["top"].set_visible(False)

    # threshold vertical lines
    ymax = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else total * 0.6

    for xval, col, label_text, ha in [
        (SIG_CONSERVED,   pal["Conserved"],   f"phyloP ≥ {SIG_CONSERVED}",  "left"),
        (SIG_ACCELERATED, pal["Accelerated"], f"phyloP ≤ {SIG_ACCELERATED}", "right"),
    ]:
        ax.axvline(xval, color=col, lw=1.5, ls=":", zorder=6)
        offset = 0.12 if ha == "left" else -0.12
        ax.text(xval + offset, ymax * 0.97, label_text,
                color=col, fontsize=7, va="top", ha=ha)

    # vertical line at 0 (neutral reference)
    ax.axvline(0, color="#BBBBBB", lw=0.8, ls="-", zorder=2)

    # summary stats annotation
    stats_text = (
        f"n = {total}\n"
        f"median = {scores.median():.2f}\n"
        f"mean = {scores.mean():.2f}"
    )
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
            fontsize=7, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", lw=0.6))

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=7.5, frameon=False,
              loc="upper left", title="Category  (FDR 5%)", title_fontsize=7.5)

    ax.set_xlabel("PhyloP score (241-way Cactus / Zoonomia)", fontsize=9)
    ax.set_ylabel("Variant count", fontsize=9)
    ax.set_title(
        "PhyloP 241-way conservation scores across variants",
        fontsize=9, loc="left"
    )
    ax.tick_params(labelsize=7)
    ax.set_xlim(s_min - 1.0, s_max + 1.0)
    ax.margins(y=0.06)
    ax.grid(axis="y", lw=0.4, alpha=0.35, zorder=0)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=180, bbox_inches="tight")
    print(f"\nFigure saved to: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()