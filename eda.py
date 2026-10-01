"""
eda.py
------
Exploratory plots for the full MPRA variant table (all 44,647 variants,
significant and not), after phyloP_scoring.py and conservation_level.py.

Input
-----
data/variants_all_conservation.tsv

Output (plots/)
---------------
phyloP447_distribution_all.png   phyloP histogram (phyloP_distribution.py)
phyloP447_by_clade.png           phyloP distribution within each conservation level
clade_counts.png                 variants per conservation level + % significant
clade_by_element.png             conservation-level mix of each element
genome_distribution.png          where the tested elements sit in hg38
element_counts.png               variants per element, significant vs not

Usage
-----
    python eda.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import phyloP_distribution

BASE_DIR = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction"
INPUT    = f"{BASE_DIR}/data/variants_all_conservation.tsv"
PLOT_DIR = f"{BASE_DIR}/plots"

SIG_P  = 0.05
LEVELS = ["human_only", "human_chimp", "great_apes", "apes", "primates", "mammals"]

# hg38 chromosome lengths (bp)
CHROM_LEN = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"
# ordinal blue ramp, light -> dark = less -> more deeply conserved
LEVEL_COLORS = dict(zip(LEVELS, ["#86b6ef", "#5598e7", "#2a78d6",
                                 "#1c5cab", "#104281", "#0d366b"]))


def style(ax, grid_axis="y"):
    ax.set_facecolor(SURFACE)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(colors=INK_2, labelsize=8, width=0.6)


def save(fig, name):
    out = f"{PLOT_DIR}/{name}"
    fig.savefig(out, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out}")


def load():
    df = pd.read_csv(INPUT, sep="\t")
    df["Chromosome"] = df["Chromosome"].astype(str).str.removeprefix("chr")
    df["conservation_level"] = pd.Categorical(df["conservation_level"],
                                              categories=LEVELS, ordered=True)
    df["significant"] = df["P-Value"] < SIG_P
    print(f"Loaded {len(df):,} variants ({df['significant'].sum():,} with P < {SIG_P}) "
          f"from {df['Element'].nunique()} elements")
    return df


# ── phyloP ──────────────────────────────────────────────────────────────────
def plot_phyloP_by_clade(df):
    present = [lv for lv in LEVELS if (df["conservation_level"] == lv).any()]
    data = [df.loc[df["conservation_level"] == lv, "phyloP447"].dropna() for lv in present]

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    for body, lv in zip(parts["bodies"], present):
        body.set_facecolor(LEVEL_COLORS[lv])
        body.set_edgecolor(SURFACE)
        body.set_alpha(0.9)
    ax.boxplot(data, widths=0.12, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=SURFACE, edgecolor=INK, linewidth=0.8),
               medianprops=dict(color=INK, linewidth=1.2),
               whiskerprops=dict(color=INK, linewidth=0.8), capprops=dict(linewidth=0))
    ax.axhline(0, color=INK_2, linewidth=0.6, zorder=1)
    ax.axhline(phyloP_distribution.SIG_CONSERVED, color=INK_2, linewidth=0.8,
               linestyle=":", zorder=1)
    ax.text(len(present) + 0.45, phyloP_distribution.SIG_CONSERVED,
            f"conserved\n(≥ {phyloP_distribution.SIG_CONSERVED})",
            fontsize=7, color=INK_2, va="center", ha="left")
    ax.set_xticks(range(1, len(present) + 1),
                  [f"{lv}\nn = {len(d):,}\nmedian {d.median():.2f}"
                   for lv, d in zip(present, data)], fontsize=7.5, color=INK)
    ax.set_xlim(0.4, len(present) + 0.6)
    # a few strongly accelerated sites reach -20; clip so the bulk is readable
    y_lo = -6
    n_below = int((df["phyloP447"] < y_lo).sum())
    ax.set_ylim(y_lo, 9)
    if n_below:
        ax.text(0.01, 0.01, f"{n_below} variants below {y_lo} not shown",
                transform=ax.transAxes, fontsize=7, color=INK_2)
    ax.set_ylabel("phyloP447", fontsize=9, color=INK)
    ax.set_title("phyloP447 within each conservation level", loc="left",
                 fontsize=10, color=INK)
    save(fig, "phyloP447_by_clade.png")


# ── clades ──────────────────────────────────────────────────────────────────
def plot_clade_counts(df):
    g = df.groupby("conservation_level", observed=False)
    counts = g.size()
    pct_sig = g["significant"].mean() * 100

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True,
                                 gridspec_kw=dict(width_ratios=[3, 2]))
    fig.patch.set_facecolor(SURFACE)
    y = np.arange(len(LEVELS))[::-1]
    colors = [LEVEL_COLORS[lv] for lv in LEVELS]

    style(a1, "x")
    a1.barh(y, counts.values, color=colors, height=0.7, edgecolor=SURFACE, linewidth=2)
    for yi, n in zip(y, counts.values):
        a1.text(n, yi, f"  {n:,} ({n / counts.sum():.0%})", va="center",
                fontsize=8, color=INK)
    a1.set_xlim(0, counts.max() * 1.3)
    a1.set_yticks(y, LEVELS, fontsize=8.5, color=INK)
    a1.set_xlabel("variants", fontsize=9, color=INK)
    a1.set_title("Variants per conservation level", loc="left", fontsize=10, color=INK)

    style(a2, "x")
    a2.barh(y, pct_sig.values, color=colors, height=0.7, edgecolor=SURFACE, linewidth=2)
    for yi, p, n in zip(y, pct_sig.values, counts.values):
        if n:
            a2.text(p, yi, f"  {p:.0f}%", va="center", fontsize=8, color=INK)
    overall = df["significant"].mean() * 100
    a2.axvline(overall, color=INK_2, linewidth=0.8, linestyle="--")
    a2.text(overall, len(LEVELS) - 0.45, f" all: {overall:.0f}%", fontsize=7,
            color=INK_2, va="bottom")
    a2.set_xlim(0, 100)
    a2.set_xlabel(f"% with MPRA P < {SIG_P}", fontsize=9, color=INK)
    a2.set_title("Share significant", loc="left", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, "clade_counts.png")


def plot_clade_by_element(df):
    tab = pd.crosstab(df["Element"], df["conservation_level"], normalize="index") \
            .reindex(columns=LEVELS, fill_value=0)
    n = df["Element"].value_counts()
    # order: most deeply conserved elements at the top
    score = (tab * np.arange(len(LEVELS))).sum(axis=1)
    tab = tab.loc[score.sort_values().index]

    fig, ax = plt.subplots(figsize=(8, 0.26 * len(tab) + 1.4))
    fig.patch.set_facecolor(SURFACE)
    style(ax, "x")
    left = np.zeros(len(tab))
    y = np.arange(len(tab))
    for lv in LEVELS:
        ax.barh(y, tab[lv].values * 100, left=left, color=LEVEL_COLORS[lv],
                height=0.75, edgecolor=SURFACE, linewidth=1, label=lv)
        left += tab[lv].values * 100
    ax.set_yticks(y, [f"{e} (n={n[e]:,})" for e in tab.index], fontsize=7.5, color=INK)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of element's variants", fontsize=9, color=INK)
    ax.legend(ncol=6, loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False,
              fontsize=7.5, labelcolor=INK, handlelength=1, columnspacing=1)
    ax.set_title("Conservation level mix per element", loc="left", fontsize=10,
                 color=INK, pad=22)
    save(fig, "clade_by_element.png")


# ── genome ──────────────────────────────────────────────────────────────────
def element_loci(df, merge_bp=1_000_000):
    """One row per locus: elements within `merge_bp` on the same chromosome
    (e.g. SORT1 / SORT1.2 / SORT1-flip, the four TERT runs) are merged."""
    el = (df.groupby("Element")
            .agg(chrom=("Chromosome", "first"), start=("Position", "min"),
                 end=("Position", "max"), n=("Position", "size"))
            .reset_index().sort_values(["chrom", "start"]))
    loci = []
    for chrom, g in el.groupby("chrom"):
        cur = None
        for r in g.itertuples():
            if cur and r.start - cur["end"] <= merge_bp:
                cur["elements"].append(r.Element)
                cur["end"] = max(cur["end"], r.end)
                cur["n"] += r.n
            else:
                cur = dict(chrom=chrom, start=r.start, end=r.end, n=r.n,
                           elements=[r.Element])
                loci.append(cur)
    return pd.DataFrame(loci)


def locus_label(elements):
    names = sorted(set(elements))
    # collapse replicate runs of one element (TERT-HEK, TERT-GBM, ... -> TERT ×4)
    stems = {e.split("-")[0].split(".")[0] for e in names}
    if len(stems) == 1 and len(names) > 1:
        return f"{stems.pop()} ×{len(names)}"
    return " / ".join(names)


def plot_genome(df):
    loci = element_loci(df)
    chroms = [c for c in CHROM_LEN if c != "Y"]
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.patch.set_facecolor(SURFACE)
    style(ax, None)
    ax.spines["left"].set_visible(False)

    ypos = {c: i for i, c in enumerate(chroms[::-1])}
    for c in chroms:
        ax.plot([0, CHROM_LEN[c] / 1e6], [ypos[c]] * 2, color=GRID,
                linewidth=6, solid_capstyle="round", zorder=1)

    size = lambda n: 20 + 180 * n / loci["n"].max()
    for chrom, g in loci.groupby("chrom"):
        g = g.sort_values("start")
        xs = ((g["start"] + g["end"]) / 2 / 1e6).to_numpy()
        for i, r in enumerate(g.itertuples()):
            x = xs[i]
            ax.scatter(x, ypos[chrom], s=size(r.n), color=BLUE, edgecolor=SURFACE,
                       linewidth=1, zorder=3)
            # labels sit above the marker; when two loci are close, push them
            # apart (left one right-aligned, right one left-aligned), and keep
            # labels near the chromosome start clear of the axis labels
            near_next = i + 1 < len(xs) and xs[i + 1] - x < 25
            near_prev = i > 0 and x - xs[i - 1] < 25
            ha, dx = ("right", 2) if near_next else ("left", -2) if near_prev or x < 15 \
                     else ("center", 0)
            ax.annotate(f"{locus_label(r.elements)} ({r.n:,})", (x, ypos[chrom]),
                        xytext=(dx, 7), textcoords="offset points",
                        ha=ha, va="bottom", fontsize=6.5, color=INK)

    ax.set_yticks(list(ypos.values()), [f"chr{c}" for c in ypos], fontsize=8, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("position (Mb, hg38)", fontsize=9, color=INK)
    ax.set_xlim(-3, max(CHROM_LEN.values()) / 1e6 + 5)
    ax.set_ylim(-0.8, len(chroms) - 0.2)
    ax.set_title(f"Genomic location of the {df['Element'].nunique()} MPRA elements "
                 f"({len(loci)} loci; marker size and label = variants)",
                 loc="left", fontsize=10, color=INK)
    save(fig, "genome_distribution.png")


def plot_element_counts(df):
    g = df.groupby("Element")["significant"].agg(["sum", "size"])
    g["not"] = g["size"] - g["sum"]
    g = g.sort_values("size")
    y = np.arange(len(g))

    fig, ax = plt.subplots(figsize=(7.5, 0.26 * len(g) + 1.4))
    fig.patch.set_facecolor(SURFACE)
    style(ax, "x")
    ax.barh(y, g["sum"], color=BLUE, height=0.75, edgecolor=SURFACE, linewidth=1,
            label=f"P < {SIG_P}")
    ax.barh(y, g["not"], left=g["sum"], color="#b7d3f6", height=0.75,
            edgecolor=SURFACE, linewidth=1, label=f"P ≥ {SIG_P}")
    for yi, (s, n) in enumerate(zip(g["sum"], g["size"])):
        ax.text(n, yi, f"  {n:,} ({s / n:.0%} sig.)", va="center", fontsize=7, color=INK)
    ax.set_yticks(y, g.index, fontsize=7.5, color=INK)
    ax.set_xlim(0, g["size"].max() * 1.3)
    ax.set_xlabel("variants", fontsize=9, color=INK)
    ax.legend(loc="lower right", frameon=False, fontsize=8, labelcolor=INK)
    ax.set_title("Variants per element", loc="left", fontsize=10, color=INK)
    save(fig, "element_counts.png")


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    df = load()

    print("\nphyloP447...")
    phyloP_distribution.main(INPUT, f"{PLOT_DIR}/phyloP447_distribution_all.png",
                             title_suffix=f" (all {len(df):,})")
    plot_phyloP_by_clade(df)

    print("\nConservation levels...")
    print(df["conservation_level"].value_counts().reindex(LEVELS).to_string())
    plot_clade_counts(df)
    plot_clade_by_element(df)

    print("\nGenomic distribution...")
    plot_genome(df)
    plot_element_counts(df)


if __name__ == "__main__":
    main()
