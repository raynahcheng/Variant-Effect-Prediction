"""
compare_ag_vs_base.py
---------------------
Compare the MPRA ranking (GRCh38.csv, ranked by P-Value) with the AlphaGenome
ranking (variants_alphagenome.tsv, ranked by |raw_score|), stratified by
conservation level, to see whether AlphaGenome agrees better with the MPRA on
more conserved variants.

Only variants AlphaGenome has scored are compared. Within each conservation
group the variants are re-ranked, then:

  kendall_tau   Kendall's tau-b between MPRA rank and AlphaGenome rank
                (95% CI from bootstrap; 0 = no agreement)
  auroc         how well |raw_score| separates MPRA-significant
                (P < 0.05) from non-significant variants (0.5 = chance)

These pooled numbers are confounded by element (each conservation level is a
different mix of elements, whose P-Values and AlphaGenome scores sit on
different scales), so tau is also computed within each element and
conservation level, then combined across elements (n-weighted), and levels
are compared pairwise using only elements that have both.

Outputs
-------
  data/ranking_comparison.tsv          per-variant ranks
  plots/ag_vs_base/by_conservation.tsv summary table (pooled)
  plots/ag_vs_base/kendall_tau_by_conservation.png
  plots/ag_vs_base/auroc_by_conservation.png
  plots/ag_vs_base/within_element_tau.tsv               tau per element x level
  plots/ag_vs_base/within_element_by_conservation.tsv   weighted tau per level
  plots/ag_vs_base/within_element_paired.tsv            paired level comparisons
  plots/ag_vs_base/within_element_tau_by_conservation.png
  plots/ag_vs_base/within_element_lines.png

Usage
-----
    python compare_ag_vs_base.py
    python compare_ag_vs_base.py \
        --alphagenome data/old_conservation_levels/variants_alphagenome.tsv \
        --plot-dir plots/ag_vs_base_old_levels \
        --ranking-out data/old_conservation_levels/ranking_comparison.tsv
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, wilcoxon
from sklearn.metrics import roc_auc_score, roc_curve

BASE = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction"
DATA_DIR = f"{BASE}/data"

parser = argparse.ArgumentParser(description="Compare MPRA and AlphaGenome rankings by conservation.")
parser.add_argument("--alphagenome", default=f"{DATA_DIR}/variants_alphagenome.tsv",
                    help="Table with conservation_level and raw_score columns.")
parser.add_argument("--plot-dir", default=f"{BASE}/plots/ag_vs_base")
parser.add_argument("--ranking-out", default=f"{DATA_DIR}/ranking_comparison.tsv")
args = parser.parse_args()
PLOT_DIR = args.plot_dir

SIG_P = 0.05
N_BOOT = 1000
MIN_N = 30   # groups smaller than this are reported but flagged

# Most to least conserved.
LEVELS = ["mammals", "primates", "apes", "great_apes", "human_chimp", "human_only"]
KEY_COLS = ["Chromosome", "Position", "Ref", "Alt", "Element"]

os.makedirs(PLOT_DIR, exist_ok=True)
rng = np.random.default_rng(0)

# ── Plot style ──────────────────────────────────────────────────────────────
# One fixed colour per conservation level, used in every plot.
LEVEL_COLORS = {
    "mammals":     "#2a78d6",   # blue
    "primates":    "#eb6834",   # orange
    "apes":        "#1baf7a",   # aqua
    "great_apes":  "#e87ba4",   # magenta
    "human_chimp": "#4a3aa7",   # violet
    "human_only":  "#e34948",   # red
    "all":         "#898781",   # grey
}
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
FADED = 0.35   # alpha for groups too small to trust

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.titlepad": 22,
    "axes.labelcolor": INK_2,
    "axes.edgecolor": AXIS,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
})


def save(fig, name):
    path = f"{PLOT_DIR}/{name}"
    fig.savefig(path)
    plt.close(fig)
    print(f"Saved {path}")


def bar_by_level(ax, levels, values, lo, hi, faded, labels):
    """Bars coloured by conservation level, with CIs and value labels."""
    x = np.arange(len(levels))
    for xi, lv, v, f in zip(x, levels, values, faded):
        ax.bar(xi, v, width=0.62, color=LEVEL_COLORS[lv], alpha=FADED if f else 1)
    ax.errorbar(x, values, yerr=[values - lo, hi - values], fmt="none",
                ecolor=INK_2, elinewidth=1.2, capsize=4)
    for xi, v, h, l in zip(x, values, hi, lo):
        y = max(h, v) if v >= 0 else min(l, v)
        ax.annotate(f"{v:.2f}", (xi, y), xytext=(0, 6 if v >= 0 else -6),
                    textcoords="offset points", ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=9, color=INK)
    ax.axhline(0, color=AXIS, lw=1)
    ax.margins(y=0.15)   # room for the value labels
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)


# ── Load and join ───────────────────────────────────────────────────────────
base = pd.read_csv(f"{DATA_DIR}/GRCh38.csv", dtype={"Chromosome": str})
ag = pd.read_csv(args.alphagenome, sep="\t",
                 dtype={"Chromosome": str}, usecols=KEY_COLS + ["conservation_level", "raw_score"])
df = base.merge(ag, on=KEY_COLS, how="left", validate="one_to_one")

df = df.dropna(subset=["raw_score", "P-Value", "conservation_level"])
print(f"{len(df):,}/{len(base):,} variants have AlphaGenome scores")

df["abs_raw_score"] = df["raw_score"].abs()
df["significant"] = df["P-Value"] < SIG_P

# Rank 1 = strongest: smallest P-Value / largest |raw_score|.
df["base_rank"] = df["P-Value"].rank(method="average")
df["ag_rank"] = df["abs_raw_score"].rank(method="average", ascending=False)
df["base_rank_in_group"] = df.groupby("conservation_level")["P-Value"].rank(method="average")
df["ag_rank_in_group"] = df.groupby("conservation_level")["abs_raw_score"].rank(
    method="average", ascending=False)


# ── Metrics per conservation group ──────────────────────────────────────────
def tau(g):
    return kendalltau(g["base_rank_in_group"], g["ag_rank_in_group"])


def metrics(g):
    t, t_p = tau(g)
    idx = np.arange(len(g))
    boot = [tau(g.iloc[rng.choice(idx, len(g))])[0] for _ in range(N_BOOT)]
    has_both = g["significant"].nunique() == 2
    return {
        "n": len(g),
        "n_significant": int(g["significant"].sum()),
        "kendall_tau": t,
        "tau_lo": np.nanpercentile(boot, 2.5),
        "tau_hi": np.nanpercentile(boot, 97.5),
        "tau_p": t_p,
        "auroc": roc_auc_score(g["significant"], g["abs_raw_score"]) if has_both else np.nan,
        "low_n": len(g) < MIN_N,
    }


rows = []
for level in LEVELS + ["all"]:
    g = df if level == "all" else df[df["conservation_level"] == level]
    if len(g) < 3:
        continue
    if level == "all":
        # Rank across all groups together rather than within each group.
        g = g.assign(base_rank_in_group=g["base_rank"], ag_rank_in_group=g["ag_rank"])
    rows.append({"conservation_level": level, **metrics(g)})
summary = pd.DataFrame(rows)


# ── Save ────────────────────────────────────────────────────────────────────
df.sort_values("base_rank").to_csv(args.ranking_out, sep="\t", index=False)
summary.to_csv(f"{PLOT_DIR}/by_conservation.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
print(summary.round(4).to_string(index=False))

faded = summary["low_n"].to_numpy()
labels = [f"{r.conservation_level}\nn={r.n:,}" for r in summary.itertuples()]

fig, ax = plt.subplots(figsize=(7.5, 4.5))
bar_by_level(ax, summary["conservation_level"], summary["kendall_tau"].to_numpy(),
             summary["tau_lo"].to_numpy(), summary["tau_hi"].to_numpy(), faded, labels)
ax.set_ylabel("Kendall's tau (95% CI)")
ax.set_title("MPRA rank vs AlphaGenome rank, pooled")
ax.text(0, 1.02, f"Most → least conserved. Faded = fewer than {MIN_N} variants.",
        transform=ax.transAxes, fontsize=8.5, color=MUTED)
save(fig, "kendall_tau_by_conservation.png")

# ROC curves: |raw_score| as a classifier for MPRA P < SIG_P.
fig, ax = plt.subplots(figsize=(6, 5.5))
ax.grid(True, axis="both")
ax.plot([0, 1], [0, 1], ls="--", color=AXIS, lw=1)
ax.text(0.62, 0.55, "chance", color=MUTED, fontsize=8.5, rotation=40)
for r in summary.itertuples():
    lv = r.conservation_level
    g = df if lv == "all" else df[df["conservation_level"] == lv]
    if g["significant"].nunique() < 2:
        continue
    fpr, tpr, _ = roc_curve(g["significant"], g["abs_raw_score"])
    ax.plot(fpr, tpr, color=LEVEL_COLORS[lv], lw=2,
            ls="--" if lv == "all" else "-", alpha=FADED if r.low_n else 1,
            label=f"{lv}  AUC {r.auroc:.2f}  (n={r.n:,})")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")
ax.set_xlabel("False positive rate")
ax.set_ylabel("True positive rate")
ax.set_title(f"|raw_score| separating MPRA P < {SIG_P} variants", pad=10)
ax.legend(loc="lower right", fontsize=8.5)
save(fig, "auroc_by_conservation.png")


# ── Within-element comparison ───────────────────────────────────────────────
# Elements differ in both their MPRA P-Value distribution and the scale of
# their AlphaGenome scores, and each conservation level is drawn from a
# different mix of elements, so the pooled numbers above are confounded by
# element. Here tau is computed within each (element, conservation level)
# cell, so variants are only ever ranked against others in the same element.
cells = []
for (element, level), g in df.groupby(["Element", "conservation_level"]):
    if len(g) < MIN_N:
        continue
    t, t_p = kendalltau(g["P-Value"], -g["abs_raw_score"])
    cells.append({"Element": element, "conservation_level": level,
                  "n": len(g), "kendall_tau": t, "tau_p": t_p})
cells = pd.DataFrame(cells)
cells.to_csv(f"{PLOT_DIR}/within_element_tau.tsv", sep="\t", index=False)

tau_table = cells.pivot(index="Element", columns="conservation_level", values="kendall_tau")
tau_table = tau_table[[lv for lv in LEVELS if lv in tau_table.columns]]
print("\nKendall's tau within each element (cells with n >= %d):" % MIN_N)
print(tau_table.round(2).to_string())


def weighted_tau(c):
    return np.average(c["kendall_tau"], weights=c["n"])


# Per conservation level: n-weighted mean tau across elements. CI from
# resampling elements, so it reflects element-to-element variability.
within_rows = []
for level in LEVELS:
    c = cells[cells["conservation_level"] == level]
    if c.empty:
        continue
    boot = [weighted_tau(c.iloc[rng.choice(len(c), len(c))]) for _ in range(N_BOOT)]
    within_rows.append({
        "conservation_level": level,
        "n_elements": len(c),
        "n": int(c["n"].sum()),
        "weighted_tau": weighted_tau(c),
        "tau_lo": np.percentile(boot, 2.5),
        "tau_hi": np.percentile(boot, 97.5),
        "elements": ",".join(c["Element"]),
    })
within = pd.DataFrame(within_rows)
within.to_csv(f"{PLOT_DIR}/within_element_by_conservation.tsv", sep="\t", index=False)
print("\nWithin-element tau by conservation level:")
print(within.drop(columns="elements").round(4).to_string(index=False))

# Paired comparisons: only elements that have both levels, so each element is
# its own control. Wilcoxon signed-rank test on the per-element differences.
paired_rows = []
present = [lv for lv in LEVELS if lv in tau_table.columns]
for i, a in enumerate(present):
    for b in present[i + 1:]:
        pair = tau_table[[a, b]].dropna()
        if len(pair) < 3:
            continue
        diff = pair[a] - pair[b]
        paired_rows.append({
            "more_conserved": a,
            "less_conserved": b,
            "n_elements": len(pair),
            "frac_more_conserved_higher": (diff > 0).mean(),
            "mean_tau_diff": diff.mean(),
            "wilcoxon_p": wilcoxon(diff).pvalue if (diff != 0).any() else np.nan,
        })
paired = pd.DataFrame(paired_rows)
paired.to_csv(f"{PLOT_DIR}/within_element_paired.tsv", sep="\t", index=False)
print("\nPaired within-element comparisons:")
print(paired.round(4).to_string(index=False))

# Plot: weighted within-element tau per conservation level.
labels = [f"{r.conservation_level}\n{r.n_elements} element{'s' * (r.n_elements != 1)}\nn={r.n:,}"
          for r in within.itertuples()]
fig, ax = plt.subplots(figsize=(7.5, 4.8))
bar_by_level(ax, within["conservation_level"], within["weighted_tau"].to_numpy(),
             within["tau_lo"].to_numpy(), within["tau_hi"].to_numpy(),
             (within["n_elements"] < 3).to_numpy(), labels)
ax.set_ylabel("Kendall's tau within element\n(n-weighted mean, 95% CI over elements)")
ax.set_title("MPRA rank vs AlphaGenome rank, within element")
ax.text(0, 1.02, "Most → least conserved. Faded = fewer than 3 elements.",
        transform=ax.transAxes, fontsize=8.5, color=MUTED)
save(fig, "within_element_tau_by_conservation.png")

# Plot: each element's tau across conservation levels.
fig, ax = plt.subplots(figsize=(7.5, 5))
xpos = {lv: i for i, lv in enumerate(present)}
for element, row in tau_table.iterrows():
    row = row.dropna()
    if len(row) < 2:
        continue
    ax.plot([xpos[lv] for lv in row.index], row.values, color=AXIS, lw=1, zorder=1)
for lv in present:
    vals = tau_table[lv].dropna()
    ax.scatter(np.full(len(vals), xpos[lv]), vals, s=28, color=LEVEL_COLORS[lv],
               edgecolor="white", linewidth=0.8, zorder=2)
w = within[within["n_elements"] >= 3]
ax.plot([xpos[lv] for lv in w["conservation_level"]], w["weighted_tau"], color=INK,
        lw=2, marker="D", ms=7, zorder=3, label="n-weighted mean (levels with ≥ 3 elements)")
ax.axhline(0, color=AXIS, lw=1)
ax.set_xticks(range(len(present)))
ax.set_xticklabels(present, fontsize=9)
ax.set_ylabel("Kendall's tau within element")
ax.set_title("Each element's agreement by conservation level")
ax.text(0, 1.02, "One dot per element; grey lines join the same element. Most → least conserved.",
        transform=ax.transAxes, fontsize=8.5, color=MUTED)
ax.legend(loc="upper right", fontsize=8.5)
save(fig, "within_element_lines.png")
