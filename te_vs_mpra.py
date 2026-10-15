"""
te_vs_mpra.py
-------------
Within each conservation level, compare transposable-element (TE) variants
with non-TE variants:

  1. how much of each level is TE (unique positions)
  2. MPRA effect size |Value| for TE vs non-TE
  3. agreement between MPRA |Value| and AlphaGenome |raw_score|
     (Kendall's tau), for TE and non-TE separately, both pooled and within
     element (elements differ in effect scale and TE content, so the pooled
     numbers are confounded by element)

non-TE = te_status non_repeat or other_repeat (simple / low-complexity
repeats, which are not transposable elements).

Outputs (plots/te_vs_mpra/)
---------------------------
  te_fraction_by_conservation.png
  abs_value_te_vs_nonte.png
  tau_te_vs_nonte.png
  scatter_abs_value_vs_abs_raw_score.png
  summary.tsv              per level x TE: n, median |Value|, pooled tau + CI
  within_element_tau.tsv   tau per element x level x TE
  mannwhitney.tsv          |Value| TE vs non-TE test per level

Usage
-----
    python te_vs_mpra.py
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kendalltau, mannwhitneyu

BASE = "/data/bentonm_shared/variant_effect_prediction/Variant-Effect-Prediction"
DATA_DIR = f"{BASE}/data"

parser = argparse.ArgumentParser(description="TE vs non-TE variants within conservation levels.")
parser.add_argument("--alphagenome", default=f"{DATA_DIR}/variants_alphagenome.tsv")
parser.add_argument("--plot-dir", default=f"{BASE}/plots/te_vs_mpra")
args = parser.parse_args()
PLOT_DIR = args.plot_dir

N_BOOT = 1000
MIN_N = 30   # groups / cells smaller than this are flagged or skipped

# Most to least conserved.
LEVELS = ["mammals", "primates", "apes", "great_apes", "human_chimp", "human_only"]
TE_GROUPS = ["non-TE", "TE"]

os.makedirs(PLOT_DIR, exist_ok=True)
rng = np.random.default_rng(0)

# ── Plot style (matches compare_ag_vs_base.py) ──────────────────────────────
TE_COLORS = {"non-TE": "#898781", "TE": "#4a3aa7"}   # neutral grey vs violet
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
FADED = 0.35

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


def fmt_p(p):
    return "p < 1e-4" if p < 1e-4 else f"p = {p:.2g}"


# ── Load ────────────────────────────────────────────────────────────────────
df = pd.read_csv(args.alphagenome, sep="\t", dtype={"Chromosome": str})
df = df.dropna(subset=["Value", "raw_score", "conservation_level", "te_status"])
df["te"] = np.where(df["te_status"] == "TE", "TE", "non-TE")
df["abs_value"] = df["Value"].abs()
df["abs_raw_score"] = df["raw_score"].abs()
levels = [lv for lv in LEVELS if lv in set(df["conservation_level"])]
print(f"{len(df):,} variants; levels: {', '.join(levels)}")


# ── 1. TE fraction per level (unique positions) ─────────────────────────────
pos = df.drop_duplicates(["Chromosome", "Position"])
counts = pd.crosstab(pos["conservation_level"], pos["te"]).reindex(levels).fillna(0)
chi2, chi_p, _, _ = chi2_contingency(counts[TE_GROUPS])
frac = counts.div(counts.sum(axis=1), axis=0)

fig, ax = plt.subplots(figsize=(7.5, 3.6))
y = np.arange(len(levels))[::-1]
left = np.zeros(len(levels))
for g in TE_GROUPS:
    ax.barh(y, frac[g] * 100, left=left, color=TE_COLORS[g], height=0.62,
            edgecolor="white", linewidth=2, label=g)
    left += frac[g].to_numpy() * 100
for yi, lv in zip(y, levels):
    ax.text(101, yi, f"{frac.loc[lv, 'TE']:.0%} TE   n={int(counts.loc[lv].sum()):,}",
            va="center", fontsize=8.5, color=INK_2)
ax.set_yticks(y, levels)
ax.set_xlim(0, 100)
ax.set_xlabel("% of unique variant positions")
ax.grid(False)
ax.set_title("TE share of variant positions by conservation level")
ax.text(0, 1.04, f"Most → least conserved, top to bottom. χ² {fmt_p(chi_p)}.",
        transform=ax.transAxes, fontsize=8.5, color=MUTED)
ax.legend(loc="upper left", bbox_to_anchor=(0, -0.18), ncol=2, fontsize=8.5)
save(fig, "te_fraction_by_conservation.png")


# ── 2. |Value| TE vs non-TE within each level ───────────────────────────────
mw_rows = []
fig, ax = plt.subplots(figsize=(7.5, 4.6))
width = 0.34
for i, lv in enumerate(levels):
    sub = df[df["conservation_level"] == lv]
    a = sub.loc[sub["te"] == "TE", "abs_value"]
    b = sub.loc[sub["te"] == "non-TE", "abs_value"]
    for j, g in enumerate(TE_GROUPS):
        vals = sub.loc[sub["te"] == g, "abs_value"]
        if vals.empty:
            continue
        x = i + (j - 0.5) * width
        bp = ax.boxplot(vals, positions=[x], widths=width * 0.85, showfliers=False,
                        patch_artist=True, medianprops={"color": "white", "lw": 2},
                        whiskerprops={"color": TE_COLORS[g]}, capprops={"color": TE_COLORS[g]})
        bp["boxes"][0].set(facecolor=TE_COLORS[g], edgecolor=TE_COLORS[g],
                           alpha=FADED if len(vals) < MIN_N else 1)
    if len(a) and len(b):
        p = mannwhitneyu(a, b).pvalue
        mw_rows.append({"conservation_level": lv, "n_TE": len(a), "n_nonTE": len(b),
                        "median_abs_value_TE": a.median(), "median_abs_value_nonTE": b.median(),
                        "mannwhitney_p": p})
        ax.annotate(fmt_p(p), (i, 1), xycoords=("data", "axes fraction"), ha="center",
                    va="bottom", fontsize=8.5, color=INK_2)
mw = pd.DataFrame(mw_rows)
mw.to_csv(f"{PLOT_DIR}/mannwhitney.tsv", sep="\t", index=False)
n_by = df.groupby(["conservation_level", "te"]).size()
ax.set_xticks(range(len(levels)),
              [f"{lv}\nTE {n_by.get((lv, 'TE'), 0):,}\nnon-TE {n_by.get((lv, 'non-TE'), 0):,}"
               for lv in levels], fontsize=8.5)
ax.set_ylabel("MPRA |Value|  (|log2 effect|)")
ax.set_title("MPRA effect size, TE vs non-TE variants", pad=34)
ax.text(0, 1.07, f"Boxes = IQR, white line = median, whiskers 1.5×IQR. Mann–Whitney p above. "
        f"Faded = n < {MIN_N}.", transform=ax.transAxes, fontsize=8.5, color=MUTED)
for g in TE_GROUPS:
    ax.bar(0, 0, color=TE_COLORS[g], label=g)
ax.legend(loc="upper right", fontsize=8.5)
save(fig, "abs_value_te_vs_nonte.png")


# ── 3. Kendall's tau |Value| vs |raw_score|, TE vs non-TE ───────────────────
def tau(g):
    return kendalltau(g["abs_value"], g["abs_raw_score"])[0]


rows = []
for lv in levels:
    for g in TE_GROUPS:
        sub = df[(df["conservation_level"] == lv) & (df["te"] == g)]
        if len(sub) < 3:
            continue
        idx = np.arange(len(sub))
        boot = [tau(sub.iloc[rng.choice(idx, len(sub))]) for _ in range(N_BOOT)]
        rows.append({"conservation_level": lv, "te": g, "n": len(sub),
                     "median_abs_value": sub["abs_value"].median(),
                     "kendall_tau": tau(sub),
                     "tau_p": kendalltau(sub["abs_value"], sub["abs_raw_score"])[1],
                     "tau_lo": np.nanpercentile(boot, 2.5),
                     "tau_hi": np.nanpercentile(boot, 97.5),
                     "low_n": len(sub) < MIN_N})
summary = pd.DataFrame(rows)
summary.to_csv(f"{PLOT_DIR}/summary.tsv", sep="\t", index=False)

# Within element: tau per (element, level, TE) cell, so variants are only
# ranked against others in the same element; n-weighted mean per level x TE,
# CI from resampling elements.
cells = []
for (el, lv, g), sub in df.groupby(["Element", "conservation_level", "te"]):
    if len(sub) < MIN_N:
        continue
    cells.append({"Element": el, "conservation_level": lv, "te": g,
                  "n": len(sub), "kendall_tau": tau(sub)})
cells = pd.DataFrame(cells)
cells.to_csv(f"{PLOT_DIR}/within_element_tau.tsv", sep="\t", index=False)


def wtau(c):
    return np.average(c["kendall_tau"], weights=c["n"])


within_rows = []
for (lv, g), c in cells.groupby(["conservation_level", "te"]):
    boot = [wtau(c.iloc[rng.choice(len(c), len(c))]) for _ in range(N_BOOT)]
    within_rows.append({"conservation_level": lv, "te": g, "n_elements": len(c),
                        "n": int(c["n"].sum()), "weighted_tau": wtau(c),
                        "tau_lo": np.percentile(boot, 2.5), "tau_hi": np.percentile(boot, 97.5)})
within = pd.DataFrame(within_rows)

pd.set_option("display.width", 200)
print("\nPooled tau (|Value| vs |raw_score|):")
print(summary.round(4).to_string(index=False))
print("\nWithin-element tau:")
print(within.round(4).to_string(index=False))
print("\n|Value| TE vs non-TE:")
print(mw.round(4).to_string(index=False))


def grouped_bars(ax, table, val, faded_col, label_fn):
    width = 0.36
    for j, g in enumerate(TE_GROUPS):
        t = table[table["te"] == g].set_index("conservation_level")
        for i, lv in enumerate(levels):
            if lv not in t.index:
                continue
            r = t.loc[lv]
            x = i + (j - 0.5) * width
            ax.bar(x, r[val], width=width * 0.92, color=TE_COLORS[g],
                   alpha=FADED if faded_col(r) else 1,
                   label=g if lv == t.index[0] else None)
            ax.errorbar(x, r[val], yerr=[[r[val] - r["tau_lo"]], [r["tau_hi"] - r[val]]],
                        fmt="none", ecolor=INK_2, elinewidth=1.1, capsize=3)
            ax.annotate(label_fn(r), (x, max(r["tau_hi"], r[val])), xytext=(0, 4),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=7, color=INK_2)
    ax.axhline(0, color=AXIS, lw=1)
    ax.margins(y=0.18)
    ax.set_xticks(range(len(levels)), levels, fontsize=9)


fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
grouped_bars(a1, summary, "kendall_tau", lambda r: r["low_n"], lambda r: f"n={r['n']:,}")
a1.set_ylabel("Kendall's tau, |Value| vs |raw_score| (95% CI)")
a1.set_title("Pooled within each level")
a1.text(0, 1.04, f"Faded = n < {MIN_N}.", transform=a1.transAxes, fontsize=8.5, color=MUTED)
grouped_bars(a2, within, "weighted_tau", lambda r: r["n_elements"] < 3,
             lambda r: f"{r['n_elements']} el")
a2.set_title("Within element (n-weighted mean)")
a2.text(0, 1.04, f"Cells with n ≥ {MIN_N}; CI over elements. Faded = fewer than 3 elements.",
        transform=a2.transAxes, fontsize=8.5, color=MUTED)
a2.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=TE_COLORS[g]) for g in TE_GROUPS],
           labels=TE_GROUPS, loc="upper right", fontsize=8.5)
fig.suptitle("Does AlphaGenome track MPRA effect size equally for TE and non-TE variants?",
             x=0.06, ha="left", fontweight="bold", fontsize=11, y=1.03)
save(fig, "tau_te_vs_nonte.png")


# ── 4. Scatter |Value| vs |raw_score|, rows = TE status, cols = level ───────
fig, axes = plt.subplots(2, len(levels), figsize=(3.2 * len(levels), 6.4),
                         sharex=True, sharey=True, squeeze=False)
xmax = df["abs_value"].quantile(0.995)
ymax = df["abs_raw_score"].quantile(0.995)
for r, g in enumerate(TE_GROUPS):
    for c, lv in enumerate(levels):
        ax = axes[r, c]
        ax.grid(True, axis="both")
        sub = df[(df["conservation_level"] == lv) & (df["te"] == g)]
        alpha = float(np.clip(30 / max(len(sub), 1), 0.04, 0.8))
        ax.scatter(sub["abs_value"], sub["abs_raw_score"], s=6, color=TE_COLORS[g],
                   alpha=alpha, linewidths=0, rasterized=True)
        if len(sub) >= 3:
            ax.text(0.97, 0.95, f"τ = {tau(sub):.2f}\nn = {len(sub):,}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK)
        if r == 0:
            ax.set_title(lv, pad=8)
        if c == 0:
            ax.set_ylabel(f"{g}\n|raw_score|")
        if r == 1:
            ax.set_xlabel("MPRA |Value|")
# Both scores are heavily right-skewed: symlog keeps zeros while spreading
# out the dense low end.
axes[0, 0].set_xscale("symlog", linthresh=0.01)
axes[0, 0].set_yscale("symlog", linthresh=0.001)
axes[0, 0].set_xlim(0, xmax)
axes[0, 0].set_ylim(0, ymax)
fig.suptitle("MPRA |Value| vs AlphaGenome |raw_score|", x=0.02, ha="left",
             fontweight="bold", fontsize=11, y=1.04)
fig.text(0.02, 1.0, "Symlog axes (linear below 0.01 / 0.001), clipped at the 99.5th percentile.",
         fontsize=8.5, color=MUTED)
save(fig, "scatter_abs_value_vs_abs_raw_score.png")
