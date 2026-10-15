"""Figures for round 2 item 02. Plots recorded values only (data/frozen_K.json,
results.json, data/census/*.json); computes nothing new.
Usage: python3 -I plots.py <item_dir>"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

D = Path(sys.argv[1])
fro = json.loads((D / "data" / "frozen_K.json").read_text())
res = json.loads((D / "results.json").read_text())
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"inc": "#2a78d6", "free": "#1baf7a", "def": "#eda100", "b": "#eb6834", "c": "#7a5bd6"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "svg.fonttype": "none"})
SIZES = ["n41", "n53", "n61"]


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# Figure 1: training sweep, mean counted cost vs K per arm (seeds 101-104).
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
for ax, n in zip(axes, SIZES):
    fz = fro["frozen"][n]
    for arm, lab in (("inc", "incumbent"), ("free", "+ free relations")):
        rows = sorted(((int(k), v) for k, v in fro["training"][n][arm].items()), key=lambda x: x[0])
        ks = [k for k, _ in rows]
        m = [v["c_total_mean"] / 1e6 for _, v in rows]
        sd = [v["c_total_sd"] / 1e6 for _, v in rows]
        off = -0.006 if arm == "inc" else 0.006
        xs = [k * (1 + off) for k in ks]
        ax.errorbar(xs, m, yerr=sd, fmt="-o", ms=4, lw=1.8, capsize=0, color=COL[arm], label=lab)
        kb = fz["K_inc" if arm == "inc" else "K_free"]
        vb = fro["training"][n][arm][str(kb)]["c_total_mean"] / 1e6
        ax.plot([kb * (1 + off)], [vb], marker="*", ms=13, color=COL[arm], mec=INK, mew=0.6)
        if arm == "free":
            for (k, v), x, y in zip(rows, xs, m):
                ax.annotate(f"{int(v['free_rank_mean'])}", (x, y), textcoords="offset points",
                            xytext=(3, -11), fontsize=7, color=COL["free"])
    ax.axvline(fz["K_def"], color=COL["def"], lw=1, ls="--", label="repository K")
    ax.set_title(f"n = {n[1:]}: K*_inc = {fz['K_inc']}, K*_free = {fz['K_free']}", color=INK, fontsize=9.5)
    ax.set_xlabel("orbit columns K")
    ax.set_ylabel("complete solve per target (M S3 root solves)")
    style(ax)
axes[0].legend(frameon=False, fontsize=8)
fig.suptitle("Training sweep (seeds 101-104, 8 targets; mean ± sd; star = frozen K*; green labels = free rank)",
             color=INK, fontsize=10)
fig.tight_layout()
fig.savefig(D / "fig_sweep.svg")
fig.savefig(D / "fig_sweep.png", dpi=150)

# Figure 2: held-out paired ratios.
fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
groups = SIZES + ["pooled"]
for ax, metric, title in zip(axes, ["counted_total", "wall_total"],
                             ["complete solve, counted unit", "complete solve, wall ms"]):
    for i, g in enumerate(groups):
        if g == "pooled":
            src = res["pooled"]
            items = [("a: free@K / inc@K at K*_inc", src["a_same_K_K_inc"], COL["inc"]),
                     ("a: same K at K*_free", src["a_same_K_K_free"], COL["free"]),
                     ("a: same K at repository K", src["a_same_K_K_def"], COL["def"]),
                     ("b: free@K*_free / inc@K*_inc", src["b_tuned_vs_tuned"], COL["b"]),
                     ("c: free@K*_free / inc@repository K", src["c_tuned_free_vs_inc_K_def"], COL["c"])]
        else:
            s = res["sizes"][g]
            K = s["K"]
            a = s["a_same_K"]
            items = [("a: free@K / inc@K at K*_inc", a[f"K_inc={K['K_inc']}"], COL["inc"]),
                     ("a: same K at K*_free", a[f"K_free={K['K_free']}"], COL["free"]),
                     ("a: same K at repository K", a[f"K_def={K['K_def']}"], COL["def"]),
                     ("b: free@K*_free / inc@K*_inc", s["b_tuned_vs_tuned"], COL["b"]),
                     ("c: free@K*_free / inc@repository K", s["c_tuned_free_vs_inc_K_def"], COL["c"])]
        for j, (lab, v, c) in enumerate(items):
            m = v[metric]
            x = i + (j - 2) * 0.14
            ax.errorbar(x, m["ratio"], yerr=[[m["ratio"] - m["ci95"][0]], [m["ci95"][1] - m["ratio"]]],
                        fmt="o", ms=5, color=c, elinewidth=2, capsize=0, label=lab if i == 0 else None)
    ax.axhline(1, color=INK2, lw=1)
    ax.set_xticks(range(len(groups)))
    lbl = []
    for g in groups:
        if g == "pooled":
            lbl.append("pooled")
        else:
            K = res["sizes"][g]["K"]
            lbl.append(f"n={g[1:]}\nK*inc={K['K_inc']} K*free={K['K_free']}\nrepo K={K['K_def']}")
    ax.set_xticklabels(lbl, fontsize=7.5)
    ax.set_title(title, color=INK, fontsize=9.5)
    ax.set_ylabel("ratio of means (95% bootstrap CI)")
    style(ax)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8, loc="lower center", ncol=3)
fig.suptitle("Held-out test, paired by seed (n=41, 53: 16 seeds; n=61: 12 seeds; 8 targets each)",
             color=INK, fontsize=10)
fig.tight_layout(rect=(0, 0.12, 1, 1))
fig.savefig(D / "fig_paired.svg")
fig.savefig(D / "fig_paired.png", dpi=150)
