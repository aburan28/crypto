#!/usr/bin/env python3
"""Figure for item 05.  Plots only values the Rust fit tool computed
(runs/fit_test.jsonl): per-cell geometric-mean paired ratios with 95 %
bootstrap CIs.

python3 -I plot.py  ->  fig_harvest.svg, fig_harvest.png
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
fit = [json.loads(l) for l in open(os.path.join(HERE, "runs", "fit_test.jsonl"))]
rat = [r for r in fit if "ratio" in r]
med = {r["medians"]: r for r in fit if "medians" in r}
ORDER = ["pipe-n13-l4-m3-a0", "pipe-n17-l6-m3-a1", "pipe-n19-l7-m3-a0", "pipe-n23-l8-m3-a0",
         "pipe-n29-l10-m3-a1", "pipe-n31-l10-m3-a0", "pipe-n17-l9-m2-a1", "pipe-n23-l12-m2-a0",
         "pipe-n29-l15-m2-a1"]


def short(c):
    p = c.split("-")
    return f"({p[1][1:]},{p[2][1:]}) {p[3]}"


def get(key, cell):
    return next(r for r in rat if r["ratio"] == key and r["cell"] == cell)


STY = {"I/G": ("#1f5fbf", "o", "I / G  frozen candidate: root, gpv < 0.9607 (verdict)"),
       "I/R": ("#2a9d6f", "^", "I / R  root, no gate (reference)"),
       "I/H": ("#d9822b", "s", "I / H  every node (reference)"),
       "AA": ("#888888", "x", "I / I′  A/A")}
fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [2.0, 1]})
for k, key in enumerate(STY):
    c, m, lab = STY[key]
    xs, ys, lo, hi = [], [], [], []
    for i, cell in enumerate(ORDER):
        r = get(key, cell)
        xs.append(i + (k - 1.5) * 0.17); ys.append(r["gm"])
        lo.append(r["gm"] - r["lo"]); hi.append(r["hi"] - r["gm"])
    a.errorbar(xs, ys, yerr=[lo, hi], color=c, marker=m, ms=5, lw=0, elinewidth=1.2, capsize=2.5, label=lab)
for i, cell in enumerate(ORDER):
    if med[cell]["g_harvests"] > 0:
        a.axvspan(i - 0.45, i + 0.45, color="#1f5fbf", alpha=0.07, lw=0)
a.axhline(1, color="k", lw=0.8, ls="--")
a.set_yscale("log")
a.set_xticks(range(len(ORDER)))
a.set_xticklabels([short(c) for c in ORDER], rotation=30, ha="right", fontsize=8.5)
a.set_ylabel("paired thread-CPU ratio  I / X   (> 1: X cheaper)")
a.set_title("Held-out targets (t = 300..): oracle solve cost per cell\n(shaded: G's gate fires; geometric mean, 95 % bootstrap CI)", fontsize=10)
a.legend(fontsize=8, loc="lower left")

OFF = {"pipe-n19-l7-m3-a0": (6, 9), "pipe-n29-l10-m3-a1": (6, 0), "pipe-n31-l10-m3-a0": (6, -9)}
for cell in ORDER:
    s = get("S_I/R", cell); w = get("W_I/R", cell)
    col = "#2a9d6f" if "-m2-" in cell else "#d9822b"
    b.errorbar([s["gm"]], [w["gm"]], yerr=[[w["gm"] - w["lo"]], [w["hi"] - w["gm"]]], color=col, marker="o", capsize=2.5)
    b.annotate(short(cell), (s["gm"], w["gm"]), fontsize=7.5, xytext=OFF.get(cell, (4, 3)), textcoords="offset points", va="center")
b.axhline(1, color="k", lw=0.8, ls="--"); b.axvline(1, color="#bbb", lw=0.8, ls=":")
b.set_xscale("log"); b.set_yscale("log")
b.set_xlabel("split reduction  splits(I) / splits(R)")
b.set_ylabel("word-operation ratio  I / R")
b.set_title("Mechanism of root harvest R (orange m=3, green m=2)\nm=2: fewer splits; m=3: same splits, smaller matrices", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_harvest.svg"))
fig.savefig(os.path.join(HERE, "fig_harvest.png"), dpi=150)
