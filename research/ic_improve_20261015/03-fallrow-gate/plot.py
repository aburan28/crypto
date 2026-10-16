#!/usr/bin/env python3
"""Figure for item 03. Plots only values computed by the Rust fit tool
(runs/fit_test.jsonl, runs/fit_diag.jsonl): per-cell geometric-mean paired
CPU ratios with 95 % bootstrap CIs.

python3 -I plot.py  ->  fig_gate.svg, fig_gate.png
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = lambda f: [json.loads(l) for l in open(os.path.join(HERE, "runs", f))]
test = [r for r in rows("fit_test.jsonl")]
diag = [r for r in rows("fit_diag.jsonl") if "ratio" in r]

ORDER = ["pipe-n13-l4-m3-a0", "pipe-n17-l6-m3-a1", "pipe-n19-l7-m3-a0", "pipe-n23-l8-m3-a0",
         "pipe-n29-l10-m3-a1", "pipe-n31-l10-m3-a0", "pipe-n17-l9-m2-a1", "pipe-n23-l12-m2-a0",
         "pipe-n29-l15-m2-a1", "barrel-n5-l3", "barrel-n7-l3", "barrel-n9-l3", "barrel-n11-l3", "barrel-n13-l3"]


def short(c):
    p = c.split("-")
    if p[0] == "barrel":
        return f"barrel n={p[1][1:]}"
    return f"({p[1][1:]},{p[2][1:]}) {p[3]}"


STY = {"I/G": ("#1f5fbf", "o", "I / G  gated candidate (verdict)"),
       "I/F": ("#d9822b", "s", "I / F  always augment (reference)"),
       "I/LR": ("#2a9d6f", "^", "I / LR  root linear tail only (reference)"),
       "AA": ("#888888", "x", "I / I′  A/A")}

fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [2.2, 1]})
for k, key in enumerate(["I/G", "I/F", "I/LR", "AA"]):
    c, m, lab = STY[key]
    xs, ys, lo, hi = [], [], [], []
    for i, cell in enumerate(ORDER):
        r = next((r for r in test if r["ratio"] == key and r["cell"] == cell), None)
        if r is None:
            continue
        xs.append(i + (k - 1.5) * 0.17); ys.append(r["gm"])
        lo.append(r["gm"] - r["lo"]); hi.append(r["hi"] - r["gm"])
    a.errorbar(xs, ys, yerr=[lo, hi], color=c, marker=m, ms=5, lw=0, elinewidth=1.2, capsize=2.5, label=lab)
for i, cell in enumerate(ORDER):
    r = next(r for r in test if r["ratio"] == "I/G" and r["cell"] == cell)
    if r["fires"] > 0:
        a.axvspan(i - 0.45, i + 0.45, color="#1f5fbf", alpha=0.07, lw=0)
a.axhline(1, color="k", lw=0.8, ls="--")
a.set_yscale("log")
a.set_xticks(range(len(ORDER)))
a.set_xticklabels([short(c) for c in ORDER], rotation=40, ha="right", fontsize=8)
a.set_ylabel("paired CPU ratio  I / X   (> 1: X cheaper)")
a.set_title("Test phase, held-out targets: whole-oracle cost per cell\n(shaded: gate fires, gpv ≥ 2.158; geometric mean, 95 % bootstrap CI)", fontsize=10)
a.legend(fontsize=8, loc="upper left")

DC = ["pipe-n13-l4-m3-a0", "pipe-n17-l6-m3-a1", "pipe-n19-l7-m3-a0", "pipe-n23-l12-m2-a0",
      "barrel-n7-l3", "barrel-n9-l3", "barrel-n11-l3"]
for cell in DC:
    s = next(r for r in diag if r["ratio"] == "Isplits/Fsplits" and r["cell"] == cell)
    f = next(r for r in diag if r["ratio"] == "I/Fsolve" and r["cell"] == cell)
    col = "#1f5fbf" if cell.startswith("barrel") else ("#2a9d6f" if "-m2-" in cell else "#d9822b")
    b.errorbar([s["gm"]], [f["gm"]], yerr=[[f["gm"] - f["lo"]], [f["hi"] - f["gm"]]], color=col, marker="o", capsize=2.5)
    b.annotate(short(cell), (s["gm"], f["gm"]), fontsize=7.5, xytext=(4, 4), textcoords="offset points")
b.axhline(1, color="k", lw=0.8, ls="--")
b.plot([1, 200], [1, 200], color="#bbbbbb", lw=0.8, ls=":")
b.text(30, 45, "cost falls as fast as splits", fontsize=7, color="#888", rotation=33)
b.set_xscale("log"); b.set_yscale("log")
b.set_xlabel("split reduction  splits(I) / splits(F)")
b.set_ylabel("solve-cost ratio  I / F_solve  (pass excluded)")
b.set_title("Diagnostic: why F loses on pipeline systems\n(orange m=3, green m=2 pipeline, blue barrel; 4 targets/cell)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_gate.svg"))
fig.savefig(os.path.join(HERE, "fig_gate.png"), dpi=150)
