#!/usr/bin/env python3
"""Figure for round 3 item 02.  Plots only values the Rust fit tool computed
(runs/fit_test.jsonl, runs/train_diag.jsonl, runs/e2e_fit.jsonl).

python3 -I plot.py  ->  fig_gate.svg, fig_gate.png
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "runs")
load = lambda f: [json.loads(l) for l in open(os.path.join(R, f))]
test, train, e2e = load("fit_test.jsonl"), load("train_diag.jsonl"), load("e2e_fit.jsonl")
THETA = 331


def short(c):
    p = c.split("/")[0].split("-")
    return f"({p[1][1:]},{p[2][1:]}) {p[3]}"


trat = {(r["ratio"], r["cell"]): r for r in test if "ratio" in r}
tmed = {r["medians"]: r for r in test if "medians" in r}
cells = sorted({c.split("/")[0] for (_, c) in trat}, key=lambda c: tmed[c + "/bench"]["rank"])

fig, (a, b, c) = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw={"width_ratios": [2.2, 1.3, 0.8]})
STY = [("G0/C", "bench", "#1f5fbf", "o", "G0 / C, oracle bench (verdict)"),
       ("G0/C", "prod", "#7b3fbf", "D", "G0 / C, production call (first root)"),
       ("G0/R", "bench", "#2a9d6f", "^", "G0 / R, ungated root (reference)"),
       ("AA", "bench", "#888888", "x", "I / I′, A/A")]
for k, (key, mode, col, mk, lab) in enumerate(STY):
    xs, ys, lo, hi = [], [], [], []
    for i, cell in enumerate(cells):
        r = trat.get((key, f"{cell}/{mode}"))
        if r is None:
            continue
        xs.append(i + (k - 1.5) * 0.17); ys.append(r["gm"])
        lo.append(r["gm"] - r["lo"]); hi.append(r["hi"] - r["gm"])
    a.errorbar(xs, ys, yerr=[lo, hi], color=col, marker=mk, ms=5, lw=0, elinewidth=1.2, capsize=2.5, label=lab)
for i, cell in enumerate(cells):
    m = tmed[cell + "/bench"]
    if m["gpv"] >= 0.9607 and m["rank"] >= THETA:
        a.axvspan(i - 0.45, i + 0.45, color="#1f5fbf", alpha=0.08, lw=0)
a.axhline(1, color="k", lw=0.8, ls="--")
a.set_yscale("log")
from matplotlib.ticker import ScalarFormatter, NullFormatter
for ax_ in (a,):
    ax_.yaxis.set_major_formatter(ScalarFormatter()); ax_.yaxis.set_minor_formatter(NullFormatter())
a.set_yticks([0.9, 1.0, 1.1, 1.25, 1.5, 1.75, 2.0])
a.set_xticks(range(len(cells)))
a.set_xticklabels([f"{short(c)}\nrank {tmed[c + '/bench']['rank']:.0f}" for c in cells], rotation=35, ha="right", fontsize=8)
a.set_ylabel("paired CPU ratio  (> 1: arm cheaper than G0)")
a.set_title("Held-out test (t = 2000..), 8–12 targets per cell, 95 % bootstrap CI\nshaded: cells where C adds coverage over G0", fontsize=10)
a.legend(fontsize=8, loc="upper left")

tr = {(r["ratio"], r["cell"]): r for r in train if "ratio" in r}
trm = {r["medians"]: r for r in train if "medians" in r}
for cell, m in trm.items():
    w = tr[("W_I/R", cell)]
    col = "#d9822b" if "-m2-" in cell else "#1f5fbf"
    mk = "s" if m["gpv"] < 0.9607 else "o"
    b.errorbar([m["rank"]], [w["gm"]], yerr=[[w["gm"] - w["lo"]], [w["hi"] - w["gm"]]], color=col, marker=mk,
               ms=6, lw=0, elinewidth=1.2, capsize=2.5, mfc=col if mk == "o" else "white")
    b.annotate(short(cell), (m["rank"], w["gm"]), fontsize=7, xytext=(4, 3), textcoords="offset points")
b.axvline(THETA, color="#1f5fbf", ls=":", lw=1.2)
b.text(THETA * 1.03, 0.9, "θ = 331\n(frozen)", fontsize=8, color="#1f5fbf")
b.axhline(1, color="k", lw=0.8, ls="--")
b.set_xscale("log"); b.set_yscale("log")
b.yaxis.set_major_formatter(ScalarFormatter()); b.yaxis.set_minor_formatter(NullFormatter()); b.set_yticks([0.9, 1.0, 1.1, 1.25, 1.5])
b.xaxis.set_major_formatter(ScalarFormatter()); b.xaxis.set_minor_formatter(NullFormatter()); b.set_xticks([200, 300, 500, 1000, 1500])
b.set_xlabel("root basis rank (the gate's feature)")
b.set_ylabel("word ops, harvest off / ungated root harvest")
b.set_title("Training (t = 0..11): gain vs root rank\norange m = 2, blue m = 3; hollow square = G0 fires", fontsize=10)

pairs = [r for r in e2e if r["e2e"] == "gpv/unset"]
for i, r in enumerate(pairs):
    c.scatter([i] * len(r["pairs"]), r["pairs"], color="#7b3fbf", s=18)
    c.plot([i - 0.25, i + 0.25], [r["gm"]] * 2, color="#7b3fbf", lw=2)
    c.text(i, r["gm"] * 1.015, f"{r['gm']:.3f}", ha="center", fontsize=8)
c.axhline(1, color="k", lw=0.8, ls="--")
c.set_xticks(range(len(pairs)))
c.set_xticklabels([f"({r['n']},{r['ell']}) m2" for r in pairs], fontsize=8.5)
c.set_xlim(-0.6, len(pairs) - 0.4)
c.set_ylabel("end-to-end CPU, G0 (gpv) / C (new default)")
c.set_title("End to end, 3 replicate pairs\n(logs verified, identical relations)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_gate.svg"))
fig.savefig(os.path.join(HERE, "fig_gate.png"), dpi=150)
