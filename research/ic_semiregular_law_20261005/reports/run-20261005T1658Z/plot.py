#!/usr/bin/env python3
"""Plot recorded refutation degrees against n, with the registered windows and the
fixed-l fit. Layout only: every plotted value is read from the run records and from
runs/registered/predictions.jsonl; nothing is computed or re-run.

    plot.py OUT.svg   (run from research/ic_semiregular_law_20261005)
"""
import json, os, sys, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = ["registered", "registered-remote-c", "registered-remote-d", "registered-remote-e"]
FIT = {("x4", 4): 9, ("x4", 5): 10, ("x4", 6): 11, ("rr", 4): 6, ("rr", 5): 8, ("rr", 6): 8}
pred = {}
for l in open("runs/registered/predictions.jsonl"):
    r = json.loads(l); pred[(r["cell"], r["draw"], r["arm"])] = r["d_reg"]
read = {}
for s in SRC:
    d = f"runs/{s}"
    for f in os.listdir(d):
        if f.startswith("K") and f.endswith(".jsonl"):
            for l in open(os.path.join(d, f)):
                r = json.loads(l)
                if r.get("refuted_at", 0) > 0:
                    read.setdefault((r["cell"], r["draw"], r["arm"]), r["refuted_at"])
fig, axes = plt.subplots(2, 1, figsize=(10, 7.5))
colors = {4: "#1b7837", 5: "#2166ac", 6: "#b2182b"}
for ax, arm in zip(axes, ["x4", "rr"]):
    pts = collections.defaultdict(list)
    for (cell, draw, a), v in read.items():
        if a != arm or cell == "K1n61l4" or pred[(cell, draw, a)] == 0:
            continue
        pts[cell].append((draw, v, pred[(cell, draw, a)]))
    order = sorted(pts, key=lambda c: (int(c.split("l")[1]), int(c.split("n")[1].split("l")[0]), c[:2]))
    for x, cell in enumerate(order):
        ell = int(cell.split("l")[1])
        vs = sorted(pts[cell])
        dr = vs[0][2]
        ax.add_patch(plt.Rectangle((x - 0.38, dr - 0.15), 0.76, 1.3, color=colors[ell], alpha=0.18, lw=0))
        ax.plot([x - 0.38, x + 0.38], [FIT[(arm, ell)]] * 2, "--", color="black", lw=1.1)
        for i, (_, v, _) in enumerate(vs):
            ax.plot(x + (i - 1.5) * 0.14, v, "o", color=colors[ell], ms=6, mec="black", mew=0.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{c[:2]}\nn={c.split('n')[1].split('l')[0]}\nℓ={c.split('l')[1]}" for c in order], fontsize=8)
    for ell in (4, 5, 6):
        ax.plot([], [], "o", color=colors[ell], mec="black", label=f"ℓ = {ell}: one dot per draw")
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, color="grey", alpha=0.3, label="registered window [d_reg, d_reg + 1]"))
    ax.plot([], [], "--", color="black", label="fixed-ℓ fit")
    ax.set_title({"x4": "direct S₄ descent (x4)", "rr": "RR norm form (rr)"}[arm], fontsize=10)
    ax.set_ylabel("refutation degree D")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", ncol=2)
fig.suptitle("ic_semiregular_law_20261005 round 1: held-out draws (4 per cell), macaulay_dense, seed 20261005", fontsize=10)
fig.tight_layout()
fig.savefig(sys.argv[1])
