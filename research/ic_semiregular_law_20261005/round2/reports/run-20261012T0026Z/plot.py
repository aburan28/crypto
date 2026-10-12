#!/usr/bin/env python3
"""Plot round 2's recorded refutation degrees per cell against the registered prediction
and the fixed-l fit. Layout only: values come from runs/registered/K*.jsonl and
runs/registered/predictions.jsonl; nothing is computed or re-run.

    plot.py OUT.svg   (run from the round2 directory)
"""
import json, os, sys, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXCL = {("K1n59l4", 0, "rr"), ("K1n59l4", 3, "rr")}
FIT = {3: 8, 4: 9, 5: 10, 6: 11}
pred = {}
for l in open("runs/registered/predictions.jsonl"):
    r = json.loads(l); pred[(r["cell"], r["draw"], r["arm"])] = r["predicted"]
read = {}
for f in os.listdir("runs/registered"):
    if f.startswith("K") and f.endswith(".jsonl"):
        for l in open(os.path.join("runs/registered", f)):
            r = json.loads(l)
            if r.get("refuted_at", 0) > 0:
                read.setdefault((r["cell"], r["draw"], r["arm"]), r["refuted_at"])
colors = {3: "#762a83", 4: "#1b7837", 5: "#2166ac", 6: "#b2182b"}
fig, axes = plt.subplots(2, 1, figsize=(10, 7.5))
for ax, arm in zip(axes, ["x4", "rr"]):
    pts = collections.defaultdict(list)
    for k, v in read.items():
        if k[2] == arm and k not in EXCL:
            pts[k[0]].append((k[1], v, pred[k]))
    order = sorted(pts, key=lambda c: (int(c.split("l")[1]), int(c.split("n")[1].split("l")[0]), c[:2]))
    for x, cell in enumerate(order):
        ell = int(cell.split("l")[1]); vs = sorted(pts[cell]); p = vs[0][2]
        ax.plot([x - 0.38, x + 0.38], [p, p], "-", color=colors[ell], lw=3, alpha=0.45)
        if arm == "x4":
            ax.plot([x - 0.38, x + 0.38], [FIT[ell]] * 2, "--", color="black", lw=1.1)
        for i, (_, v, _) in enumerate(vs):
            ax.plot(x + (i - 1.5) * 0.14, v, "o", color=colors[ell], ms=6, mec="black", mew=0.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{c[:2]}\nn={c.split('n')[1].split('l')[0]}\nℓ={c.split('l')[1]}" for c in order], fontsize=8)
    for ell in sorted(colors):
        ax.plot([], [], "o", color=colors[ell], mec="black", label=f"ℓ = {ell}: one dot per draw")
    ax.plot([], [], "-", color="grey", lw=3, alpha=0.45, label="registered prediction" + (" (d_aff)" if arm == "x4" else " (d_red rule)"))
    if arm == "x4":
        ax.plot([], [], "--", color="black", label="fixed-ℓ fit (ℓ + 5)")
    ax.set_title({"x4": "direct S₄ descent (x4): prediction A", "rr": "RR norm form (rr): prediction B"}[arm], fontsize=10)
    ax.set_ylabel("refutation degree D"); ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=7, loc="upper left", ncol=2)
fig.suptitle("ic_semiregular_law_20261005 round 2: held-out draws (4 per cell), macaulay_dense, seed 20261012", fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[1])
