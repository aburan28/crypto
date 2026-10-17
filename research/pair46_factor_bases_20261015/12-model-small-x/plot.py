#!/usr/bin/env python3
"""Plot recorded summary values from results.json (no computation beyond layout)."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("results.json"))
S = d["summary"]
pops = [("M", ["W-x", "M-u", "E-y", "E-x", "R-sym2"]), ("H", ["W-x", "H-x", "H-s", "R-sym2"])]
sizes = sorted({s["bits"] for s in S})
colors = {16: "#4C78A8", 20: "#F58518", 22: "#54A24B"}
fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), gridspec_kw={"width_ratios": [5, 4, 3]})
for ax, (pop, arms) in zip(axes[:2], pops):
    for si, b in enumerate(sizes):
        for ai, arm in enumerate(arms):
            r = [s for s in S if s["pop"] == pop and s["bits"] == b and s["arm"] == arm]
            if not r:
                continue
            r = r[0]
            y = ai + (si - 1) * 0.22
            m = 100 * r["mean_rel"]
            lo, hi = 100 * r["ci_bonf"][0], 100 * r["ci_bonf"][1]
            ax.errorbar(m, y, xerr=[[m - lo], [hi - m]], fmt="o", color=colors[b], ms=4, capsize=2,
                        label=f"b={b} (n={r['n_curves']})" if ai == 0 else None)
    ax.axvspan(-1, 1, color="0.9", zorder=0)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(range(len(arms)))
    ax.set_yticklabels([f"{a} vs {'R-asym' if a in ('E-x','H-x') else 'R-sym'}" for a in arms])
    ax.invert_yaxis()
    ax.set_xlim(-2, 2)
    ax.set_xlabel("relative coverage difference C_arm/C_ctrl - 1  (%)")
    ax.set_title(f"Pop {pop}: {'Montgomery/Edwards' if pop=='M' else 'Hessian'} curves")
    ax.legend(fontsize=7, loc="lower left")
ax = axes[2]
for si, b in enumerate(sizes):
    for ai, (pop, arm) in enumerate([("M", "PC1"), ("M", "PC2"), ("H", "PC1"), ("H", "PC2")]):
        r = [s for s in S if s["pop"] == pop and s["bits"] == b and s["arm"] == arm]
        if not r:
            continue
        r = r[0]
        y = ai + (si - 1) * 0.22
        m = 100 * r["mean_rel"]
        lo, hi = 100 * r["ci_bonf"][0], 100 * r["ci_bonf"][1]
        ax.errorbar(m, y, xerr=[[m - lo], [hi - m]], fmt="s", color=colors[b], ms=4, capsize=2)
ax.axvspan(-1, 1, color="0.9", zorder=0)
ax.set_yticks(range(4)); ax.set_yticklabels(["M PC1", "M PC2", "H PC1", "H PC2"]); ax.invert_yaxis()
ax.set_xlim(-105, 5); ax.set_xlabel("relative coverage difference (%)")
ax.set_title("Positive controls (vs R-sym)")
fig.suptitle("12-model-small-x: m=2 coverage |F+F|/N of small-coordinate factor bases, mean over curves, "
             "Bonferroni (21 tests) bootstrap CI; grey = +-1% equivalence margin", fontsize=9)
fig.tight_layout()
fig.savefig("fig_coverage.svg"); fig.savefig("fig_coverage.png", dpi=150)
