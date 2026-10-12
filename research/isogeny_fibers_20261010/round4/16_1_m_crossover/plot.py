"""Plot recorded values only (results/panel_*.json written by
examples/r4_16_1_m_crossover.rs).  No computation beyond reading fields.

Writes crossover.svg (and crossover_data.csv, the plotted rows)."""
import csv
import glob
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for path in sorted(glob.glob(os.path.join(HERE, "results", "panel_*.json"))):
    doc = json.load(open(path))
    for c in doc["curves"]:
        for cell in c["cells"]:
            if cell["arm"] == "_base" or cell.get("skipped"):
                continue
            rows.append({
                "n": c["n"], "m": cell["m"], "ell": cell["ell"],
                "inv": cell.get("v_invariant", False), "arm": cell["arm"],
                "attempts": cell["attempts"], "found": cell["found"],
                "refuted": cell["refuted"], "budget": cell["budget_exhausted"],
                "unbuildable": cell["unbuildable"],
                "solved_fraction": cell["solved_fraction"],
                "xor_per_relation": cell["xor_per_relation"],
                "feasible": cell["feasible"], "wall_truncated": cell["wall_truncated"],
            })

with open(os.path.join(HERE, "crossover_data.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

ns = sorted({r["n"] for r in rows})
COL = {"x-direct": "#4269d0", "x-chained": "#97bbf5", "sym": "#ff725c", "sym-chained": "#ffb3a8"}
fig, axes = plt.subplots(2, len(ns), figsize=(4.2 * len(ns), 6.4), sharey="row")
for j, n in enumerate(ns):
    sub = [r for r in rows if r["n"] == n]
    groups = []
    for r in sub:
        key = (r["m"], r["ell"], r["inv"])
        if key not in groups:
            groups.append(key)
    groups.sort()
    for k, key in enumerate(groups):
        cell = [r for r in sub if (r["m"], r["ell"], r["inv"]) == key]
        for i, r in enumerate(cell):
            x = k + (i - (len(cell) - 1) / 2) * 0.22
            y = r["xor_per_relation"]
            ok = isinstance(y, (int, float)) and math.isfinite(y)
            ax = axes[0][j]
            if ok:
                ax.bar(x, math.log10(y), 0.2, color=COL[r["arm"]],
                       edgecolor="black" if r["feasible"] else "none",
                       hatch=None if r["feasible"] else "//")
            else:
                ax.text(x, 0.3, "no rel." if r["unbuildable"] < r["attempts"] else ">64 vars",
                        rotation=90, ha="center", va="bottom", fontsize=7)
            axes[1][j].bar(x, r["solved_fraction"], 0.2, color=COL[r["arm"]])
    labels = [f"m={m}\nℓ={l}{' inv' if inv else ''}" for (m, l, inv) in groups]
    for row in (0, 1):
        axes[row][j].set_xticks(range(len(groups)))
        axes[row][j].set_xticklabels(labels, fontsize=8)
    axes[0][j].set_title(f"n = {n}")
axes[0][0].set_ylabel("log10 word XORs per relation\n(hatched: not feasible)")
axes[1][0].set_ylabel("solved fraction\n(found + refuted) / attempts")
axes[1][0].axhline(0.5, color="gray", lw=0.8, ls=":")
for j in range(1, len(ns)):
    axes[1][j].axhline(0.5, color="gray", lw=0.8, ls=":")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COL.values()]
fig.legend(handles, list(COL.keys()), loc="lower center", ncol=4, frameon=False)
fig.suptitle("Lead 16.1: plain x-frame vs kernel-saturated u-frame, InheritedF4, node budget 4096, 24 targets per n",
             fontsize=10)
fig.tight_layout(rect=(0, 0.05, 1, 0.95))
fig.savefig(os.path.join(HERE, "crossover.svg"))
fig.savefig(os.path.join(HERE, "crossover.png"), dpi=130)
print("wrote crossover.svg", len(rows), "rows")
