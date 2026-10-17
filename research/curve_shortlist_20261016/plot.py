"""Lay out recorded values from curves.json as volcano_depths.svg.

Plotting only: every number comes from curves.json, which shortlist.gp writes.
Run: python3 research/curve_shortlist_20261016/plot.py  (needs matplotlib)
"""
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
curves = {c["name"]: c for c in json.loads((HERE / "curves.json").read_text())["curves"]}

deep = [
    ("shortlist-deep-tall", "tall", "#3b6fb6"),
    ("shortlist-deep-broad", "broad", "#d9822b"),
    ("shortlist-deep-cycle", "cycle", "#2f9e6e"),
    ("shortlist-deep-binary", "binary", "#a14ac9"),
]
targets = [
    ("shortlist-prime-56", "prime-56 target"),
    ("shortlist-koblitz-73", "Koblitz-73 target (one prime ℓ)"),
    ("shortlist-binary-59", "binary-59 target"),
]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
x, ticks, labels = 0, [], []
for name, short, colour in deep:
    c = curves[name]
    walked = c["construction"]["descent_steps"]
    for ell, depth in c["v_factorization"].items():
        ax[0].bar(x, depth, color=colour, alpha=1.0 if ell in walked else 0.35,
                  hatch=None if ell in walked else "//")
        ticks.append(x)
        labels.append(f"{short}\nℓ={ell}")
        x += 1
    x += 0.8
ax[0].set_xticks(ticks)
ax[0].set_xticklabels(labels, fontsize=6.5, rotation=60)
ax[0].set_ylabel("volcano depth v_ℓ(v)  (levels below the crater)")
ax[0].set_title("Depth of each ℓ-volcano (hatched: not walked)", fontsize=10)

rows = [(short, math.log2(int(curves[n]["conductor_gap_v"])), col) for n, short, col in deep]
rows += [(label, math.log2(int(curves[n]["conductor_gap_v"])), "#888888") for n, label in targets]
ax[1].barh([r[0] for r in rows], [r[1] for r in rows], color=[r[2] for r in rows])
for i, r in enumerate(rows):
    ax[1].text(r[1] + 0.3, i, f"{r[1]:.1f}", va="center", fontsize=8, clip_on=False)
ax[1].set_xlabel("log2 of the conductor gap v = [O_K : Z[π]]")
ax[1].set_title("Conductor gap, crater to Z[π]", fontsize=10)
fig.suptitle("research/curve_shortlist_20261016/curves.json (7 curves, exact values)", fontsize=8, y=0.995)
plt.tight_layout()
plt.savefig(HERE / "volcano_depths.svg")

