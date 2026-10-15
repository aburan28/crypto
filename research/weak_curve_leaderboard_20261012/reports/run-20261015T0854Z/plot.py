"""Render recorded GHS model prices; no research arithmetic is done here."""

import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
receipt = json.loads((HERE / "ghs_genus128_cost_audit.json").read_text())
fig, ax = plt.subplots(figsize=(8.2, 4.6))
for curve, color in zip(receipt["curves"], ["#176a9c", "#b45625"]):
    rows = [row for row in curve["audit"]["sweep"] if row["degree_bound"] <= 30]
    ax.plot(
        [row["degree_bound"] for row in rows],
        [row["log2_total"] for row in rows],
        marker="o", markersize=2.8, linewidth=1.7,
        label=curve["name"], color=color,
    )
ax.axvline(15, color="#56616b", linestyle="--", linewidth=1)
ax.annotate("model minimum: b = 15", (15, 74), (18, 105),
            arrowprops={"arrowstyle": "->", "color": "#56616b"}, fontsize=9)
ax.set(xlabel="maximum irreducible degree b over GF(4)",
       ylabel="modeled total cost (log2 bit operations)",
       title="Conditional genus-128 GHS cost sweep (84 bounds tested)")
ax.set_xlim(1, 30)
ax.set_ylim(60, 265)
ax.grid(alpha=0.2)
ax.legend()
fig.text(0.5, 0.01,
         "run-20261015T0854Z · two named source curves · computed model, no measured relation yield",
         ha="center", fontsize=8)
fig.tight_layout(rect=(0, 0.035, 1, 1))
fig.savefig(HERE / "ghs_cost_sweep.svg")
plt.close(fig)
