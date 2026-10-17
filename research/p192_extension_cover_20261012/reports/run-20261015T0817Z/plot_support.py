"""Render only the already recorded native support counts."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
receipt = json.loads((HERE.parent.parent / "type_i_conditioned_support_native.json").read_text())
cells = ["++", "+-", "-+", "--"]
counts = [receipt["counts_by_character_pair"][cell] for cell in cells]
assert counts == [224, 211, 203, 178]
fig, ax = plt.subplots(figsize=(7.5, 4.1))
bars = ax.bar(cells, counts, color=["#167d8d", "#a7b8c1", "#a7b8c1", "#a7b8c1"])
ax.bar_label(bars, padding=5, fontsize=11)
ax.set_ylim(0, 255)
ax.set_ylabel("Distinct fixed-base three-sum targets")
ax.set_xlabel("Character of P and 2P; q = -2P")
ax.set_title("P-192 Type I: exact partition of 816 support targets")
ax.text(0.98, 0.98, "Accepted cell: 224\nNative finite census; n=816",
        transform=ax.transAxes, va="top", ha="right", fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "#cccccc"})
ax.grid(axis="y", alpha=0.2)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(HERE / "support.svg")
fig.savefig(HERE / "support.png", dpi=150)
plt.close(fig)
svg = HERE / "support.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
