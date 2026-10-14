"""Plot S from quartic_ic trial JSONL.
Usage: python3 plot.py <results.jsonl> <outdir>"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

src, out = sys.argv[1], sys.argv[2]
trials = [json.loads(line) for line in open(src) if '"rec":"trial"' in line]
by_key = {}
for t in trials:
    key = (t["q"], t["variant"])
    by_key.setdefault(key, []).append(t["S"])
keys = sorted(by_key.keys())
labels = [f"q={q}\n{v}" for q, v in keys]
means = [sum(by_key[k]) / len(by_key[k]) for k in keys]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#6b6a63", "#e4e3dc"
plt.rcParams.update(
    {
        "font.size": 9,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
    }
)
fig, ax = plt.subplots(figsize=(7.2, 3.4))
xs = range(len(keys))
ax.bar(xs, means, color=BLUE, label="mean S (ops/√ℓ)", zorder=2)
for i, k in enumerate(keys):
    ys = by_key[k]
    ax.scatter([i] * len(ys), ys, color=ORANGE, s=28, zorder=3)
ax.set_ylabel("S")
ax.set_xticks(list(xs))
ax.set_xticklabels(labels, fontsize=8)
ax.set_yscale("log")
ax.grid(axis="y", color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left")
ax.set_title("Toy smokes — not a security claim")
fig.savefig(f"{out}/S_mean.svg", bbox_inches="tight")
fig.savefig(f"{out}/S_mean.png", dpi=150, bbox_inches="tight")
