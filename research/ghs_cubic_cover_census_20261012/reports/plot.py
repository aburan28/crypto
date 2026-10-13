"""Render recorded values from results/run-*.jsonl. Rendering only: every plotted number (curve counts n,
genus-3 fractions) is a `mod4` record written by cover_census.rs; this script derives nothing.
Usage: python3 plot.py <results.jsonl> <outdir>"""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

src, out = sys.argv[1], sys.argv[2]
rows = [json.loads(l) for l in open(src)]
mod4 = [r for r in rows if r["rec"] == "mod4"]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#6b6a63", "#e4e3dc"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
GROUPS = [("odd", "odd"), ("2mod4", "2 mod 4"), ("0mod4", "0 mod 4")]
qs = sorted({r["q"] for r in mod4})
fig, axes = plt.subplots(1, len(qs), figsize=(7.2, 3.3), sharey=True)
for ax, q in zip(axes, qs):
    for k, (key, label) in enumerate(GROUPS):
        for j, (scope, col) in enumerate((("self", BLUE), ("2isog", ORANGE))):
            rec = next(r for r in mod4 if r["q"] == q and r["scope"] == scope and r["order_mod4"] == key)
            frac = rec["genus3_fraction"]
            x = k + (j - 0.5) * 0.34
            ax.bar(x, frac, width=0.30, color=col, edgecolor="white", linewidth=1)
            ax.text(x, frac + 0.02, f"{frac:.2f}" if frac > 0 else "0", ha="center", va="bottom", color=INK, fontsize=7.5)
            if j == 0:
                ax.text(k, -0.19, f"n={rec['curves']}", ha="center", va="top", color=MUTED, fontsize=7.5, transform=ax.get_xaxis_transform())
    mode = next(r["mode"] for r in rows if r["rec"] == "field" and r["q"] == q)
    ax.set_title(f"q={q} ({'all classes' if mode == 'exhaustive' else 'sample'})", color=INK, fontsize=9.5)
    ax.set_xticks(range(3)); ax.set_xticklabels([g[1] for g in GROUPS]); ax.set_ylim(0, 1.12)
    ax.grid(axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
axes[0].set_ylabel("fraction of curves with a genus-3 cover")
fig.text(0.5, -0.02, "#E(K) modulo 4", ha="center", color=INK)
fig.legend(handles=[Patch(color=BLUE, label="on the curve itself"), Patch(color=ORANGE, label="on the curve or a 2-isogenous neighbour")],
           frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=2)
fig.savefig(f"{out}/genus3_fraction.svg", bbox_inches="tight"); fig.savefig(f"{out}/genus3_fraction.png", dpi=150, bbox_inches="tight")
