"""Plot recorded values from results/run-*.jsonl. Layout only: it groups the per-curve records the tool wrote.
Usage: python3 plot.py <results.jsonl> <outdir>"""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src, out = sys.argv[1], sys.argv[2]
rows = [json.loads(l) for l in open(src)]
curves = [r for r in rows if r["rec"] == "curve"]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#6b6a63", "#e4e3dc"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})

def v2(n):
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return min(k, 2)

qs = sorted({r["q"] for r in curves})
labels = ["odd", "2 mod 4", "0 mod 4"]
fig, axes = plt.subplots(1, len(qs), figsize=(7.2, 3.3), sharey=True)
for ax, q in zip(axes, qs):
    cs = [r for r in curves if r["q"] == q]
    for k in range(3):
        grp = [r for r in cs if v2(r["order"]) == k]
        if not grp: continue
        for j, (key, col) in enumerate((("min_genus", BLUE), ("min_genus_2isog", ORANGE))):
            frac = sum(1 for r in grp if r[key] == 3) / len(grp)
            x = k + (j - 0.5) * 0.34
            ax.bar(x, frac, width=0.30, color=col, edgecolor="white", linewidth=1)
            if frac > 0: ax.text(x, frac + 0.02, f"{frac:.2f}", ha="center", va="bottom", color=INK, fontsize=7.5)
        if all(r["min_genus_2isog"] != 3 for r in grp): ax.text(k, 0.02, "0", ha="center", va="bottom", color=INK, fontsize=7.5)
        ax.text(k, -0.19, f"n={len(grp)}", ha="center", va="top", color=MUTED, fontsize=7.5, transform=ax.get_xaxis_transform())
    mode = next(r["mode"] for r in rows if r["rec"] == "field" and r["q"] == q)
    ax.set_title(f"q={q} ({'all classes' if mode == 'exhaustive' else 'sample'})", color=INK, fontsize=9.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(labels); ax.set_ylim(0, 1.12)
    ax.grid(axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
axes[0].set_ylabel("fraction of curves with a genus-3 cover")
fig.text(0.5, -0.02, "#E(K) modulo 4", ha="center", color=INK)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=BLUE, label="on the curve itself"), Patch(color=ORANGE, label="on the curve or a 2-isogenous neighbour")],
           frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=2)
fig.savefig(f"{out}/genus3_fraction.svg", bbox_inches="tight"); fig.savefig(f"{out}/genus3_fraction.png", dpi=150, bbox_inches="tight")
