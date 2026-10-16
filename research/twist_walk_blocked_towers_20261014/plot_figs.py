"""Lay out the recorded rows of summary.csv and class_ceiling.csv (layout only).

usage: python3 plot_figs.py summary.csv class_ceiling.csv OUT_DIR
Writes escape_by_config.svg and class_ceiling.svg.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#c9c8c2"
summary, ceiling, out = sys.argv[1], sys.argv[2], sys.argv[3]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="x", color=GRID, linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)


rows = list(csv.DictReader(open(summary)))
null = [r for r in rows if r["kind"] == "null"]
cat = [r for r in rows if r["kind"] == "catalog"]
frac = lambda s: int(s.split("/")[0]) / int(s.split("/")[1])  # noqa: E731

fig, ax = plt.subplots(figsize=(7.4, 4.6), dpi=150)
labels = []
for i, r in enumerate(null):
    y = len(null) - 1 - i
    key = (r["N"], r["l"])
    labels.append(f"N={r['N']}, l={r['l']} (n={r['n']})")
    p, lo, hi = frac(r["escape_le4"]), float(r["escape_le4_ci_lo"]), float(r["escape_le4_ci_hi"])
    ax.errorbar([p], [y], xerr=[[p - lo], [hi - p]], fmt="o", color=BLUE, ms=5, lw=1.2, capsize=2.5,
                label="random blocked starts: escaped within 4 steps (95% CI)" if i == 0 else None)
    for c in [c for c in cat if (c["N"], c["l"]) == key]:
        escaped = c["first_unblocked_depth"] not in ("", "none")
        has_nb = c["has_neighbour"] == "true"
        if escaped:
            ax.plot([1.0], [y + 0.22], "D", color=ORANGE, ms=6, label="catalog tower: escaped")
        elif has_nb:
            ax.plot([0.0], [y + 0.22], "D", mec=ORANGE, mfc="#fcfcfb", ms=6, mew=1.4, label="catalog tower: has 3-isogenies, none unblocked")
        else:
            ax.plot([0.0], [y + 0.22], "x", color=ORANGE, ms=7, mew=1.6, label="catalog tower: no rational 3-isogeny")
ax.set_yticks(range(len(null)))
ax.set_yticklabels(labels[::-1], fontsize=8, color=INK)
ax.set_xlim(-0.04, 1.04)
ax.set_xlabel("share of starts that reach an unblocked twist within 4 steps of the 3-isogeny graph", color=MUTED, fontsize=8.5)
style(ax)
h, l = ax.get_legend_handles_labels()
seen = {}
for hh, ll in zip(h, l):
    seen.setdefault(ll, hh)
ax.legend(seen.values(), seen.keys(), loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=7.5)
fig.tight_layout()
fig.savefig(os.path.join(out, "escape_by_config.svg"), format="svg", facecolor="#fcfcfb")

cr = list(csv.DictReader(open(ceiling)))
fig, ax = plt.subplots(figsize=(7.4, 3.2), dpi=150)
for r in cr:
    bl, bm = int(r["curves_blocked"]), int(r["curves_blocked_in_mixed_classes"])
    ax.plot([int(r["l"])], [bm / bl], "o", color=BLUE, ms=5, alpha=0.85)
ax.set_xlabel("base degree l (one point per field degree N = 8 ... 22 with n = N/l even)", color=MUTED, fontsize=8.5)
ax.set_ylabel("blocked curves whose class also\ncontains an unblocked curve", color=MUTED, fontsize=8.5)
ax.set_ylim(-0.04, 1.04)
ax.set_xticks(sorted({int(r["l"]) for r in cr}))
ax.annotate("l = 1: no mixed class in any of\nthe 8 field sizes, so no isogeny\ncan escape", xy=(1, 0.0), xytext=(1.6, 0.2), fontsize=8, color=INK,
            arrowprops=dict(arrowstyle="-", color=GRID))
style(ax)
ax.grid(axis="y", color=GRID, linewidth=0.5, alpha=0.6)
fig.tight_layout()
fig.savefig(os.path.join(out, "class_ceiling.svg"), format="svg", facecolor="#fcfcfb")
print("wrote figures")
