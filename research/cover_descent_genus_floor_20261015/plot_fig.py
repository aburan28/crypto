"""Lay out route_changes.csv as a stacked-bar figure (layout only).

usage: python3 plot_fig.py route_changes.csv OUT.svg
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#c9c8c2"
COLOR = {"formula_advantage": BLUE, "formula_marginal": ORANGE, "boundary": GRID}
LABEL = {"formula_advantage": "formula advantage (≥ 8 bits below rho)",
         "formula_marginal": "formula marginal", "boundary": "boundary (at or above rho)"}
rows = list(csv.DictReader(open(sys.argv[1])))
ks = sorted({int(r["k"]) for r in rows})
# only routes that were a lead at genus k
leads = [r for r in rows if r["class_before"] in ("formula_advantage", "formula_marginal")]

fig, ax = plt.subplots(figsize=(7.4, 3.6), dpi=150)
ylabels, y = [], 0
for k in reversed(ks):
    sel = [r for r in leads if int(r["k"]) == k]
    if not sel:
        continue
    g0, g1 = sel[0]["genus_before"], sel[0]["genus_after"]
    for tag, key in ((f"priced at genus {g0}", "class_before"), (f"priced at GHS floor, genus {g1}", "class_after")):
        left = 0
        for cls in ("formula_advantage", "formula_marginal", "boundary"):
            n = sum(int(r["routes"]) for r in sel if r[key] == cls)
            if n:
                ax.barh(y, n, left=left, color=COLOR[cls], height=0.7, label=LABEL[cls])
                if n >= 15:
                    ax.text(left + n / 2, y, str(n), ha="center", va="center", fontsize=7.5,
                            color="white" if cls == "formula_advantage" else INK)
                left += n
        ylabels.append(f"k = {k}: {tag}")
        y += 1
    y += 0.6
pos, yy = [], 0
for k in reversed(ks):
    if not [r for r in leads if int(r["k"]) == k]:
        continue
    pos += [yy, yy + 1]
    yy += 2.6
ax.set_yticks(pos)
ax.set_yticklabels(ylabels, fontsize=8, color=INK)
ax.set_xlabel("cover_descent routes that were a lead at genus k (count)", color=MUTED, fontsize=8.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=MUTED, labelsize=8)
ax.grid(axis="x", color=GRID, linewidth=0.5, alpha=0.6)
ax.set_axisbelow(True)
h, l = ax.get_legend_handles_labels()
seen = {}
for hh, ll in zip(h, l):
    seen.setdefault(ll, hh)
ax.legend(seen.values(), seen.keys(), fontsize=7.5, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(sys.argv[2], format="svg")
fig.savefig(sys.argv[2].replace(".svg", ".png"), format="png", dpi=150)
