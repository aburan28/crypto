"""Lay out summary.csv (layout only).  usage: python3 plot_figs.py summary.csv OUT_DIR"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#c9c8c2"
rows = list(csv.DictReader(open(sys.argv[1])))
out = sys.argv[2]
g_rows = [r for r in rows if r["kind"] == "genus"]
sens = [r for r in rows if r["kind"] == "sensitivity"][0]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)


fig, ax = plt.subplots(figsize=(7.4, 3.6), dpi=150)
ax.axhspan(-1, 1, color=GRID, alpha=0.35, lw=0, label="within one bit of the model (decision rule)")
for op, colour in (("add", BLUE), ("double", ORANGE)):
    for p, filled in (("1009", True), ("65521", False)):
        sel = [r for r in g_rows if r["op"] == op and r["p"] == p]
        xs = [int(r["g"]) for r in sel]
        ys = [float(r["log2_ratio_to_model"]) for r in sel]
        ax.plot(xs, ys, "-o" if filled else "o", color=colour, ms=4.5, lw=1.2 if filled else 0,
                mfc=colour if filled else "#fcfcfb", mew=1.2, label=f"{op}, p = {p}")
ax.axhline(0, color=INK, lw=0.8)
ax.set_xscale("log", base=2)
ax.set_xticks([2, 4, 6, 8, 10, 14, 16, 24, 32])
ax.set_xticklabels([2, 4, 6, 8, 10, 14, 16, 24, 32])
ax.set_xlabel("genus g", color=MUTED, fontsize=8.5)
ax.set_ylabel("log2 (counted cost / 10 g²)", color=MUTED, fontsize=8.5)
ax.set_ylim(-1.5, 4)
for g in (6, 10, 14):
    ax.axvline(g, color=GRID, lw=0.8, ls=":")
style(ax)
ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(out, "ratio_by_genus.svg"), format="svg", facecolor="#fcfcfb")

field = next(v for v in sens.values() if v and "advantage_curves=" in v)
detail = dict(kv.split("=") for kv in field.split(";"))
vals = [
    ("slack ≤ 1 bit", int(detail["slack_le_1bit"])),
    ("slack ≤ 2 bits", int(detail["slack_le_2bits"])),
    ("slack ≤ 3 bits", int(detail["slack_le_3bits"])),
    ("would drop below the margin\n(measured cost, whole route)", int(detail["would_drop_below_margin_if_measured_cost_applied_to_the_whole_route"])),
]
fig, ax = plt.subplots(figsize=(7.4, 2.6), dpi=150)
ys = list(range(len(vals)))[::-1]
ax.barh(ys, [v for _, v in vals], color=[BLUE, BLUE, BLUE, ORANGE], height=0.55)
for y, (_, v) in zip(ys, vals):
    ax.text(v + 0.6, y, str(v), va="center", fontsize=8.5, color=INK)
ax.set_yticks(ys)
ax.set_yticklabels([l for l, _ in vals], fontsize=8, color=INK)
ax.set_xlim(0, 45)
ax.set_xlabel(f"of {detail['advantage_curves']} formula_advantage verdicts ({detail['marginal_curves']} more are formula_marginal)", color=MUTED, fontsize=8.5)
style(ax)
ax.grid(axis="y", visible=False)
fig.tight_layout()
fig.savefig(os.path.join(out, "verdict_fragility.svg"), format="svg", facecolor="#fcfcfb")
print("wrote figures")
