#!/usr/bin/env python3
"""Lay out summary.json (written by examples/generic_dlp_correctness_summary.rs)
as figures; every plotted number is a count already in summary.json.

    plot.py <study_dir>

Writes fig_correctness.svg and a 150 dpi PNG of it into <study_dir>.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Reference palette (dataviz skill, references/palette.md), light mode.
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BEFORE, AFTER = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "font.size": 9, "axes.titlesize": 10, "legend.frameon": False, "svg.fonttype": "none",
})

d = Path(sys.argv[1])
s = json.loads((d / "summary.json").read_text())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Left: dead trails per curve and mode, before and after.
rows = s["collab"]
labels = [f"{r['curve']}\n{'negation' if r['negation_map'] else 'plain'}" for r in rows]
x = range(len(rows))
w = 0.38
b1 = ax1.bar([i - w / 2 for i in x], [r["before"]["dead_trails"] for r in rows], w, color=BEFORE, label="before (main)")
b2 = ax1.bar([i + w / 2 for i in x], [r["after"]["dead_trails"] for r in rows], w, color=AFTER, label="after (fix)")
for bars in (b1, b2):
    for bar in bars:
        ax1.annotate(str(int(bar.get_height())), (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     ha="center", va="bottom", fontsize=7, color=INK2)
ax1.set_xticks(list(x))
ax1.set_xticklabels(labels, fontsize=7.5)
ax1.set_ylabel("walkers abandoned at the step cap (dead trails)")
ax1.set_title(f"pollard_collab walkers ({rows[0]['walkers']} or {rows[-1]['walkers']} per arm)")
ax1.legend(fontsize=8)

# Right: Pohlig–Hellman total_steps per solve, before against after.
p = s["ph_panel"]["rows"]
for bound, filled in ((16384, True), (131072, False)):
    pts = [r for r in p if r["bound"] == bound]
    ax2.scatter([r["steps_before"] for r in pts], [r["steps_after"] for r in pts], s=34 if filled else 70,
                marker="o", facecolors=AFTER if filled else "none", edgecolors=AFTER,
                label=f"bound 2^{bound.bit_length() - 1}" + ("" if filled else " (same counts)"))
lo = min(min(r["steps_after"] for r in p), min(r["steps_before"] for r in p)) * 0.8
hi = max(r["steps_before"] for r in p) * 1.2
ax2.plot([lo, hi], [lo, hi], color=INK2, lw=0.8, ls="--", label="no change")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("total_steps before (linear digits)")
ax2.set_ylabel("total_steps after (baby-step giant-step digits)")
ax2.set_title(f"Pohlig–Hellman, {s['ph_panel']['solves']} solves, same answers")
ax2.legend(fontsize=8)

fig.suptitle("generic-DLP correctness fixes: counts before and after (deterministic, nothing timed)", fontsize=9)
fig.tight_layout()
fig.savefig(d / "fig_correctness.svg")
fig.savefig(d / "fig_correctness.png", dpi=150)
print("figure written to", d)
