"""Campaign summary graph for the pair-46 factor-base campaign.

Plots only values recorded in each item's RESULTS.md / results.json (as
quoted in summary_data.json). No computation beyond layout.
"""

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
data = json.loads((HERE / "summary_data.json").read_text())

COLORS = {
    "negative": "#c0392b",
    "exact_obstruction": "#2c3e50",
    "null": "#7f8c8d",
    "constant_only": "#d68910",
    "empirical_negative": "#8e44ad",
}

rows = data["ratios"]
fig, (ax, az) = plt.subplots(
    1, 2, figsize=(10.5, 5.6), sharey=True, gridspec_kw={"width_ratios": [1.25, 1]}
)
for i, r in enumerate(rows):
    y = len(rows) - 1 - i
    c = COLORS[r["verdict"]]
    lo, hi = r["lo"], r["hi"]
    pt = r["point"]
    for a in (ax, az):
        a.plot([lo, hi], [y, y], color=c, lw=2.4, solid_capstyle="butt")
        if pt is not None:
            a.plot(pt, y, "o", color=c, ms=5.5)
        else:
            a.plot([lo, hi], [y, y], "|", color=c, ms=9, mew=1.6)
    if not (0.94 <= lo and hi <= 1.06):
        az.annotate(
            "off scale" + (" (>)" if lo > 1 else " (<)"),
            xy=(1.055 if lo > 1 else 0.945, y),
            ha="right" if lo > 1 else "left",
            va="center",
            fontsize=7,
            color=c,
        )
for a in (ax, az):
    a.axvline(1.0, color="#555", lw=0.8, ls="--")
    a.grid(axis="x", color="#ddd", lw=0.5)
ax.set_xscale("log")
ax.set_xlim(0.005, 3.0)
az.set_xlim(0.94, 1.06)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r["label"] for r in reversed(rows)], fontsize=8)
ax.set_xlabel("ratio to matched reference (log scale)")
az.set_xlabel("same, zoom 0.94 to 1.06 (linear)")
fig.suptitle(
    "pair-46 factor-base campaign (2026-10-15): headline ratios with 95% intervals "
    "(1 = no effect)",
    fontsize=9.5,
)
handles = [
    plt.Line2D([], [], color=c, lw=3, label=k.replace("_", " "))
    for k, c in COLORS.items()
    if any(r["verdict"] == k for r in rows)
]
ax.legend(handles=handles, loc="lower left", fontsize=7.5, frameon=True)
fig.text(
    0.01,
    0.01,
    "Recorded values from each item's RESULTS.md (bootstrap/test 95% CIs as reported); "
    "item 01 is a cost ratio (>1 = worse), the rest are yield/coverage/density ratios.\n"
    "Item 03 is drawn as its CI-end range over 8 toys (no pooled point). Items 02, 06, 07, "
    "10, 11, 12, 14 have no single ratio of this kind; see README table.",
    fontsize=6.8,
    color="#444",
)
fig.subplots_adjust(left=0.21, right=0.98, bottom=0.17, top=0.92, wspace=0.06)
fig.savefig(HERE / "fig_campaign_summary.svg")
fig.savefig(HERE / "fig_campaign_summary.png", dpi=150)
