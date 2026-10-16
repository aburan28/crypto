"""Plot gf4_rows.csv and binary_curves.csv (recorded rows) as SVG.

Draws recorded ratio_log2 values only.  usage: python3 plot_gf4_rows.py DIR
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

d = Path(sys.argv[1])
rows = list(csv.DictReader(open(d / "gf4_rows.csv")))
curves = ["c2pnb176w1", "c2pnb272w1", "c2pnb304w1"]
y = {c: i for i, c in enumerate(reversed(curves))}

INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dd"
style = {
    ("before", "lower"): ("o", "#2a78d6", "before: m = 8 lower-bound row (heuristic)"),
    ("before", "upper"): ("s", "#eb6834", "before: m = 9 upper-bound row (heuristic)"),
    ("after", "gf4"): ("D", "#1baf7a", "after: m = 8 exact row (condition (8) holds)"),
    ("after", "other"): ("^", "#1baf7a", "after: best route on another tower"),
}
offset = {("before", "lower"): 0.18, ("before", "upper"): 0.18, ("after", "gf4"): -0.18, ("after", "other"): -0.18}

fig, ax = plt.subplots(figsize=(8.6, 3.9))
ax.axvline(0, color=MUTED, lw=1)
ax.axvline(8, color=MUTED, lw=1, ls=(0, (4, 3)))
ax.text(0, 2.62, "rho", ha="center", va="bottom", fontsize=8, color=MUTED)
ax.text(8, 2.62, "8-bit margin", ha="center", va="bottom", fontsize=8, color=MUTED)
seen = set()
for r in rows:
    gf4 = r["role"].startswith("gf4")
    key = (r["state"], r["bound"]) if r["state"] == "before" else ("after", "gf4" if gf4 else "other")
    marker, color, label = style[key]
    x = float(r["ratio_log2"])
    yy = y[r["curve"]] + offset[key]
    carried = r["class"] != "not_carried"
    ax.plot(
        x, yy, marker, ms=9, mew=2, color=color,
        mfc=color if carried else "white", mec=color,
        label=label if key not in seen else None, zorder=3,
    )
    seen.add(key)
    tag = r["class"].replace("formula_", "")
    if not gf4:
        tag += " (GF(2^76), m = 4)"
    if "best" in r["role"]:
        tag += " *"
    ax.annotate(tag, (x, yy), xytext=(7, 0), textcoords="offset points",
                va="center", fontsize=7.5, color=INK)
ax.set_yticks([y[c] for c in curves], curves, fontsize=9, color=INK)
ax.set_ylim(-0.6, 2.75)
ax.set_xlim(-12, 112)
ax.set_xlabel("ratio_log2 = log2(rho) - log2(route), bits (formula, not a measurement)", fontsize=8.5, color=INK)
ax.grid(axis="x", color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(MUTED)
ax.tick_params(axis="x", colors=MUTED, labelsize=8)
ax.tick_params(axis="y", length=0)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=7.5, frameon=False)
ax.set_title(
    "GHS descent rows over GF(2^2), before (b8ce3f697) and after the MMT Lemma 6 change\n"
    "hollow = not carried (Jacobian smaller than the subgroup); * = curve's best route",
    fontsize=9, color=INK, loc="left",
)
fig.tight_layout()
fig.savefig(d / "gf4_rows.svg")

# Second figure: every binary-field curve's best carried route (binary_curves.csv).
rows = list(csv.DictReader(open(d / "binary_curves.csv")))
changed = {"c2pnb176w1", "c2pnb304w1"}
cls = {
    "boundary": ("o", "#a3a29a", "boundary"),
    "formula_marginal": ("s", "#2a78d6", "formula_marginal"),
    "formula_advantage": ("D", "#eb6834", "formula_advantage"),
}
fig, ax = plt.subplots(figsize=(8.6, 4.4))
ax.axhline(0, color=MUTED, lw=1)
ax.axhline(8, color=MUTED, lw=1, ls=(0, (4, 3)))
ax.text(575, 0.8, "rho", fontsize=8, color=MUTED, ha="right")
ax.text(575, 8.8, "8-bit margin", fontsize=8, color=MUTED, ha="right")
# label placement: (dx, dy, ha) in points, chosen so labels clear markers
place = {"c2pnb304w1": (10, -4, "left"), "c2pnb368w1": (10, -4, "left"),
         "Oakley Group 4": (-8, 6, "right")}
seen = set()
for r in rows:
    if not r["best_ratio_log2"]:
        continue
    marker, color, label = cls[r["verdict"]]
    x, yv = float(r["subgroup_bits"]), float(r["best_ratio_log2"])
    big = r["verdict"] != "boundary"
    ax.plot(x, yv, marker, ms=8 if big else 5, color=color, mec="white", mew=1,
            label=label if label not in seen else None, zorder=3 if big else 2)
    seen.add(label)
    if big:
        tag = r["name"] + ("  (changed by Lemma 6)" if r["name"] in changed else "")
        dx, dy, ha = place.get(r["name"], (-8, 6, "right"))
        ax.annotate(tag, (x, yv), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va="center" if dy < 0 else "baseline", fontsize=7.5, color=INK)
ax.set_xlabel("subgroup order r, bits", fontsize=8.5, color=INK)
ax.set_ylabel("ratio_log2 of best carried route, bits", fontsize=8.5, color=INK)
ax.grid(color=GRID, lw=0.8)
ax.set_axisbelow(True)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
for s in ["left", "bottom"]:
    ax.spines[s].set_color(MUTED)
ax.tick_params(colors=MUTED, labelsize=8)
ax.legend(loc="upper left", fontsize=7.5, frameon=False)
ax.set_title(
    f"Every binary-field curve in docs/curves/jacobian-transfers.json ({len(rows)} curves), "
    "best carried route\nformula ratios under the frozen cost model, not measurements",
    fontsize=9, color=INK, loc="left",
)
fig.tight_layout()
fig.savefig(d / "binary_curves.svg")
