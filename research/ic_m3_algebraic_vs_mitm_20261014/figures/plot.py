"""Layout only: draws the values `ic_algebraic_floor analyze` recorded
(cells.csv, ANALYSIS.json). No value is computed here beyond axis
placement. Usage: python plot.py <study dir> [<reference dir>]: arms of the
optional reference dir (other than its floor) are overlaid, faded, labelled
"original run"."""
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

study = Path(sys.argv[1])
rows = list(csv.DictReader(open(study / "cells.csv")))
if len(sys.argv) > 2:
    for r in csv.DictReader(open(Path(sys.argv[2]) / "cells.csv")):
        if r["arm"] != "floor-mitm":
            r["arm"] = "orig:" + r["arm"]
            rows.append(r)
ana = json.load(open(study / "ANALYSIS.json"))
out = study / "figures"
out.mkdir(exist_ok=True)

STYLE = {
    "floor-mitm": ("#1b1b1b", "o", "MITM floor (pair table, amortized)"),
    "f4-chain": ("#2a6fdb", "s", "inherited F4, S3 chain"),
    "f6-ic": ("#13a37f", "D", "F6-IC (F4 + support/closure)"),
    "sat": ("#d1495b", "^", "CDCL+XOR SAT, chain"),
    "s4-hybrid-f10": ("#8f5bd6", "v", "unsplit S4 hybrid, 10 free"),
    "s4-hybrid-f12": ("#c27c0e", "P", "unsplit S4 hybrid, 12 free"),
    "s22-hybrid-g8": ("#8f5bd6", "v", "S3+S4 split hybrid, g=8"),
    "wide-chain": ("#e36414", "X", "wide u128 chain solver (amendment 1)"),
}


def save(fig, name):
    fig.savefig(out / f"{name}.svg")
    fig.savefig(out / f"{name}.png", dpi=150)
    plt.close(fig)


for m in sorted({r["m"] for r in rows if r["arm"] == "floor-mitm"}):
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    arms = []
    for r in rows:
        if r["m"] == m and r["arm"] not in arms:
            arms.append(r["arm"])
    for arm in arms:
        orig = arm.startswith("orig:")
        color, marker, label = STYLE.get(arm.removeprefix("orig:"), ("grey", "x", arm))
        if orig:
            label = f"{label} (original run)"
        alpha = 0.35 if orig else 1.0
        pts = [r for r in rows if r["m"] == m and r["arm"] == arm]
        meas = [r for r in pts if r["status"] == "measured"]
        lb = [r for r in pts if r["status"] == "lower_bound"]
        if meas:
            xs = [int(r["n"]) for r in meas]
            ys = [float(r["log2_cost_ns"]) for r in meas]
            lo = [y - float(r["lo"]) for y, r in zip(ys, meas)]
            hi = [float(r["hi"]) - y for y, r in zip(ys, meas)]
            ax.errorbar(xs, ys, yerr=[lo, hi], color=color, marker=marker,
                        capsize=3, lw=1.4, ms=6, label=label, alpha=alpha)
        if lb:
            ax.scatter([int(r["n"]) for r in lb], [float(r["log2_cost_ns"]) for r in lb],
                       facecolors="none", edgecolors=color, marker=marker, s=60, alpha=alpha,
                       label=f"{label}: lower bound (censored)" if not meas else None)
            for r in lb:
                ax.annotate("", xy=(int(r["n"]), float(r["log2_cost_ns"]) + 1.6),
                            xytext=(int(r["n"]), float(r["log2_cost_ns"])),
                            arrowprops=dict(arrowstyle="->", color=color, lw=1))
    lo_y, hi_y = ax.get_ylim()
    ax.set_ylim(lo_y, hi_y + 2.5)  # room for the lower-bound arrows
    ax.set_xlabel("field degree n (cells a:n:d with d = ceil(n/m))")
    ax.set_ylabel("log2 wall ns per relation (one core)")
    ax.set_title(f"m = {m}: cost per relation, algebraic/SAT arms vs MITM floor\n"
                 "bars: 95% bootstrap; open markers with arrows: lower bounds")
    ax.grid(alpha=0.3)
    if m == "4" and not any(
        r["m"] == "4" and r["arm"] == "wide-chain" for r in rows
    ):
        ax.axvspan(21.5, 32, color="#999999", alpha=0.12)
        ax.text(26.5, 31, "F4 / F6-IC / SAT unsupported here:\nchain has 4d+2n > 64 unknowns\n(MAX_VARS); floor only", ha="center", fontsize=7)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    save(fig, f"cost_per_relation_m{m}")

# Slopes with intervals (from ANALYSIS.json).
fig, ax = plt.subplots(figsize=(7.2, 3.8))
labels, vals, los, his, cols = [], [], [], [], []
for m, blk in sorted(ana["by_m"].items()):
    labels.append(f"m={m} floor")
    vals.append(blk["floor_slope_amortized"])
    lo, hi = blk["floor_slope_amortized_ci"]
    los.append(vals[-1] - lo)
    his.append(hi - vals[-1])
    cols.append(STYLE["floor-mitm"][0])
    for arm, a in blk["arms"].items():
        if a["slope_measured"] is not None and a["slope_measured_ci"] is not None and len(a["measured_cells"]) >= 3:
            labels.append(f"m={m} {arm}")
            vals.append(a["slope_measured"])
            lo, hi = a["slope_measured_ci"]
            los.append(vals[-1] - lo)
            his.append(hi - vals[-1])
            cols.append(STYLE.get(arm, ("grey",))[0])
ys = list(range(len(labels)))
for v, y, lo, hi, c in zip(vals, ys, los, his, cols):
    ax.errorbar([v], [y], xerr=[[lo], [hi]], fmt="o", color=c, capsize=3)
ax.set_yticks(ys)
ax.set_yticklabels(labels, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("slope: log2 ns per relation per unit n (95% bootstrap)")
ax.set_title("Slope on each arm's measured cells (>= 3); floor rows: slope over all cells.\nThe decision uses arm − floor fitted on the same cells (ANALYSIS.md)", fontsize=9)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
save(fig, "slopes")
print("figures written to", out)
