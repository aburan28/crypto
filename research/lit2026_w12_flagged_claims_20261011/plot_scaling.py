"""Plot recorded values of a lit2026 W12 run (layout only; computes nothing).

Usage: python -I plot_scaling.py <run-dir>
Reads <run-dir>/scaling.json, <run-dir>/fits.json and, if present,
<run-dir>/dsweep.json; writes scaling.svg and dsweep.svg there.  Every number
drawn comes from the JSON the Rust example wrote (runs, fitted slopes and
intercepts with their intervals, paired ratios with their intervals); the only
arithmetic here is log2 of recorded counts for the axes.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

run = Path(sys.argv[1])
doc = json.loads((run / "scaling.json").read_text())
fits = json.loads((run / "fits.json").read_text())["fits_log2_vs_log2r"]
runs = [r for r in doc["runs"] if r["verified"]]
rho = [r for r in runs if r["method"] == "rho.negation"]
zm = [r for r in runs if r["method"].startswith("zero_minor")]


def line(ax, f, xs, color, ls, label):
    x0, x1 = min(xs), max(xs)
    ax.plot([x0, x1], [f["intercept"] + f["slope"] * x0, f["intercept"] + f["slope"] * x1],
            color=color, ls=ls, lw=1.6, label=label)


def lab(name, f):
    lo, hi = f["slope_ci95_bootstrap"]
    return f"{name}: slope {f['slope']:.3f} [{lo:.3f}, {hi:.3f}], n={f['n']}"


fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
ax = axes[0]
ax.scatter([r["log2_r"] for r in rho], [math.log2(r["ops"]) for r in rho], s=10, alpha=0.5, color="#1f77b4")
line(ax, fits["rho.negation [field ops]"], [r["log2_r"] for r in rho], "#1f77b4", "-",
     lab("rho.negation", fits["rho.negation [field ops]"]))
for d, mk in [(3, "o"), (4, "s")]:
    sel = [r for r in zm if r["detail"]["d"] == d]
    if sel:
        ax.scatter([r["log2_r"] for r in sel], [math.log2(r["ops"]) for r in sel], s=14, alpha=0.6,
                   color="#d62728", marker=mk, label=f"zero-minor runs, d = {d}")
k = "zero_minor, d by Eq.(4) rule [field ops, all phases]"
line(ax, fits[k], [r["log2_r"] for r in zm], "#d62728", "-", lab("zero-minor, all phases", fits[k]))
k = "zero_minor, d by Eq.(4) rule [field ops, search phases only]"
line(ax, fits[k], [r["log2_r"] for r in zm], "#ff7f0e", "--", lab("zero-minor, search phases only", fits[k]))
ax.set_xlabel("log2 r (prime group order)")
ax.set_ylabel("log2 field ops (F_p mul + sqr + inv)")
ax.set_title("(a) Cost per solved target", fontsize=10)
ax.legend(fontsize=6.5, loc="upper left")
ax.grid(alpha=0.3)

ax = axes[1]
sel = [r for r in zm if r["detail"]["counters"].get("covered_d_subsets", 0) > 0]
xs = [r["log2_r"] for r in sel]
ax.scatter(xs, [math.log2(r["detail"]["counters"]["covered_d_subsets"]) for r in sel], s=14, alpha=0.6,
           color="#2ca02c", label="d-subsets covered until the zero minor")
k = "zero_minor, d by Eq.(4) rule [d-subsets covered]"
line(ax, fits[k], xs, "#2ca02c", "-", lab("fit", fits[k]))
ax.plot([min(xs), max(xs)], [min(xs), max(xs)], color="gray", ls=":", label="reference line: covered = r")
ax.set_xlabel("log2 r")
ax.set_ylabel("log2 candidate minors covered")
ax.set_title("(b) Determine step: minors examined", fontsize=10)
ax.legend(fontsize=6.5, loc="upper left")
ax.grid(alpha=0.3)

ax = axes[2]
pr = doc["paired_ratio"]
xs = [p["bits"] for p in pr]
ys = [p["ratio_zero_minor_over_rho_geomean"] for p in pr]
lo = [p["ratio_zero_minor_over_rho_geomean"] - p["ci95_bootstrap"][0] for p in pr]
hi = [p["ci95_bootstrap"][1] - p["ratio_zero_minor_over_rho_geomean"] for p in pr]
ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o-", color="#9467bd", capsize=3,
            label="geometric mean over 8 targets, 95% bootstrap CI")
ax.axhline(1.0, color="gray", ls=":", label="parity with rho.negation")
ax.set_yscale("log")
ax.set_xlabel("curve (registry ladder bits)")
ax.set_ylabel("zero-minor / rho.negation field ops")
ax.set_title("(c) Paired ratio per target", fontsize=10)
ax.legend(fontsize=6.5, loc="upper left")
ax.grid(alpha=0.3, which="both")
fig.suptitle("lit2026 W12: zero-minor guess-and-determine (arXiv 2607.09814) vs rho.negation on the registered "
             f"prime ladder, run {run.name}. Intervals: 95% bootstrap.", fontsize=9)
fig.tight_layout()
fig.savefig(run / "scaling.svg")

ds = run / "dsweep.json"
if ds.exists():
    d = json.loads(ds.read_text())
    pcs = d["per_curve"]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    labels = [p["method"].split("(")[1].rstrip(")") for p in pcs]
    ax.bar(range(len(pcs)), [math.log2(p["ops_median"]) for p in pcs], color="#d62728", alpha=0.7)
    for i, p in enumerate(pcs):
        ax.vlines(i, math.log2(p["ops_min"]), math.log2(p["ops_max"]), color="k")
    ax.set_xticks(range(len(pcs)), labels)
    ax.set_ylabel("log2 field ops (median; bar = min-max)")
    ax.set_title(f"Defect sweep on {pcs[0]['slug']} ({pcs[0]['runs']} targets per d), run {run.name}",
                 fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(run / "dsweep.svg")
