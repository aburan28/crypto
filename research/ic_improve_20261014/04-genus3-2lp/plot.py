#!/usr/bin/env python3
"""Figures for item 04. Plots recorded values from summary.json only (no computation
beyond reading the analyzer's numbers). Usage: python3 -I plot.py results/run-20261014"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = sys.argv[1]
S = json.load(open(f"{run}/summary.json"))
ARMS = [("basic", "IC-basic (incumbent)", "#c0392b", "s"),
        ("ic1lp", "IC-1LP, kappa1=2 (middle)", "#8e44ad", "^"),
        ("ic2lp", "IC-2LP, kappa=2 (candidate)", "#1f6fb2", "o"),
        ("aa", "IC-2LP A/A (null)", "#7fb3d5", "x")]

fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
ax = axs[0]
for arm, lab, col, mk in ARMS:
    xs, ys = [], []
    for row in S["per_size"]:
        for r in row.get(f"{arm}_rows", []):
            xs.append(2 ** row["bits"]); ys.append(r["cost"])
    if xs:
        ax.scatter(xs, ys, s=18, marker=mk, color=col, label=lab, alpha=0.85)
xs, ys = [], []
for row in S["per_size"]:
    if "rho_model_neg_mean" in row:
        xs.append(2 ** row["bits"]); ys.append(row["rho_model_neg_mean"])
ax.plot(xs, ys, "--", color="#555", label="rho model, negation map (reference)")
ax.set_xscale("log", base=2); ax.set_yscale("log")
ax.set_xlabel("q"); ax.set_ylabel("cost per solve (F_q-mult equivalents)")
ax.set_title("End-to-end cost, every solve")
ax.legend(fontsize=7.5)
ax.grid(alpha=0.3, which="both")

ax = axs[1]
off = {"basic": -0.12, "ic1lp": -0.04, "ic2lp": 0.04, "aa": 0.12}
for arm, lab, col, mk in ARMS:
    le = S["local_exponents"].get(arm)
    if not le:
        continue
    xs = [(iv["from"] + iv["to"]) / 2 + off[arm] for iv in le["intervals"]]
    ys = [iv["local"] for iv in le["intervals"]]
    lo = [iv["local"] - iv["ci"][0] for iv in le["intervals"]]
    hi = [iv["ci"][1] - iv["local"] for iv in le["intervals"]]
    ax.errorbar(xs, ys, yerr=[lo, hi], fmt=mk + "-", color=col, ms=4, capsize=2, lw=1, label=lab)
for v, t in [(4 / 3, "4/3 (GTTD, g=3)"), (10 / 7, "10/7 (1LP)"), (2.0, "2")]:
    ax.axhline(v, color="#999", ls=":", lw=1)
    ax.text(10.1, v + 0.01, t, fontsize=7, color="#555")
ax.set_xlabel("log2 q (interval midpoint)"); ax.set_ylabel("local exponent d ln cost / d ln q")
ax.set_title("Local exponents, 95% bootstrap CI")
ax.grid(alpha=0.3)
ax.legend(fontsize=7.5, loc="lower right")

ax = axs[2]
for key, col, lab in [("ic2lp/basic", "#1f6fb2", "IC-2LP / IC-basic"),
                      ("ic1lp/basic", "#8e44ad", "IC-1LP / IC-basic"),
                      ("aa/ic2lp", "#7fb3d5", "A/A: IC-2LP' / IC-2LP")]:
    rows = S["paired_ratios"].get(key, [])
    xs = [r["bits"] for r in rows]; ys = [r["ratio"] for r in rows]
    lo = [r["ratio"] - r["ci"][0] for r in rows]; hi = [r["ci"][1] - r["ratio"] for r in rows]
    ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o-", color=col, ms=4, capsize=2, lw=1, label=lab)
ax.axhline(1, color="#999", ls=":")
ax.set_yscale("log")
ax.set_xlabel("log2 q"); ax.set_ylabel("paired cost ratio (geometric mean, 95% CI)")
ax.set_title("Paired ratio against the incumbent")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=7.5)
fig.suptitle("Item 04: genus-3 IC after cubic GHS, run-20261014 (3 curves per size, paired on instance, target and walk seed)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{run}/../../cost_vs_q.svg")
fig.savefig(f"{run}/../../cost_vs_q.png", dpi=150)
