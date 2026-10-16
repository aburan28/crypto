#!/usr/bin/env python3
"""Figure for item 02. Plots only values computed by the Rust fit tool
(runs/fit_oracle.jsonl, runs/fit_step.jsonl): geometric-mean paired ratios
per n with 95 % bootstrap CIs.

python3 -I plot.py  ->  fig_ratios.svg, fig_ratios.png
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    return {r["ratio"]: r for r in (json.loads(l) for l in open(os.path.join(HERE, "runs", name)))}


oracle = load("fit_oracle.jsonl")
step = load("fit_step.jsonl")

STYLE = {
    "I/C": ("#1f5fbf", "o", "I / C  (incumbent ÷ equivariant candidate)"),
    "I/F": ("#d9822b", "s", "I / F  (incumbent ÷ plain fall-row augmentation)"),
    "F/C": ("#2a9d6f", "^", "F / C  (what equivariance adds)"),
    "AA": ("#888888", "x", "I / I′  (A/A null)"),
    "POS": ("#9b3fbf", "D", "Sym-full ÷ C  (positive control)"),
}


def panel(ax, fits, keys, title, loc="upper left"):
    for i, k in enumerate(keys):
        if k not in fits:
            continue
        c, m, lab = STYLE[k]
        pts = fits[k]["per_n"]
        xs = [p["n"] + (i - len(keys) / 2) * 0.12 for p in pts]
        ys = [p["geomean"] for p in pts]
        lo = [p["geomean"] - p["lo"] for p in pts]
        hi = [p["hi"] - p["geomean"] for p in pts]
        s = fits[k].get("slope_vs_ln_n")
        if s:
            lab += f"; slope {s['b']:.2f} [{s['lo']:.2f}, {s['hi']:.2f}]"
        ax.errorbar(xs, ys, yerr=[lo, hi], color=c, marker=m, ms=5, lw=1.2, capsize=3, label=lab)
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel("n  (field degree of F_{2^n})")
    ax.set_ylabel("paired CPU ratio (> 1: right-hand arm cheaper)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc=loc)
    ax.grid(alpha=0.3, which="both")


fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
panel(axes[0], oracle, ["I/C", "I/F", "F/C", "AA"],
      "Whole oracle (barrel, m = 2, ℓ = 3): 12 targets per n, all gates pass")
panel(axes[1], step, ["I/C", "POS", "AA"],
      "Degree-3 Macaulay step: step gate FAILS (C misses a linear row)", loc="lower right")
fig.suptitle("Item 02, G2 equivariant pass vs incumbent: geometric means, 95 % bootstrap CI "
             "(runs/oracle.jsonl, runs/step_*.jsonl)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_ratios.svg"))
fig.savefig(os.path.join(HERE, "fig_ratios.png"), dpi=150)
print("wrote fig_ratios.svg/.png")
