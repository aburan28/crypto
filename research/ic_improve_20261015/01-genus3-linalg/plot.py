"""Figure for round 2 item 01. Plots summary.json only (written by ic_improve2_01_analyze).

usage: python3 -I plot.py <run_dir> [arm ...]
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = sys.argv[1]
S = json.load(open(f"{run}/summary.json"))
pa = S["per_arm"]
arms = sys.argv[2:] or [a for a in S["arms"] if "sizes" in pa.get(a, {})]
INC = "main:inc"
palette = ["#4c6ef5", "#e8590c", "#2b8a3e", "#ae3ec9", "#1098ad", "#c92a2a", "#5c940d", "#862e9c"]
col = {INC: "#343a40"}
for i, a in enumerate([a for a in arms if a != INC]):
    col[a] = palette[i % len(palette)]

fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))

# left: whole cost and LA cost vs q
for a in arms:
    sz = pa[a]["sizes"]
    q = [s["q"] for s in sz]
    ax[0].plot(q, [s["gm_cost"] for s in sz], "-o", ms=3.5, color=col[a], label=f"{a} whole")
    ax[0].plot(q, [s["gm_la"] for s in sz], "--", lw=1.1, color=col[a], alpha=0.8)
sz = pa[INC]["sizes"]
q0, c0 = sz[-1]["q"], sz[-1]["gm_cost"]
qs = [s["q"] for s in sz]
ax[0].plot(qs, [c0 * (x / q0) ** (4 / 3) for x in qs], ":", color="#868e96", label="slope 4/3 (reference)")
ax[0].set_xscale("log", base=2)
ax[0].set_yscale("log")
ax[0].set_xlabel("q")
ax[0].set_ylabel("cost (F_q-mult. equivalents, c_ell = 1)")
ax[0].set_title("Whole solve (solid) and linear algebra (dashed)")
ax[0].legend(fontsize=7)
ax[0].grid(alpha=0.3, which="both")

# middle: local exponents (whole)
for i, a in enumerate(arms):
    loc = pa[a].get("local_cost", [])
    if not loc:
        continue
    x = [l[0] + 0.5 + 0.04 * (i - len(arms) / 2) for l in loc]
    y = [l[1] for l in loc]
    lo = [l[1] - l[2] for l in loc]
    hi = [l[3] - l[1] for l in loc]
    ax[1].errorbar(x, y, yerr=[lo, hi], fmt="-o", ms=3, capsize=2, color=col[a], label=a)
ax[1].axhline(4 / 3, ls=":", color="#868e96", label="4/3")
ax[1].set_xlabel("interval midpoint, log2 q")
ax[1].set_ylabel("local exponent d ln cost / d ln q")
ax[1].set_title("Local exponents of the whole solve (95% CI)")
ax[1].legend(fontsize=7)
ax[1].grid(alpha=0.3)

# right: paired ratios vs incumbent
for i, a in enumerate(arms):
    if a == INC:
        continue
    sz = pa[a]["sizes"]
    x = [s["bits"] + 0.05 * (i - len(arms) / 2) for s in sz]
    r = [s["ratio"] for s in sz]
    ax[2].errorbar(x, [v[0] for v in r], yerr=[[v[0] - v[1] for v in r], [v[2] - v[0] for v in r]],
                   fmt="-o", ms=3, capsize=2, color=col[a], label=f"{a} whole")
    rl = [s["ratio_la"] for s in sz]
    ax[2].plot(x, [v[0] for v in rl], "--", lw=1, color=col[a], alpha=0.8)
ax[2].axhline(1, color="#343a40", lw=0.8)
ax[2].set_yscale("log")
ax[2].set_xlabel("log2 q")
ax[2].set_ylabel("paired ratio vs main:inc")
ax[2].set_title("Paired ratios: whole (solid, CI) and LA (dashed)")
ax[2].legend(fontsize=7)
ax[2].grid(alpha=0.3, which="both")

fig.tight_layout()
fig.savefig("cost_vs_q.svg")
fig.savefig("cost_vs_q.png", dpi=130)
