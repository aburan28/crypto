"""Figure for round 2 item 06. Plots recorded values only: summary.json (ic_improve2_06_analyze analyze)
and train.json (ic_improve2_06_analyze train).

usage: python3 -I plot.py <run_dir> <train.json> <incumbent> <arm> [arm ...]
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

run, trainf, INC = sys.argv[1], sys.argv[2], sys.argv[3]
arms = [INC] + sys.argv[4:]
S = json.load(open(f"{run}/summary.json"))
T = json.load(open(trainf))
pa = S["per_arm"]
palette = ["#e8590c", "#2b8a3e", "#4c6ef5", "#ae3ec9", "#1098ad", "#c92a2a"]
col = {INC: "#343a40"}
for i, a in enumerate(arms[1:]):
    col[a] = palette[i % len(palette)]

fig, ax = plt.subplots(2, 2, figsize=(14, 10))
ax = ax.ravel()

# (a) whole cost and LA cost vs q
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
ax[0].set_title("(a) Whole solve (solid) and linear algebra (dashed), geometric means")
ax[0].legend(fontsize=7)
ax[0].grid(alpha=0.3, which="both")

# (b) local exponents of the whole solve
for i, a in enumerate(arms):
    loc = pa[a].get("local_cost", [])
    x = [l[0] + 0.5 + 0.05 * (i - len(arms) / 2) for l in loc]
    ax[1].errorbar(x, [l[1] for l in loc], yerr=[[l[1] - l[2] for l in loc], [l[3] - l[1] for l in loc]],
                   fmt="-o", ms=3, capsize=2, color=col[a], label=a)
ax[1].axhline(4 / 3, ls=":", color="#868e96", label="4/3")
ax[1].set_xlabel("interval midpoint, log2 q")
ax[1].set_ylabel("local exponent d ln cost / d ln q")
ax[1].set_title("(b) Local exponents of the whole solve (95% bootstrap CI)")
ax[1].legend(fontsize=7)
ax[1].grid(alpha=0.3)

# (c) paired ratios against the incumbent
for i, a in enumerate(arms[1:]):
    sz = pa[a]["sizes"]
    x = [s["bits"] + 0.06 * (i - len(arms) / 2) for s in sz]
    r = [s["ratio"] for s in sz]
    ax[2].errorbar(x, [v[0] for v in r], yerr=[[v[0] - v[1] for v in r], [v[2] - v[0] for v in r]],
                   fmt="-o", ms=3, capsize=2, color=col[a], label=f"{a} whole")
    ax[2].plot(x, [s["ratio_la"][0] for s in sz], "--", lw=1, color=col[a], alpha=0.8)
ax[2].axhline(1, color="#343a40", lw=0.8)
ax[2].set_xlabel("log2 q")
ax[2].set_ylabel(f"paired ratio vs {INC}")
ax[2].set_title("(c) Paired ratios: whole (solid, 95% CI) and linear algebra (dashed)")
ax[2].legend(fontsize=7)
ax[2].grid(alpha=0.3)

# (d) training: cost vs kappa per size, kappa*(b) and the frozen schedule
ax3 = ax[3]
cm = plt.get_cmap("viridis")
for i, s in enumerate(T["sizes"]):
    g = sorted(s["grid"])
    c = cm(i / max(1, len(T["sizes"]) - 1))
    base = min(v for _, v in g)
    ax3.plot([k for k, _ in g], [v / base for _, v in g], "-o", ms=3, color=c,
             label=f"b = {s['bits']} (cx {s['cx']}), kappa* = {s['kappa_star']:.2f}")
    ax3.axvline(s["kappa_star"], color=c, lw=0.8, ls=":")
ax3.set_xscale("log", base=2)
ax3.set_xlabel("kappa (factor base = kappa q^(2/3))")
ax3.set_ylabel("training whole cost / minimum at that size")
import math
sched = ", ".join(f"{b}: {math.exp(T['alpha'] + T['beta'] * (b - 20)):.2f}" for b in (12, 16, 20, 21))
ax3.set_title(f"(d) Training curves 5-7: kappa(b) = {math.exp(T['alpha']):.3f} x {math.exp(T['beta']):.3f}^(b-20)\n(b {sched})", fontsize=10)
ax3.legend(fontsize=7)
ax3.grid(alpha=0.3, which="both")

fig.tight_layout()
fig.savefig("cost_vs_q.svg")
fig.savefig("cost_vs_q.png", dpi=150)
