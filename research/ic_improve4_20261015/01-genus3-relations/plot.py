"""Figure for round 4 item 01. Plots recorded values only: summary.json (ic_improve4_01_analyze
analyze), pertrial.json (ic_improve4_01_analyze pertrial) and tune.json (ic_improve4_01_analyze tune).
No value is computed here beyond axis scaling and the schedule lines from the frozen coefficients.

usage: python3 -I plot.py <summary.json> <pertrial.json> <tune.json> <out_stem> <incumbent> <arm> [arm ...]
"""
import json
import math
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sumf, ptf, tunef, stem, INC = sys.argv[1:6]
arms = [INC] + sys.argv[6:]
S = json.load(open(sumf))
PT = json.load(open(ptf))
T = json.load(open(tunef))
pa = S["per_arm"]
palette = ["#2b8a3e", "#4c6ef5", "#e8590c", "#ae3ec9", "#1098ad", "#c92a2a"]
col = {INC: "#343a40"}
for i, a in enumerate(arms[1:]):
    col[a] = palette[i % len(palette)]

fig, ax = plt.subplots(2, 3, figsize=(19, 10.5))
ax = ax.ravel()

# (a) whole cost (solid) and relation phase (dashed) vs q
for a in arms:
    sz = pa[a]["sizes"]
    q = [s["q"] for s in sz]
    ax[0].plot(q, [s["gm_cost"] for s in sz], "-o", ms=3.5, color=col[a], label=f"{a} whole")
    ax[0].plot(q, [s["gm_rel"] for s in sz], "--", lw=1.1, color=col[a], alpha=0.85)
sz = pa[INC]["sizes"]
q0, c0 = sz[-1]["q"], sz[-1]["gm_cost"]
qs = [s["q"] for s in sz]
ax[0].plot(qs, [c0 * (x / q0) ** (4 / 3) for x in qs], ":", color="#868e96", label="slope 4/3 (reference)")
ax[0].set_xscale("log", base=2)
ax[0].set_yscale("log")
ax[0].set_xlabel("q")
ax[0].set_ylabel("cost (F_q-mult. equivalents, c_ell = 1)")
ax[0].set_title("(a) Whole solve (solid) and relation phase (dashed)")
ax[0].legend(fontsize=7)
ax[0].grid(alpha=0.3, which="both")

# (b) local exponents of the whole solve
for i, a in enumerate(arms):
    loc = pa[a].get("local_cost", [])
    x = [l[0] + 0.5 + 0.06 * (i - len(arms) / 2) for l in loc]
    ax[1].errorbar(x, [l[1] for l in loc], yerr=[[l[1] - l[2] for l in loc], [l[3] - l[1] for l in loc]],
                   fmt="-o", ms=3, capsize=2, color=col[a], label=a)
ax[1].axhline(4 / 3, ls=":", color="#868e96", label="4/3")
ax[1].set_xlabel("interval midpoint, log2 q")
ax[1].set_ylabel("local exponent d ln cost / d ln q")
ax[1].set_title("(b) Local exponents of the whole solve, 95% CI")
ax[1].legend(fontsize=7)
ax[1].grid(alpha=0.3)

# (c) paired ratios against the incumbent (whole, CI; relation dashed; LA dotted)
for a in arms[1:]:
    sz = pa[a]["sizes"]
    b = [s["bits"] for s in sz]
    r = [s["ratio"] for s in sz]
    ax[2].errorbar(b, [x[0] for x in r], yerr=[[x[0] - x[1] for x in r], [x[2] - x[0] for x in r]],
                   fmt="-o", ms=3.5, capsize=2, color=col[a], label=f"{a} whole")
    ax[2].plot(b, [s["ratio_rel"][0] for s in sz], "--", lw=1, color=col[a], alpha=0.8)
    ax[2].plot(b, [s["ratio_la"][0] for s in sz], ":", lw=1.2, color=col[a], alpha=0.8)
ax[2].axhline(1.0, color="#343a40", lw=0.8)
ax[2].set_xlabel("log2 q")
ax[2].set_ylabel(f"paired ratio to {INC}")
ax[2].set_title("(c) Paired ratios: whole (solid, CI), relation (dashed), LA (dotted)")
ax[2].legend(fontsize=7)
ax[2].grid(alpha=0.3)

# (d) cost per trial by primitive: R1 | C1 | C2 | C1+C2 (stacked)
rows = PT["rows"]
w = 0.2
parts = [("walk", "#495057"), ("split", "#4c6ef5"), ("roots", "#e8590c")]
vkeys = [("R1", "r1_walk", "r1_split", "r1_roots"), ("C1", "r1_walk", "c1_split", "c1_roots"),
         ("C2", "c2_walk", "r1_split", "r1_roots"), ("C1+C2", "c2_walk", "c1_split", "c1_roots")]
for vi, (vn, kw, ks, kr) in enumerate(vkeys):
    xs = [r["bits"] + (vi - 1.5) * w for r in rows]
    bottom = [0.0] * len(rows)
    for (pn, pc), k in zip(parts, [kw, ks, kr]):
        h = [r[k] for r in rows]
        ax[3].bar(xs, h, width=w * 0.95, bottom=bottom, color=pc, alpha=[1.0, 0.75, 0.55, 0.35][vi],
                  label=pn if vi == 0 else None, edgecolor="white", lw=0.3)
        bottom = [b0 + h0 for b0, h0 in zip(bottom, h)]
ax[3].set_xlabel("log2 q (bars: R1, C1, C2, C1+C2)")
ax[3].set_ylabel("counted units per trial")
fits = PT["fits"]
ax[3].set_title("(d) Cost per trial; A + B log2 q: " + ", ".join(
    f"{k} B={fits[k]['B']:.2f}" for k in ["R1", "C1", "C2", "C1+C2"]), fontsize=9)
ax[3].legend(fontsize=7)
ax[3].grid(alpha=0.3, axis="y")

# (e) solver dimension per factor-base element and nnz per column
for a in arms:
    sz = pa[a]["sizes"]
    b = [s["bits"] for s in sz]
    ax[4].plot(b, [s["la_cols"] / s["nfb"] for s in sz], "-o", ms=3, color=col[a], label=f"{a} cols/nfb")
    ax4b = ax[4]
    ax4b.plot(b, [s["la_nnz"] / s["la_cols"] / 20 for s in sz], "--", lw=1, color=col[a], alpha=0.7)
ax[4].set_xlabel("log2 q")
ax[4].set_ylabel("solver columns / factor base (solid); nnz per column / 20 (dashed)")
ax[4].set_title("(e) Matrix handed to block Wiedemann")
ax[4].legend(fontsize=7)
ax[4].grid(alpha=0.3)

# (f) kappa schedules
bb = [12 + 0.1 * i for i in range(101)]
ax[5].plot(bb, [math.exp(0.2669 - 0.0941 * (b - 20)) for b in bb], "-", color="#343a40", label="kappa_R2 (incumbent)")
c4 = T["c4"]
ax[5].plot(bb, [math.exp(c4["alpha"] + c4["beta"] * (b - 20)) for b in bb], "-", color="#c92a2a",
           label=f"kappa_C4 (frozen, {c4['arm']})")
ax[5].plot([p[0] for p in c4["points"]], [math.exp(p[1]) for p in c4["points"]], "o", color="#c92a2a",
           label="training minimisers")
ax[5].set_yscale("log")
ax[5].set_xlabel("log2 q")
ax[5].set_ylabel("kappa (factor base = kappa q^(2/3))")
ax[5].set_title("(f) kappa schedules")
ax[5].legend(fontsize=7)
ax[5].grid(alpha=0.3, which="both")

fig.tight_layout()
fig.savefig(stem + ".svg")
fig.savefig(stem + ".png", dpi=110)
print("wrote", stem + ".svg", stem + ".png")
