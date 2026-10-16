"""Figure for round 3 item 01. Plots recorded values only: summary.json (ic_improve3_01_analyze
analyze) and diag.json (ic_improve3_01_analyze diag). No value is computed here beyond axis scaling.

usage: python3 -I plot.py <summary.json> <diag.json> <r2.json> <out_stem> <incumbent> <arm> [arm ...]
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sumf, diagf, r2f, stem, INC = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
arms = [INC] + sys.argv[6:]
S = json.load(open(sumf))
D = json.load(open(diagf))
D["r2"] = json.load(open(r2f))
pa = S["per_arm"]
palette = ["#2b8a3e", "#4c6ef5", "#e8590c", "#ae3ec9", "#1098ad"]
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

# (b) local exponents of the whole solve and the relation phase
for i, a in enumerate(arms):
    for key, fmt, lab in [("local_cost", "-o", "whole"), ("local_rel_cost", "--s", "relation")]:
        loc = pa[a].get(key, [])
        x = [l[0] + 0.5 + 0.05 * (i - len(arms) / 2) for l in loc]
        ax[1].errorbar(x, [l[1] for l in loc], yerr=[[l[1] - l[2] for l in loc], [l[3] - l[1] for l in loc]],
                       fmt=fmt, ms=3, capsize=2, color=col[a], alpha=1.0 if lab == "whole" else 0.55,
                       label=f"{a} {lab}")
ax[1].axhline(4 / 3, ls=":", color="#868e96", label="4/3")
ax[1].set_xlabel("interval midpoint, log2 q")
ax[1].set_ylabel("local exponent d ln cost / d ln q")
ax[1].set_title("(b) Local exponents, whole (solid) and relation (dashed), 95% CI")
ax[1].legend(fontsize=6.5, ncol=2)
ax[1].grid(alpha=0.3)

# (c) paired ratios against the incumbent
for i, a in enumerate(arms[1:]):
    sz = pa[a]["sizes"]
    x = [s["bits"] + 0.06 * (i - len(arms) / 2) for s in sz]
    r = [s["ratio"] for s in sz]
    ax[2].errorbar(x, [v[0] for v in r], yerr=[[v[0] - v[1] for v in r], [v[2] - v[0] for v in r]],
                   fmt="-o", ms=3, capsize=2, color=col[a], label=f"{a} whole")
    ax[2].plot(x, [s["ratio_rel"][0] for s in sz], "--", lw=1, color=col[a], alpha=0.85, label=f"{a} relation")
ax[2].axhline(1, color="#343a40", lw=0.8)
ax[2].set_xlabel("log2 q")
ax[2].set_ylabel(f"paired ratio vs {INC}")
ax[2].set_title("(c) Paired ratios: whole (solid, 95% CI) and relation phase (dashed)")
ax[2].legend(fontsize=7)
ax[2].grid(alpha=0.3)

# (d) diagnosis: cost per trial by primitive, incumbent vs R1 (training curves, kappa 2)
ds = D["sizes"]
b = [d["bits"] for d in ds]
w = 0.38
parts = [("walk", "cand_walk", "#4c6ef5", "walk addition (Cantor)"),
         ("split", "cand_split", "#e8590c", "splitting test"),
         ("roots", "cand_roots", "#2b8a3e", "root finding")]
bot_i = [0.0] * len(b)
bot_c = [0.0] * len(b)
for ki, kc, c, lab in parts:
    vi = [d[ki] for d in ds]
    vc = [d[kc] for d in ds]
    ax[3].bar([x - w / 2 for x in b], vi, w, bottom=bot_i, color=c, alpha=0.9, label=f"{lab}")
    ax[3].bar([x + w / 2 for x in b], vc, w, bottom=bot_c, color=c, alpha=0.45, hatch="//")
    bot_i = [x + y for x, y in zip(bot_i, vi)]
    bot_c = [x + y for x, y in zip(bot_c, vc)]
ax[3].set_xlabel("log2 q (left bar: incumbent; right, hatched: R1)")
ax[3].set_ylabel("counted units per trial")
ax[3].set_title("(d) Cost per trial by primitive (training curves 5-7, kappa 2)")
ax[3].legend(fontsize=7)
ax[3].grid(alpha=0.3, axis="y")

# (e) diagnosis: trials needed
ax[4].plot(b, [d["t_kappa_q43"] for d in ds], "-o", color="#343a40", label="trials x kappa / q^(4/3)")
ax[4].plot(b, [d["mean_degree"] / 4 for d in ds], "--s", color="#ae3ec9", ms=3,
           label="(2LP-graph mean degree 2E/V) / 4")
ax[4].plot(b, [d["t_model_over_t"] for d in ds], ":^", color="#1098ad", ms=3,
           label="independent-trial model / measured trials")
ax[4].set_xlabel("log2 q")
ax[4].set_title("(e) Trials needed (training curves 5-7, kappa 2)")
ax[4].legend(fontsize=7)
ax[4].grid(alpha=0.3)

# (f) R2 rule: kappa schedules
if "r2" in D:
    r2 = D["r2"]
    xb = list(range(12, 23))
    ax[5].plot(xb, [2.718281828 ** (0.5533 - 0.0806 * (x - 20)) for x in xb], "-", color="#343a40",
               label="item 06 kappa(b) (incumbent)")
    ax[5].plot(xb, [2.718281828 ** (r2["alpha"] + r2["beta"] * (x - 20)) for x in xb], "-", color=palette[1],
               label="kappa_R2(b) (frozen, amendment C1)")
    ax[5].plot([s["bits"] for s in r2["sizes"]], [s["kappa_star"] for s in r2["sizes"]], "o", color=palette[1],
               label="kappa*_R1(b), training")
    ax[5].set_xlabel("log2 q")
    ax[5].set_ylabel("kappa")
    ax[5].set_title("(f) Factor-base constant: incumbent and R2 schedules")
    ax[5].legend(fontsize=7)
    ax[5].grid(alpha=0.3)
else:
    ax[5].axis("off")

fig.tight_layout()
fig.savefig(f"{stem}.svg")
fig.savefig(f"{stem}.png", dpi=150)
