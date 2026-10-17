"""Graphs for E19 (RESEARCH_GLV_INVARIANT_FACTOR_BASES.md §9.5–§9.6).

Lays out recorded values from experiments/27_quartic_pipeline_d4/p*s*.json;
computes nothing the Rust run did not record, except the registered
prediction lines (§9.5's formulas, drawn as references).

    python3 research/glv_k4_pipeline_d4_20261017/figures/plot.py
"""

import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
rows = []
for path in sorted(glob.glob(os.path.join(ROOT, "experiments/27_quartic_pipeline_d4/p*s*.json"))):
    rows.extend(json.load(open(path)))
rows.sort(key=lambda r: (r["p"], r["seed"]))
ids = ", ".join(f"p{r['p']}s{r['seed']}" for r in rows)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, name + ".svg"))
    fig.savefig(os.path.join(HERE, name + ".png"), dpi=150)
    plt.close(fig)


# 1. The distance to rho: S / rho per curve against log2 n, with the
#    registered prediction 0.584·C/p at d4-s5's C = 1.3e9 and its ×2 band.
fig, ax = plt.subplots(figsize=(6.4, 4.2))
bits = [r["bits"] for r in rows]
ax.scatter(bits, [r["ratio"] for r in rows], color="C0", zorder=3,
           label="measured S / rho (negation walks, 16 per curve)")
ps = sorted({r["p"] for r in rows})
pb = {r["p"]: r["bits"] for r in rows}
pred = [0.584 * 1.3e9 / p for p in ps]
ax.plot([pb[p] for p in ps], pred, "k--", label="G2.12 prediction 0.584·C/p, C = 1.3e9")
ax.fill_between([pb[p] for p in ps], [v / 2 for v in pred], [v * 2 for v in pred],
                color="0.85", label="×2 band (falsification bound)")
ax.set_yscale("log")
ax.set_xlabel("log2 n (prime subgroup order)")
ax.set_ylabel("S_IC / S_rho  (F_p multiplications)")
ax.set_title("E19: k = 4 index calculus, end to end, against matched rho")
ax.legend(fontsize=7)
ax.text(0.01, 0.01, f"curves: {ids}", transform=ax.transAxes, fontsize=5, color="0.4")
save(fig, "distance_to_rho")

# 2. Per-call solve cost: p10 / median / p90 per curve, d4-s5's band ±30 %.
fig, ax = plt.subplots(figsize=(6.4, 4.0))
x = list(range(len(rows)))
med = [r["call"]["median"] for r in rows]
lo = [r["call"]["median"] - r["call"]["p10"] for r in rows]
hi = [r["call"]["p90"] - r["call"]["median"] for r in rows]
ax.errorbar(x, med, yerr=[lo, hi], fmt="o", capsize=3, label="median, bars p10–p90")
ax.scatter(x, [r["call"]["mean"] for r in rows], marker="x", color="C3", label="mean", zorder=3)
ax.axhspan(1.2e9 * 0.7, 1.44e9 * 1.3, color="0.88", label="d4-s5 medians 1.2–1.44e9, ±30 % (G2.11)")
ax.set_xticks(x)
ax.set_xticklabels([f"{r['p']}/{r['seed']}" for r in rows], rotation=45, fontsize=7)
ax.set_xlabel("p / seed")
ax.set_ylabel("C per call (F_p multiplications)")
ax.set_title("E19: one weighted D4 solve inside the pipeline")
ax.legend(fontsize=7)
save(fig, "solve_cost")

# 3. Residuals per relation y, with G2.10's band.
fig, ax = plt.subplots(figsize=(6.4, 3.6))
y = [r["residuals"] / r["folded"]["relations"] for r in rows]
ax.bar(x, y, color="C2")
ax.axhspan(70, 160, color="0.88", zorder=0, label="G2.10 band [70, 160] (pooled)")
ax.axhline(106, color="k", ls="--", lw=1, label="E18 stand-in, 106")
ax.set_xticks(x)
ax.set_xticklabels([f"{r['p']}/{r['seed']}" for r in rows], rotation=45, fontsize=7)
ax.set_xlabel("p / seed")
ax.set_ylabel("targets per kept relation")
ax.set_title("E19: residuals per relation from the exact D4 oracle")
ax.legend(fontsize=7)
save(fig, "residuals_per_relation")
