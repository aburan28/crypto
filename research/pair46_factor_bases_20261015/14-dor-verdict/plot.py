"""Layout only: plots values that dorv wrote to out/results.json."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(here, "out", "results.json")))
COL = {"standard": "#1f77b4", "random": "#7f7f7f", "nb_window": "#2ca02c",
       "nb_break": "#d62728", "subfield": "#9467bd"}

def save(fig, name):
    fig.savefig(os.path.join(here, name + ".svg"))
    fig.savefig(os.path.join(here, name + ".png"), dpi=150)
    plt.close(fig)

# Fig 1: D* vs dim per n, by structure, 95% CI
ns = sorted({c["n"] for c in R["cells"]})
fig, axes = plt.subplots(2, 4, figsize=(14, 6.5), sharey=True)
for ax, n in zip(axes.flat, ns):
    for s, col in COL.items():
        cs = sorted([c for c in R["cells"] if c["n"] == n and c["structure"] == s], key=lambda c: c["dim"])
        if not cs:
            continue
        x = [c["dim"] + {"standard": -0.15, "random": -0.05, "nb_window": 0.05, "nb_break": 0.15, "subfield": 0}[s] for c in cs]
        y = [c["d_star_mean"] for c in cs]
        lo = [c["d_star_mean"] - c["ci95"][0] for c in cs]
        hi = [c["ci95"][1] - c["d_star_mean"] for c in cs]
        ax.errorbar(x, y, yerr=[lo, hi], fmt="o-", ms=3, lw=1, capsize=2, color=col, label=s)
    ax.set_title(f"n = {n}")
    ax.set_xlabel("factor-base dim n'")
    ax.grid(alpha=0.3)
axes.flat[0].set_ylabel("mean D* (95% cluster bootstrap)")
axes.flat[4].set_ylabel("mean D* (95% cluster bootstrap)")
axes.flat[-1].axis("off")
axes.flat[1].legend(fontsize=7)
fig.suptitle("Recorded campaign D* (m=2, unsat targets) by structure and dimension; data: factor_base_dor_campaign_20261013/raw")
fig.tight_layout()
save(fig, "fig_dstar_vs_dim")

# Fig 2: D* against the exact invariant slack (l + r)/n
fig, ax = plt.subplots(figsize=(7, 4.8))
for s, col in COL.items():
    cs = [c for c in R["cells"] if c["structure"] == s]
    ax.scatter([(c["l_dim_V1_plus_V2"] + c["r_dim_V1V2_max_over_seeds"]) / c["n"] for c in cs],
               [c["d_star_mean"] for c in cs], s=16, color=col, label=s, alpha=0.8)
ax.axvline(1.0, color="k", lw=0.8, ls="--")
ax.text(1.01, 2.1, "l + r = n", fontsize=8)
ax.set_xlabel("(dim(V1+V2) + dim(V1·V2)) / n   [exact invariants of the factor base]")
ax.set_ylabel("cell mean D*")
ax.set_title("D* tracks the degree-2 rank budget l + r, not the structure label")
ax.legend(fontsize=7); ax.grid(alpha=0.3)
fig.tight_layout()
save(fig, "fig_dstar_vs_invariant")

# Fig 3: Semaev vs invariant-matched null
if "n2_invariant_null" in R:
    rows = R["n2_invariant_null"]["rows"]
    fig, ax = plt.subplots(figsize=(6, 5.5))
    for s, col in COL.items():
        rs = [r for r in rows if r["structure"] == s]
        ax.errorbar([r["null_mean"] for r in rs], [r["semaev_mean"] for r in rs],
                    yerr=[[r["diff"] - r["ci95"][0] for r in rs], [r["ci95"][1] - r["diff"] for r in rs]],
                    fmt="o", ms=4, capsize=2, color=col, label=s)
    lim = [1.9, 5.2]
    ax.plot(lim, lim, "k-", lw=0.8)
    ax.fill_between(lim, [v - 0.25 for v in lim], [v + 0.25 for v in lim], color="k", alpha=0.08, label="±0.25 margin")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("N2 null mean D* (random bilinear systems, matched n, l, q; 200 each)")
    ax.set_ylabel("Semaev cell mean D* (bars: 95% CI of difference)")
    ax.set_title("Invariant null N2 vs recorded Semaev D*")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_null_vs_semaev")
