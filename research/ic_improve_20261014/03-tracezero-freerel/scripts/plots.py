"""Figures for item 03. Plots recorded values only (results.json, census_*.json,
trace_theory.jsonl); computes nothing new. Usage: python3 -I plots.py <item_dir>"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

D = Path(sys.argv[1])
res = json.loads((D / "results.json").read_text())
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
COL = {"inc": "#2a78d6", "F": "#eb6834", "R": "#1baf7a", "FR": "#eda100"}
LAB = {"inc": "incumbent binary (A/A)", "F": "F: trace filter", "R": "R: free relations", "FR": "FR: both"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "svg.fonttype": "none"})


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# Figure 1: paired ratios vs arm A, counted total / counted rank / wall total.
sizes = ["n41", "n53", "n61"]
panels = [("counted_total", "complete solve, counted unit"),
          ("counted_rank", "rank stage, counted unit"),
          ("wall_total", "complete solve, wall ms")]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=False)
for ax, (metric, title) in zip(axes, panels):
    for i, size in enumerate(sizes):
        for j, arm in enumerate(["inc", "F", "R", "FR"]):
            a = res[size].get(arm)
            if not a:
                continue
            m = a[metric]
            x = i + (j - 1.5) * 0.17
            ax.errorbar(x, m["ratio"], yerr=[[m["ratio"] - m["ci95"][0]], [m["ci95"][1] - m["ratio"]]],
                        fmt="o", ms=5, color=COL[arm], ecolor=COL[arm], elinewidth=2, capsize=0,
                        label=LAB[arm] if i == 0 else None)
            if metric == "counted_rank" and arm in ("R", "FR") and arm == "R":
                ax.plot([x - 0.12, x + 0.12 + 0.17], [a["expected_rank_ratio_K_minus_F_over_K"]] * 2,
                        color=INK, lw=1, ls="--", label="expected (K - free rank)/K" if i == 0 else None)
    ax.axhline(1, color=INK2, lw=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels([f"n={s[1:]}\nK={res[s]['K']}" for s in sizes])
    ax.set_title(title, color=INK, fontsize=9.5)
    ax.set_ylabel("arm / A (ratio of means, 95% CI)")
    style(ax)
axes[0].legend(frameon=False, fontsize=8, loc="lower left")
axes[1].legend(frameon=False, fontsize=8, loc="lower left")
reps = ", ".join(f"{s}: {res[s]['A_replicates']} seeds" for s in sizes)
fig.suptitle(f"Item 03 paired arms against the incumbent copy A ({reps}; 8 targets each)", color=INK, fontsize=10)
fig.tight_layout()
fig.savefig(D / "fig_paired.svg")
fig.savefig(D / "fig_paired.png", dpi=150)

# Figure 2: (left) filter firing rate by target law; (right) free rank vs K.
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.4))
rows = [json.loads(l) for l in (D / "data" / "trace_theory.jsonl").read_text().splitlines()]
labels, vals = [], []
for r in rows:
    if r["check"] == "control_uniform_abscissae":
        labels.append(f"uniform curve x\nn={r['n']} a={r['a']}")
        vals.append(r["fraction"])
for r in rows:
    if r["check"] == "control_subgroup" and r["n"] == 61:
        labels.append("subgroup [k]G\nn=61")
        vals.append(r["filter_fires_kG"] / r["samples"])
inc_rej = sum(res[s]["F"]["filter_rejections_total"] for s in sizes if "F" in res[s])
labels.append("incumbent pipeline\n(all F-arm attempts)")
vals.append(0.0 if inc_rej == 0 else float("nan"))
a1.bar(range(len(vals)), vals, color=COL["F"], width=0.6)
for i, v in enumerate(vals):
    a1.text(i, v + 0.01, f"{v:.3f}", ha="center", color=INK, fontsize=8)
a1.set_xticks(range(len(vals)))
a1.set_xticklabels(labels, fontsize=7)
a1.set_ylabel("fraction with Tr(sqrt(b)/X) = 1")
a1.set_ylim(0, 0.6)
a1.set_title("Trace filter fires only outside 2E(F)", color=INK, fontsize=9.5)
style(a1)
cens = {}
for f in sorted((D / "data").glob("census_n*_K*.json")):
    c = json.loads(f.read_text())
    cens.setdefault(c["n"], []).append((c["K"], c["free_rank"]))
mk = {41: "o", 53: "s", 61: "^"}
for j, (n, pts) in enumerate(sorted(cens.items())):
    pts.sort()
    a2.plot([p[0] for p in pts], [100 * p[1] / p[0] for p in pts], marker=mk[n], ms=6, lw=2,
            color=list(COL.values())[j], label=f"n={n}")
    for k, fr in pts:
        a2.annotate(f"{fr}", (k, 100 * fr / k), textcoords="offset points", xytext=(4, 4), fontsize=7, color=INK2)
a2.set_yscale("log")
a2.set_xlabel("orbit columns K")
a2.set_ylabel("free rank / K (%)  [label: free rank]")
a2.set_title("Free relations harvested from the root table", color=INK, fontsize=9.5)
a2.legend(frameon=False, fontsize=8)
style(a2)
fig.tight_layout()
fig.savefig(D / "fig_trace_free.svg")
fig.savefig(D / "fig_trace_free.png", dpi=150)
