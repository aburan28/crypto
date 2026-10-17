"""Plot recorded values from results/run-*.jsonl. Layout only.
Usage: python3 plot.py <results.jsonl> <outdir>"""
import json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]
cases = [r for r in map(json.loads, open(src)) if r["rec"] == "case"]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1a1a19", "#6b6a63", "#e4e3dc"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(7.2, 3.4))
xs = range(len(cases))
ax.scatter(xs, [r["order_E0"] for r in cases], s=150, facecolors="none", edgecolors=ORANGE, linewidths=2, label="#E0(K), counted on the elliptic curve over K", zorder=2)
ax.scatter(xs, [r["jac_order"] for r in cases], s=36, color=BLUE, label="#Jac(C)(F_q), from point counts of the F_q quartic", zorder=3)
for x, r in zip(xs, cases): ax.text(x, r["jac_order"] * 1.22, f'{r["jac_order"]}\n= 4·{r["ell"]}', ha="center", va="bottom", color=INK, fontsize=7.5)
ax.set_yscale("log"); ax.set_ylim(80, 4500); ax.set_ylabel("group order")
ax.set_xticks(list(xs)); ax.set_xticklabels([f'q={r["q"]}\n{r["kind"][:5]}\n{r["u_idx"][1]}' for r in cases], fontsize=7.5)
ax.grid(axis="y", color=GRID, linewidth=0.8); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left")
ax.set_xlabel("case: field, type of u, index of its constant term")
fig.savefig(f"{out}/orders.svg", bbox_inches="tight"); fig.savefig(f"{out}/orders.png", dpi=150, bbox_inches="tight")
