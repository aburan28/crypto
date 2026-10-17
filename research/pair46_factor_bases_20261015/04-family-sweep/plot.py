"""Layout only: plots recorded values from curves.csv / results.json (no new computation of results)."""
import csv, json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1] if len(sys.argv) > 1 else "."
rows = list(csv.DictReader(open(f"{d}/curves.csv")))
res = json.load(open(f"{d}/results.json"))
arms = ["F", "N1", "N2", "PC"]
col = {"F": "#c0392b", "N1": "#2e86c1", "N2": "#7f8c8d", "PC": "#27ae60"}
lab = {"F": "F family\np=5t²+19t+19", "N1": "N1 D=-19\nnon-family", "N2": "N2 random\nprime order", "PC": "PC planted\n(+6.1% M1)"}
mets = [("M1", "M1_norm_decomp_yield", "M1 normalized Semaev m=2 yield"), ("M2", "M2_relation_rate", "M2 relation rate (>=1 decomp.)"), ("M3", "M3_smallx_yield", "M3 small-x yield |F_x|/(B/2)")]

fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, (k, key, title) in zip(axs, mets):
    for i, a in enumerate(arms):
        v = [float(r[k]) for r in rows if r["arm"] == a]
        xs = [i + 0.25 * ((j * 0.6180339) % 1 - 0.5) for j in range(len(v))]
        ax.scatter(xs, v, s=6, alpha=0.35, color=col[a])
        m = sum(v) / len(v)
        # recorded bootstrap CI vs N2 (for F/N1/PC) shown as diff interval added to N2 mean
        ax.hlines(m, i - 0.3, i + 0.3, color="k", lw=2)
    n2 = [t for t in res["tests"] if t["band"] == "all" and t["metric"] == key and t["y"] == "N2"]
    for t in n2:
        i = arms.index(t["x"])
        lo, hi = t["mean_y"] + t["diff_ci95"][0], t["mean_y"] + t["diff_ci95"][1]
        ax.vlines(i + 0.36, lo, hi, color="k", lw=1.2)
        ax.text(i, ax.get_ylim()[1] if False else max(float(r[k]) for r in rows), f"p={t['p_perm']:.3g}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(4)); ax.set_xticklabels([lab[a] for a in arms], fontsize=7)
    ax.set_title(title, fontsize=9)
fig.suptitle(f"04-family-sweep: per-curve metrics, T={res['targets_per_curve']} targets/curve, bits 20-30, seed {res['seed']} "
             f"(n F={res['counts']['F']}, N1={res['counts']['N1']}, N2={res['counts']['N2']}, PC={res['counts']['PC']}); "
             "bar = arm mean; whisker = N2 mean + 95% bootstrap CI of (arm − N2); p = stratified permutation vs N2", fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(f"{d}/fig_metrics.svg"); fig.savefig(f"{d}/fig_metrics.png", dpi=150)

# effect sizes (ratio of means with CI) for primary + PC tests
fig, ax = plt.subplots(figsize=(8, 4.2))
sel = [t for t in res["tests"] if t["band"] == "all" and t["role"] in ("primary", "positive_control", "secondary")]
for i, t in enumerate(sel):
    c = col[t["x"]] if t["x"] != "N1" or t["y"] != "N2" else "#2e86c1"
    ax.errorbar(t["ratio"], i, xerr=[[t["ratio"] - t["ratio_ci95"][0]], [t["ratio_ci95"][1] - t["ratio"]]], fmt="o", color=c, capsize=3)
    ax.text(t["ratio_ci95"][1] + 0.003, i, f"p={t['p_perm']:.3g}" + (" *" if t["significant"] else ""), va="center", fontsize=7)
ax.axvline(1, color="k", lw=0.8, ls="--")
ax.axvline(1.0609, color=col["PC"], lw=0.8, ls=":")
ax.set_yticks(range(len(sel))); ax.set_yticklabels([f"{t['metric'].split('_')[0]}  {t['x']} / {t['y']} ({t['role'][:9]})" for t in sel], fontsize=7)
ax.set_xlabel("ratio of arm means (95% stratified bootstrap CI); dotted = planted M1 effect 1.061; * = p < 0.05/6")
ax.set_title("Effect sizes, all bits pooled (stratified by bit length)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{d}/fig_effects.svg"); fig.savefig(f"{d}/fig_effects.png", dpi=150)

# pair46-scale small-x add-on
s = res["pair46_smallx"]
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.hist(s["ctrl_models"], bins=25, alpha=0.5, label="200 random crater models (a c⁴, b c⁶)", color="#c0392b")
ax.hist(s["ctrl_random"], bins=25, alpha=0.5, label="200 random (a,b) over same p", color="#7f8c8d")
ax.axvline(s["crater"], color="#c0392b", lw=2, label=f"given crater model (rank {s['crater_rank_models']}/200)")
ax.axvline(s["floor"], color="k", lw=2, ls="--", label=f"given floor model (rank {s['floor_rank_random']}/200 vs random)")
ax.set_xlabel("M3 small-x yield at p46 (B = 2^16)"); ax.legend(fontsize=7)
ax.set_title("Pair-46 scale: small-x yield of the given models vs controls", fontsize=9)
fig.tight_layout(); fig.savefig(f"{d}/fig_pair46_smallx.svg"); fig.savefig(f"{d}/fig_pair46_smallx.png", dpi=150)
print("plots written")
