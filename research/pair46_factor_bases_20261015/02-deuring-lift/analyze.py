#!/usr/bin/env python3
"""Summary statistics and graphs from results.json (recorded values only)."""
import json, random, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

r = json.load(open("results.json"))
toy = [t for t in r["toy"] if "A2_rho" in t]
rng = random.Random(20261015)

def boot_diff_median(a, b, n=2000):
    ds = []
    for _ in range(n):
        aa = [rng.choice(a) for _ in a]; bb = [rng.choice(b) for _ in b]
        ds.append(st.median(aa) - st.median(bb))
    ds.sort(); return st.median(a) - st.median(b), ds[int(0.025 * n)], ds[int(0.975 * n)]

def boot_diff_mean(a, b, n=2000):
    ds = []
    for _ in range(n):
        aa = [rng.choice(a) for _ in a]; bb = [rng.choice(b) for _ in b]
        ds.append(st.mean(aa) - st.mean(bb))
    ds.sort(); return st.mean(a) - st.mean(b), ds[int(0.025 * n)], ds[int(0.975 * n)]

def mann_whitney(a, b):
    # two-sided, normal approximation with tie correction ignored (continuous data)
    allv = sorted([(x, 0) for x in a] + [(x, 1) for x in b])
    ranks = {}
    R1 = sum(i + 1 for i, (x, g) in enumerate(allv) if g == 0)
    n1, n2 = len(a), len(b)
    U = R1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2; sd = (n1 * n2 * (n1 + n2 + 1) / 12) ** 0.5
    z = (U - mu) / sd
    import math
    p = math.erfc(abs(z) / 2 ** 0.5)
    return U, z, p

A2 = [t["A2_rho"] for t in toy]; NC = [t["NC_rho"] for t in toy]; PC = [t["PC_rho"] for t in toy]
cA = [t["A2_cov"] for t in toy]; cN = [t["NC_cov"] for t in toy]
summ = {
    "n_primes": len(toy),
    "rho_median": {"A2": st.median(A2), "NC": st.median(NC), "PC": st.median(PC)},
    "rho_A2_minus_NC_median": boot_diff_median(A2, NC),
    "mann_whitney_A2_vs_NC": mann_whitney(A2, NC),
    "PC_all_below_0.5": all(x < 0.5 for x in PC),
    "PC_max": max(PC),
    "cov_A2_minus_NC_mean": boot_diff_mean(cA, cN),
    "lifted_relations_total": {
        "A2": sum(t["A2_lifted_relations"] for t in toy),
        "PC": sum(t["PC_lifted_relations"] for t in toy),
        "A3_xedni": sum(x["lifted_relations"] for t in r["toy"] for x in t["A3_xedni"]),
        "A3_trials": sum(len(t["A3_xedni"]) for t in r["toy"]),
        "A3_all_p_integral": all(x["p_integral"] and x["on_curve"] for t in r["toy"] for x in t["A3_xedni"]),
        "PC_xedni": sum(t["PC_xedni"]["lifted_relations"] for t in r["toy"]),
        "PC_xedni_trials": len(r["toy"]),
    },
    "relations_needed_per_prime": "|F| = 6 (plus 1 for the target)",
}
from collections import Counter
summ["twist_rank_census"] = dict(Counter(str(t["rank"]) for t in r["twists"]))
summ["twists_total"] = len(r["twists"])
json.dump(summ, open("summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))

# Figure 1: rho by arm and band
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
bands = [16, 20, 24]; arms = [("A2 canonical twists", "A2_rho", "#2a6fb0"), ("NC random points", "NC_rho", "#999999"), ("PC planted G,2G,3G", "PC_rho", "#c0392b")]
for k, (lab, key, col) in enumerate(arms):
    for i, b in enumerate(bands):
        ys = [t[key] for t in toy if t["bits"] == b]
        xs = [i + (k - 1) * 0.22 + rng.uniform(-0.05, 0.05) for _ in ys]
        ax.scatter(xs, ys, s=18, color=col, label=lab if i == 0 else None, alpha=0.85)
ax.axhline(1.0, color="k", lw=0.6, ls=":")
ax.set_xticks(range(3)); ax.set_xticklabels([f"{b}-bit p" for b in bands])
ax.set_ylabel(r"$\lambda_1(L)/\mathrm{GH}$ (mod-p relation lattice, R=6)")
ax.set_title("Shortest mod-p relation among factor-base points\n(8 primes per band; 1 = Gaussian heuristic)", fontsize=10)
ax.legend(fontsize=8, loc="upper right"); ax.set_ylim(0, 2.0)

ax = axes[1]
labels = ["A1 E(Q)\n(rank 1)", "A2 6 twists\n(rank 6)", "A3 Xedni\n(r=2..5)", "PC planted\nG,2G,3G", "PC Xedni\nP,-P"]
need = [1, 6, 14 / 4, 6, 3]
lt = summ["lifted_relations_total"]
got = [0, lt["A2"] / len(toy), lt["A3_xedni"] / lt["A3_trials"], lt["PC"] / len(toy), lt["PC_xedni"] / lt["PC_xedni_trials"]]
x = range(len(labels))
ax.bar([i - 0.2 for i in x], need, 0.4, color="#dddddd", edgecolor="k", label="relations needed (= |F|)")
ax.bar([i + 0.2 for i in x], got, 0.4, color=["#2a6fb0", "#2a6fb0", "#e08e0b", "#c0392b", "#c0392b"], label="lifted (Q-) relations found, mean/prime")
for i, g in enumerate(got):
    ax.text(i + 0.2, g + 0.1, f"{g:g}", ha="center", fontsize=8)
ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("relations per factor base")
ax.set_title("Characteristic-0 relations: |F| - rank(height Gram)\n(24 toy primes; exact rank test)", fontsize=10)
ax.legend(fontsize=8, loc="upper left"); ax.set_ylim(0, 8)
fig.suptitle("02-deuring-lift — canonical CM(-19) lift (361a1) of the pair-46 crater", fontsize=11)
fig.tight_layout()
fig.savefig("fig_relations.svg"); fig.savefig("fig_relations.png", dpi=150)

# Figure 2: twist rank census
fig, ax = plt.subplots(figsize=(8, 3.2))
tw = r["twists"]
for t in tw:
    c = {0: "#bbbbbb", 1: "#2a6fb0", 2: "#e08e0b", 3: "#c0392b"}.get(t["rank"], "k")
    ax.scatter(t["d"], t["rank"], color=c, s=14)
ax.set_xlabel("fundamental discriminant d of the twist E_d (one per {d, -19d} isogeny pair)")
ax.set_ylabel("rank E_d(Q)")
ax.set_yticks([0, 1, 2, 3])
ax.set_title(f"Mordell-Weil ranks of quadratic twists of 361a1, 1<|d|<=200 (n={len(tw)}, PARI ellrank, all bounds closed)", fontsize=9)
fig.tight_layout(); fig.savefig("fig_twist_ranks.svg"); fig.savefig("fig_twist_ranks.png", dpi=150)
