# Layout only: plots recorded values from results.json. Computes nothing new beyond means of recorded rows.
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = json.load(open("results.json"))
S = {(s["bits"], s["quantity"]): s for s in d["summary"]}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
bits = [12, 14, 16]
arms = [("SX", "SX_coverage", "#333333"), ("RND (null)", "RND_minus_SX", "#888888"), ("KL (pos. ctrl)", "KL_minus_SX", "#2ca02c")] + \
       [(f"PRE_{m}", f"PRE_{m}_minus_SX", c) for m, c in zip([2, 3, 5, 7], ["#c6dbef", "#6baed6", "#2171b5", "#08306b"])]
w = 0.11
for i, (lab, q, col) in enumerate(arms):
    ys, lo, hi = [], [], []
    for b in bits:
        s = S[(b, q)]; base = 0 if q == "SX_coverage" else S[(b, "SX_coverage")]["mean"]
        ys.append(s["mean"] + base); lo.append(s["mean"] - s["ci95"][0]); hi.append(s["ci95"][1] - s["mean"])
    ax[0].bar([j + (i - 3) * w for j in range(3)], ys, w, yerr=[lo, hi], capsize=2, color=col, label=lab)
ax[0].set_xticks(range(3)); ax[0].set_xticklabels([f"{b}-bit" for b in bits])
ax[0].set_ylabel("coverage C = |{±P1±P2}|/N  (mean of 5 curves, 95% bootstrap CI)")
ax[0].set_title("M1: 2-term coverage, |V| = ceil(sqrt(N)/2)\nPRE_m = [m]^-1(SX) equals SX exactly (E4)", fontsize=9)
ax[0].legend(fontsize=7, ncol=2)
# rho panel
labs, vals, cols = ["SX/RND (list)"], [1.0], ["#333333"]
for m in [2, 3, 5, 7]: labs.append(f"PRE_{m}"); vals.append(1 / m**2); cols.append("#2171b5")
for m in [3, 5, 7]:
    labs.append(f"LVL_{m}"); vals.append(sum(S[(b, f"LVL_{m}_rho_mean")]["mean"] for b in bits) / 3); cols.append("#d62728")
for l in [3, 5, 7]:
    rows = [t for t in d["torsion"] if t["l"] == l]
    labs.append(f"TOR_{l}"); vals.append(sum(t["roots_on_E"] / t["psi_deg"] for t in rows) / len(rows)); cols.append("#ff7f0e")
ax[1].bar(range(len(vals)), vals, color=cols); ax[1].set_yscale("log")
ax[1].set_xticks(range(len(vals))); ax[1].set_xticklabels(labs, rotation=60, fontsize=7)
ax[1].axhline(1, color="gray", ls=":", lw=1)
ax[1].set_ylabel("rho = |V| / deg(membership polynomial)")
ax[1].set_title("M2: size per unit membership degree (all <= 1)\nLVL: 50 c x 15 curves; TOR: 9 curves per l", fontsize=9)
# pair46 panel
for k, c in enumerate(d["pair46"]):
    ms = [r["m"] for r in c["psi_roots"]]; onE = [r["on_E"] for r in c["psi_roots"]]; onT = [r["on_twist"] for r in c["psi_roots"]]
    xs = [j + (k - 0.5) * 0.35 for j in range(len(ms))]
    ax[2].bar(xs, onT, 0.35, color=["#9467bd", "#c5b0d5"][k], label=f"{c['curve']}: roots on twist")
    ax[2].bar(xs, onE, 0.35, bottom=onT, color="#d62728", label=f"{c['curve']}: roots on E" if k == 0 else None)
ax[2].set_xticks(range(len(ms))); ax[2].set_xticklabels([f"psi_{m}" if m > 2 else "f (m=2)" for m in ms])
ax[2].set_ylabel("number of F_p-roots (exact, gcd with x^p - x)")
ax[2].set_title("Pair 46: F_p-roots of division polynomials\n0 on E for every m (N prime); psi_5: 2 twist roots", fontsize=9)
ax[2].legend(fontsize=7)
fig.suptitle("07-torsion-division (results.json, protocol sha256 97f57630...)", fontsize=10)
fig.tight_layout()
fig.savefig("fig_torsion_division.svg"); fig.savefig("fig_torsion_division.png", dpi=150)
