# Layout only: plots recorded values from results.json / posthoc_v2_randomwalk.json. Computes nothing new.
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
r = json.load(open("results.json")); ph = json.load(open("posthoc_v2_randomwalk.json"))
bits = [t["bits"] for t in r["toys"]]
arms = ["rho", "rho-neg", "bsgs-pm", "abscissa-k1", "endo-chain-k4", "endo-chain-k16", "endo-chain-k64",
        "endo-chain-kF", "endo-chain-known", "endo-chain-k16-nofold"]
cols = {"rho": "#555555", "rho-neg": "#999999", "bsgs-pm": "#000000", "abscissa-k1": "#d62728",
        "endo-chain-k4": "#9ecae1", "endo-chain-k16": "#4292c6", "endo-chain-k64": "#2171b5",
        "endo-chain-kF": "#08306b", "endo-chain-known": "#2ca02c", "endo-chain-k16-nofold": "#ff7f0e"}
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for a in arms:
    m = [next(x for x in t["arms"] if x["arm"] == a) for t in r["toys"]]
    y = [x["S_mean"] for x in m]; lo = [x["S_mean"] - x["S_ci95"][0] for x in m]; hi = [x["S_ci95"][1] - x["S_mean"] for x in m]
    ax[0].errorbar(bits, y, yerr=[lo, hi], marker="o", capsize=3, label=a, color=cols[a])
for v, lab, ls in [(1.2533, "theory: rho walk 1.253", ":"), (0.886, "theory: rho-neg walk 0.886", "--"), (1.0, "theory: +-BSGS 1.0", "-.")]:
    ax[0].axhline(v, ls=ls, lw=1, color="gray", label=lab)
ax[0].set_yscale("log"); ax[0].set_xlabel("log2 p (toy CM(-19) curves, p = 5t^2+19t+19)")
ax[0].set_ylabel("S = group-op equivalents / sqrt(N)  (mean, bootstrap 95% CI, n=40)")
ax[0].set_title("Phase 1 (preregistered): full DLP solves, fixed-b walk"); ax[0].legend(fontsize=7, ncol=2); ax[0].set_xlim(15.5, 22.5)
for t in ph["toys"]:
    ks = []; ys = []; e1 = []; e2 = []
    for x in t["arms"]:
        k = x["k"]; ks.append(k); ys.append(x["S_mean"]); e1.append(x["S_mean"] - x["S_ci95"][0]); e2.append(x["S_ci95"][1] - x["S_mean"])
    ax[1].errorbar(ks, ys, yerr=[e1, e2], marker="s", capsize=3, label=f"{t['bits']}-bit (F={t['F']})")
ax[1].set_xscale("log", base=2); ax[1].set_xlabel("chain length k (k=1: abscissa/small-x; k=F: single chain)")
ax[1].set_ylabel("S (mean, 95% CI, n=40)"); ax[1].set_title("POST-HOC: random-relation walk, F = ceil(sqrt N) fixed")
ax[1].legend(fontsize=8)
fig.suptitle("01-endo-chain: degree-5 endomorphism chains vs rho / BSGS / abscissa (results.json, posthoc_v2_randomwalk.json)", fontsize=10)
fig.tight_layout(); fig.savefig("fig_cost.svg"); fig.savefig("fig_cost.png", dpi=150)

p2 = r["phase2_m2_rate"]
fig, ax = plt.subplots(figsize=(8, 4))
names = [x["fb"] for x in p2["rows"]]; rates = [x["rate"] for x in p2["rows"]]
err = [[x["rate"] - x["wilson95"][0] for x in p2["rows"]], [x["wilson95"][1] - x["rate"] for x in p2["rows"]]]
ax.bar(range(len(names)), rates, yerr=err, capsize=4, color=["#999999", "#d62728", "#4292c6", "#08306b", "#ff7f0e"])
ax.axhline(p2["rows"][0]["random_set_heuristic"], ls=":", color="k", lw=1); ax.text(3.6, p2["rows"][0]["random_set_heuristic"] + 0.004, "random-set heuristic", fontsize=7)
ax.set_xticks(range(len(names))); ax.set_xticklabels([n.replace(" (positive control)", "\n(positive ctrl)") for n in names], fontsize=8)
ax.set_ylabel("P(R = +-Pi +-Pj), Wilson 95% CI, T=3000")
ax.set_title(f"Phase 2: m=2 decomposition rate, 20-bit toy N={p2['N']}, |FB|={p2['F2']}", fontsize=10)
fig.tight_layout(); fig.savefig("fig_m2_rate.svg"); fig.savefig("fig_m2_rate.png", dpi=150)
print("ok")
