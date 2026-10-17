"""Summarise results.json (recorded per-unit values): log-ratio effect sizes
with bootstrap CIs, control checks, plots. Usage: python3 analyze.py"""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = json.load(open("results.json"))
U = R["units"]
ARMS = ["random2", "greedy", "learned_greedy", "learned_success", "known_log"]
rng = np.random.default_rng(0)

def boot(v):
    v = np.asarray(v); bs = rng.choice(v, (2000, len(v))).mean(1)
    return [float(v.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]

S = {"cells": [], "controls": {}, "identity_sum_count_eq_S2": all(
    m["sum_count_all_r"] == m["S_size"] ** 2 for u in U for m in u["arms"].values())}
for mu in sorted({u["mu"] for u in U}):
    for bits in sorted({u["bits"] for u in U}):
        us = [u for u in U if u["mu"] == mu and u["bits"] == bits]
        cell = {"mu": mu, "bits": bits, "n_units": len(us),
                "ceiling_ratio_bound": None, "arms": {}}
        L = np.mean([u["L"] for u in us]); N = np.mean([u["N"] for u in us])
        cell["ceiling_ratio_bound"] = float(np.mean([(u["L"] * (2 * u["L"] + 1) / u["N"]) / u["arms"]["random"]["held_success"] for u in us]))
        for a in ARMS + ["random"]:
            d = {"held_success": boot([u["arms"][a]["held_success"] for u in us]),
                 "train_success": boot([u["arms"][a]["train_success"] for u in us]),
                 "held_mean_count_ratio": boot([u["arms"][a]["held_mean_count"] / u["arms"]["random"]["held_mean_count"] for u in us])}
            if a != "random":
                d["logratio_vs_random"] = boot([math.log(u["arms"][a]["held_success"] / u["arms"]["random"]["held_success"]) for u in us])
            if a.startswith("learned"):
                d["logratio_vs_greedy"] = boot([math.log(u["arms"][a]["held_success"] / u["arms"]["greedy"]["held_success"]) for u in us])
            cell["arms"][a] = d
        cell["auc_learned_greedy"] = boot([u["auc_learned_greedy"] for u in us])
        cell["auc_learned_success"] = boot([u["auc_learned_success"] for u in us])
        cell["auc_pc2_mod7"] = boot([u["auc_pc2_mod7"] for u in us])
        S["cells"].append(cell)

C = S["cells"]
S["controls"] = {
    "PC1_known_log_upper_lt_0": all(c["arms"]["known_log"]["logratio_vs_random"][2] < 0 for c in C),
    "random2_covers_0": [ (c["mu"], c["bits"], c["arms"]["random2"]["logratio_vs_random"][1] <= 0 <= c["arms"]["random2"]["logratio_vs_random"][2]) for c in C],
    "PC2_auc_min": min(u["auc_pc2_mod7"] for u in U),
}
win = {}
for a in ["learned_greedy", "learned_success"]:
    win[a] = all(c["arms"][a]["logratio_vs_random"][1] > 0 and c["arms"][a]["logratio_vs_greedy"][1] > 0 for c in C)
S["decision_learned_advances"] = win
json.dump(S, open("summary.json", "w"), indent=1)

# Figure 1: held-out log-ratio vs random per arm, per bits, per mu
cols = {"random2": "#888888", "greedy": "#d95f02", "learned_greedy": "#1b9e77", "learned_success": "#7570b3"}
fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, mu in zip(axs, [0.25, 1.0]):
    cs = [c for c in C if c["mu"] == mu]
    for j, a in enumerate(cols):
        xs = [c["bits"] + (j - 1.5) * 0.25 for c in cs]
        m = [c["arms"][a]["logratio_vs_random"] for c in cs]
        ax.errorbar(xs, [v[0] for v in m], yerr=[[v[0] - v[1] for v in m], [v[2] - v[0] for v in m]], fmt="o", color=cols[a], label=a, capsize=3)
    ax.axhline(0, color="k", lw=0.6)
    cb = [math.log(c["ceiling_ratio_bound"]) for c in cs]
    ax.plot([c["bits"] for c in cs], cb, "k--", lw=0.8, label="exact ceiling log(L(2L+1)/N / H_rand)")
    ax.set_title(f"mu = 2L^2/N = {mu}"); ax.set_xlabel("bits of p"); ax.set_xticks([14, 16, 18])
axs[0].set_ylabel("held-out log(H_arm / H_random)  (mean, 95% bootstrap CI)")
axs[1].legend(fontsize=7, loc="upper left")
fig.suptitle("08-learned: held-out m=2 success vs matched random (12 units per point = 4 curves x 3 seeds)", fontsize=9)
fig.tight_layout(); fig.savefig("fig_heldout.svg"); fig.savefig("fig_heldout.png", dpi=150)

# Figure 2: train vs held-out success (overfit picture), mu=1.0 and 0.25 pooled over bits
fig, ax = plt.subplots(figsize=(6, 4))
for a, cc in list(cols.items()) + [("random", "#000000"), ("known_log", "#e7298a")]:
    tr = [u["arms"][a]["train_success"] for u in U]; he = [u["arms"][a]["held_success"] for u in U]
    ax.scatter(tr, he, s=10, color=cc, label=a, alpha=0.7)
ax.plot([0, 1], [0, 1], "k:", lw=0.7)
ax.set_xlabel("training-target success"); ax.set_ylabel("held-out success")
ax.set_title("in-sample vs held-out (all 72 units, both mu)", fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig("fig_overfit.svg"); fig.savefig("fig_overfit.png", dpi=150)
print(json.dumps({"controls": S["controls"], "decision": win, "identity": S["identity_sum_count_eq_S2"]}, indent=1))
for c in C:
    print(c["mu"], c["bits"], {a: [round(x, 4) for x in c["arms"][a].get("logratio_vs_random", [0])] for a in ARMS},
          "G:", {a: [round(x, 4) for x in c["arms"][a]["logratio_vs_greedy"]] for a in ["learned_greedy", "learned_success"]},
          "ceil", round(c["ceiling_ratio_bound"], 3), "auc", [round(c[k][0], 3) for k in ["auc_learned_greedy", "auc_learned_success", "auc_pc2_mod7"]],
          "train", {a: round(c["arms"][a]["train_success"][0], 3) for a in ["random", "greedy", "learned_greedy"]})
