"""Lay out a recorded `volcano-rho --instance-study` run (instance_study.json,
instance_trials.csv): no computation beyond arranging recorded values
(AGENTS.md "Implementation language").

usage: python plot_instance.py RUN_DIR
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

run = Path(sys.argv[1])
R = json.loads((run / "instance_study.json").read_text())
rid = run.name
rows = defaultdict(list)
for r in csv.DictReader(open(run / "instance_trials.csv")):
    rows[r["arm"]].append(r)
arms = [a["arm"] for a in R["arms"]]
label = {
    "walk": "walk seed varies\n(P, k fixed)",
    "k": "k varies\n(P, walk fixed)",
    "generator": "generator varies\n(k, walk fixed)",
    "all": "all three vary",
    "small_k": "k in [1, 2^16]\n(P, walk fixed)",
}
colour = {"walk": "#2a6fdb", "k": "#d1495b", "generator": "#2e933c", "all": "#7b2cbf", "small_k": "#e07a00"}


def save(fig, name):
    fig.savefig(run / f"{name}.svg", bbox_inches="tight")
    fig.savefig(run / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


# 1. Every solve per arm, with mean ± 95% CI and the reference.
fig, ax = plt.subplots(figsize=(8, 4))
for i, a in enumerate(R["arms"]):
    xs = [int(r["iters"]) for r in rows[a["arm"]]]
    ax.scatter([i + 0.25 * ((j * 0.618) % 1 - 0.5) for j in range(len(xs))], xs, s=3, alpha=0.25,
               color=colour[a["arm"]], linewidths=0)
    s = a["iters"]
    ax.errorbar(i, s["mean"], yerr=s["ci95_half"], fmt="o", color="black", capsize=4, zorder=3)
    note = f"{s['mean']:,.0f}"
    if "vs_walk" in a:
        v = a["vs_walk"]
        note += f"\nz = {v['welch_z']:+.2f}, KS p = {v['ks_p']:.2f}"
    ax.annotate(note, (i, s["mean"]), xytext=(8, 4), textcoords="offset points", fontsize=7)
ref = R["rho"]["expected_iters_r32"]
ax.axhline(ref, color="grey", ls="--", lw=1, label=f"sqrt(pi n/2)/sqrt(1-1/32) = {ref:,.0f}")
ax.set_xticks(range(len(arms)), [label[a] for a in arms], fontsize=7)
ax.set_ylabel("rho iterations to solve")
ax.set_yscale("log")
ax.set_title(f"Does rho's cost depend on the generator or on k?  {R['trials_per_arm']} solves per arm, "
             f"curve C (n = {R['curve']['n']:,})\n{rid}; every solve returned its own k", fontsize=8)
ax.legend(fontsize=7, loc="lower right")
save(fig, "instance_arms")

# 2. Cumulative distributions of every arm.
fig, ax = plt.subplots(figsize=(7, 3.8))
for a in arms:
    xs = sorted(int(r["iters"]) for r in rows[a])
    ax.step(xs, [(i + 1) / len(xs) for i in range(len(xs))], where="post", color=colour[a], lw=1.2,
            label=label[a].replace("\n", " "))
ax.set_xlabel("rho iterations to solve")
ax.set_ylabel("fraction of solves finished")
ax.set_title(f"Cumulative distribution per arm ({rid})", fontsize=8)
ax.legend(fontsize=7, loc="lower right")
save(fig, "instance_cdf")

# 3. Cost against k itself: the k and small-k arms.
fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
for ax, a, scale in [(axs[0], "k", R["curve"]["n"]), (axs[1], "small_k", 1)]:
    ks = [int(r["k"]) / scale for r in rows[a]]
    it = [int(r["iters"]) for r in rows[a]]
    ax.scatter(ks, it, s=4, alpha=0.4, color=colour[a], linewidths=0)
    ax.axhline(ref, color="grey", ls="--", lw=1)
    ax.set_xlabel("k / n" if scale > 1 else "k")
    ax.set_title(label[a].replace("\n", " "), fontsize=8)
axs[0].set_ylabel("rho iterations to solve")
axs[0].set_yscale("log")
fig.suptitle(f"Iterations against the secret k ({rid})", fontsize=8)
save(fig, "instance_vs_k")
print("wrote", run)
