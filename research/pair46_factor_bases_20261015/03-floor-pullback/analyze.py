#!/usr/bin/env python3
"""Plot recorded values from results.json + results_ext.json; write summary.json."""
import json, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

runs = []
for f, tag in (("results.json", "main"), ("results_ext.json", "ext E1")):
    d = json.load(open(f))
    for r in d["instances"]:
        r["_run"] = tag
        runs.append(r)
runs.sort(key=lambda r: r["v"])
ARMS = [("A1_pull_smallx", "A1 phi(crater small-x)", "#2a6fdb"),
        ("A2_pull_coset", "A2 phi(crater coset)", "#e07b00"),
        ("C1_floor_smallx", "C1 floor small-x", "#3b9c5a")]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
xs = list(range(len(runs)))
for k, (key, lab, col) in enumerate(ARMS):
    ax.plot([x + (k - 1) * 0.15 for x in xs], [r["stats"][key]["z"] for r in runs],
            "o", color=col, label=lab)
ax.axhspan(-3, 3, color="#bbbbbb", alpha=0.3, label="|z|<3 (null band)")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(xs); ax.set_xticklabels([f"v={r['v']}\n({r['_run']})" for r in runs], fontsize=7)
ax.set_ylabel("z vs random-x floor control (10 seeds)")
ax.set_title("2-term decomposition yield, z-score")
ax.set_ylim(-4, 4)
ax.legend(fontsize=7, loc="lower left")

ax = axes[1]
for k, (key, lab, col) in enumerate(ARMS):
    y = [r["stats"][key]["ratio"] for r in runs]
    lo = [r["stats"][key]["ratio"] - r["stats"][key]["ratio_ci95"][0] for r in runs]
    hi = [r["stats"][key]["ratio_ci95"][1] - r["stats"][key]["ratio"] for r in runs]
    ax.errorbar([x + (k - 1) * 0.15 for x in xs], y, yerr=[lo, hi], fmt="o", color=col,
                capsize=2, label=lab)
pp = [r["stats"]["Pplus_interval"]["ratio"] for r in runs]
ax.axhline(1, color="k", lw=0.6)
ax.set_xticks(xs); ax.set_xticklabels([f"v={r['v']}" for r in runs], fontsize=7)
ax.set_ylabel("Y(arm) / mean Y(random-x)  [95% CI]")
ax.set_title(f"yield ratio (positive control P+ ratio: {min(pp):.4f}-{max(pp):.4f}, off-scale)")
ax.title.set_fontsize(9)
ax.legend(fontsize=7)
fig.suptitle("Item 03 floor-pullback: toy pair-46 analogues p=5v^2-v+1, N prime, m=2*round(sqrt(N)/2)",
             fontsize=9)
fig.tight_layout()
fig.savefig("fig_yield.svg"); fig.savefig("fig_yield.png", dpi=150)

summary = {"instances": len(runs),
           "v": [r["v"] for r in runs],
           "X0_identity_all_exact": all(r["arms"][a]["X0_exact_equal"] for r in runs
                                        for a in ("A1_pull_smallx", "A2_pull_coset")),
           "max_abs_z": {k: max(abs(r["stats"][k]["z"]) for r in runs) for k, _, _ in ARMS},
           "mean_z": {k: sum(r["stats"][k]["z"] for r in runs) / len(runs) for k, _, _ in ARMS},
           "Pplus_min_abs_z": min(abs(r["stats"]["Pplus_interval"]["z"]) for r in runs),
           "ratio_ranges": {k: [min(r["stats"][k]["ratio_ci95"][0] for r in runs),
                                max(r["stats"][k]["ratio_ci95"][1] for r in runs)] for k, _, _ in ARMS},
           "xfrac_lt_p16": {k: [r["arms"][k]["xfrac_lt_p16"] for r in runs]
                            for k in ("A1_pull_smallx", "A2_pull_coset", "C1_floor_smallx")},
           "xfrac_C0_mean": [sum(c["xfrac_lt_p16"] for c in r["arms"]["C0_random_x"]) / 10 for r in runs],
           "horizontal_subgroups_skipped": [sum(t["horizontal"] for t in r["subgroups_tried"]) for r in runs],
           }
json.dump(summary, open("summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
