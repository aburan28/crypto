#!/usr/bin/env python3
"""Render recorded ECBS1hb77d7527b579 values without recomputing the run."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
with (HERE / "data.json").open(encoding="utf-8") as handle:
    data = json.load(handle)

arms = data["arms"]
labels = [arm["label"] for arm in arms]
means = np.array([arm["mean_solve_seconds"] for arm in arms])
mins = np.array([arm["min_solve_seconds"] for arm in arms])
maxs = np.array([arm["max_solve_seconds"] for arm in arms])

fig, (ax_time, ax_phase) = plt.subplots(
    1, 2, figsize=(11.2, 4.8), gridspec_kw={"width_ratios": [1.08, 1.35]}
)

colors = ["#355c7d", "#2a9d8f", "#e76f51"]
x = np.arange(len(arms))
ax_time.bar(x, means, color=colors, width=0.63)
ax_time.errorbar(
    x,
    means,
    yerr=np.vstack((means - mins, maxs - means)),
    fmt="none",
    ecolor="#1f2933",
    capsize=4,
    linewidth=1.2,
)
ax_time.set_yscale("log")
ax_time.set_ylabel("cold solve wall time (seconds, log scale)")
ax_time.set_xticks(x, labels, rotation=14, ha="right")
ax_time.set_title("Mean of five measured solves; whiskers are min–max")
ax_time.grid(axis="y", which="both", alpha=0.25)
for xpos, mean in zip(x, means):
    shown = f"{mean * 1000:.1f} ms" if mean < 1 else f"{mean:.3f} s"
    ax_time.text(xpos, mean * 1.15, shown, ha="center", va="bottom", fontsize=9)

ic_arms = [next(arm for arm in arms if arm["arm"] == key) for key in ("ic-b3-chain", "ic-b2-quartic")]
phase_names = ["relations", "oracle setup", "factor base", "linear algebra", "verification"]
phase_colors = {
    "relations": "#d1495b",
    "oracle setup": "#edae49",
    "factor base": "#00798c",
    "linear algebra": "#30638e",
    "verification": "#8d99ae",
}
left = np.zeros(len(ic_arms))
for phase in phase_names:
    shares = np.array(
        [arm["phase_mean_seconds"][phase] / arm["mean_solve_seconds"] * 100 for arm in ic_arms]
    )
    ax_phase.barh(
        [arm["label"] for arm in ic_arms],
        shares,
        left=left,
        color=phase_colors[phase],
        label=phase,
    )
    left += shares
ax_phase.set_xlim(0, 100)
ax_phase.set_xlabel("share of cold solve wall time (%)")
ax_phase.set_title("Relation extraction consumes the cold solve")
ax_phase.grid(axis="x", alpha=0.25)
for ypos, arm in enumerate(ic_arms):
    share = arm["relation_share_of_solve"] * 100
    ax_phase.text(50, ypos, f"relations {share:.2f}%", ha="center", va="center", color="white", fontweight="bold")
handles, legend_labels = ax_phase.get_legend_handles_labels()
fig.legend(
    handles,
    legend_labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.055),
    ncol=5,
    frameon=False,
    fontsize=8,
)

fig.suptitle(
    "n=41 Koblitz end-to-end timing — ECBS1hb77d7527b579 (L0, descriptive)",
    fontsize=13,
    fontweight="bold",
)
fig.text(
    0.5,
    0.012,
    "All 20 measured executions verified; wall results are descriptive because PSI telemetry was unavailable.",
    ha="center",
    fontsize=9,
)
fig.tight_layout(rect=(0, 0.14, 1, 0.93))
fig.savefig(HERE / "e2e-timing.svg")
fig.savefig(HERE / "e2e-timing.png", dpi=150)
