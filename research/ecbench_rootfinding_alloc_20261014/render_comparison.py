#!/usr/bin/env python3
"""Render descriptive paired timing diagnostics from sealed comparison JSON."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

study = Path(__file__).parent
candidate = json.loads((study / "candidate-vs-baseline-a0-final.json").read_text())
drift = json.loads((study / "baseline-drift-final.json").read_text())
rows = [
    ("Batch 8, relations", "b48-batch-8", "relations_wall"),
    ("Batch 8, full solve", "b48-batch-8", "solve_wall"),
    ("Serial 8, relations", "b48-serial-8", "relations_wall"),
    ("Rho, full solve", "rho-neg", "solve_wall"),
    ("Rho A/A, full solve", "rho-neg-aa", "solve_wall"),
]

fig, ax = plt.subplots(figsize=(9.5, 4.1), layout="constrained")
for offset, (label, arm, metric) in enumerate(rows):
    y = len(rows) - 1 - offset
    for panel, color, dy, legend in (
        (candidate, "#2468a0", 0.14, "Candidate / first baseline"),
        (drift, "#c56828", -0.14, "Repeat baseline / first baseline"),
    ):
        entry = panel["l0_timing"][arm][metric]
        center = entry["candidate_over_baseline_geomean"]
        lo, hi = entry["target_cluster_bootstrap_95_percent"]
        ax.errorbar(
            center, y + dy, xerr=[[center - lo], [hi - center]], fmt="o",
            color=color, capsize=2.5, label=legend if offset == 0 else None,
        )
ax.axvline(1.0, color="#444444", linewidth=1, linestyle="--")
ax.set_yticks(range(len(rows)), [r[0] for r in reversed(rows)])
ax.set_xlabel("Paired wall-time ratio; target-cluster bootstrap interval")
ax.set_xlim(0.5, 1.8)
ax.set_title("Prime relation collector buffer reuse: L0 shared-host diagnostics")
ax.grid(axis="x", alpha=0.2)
ax.legend(loc="lower right", fontsize=8)
fig.savefig(study / "timing-ratios.svg")
fig.savefig(study / "timing-ratios.png", dpi=180)
plt.close(fig)
svg = study / "timing-ratios.svg"
svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
