#!/usr/bin/env python3
"""Render an immutable ecbench session's measured cold GAE decomposition."""

import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


session = Path(sys.argv[1])
out = Path(sys.argv[2])
arms = sys.argv[3:]
out.mkdir(parents=True, exist_ok=True)
by_arm = defaultdict(list)
for line in (session / "records.jsonl").read_text().splitlines():
    record = json.loads(line)
    if not record["warmup"]:
        by_arm[record["arm"]].append(record)


def mean_phase(records, phase):
    return statistics.mean(
        next((entry["gae"] for entry in row["phases"] if entry["name"] == phase), 0)
        for row in records
    )


summary = {}
for arm in arms:
    rows = by_arm[arm]
    summary[arm] = {
        "measured": len(rows),
        "verified": sum(row["outcome"]["status"] == "verified" for row in rows),
        "cold_gae_mean": statistics.mean(row["cost"]["total_gae"] for row in rows),
        "setup_gae_mean": mean_phase(rows, "oracle_setup"),
        "relations_gae_mean": mean_phase(rows, "relations"),
        "verify_gae_mean": mean_phase(rows, "verify"),
        "trials_mean": statistics.mean(row["counters"].get("targets_tried", 0) for row in rows),
        "relations_mean": statistics.mean(
            row["counters"].get("relations_found", 0) for row in rows
        ),
    }
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

fig, ax = plt.subplots(figsize=(8.0, 4.6), layout="constrained")
x = list(range(len(arms)))
bottom = [0.0] * len(arms)
for name, color in [
    ("setup_gae_mean", "#4c78a8"),
    ("relations_gae_mean", "#f58518"),
    ("verify_gae_mean", "#54a24b"),
]:
    values = [summary[arm][name] for arm in arms]
    ax.bar(x, values, bottom=bottom, color=color, label=name.removesuffix("_gae_mean"))
    bottom = [a + b for a, b in zip(bottom, values)]
for i, arm in enumerate(arms):
    ax.text(i, bottom[i] * 1.02, f"{summary[arm]['cold_gae_mean']:,.0f}", ha="center", fontsize=9)
ax.set_xticks(x, arms)
ax.set_ylim(0, max(bottom) * 1.13)
ax.set_ylabel("Mean cold counted group-addition equivalents")
ax.set_title("One public target per execution; setup charged each time")
ax.legend(ncol=3, loc="upper center")
ax.grid(axis="y", alpha=0.2)
fig.savefig(out / "counts.svg")
fig.savefig(out / "counts.png", dpi=150)
