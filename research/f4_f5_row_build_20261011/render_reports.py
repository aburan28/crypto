#!/usr/bin/env python3
"""Render recorded native audit results only; no research or benchmark execution."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import typst

ROOT = Path(__file__).resolve().parent
for report in sorted((ROOT / "reports").glob("run-*")):
    data = json.loads((report / "graph-data.json").read_text())
    arms = {row["arm"]: row for row in data["arms"]}
    mode = data["graph_mode"]
    labels = ["F4 complete call", "F5 complete call", "F5 criterion"]
    keys = ["f4", "f5", "criterion"]
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 3.65), constrained_layout=True)
    for ax, field, unit in zip(
        axes,
        ["allocations", "allocated_bytes"],
        ["Allocation / reallocation requests", "Requested bytes (cumulative)"],
    ):
        baseline = [arms[key]["baseline_" + field] for key in keys]
        if mode == "baseline":
            bars = ax.bar(range(3), baseline, width=0.48, color="#44718d", label="Baseline")
            ax.bar_label(bars, labels=[f"{v:,}" for v in baseline], fontsize=8, padding=3)
        else:
            candidate = [arms[key]["candidate_" + field] for key in keys]
            left = ax.bar([x - 0.18 for x in range(3)], baseline, width=0.34, color="#44718d", label="Baseline")
            right = ax.bar([x + 0.18 for x in range(3)], candidate, width=0.34, color="#ca7854",
                           label="Baseline copy" if mode == "aa" else "Candidate")
            ax.bar_label(left, labels=[f"{v:,}" for v in baseline], fontsize=7, padding=3, rotation=90)
            ax.bar_label(right, labels=[f"{v:,}" for v in candidate], fontsize=7, padding=3, rotation=90)
        ax.set_xticks(range(3), labels, fontsize=9)
        ax.set_ylabel(unit, fontsize=9)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        ax.set_ylim(0, ax.get_ylim()[1] * 1.30)
        ax.grid(axis="y", alpha=0.2)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle(data["caption"], fontsize=11)
    fig.savefig(report / "allocation-audit.svg", metadata={"Description": data["caption"]})
    plt.close(fig)
    typst.compile(str(report / "report.typ"), output=str(report / "report.pdf"))
result = ROOT / "report.typ"
if result.exists():
    typst.compile(str(result), output=str(ROOT / "report.pdf"))
