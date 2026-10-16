#!/usr/bin/env python3
"""Render recorded p61 adjacency calibration metrics without rerunning work."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "metrics.json").read_text(encoding="utf-8"))
ROWS = DATA["calibrations"]


def main() -> None:
    labels = [
        "1 vertex" if row["vertices"] == 1 else f"{row['vertices']:,} vertices"
        for row in ROWS
    ]
    generation = [row["generation_elapsed_ms"] / 1000 for row in ROWS]
    replay = [row["verification_elapsed_ms"] / 1000 for row in ROWS]
    logical_per_vertex = [
        row["logical_uncompressed_bytes"] / row["vertices"] for row in ROWS
    ]
    stored_per_vertex = [row["stored_gzip_bytes"] / row["vertices"] for row in ROWS]
    x = np.arange(len(labels))
    width = 0.34

    plt.rcParams.update(
        {
            "font.size": 10.5,
            "axes.titlesize": 13,
            "axes.labelsize": 10.5,
            "figure.facecolor": "white",
            "axes.facecolor": "#f7f9fc",
        }
    )
    fig, (timing_axis, byte_axis) = plt.subplots(1, 2, figsize=(12.4, 5.4))
    fig.suptitle(
        "p61 selected-degree adjacency: replayed calibration windows",
        y=0.98,
        fontsize=15,
        weight="bold",
    )

    timing_axis.bar(x - width / 2, generation, width, label="Generate", color="#2367a8")
    timing_axis.bar(x + width / 2, replay, width, label="Independent replay", color="#b45309")
    timing_axis.set_yscale("log")
    timing_axis.set_ylabel("Elapsed seconds (log scale)")
    timing_axis.set_title("Startup-inclusive elapsed time")
    timing_axis.set_xticks(x, labels)
    timing_axis.grid(axis="y", which="both", color="#d7dde5", linewidth=0.8)
    timing_axis.set_axisbelow(True)
    timing_axis.legend(frameon=False, loc="upper left")
    for position, value in zip(x - width / 2, generation):
        timing_axis.text(position, value * 1.07, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    for position, value in zip(x + width / 2, replay):
        timing_axis.text(position, value * 1.07, f"{value:.3f}", ha="center", va="bottom", fontsize=8)

    byte_axis.bar(x - width / 2, logical_per_vertex, width, label="Logical JSONL", color="#15803d")
    byte_axis.bar(x + width / 2, stored_per_vertex, width, label="Stored gzip", color="#7c3aed")
    byte_axis.set_ylabel("Bytes per source vertex")
    byte_axis.set_title("Observed representation density")
    byte_axis.set_xticks(x, labels)
    byte_axis.grid(axis="y", color="#d7dde5", linewidth=0.8)
    byte_axis.set_axisbelow(True)
    byte_axis.legend(frameon=False, loc="upper right")
    for position, value in zip(x - width / 2, logical_per_vertex):
        byte_axis.text(position, value + 180, f"{value:,.0f}", ha="center", va="bottom", fontsize=8)
    for position, value in zip(x + width / 2, stored_per_vertex):
        byte_axis.text(position, value + 180, f"{value:,.0f}", ha="center", va="bottom", fontsize=8)

    fig.text(
        0.5,
        0.018,
        "Single observations on a shared host; every window passed exact independent replay. Timings are descriptive.",
        ha="center",
        fontsize=8.5,
        color="#475569",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.17, wspace=0.25)
    fig.savefig(ROOT / "adjacency-calibration.svg")
    fig.savefig(ROOT / "adjacency-calibration.png", dpi=160)


if __name__ == "__main__":
    main()
