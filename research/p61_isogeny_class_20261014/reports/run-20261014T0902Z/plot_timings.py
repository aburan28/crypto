#!/usr/bin/env python3
"""Render recorded p61 pilot timings without recomputing experiment results."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "timings.json").read_text(encoding="utf-8"))
SEGMENTS = DATA["segments"]


def main() -> None:
    labels = [segment["label"] for segment in SEGMENTS]
    generation = [segment["generation_elapsed_ms"] / 1000 for segment in SEGMENTS]
    verification = [
        None
        if segment["verification_elapsed_ms"] is None
        else segment["verification_elapsed_ms"] / 1000
        for segment in SEGMENTS
    ]
    verification_plot = [value if value is not None else 1.0 for value in verification]
    parallel_verification = [
        segment["parallel_verification_elapsed_ms"] / 1000 for segment in SEGMENTS
    ]
    positions = range(len(labels))
    width = 0.24

    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 11,
        "figure.facecolor": "white",
        "axes.facecolor": "#f7f9fc",
    })
    fig, axis = plt.subplots(figsize=(10.2, 5.5))
    generation_bars = axis.bar(
        [position - width for position in positions],
        [value - 1.0 for value in generation],
        width,
        bottom=1.0,
        label="Generation",
        color="#2367a8",
    )
    verification_bars = axis.bar(
        list(positions),
        [max(value - 1.0, 0.0) for value in verification_plot],
        width,
        bottom=1.0,
        label="c9 replay",
        color="#d07a2d",
    )
    parallel_bars = axis.bar(
        [position + width for position in positions],
        [value - 1.0 for value in parallel_verification],
        width,
        bottom=1.0,
        label="d064 parallel replay",
        color="#9b3a4a",
    )
    fig.suptitle("p61 certified-prefix timings and parallel replay candidate", y=0.97, fontsize=14)
    axis.set_ylabel("Recorded elapsed time (seconds, log scale)")
    axis.set_xticks(
        list(positions),
        [f"{label}\n{segment['vertices']:,} vertices" for label, segment in zip(labels, SEGMENTS)],
    )
    axis.set_yscale("log")
    axis.set_ylim(
        1,
        max(generation + [value for value in verification if value is not None] + parallel_verification)
        * 1.45,
    )
    axis.grid(axis="y", which="both", color="#d7dde5", linewidth=0.8)
    axis.set_axisbelow(True)
    handles, legend_labels = axis.get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        frameon=False,
        fontsize=9.5,
    )
    for position, value in zip(positions, generation):
        axis.text(position - width, value * 1.06, f"{value:.3f} s", ha="center", va="bottom", fontsize=8, rotation=90)
    for position, value in zip(positions, verification):
        if value is None:
            axis.text(position, 1.18, "not run", ha="center", va="bottom", fontsize=8, rotation=90)
        else:
            axis.text(position, value * 1.06, f"{value:.3f} s", ha="center", va="bottom", fontsize=8, rotation=90)
    for position, value in zip(positions, parallel_verification):
        axis.text(position + width, value * 1.06, f"{value:.3f} s", ha="center", va="bottom", fontsize=8, rotation=90)
    fig.text(
        0.5,
        0.025,
        "Single observations on a contended host. All replays verify; timing is diagnostic.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#3f4b5b",
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.79, bottom=0.22)
    fig.savefig(ROOT / "pilot-timings.svg")
    fig.savefig(ROOT / "pilot-timings.png", dpi=150)


if __name__ == "__main__":
    main()
