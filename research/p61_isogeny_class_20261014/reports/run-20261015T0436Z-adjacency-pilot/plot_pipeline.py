#!/usr/bin/env python3
"""Render the p61 complete-graph evidence pipeline and current boundary."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent


def stage(axis, x, title, detail, color, status):
    width = 0.17
    height = 0.42
    patch = FancyBboxPatch(
        (x, 0.31),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=1.5,
        edgecolor=color,
        facecolor="#ffffff",
    )
    axis.add_patch(patch)
    axis.text(x + 0.015, 0.68, title, fontsize=10.5, weight="bold", va="top")
    axis.text(x + 0.015, 0.60, detail, fontsize=8.3, va="top", color="#334155", linespacing=1.35)
    axis.text(x + width / 2, 0.275, status, fontsize=8.2, weight="bold", ha="center", color=color)


def arrow(axis, start, end):
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.4,
            color="#64748b",
            shrinkA=3,
            shrinkB=3,
        )
    )


def main() -> None:
    fig, axis = plt.subplots(figsize=(13.2, 5.5))
    fig.patch.set_facecolor("white")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    fig.suptitle(
        "p61 selected-class graph: evidence pipeline and current boundary",
        y=0.96,
        fontsize=16,
        weight="bold",
    )
    axis.text(
        0.5,
        0.86,
        "Each downstream stage requires immutable completion evidence from the stage before it.",
        ha="center",
        fontsize=10,
        color="#475569",
    )

    stage(axis, 0.02, "1. Class order", "358,094,304 vertices\nPARI flag-0 certificate\nclass lattice frozen", "#15803d", "CERTIFIED")
    stage(axis, 0.22, "2. Vertex census", "2,659,968 records replayed\nx = 0..147,776\nglobal uniqueness: x = 0..33,088", "#b45309", "ACTIVE / PARTIAL")
    stage(axis, 0.42, "3. Adjacency", "4,096 production vertices\n28,672 edges replayed\nfull target: 2,506,660,128", "#b45309", "ACTIVE / PARTIAL")
    stage(axis, 0.62, "4. Dense join", "Resolve every target j\nreciprocity + exact counts\npartitioned immutable index", "#64748b", "PENDING")
    stage(axis, 0.82, "5. ML tables", "hardness contexts + masks\nfield/class-held-out splits\nnull controls + abstention", "#64748b", "PENDING")
    for x in (0.19, 0.39, 0.59, 0.79):
        arrow(axis, (x, 0.52), (x + 0.03, 0.52))

    axis.text(
        0.5,
        0.11,
        "Storage contract: 6.0 TB logical corpus cap, 1.5 TB free-space floor, create-only shards, distinct producer and verifier identities.",
        ha="center",
        fontsize=9.2,
        color="#334155",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f1f5f9", "edgecolor": "#cbd5e1"},
    )
    fig.subplots_adjust(left=0.015, right=0.985, top=0.89, bottom=0.02)
    fig.savefig(ROOT / "evidence-pipeline.svg")
    fig.savefig(ROOT / "evidence-pipeline.png", dpi=160)


if __name__ == "__main__":
    main()
