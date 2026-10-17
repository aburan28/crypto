#!/usr/bin/env python3
"""Render the p61 selected-class evidence gates from recorded receipts."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent


def box(axis, x, y, width, height, title, detail, color, status):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.6,
        edgecolor=color,
        facecolor="#ffffff",
    )
    axis.add_patch(patch)
    axis.text(x + 0.025, y + height - 0.055, title, fontsize=11.5, weight="bold", va="top")
    axis.text(x + 0.025, y + height - 0.125, detail, fontsize=8.9, va="top", color="#334155", linespacing=1.3)
    axis.text(
        x + width,
        y + height + 0.012,
        status,
        fontsize=8.5,
        weight="bold",
        ha="right",
        va="bottom",
        color=color,
    )


def arrow(axis, start, end, color="#64748b"):
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.5,
            color=color,
            shrinkA=3,
            shrinkB=3,
        )
    )


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig, axis = plt.subplots(figsize=(13.5, 6.6))
    fig.patch.set_facecolor("white")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")

    fig.suptitle("p61 selected-class graph: evidence gates and current boundary", y=0.97, fontsize=17, weight="bold")
    axis.text(
        0.5,
        0.91,
        "Certified class order and a replayed prefix are inputs to - not substitutes for - complete graph coverage.",
        ha="center",
        fontsize=10.5,
        color="#475569",
    )

    green = "#15803d"
    amber = "#b45309"
    slate = "#64748b"
    blue = "#1d4ed8"

    box(
        axis,
        0.04,
        0.58,
        0.25,
        0.23,
        "1. Exact class order",
        "PARI flag-0 certificate = 1\n|Cl(D)| = 358,094,304\nCl(D) = C179047152 x C2",
        green,
        "CERTIFIED",
    )
    box(
        axis,
        0.375,
        0.58,
        0.25,
        0.23,
        "2. Immutable census prefix",
        "x = 0..16,704 across 18 lanes\n300,672 records in four segments\ncheckpoint + hash-chain bound",
        green,
        "GENERATED + REPLAYED",
    )
    box(
        axis,
        0.71,
        0.58,
        0.25,
        0.23,
        "3. External uniqueness",
        "300,672 unique j-invariants\n0 duplicates\ncomplete_class_coverage = false",
        amber,
        "PARTIAL CENSUS PASS",
    )

    arrow(axis, (0.29, 0.695), (0.375, 0.695))
    arrow(axis, (0.625, 0.695), (0.71, 0.695))

    box(
        axis,
        0.12,
        0.20,
        0.31,
        0.22,
        "4. Full selected-class census",
        "Generate and replay all 358,094,304 vertices\nthen externally sort and prove global uniqueness\ncurrent coverage: 0.0839645%",
        slate,
        "NOT COMPLETE",
    )
    box(
        axis,
        0.57,
        0.20,
        0.31,
        0.22,
        "5. Selected-degree adjacency",
        "Join every vertex and emit degree 5/13/37/59\nadjacency: 2,506,660,128 directed records\nthen attach frozen hardness labels for M7",
        slate,
        "SEPARATE PASS PENDING",
    )

    arrow(axis, (0.835, 0.58), (0.35, 0.42), amber)
    arrow(axis, (0.43, 0.31), (0.57, 0.31), blue)

    axis.text(
        0.5,
        0.075,
        "Scope: one selected ordinary isogeny class over p = 2^61 - 1. This is not an enumeration of every curve or every isogeny class over the field.",
        ha="center",
        fontsize=9.5,
        color="#334155",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f1f5f9", "edgecolor": "#cbd5e1"},
    )

    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.04)
    fig.savefig(ROOT / "evidence-pipeline.svg")
    fig.savefig(ROOT / "evidence-pipeline.png", dpi=160)


if __name__ == "__main__":
    main()
