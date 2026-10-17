#!/usr/bin/env python3
"""Render recorded ECBench/Callgrind values; no measurement is done here."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "figure-data.json").read_text())

INK = "#17233c"
BLUE = "#1769aa"
TEAL = "#138a72"
AMBER = "#d08a00"
RED = "#b5443c"
GRID = "#d7dfeb"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.labelsize": 10,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(ROOT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(ROOT / f"{stem}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


callgrind = DATA["callgrind"]
rows = callgrind["iterations"]
values = [row["solve_ir"] for row in rows]
labels = [row["label"] for row in rows]
colors = [INK, BLUE, AMBER, TEAL]
threshold = callgrind["baseline_ir"] * callgrind["acceptance_ratio"]

fig, ax = plt.subplots(figsize=(9.0, 4.8))
bars = ax.bar(range(len(rows)), values, color=colors, width=0.68)
ax.axhline(threshold, color=RED, linestyle="--", linewidth=1.6,
           label="preregistered gate: 0.85× baseline")
ax.set_xticks(range(len(labels)), labels)
ax.set_ylabel("complete-solve instructions (Callgrind Ir)")
ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e9:.1f}B"))
ax.set_ylim(0, max(values) * 1.18)
ax.grid(axis="y", color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Strong-rho optimization converged below the frozen instruction gate", loc="left")
for bar, row in zip(bars, rows):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() - max(values) * 0.065,
        f"{row['ratio']:.3f}×",
        ha="center",
        va="top",
        fontweight="bold",
        color="white",
    )
ax.legend(frameon=False, loc="upper right")
fig.text(
    0.01,
    0.01,
    "n=61 public target W3a517768b166 · one deterministic profile per binary · "
    "Valgrind 3.19.0 · session ECBS1hffaa276d8a5c",
    fontsize=8.6,
    color="#4f6078",
)
fig.tight_layout(rect=(0, 0.055, 1, 1))
save(fig, "callgrind-progress")


controls = DATA["session_controls"]
fig, ax = plt.subplots(figsize=(9.0, 4.5))
y = list(range(len(controls)))
control_colors = [TEAL, AMBER, BLUE]
bars = ax.barh(y, [row["ratio"] for row in controls], color=control_colors, height=0.58)
bars[1].set_hatch("///")
bars[2].set_hatch("///")
ax.axvline(1.0, color=INK, linestyle="--", linewidth=1.4)
ax.set_yticks(y, [row["label"] for row in controls])
ax.invert_yaxis()
ax.set_xlim(0.9, 1.015)
ax.set_xlabel("ratio")
ax.grid(axis="x", color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Matched counts are identical; native wall samples stayed below L2", loc="left")
for bar, row in zip(bars, controls):
    ax.text(
        min(row["ratio"] + 0.002, 1.008),
        bar.get_y() + bar.get_height() / 2,
        f"{row['ratio']:.4f}×",
        va="center",
        ha="left",
        fontweight="bold",
    )
ax.text(0.999, -0.43, "parity = 1.0", ha="right", va="center", color=INK)
fig.text(
    0.01,
    0.01,
    "48 workloads × 5 measured rounds = 240 pairs · hatched wall rows are L0 descriptive, "
    "have no admitted interval, and are not performance claims · ECBS1hffaa276d8a5c",
    fontsize=8.6,
    color="#4f6078",
)
fig.tight_layout(rect=(0, 0.065, 1, 1))
save(fig, "session-controls")
