#!/usr/bin/env python3
"""Lay out a run's summary.json as figures (no computation: every plotted
number is a median, minimum or ratio already in summary.json, written by
examples/generic_dlp_analyze.rs).

    plot.py <run_dir>

Writes fig_bsgs.svg, fig_rho.svg, fig_external.svg and a 150 dpi PNG of each
into <run_dir>.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Reference palette (dataviz skill, references/palette.md), light mode.
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
ARM_A, ARM_B, PARI, BAND = "#2a78d6", "#eb6834", "#1baf7a", "#c9c8c2"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "font.size": 9, "axes.titlesize": 10, "legend.frameon": False, "lines.linewidth": 2,
    "lines.markersize": 6, "svg.fonttype": "none",
})

run = Path(sys.argv[1])
s = json.loads((run / "summary.json").read_text())
rid = s["run"]


def rows(part, kernel=None):
    out = [c for c in s["configs"] if c["part"] == part and (kernel is None or c["kernel"] == kernel)]
    return sorted(out, key=lambda c: c["bits"])


def save(fig, name):
    fig.savefig(run / f"{name}.svg")
    fig.savefig(run / f"{name}.png", dpi=150)
    plt.close(fig)


def ratio_panel(ax, series, timer, title):
    """Ratio of medians A/B with the per-round [min, max] as a bar, and the
    A/A2 per-round [min, max] as a grey band per size."""
    for i, (label, cs, color, marker) in enumerate(series):
        x = [c["bits"] + (i - (len(series) - 1) / 2) * 0.5 for c in cs]
        t = [c["timers"][timer] for c in cs]
        y = [r["ratio_of_medians_A_over_B"] for r in t]
        lo = [r["ratio_of_medians_A_over_B"] - r["per_round_A_over_B"]["min"] for r in t]
        hi = [r["per_round_A_over_B"]["max"] - r["ratio_of_medians_A_over_B"] for r in t]
        ax.errorbar(x, y, yerr=[lo, hi], color=color, marker=marker, capsize=3, label=f"{label}: A/B")
        for xi, r in zip(x, t):
            aa = r["per_round_A_over_A2"]
            ax.fill_between([xi - 0.25, xi + 0.25], aa["min"], aa["max"], color=BAND, alpha=0.7, lw=0)
    ax.set_xticks(sorted({c["bits"] for _, cs, _, _ in series for c in cs}))
    ax.axhline(1.0, color=INK2, lw=0.8)
    ax.fill_between([], [], [], color=BAND, label="A/A′ per-round range (noise)")
    ax.set_xlabel("group order (bits)")
    ax.set_ylabel(f"baseline / candidate, {timer}")
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)


# ── BSGS ────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
for kernel, ls in (("chain", "--"), ("paired", "-")):
    cs = rows("bsgs1", kernel)
    bits = [c["bits"] for c in cs]
    for arm, color, name in (("A", ARM_A, "baseline A"), ("B", ARM_B, "candidate B")):
        ax1.plot(bits, [c["timers"]["total_s"][arm]["median"] for c in cs], ls=ls, color=color,
                 marker="o", label=f"{name}, {kernel}")
ax1.set_yscale("log")
ax1.set_xlabel("group order (bits)")
ax1.set_ylabel("one-target solve, total_s (s, median)")
ax1.set_title("bsgs_fast, 1 thread, full group, negation map")
ax1.set_xticks(bits)
ax1.legend(fontsize=8)
ratio_panel(ax2, [("chain", rows("bsgs1", "chain"), ARM_A, "s"),
                  ("paired", rows("bsgs1", "paired"), ARM_B, "o")],
            "total_s", "speed ratio with per-round range")
fig.suptitle(f"BSGS table prefetch: {rid}; 7 rounds (5 at 52 bits), identical counted outputs", fontsize=9)
fig.tight_layout()
save(fig, "fig_bsgs")

# ── rho ─────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
cs = rows("rho1")
bits = [c["bits"] for c in cs]
for arm, color, name in (("A", ARM_A, "baseline A"), ("B", ARM_B, "candidate B")):
    ax1.plot(bits, [c["timers"]["solve_s"][arm]["median"] for c in cs], color=color, marker="o",
             label=f"{name}: whole solve")
    ax1.plot(bits, [c["timers"]["setup_s"][arm]["median"] for c in cs], color=color, marker="^",
             ls=":", label=f"{name}: start-up combine calls alone")
ax1.set_yscale("log")
ax1.set_xlabel("group order (bits)")
ax1.set_ylabel("seconds (median)")
ax1.set_title("prime_fast::rho_parallel, 1 thread, negation map")
ax1.set_xticks(bits)
ax1.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
ratio_panel(ax2, [("1 thread", rows("rho1"), ARM_A, "o"), ("3 threads", rows("rho3"), ARM_B, "s")],
            "solve_s", "speed ratio with per-round range")
fig.suptitle(f"rho start-up in Jacobian coordinates: {rid}; 7 rounds (5 at 52 bits)", fontsize=9)
fig.tight_layout()
save(fig, "fig_rho")

# ── external reference ──────────────────────────────────────────────────
ext = s.get("external", [])
if ext:
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for ax, role, title in ((axs[0], "bsgs", "BSGS curves: bsgs_fast paired vs PARI/GP elllog"),
                            (axs[1], "rho", "rho curves: rho_parallel vs PARI/GP elllog")):
        es = sorted((e for e in ext if e["role"] == role), key=lambda e: e["bits"])
        if not es:
            continue
        bits = [e["bits"] for e in es]
        ax.plot(bits, [e["pari_elllog_s"]["median"] for e in es], color=PARI, marker="D",
                label=f"{es[0]['tool']} elllog (median)")
        ax.plot(bits, [e["native_A_median_s"] for e in es], color=ARM_A, marker="o", label="baseline A (median)")
        ax.plot(bits, [e["native_B_median_s"] for e in es], color=ARM_B, marker="o", label="candidate B (median)")
        ax.set_yscale("log")
        ax.set_xticks(bits)
        ax.set_xlabel("group order (bits)")
        ax.set_ylabel("one-target solve (s), single thread")
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=8)
    fig.suptitle(f"External reference (AGENTS.md §15), same host, same curves and targets: {rid}", fontsize=9)
    fig.tight_layout()
    save(fig, "fig_external")
print("figures written to", run)
