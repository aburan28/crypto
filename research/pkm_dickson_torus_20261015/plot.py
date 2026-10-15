"""Lay out the recorded D and width values of cells.tsv and recorded_split.tsv.

No computation beyond grouping recorded rows (AGENTS.md "no Python": this
only plots). Writes d_vs_n.svg/.png and width_vs_n.svg/.png.
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

here = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(here / "cells.tsv"), delimiter="\t"))
rec = list(csv.DictReader(open(here / "recorded_split.tsv"), delimiter="\t"))

T24, T31, P1, P2, P0 = "14680063", "2046820351", "2013265921", "3221225473", "786433"
ARMS = [
    # label, kind, p, control, colour, marker, x-offset
    ("torus T31 (61·2^25−1)", "dickson-torus", T31, "tower", "#1f6fb4", "o", -0.24),
    ("torus T24 (7·2^21−1)", "dickson-torus", T24, "tower", "#2ca0d6", "s", -0.12),
    ("split Dickson p1 (this build)", "dickson", P1, "tower", "#d1495b", "^", 0.0),
    ("torus null T31", "dickson-torus", T31, "null", "#7d7d7d", "x", 0.12),
    ("torus null T24", "dickson-torus", T24, "null", "#a9a9a9", "+", 0.12),
    ("deep Kummer p2 (v2=30)", "kummer", P2, "tower", "#2e8b57", "D", 0.24),
    ("isogeny p0 (positive control)", "isogeny", P0, "tower", "#9467bd", "v", 0.36),
]


def finished(r):
    return r["timed_out"] == "false" and r["oversize"] == "false"


def series(kind, p, control, m, key):
    """Per N: the set of recorded values over finished targets."""
    by_n = defaultdict(set)
    for r in rows:
        if (r["kind"], r["p"], r["control"], r["m"]) == (kind, p, control, str(m)) and finished(r):
            by_n[int(r["N"])].add(int(r[key]))
    return dict(sorted(by_n.items()))


def recorded(m):
    by_n = {}
    for r in rec:
        if r["m"] == str(m) and r["p"] in (P0, P1):
            by_n[int(r["N"])] = int(r["D"])
    return dict(sorted(by_n.items()))


def d_plot(fmt):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), gridspec_kw={"width_ratios": [2.2, 1, 0.8]})
    for ax, m in zip(axes, (2, 3, 4)):
        rs = recorded(m)
        ax.step(list(rs), list(rs.values()), where="mid", color="#000000", lw=1.0, alpha=0.35,
                label="split, recorded (Kummer p1/p0, round 2)")
        for label, kind, p, control, colour, marker, dx in ARMS:
            s = series(kind, p, control, m, "D")
            if not s:
                continue
            xs, ys = [], []
            for n, vals in s.items():
                for v in vals:
                    xs.append(n + dx)
                    ys.append(v)
            ax.plot(xs, ys, marker=marker, ls="none", color=colour, ms=6, label=label)
        ax.set_title(f"m = {m}" + ("  (N = 15: one target per prime)" if m == 3 else ""), fontsize=9)
        ax.set_xticks({2: list(range(2, 21, 2)), 3: [6, 9, 12, 15], 4: [8, 12]}[m])
        ax.set_xlabel("N = m·t (tower variables)")
        ax.grid(alpha=0.3)
        ax.set_ylim(2.5, 7.5)
    axes[0].set_ylabel("F4 solving degree D (f4_fp_tower)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.suptitle("PKM towers: D vs N, torus vs split vs null, and deep-2-adic Kummer "
                 "(research/pkm_dickson_torus_20261015, runs R0–R15; 2 random targets per cell unless noted)",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0.12, 1, 0.95))
    out = here / f"d_vs_n.{fmt}"
    fig.savefig(out, dpi=150 if fmt == "png" else None)
    plt.close(fig)


def width_plot(fmt):
    fig, ax = plt.subplots(figsize=(7, 4.4))
    for m, ls in ((2, "-"), (3, "--"), (4, ":")):
        for label, kind, p, control, colour, marker, dx in ARMS[:3] + ARMS[5:6]:
            s = series(kind, p, control, m, "width_to_solution")
            if not s:
                continue
            ax.plot(list(s), [max(v) for v in s.values()], ls=ls, marker=marker, color=colour,
                    ms=4, label=f"{label}, m={m}")
    ax.set_yscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel("widest F4 matrix up to the solution (columns)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=6, ncol=2)
    ax.set_title("Width to solution (recorded rows, research/pkm_dickson_torus_20261015)", fontsize=9)
    fig.tight_layout()
    fig.savefig(here / f"width_vs_n.{fmt}", dpi=150 if fmt == "png" else None)
    plt.close(fig)


for fmt in ("svg", "png"):
    d_plot(fmt)
    width_plot(fmt)
print("ok", file=sys.stderr)
