#!/usr/bin/env python3
"""Scheduling estimate for grid cells: an EXTRAPOLATION, never a measurement.

    python3 research/dreg_ell_grid_20260925/cost_model.py

Macaulay shape of the chained S3 system (m = 3, N = n + 3*ell unknowns,
n equations of degree 3 and n of degree 2) at degree D:
    rows = n * (C(N, <=D-2) + C(N, <=D-3))     (matches every recorded row count)
    cols <= C(N, <=D)                            (the builder keeps occurring monomials only)
Time is fitted as t = c * (rows * cols)^alpha against the ladder's measured
single-draw wall times on the four-core container.  Two fits bracket the
estimate: alpha from the two degree-6 cells alone, and alpha from those two
plus the degree-5 probe.  Neither is evidence about any degree.
"""
import json
import statistics
from math import comb, exp, log
from pathlib import Path

LADDER = Path(__file__).resolve().parent.parent / "dreg_fixed_surplus_20260923" / "runs"


def cum(N, k):
    return sum(comb(N, i) for i in range(k + 1))


def shape(n, ell, D):
    N = n + 3 * ell
    return n * (cum(N, D - 2) + cum(N, D - 3)), cum(N, D), N


def median_secs(name):
    rows = [json.loads(l) for l in (LADDER / name).read_text().splitlines() if l.strip()]
    return statistics.median(r["secs"] for r in rows if r.get("outcome", {}).get("kind") == "resolved")


POINTS = [
    ((11, 4, 6), median_secs("cell-11-4-6.jsonl")),
    ((13, 4, 6), median_secs("cell-13-4-6.jsonl")),
    ((15, 4, 5), 60.0),  # degree-5 probe, PREREGISTRATION.md of the ladder, "the cost probe" table
]


def fit(points):
    xs = [log(shape(*c)[0] * shape(*c)[1]) for c, _ in points]
    ys = [log(t) for _, t in points]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return a, my - a * mx


def fmt(t):
    return f"{t / 3600:.1f} h" if t >= 3600 else f"{t / 60:.1f} min" if t >= 60 else f"{t:.0f} s"


def main():
    fits = {"two degree-6 cells": fit(POINTS[:2]), "plus the degree-5 probe": fit(POINTS)}
    print("timings:", {f"{c}": round(t) for c, t in POINTS})
    for name, (a, _) in fits.items():
        print(f"alpha ({name}): {a:.2f}")
    print("\n| cell (n, ell, D) | N | S | rows | cols (upper) | estimate a draw |")
    print("|---|--:|--:|--:|--:|---|")
    cells = [(4, 3, 7), (5, 3, 7), (8, 2, 7), (10, 2, 7), (12, 2, 7), (10, 5, 6),
             (11, 4, 6), (13, 4, 6), (13, 5, 6), (15, 5, 6)]
    for c in cells:
        r, k, N = shape(*c)
        est = [exp(b + a * log(r * k)) for a, b in fits.values()]
        print(f"| {c} | {N} | {c[0] - 3 * c[1]:+d} | {r:,} | {k:,} | {fmt(min(est))} – {fmt(max(est))} |")
    print("\nDegree-7 rows are the cost if the cell reaches degree 7; most resolve lower.")


if __name__ == "__main__":
    main()
