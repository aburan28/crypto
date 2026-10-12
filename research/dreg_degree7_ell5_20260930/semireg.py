#!/usr/bin/env python3
"""The semi-regular degree of regularity of each measured cell's shape.

    python3 research/dreg_degree7_ell5_20260930/semireg.py

A derived reference, not a measurement.  Over F_2 with the field equations,
a semi-regular sequence of polynomials of degrees d_1..d_k in N unknowns
has Hilbert series (1+z)^N / prod(1+z^{d_i}), truncated before its first
non-positive coefficient (Bardet, Faugère, Salvy and Yang), and D_reg is
the index of that coefficient.  The chained S_3 system of a cell (n, l) has
N = n + 3l unknowns and n cubic plus n quadratic equations.

It prints the measured refutation degrees (committed RESULTS files, as
quoted in RESEARCH_DREG_MEASUREMENT.md Results 4-8) next to D_reg.
"""
from math import comb


def dreg(n_vars, degrees, cap=400):
    c = [comb(n_vars, i) for i in range(cap + 1)]
    for d in degrees:
        out = [0] * (cap + 1)
        for i in range(cap + 1):
            out[i] = c[i] - (out[i - d] if i >= d else 0)
        c = out
    return next(i for i, x in enumerate(c) if x <= 0)


def cell_dreg(n, ell):
    return dreg(n + 3 * ell, [2] * n + [3] * n)


# (n, l) -> the four committed values
MEASURED = {
    (5, 2): "5 5 5 5", (7, 2): "5 5 5 5", (8, 2): "5 5 5 5", (10, 2): "5 5 5 5", (12, 2): "5 5 5 5",
    (4, 3): "5 5 5 5", (5, 3): "6 6 6 6", (7, 3): "6 6 6 6", (9, 3): "6 6 6 6",
    (7, 4): "7 7 6 7", (8, 4): "6 7 6 7", (11, 4): "6 6 6 6", (13, 4): "6 6 6 6",
    (10, 5): "≥7 ≥7 ≥7 ≥7", (11, 5): "≥7 ≥7 ≥7 ≥7", (13, 5): "≥7 ≥7 ≥7 ≥7", (15, 5): "≥7 ≥7 ≥7 ≥7",
}

if __name__ == "__main__":
    print("| cell | N | S | measured | semi-regular D_reg |")
    print("|---|--:|--:|---|--:|")
    for (n, ell), v in sorted(MEASURED.items(), key=lambda t: (t[0][1], t[0][0] - 3 * t[0][1])):
        print(f"| ({n}, {ell}) | {n + 3 * ell} | {n - 3 * ell:+d} | {v} | {cell_dreg(n, ell)} |")
