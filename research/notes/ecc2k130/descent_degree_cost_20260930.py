#!/usr/bin/env python3
"""What the Weil-descent route costs at ECC2K-130 as a function of the degree.

    python3 research/notes/ecc2k130/descent_degree_cost_20260930.py

Derived, not measured.  Every row is an evaluation of the cost model of
RESEARCH_ECC2K130_DECOMPOSITION.md §5.1 with one change: the oracle is a
Macaulay elimination of the Weil-descended chained S_3 system, priced per
target from below by the size of its degree-D matrix.

    relations      |F| = 2^l                      (l = dim V, an integer)
    targets        m! * 2^{n-(m-1)l}, at least one per relation
    oracle         C(N, <=D)^w per target, N = (m-2) n + m l unknowns
    linear algebra m * 2^{2l}

C(N, <=D) is the number of square-free monomials of degree <= D in N
unknowns: the column count of the degree-D Macaulay matrix.  w = 1 prices
one touch per column, which no elimination beats; w = 2 is the usual
estimate at omega = 2.  The unit is a column touch, and rho's 2^60.81 is
in group operations: the two are not converted here (see the note).

The degree scenarios:
  * D = 2: any Macaulay matrix at all -- a floor independent of the degree;
  * D = 6, 7: constant, at the values measured for m = 3 at l <= 5;
  * semi-regular D_reg(N) of the m = 3 shape (n cubic + n quadratic
    equations), which the measured degrees track within one at n <= 15.
Each is minimised over integer l.
"""
from math import comb, factorial, log2

N_FIELD = 131
RHO = 60.81  # RESEARCH_ECC2K130_DECOMPOSITION.md §5.2, <-1> x <pi>


def dreg(n_vars, degrees, cap=400):
    """Semi-regular degree of regularity over F_2 with field equations."""
    c = [comb(n_vars, i) for i in range(cap + 1)]
    for d in degrees:
        out = [0] * (cap + 1)
        for i in range(cap + 1):
            out[i] = c[i] - (out[i - d] if i >= d else 0)
        c = out
    return next(i for i, x in enumerate(c) if x <= 0)


def log2_add(a, b):
    return max(a, b) + log2(1 + 2 ** (-abs(a - b)))


def cost(m, l, degree, w):
    n = N_FIELD
    unknowns = (m - 2) * n + m * l
    per_relation = max(0.0, log2(factorial(m)) + n - m * l)
    targets = l + per_relation
    oracle = w * log2(sum(comb(unknowns, i) for i in range(degree + 1))) if degree else 0.0
    algebra = log2(m) + 2 * l
    return log2_add(targets + oracle, algebra), targets, oracle, algebra, unknowns


def best(m, degree_of, w):
    return min((cost(m, l, degree_of(l), w)[0], l) for l in range(4, 70))


def main():
    n = N_FIELD
    print("| m | degree D | w | best l | unknowns N | D | targets | oracle per target | linear algebra | total | vs rho |")
    print("|--:|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for m in (3, 4):
        scenarios = [("free oracle (floor)", lambda l: 0), ("2 (any Macaulay matrix)", lambda l: 2),
                     ("6, constant", lambda l: 6), ("7, constant", lambda l: 7)]
        if m == 3:
            scenarios.append(("semi-regular D_reg(N)", lambda l: dreg(n + 3 * l, [2] * n + [3] * n)))
        for name, degree_of in scenarios:
            for w in ((1,) if name.startswith("free") else (1, 2)):
                _, l = best(m, degree_of, w)
                total, targets, oracle, algebra, unknowns = cost(m, l, degree_of(l), w)
                print(f"| {m} | {name} | {w} | {l} | {unknowns} | {degree_of(l) or '—'} | 2^{targets:.2f} | "
                      f"2^{oracle:.2f} | 2^{algebra:.2f} | **2^{total:.2f}** | 2^{total - RHO:+.2f} |")
    print()
    print("semi-regular D_reg of the m = 3 shape at n = 131:",
          ", ".join(f"l={l}: {dreg(n + 3 * l, [2] * n + [3] * n)}" for l in (20, 26, 33, 43)))


if __name__ == "__main__":
    main()
