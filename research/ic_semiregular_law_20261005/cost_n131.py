#!/usr/bin/env python3
"""Index-calculus cost at the ECC2K-130 size, m = 3, IF the refutation degree follows the
semi-regular law (window [d_reg, d_reg + 1]).  Registered before any held-out reading.

Model (all logs base 2; every term is a stated assumption, not a measurement):
  r          subgroup order ~ 2^129 (n = 131 Koblitz curve, cofactor 4)
  F          factor base size ~ 2^l points (x in a random l-dimensional subspace V)
  p          P(random point = sum of 3 factor-base points) ~ F^3 / (3! r)
  relations  ~ F;  attempts ~ F / p
  per attempt one refutation-sized Macaulay solve at degree D in v unknowns:
             cost ~ N(D)^omega bit operations, N(D) = sum_{k<=D} C(v, k) columns,
             omega in {2 (sparse, optimistic), 2.807 (Strassen-dense)}
  linear algebra ~ F^2 (sparse, optimistic)
  arms: x4 (v = 3l, n equations of degree 6); rr (v = 3l + 1, 2n equations of degree 4 and 1 linear)
  rho baseline: 2^60.9 group operations (published ECC2K-130 estimate), counted as bit
  operations x 1 (generous to IC: a group operation costs >= 2^12 bit operations).
Prints, for each arm, omega and window end, the best l and its log2 total.
"""
import math

N_BITS, LOG_R, RHO = 131, 129.0, 60.9


def d_reg(v, degs):
    h = [math.comb(v, k) for k in range(v + 2)]
    for d in degs:
        out = [0] * len(h)
        for k in range(len(h)):
            out[k] = h[k] - (out[k - d] if k >= d else 0)
        h = out
    return next((k for k, c in enumerate(h) if c <= 0), v + 1)


def cols(v, D):
    return sum(math.comb(v, k) for k in range(min(D, v) + 1))


def shape(arm, l):
    if arm == "x4":
        return 3 * l, [6] * N_BITS
    return 3 * l + 1, [4] * (2 * N_BITS) + [1]


def total(arm, l, omega, plus):
    v, degs = shape(arm, l)
    D = d_reg(v, degs) + plus
    log_attempts = l - (3 * l - math.log2(6) - LOG_R)  # F / p
    log_solve = omega * math.log2(cols(v, D))
    log_la = 2 * l
    t = max(log_attempts + log_solve, log_la)
    return t, D, log_attempts, log_solve


if __name__ == "__main__":
    print(f"rho baseline 2^{RHO}")
    for arm in ("x4", "rr"):
        for omega in (2.0, 2.807):
            for plus in (0, 1):
                best = min((total(arm, l, omega, plus) + (l,) for l in range(1, 44)), key=lambda t: t[0])
                t, D, la, ls, l = best
                print(f"{arm} omega={omega} D=d_reg+{plus}: best l={l} D={D} "
                      f"attempts 2^{la:.1f} x solve 2^{ls:.1f} = 2^{t:.1f}  ({'below' if t < RHO else 'above'} rho by 2^{abs(t - RHO):.1f})")
