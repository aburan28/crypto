#!/usr/bin/env python3
"""The descent_degree_cost model (examples/descent_degree_cost.rs, same
formulas) evaluated at the constant degrees D = 3, 4, 5 that the F4 arm
measures, n = 131, l in 4..70.  Derived, not measured: extrapolation."""
from math import comb, log2, factorial
RHO = 60.81
def la(a, b): return max(a, b) + log2(1 + 2 ** (-abs(a - b)))
def mono(N, d): return log2(sum(comb(N, i) for i in range(d + 1)))
n = 131
print("| m | D | w | best l | N | total | vs rho |")
print("|--:|--:|--:|--:|--:|--:|--:|")
for m in (3, 4):
    for D in (3, 4, 5, 6):
        for w in (1, 2):
            best = None
            for l in range(4, 71):
                per = max(0.0, log2(factorial(m)) + n - m * l)
                targets = l + per
                N = (m - 2) * n + m * l
                tot = la(targets + w * mono(N, D), log2(m) + 2 * l)
                if best is None or tot < best[0]:
                    best = (tot, l, N)
            print(f"| {m} | {D} | {w} | {best[1]} | {best[2]} | 2^{best[0]:.2f} | 2^{best[0]-RHO:+.2f} |")
