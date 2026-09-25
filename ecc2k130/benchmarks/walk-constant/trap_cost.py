"""What fruitless cycles cost the table walk at ECC2K-130, WALK-CONSTANT.md section 6.

A walk steps until its point is distinguished (weight <= dpWeight, probability
theta per step) or until the guard restarts it after maxIters = M steps.  A
walk that enters a fruitless cycle (rate r per step, fruitless_patterns.py and
the measured rates) never becomes distinguished and runs to the guard.  Only
the points of a trail that ends distinguished are in the corpus, so the cost
is total steps over steps on distinguished trails, relative to a walk with
r = 0 and no guard:

    steps  = integral_0^M P(running at t) dt,
    P(running at t) = e^{-(theta+r)t} + r/(theta+r) (1 - e^{-(theta+r)t}),
    useful = integral_0^M t theta e^{-(theta+r)t} dt.

A guard below a few trail lengths discards honest trails; one far above them
lets a trapped walk run long.  best_guard() finds the balance.

This prices the loss only; it is a model built on the rates the harnesses
measure, and every number it prints is marked as such in the note.
"""
from math import comb, exp, log2

N_DEG = 131
M_DEFAULT = 2 ** 30          # aws/campaign.json "maxIters"
DP_WEIGHTS = (32, 34)        # aws/campaign.json "dpWeight" (the live bucket); the benchmarks' 34


def branch_probabilities(n, h):
    p = [0.0] * h
    for k in range(0, n + 1, 2):
        p[(k // 2) % h] += comb(n, k)
    s = sum(p)
    return [x / s for x in p]


def theta(n, w):
    even = sum(comb(n, k) for k in range(0, n + 1, 2))
    return sum(comb(n, k) for k in range(0, w + 1, 2)) / even


def overhead(r, th, M):
    a = th + r
    e = exp(-a * M)
    run = (1 - e) / a + (r / a) * (M - (1 - e) / a)
    useful = th * (1 - e * (1 + a * M)) / (a * a)   # integral of t th e^{-at}
    return run / useful                              # r = 0, M = inf gives 1


def best_guard(r, th):
    grid = [2 ** (e / 8) for e in range(8 * 24, 8 * 32 + 1)]
    return min(((overhead(r, th, M), M) for M in grid))


if __name__ == "__main__":
    two_m = 2 * N_DEG
    for dp_weight in DP_WEIGHTS:
        th = theta(N_DEG, dp_weight)
        print(f"n = {N_DEG}, dpWeight {dp_weight}: distinguished with probability 2^{log2(th):.2f} per step, "
              f"trail 2^{-log2(th):.2f}; guard maxIters = 2^{log2(M_DEFAULT):.0f}")
        for H in (8, 16):
            p = branch_probabilities(N_DEG, H)
            s2, s4 = sum(x * x for x in p), sum(x ** 4 for x in p)
            pair = 4 * (s2 / two_m) ** 3
            rel = 24 * s4 / two_m ** 3
            print(f"  H = {H}: sum p^2 = {s2:.5f}, sum p^4 = {s4:.6f}; pairwise 6-cycles {pair:.3e}/step, "
                  f"tau-relation 4-cycles {rel:.3e}/step")
            for label, r in (("as built (rule: 2- and 4-step pairwise)", pair + rel),
                             ("rule also refuses the tau-relation 4-cycles", pair),
                             ("no fruitless cycles (sigma walk)", 0.0)):
                f = r / (r + th)
                ov, M = best_guard(r, th)
                print(f"    {label}: {100 * f:.2f}% of trails trapped; cost x{overhead(r, th, M_DEFAULT):.3f} "
                      f"at maxIters 2^30, x{ov:.3f} at the best guard 2^{log2(M):.2f}")
        print()
