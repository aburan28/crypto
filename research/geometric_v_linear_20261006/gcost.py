"""Corrected n = 131 cost map with the geometric-V linear (A, P) oracle.

Same heuristics and budget convention as ../linearization_reach_20260930/costmap_a1.py:
Pr[decompose] = min(1, X^m/(m! 2^n)), relations = columns 2^c, LA = m * 2^(2c),
per-attempt budget 2^(RHO - c) once Pr ~ 1.  Per trial the oracle guesses the first
m - 2 summands and resolves the last pair fully (b = c) by one n x (3c - 1) solve,
valid while 3c - 1 <= n; past that the solution family 2^(3c-1-n) is enumerated and
charged.  Per relation: 2^(n - 2c) trials (pair success 2^(2c - n - 1), two signs).
    gap_bits    = log2(trials per relation) - log2(budget)
    gap_charged = gap_bits + omega * log2(3c - 1) + max(0, 3c - 1 - n)
Random V (the frozen thread's base) is printed alongside for comparison.
"""
import json
import math
import sys

N, RHO = 131, 60.8090


def lf(m):
    return math.lgamma(m + 1) / math.log(2)


def rows(omega):
    out = []
    for orbit in (False, True):
        for m in (3, 4, 5, 6, 8, 12):
            best = None
            for c in range(6, 80):
                X = c + (math.log2(N) if orbit else 0.0)
                p = min(0.0, m * X - lf(m) - N)
                attempts = c - p
                la = math.log2(m) + 2 * c
                if la >= RHO:
                    break
                budget = math.log2(2 ** RHO - 2 ** la) - attempts
                trials = N - 2 * c          # per relation, both summands of the pair resolved
                U = 3 * c - 1
                gb = trials - budget
                gc = gb + omega * math.log2(U) + max(0, U - N)
                row = dict(orbit=orbit, m=m, c=c, log2_budget=round(budget, 2), log2_trials=trials,
                           unknowns=U, gap_bits=round(gb, 2), gap_charged=round(gc, 2))
                if best is None or row["gap_charged"] < best["gap_charged"]:
                    best = row
            out.append(best)
    return out


if __name__ == "__main__":
    res = {}
    for omega in (2.0, 2.81):
        rs = rows(omega)
        res[str(omega)] = rs
        print(f"omega={omega}")
        print("| base | m | c | budget | trials/relation | unknowns | gap bits | gap charged |")
        print("|---|---:|---:|---:|---:|---:|---:|---:|")
        for r in rs:
            print(f"| {'orbit union' if r['orbit'] else 'subspace'} | {r['m']} | {r['c']} | 2^{r['log2_budget']} | "
                  f"2^{r['log2_trials']} | {r['unknowns']} | {r['gap_bits']} | {r['gap_charged']} |")
        bb = min(rs, key=lambda r: r["gap_bits"])
        bc = min(rs, key=lambda r: r["gap_charged"])
        print(f"best bits-only {bb['gap_bits']} (m={bb['m']}, c={bb['c']}); best charged {bc['gap_charged']} (m={bc['m']}, c={bc['c']})\n")
    if len(sys.argv) > 1:
        json.dump(res, open(sys.argv[1], "w"), indent=1)
