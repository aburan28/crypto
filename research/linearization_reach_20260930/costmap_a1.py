"""(Amendment 1 b_max) ECC2K-130 decomposition cost map and the linearization gap (arithmetic only).

For arity m and factor-base dimension c (optionally a 131-orbit union), print the
best cell under rho, the per-attempt oracle budget it leaves, and the per-relation
cost of the guess-and-linearize oracle, 2^(n - c - bmax(c)), measured law of lr.py.
Heuristics (stated, not proved): Pr[decompose] = min(1, X^m / (m! 2^n)) with X the
number of base abscissae; relations = columns; linear algebra = m * columns^2.
"""
import json
import math

RHO, N = 60.8090, 131


def lfact(m):
    return math.lgamma(m + 1) / math.log(2)


def bmax(c):
    # AMENDMENT-1: largest b <= c with c(b+1) + b - b(b+1)/2 <= N
    return max(b for b in range(0, c + 1) if c * (b + 1) + b - b * (b + 1) // 2 <= N)


rows = []
for orbit in (False, True):
    for m in (3, 4, 5, 6, 8, 10, 12, 16):
        best = None
        for c in range(8, 80):
            X = c + (math.log2(N) if orbit else 0.0)
            p = min(0.0, m * X - lfact(m) - N)
            attempts = c - p
            la = math.log2(m) + 2 * c
            if la >= RHO:
                break
            budget = math.log2(2 ** RHO - 2 ** la) - attempts
            lin = N - c - bmax(c)          # log2 cost per relation, measured law
            cell = dict(orbit=orbit, m=m, c=c, log2_pr=round(p, 2), log2_attempts=round(attempts, 2),
                        log2_la=round(la, 2), log2_budget=round(budget, 2), bmax=bmax(c),
                        log2_lin_per_relation=lin, gap_bits=round(lin - budget, 2),
                        total_lin=round(c + lin, 2))
            if best is None or cell["gap_bits"] < best["gap_bits"]:
                best = cell
        rows.append(best)
print("| base | m | c | Pr[dec] | attempts | LA | oracle budget | bmax | lin oracle / relation | gap (bits) |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for r in rows:
    print(f"| {'orbit union' if r['orbit'] else 'subspace'} | {r['m']} | {r['c']} | 2^{r['log2_pr']} | 2^{r['log2_attempts']} "
          f"| 2^{r['log2_la']} | 2^{r['log2_budget']} | {r['bmax']} | 2^{r['log2_lin_per_relation']} | {r['gap_bits']} |")
print()
print("bmax(l) at n=131 around n/3:", {l: bmax(l) for l in range(38, 50)})
print("max_c bmax(c) =", max(bmax(c) for c in range(1, N)), "at c =", [c for c in range(1, N) if bmax(c) == max(bmax(k) for k in range(1, N))])
print("gap floor when Pr[dec] = 1: n - log2(rho) - max bmax =", round(N - RHO - max(bmax(c) for c in range(1, N)), 2))
json.dump(rows, open("results/costmap_a1.json", "w"), indent=1)
