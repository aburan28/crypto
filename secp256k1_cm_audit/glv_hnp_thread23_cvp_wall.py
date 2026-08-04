"""
Thread 23, part 2 — separate the *information-theoretic* wall from the
*reduction-quality* wall in the GLV-HNP Phase-2 lattice.

Part 1 (glv_hnp_thread23_bdd.py, EXP U3) found that on the historical failure
curve 12-bit/2677 (lam*=0.070) exact CVP by enumeration recovers d at
K1 = 6 (4/5) and K1 = 8 (1/5), where Kannan+LLL, Kannan+BKZ-20, LLL+Babai and
BKZ20+Babai all score 0/5.  So the K1~4-6 wall reported on 2026-07-26 and
re-measured on 2026-07-29 (EXP T4) is NOT information-theoretic.

This script pins down the two walls:

  U5  exact-CVP wall.  For each (curve, K1) run exact CVP and record
        - recovery,
        - dist(closest) vs ||e_true||.
      dist < ||e_true||  =>  a decoy is genuinely closer than the planted
      error: the instance is information-theoretically lost, no algorithm
      that decodes to the closest vector can win.
      dist == ||e_true|| but no recovery => tie (measure-zero, report it).

  U6  reduction-quality ladder at the critical K1: Kannan+LLL, Kannan+BKZ-b
      for b in {20,30,40,50,60}, LLL+Babai, BKZ20+Babai, exact CVP.
      Shows how much blocksize is needed to reach the CVP wall.

  U7  the exact nearest-plane criterion.  Babai returns the correct vector iff
      |<e, b*_i>| / ||b*_i||^2 < 1/2 for every i (given the reduced basis).
      Report max_i of that quantity; compare with the loose global bound
      ||e|| < min_i||b*_i||/2 used in U4 (which was useless: ratio 3-6 while
      winning 5/5, because min||b*_i|| is pinned to n by the trivial vector).

Run: python3 glv_hnp_thread23_cvp_wall.py
"""

import math
import random

from fpylll import IntegerMatrix, LLL, BKZ, CVP

from glv_hnp_thread23_bdd import (
    find_generator, lam_star, gen_signatures, scales,
    build_kannan, planted_kannan, recover_kannan,
    build_bdd, error_vector, recover_bdd,
    babai_nearest_plane, gram_schmidt, norm,
)

SEEDS = [42, 1234, 9999, 555, 31337]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

HIST = [
    ("8-bit/199",     211,  2, 199,  106,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8),
    ("12-bit/2677",   2677, 2, 2647, 185,  10),
]

print("=" * 78)
print("Thread 23 part 2 — information wall vs reduction wall")
print("=" * 78)
CURVES = []
for label, p, b, n, lam, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    CURVES.append((label, (p, b, n, lam, G), m))
    print(f"  {label:14s} n={n:5d}  lam={lam:5d}  lam*={lam_star(lam,n):.4f}  m={m}")


def cvp_trial(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M, t = build_bdd(sigs, n, lam, k1_bound, k2_bound, 1)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    w = list(CVP.closest_vector(A, tuple(t)))
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, 1)
    ok = recover_bdd(w, t, m, n, S_D, d_secret) is not None
    dist = norm([w[j] - t[j] for j in range(dim)])
    e_true = error_vector(sigs, d_secret, n, k1_bound, k2_bound, 1, reduced=True)
    e_raw = error_vector(sigs, d_secret, n, k1_bound, k2_bound, 1, reduced=False)
    en = min(norm(e_true), norm(e_raw))
    return dict(ok=ok, dist=dist, en=en, n=n)


print("\n" + "=" * 78)
print("EXP U5 — exact-CVP wall.  'closer' = a decoy is strictly closer to the")
print("         target than the planted error.  Necessary but NOT sufficient")
print("         for failure: a closer decoy sharing the d-coordinate mod n")
print("         still decodes correctly.  Verdict follows the CVP column.")
print("=" * 78)
print("curve            K1    CVP    dist/n   ||e||/n   dist<||e||  verdict")
print("-" * 78)
U5 = {}
for label, curve, m in CURVES:
    for k1 in K1_GRID:
        wins, dists, ens, closer = 0, [], [], 0
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, curve[2] - 1)
            r = cvp_trial(curve, m, d, k1, s)
            if r is None:
                continue
            wins += bool(r['ok'])
            dists.append(r['dist'] / r['n'])
            ens.append(r['en'] / r['n'])
            if r['dist'] < r['en'] - 1e-9:
                closer += 1
        md = sum(dists) / len(dists)
        me = sum(ens) / len(ens)
        # NB: "a strictly closer decoy exists" is necessary but NOT sufficient
        # for failure — a closer decoy that happens to share the d-coordinate
        # mod n still decodes correctly (e.g. any coset shift by the trivial
        # vector).  So the verdict is driven by the CVP recovery count, and the
        # 'closer' column is reported as the diagnostic it is.
        if wins == len(SEEDS):
            verdict = "CVP OK"
        elif closer == len(SEEDS) and wins == 0:
            verdict = "CVP LOST"
        else:
            verdict = "mixed"
        U5[(label, k1)] = wins
        print(f"{label:16s} {k1:<4d} {wins}/5    {md:6.3f}   {me:7.3f}   "
              f"{closer}/5        {verdict}")
    print()


print("=" * 78)
print("EXP U6 — reduction ladder at the critical K1")
print("=" * 78)


def kannan_rate(curve, m, k1, seeds, beta=None):
    n = curve[2]
    wins = 0
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        p, b, nn, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2_bound, s)
        if len(sigs) < m:
            continue
        M = build_kannan(sigs, n, lam, k1, k2_bound, 1)
        dim = 2 * m + 2
        A = IntegerMatrix.from_matrix(M)
        if beta is None:
            LLL.reduction(A)
        else:
            BKZ.reduction(A, BKZ.Param(beta))
        rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2_bound, 1)
        wins += recover_kannan(rows, m, n, S_KAN, S_D, d) is not None
    return wins


def babai_rate(curve, m, k1, seeds, beta=None):
    n = curve[2]
    wins = 0
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        p, b, nn, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2_bound, s)
        if len(sigs) < m:
            continue
        M, t = build_bdd(sigs, n, lam, k1, k2_bound, 1)
        dim = 2 * m + 1
        A = IntegerMatrix.from_matrix(M)
        if beta is None:
            LLL.reduction(A)
        else:
            BKZ.reduction(A, BKZ.Param(beta))
        rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
        w = babai_nearest_plane(rows, t)
        S_K1, S_D, S_K2, _ = scales(n, k1, k2_bound, 1)
        wins += recover_bdd(w, t, m, n, S_D, d) is not None
    return wins


CRIT = {"8-bit/199": [4, 6, 8], "12-bit/2557": [12, 16], "12-bit/2677": [4, 6, 8]}
for label, curve, m in CURVES:
    print(f"\n{label} (m={m})")
    print("  K1  | K+LLL  K+B20  K+B30  K+B40  K+B50  K+B60 | Bab   B20+Bab | CVP")
    print("  " + "-" * 72)
    for k1 in CRIT[label]:
        row = [kannan_rate(curve, m, k1, SEEDS)]
        for beta in (20, 30, 40, 50, 60):
            row.append(kannan_rate(curve, m, k1, SEEDS, beta=beta))
        bab = babai_rate(curve, m, k1, SEEDS)
        bab20 = babai_rate(curve, m, k1, SEEDS, beta=20)
        print(f"  {k1:<3d} | " + "  ".join(f"{v}/5  " for v in row) +
              f"| {bab}/5   {bab20}/5     | {U5[(label,k1)]}/5")


print("\n" + "=" * 78)
print("EXP U7 — exact nearest-plane criterion  max_i |<e,b*_i>|/||b*_i||^2 < 1/2")
print("=" * 78)
print("curve            K1   max_i coef   #i over 1/2   Babai   (U4 loose ratio)")
print("-" * 78)
for label, curve, m in CURVES:
    p, b, n, lam, G = curve
    for k1 in K1_GRID:
        maxes, overs, wins = [], [], 0
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, n - 1)
            k2_bound = math.isqrt(n) + 1
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2_bound, s)
            if len(sigs) < m:
                continue
            M, t = build_bdd(sigs, n, lam, k1, k2_bound, 1)
            dim = 2 * m + 1
            A = IntegerMatrix.from_matrix(M)
            LLL.reduction(A)
            rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
            Bs, _ = gram_schmidt(rows)
            e = error_vector(sigs, d, n, k1, k2_bound, 1, reduced=True)
            e2 = error_vector(sigs, d, n, k1, k2_bound, 1, reduced=False)
            if norm(e2) < norm(e):
                e = e2
            coefs = []
            for i in range(dim):
                nb = sum(y * y for y in Bs[i])
                if nb == 0:
                    continue
                coefs.append(abs(sum(float(e[j]) * Bs[i][j] for j in range(dim))) / nb)
            maxes.append(max(coefs))
            overs.append(sum(1 for c in coefs if c > 0.5))
            w = babai_nearest_plane(rows, t)
            S_K1, S_D, S_K2, _ = scales(n, k1, k2_bound, 1)
            wins += recover_bdd(w, t, m, n, S_D, d) is not None
        print(f"{label:16s} {k1:<4d} {sum(maxes)/len(maxes):10.3f}   "
              f"{sum(overs)/len(overs):11.1f}   {wins}/5")
    print()

print("=" * 78)
print("done")
print("=" * 78)
