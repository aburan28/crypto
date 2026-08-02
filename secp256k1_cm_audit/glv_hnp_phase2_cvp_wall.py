"""
GLV-HNP Phase 2 — Thread 23c: confirmation run.  The Kannan embedding is the lossy step.

What Thread 23a/23b found today
-------------------------------
23a: quotienting out the trivial vector n*S_D*e_m does NOT move the K1 wall.
23b (m=8, 5 seeds): at the lam*=0.07 failure point, EXACT CVP on the rank-2m
     homogeneous lattice L0 recovered d where LLL+Kannan recovered nothing:

        curve         K1  eff    exactCVP d ok   LLL(B) d ok
        12-bit/2677    6  0.118      3/5             0/5
        12-bit/2677    8  0.157      2/5             0/5

     and the BKZ ladder beta = 2/10/20/30/40 on the Kannan lattice was completely
     FLAT -- stronger reduction buys nothing.

Taken together those say the wall is not reduction strength and not information:
it is the Kannan embedding.  Embedding the inhomogeneous problem as an SVP instance
turns "closest lattice point to -t" into "shortest vector in a coset", and the planted
vector is not the shortest vector of that lattice (2026-07-29 T5), so the embedding
loses the instance even when the underlying CVP is solvable.

This script re-runs the head-to-head at m=12 (matching the T4/X3 wall grid exactly)
with 10 seeds instead of 5, so the comparison is apples-to-apples with the wall table
and has enough samples to be worth quoting.

Also recorded: the non-uniqueness that makes "exact == planted" the wrong success
metric.  The GLV decomposition of a nonce is not unique -- (k1, k2) and
(k1 - lam + c*n, k2 + 1) give the same k_full -- so exact CVP legitimately returns a
SHORTER decomposition of the SAME d.  In 23b's table that showed up as
"exact=planted 0/5" together with "exactCVP d ok 5/5" at 12-bit/2557, K1=6.
Success must be measured on the recovered d, never on the recovered vector.

Run: python3 glv_hnp_phase2_cvp_wall.py
"""

import math
import os
import random
import time

from fpylll import CVP, IntegerMatrix, LLL, BKZ

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "glv_hnp_phase2_projected.py")
with open(_SRC) as _f:
    _text = _f.read()
_P = {'__file__': _SRC, '__name__': '_t23a'}
exec(compile(_text[:_text.index("# EXPERIMENT X1")].rsplit("# ====", 1)[0],
             _SRC, "exec"), _P)

scales = _P['scales']
gen_signatures = _P['gen_signatures']
find_generator = _P['find_generator']
lam_star = _P['lam_star']
norm = _P['norm']
d_from_k = _P['d_from_k']
build_homogeneous_lattice = _P['build_homogeneous_lattice']
build_projected_lattice = _P['build_projected_lattice']
recover_d_projected = _P['recover_d_projected']

SEEDS = [42, 1234, 9999, 555, 31337, 7, 20260802, 111, 65537, 314159]
K1_GRID = [4, 6, 8, 12, 16]
M = 12

HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]
curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    curves.append((label, (p, b, n, lam, G)))


def instance(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    return sigs, d_secret, k2_bound


def solve_cvp(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    Mh = build_homogeneous_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(Mh)
    LLL.reduction(A)
    rows = [r for r in ([A[i][j] for j in range(2 * m)] for i in range(A.nrows))
            if any(r)]
    B = IntegerMatrix.from_matrix(rows)
    t_pos = [0] * (2 * m)
    for i in range(m):
        t_pos[i] = (sigs[i]['A'] % n) * S_K1
    u = list(CVP.closest_vector(B, tuple(-x for x in t_pos)))
    v = [t_pos[j] + u[j] for j in range(2 * m)]
    if any(v[i] % S_K1 or v[m + i] % S_K2 for i in range(m)):
        return None
    cands = {d_from_k(v[i] // S_K1, v[m + i] // S_K2, sigs[i], n, lam)
             for i in range(m)}
    cands.discard(None)
    return cands.pop() if len(cands) == 1 else None


def solve_lattice(sigs, n, lam, k1_bound, k2_bound, beta=0):
    m = len(sigs)
    Mb = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(Mb)
    LLL.reduction(A)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    rows = [[A[i][j] for j in range(2 * m + 1)] for i in range(A.nrows)]
    return recover_d_projected(rows, sigs, n, lam, k1_bound, k2_bound)


print("=" * 78)
print("Thread 23c — exact CVP vs LLL/BKZ + Kannan, matched at m=12, 10 seeds")
print("=" * 78)
print(f"seeds = {SEEDS}")
print("Success is measured on the RECOVERED d, not on the recovered vector:")
print("the GLV decomposition of a nonce is not unique, so a correct attack may")
print("legitimately return a shorter (k1,k2) pair than the planted one.\n")

print(f"{'curve':<14} {'lam*':>7} {'K1':>3} {'eff':>6} "
      f"{'LLL+Kan':>9} {'BKZ40+Kan':>10} {'exact CVP':>10} {'cvp ms':>8}")
summary = {}
for label, curve in curves:
    p, b, n, lam, G = curve
    for k1 in K1_GRID:
        K2 = math.isqrt(n) + 1
        w_lll = w_bkz = w_cvp = tot = 0
        t_cvp = 0.0
        for s in SEEDS:
            inst = instance(curve, M, k1, s)
            if inst is None:
                continue
            sigs, d_secret, k2b = inst
            tot += 1
            w_lll += (solve_lattice(sigs, n, lam, k1, k2b, 0) == d_secret)
            w_bkz += (solve_lattice(sigs, n, lam, k1, k2b, 40) == d_secret)
            t0 = time.time()
            w_cvp += (solve_cvp(sigs, n, lam, k1, k2b) == d_secret)
            t_cvp += (time.time() - t0) * 1000
        summary[(label, k1)] = (w_lll, w_bkz, w_cvp, tot)
        print(f"{label:<14} {lam_star(lam,n):>7.4f} {k1:>3} {k1*K2/n:>6.3f} "
              f"{str(w_lll)+'/'+str(tot):>9} {str(w_bkz)+'/'+str(tot):>10} "
              f"{str(w_cvp)+'/'+str(tot):>10} {t_cvp/max(tot,1):>8.1f}")

print("\n" + "-" * 78)
print("Wall (largest K1 with >= 80% recovery), and the eff gained by exact CVP:")
print("-" * 78)
print(f"{'curve':<14} {'LLL+Kan':>9} {'BKZ40+Kan':>10} {'exact CVP':>10} {'gain':>7}")
for label, curve in curves:
    n = curve[2]
    K2 = math.isqrt(n) + 1
    walls = []
    for idx in range(3):
        w = 0
        for k1 in K1_GRID:
            a, bq, c, tot = summary[(label, k1)]
            got = (a, bq, c)[idx]
            if tot and got / tot >= 0.8:
                w = k1
        walls.append(w)
    gain = (walls[2] / walls[0]) if walls[0] else float('nan')
    print(f"{label:<14} {walls[0]:>9} {walls[1]:>10} {walls[2]:>10} {gain:>6.2f}x")
print(f"\nunicity bound at m={M}: eff < n^(-1/{M}) ~ {2647 ** (-1.0/M):.4f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
