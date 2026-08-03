"""
Thread 23 (part 4) — per-seed test of the Babai/GS-profile criterion.

Part 3 (glv_hnp_phase2_babai.py) compared SEED-AVERAGED beta against
seed-averaged success, and found the flip bracketed in (0.666, 0.851) for
12-bit/2557 but (0.435, 0.471) for 12-bit/2677 -- much tighter than the
lambda_1 predictor gamma (2.13/2.31 vs 0.71/0.86, a factor ~3 apart) but still
not a single threshold.

Averaging beta over seeds is the wrong comparison: the Babai condition

    beta = max_i |<v_planted, b*_i>| / ||b*_i||^2  <  1/2

is a per-INSTANCE statement, so it must be tested per instance.  This script
pairs (beta, recovered) seed by seed and asks:

  E10  does beta < 1/2 imply recovery, instance by instance, pooled over both
       curves and the whole K1 grid?  Report violations in both directions and
       the separation quality.

Run: python3 glv_hnp_phase2_beta_perseed.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

from glv_hnp_phase2_projected import (
    find_generator, lam_star, scales, gen_signatures,
    build_projected_lattice, planted_vector_projected, recover_d_projected,
    norm,
)
from glv_hnp_phase2_babai import gs_profile

SEEDS = [42, 1234, 9999, 555, 31337, 8675309, 271828, 141421, 1618033, 6022140]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755, 8),
    ("12-bit/2677", 2677, 2, 2647, 185,  10),
]
curves = []
for label, p, b, n, lam, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0
    curves.append((label, (p, b, n, lam, G), m))


def instance(curve, m, d_secret, k1_bound, seed, beta_bkz=30):
    """One instance: returns (beta, recovered) on the SAME signature set."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m:
        return None
    rows = build_projected_lattice(sigs, n, lam, k1_bound, k2b)
    cols = 2 * m + 1

    # (a) beta on the Kannan-free sublattice L'_0
    A0 = IntegerMatrix.from_matrix(rows[:-1])
    BKZ.reduction(A0, BKZ.Param(beta_bkz))
    red0 = [r for r in ([A0[i][j] for j in range(cols)] for i in range(A0.nrows))
            if any(r)]
    t = planted_vector_projected(sigs, n, k1_bound, k2b)
    beta = 0.0
    for bs in gs_profile(red0):
        nb = sum(x * x for x in bs)
        if nb > 0:
            beta = max(beta, abs(sum(float(a) * c for a, c in zip(t, bs))) / nb)

    # (b) the actual attack on the full projected lattice
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(cols)] for i in range(A.nrows)]
    ok = recover_d_projected(red, sigs, n, lam, k1_bound, k2b) == d_secret
    return beta, ok


print("=" * 78)
print("Thread 23 part 4 — per-instance test of the Babai criterion beta < 1/2")
print("=" * 78)

data = []   # (label, K1, seed, beta, ok)
for label, cur, m in curves:
    for k1 in K1_GRID:
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, cur[2] - 1)
            r = instance(cur, m, d, k1, s)
            if r is None:
                continue
            data.append((label, k1, s, r[0], r[1]))

print(f"\ncollected {len(data)} instances "
      f"({len(curves)} curves x {len(K1_GRID)} K1 x {len(SEEDS)} seeds)")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("E10 — does beta < 1/2 imply recovery, per instance?")
print("=" * 78)

viol_a = [d for d in data if d[3] < 0.5 and not d[4]]      # beta<0.5 but FAILED
viol_b = [d for d in data if d[3] >= 0.5 and d[4]]         # beta>=0.5 but WORKED
n_lo = sum(1 for d in data if d[3] < 0.5)
n_hi = len(data) - n_lo
print(f"  instances with beta <  0.5 : {n_lo:4d}   of which recovered: "
      f"{n_lo - len(viol_a):4d}   ({100*(n_lo-len(viol_a))/max(n_lo,1):.1f}%)")
print(f"  instances with beta >= 0.5 : {n_hi:4d}   of which recovered: "
      f"{len(viol_b):4d}   ({100*len(viol_b)/max(n_hi,1):.1f}%)")
print(f"\n  sufficiency violations (beta<0.5 yet failed): {len(viol_a)}")
for d in viol_a[:10]:
    print(f"    {d[0]}  K1={d[1]:<3} seed={d[2]:<8} beta={d[3]:.4f}")
print(f"\n  successes beyond the Babai guarantee (beta>=0.5 yet worked): {len(viol_b)}")
for d in sorted(viol_b, key=lambda x: -x[3])[:10]:
    print(f"    {d[0]}  K1={d[1]:<3} seed={d[2]:<8} beta={d[3]:.4f}")

# separation quality
succ = sorted(d[3] for d in data if d[4])
fail = sorted(d[3] for d in data if not d[4])
if succ and fail:
    print(f"\n  beta | recovered  : n={len(succ):3d}  min={succ[0]:.3f} "
          f"median={succ[len(succ)//2]:.3f}  max={succ[-1]:.3f}")
    print(f"  beta | failed     : n={len(fail):3d}  min={fail[0]:.3f} "
          f"median={fail[len(fail)//2]:.3f}  max={fail[-1]:.3f}")
    pairs = sum(1 for a in succ for b in fail if a < b)
    ties = sum(1 for a in succ for b in fail if a == b)
    auc = (pairs + 0.5 * ties) / (len(succ) * len(fail))
    print(f"  AUC (beta separates failure from success): {auc:.4f}")
    # best single threshold
    cands = sorted(set(succ + fail))
    best = max(((sum(1 for x in succ if x < c) + sum(1 for x in fail if x >= c)) / len(data), c)
               for c in cands)
    print(f"  best empirical threshold: beta* = {best[1]:.4f}  "
          f"(accuracy {100*best[0]:.1f}%)")
    base = max(len(succ), len(fail)) / len(data)
    print(f"  majority baseline accuracy: {100*base:.1f}%")

# per-curve view
print("\n  --- per-curve breakdown at the empirical threshold ---")
for label, _cur, _m in curves:
    sub = [d for d in data if d[0] == label]
    s = sorted(d[3] for d in sub if d[4]); f = sorted(d[3] for d in sub if not d[4])
    print(f"  {label}: max beta among successes = {max(s) if s else float('nan'):.3f}, "
          f"min beta among failures = {min(f) if f else float('nan'):.3f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
