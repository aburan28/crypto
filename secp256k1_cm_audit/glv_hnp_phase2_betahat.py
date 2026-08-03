"""
Thread 23 (part 5) — a SECRET-FREE proxy for the Babai criterion.

Part 4 (glv_hnp_phase2_beta_perseed.py) established, over 160 instances:

    beta = max_i |<v_planted, b*_i>| / ||b*_i||^2      (AUC 0.961)
    beta < 0.466 -> 50/50 recovered;  beta > 0.934 -> 0/N recovered.

That is the best instance-level predictor this project has found (previous best:
nu_hat, AUC 0.935, 2026-07-29).  But beta is a DIAGNOSTIC, not an attack tool:
it needs v_planted, i.e. the secret.

The isotropic approximation removes that dependence.  For an error e of norm R
in dimension D, |<e, b*_i>| ~ R*||b*_i||/sqrt(D), so the Babai condition
|<e,b*_i>|/||b*_i||^2 < 1/2 becomes R < sqrt(D) * min_i ||b*_i|| / 2, i.e.

    beta_hat = E||v_planted|| / ( sqrt(D) * min_i ||b*_i|| )   <  1/2

where E||v_planted|| comes from the PUBLIC bounds K1, K2 and the scale factors,
and the GS profile of L'_0 is public too.  Nothing secret enters.

  E11  beta_hat vs observed recovery, per instance, both curves, full K1 grid.
       Compare AUC against beta (0.961) and nu_hat (0.935).

Run: python3 glv_hnp_phase2_betahat.py
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


def expected_planted_norm(m, n, k1_bound, k2_bound):
    """E||v_planted|| in the PROJECTED lattice, from public bounds only."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                     + m * (k2_bound * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)


def instance(curve, m, d_secret, k1_bound, seed, beta_bkz=30):
    """Returns (beta, beta_hat, recovered) on one signature set."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m:
        return None
    rows = build_projected_lattice(sigs, n, lam, k1_bound, k2b)
    cols = 2 * m + 1

    A0 = IntegerMatrix.from_matrix(rows[:-1])
    BKZ.reduction(A0, BKZ.Param(beta_bkz))
    red0 = [r for r in ([A0[i][j] for j in range(cols)] for i in range(A0.nrows))
            if any(r)]
    bstar = gs_profile(red0)

    # secret-dependent beta
    t = planted_vector_projected(sigs, n, k1_bound, k2b)
    beta = 0.0
    for bs in bstar:
        nb = sum(x * x for x in bs)
        if nb > 0:
            beta = max(beta, abs(sum(float(a) * c for a, c in zip(t, bs))) / nb)

    # secret-FREE beta_hat
    gsn = [math.sqrt(sum(x * x for x in bs)) for bs in bstar]
    D = len(bstar)
    beta_hat = expected_planted_norm(m, n, k1_bound, k2b) / (math.sqrt(D) * min(gsn))

    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(cols)] for i in range(A.nrows)]
    ok = recover_d_projected(red, sigs, n, lam, k1_bound, k2b) == d_secret
    return beta, beta_hat, ok


def auc_and_threshold(vals_ok, vals_bad):
    if not vals_ok or not vals_bad:
        return float('nan'), float('nan'), float('nan')
    pairs = sum(1 for a in vals_ok for b in vals_bad if a < b)
    ties = sum(1 for a in vals_ok for b in vals_bad if a == b)
    auc = (pairs + 0.5 * ties) / (len(vals_ok) * len(vals_bad))
    tot = len(vals_ok) + len(vals_bad)
    acc, thr = max(((sum(1 for x in vals_ok if x < c) +
                     sum(1 for x in vals_bad if x >= c)) / tot, c)
                   for c in sorted(set(vals_ok + vals_bad)))
    return auc, thr, acc


print("=" * 78)
print("Thread 23 part 5 — secret-free proxy beta_hat for the Babai criterion")
print("=" * 78)

data = []
for label, cur, m in curves:
    for k1 in K1_GRID:
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, cur[2] - 1)
            r = instance(cur, m, d, k1, s)
            if r is not None:
                data.append((label, k1, s, r[0], r[1], r[2]))

print(f"\ncollected {len(data)} instances")

print("\n" + "=" * 78)
print("E11 — beta_hat (public) vs beta (secret) as predictors of recovery")
print("=" * 78)
for name, idx in (("beta      (secret-dependent)", 3),
                  ("beta_hat  (public)          ", 4)):
    ok = [d[idx] for d in data if d[5]]
    bad = [d[idx] for d in data if not d[5]]
    auc, thr, acc = auc_and_threshold(ok, bad)
    print(f"\n  {name}")
    print(f"    recovered : n={len(ok):3d} min={min(ok):.3f} "
          f"median={sorted(ok)[len(ok)//2]:.3f} max={max(ok):.3f}")
    print(f"    failed    : n={len(bad):3d} min={min(bad):.3f} "
          f"median={sorted(bad)[len(bad)//2]:.3f} max={max(bad):.3f}")
    print(f"    AUC={auc:.4f}   best threshold={thr:.4f} (accuracy {100*acc:.1f}%)")
    print(f"    certified-success region: x < {min(bad):.4f}  "
          f"({sum(1 for x in ok if x < min(bad))}/{len(ok)} of all successes)")
    print(f"    certified-failure region: x > {max(ok):.4f}  "
          f"({sum(1 for x in bad if x > max(ok))}/{len(bad)} of all failures)")

# how well does beta_hat track beta?
import statistics
b = [d[3] for d in data]; bh = [d[4] for d in data]
mb, mbh = statistics.mean(b), statistics.mean(bh)
cov = sum((x - mb) * (y - mbh) for x, y in zip(b, bh))
sb = math.sqrt(sum((x - mb) ** 2 for x in b))
sbh = math.sqrt(sum((y - mbh) ** 2 for y in bh))
print(f"\n  Pearson corr(beta, beta_hat) = {cov / (sb * sbh):.4f}")
print(f"  mean ratio beta/beta_hat     = "
      f"{statistics.mean(x / y for x, y in zip(b, bh) if y):.4f}")

print("\n  --- per-curve, per-K1 mean beta_hat vs success ---")
for label, cur, m in curves:
    print(f"\n  {label}  m={m}   lam*={lam_star(cur[3], cur[2]):.4f}")
    print("    K1        " + " ".join(f"{k:>7}" for k in K1_GRID))
    row_bh, row_w = [], []
    for k1 in K1_GRID:
        sub = [d for d in data if d[0] == label and d[1] == k1]
        row_bh.append(statistics.mean(d[4] for d in sub))
        row_w.append(sum(1 for d in sub if d[5]))
    print("    beta_hat  " + " ".join(f"{x:>7.3f}" for x in row_bh))
    print("    wins      " + " ".join(f"{str(w)+'/10':>7}" for w in row_w))

print("\n" + "=" * 78)
print("done")
print("=" * 78)
