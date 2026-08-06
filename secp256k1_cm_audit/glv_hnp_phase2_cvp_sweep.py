"""
GLV-HNP Phase 2, Thread 23 part 2: does target-centring generalise?

The 2x2 design in `glv_hnp_phase2_cvp.py` was run on the three historical
curves only.  This script re-runs the 2026-07-29 experiment T3 sweep --
20 fresh j=0 GLV curves with p in [2^16, 2^17), lam* spread over (0, 0.5)
in 10 bins -- with the centred methods added, so the effect is measured
across curves rather than on the two curves that motivated it.

T3 baseline to beat (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, uncentred
KAN-LLL, m = 12, 5 seeds, "curves recovering 5/5"):

    eff = K1*K2/n = 0.05  ->  19/20
    eff = 0.15            ->   3/20
    eff = 0.25            ->   0/20

The counting bound [IT] for these curves is eff < n^(-1/m) = 0.383 at
m = 12, so eff = 0.25 is still well inside the feasible region: the T3
zeroes were a lattice failure, not an information failure.

Curve-search helpers are verbatim from
glv_hnp_phase2_lambda_threshold.py:146-414; the lattice methods are
imported from glv_hnp_phase2_cvp.py.

Run: python3 glv_hnp_phase2_cvp_sweep.py
"""

import math
import random

import sympy

from glv_hnp_phase2_cvp import (
    ec_mul, tonelli_shanks, modinv, find_generator, lam_star,
    rates, it_bound_k1, fmt,
)

# ---------------------------------------------------------------------------
# Eisenstein / j=0 curve search -- verbatim from
# glv_hnp_phase2_lambda_threshold.py:146
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p)."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n.  Requires n = 1 mod 3."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def build_curve(p, n, seed=12345):
    """Find the sextic twist b with #E = n, plus a generator."""
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by
    lam* into `nbins` bins over [0, 0.5]."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    ls = lam_star(lam, n_cand)
                    idx = min(nbins - 1, int(ls / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 12
LO, HI = 2 ** 16, 2 ** 17
SWEEP_METHODS = ("KAN-LLL", "KAN-LLL-C", "BAB-LLL-C")
EFFS = (0.05, 0.15, 0.25, 0.32, 0.38, 0.45)

print("=" * 92)
print("Thread 23 part 2 -- does centring generalise across curves?  (T3 re-run)")
print("=" * 92)
print(f"\nSearching j=0 GLV curves with p in [{LO}, {HI}) ...")
curves = search_curves(LO, HI, per_bin=2, nbins=10)
print(f"Found {len(curves)} curves.")

n_med = sorted(c[2] for c in curves)[len(curves) // 2]
print(f"median n = {n_med}   [IT] eff bound at m={M_SIGS} is "
      f"n^(-1/m) = {n_med ** (-1.0 / M_SIGS):.3f}\n")

summary = {}
for eff_t in EFFS:
    print("\n" + "-" * 92)
    print(f"=== eff target = {eff_t}   (m = {M_SIGS}, {len(SEEDS)} seeds) ===")
    print("-" * 92)
    print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>5} {'eff':>6} {'ITK1':>6} "
          + " ".join(f"{k:>9}" for k in SWEEP_METHODS))
    full = {k: 0 for k in SWEEP_METHODS}
    anyw = {k: 0 for k in SWEEP_METHODS}
    tot_curves = 0
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        acc, tot, _ = rates((p, b, n, lam, G), M_SIGS, k1, SEEDS,
                            methods=SWEEP_METHODS)
        if tot == 0:
            continue
        tot_curves += 1
        for k in SWEEP_METHODS:
            full[k] += (acc[k] == tot)
            anyw[k] += (acc[k] > 0)
        print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1:>5} "
              f"{k1 * k2 / n:>6.3f} {it_bound_k1(n, M_SIGS, k2):>6.1f} "
              + fmt(acc, tot, SWEEP_METHODS))
    summary[eff_t] = (full, anyw, tot_curves)
    print(f"{'FULL 5/5':>34} "
          + " ".join(f"{str(full[k]) + '/' + str(tot_curves):>9}"
                     for k in SWEEP_METHODS))
    print(f"{'ANY >=1/5':>34} "
          + " ".join(f"{str(anyw[k]) + '/' + str(tot_curves):>9}"
                     for k in SWEEP_METHODS))

print("\n" + "=" * 92)
print("SUMMARY -- curves recovering 5/5 out of all curves, by eff")
print("=" * 92)
print(f"{'eff':>6} {'IT feasible?':>13} " + " ".join(f"{k:>11}" for k in SWEEP_METHODS))
for eff_t in EFFS:
    full, _anyw, tc = summary[eff_t]
    feas = "yes" if eff_t < n_med ** (-1.0 / M_SIGS) else "NO (over [IT])"
    print(f"{eff_t:>6.2f} {feas:>13} "
          + " ".join(f"{str(full[k]) + '/' + str(tc):>11}" for k in SWEEP_METHODS))

print("\nT3 baseline (2026-07-29, uncentred KAN-LLL): 0.05 -> 19/20, "
      "0.15 -> 3/20, 0.25 -> 0/20")
print("=" * 92)
