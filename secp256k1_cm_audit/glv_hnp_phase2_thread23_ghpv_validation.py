"""
GLV-HNP Phase 2, Thread 23b: validate the GH/pv sufficient condition.

Thread 23 (glv_hnp_phase2_thread23_reformulation.py, run 2026-08-03) found
that in the d-eliminated lattice V3 (no Kannan row, no d column, dim 2m) the
ratio

    GH/pv  =  Gaussian-heuristic lambda_1  /  E||planted vector||

separated success from failure with ZERO false positives on the three
historical curves (6/6 rows with GH/pv > 1 recovered d at 5/5), while being
conservative in the other direction (LLL also succeeded on several rows with
GH/pv < 1).  Closed form, with S_K1 = n/K1, S_K2 = n/K2, det = n^(m-1) *
S_K1^m * S_K2^m and dim = 2m:

    GH/pv  =  sqrt(3 / (2*pi*e)) * eff^(-1/2) * n^(-1/(2m)),   eff = K1*K2/n

so the predicted sufficient condition for recovery is

    (*)   eff * n^(1/m)  <  3/(2*pi*e)  =  0.17567...

This script tests (*) on fresh curves.  The claim under test is one-sided:

    H23b:  eff * n^(1/m) < 0.1757  =>  the V2/V3 attack recovers d.

A single fresh curve/K1 row satisfying (*) that fails to recover d falsifies
H23b.  Rows violating (*) that still succeed do NOT falsify it (the condition
is sufficient, not necessary) but are counted, since a large excess means the
bound is loose.

Run: python3 glv_hnp_phase2_thread23_ghpv_validation.py
"""

import math
import random
import sympy

from glv_hnp_phase2_thread23_reformulation import (
    ec_mul, tonelli_shanks, find_generator, lam_star, norm, scales,
    gen_signatures, planted_v2, attack_v2, attack_v3,
)

GH_CONST = 3.0 / (2.0 * math.pi * math.e)          # 0.175666...

# ---------------------------------------------------------------------------
# Curve search (verbatim logic from glv_hnp_phase2_lambda_threshold.py:373)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
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
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = pow(2, -1, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=7):
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
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

# ---------------------------------------------------------------------------

def gh_over_pv_closed(n, m, k1_bound, k2_bound):
    """sqrt(3/(2 pi e)) * eff^(-1/2) * n^(-1/(2m)), the closed form of (*)."""
    eff = k1_bound * k2_bound / n
    return math.sqrt(GH_CONST / eff) * n ** (-1.0 / (2 * m))

def gh_over_pv_exact(curve, m, k1_bound, seeds):
    """The same ratio computed from the actual lattice determinant and the
    actual (not expected) planted-vector norms — checks the closed form."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    dim = 2 * m
    logdet = (m - 1) * math.log(n) + m * math.log(S_K1) + m * math.log(S_K2)
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)
    pvs = []
    for s in seeds:
        d_s = random.Random(s + 7777).randint(1, n - 1)
        sg = gen_signatures(G, d_s, m, n, lam, p, k1_bound, k2_bound, s)
        if len(sg) == m:
            pvs.append(norm(planted_v2(sg, n, k1_bound, k2_bound, kannan=False)))
    return gh / (sum(pvs) / len(pvs)) if pvs else float('nan')

def run_rate(curve, m, k1_bound, seeds, variant='V2'):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    wins = 0
    for s in seeds:
        d_s = random.Random(s + 7777).randint(1, n - 1)
        sg = gen_signatures(G, d_s, m, n, lam, p, k1_bound, k2_bound, s)
        if len(sg) < m:
            continue
        if variant == 'V2':
            ok = attack_v2(sg, n, lam, k1_bound, k2_bound, d_s)[0]
        else:
            ok = attack_v3(sg, n, lam, k1_bound, k2_bound, d_s)[0]
        wins += bool(ok)
    return wins, len(seeds)

# ===========================================================================

print("=" * 84)
print("Thread 23b — validating  eff * n^(1/m) < 3/(2*pi*e) = 0.17567  as a")
print("             SUFFICIENT condition for GLV-HNP Phase 2 recovery")
print("=" * 84)

SEEDS = [42, 1234, 9999, 555, 31337]
M_FIXED = 8
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]

print("\nSearching fresh j=0 GLV curves in [2^13, 2^15), lam* bucketed into 7 bins...")
curves = search_curves(2 ** 13, 2 ** 15, per_bin=2, nbins=7)
print(f"  found {len(curves)} curves\n")

print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>3} {'eff':>7} {'eff*n^(1/m)':>12} "
      f"{'GH/pv cf':>9} {'GH/pv ex':>9} {'V2':>5} {'V3':>5} {'(*)':>4} {'verdict':>9}")
print("-" * 84)

rows = []
for (p, b, n, lam, G) in curves:
    curve = (p, b, n, lam, G)
    k2_bound = math.isqrt(n) + 1
    for k1 in K1_GRID:
        eff = k1 * k2_bound / n
        if eff >= 2.0:
            continue
        stat = eff * n ** (1.0 / M_FIXED)
        cf = gh_over_pv_closed(n, M_FIXED, k1, k2_bound)
        ex = gh_over_pv_exact(curve, M_FIXED, k1, SEEDS)
        w2, t2 = run_rate(curve, M_FIXED, k1, SEEDS, 'V2')
        w3, t3 = run_rate(curve, M_FIXED, k1, SEEDS, 'V3')
        sat = stat < GH_CONST
        full = (w2 == t2)
        if sat and not full:
            verdict = "FALSIFY"
        elif sat and full:
            verdict = "ok"
        elif (not sat) and full:
            verdict = "loose"
        else:
            verdict = "-"
        rows.append((p, n, lam, k1, eff, stat, cf, ex, w2, t2, w3, sat, full, verdict))
        print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1:>3} {eff:>7.4f} "
              f"{stat:>12.4f} {cf:>9.3f} {ex:>9.3f} {str(w2)+'/'+str(t2):>5} "
              f"{str(w3)+'/'+str(t3):>5} {'Y' if sat else 'n':>4} {verdict:>9}")

print("-" * 84)
n_sat = sum(1 for r in rows if r[11])
n_fal = sum(1 for r in rows if r[13] == "FALSIFY")
n_loose = sum(1 for r in rows if r[13] == "loose")
print(f"\nRows total          : {len(rows)}")
print(f"Rows satisfying (*) : {n_sat}")
print(f"  of which 5/5      : {n_sat - n_fal}")
print(f"  of which FALSIFY  : {n_fal}")
print(f"Rows violating (*) that still recovered 5/5 (bound is loose): {n_loose}")

if rows:
    succ = [r[5] for r in rows if r[12]]
    fail = [r[5] for r in rows if not r[12]]
    if succ and fail:
        print(f"\nempirical eff*n^(1/m):  max over 5/5 rows = {max(succ):.4f}")
        print(f"                        min over non-5/5   = {min(fail):.4f}")
        print(f"predicted threshold  :  {GH_CONST:.5f}")

print("\nVERDICT: " + ("H23b FALSIFIED" if n_fal else
                       "H23b SURVIVES — no counterexample among rows satisfying (*)"))
print("=" * 84)
