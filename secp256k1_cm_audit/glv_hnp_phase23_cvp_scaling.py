"""
GLV-HNP Thread 23, part 2: does the CVP wall shift replicate on fresh curves?

glv_hnp_phase23_cvp.py showed, on the two historical 12-bit curves, that
replacing the Kannan/SVP embedding with an explicit CVP on L23 (dim 2m+1)
moves the K1 wall from 4-6 out to 12.  Two curves is not a result.  This
script replicates on fresh curves at two bit-lengths, with lam* spread over
(0, 0.5), and reports the wall in terms of the bias strength

    eff = K1 * K2 / n           (K2 = isqrt(n)+1, so eff ~ K1/sqrt(n))

which the 2026-07-29 run identified as the real controlling variable.

Only babai_c (Babai nearest-plane on the centered target) is run against the
Kannan baseline: part 1 found exact CVP and Babai agree on the wall location
at every grid point but one, and Babai is O(dim^2) after LLL rather than an
enumeration, so it is what would actually scale.

Run: python3 glv_hnp_phase23_cvp_scaling.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, GSO

from glv_hnp_phase23_cvp_lib import (
    ec_mul, tonelli_shanks, find_generator, gen_signatures,
    attack_kannan, attack_cvp, lam_star,
)

# ---------------------------------------------------------------------------
# Curve search (adapted from glv_hnp_phase2_lambda_threshold.py:374)
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

def search_curves(lo, hi, per_bin=1, nbins=6, max_primes=200000):
    """j=0 GLV curves with n prime, bucketed by lam* over [0, 0.5]."""
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

# ---------------------------------------------------------------------------
# Wall measurement
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]
# Grid extended past 48 because the first run right-censored 5 of 6 20-bit
# CVP walls at the old grid maximum, which understates the shift.
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192]

def trials(curve, k1, m, method):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    wins = tot = 0
    for seed in SEEDS:
        d_secret = random.Random(seed ^ 0xABCDEF).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1, k2, seed)
        if len(sigs) < m:
            continue
        tot += 1
        if method == 'kannan':
            wins += bool(attack_kannan(sigs, n, lam, k1, k2, d_secret))
        else:
            ok, _, _, _ = attack_cvp(sigs, n, lam, k1, k2, d_secret,
                                     centered=True, exact=False)
            wins += bool(ok)
    return wins, tot

def wall(curve, m, method):
    """Largest K1 in the grid at which the method scores 5/5."""
    best = None
    for k1 in K1_GRID:
        w, t = trials(curve, k1, m, method)
        if t and w == t:
            best = k1
        elif best is not None and k1 >= 2 * best:
            break          # two grid steps past the last success: stop
    return best

print("=" * 78)
print("Thread 23 part 2 - does the CVP wall shift replicate on fresh curves?")
print("=" * 78)

for lo, hi, m, tag in ((2**16, 2**17, 10, "17-bit"),
                       (2**19, 2**20, 12, "20-bit")):
    print(f"\n{'-'*78}\n{tag} curves, m={m}, seeds={len(SEEDS)}\n{'-'*78}")
    found = search_curves(lo, hi, per_bin=1, nbins=6)
    hdr = (f"{'p':>8s} {'n':>8s} {'lam*':>7s} {'K2':>5s} | "
           f"{'K1 wall kannan':>14s} {'K1 wall babai_c':>15s} "
           f"{'shift':>6s} | {'eff kan':>8s} {'eff cvp':>8s}")
    print(hdr)
    print("-" * len(hdr))
    shifts = []
    for (p, b, n, lam, G) in found:
        curve = (p, b, n, lam, G)
        k2 = math.isqrt(n) + 1
        wk = wall(curve, m, 'kannan')
        wc = wall(curve, m, 'babai_c')
        if wk and wc:
            shifts.append(wc / wk)
        sh = f"{wc/wk:.2f}x" if (wk and wc) else "-"
        ek = f"{wk*k2/n:.4f}" if wk else "-"
        ec = f"{wc*k2/n:.4f}" if wc else "-"
        print(f"{p:8d} {n:8d} {lam_star(lam,n):7.4f} {k2:5d} | "
              f"{str(wk):>14s} {str(wc):>15s} {sh:>6s} | {ek:>8s} {ec:>8s}")
    if shifts:
        print(f"\n  {tag}: CVP wall / Kannan wall = "
              f"{min(shifts):.2f}x .. {max(shifts):.2f}x "
              f"(mean {sum(shifts)/len(shifts):.2f}x) over {len(shifts)} curves")
        print(f"  replication: {sum(1 for s in shifts if s > 1)}/{len(shifts)} "
              f"curves show a strictly wider CVP wall")

print("\nDone.")
