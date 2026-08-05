"""
GLV-HNP Thread 23, part B — sharpening the U2/U3 boundary.

Follows glv_hnp_phase2_thread23.py.  Three additions:

  U4  fine K1 grid, 10 seeds, on the hard curve (lam*=0.07) — locate the wall
      for V0 (baseline Kannan) vs V2_exact (exact CVP) precisely.
  U5  per-seed contingency of  [true error is the closest vector]  against
      [exact CVP recovers d].  Tests the claim: recovery fails exactly when
      the planted error stops being the closest lattice point.
  U6  replication on fresh 13-16 bit j=0 GLV curves, so the U2 conclusion does
      not rest on the two historical curves.

Run: python3 glv_hnp_phase2_thread23b.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, CVP

from glv_hnp_phase2_thread23 import (
    ec_mul, tonelli_shanks, find_generator, modinv, lam_star, norm, scales,
    gen_signatures, build_V0, planted_V0, recover_V0,
    build_V2, target_V2, true_error_V2, d_from_error, trial,
)

# ---------------------------------------------------------------------------
# CM curve search (verbatim from glv_hnp_phase2_lambda_threshold.py:146-178,356)
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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
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

def fresh_curves(lo, hi, want, lam_star_lo=0.0, lam_star_hi=0.5):
    """Collect `want` j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = p + 1 - t
                    if nc < 8 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    lam = roots[0]
                    ls = lam_star(lam, nc)
                    if not (lam_star_lo <= ls <= lam_star_hi):
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    out.append((f"{nc.bit_length()}b/{p}", (p, cur[1], nc, lam, cur[4])))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Per-seed record
# ---------------------------------------------------------------------------

def per_seed(curve, m, k1_bound, seeds):
    rows = []
    p, b, n, lam, G = curve
    for s in seeds:
        ds = random.Random(s).randrange(1, n)
        r = trial(curve, m, ds, k1_bound, s)
        if r is None:
            continue
        rows.append(r)
    return rows


# ===========================================================================
if __name__ == '__main__':
    print("=" * 78)
    print("Thread 23b — sharpening the wall")
    print("=" * 78)

    SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 20260805, 88, 4242, 123456]

    HARD = ("12-bit/2677", (2677, 2, 2647, 185, None), 10)
    EASY = ("12-bit/2557", (2557, 2, 2659, 1755, None), 8)
    hist = []
    for label, (p, b, n, lam, _), m in (EASY, HARD):
        G = find_generator(p, b, n)
        assert (lam * lam + lam + 1) % n == 0
        hist.append((label, (p, b, n, lam, G), m))

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U4: fine K1 grid, 10 seeds — where exactly is the wall?")
    print("-" * 78)
    FINE = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
    for label, curve, m in hist:
        p, b, n, lam, G = curve
        print(f"\n{label}: n={n}, lam*={lam_star(lam,n):.4f}, m={m}, "
              f"eff=K1*K2/n at K1=1 is {(math.isqrt(n)+1)/n:.4f}")
        print(f"  {'solver':<14}" + "".join(f"{('K1=%d'%k):>8}" for k in FINE))
        tab = {k: [] for k in ('V0', 'V1', 'V2_babai', 'V2_exact', 'closest')}
        for k1b in FINE:
            rs = per_seed(curve, m, k1b, SEEDS10)
            t = len(rs)
            for k in ('V0', 'V1', 'V2_babai', 'V2_exact'):
                tab[k].append(f"{sum(1 for r in rs if r.get(k))}/{t}")
            tab['closest'].append(
                f"{sum(1 for r in rs if r.get('V2_true_is_closest'))}/{t}")
        for k in ('V0', 'V1', 'V2_babai', 'V2_exact', 'closest'):
            nm = 'true=closest' if k == 'closest' else k
            print(f"  {nm:<14}" + "".join(f"{v:>8}" for v in tab[k]))

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U5: contingency  [true error is closest]  vs  [exact CVP wins]")
    print("        pooled over both curves, K1 in 2..12, 10 seeds")
    print("-" * 78)
    cont = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    detail = []
    for label, curve, m in hist:
        for k1b in FINE:
            for r in per_seed(curve, m, k1b, SEEDS10):
                c = bool(r.get('V2_true_is_closest'))
                w = bool(r.get('V2_exact'))
                cont[(c, w)] += 1
                if c and not w:
                    detail.append((label, k1b, r['seed'], r['V2_etrue'],
                                   r['V2_exact_err']))
    tot = sum(cont.values())
    print(f"  total trials: {tot}")
    print(f"  {'':<22}{'CVP wins':>12}{'CVP fails':>12}")
    print(f"  {'true IS closest':<22}{cont[(True,True)]:>12}{cont[(True,False)]:>12}")
    print(f"  {'true is NOT closest':<22}{cont[(False,True)]:>12}{cont[(False,False)]:>12}")
    agree = cont[(True, True)] + cont[(False, False)]
    print(f"  agreement: {agree}/{tot} = {agree/tot:.3f}")
    print(f"  false negatives (closest but no d): {len(detail)}")
    for d in detail[:8]:
        print(f"    {d[0]} K1={d[1]} seed={d[2]} etrue={d[3]:.4e} cvp={d[4]:.4e}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U6: replication on fresh curves (does V2_exact beat V0 in general?)")
    print("-" * 78)
    fresh = fresh_curves(2 ** 13, 2 ** 14, 4)
    fresh += fresh_curves(2 ** 15, 2 ** 16, 2)
    SEEDS5 = SEEDS10[:5]
    print(f"  {'curve':<14}{'n':>8}{'lam*':>8}{'m':>4}"
          + "".join(f"{('K1=%d'%k):>10}" for k in (4, 6, 8, 12)))
    for label, curve in fresh:
        p, b, n, lam, G = curve
        m = 10
        cells = []
        for k1b in (4, 6, 8, 12):
            rs = per_seed(curve, m, k1b, SEEDS5)
            t = len(rs)
            v0 = sum(1 for r in rs if r.get('V0'))
            ve = sum(1 for r in rs if r.get('V2_exact'))
            cells.append(f"{v0}/{t}|{ve}/{t}")
        print(f"  {label:<14}{n:>8}{lam_star(lam,n):>8.4f}{m:>4}"
              + "".join(f"{c:>10}" for c in cells))
    print("  (cells are V0|V2_exact)")

    print("\n" + "=" * 78)
    print("done")
