#!/usr/bin/env python3
"""
GLV-HNP Thread 23 — closed form for lambda_1(L2), the quantity behind nu_hat.

Background (2026-07-29, commit e845207).  The Phase 2 lattice contains m copies
of a non-planted 2-dimensional sublattice

    L2 = < (n*S_K1, 0),  (-lam*S_K1, S_K2) >,      det L2 = n*S_K1*S_K2,

and nu_hat = lambda_1(L2)/sqrt(det L2) separates the June C1/C2 classes with
AUC 0.935 where eight curve-level invariants sat at the trivial baseline.
nu_hat was computed (Lagrange-Gauss) but not understood.

This script proves and machine-checks the closed form.  Every vector of L2 is

    w(a,b) = ( (a*n - b*lam)*S_K1,  b*S_K2 ),     a, b in Z,

so for fixed b the best a minimises |b*lam - a*n|.  Writing

    eps(b) = min_a |b*lam - a*n| = min(r, n-r),  r = (b*lam) mod n,

we get  ||w||^2 = S_K1^2*eps(b)^2 + S_K2^2*b^2, and

    lambda_1(L2)^2 = min( (n*S_K1)^2,  min_{b>=1} [ S_K1^2 eps(b)^2 + S_K2^2 b^2 ] ).

The objective is strictly increasing in both b and eps(b), so its minimiser is
Pareto-optimal for (b, eps(b)) -- i.e. b is a best approximation of the second
kind to lam/n -- and by the classical theorem those are exactly the continued
fraction convergent denominators q_j of lam/n.  Hence the CLAIM:

    lambda_1(L2)^2 = min over { 0 } u { q_j } of  S_K1^2 eps(q)^2 + S_K2^2 q^2 .

Tests
  T1  brute force over all b in [0,n) on small n  ==  CF formula  ==  exact Gauss
  T2  CF formula == exact Gauss at 20/24/64/128/256/384/521 bits
  T3  precision audit of the SHIPPED float gauss_reduce_2d against exact integers
  T4  scale localisation: which convergent wins, vs the balance scale
  T5  secp256k1 exact nu_hat, cross-checked against the values logged 2026-07-29
  T6  the two-term closed form and its AM-GM lower bound
  T7  root-choice sensitivity: nu_hat(lam) vs nu_hat(lam') for the two eigenvalues
  T8  the eff=1 exact zero is swap duality (lam' = lam^{-1} mod n)

Run:  python3 secp256k1_cm_audit/glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random

import sympy

# ---------------------------------------------------------------------------
# Scalings -- identical to glv_hnp_phase2_nuhat_base.py:scales()
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def k1k2(n, eff_num, eff_den):
    """K2 ~ sqrt(n), K1 ~ eff*sqrt(n); integer-only so it works at 521 bits."""
    s = math.isqrt(n)
    return max(2, (eff_num * s) // eff_den), s + 1

# ---------------------------------------------------------------------------
# Ground truth: exact-integer Lagrange-Gauss
# ---------------------------------------------------------------------------

def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss with exact integer rounding.  Returns lambda_1^2 (an int)."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) < n2(v):
        u, v = v, u
    while True:
        nv = n2(v)
        num = u[0] * v[0] + u[1] * v[1]
        # round-half-up of num/nv, exact for negative num too (Python // is floor)
        q = (2 * num + nv) // (2 * nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u)

def gauss_reduce_2d_float(u, v):
    """VERBATIM copy of the shipped float version (nuhat_base.py:216), for T3."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    while True:
        q = round((u[0] * v[0] + u[1] * v[1]) / nrm2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return math.sqrt(nrm2(u))

# ---------------------------------------------------------------------------
# The closed form
# ---------------------------------------------------------------------------

def convergent_denominators(lam, n):
    """Denominators q_j of the continued fraction convergents of lam/n."""
    a_list = []
    x, y = lam, n
    while y:
        a_list.append(x // y)
        x, y = y, x % y
    qm1, qm2 = 0, 1
    out = []
    for a in a_list:
        q = a * qm1 + qm2
        qm2, qm1 = qm1, q
        out.append(q)
    return out

def eps_centered(b, lam, n):
    """min_a |b*lam - a*n|."""
    r = (b * lam) % n
    return min(r, n - r)

def lambda1_sq_cf(n, lam, k1_bound, k2_bound, want_witness=False):
    """lambda_1(L2)^2 via the convergent enumeration.  Exact integer."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    best = (n * s1) ** 2               # the b = 0 vector
    wit = (0, n, 0)                    # (b, eps, index)
    qs = convergent_denominators(lam, n)
    for j, q in enumerate(qs):
        e = eps_centered(q, lam, n)
        val = s1 * s1 * e * e + s2 * s2 * q * q
        if val < best:
            best, wit = val, (q, e, j)
    return (best, wit, qs) if want_witness else best

def lambda1_sq_brute(n, lam, k1_bound, k2_bound):
    """Exhaustive min over every b in [0, n).  Only usable for tiny n."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    best = (n * s1) ** 2
    for b in range(1, n):
        e = eps_centered(b, lam, n)
        val = s1 * s1 * e * e + s2 * s2 * b * b
        if val < best:
            best = val
    return best

def nu_hat_exact_safe(n, lam, k1_bound, k2_bound):
    """Overflow-safe nu_hat for 256..521 bits: work with integer bit lengths."""
    s1, _, s2, _ = scales(n, k1_bound, k2_bound)
    l1sq = lambda1_sq_cf(n, lam, k1_bound, k2_bound)
    det = n * s1 * s2
    # log(l1sq/det) without float overflow
    ratio_log = _log_ratio(l1sq, det)
    return math.exp(0.5 * ratio_log)

def _log_ratio(a, b):
    """log(a/b) for arbitrarily large ints, without float overflow."""
    if a == 0:
        return float('-inf')
    sa, sb = a.bit_length(), b.bit_length()
    shift = max(sa, sb) - 900
    if shift > 0:
        a >>= shift
        b >>= shift
    return math.log(a) - math.log(b)

# ---------------------------------------------------------------------------
# Instance generation
# ---------------------------------------------------------------------------

def rand_instance(bits, rng, eff_num=5, eff_den=100):
    n = rng.randrange(1 << (bits - 1), 1 << bits) | 1
    lam = rng.randrange(1, n)
    k1, k2 = k1k2(n, eff_num, eff_den)
    return n, lam, k1, k2

# ---------------------------------------------------------------------------
# T1 -- brute force validation on small n
# ---------------------------------------------------------------------------

def t1_brute(rng, trials=400, nmax=6000):
    print("=" * 74)
    print("T1  brute force over all b  vs  CF closed form  vs  exact Gauss")
    print("=" * 74)
    bad_cf = bad_gauss = 0
    for _ in range(trials):
        n = rng.randrange(500, nmax) | 1
        lam = rng.randrange(1, n)
        k1, k2 = k1k2(n, 5, 100)
        s1, _, s2, _ = scales(n, k1, k2)
        b_val = lambda1_sq_brute(n, lam, k1, k2)
        c_val = lambda1_sq_cf(n, lam, k1, k2)
        g_val = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
        if c_val != b_val:
            bad_cf += 1
            if bad_cf <= 3:
                print(f"  CF MISMATCH n={n} lam={lam}: cf={c_val} brute={b_val}")
        if g_val != b_val:
            bad_gauss += 1
            if bad_gauss <= 3:
                print(f"  GAUSS MISMATCH n={n} lam={lam}: gauss={g_val} brute={b_val}")
    print(f"  trials={trials}  CF mismatches={bad_cf}  exact-Gauss mismatches={bad_gauss}")
    print(f"  verdict: CF closed form {'EXACT' if not bad_cf else 'WRONG'}; "
          f"exact Gauss {'EXACT' if not bad_gauss else 'WRONG'}")
    print()
    return bad_cf, bad_gauss

# ---------------------------------------------------------------------------
# T2 -- CF vs exact Gauss at cryptographic sizes
# ---------------------------------------------------------------------------

def t2_scale(rng, sizes=(20, 24, 64, 128, 256, 384, 521), trials=500):
    print("=" * 74)
    print("T2  CF closed form vs exact-integer Lagrange-Gauss, by bit length")
    print("=" * 74)
    print(f"  {'bits':>5} {'trials':>7} {'mismatch':>9}  verdict")
    total_bad = 0
    for bits in sizes:
        bad = 0
        for _ in range(trials):
            n, lam, k1, k2 = rand_instance(bits, rng)
            s1, _, s2, _ = scales(n, k1, k2)
            c = lambda1_sq_cf(n, lam, k1, k2)
            g = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
            if c != g:
                bad += 1
        total_bad += bad
        print(f"  {bits:>5} {trials:>7} {bad:>9}  {'EXACT' if not bad else 'MISMATCH'}")
    print(f"  total mismatches across all sizes: {total_bad}")
    print()
    return total_bad

# ---------------------------------------------------------------------------
# T3 -- precision audit of the shipped float reduction
# ---------------------------------------------------------------------------

def t3_precision(rng, sizes=(20, 24, 64, 128, 256, 384, 521), trials=400):
    print("=" * 74)
    print("T3  precision audit: SHIPPED float gauss_reduce_2d vs exact integers")
    print("=" * 74)
    print(f"  {'bits':>5} {'trials':>7} {'wrong':>6} {'err':>6} {'max rel dev':>12}  note")
    for bits in sizes:
        wrong = errs = 0
        maxdev = 0.0
        for _ in range(trials):
            n, lam, k1, k2 = rand_instance(bits, rng)
            s1, _, s2, _ = scales(n, k1, k2)
            exact_sq = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
            try:
                f = gauss_reduce_2d_float((n * s1, 0), (-lam * s1, s2))
            except (OverflowError, ZeroDivisionError, ValueError):
                errs += 1
                continue
            ex = math.exp(0.5 * _log_ratio(exact_sq, 1))
            if ex > 0:
                dev = abs(f - ex) / ex
                maxdev = max(maxdev, dev)
                if dev > 1e-9:
                    wrong += 1
        note = "ok" if not wrong and not errs else \
               ("raises" if errs else "loses the true lambda_1")
        print(f"  {bits:>5} {trials:>7} {wrong:>6} {errs:>6} {maxdev:>12.3e}  {note}")
    print()

# ---------------------------------------------------------------------------
# T4 -- scale localisation
# ---------------------------------------------------------------------------

def t4_scale_localisation(rng, bits=24, trials=3000, eff_list=((5, 100), (10, 100), (25, 100))):
    print("=" * 74)
    print("T4  which convergent wins?  q_j* vs the balance scale sqrt(n*S_K1/S_K2)")
    print("=" * 74)
    print("  balance scale q* solves S_K1*(n/q_{j+1}) = S_K2*q_j, i.e. q ~ sqrt(n*S_K1/S_K2)")
    print(f"  {'eff':>6} {'trials':>7} {'med q*/qbal':>12} {'IQR':>20} "
          f"{'med j*/J':>9} {'frac j*<J/4':>12} {'frac j*>3J/4':>13}")
    for en, ed in eff_list:
        ratios, fracs = [], []
        tail_lo = tail_hi = 0
        for _ in range(trials):
            n, lam, k1, k2 = rand_instance(bits, rng, en, ed)
            s1, _, s2, _ = scales(n, k1, k2)
            _, wit, qs = lambda1_sq_cf(n, lam, k1, k2, want_witness=True)
            q, _e, j = wit
            if q == 0:
                continue
            qbal = math.sqrt(n * s1 / s2)
            ratios.append(q / qbal)
            J = max(1, len(qs) - 1)
            fr = j / J
            fracs.append(fr)
            if fr < 0.25:
                tail_lo += 1
            if fr > 0.75:
                tail_hi += 1
        ratios.sort()
        fracs.sort()
        N = len(ratios)
        med = ratios[N // 2]
        q1, q3 = ratios[N // 4], ratios[(3 * N) // 4]
        print(f"  {en/ed:>6.2f} {N:>7} {med:>12.3f} [{q1:>8.3f},{q3:>8.3f}] "
              f"{fracs[N//2]:>9.3f} {tail_lo/N:>12.3f} {tail_hi/N:>13.3f}")
    print()
    print("  Reading: the winning convergent sits at a FIXED scale set by (n, K1, K2).")
    print("  Scale-free CF statistics (q_cf, max_q_cf, max_a over the whole expansion)")
    print("  average over all j and therefore cannot see it -- which is exactly why all")
    print("  six June CF/curve invariants sat at the trivial baseline.")
    print()

# ---------------------------------------------------------------------------
# T5 -- secp256k1, exact
# ---------------------------------------------------------------------------

SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

def secp_lambda():
    """The GLV eigenvalue: the root of x^2+x+1 = 0 mod n.  Returns min(r, n-1-r)."""
    # r = (-1 + sqrt(-3))/2 mod n; find sqrt(-3) by Tonelli-Shanks on -3 mod n
    n = SECP_N
    s = tonelli(n - 3, n)
    assert s is not None and (s * s - (n - 3)) % n == 0
    inv2 = pow(2, n - 2, n)
    r1 = ((n - 1) + s) % n * inv2 % n
    r2 = (n - 1 - r1) % n
    for r in (r1, r2):
        assert (r * r + r + 1) % n == 0, "not a primitive cube root of unity"
    return min(r1, r2), max(r1, r2)

def tonelli(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r

def t5_secp(rng):
    print("=" * 74)
    print("T5  secp256k1: exact nu_hat, cross-check against the 2026-07-29 log")
    print("=" * 74)
    lam, lam2 = secp_lambda()
    print(f"  n   = {SECP_N:#x}")
    print(f"  lam = {lam:#x}   (lam^2+lam+1 = 0 mod n: verified)")
    print(f"  mu  = min(lam,n-lam)/n = {min(lam, SECP_N-lam)/SECP_N:.6f}")
    print()
    logged = {(5, 100): 0.8709, (10, 100): 0.6624, (25, 100): 0.5852, (50, 100): 0.6851}
    print(f"  {'eff':>6} {'nu_hat exact':>13} {'logged 07-29':>13} {'delta':>10} "
          f"{'q_j*':>22} {'j*/J':>7}")
    for (en, ed), want in logged.items():
        k1, k2 = k1k2(SECP_N, en, ed)
        nh = nu_hat_exact_safe(SECP_N, lam, k1, k2)
        _, wit, qs = lambda1_sq_cf(SECP_N, lam, k1, k2, want_witness=True)
        q, _e, j = wit
        J = max(1, len(qs) - 1)
        print(f"  {en/ed:>6.2f} {nh:>13.4f} {want:>13.4f} {nh-want:>+10.4f} "
              f"{q:>22} {j/J:>7.3f}")
    print()
    print("  root-choice control (lam vs n-1-lam, must agree -- both are 'the' eigenvalue):")
    for (en, ed) in logged:
        k1, k2 = k1k2(SECP_N, en, ed)
        a = nu_hat_exact_safe(SECP_N, lam, k1, k2)
        b = nu_hat_exact_safe(SECP_N, lam2, k1, k2)
        print(f"    eff={en/ed:.2f}  nu_hat(lam)={a:.4f}  nu_hat(lam')={b:.4f}  "
              f"|d|={abs(a-b):.2e}")
    print()

# ---------------------------------------------------------------------------
# T6 -- the normalised two-term form and the AM-GM bound
# ---------------------------------------------------------------------------

def t6_closed_form(rng, bits=64, trials=2000):
    print("=" * 74)
    print("T6  normalised form  nu_hat^2 = min_j (x_j^2 + y_j^2),  x_j*y_j = q_j*eps_j/n")
    print("=" * 74)
    print("  with t = S_K1/S_K2,  x_j = q_j/sqrt(t*n),  y_j = eps_j*sqrt(t/n).")
    print("  AM-GM then gives the parameter-free lower bound")
    print("      nu_hat >= sqrt(2 * min_j q_j*eps_j/n)   and   q_j*eps_j/n ~ q_j/q_{j+1}.")
    print()
    bad = 0
    tight = []
    for _ in range(trials):
        n, lam, k1, k2 = rand_instance(bits, rng)
        s1, _, s2, _ = scales(n, k1, k2)
        l1sq, wit, qs = lambda1_sq_cf(n, lam, k1, k2, want_witness=True)
        det = n * s1 * s2
        nu = math.exp(0.5 * _log_ratio(l1sq, det))
        # normalised two-term form, recomputed independently
        t = s1 / s2
        best = None
        gmin = None
        for q in qs:
            e = eps_centered(q, lam, n)
            x = q / math.sqrt(t * n)
            y = e * math.sqrt(t / n)
            v = x * x + y * y
            best = v if best is None else min(best, v)
            # the last convergent has q = n, eps = 0 exactly; it is a real lattice
            # vector (0, n*S_K2) and stays in `best`, but it zeroes the AM-GM
            # product and must be excluded from the bound.
            if e > 0:
                g = q * e / n
                gmin = g if gmin is None else min(gmin, g)
        # b = 0 term in normalised coordinates: x = 0, y = n*sqrt(t/n)
        best = min(best, (n * math.sqrt(t / n)) ** 2)
        if best is None:
            continue
        if abs(math.sqrt(best) - nu) > 1e-6 * max(1.0, nu):
            bad += 1
        lb = math.sqrt(2 * gmin) if gmin is not None and gmin > 0 else 0.0
        if lb > 0:
            tight.append(nu / lb)
    tight.sort()
    N = len(tight)
    print(f"  trials={trials}  two-term-form mismatches={bad}  "
          f"{'CONSISTENT' if not bad else 'INCONSISTENT'}")
    if N:
        print(f"  nu_hat / AM-GM lower bound:  min={tight[0]:.4f}  "
              f"q25={tight[N//4]:.4f}  med={tight[N//2]:.4f}  "
              f"q75={tight[(3*N)//4]:.4f}  max={tight[-1]:.4f}")
        print(f"  bound violated (ratio < 1) in {sum(1 for r in tight if r < 1-1e-9)}/{N} cases")
    print()

# ---------------------------------------------------------------------------
# T7 -- root-choice sensitivity (checks the 2026-07-29 "<0.02" claim)
# ---------------------------------------------------------------------------

def t7_root_choice(rng, sizes=(20, 24, 64, 128, 256), trials=300, cut=0.645):
    print("=" * 74)
    print("T7  root choice: nu_hat(lam) vs nu_hat(lam') for the TWO GLV eigenvalues")
    print("=" * 74)
    print("  For j=0 curves the eigenvalues are the two roots of x^2+x+1 mod n,")
    print("  i.e. lam and lam' = n-1-lam.  nu_hat is exactly invariant under")
    print("  lam -> n-lam (eps_centered is), but that is the WRONG involution.")
    print()
    print(f"  {'bits':>5} {'curves':>7} {'med |d|':>9} {'q90 |d|':>9} {'max |d|':>9} "
          f"{'frac>0.02':>10} {'cut flips':>10}")
    for bits in sizes:
        ds, flips, mirror_bad = [], 0, 0
        got = 0
        p = rng.randrange(1 << (bits - 1), 1 << bits)
        while got < trials:
            p = int(sympy.nextprime(p))
            if p.bit_length() > bits:
                p = rng.randrange(1 << (bits - 1), 1 << bits)
                continue
            if p % 3 != 1:
                continue
            s = tonelli(p - 3, p)
            if s is None:
                continue
            inv2 = pow(2, p - 2, p)
            lam = ((p - 1) + s) % p * inv2 % p
            if (lam * lam + lam + 1) % p != 0:
                continue
            lam2 = (p - 1 - lam) % p
            k1, k2 = k1k2(p, 10, 100)
            a = nu_hat_exact_safe(p, lam, k1, k2)
            b = nu_hat_exact_safe(p, lam2, k1, k2)
            # sanity: the true symmetry lam -> n-lam must be exact
            c = nu_hat_exact_safe(p, (p - lam) % p, k1, k2)
            if abs(a - c) > 1e-12:
                mirror_bad += 1
            ds.append(abs(a - b))
            if (a <= cut) != (b <= cut):
                flips += 1
            got += 1
        ds.sort()
        N = len(ds)
        print(f"  {bits:>5} {N:>7} {ds[N//2]:>9.4f} {ds[(9*N)//10]:>9.4f} "
              f"{ds[-1]:>9.4f} {sum(1 for d in ds if d > 0.02)/N:>10.3f} "
              f"{flips/N:>10.3f}"
              + ("   MIRROR-BROKEN" if mirror_bad else ""))
    print()
    # eff sweep at fixed bit size: the split grows with eff
    print(f"  eff sweep at 256 bits ({trials} curves each):")
    print(f"  {'eff':>6} {'med |d|':>9} {'q90 |d|':>9} {'max |d|':>9} {'cut flips':>10}")
    p = rng.randrange(1 << 255, 1 << 256)
    curves = []
    while len(curves) < trials:
        p = int(sympy.nextprime(p))
        if p % 3 != 1:
            continue
        s = tonelli(p - 3, p)
        if s is None:
            continue
        inv2 = pow(2, p - 2, p)
        lam = ((p - 1) + s) % p * inv2 % p
        if (lam * lam + lam + 1) % p != 0:
            continue
        curves.append((p, lam, (p - 1 - lam) % p))
    for en, ed in ((5, 100), (10, 100), (25, 100), (50, 100), (100, 100)):
        ds, flips = [], 0
        for (nn, l1, l2) in curves:
            k1, k2 = k1k2(nn, en, ed)
            a = nu_hat_exact_safe(nn, l1, k1, k2)
            b = nu_hat_exact_safe(nn, l2, k1, k2)
            ds.append(abs(a - b))
            if (a <= cut) != (b <= cut):
                flips += 1
        ds.sort()
        N = len(ds)
        print(f"  {en/ed:>6.2f} {ds[N//2]:>9.4f} {ds[(9*N)//10]:>9.4f} {ds[-1]:>9.4f} "
              f"{flips/N:>10.3f}")
    print()
    print("  Consequence: nu_hat is a function of the ATTACK INSTANCE (n, lam, K1, K2),")
    print("  not of the curve -- the two eigenvalues of one curve give different nu_hat")
    print(f"  and land on opposite sides of the nu_hat={cut} C1/C2 cut a nonzero")
    print("  fraction of the time.  It also shows mu = min(lam, n-lam)/n symmetrises")
    print("  over the wrong involution: nu_hat is already exactly invariant under")
    print("  lam -> n-lam, so mu cannot be the coordinate that distinguishes the roots.")
    print()

# ---------------------------------------------------------------------------
# T8 -- why the eff=1 row of T7 is an exact zero: swap duality
# ---------------------------------------------------------------------------

def t8_swap_duality(rng, bits=256, trials=200):
    print("=" * 74)
    print("T8  the eff=1 exact zero in T7 is a genuine isometry, not a bug")
    print("=" * 74)
    print("  For a primitive cube root of unity, lam^2+lam+1=0 mod n gives lam^3=1,")
    print("  hence lam' = n-1-lam = lam^2 = lam^{-1} mod n.  Writing")
    print("      L_lam = { (a*n - b*lam)*S_K1, b*S_K2 ) },")
    print("  the coordinate swap (x,y) -> (y,x) carries L_lam to L_{lam^{-1}} = L_{lam'}.")
    print("  That swap is an isometry exactly when S_K1 = S_K2, giving")
    print("      nu_hat(lam) = nu_hat(lam')   and   b*(lam) = eps*(lam'),  eps*(lam) = b*(lam').")
    print()
    print("  NOTE: the codebase's eff=1 setting is K1 = isqrt(n), K2 = isqrt(n)+1, which")
    print("  gives S_K1 != S_K2 -- only S_K1/S_K2 = 1 + O(n^-1/2).  So the T7 eff=1.00 row")
    print("  is a NEAR-isometry (agreement to O(n^-1/2), i.e. below print precision), not")
    print("  an exact one.  Arm A forces S_K1 = S_K2 to test the exact statement; Arm B")
    print("  measures the codebase parameterisation.")
    print()
    p = rng.randrange(1 << (bits - 1), 1 << bits)
    curves = []
    ok_inv = 0
    while len(curves) < trials:
        p = int(sympy.nextprime(p))
        if p % 3 != 1:
            continue
        s = tonelli(p - 3, p)
        if s is None:
            continue
        lam = ((p - 1) + s) % p * pow(2, p - 2, p) % p
        if (lam * lam + lam + 1) % p != 0:
            continue
        lam2 = (p - 1 - lam) % p
        if pow(lam, -1, p) == lam2:
            ok_inv += 1
        curves.append((p, lam, lam2))
    print(f"  curves={len(curves)}   lam' == lam^-1 mod n: {ok_inv}/{len(curves)}")
    print()

    # Arm A: force S_K1 == S_K2 exactly by using the same bound twice.
    ok_nu = ok_wit = 0
    for (nn, l1, l2) in curves:
        kk = math.isqrt(nn)
        s1, _, s2, _ = scales(nn, kk, kk)
        assert s1 == s2
        a = lambda1_sq_cf(nn, l1, kk, kk)
        b = lambda1_sq_cf(nn, l2, kk, kk)
        if a == b:
            ok_nu += 1
        _, wa, _ = lambda1_sq_cf(nn, l1, kk, kk, want_witness=True)
        _, wb, _ = lambda1_sq_cf(nn, l2, kk, kk, want_witness=True)
        if wa[0] == wb[1] and wa[1] == wb[0]:
            ok_wit += 1
    print(f"  Arm A (S_K1 = S_K2 forced):")
    print(f"    lambda_1^2(lam) == lambda_1^2(lam') exactly:  {ok_nu}/{len(curves)}")
    print(f"    witness pair exchanged (b* <-> eps*):         {ok_wit}/{len(curves)}")
    print()

    # Arm B: the codebase's eff=1 parameterisation.
    ok_nu_b = ok_wit_b = 0
    worst = 0.0
    ratios = []
    for (nn, l1, l2) in curves:
        k1, k2 = k1k2(nn, 100, 100)
        s1, _, s2, _ = scales(nn, k1, k2)
        ratios.append(s1 / s2)
        a = lambda1_sq_cf(nn, l1, k1, k2)
        b = lambda1_sq_cf(nn, l2, k1, k2)
        if a == b:
            ok_nu_b += 1
        _, wa, _ = lambda1_sq_cf(nn, l1, k1, k2, want_witness=True)
        _, wb, _ = lambda1_sq_cf(nn, l2, k1, k2, want_witness=True)
        if wa[0] == wb[1] and wa[1] == wb[0]:
            ok_wit_b += 1
        worst = max(worst, abs(nu_hat_exact_safe(nn, l1, k1, k2)
                               - nu_hat_exact_safe(nn, l2, k1, k2)))
    print(f"  Arm B (codebase eff=1: K1=isqrt(n), K2=isqrt(n)+1):")
    print(f"    S_K1/S_K2 - 1:  max {max(ratios)-1:.3e}  min {min(ratios)-1:.3e}")
    print(f"    lambda_1^2 equal exactly:                     {ok_nu_b}/{len(curves)}")
    print(f"    witness pair exchanged (b* <-> eps*):         {ok_wit_b}/{len(curves)}")
    print(f"    max |nu_hat(lam) - nu_hat(lam')|:             {worst:.3e}")
    print()
    print("  So the duality is structural (witnesses exchange even in Arm B), while the")
    print("  numerical equality of nu_hat is exact only when S_K1 = S_K2 and otherwise")
    print("  holds to O(n^-1/2).  Away from eff=1 the swap rescales by t = S_K1/S_K2 and")
    print("  |d| in T7 grows with |t-1|.")
    print()

# ---------------------------------------------------------------------------

def main():
    rng = random.Random(20260802)
    print()
    print("GLV-HNP Thread 23 — closed form for lambda_1(L2) behind nu_hat")
    print("seed=20260802")
    print()
    t1_brute(rng)
    t2_scale(rng)
    t3_precision(rng)
    t4_scale_localisation(rng)
    t5_secp(rng)
    t6_closed_form(rng)
    t7_root_choice(rng)
    t8_swap_duality(rng)
    print("done.")

if __name__ == "__main__":
    main()
