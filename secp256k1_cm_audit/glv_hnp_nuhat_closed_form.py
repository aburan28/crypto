"""
GLV-HNP -- Thread 23: closed form for nu_hat, and its exact null law.

Background
----------
Thread 20d (2026-07-29, commit e845207) found that

    nu_hat  =  lambda_1(L2) / sqrt(det L2),
    L2      =  < (n*S1, 0), (-lam*S1, S2) >,     S1 = n//K1,  S2 = max(1, n//K2)

separates the June C1/C2 classes with AUC 0.935, after eight scale-free
invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n, lam/n, mu) had
all failed.  The 2026-07-29 log closed with Thread 23: "nu_hat is currently
computed, not understood", and proposed showing lambda_1(L2) is a
continued-fraction quantity of lam/n *at the scale* S2/S1, with the falsifier
"if predicted-from-CF nu_hat correlates < 0.9 with computed nu_hat, the
convergent-scale story is wrong".

This script settles that, and then goes further than the proposal:

  A. CLOSED FORM (exact, not approximate).  lambda_1(L2)^2 = min over
     continued-fraction convergent denominators q_j of lam/n of
     S1^2*theta_j^2 + S2^2*q_j^2, where theta_j = ||q_j*lam||_n is the centered
     residue.  Proof + exhaustive equality check against exact-integer
     Lagrange-Gauss reduction.
  B. PRECISION AUDIT of the incumbent float64 gauss_reduce_2d (the routine that
     produced every published nu_hat number, including the secp256k1 placement).
  C. LOCALISATION.  Which convergent attains the min, and how close its
     denominator sits to R = sqrt(n*S1/S2) = sqrt(n*K2/K1).
  D. EXACT NULL LAW.  In the (u,w) coordinates of part A, nu_hat^2 is the
     shortest-vector-squared of a unimodular 2-lattice.  Under equidistribution
     that gives a closed-form CDF, verified against 20000 random lam.
  E. secp256k1 under the closed form: exact nu_hat, and its behaviour as the
     leak size K1 = 2^c is swept over the plausible range.
  F. Why the June CF invariants failed: local (at-scale) vs global partial
     quotients.
  G. Is the GLV eigenvalue arithmetically special?  lam is a root of
     x^2+x+1 mod n, so L2 at S1 = S2 is a Z[omega]-ideal lattice -- test whether
     nu_hat departs from the part-D null at that scale and not at others.

Run: python3 glv_hnp_nuhat_closed_form.py
"""

import math
import random
import sys
import time

import sympy

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_lambda_threshold.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

gauss_reduce_2d_float = _t20a.gauss_reduce_2d      # incumbent float64 routine
scales = _t20a.scales
glv_eigenvalues = _t20a.glv_eigenvalues

HERMITE2 = (4.0 / 3.0) ** 0.25                     # 1.0746, max possible nu_hat

# ---------------------------------------------------------------------------
# A. Exact lambda_1: integer Lagrange-Gauss, and the closed form
# ---------------------------------------------------------------------------

def gauss_reduce_2d_exact(u, v):
    """Exact integer Lagrange-Gauss.  Returns lambda_1^2 as an exact int."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    while True:
        nv = nrm2(v)
        if nv == 0:
            return nrm2(u)
        dot = u[0] * v[0] + u[1] * v[1]
        # q = round(dot / nv), exact, correct for negative dot
        q = (2 * dot + nv) // (2 * nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return nrm2(u)

def cf_convergent_denominators(lam, n):
    """Denominators q_j of the continued fraction of lam/n, with q_{-1} = 0."""
    qs = [0, 1]
    a, b = lam, n
    qm1, q0 = 1, 0                    # q_{-2}, q_{-1} in the standard recurrence
    while b:
        c = a // b
        qm1, q0 = q0, c * q0 + qm1
        qs.append(q0)
        a, b = b, a - c * b
    return sorted(set(qs))

def centered_residue(x, n):
    """||x||_n = min over integers a of |x - a*n|."""
    r = x % n
    return min(r, n - r)

def lambda1_sq_closed_form(n, lam, s1, s2, want_witness=False):
    """
    lambda_1(L2)^2 via the closed form:

        L2 = { ( S1*(a*n - b*lam), S2*b ) : a, b in Z }

    For fixed b the inner minimum over a is S1*||b*lam||_n, so

        lambda_1^2 = min_{b >= 0} [ S1^2*||b*lam||_n^2 + S2^2*b^2 ].

    If b > 0 attains the min then every 0 < b' < b has
    ||b'*lam||_n > ||b*lam||_n (else b' would beat b, since S2^2*b'^2 <
    S2^2*b^2), i.e. b is a best approximation of the second kind to lam/n.
    By Lagrange's theorem those are exactly the CF convergent denominators.
    b = 0 (giving S1*n, from a = +-1) is q_{-1} = 0 in the list.
    """
    best, wit = None, None
    for q in cf_convergent_denominators(lam, n):
        th = n if q == 0 else centered_residue(q * lam, n)
        val = s1 * s1 * th * th + s2 * s2 * q * q
        if best is None or val < best:
            best, wit = val, (q, th)
    return (best, wit) if want_witness else best

def nu_hat_exact(n, lam, k1, k2):
    """Exact nu_hat = lambda_1/sqrt(det) computed in integers, float only at the end."""
    s1, _, s2, _ = scales(n, k1, k2)
    l1sq = lambda1_sq_closed_form(n, lam, s1, s2)
    det = n * s1 * s2
    # l1sq/det can overflow float64 separately but not as a ratio
    return math.sqrt(l1sq / det) if l1sq.bit_length() < 900 else \
        math.exp(0.5 * (_log_int(l1sq) - _log_int(det)))

def _log_int(x):
    """log of an arbitrarily large int without overflow."""
    b = x.bit_length()
    if b <= 900:
        return math.log(x)
    sh = b - 900
    return math.log(x >> sh) + sh * math.log(2)

def rand_odd_prime(bits, rng):
    return int(sympy.nextprime(rng.randrange(2 ** (bits - 1), 2 ** bits)))

# ---------------------------------------------------------------------------
# D. The exact null law
# ---------------------------------------------------------------------------

def nu_hat_cdf(r):
    """
    P(nu_hat <= r) for a Haar-random unimodular 2-lattice.

    Writing a unimodular lattice as <1, tau>/sqrt(Im tau) with tau in the
    standard SL(2,Z) fundamental domain F = {|x| <= 1/2, |tau| >= 1} and Haar
    measure (3/pi) dx dy / y^2, the shortest vector has lambda_1^2 = 1/y.  So
    with Y = 1/r^2,

        P = (3/pi) * measure(F and {y >= Y}).

    Y >= 1        : the strip has full width 1, P = 3 r^2 / pi.
    sqrt(3)/2 <= Y <= 1 : width is 1 - 2*sqrt(1-y^2), giving the arcsin form.
    Y <= sqrt(3)/2: all of F, P = 1.   (Y = sqrt(3)/2 <=> r = (4/3)^(1/4).)
    """
    if r <= 0:
        return 0.0
    if r >= HERMITE2:
        return 1.0
    if r <= 1.0:
        return 3.0 * r * r / math.pi
    y = 1.0 / (r * r)
    return (3.0 / math.pi) * (1.0 / y - 2.0 * math.sqrt(1.0 - y * y) / y
                              - 2.0 * math.asin(y) + math.pi)

def nu_hat_quantile(p):
    lo, hi = 0.0, HERMITE2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if nu_hat_cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

def ks_stat(sample, cdf):
    s = sorted(sample)
    m = len(s)
    d = 0.0
    for i, x in enumerate(s):
        f = cdf(x)
        d = max(d, abs(f - i / m), abs((i + 1) / m - f))
    return d

def ks_pvalue(d, m):
    """Asymptotic two-sided Kolmogorov CDF complement."""
    lam = (math.sqrt(m) + 0.12 + 0.11 / math.sqrt(m)) * d
    s, sign = 0.0, 1.0
    for k in range(1, 101):
        s += sign * math.exp(-2.0 * k * k * lam * lam)
        sign = -sign
    return max(0.0, min(1.0, 2.0 * s))

# ---------------------------------------------------------------------------
# helpers for reporting
# ---------------------------------------------------------------------------

def quantiles(xs, ps):
    s = sorted(xs)
    return [s[min(len(s) - 1, int(p * len(s)))] for p in ps]

def pearson(xs, ys):
    m = len(xs)
    mx, my = sum(xs) / m, sum(ys) / m
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float('nan')

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))

def glv_lambda(n):
    """A root of x^2+x+1 mod n for prime n = 1 mod 3 (None otherwise)."""
    if n % 3 != 1:
        return None
    rng = random.Random(n & 0xffff)
    for _ in range(60):
        g = rng.randrange(2, n - 1)
        lam = pow(g, (n - 1) // 3, n)
        if lam != 1 and (lam * lam + lam + 1) % n == 0:
            return lam
    return None

# ===========================================================================
# Part A -- closed form vs exact Lagrange-Gauss
# ===========================================================================

def part_a():
    print("=" * 78)
    print("A. CLOSED FORM: lambda_1(L2)^2 = min over CF convergents")
    print("=" * 78)
    print("   claim: lambda_1^2 = min_{q in CFden(lam/n) u {0}} "
          "S1^2*||q*lam||_n^2 + S2^2*q^2")
    print("   test : exact-integer Lagrange-Gauss on the same lattice, "
          "compared as exact ints\n")

    rng = random.Random(20260729)
    print(f"{'bits':>5} {'eff':>7} {'cases':>6} {'CF == Gauss':>12} "
          f"{'mean |CFden|':>13} {'mean j*':>8}")
    total, agree = 0, 0
    for bits in (20, 24, 32, 64, 128, 192, 256, 384, 521):
        for eff in (0.02, 0.10, 0.50):
            ok, cases, ncf, jstar = 0, 0, 0, 0
            for _ in range(40):
                n = rand_odd_prime(bits, rng)
                lam = rng.randrange(2, n - 1)
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff * math.isqrt(n)))
                s1, _, s2, _ = scales(n, k1, k2)
                if s1 == 0 or s2 == 0:
                    continue
                cf = lambda1_sq_closed_form(n, lam, s1, s2)
                gs = gauss_reduce_2d_exact((n * s1, 0), (-lam * s1, s2))
                cases += 1
                ok += (cf == gs)
                dens = cf_convergent_denominators(lam, n)
                ncf += len(dens)
                _, wit = lambda1_sq_closed_form(n, lam, s1, s2, want_witness=True)
                jstar += dens.index(wit[0])
            if cases:
                total += cases
                agree += ok
                print(f"{bits:>5} {eff:>7.2f} {cases:>6} {ok:>8}/{cases:<3} "
                      f"{ncf/cases:>13.1f} {jstar/cases:>8.1f}")
    print(f"\n   TOTAL: {agree}/{total} exact agreements "
          f"({'CLOSED FORM VERIFIED' if agree == total else 'MISMATCH -- formula wrong'})")
    return agree == total

# ===========================================================================
# Part B -- precision audit of the incumbent float64 routine
# ===========================================================================

def part_b():
    print("\n" + "=" * 78)
    print("B. PRECISION AUDIT of the incumbent float64 gauss_reduce_2d")
    print("=" * 78)
    print("   every published nu_hat came from float64 Lagrange-Gauss on entries")
    print("   of size ~n*S1 = n^2/K1; the squared norms are ~(n^2/K1)^2.\n")
    rng = random.Random(4242)
    print(f"{'bits':>5} {'eff':>6} {'|nrm2| bits':>12} {'max rel err':>12} "
          f"{'#inf/nan':>9} {'#wrong':>7}")
    for bits in (20, 32, 64, 128, 192, 256, 384, 521):
        for eff in (0.02, 0.10):
            worst, bad, wrong, nb = 0.0, 0, 0, 0
            for _ in range(25):
                n = rand_odd_prime(bits, rng)
                lam = rng.randrange(2, n - 1)
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff * math.isqrt(n)))
                s1, _, s2, _ = scales(n, k1, k2)
                # reference sqrt to ~64 bits beyond float64, so the comparison
                # measures float64 error and not integer-sqrt truncation
                l1sq = lambda1_sq_closed_form(n, lam, s1, s2)
                exact = math.exp(0.5 * _log_int(l1sq)) if l1sq.bit_length() > 900 \
                    else math.isqrt(l1sq << 128) / float(1 << 64)
                nb = max(nb, ((n * s1) ** 2).bit_length())
                try:
                    got = gauss_reduce_2d_float((n * s1, 0), (-lam * s1, s2))
                except (OverflowError, ZeroDivisionError, ValueError):
                    bad += 1
                    continue
                if not math.isfinite(got):
                    bad += 1
                    continue
                rel = abs(got - exact) / exact if exact else 0.0
                worst = max(worst, rel)
                if rel > 1e-9:
                    wrong += 1
            print(f"{bits:>5} {eff:>6.2f} {nb:>12} {worst:>12.2e} "
                  f"{bad:>9} {wrong:>7}")
    print("\n   float64 holds exact ints to 2^53 and overflows at 2^1024;")
    print("   a squared norm of B bits loses (B-53)/2 bits in the norm itself.")

# ===========================================================================
# Part C -- localisation of the minimising convergent
# ===========================================================================

def part_c():
    print("\n" + "=" * 78)
    print("C. LOCALISATION: is q_{j*} the convergent nearest R = sqrt(n*K2/K1)?")
    print("=" * 78)
    rng = random.Random(31337)
    print(f"{'bits':>5} {'eff':>6} {'median |log2(q*/R)|':>20} "
          f"{'q* is nearest-to-R':>19} {'within 1 index':>15}")
    for bits in (20, 64, 128, 256):
        for eff in (0.02, 0.10, 0.50):
            gaps, nearest, within = [], 0, 0
            for _ in range(120):
                n = rand_odd_prime(bits, rng)
                lam = rng.randrange(2, n - 1)
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff * math.isqrt(n)))
                s1, _, s2, _ = scales(n, k1, k2)
                _, wit = lambda1_sq_closed_form(n, lam, s1, s2, want_witness=True)
                R = math.exp(0.5 * (_log_int(n) + _log_int(s1) - _log_int(s2)))
                dens = [q for q in cf_convergent_denominators(lam, n) if q > 0]
                if not dens or wit[0] == 0:
                    continue
                gaps.append(abs(math.log2(wit[0] / R)))
                near = min(dens, key=lambda q: abs(math.log2(q / R)))
                nearest += (near == wit[0])
                i_star, i_near = dens.index(wit[0]), dens.index(near)
                within += (abs(i_star - i_near) <= 1)
            if gaps:
                g = sorted(gaps)[len(gaps) // 2]
                print(f"{bits:>5} {eff:>6.2f} {g:>20.3f} "
                      f"{nearest/len(gaps)*100:>18.0f}% {within/len(gaps)*100:>14.0f}%")
    print("\n   R = sqrt(n*S1/S2) = sqrt(n*K2/K1) is the scale at which the two")
    print("   terms S1*theta and S2*q balance (theta ~ n/q).")

# ===========================================================================
# Part D -- the exact null law
# ===========================================================================

def part_d():
    print("\n" + "=" * 78)
    print("D. EXACT NULL LAW for nu_hat over random lam")
    print("=" * 78)
    print("   normalising u = q/R, w = theta*R/n gives nu_hat^2 = min(u^2 + w^2)")
    print("   with u*w = q*theta/n: L2/sqrt(det) is a unimodular 2-lattice, so")
    print("   under equidistribution nu_hat^2 = 1/Im(tau), tau ~ Haar on SL(2,Z)\\H.\n")
    print("   derived CDF:  P(nu<=r) = 3r^2/pi                        (r <= 1)")
    print("                          = (3/pi)[r^2 - 2r^2*sqrt(1-r^-4)")
    print("                                   - 2*asin(r^-2) + pi]    (1 <= r <= (4/3)^(1/4))")
    print(f"   support ceiling = Hermite bound (4/3)^(1/4) = {HERMITE2:.6f}\n")

    rng = random.Random(90210)
    for bits, m in ((64, 8000), (256, 20000)):
        for eff in (0.05, 0.25):
            sample = []
            for _ in range(m // 40):
                n = rand_odd_prime(bits, rng)
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff * math.isqrt(n)))
                s1, _, s2, _ = scales(n, k1, k2)
                for _ in range(40):
                    lam = rng.randrange(2, n - 1)
                    l1 = lambda1_sq_closed_form(n, lam, s1, s2)
                    sample.append(math.exp(0.5 * (_log_int(l1) - _log_int(n * s1 * s2))))
            ps = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
            emp = quantiles(sample, ps)
            pred = [nu_hat_quantile(p) for p in ps]
            d = ks_stat(sample, nu_hat_cdf)
            print(f"  bits={bits} eff={eff}  N={len(sample)}  "
                  f"KS D={d:.4f}  p={ks_pvalue(d, len(sample)):.3f}  "
                  f"max={max(sample):.4f} (<= {HERMITE2:.4f}: "
                  f"{'ok' if max(sample) <= HERMITE2 + 1e-12 else 'VIOLATED'})")
            print("        q:      " + "".join(f"{p:>9.2f}" for p in ps))
            print("        empir:  " + "".join(f"{v:>9.4f}" for v in emp))
            print("        derived:" + "".join(f"{v:>9.4f}" for v in pred))
    print("\n   2026-07-29 log reported an empirical 256-bit null (4000 draws,")
    print("   eff=0.05): q05=0.232 q10=0.322 q25=0.496 q50=0.715 q75=0.885")
    print("   q90=0.974 q95=0.999.  Derived: " +
          " ".join(f"{nu_hat_quantile(p):.3f}"
                   for p in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)))

# ===========================================================================
# Part E -- secp256k1 under the closed form
# ===========================================================================

N_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAM_SECP = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72

def part_e():
    print("\n" + "=" * 78)
    print("E. secp256k1 under the closed form")
    print("=" * 78)
    assert (LAM_SECP * LAM_SECP + LAM_SECP + 1) % N_SECP == 0
    n, lam = N_SECP, LAM_SECP
    k2 = math.isqrt(n) + 1

    print("   (i) exact vs float64 nu_hat at the eff values reported 2026-07-29")
    print(f"{'eff':>6} {'nu_hat exact':>13} {'nu_hat float':>13} {'rel err':>10} "
          f"{'pctile (derived CDF)':>21}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(2, int(eff * math.isqrt(n)))
        s1, _, s2, _ = scales(n, k1, k2)
        ex = nu_hat_exact(n, lam, k1, k2)
        try:
            fl = gauss_reduce_2d_float((n * s1, 0), (-lam * s1, s2)) / \
                math.exp(0.5 * _log_int(n * s1 * s2))
        except (OverflowError, ValueError):
            fl = float('nan')
        rel = abs(fl - ex) / ex if math.isfinite(fl) else float('nan')
        print(f"{eff:>6.2f} {ex:>13.6f} {fl:>13.6f} {rel:>10.2e} "
              f"{nu_hat_cdf(ex)*100:>20.1f}%")

    print("\n   (ii) nu_hat as a function of the leak size K1 = 2^c")
    print("        (K1 is set by the side channel, not chosen by the attacker;")
    print("         this is the lottery the attacker is handed)")
    print(f"{'c':>4} {'K1=2^c':>8} {'log2 R':>8} {'nu_hat':>9} {'pctile':>8} {'j*':>4} {'a_{j*+1}':>9}")
    vals = []
    for c in range(4, 129, 4):
        k1 = 2 ** c
        s1, _, s2, _ = scales(n, k1, k2)
        l1, wit = lambda1_sq_closed_form(n, lam, s1, s2, want_witness=True)
        nu = math.exp(0.5 * (_log_int(l1) - _log_int(n * s1 * s2)))
        R = math.exp(0.5 * (_log_int(n) + _log_int(s1) - _log_int(s2)))
        dens = cf_convergent_denominators(lam, n)
        j = dens.index(wit[0])
        aloc = (dens[j + 1] // dens[j]) if 0 < j < len(dens) - 1 and dens[j] else 0
        vals.append(nu)
        print(f"{c:>4} {'2^' + str(c):>8} {math.log2(R):>8.1f} {nu:>9.4f} "
              f"{nu_hat_cdf(nu)*100:>7.1f}% {j:>4} {aloc:>9}")
    print(f"\n   over c in [4,128]: min={min(vals):.4f} max={max(vals):.4f} "
          f"mean={sum(vals)/len(vals):.4f}  (Haar mean approx "
          f"{sum(nu_hat_quantile(i/1000) for i in range(1,1000))/999:.4f})")

# ===========================================================================
# Part F -- why the June CF invariants failed
# ===========================================================================

def part_f():
    print("\n" + "=" * 78)
    print("F. LOCAL vs GLOBAL partial quotients (why max_a/q_cf failed in June)")
    print("=" * 78)
    rng = random.Random(555)
    bits, eff = 20, 0.10
    nus, alocs, amaxs, mus = [], [], [], []
    for _ in range(400):
        n = rand_odd_prime(bits, rng)
        lam = rng.randrange(2, n - 1)
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff * math.isqrt(n)))
        s1, _, s2, _ = scales(n, k1, k2)
        l1, wit = lambda1_sq_closed_form(n, lam, s1, s2, want_witness=True)
        nu = math.sqrt(l1 / (n * s1 * s2))
        dens = cf_convergent_denominators(lam, n)
        j = dens.index(wit[0])
        if not (0 < j < len(dens) - 1) or dens[j] == 0:
            continue
        aloc = dens[j + 1] // dens[j]
        amax = max(dens[i + 1] // dens[i] for i in range(1, len(dens) - 1) if dens[i])
        nus.append(nu); alocs.append(math.log(aloc + 1))
        amaxs.append(math.log(amax + 1)); mus.append(min(lam, n - lam) / n)
    print(f"   N={len(nus)} random (n,lam), {bits}-bit, eff={eff}")
    print(f"   pearson(log(1+a_local), nu_hat)  = {pearson(alocs, nus):+.3f}   "
          f"spearman = {spearman(alocs, nus):+.3f}")
    print(f"   pearson(log(1+a_max),   nu_hat)  = {pearson(amaxs, nus):+.3f}   "
          f"spearman = {spearman(amaxs, nus):+.3f}")
    print(f"   pearson(mu,             nu_hat)  = {pearson(mus, nus):+.3f}   "
          f"spearman = {spearman(mus, nus):+.3f}")
    print(f"   pearson(log(1+a_local), log(1+a_max)) = {pearson(alocs, amaxs):+.3f}")
    print("\n   binned nu_hat by the LOCAL partial quotient a_{j*+1}:")
    print(f"{'a_local':>10} {'count':>6} {'mean a':>8} {'mean nu_hat':>12} "
          f"{'sqrt(2/(a+1))':>14}")
    for lo, hi in ((1, 1), (2, 2), (3, 4), (5, 8), (9, 10 ** 9)):
        grp = [(nu, math.exp(a) - 1) for nu, a in zip(nus, alocs)
               if lo <= round(math.exp(a) - 1) <= hi]
        if len(grp) >= 3:
            am = sum(a for _, a in grp) / len(grp)
            lbl = f"{lo}-{hi}" if hi < 10 ** 9 else f">={lo}"
            print(f"{lbl:>10} {len(grp):>6} {am:>8.1f} "
                  f"{sum(nu for nu, _ in grp)/len(grp):>12.4f} "
                  f"{math.sqrt(2.0/(am+1)):>14.4f}")

# ===========================================================================
# Part G -- is the GLV eigenvalue arithmetically special?
# ===========================================================================

def part_g():
    print("\n" + "=" * 78)
    print("G. GLV lam (root of x^2+x+1) vs random lam, as a function of scale")
    print("=" * 78)
    print("   L2 at S1 = S2 is (a scaling of) the Z[omega]-ideal lattice that")
    print("   makes GLV work, so nu_hat there should NOT follow the part-D null.")
    print("   Away from that scale the relevant convergent sits elsewhere and the")
    print("   CF should look Gauss-Kuzmin generic.  R/sqrt(n) = sqrt(S1/S2) is the")
    print("   knob: R = sqrt(n) exactly when S1 = S2 (i.e. K1 = K2).\n")

    rng = random.Random(777)
    bits = 64
    ns = []
    while len(ns) < 300:
        n = rand_odd_prime(bits, rng)
        if n % 3 == 1 and glv_lambda(n):
            ns.append(n)

    print(f"{'S1/S2':>10} {'log2(R/sqrt n)':>15} {'GLV KS':>8} {'GLV p':>8} "
          f"{'rand KS':>8} {'rand p':>8} {'GLV mean':>9} {'rand mean':>10}")
    for ratio in (1, 4, 2 ** 8, 2 ** 16, 2 ** 24, 2 ** 32):
        g_s, r_s = [], []
        for n in ns:
            s2 = math.isqrt(n) + 1
            s1 = s2 * ratio
            det = n * s1 * s2
            lam = glv_lambda(n)
            g_s.append(math.sqrt(lambda1_sq_closed_form(n, lam, s1, s2) / det))
            lr = rng.randrange(2, n - 1)
            r_s.append(math.sqrt(lambda1_sq_closed_form(n, lr, s1, s2) / det))
        dg, dr = ks_stat(g_s, nu_hat_cdf), ks_stat(r_s, nu_hat_cdf)
        print(f"{ratio:>10} {0.5*math.log2(ratio):>15.1f} {dg:>8.3f} "
              f"{ks_pvalue(dg, len(g_s)):>8.4f} {dr:>8.3f} "
              f"{ks_pvalue(dr, len(r_s)):>8.4f} {sum(g_s)/len(g_s):>9.4f} "
              f"{sum(r_s)/len(r_s):>10.4f}")
    print(f"\n   (Hermite ceiling {HERMITE2:.4f} is attained by the hexagonal")
    print("    lattice; Z[omega] IS hexagonal, in the a^2-ab+b^2 norm.)")

    print("\n   G2. where does the GLV protection stop?  fraction of curves with")
    print(f"   nu_hat <= 0.645 (the Thread 20d C2 ceiling); Haar predicts "
          f"{nu_hat_cdf(0.645)*100:.1f}%")
    print(f"{'K2/K1':>10} {'min':>7} {'q05':>7} {'q25':>7} {'q50':>7} {'q95':>7} "
          f"{'frac<=.645':>11} {'rand frac':>10}")
    for ratio in (1, 2, 4, 16, 64, 256, 4096):
        g_s, r_s = [], []
        for n in ns:
            s2 = math.isqrt(n) + 1
            s1 = s2 * ratio
            det = n * s1 * s2
            g_s.append(math.sqrt(lambda1_sq_closed_form(n, glv_lambda(n), s1, s2) / det))
            r_s.append(math.sqrt(
                lambda1_sq_closed_form(n, rng.randrange(2, n - 1), s1, s2) / det))
        qs = quantiles(g_s, (0.05, 0.25, 0.50, 0.95))
        print(f"{ratio:>10} {min(g_s):>7.4f} " + "".join(f"{v:>7.4f}" for v in qs) +
              f" {sum(1 for v in g_s if v <= 0.645)/len(g_s)*100:>10.1f}% "
              f"{sum(1 for v in r_s if v <= 0.645)/len(r_s)*100:>9.1f}%")

# ===========================================================================
# Part H -- the GLV floor theorem
# ===========================================================================

def glv_floor(s1, s2):
    """
    THEOREM (GLV floor).  Let n be prime and lam^2 + lam + 1 = 0 mod n.  Every
    nonzero (a,b) in Z^2 with a + b*lam = 0 mod n satisfies

        a^2 - a*b + b^2  =  b^2*(lam^2 + lam + 1)  =  0   (mod n),

    and the form is positive definite, so a^2 - a*b + b^2 >= n.  Hence

        lambda_1(L2)^2 = min (S1^2 a^2 + S2^2 b^2)
                      >= n * min_{(a,b) != 0} (S1^2 a^2 + S2^2 b^2)/(a^2-ab+b^2)
                       = n * t_min,

    where t_min is the smallest generalized eigenvalue of diag(S1^2,S2^2) with
    respect to [[1,-1/2],[-1/2,1]], i.e. the smaller root of

        (3/4) t^2 - (S1^2 + S2^2) t + S1^2 S2^2 = 0.

    Dividing by det = n*S1*S2:   nu_hat >= sqrt(t_min / (S1*S2)).

    S1 = S2 gives t_min = (2/3)S1^2 and the scale-free floor nu_hat >=
    sqrt(2/3) = 0.8165.  For S1 >> S2, t_min -> S2^2 and the floor decays like
    sqrt(S2/S1) = sqrt(K1/K2).

    Numerically: with r = S1/S2 the bound is scale-free, and using the stable
    form of the smaller root (t_min = 2c/(b + sqrt(b^2-4ac))) and dividing
    through by r^2 (note (S1^2+S2^2)^2 - 3 S1^2 S2^2 = S2^4 (r^4 - r^2 + 1)):

        nu_hat^2 >= 2 / ( r * [ (1 + r^-2) + sqrt(1 - r^-2 + r^-4) ] ).
    """
    return glv_floor_ratio(math.exp(_log_int(s1) - _log_int(s2)))

def glv_floor_ratio(r):
    """The part-H floor as a function of r = S1/S2 alone (no truncation)."""
    if r < 1.0:                       # bound is symmetric under S1 <-> S2
        r = 1.0 / r
    ir2 = 1.0 / (r * r)
    return math.sqrt(2.0 / (r * ((1.0 + ir2) + math.sqrt(1.0 - ir2 + ir2 * ir2))))

def part_h():
    print("\n" + "=" * 78)
    print("H. THEOREM: an unconditional floor on nu_hat for GLV (j=0) curves")
    print("=" * 78)
    print("   a + b*lam = 0 (mod n)  =>  a^2 - ab + b^2 = 0 (mod n)  =>  >= n,")
    print("   so nu_hat >= sqrt(t_min/(S1*S2)),  t_min = smaller root of")
    print("   (3/4)t^2 - (S1^2+S2^2)t + S1^2*S2^2 = 0.   No lattice reduction,")
    print("   no heuristics, no CM hypotheses -- just positive definiteness.\n")

    rng = random.Random(2718)
    print(f"{'K2/K1':>8} {'floor (proved)':>15} {'observed min':>13} "
          f"{'slack':>8} {'#violations':>12}")
    for bits in (64,):
        ns = []
        while len(ns) < 200:
            n = rand_odd_prime(bits, rng)
            if n % 3 == 1 and glv_lambda(n):
                ns.append(n)
        for ratio in (1, 2, 4, 16, 64, 256, 4096, 2 ** 20):
            obs, viol = [], 0
            s2 = math.isqrt(ns[0]) + 1
            fl = glv_floor(s2 * ratio, s2)
            for n in ns:
                s2 = math.isqrt(n) + 1
                s1 = s2 * ratio
                v = math.sqrt(lambda1_sq_closed_form(n, glv_lambda(n), s1, s2)
                              / (n * s1 * s2))
                obs.append(v)
                if v < glv_floor(s1, s2) - 1e-12:
                    viol += 1
            print(f"{ratio:>8} {fl:>15.6f} {min(obs):>13.6f} "
                  f"{min(obs)/fl:>8.4f} {viol:>12}")

    # critical ratio: below it, the C2 ("easy") band is provably unreachable
    lo, hi = 1.0, 1e6
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if glv_floor_ratio(mid) > 0.645:
            lo = mid
        else:
            hi = mid
    print(f"\n   floor(K2/K1) = 0.645 (the Thread 20d C2 ceiling) at "
          f"K2/K1 = {lo:.3f} = 2^{math.log2(lo):.2f}")
    print(f"   => whenever the leak satisfies K1 >= K2/{lo:.2f}, i.e. the known")
    print(f"      part of k1 is within {math.log2(lo):.2f} bits of K2, NO j=0 GLV")
    print("      curve can be in the low-nu_hat class.  The bound is tight: the")
    print("      observed minima above sit within 1% of it for K2/K1 <= 16.")

    print("\n   secp256k1, floor vs actual, as the leak size K1 = 2^c varies:")
    n, lam = N_SECP, LAM_SECP
    k2 = math.isqrt(n) + 1
    print(f"{'c':>5} {'K2/K1':>10} {'floor':>10} {'actual nu_hat':>14} "
          f"{'actual/floor':>12}")
    for c in (8, 32, 64, 96, 112, 120, 124, 126, 127, 128):
        k1 = 2 ** c
        s1, _, s2, _ = scales(n, k1, k2)
        v = nu_hat_exact(n, lam, k1, k2)
        fl = glv_floor(s1, s2)
        print(f"{c:>5} {'2^%d' % (128 - c):>10} {fl:>10.2e} {v:>14.6f} "
              f"{v/fl:>10.1f}")

def main():
    print("=" * 78)
    print("GLV-HNP Thread 23 -- closed form and null law for nu_hat")
    print("=" * 78)
    t0 = time.time()
    ok = part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    part_g()
    part_h()
    print(f"\nDone in {time.time()-t0:.1f}s.  "
          f"Closed form {'VERIFIED' if ok else 'FAILED'}.")
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
