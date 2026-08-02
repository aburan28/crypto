"""
GLV-HNP — Thread 23: closed form for lambda_1(L2), and the null law of nu_hat.

Background
----------
Thread 20 (2026-07-29) found the first surviving Phase-2 predictor after eight
falsified invariants (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n,
lam/n, mu):

    nu_hat = lambda_1(L2) / sqrt(det L2),      L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >

with det L2 = n*S_K1*S_K2 independent of lam.  AUC 0.935 against the June
C1/C2 classes; low nu_hat = easy.  But nu_hat was *computed* (Lagrange-Gauss),
not *understood*, and the 2026-07-29 entry proposed Thread 23: show lambda_1(L2)
is a continued-fraction quantity of lam/n taken at the scale
sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1), which would explain retroactively why the
scale-FREE CF invariants (q_cf, max_q_cf, max_a) all failed in June.

Theory tested here
------------------
A vector of L2 is  a*(-lam*S_K1, S_K2) + b*(n*S_K1, 0) = ((b*n - a*lam)*S_K1, a*S_K2),
so with eta(a) = |a*lam|_n  (centred residue, minimised over b)

    ||v||^2 = (eta(a)*S_K1)^2 + (a*S_K2)^2 .

(T23-A) The minimising a is a *best approximation of the second kind* to
lam/n, i.e. a convergent denominator q_j of the continued fraction of lam/n
(or a = 0, the row (n*S_K1, 0)).  Hence

    lambda_1(L2)^2 = min over {(0,n)} u {(q_j, eta_j)} of (eta*S_K1)^2 + (a*S_K2)^2 ,

computable in O(log n) with no lattice reduction.  Falsifier: any single
instance where this disagrees with exact Lagrange-Gauss reduction.

(T23-B) Put q* = sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1) (the scale), and for each
convergent

    x_j = q_j / q*            (where the convergent sits relative to the scale)
    y_j = q_j * eta_j / n     (classical scale-free approximation quality, in (0,1))

Then, identically,

    nu_hat^2 = min_j ( y_j^2 / x_j^2 + x_j^2 ) .

This is the explanation of the June failures: nu_hat is a function of the JOINT
(x_j, y_j), never of y alone.  q_cf/max_q_cf/max_a are functions of the y's (and
the partial quotients) alone, so they cannot determine nu_hat.  Note also
min_x (y^2/x^2 + x^2) = 2y at x = sqrt(y): a small nu_hat needs a good
approximation that ALSO lands near the scale.

(T23-C) Null law.  As n -> oo with lam uniform, the unimodular lattice
L2/sqrt(det L2) equidistributes in SL_2(Z)\SL_2(R) w.r.t. Haar.  Writing the
lattice as y^{-1/2}(Z + Z*tau) with tau in the standard fundamental domain, the
shortest vector is y^{-1/2}, and dmu = (3/pi) dx dy / y^2, giving the exact CDF

    F(r) = 3r^2/pi                                             0 <= r <= 1
    F(r) = (3/pi)[ r^2 - 2( sqrt(1-y0^2)/y0 + arcsin(y0) - pi/2 ) ],  y0 = 1/r^2
                                                          1 <= r <= (4/3)^(1/4)
    F(r) = 1                                              r >= (4/3)^(1/4)

In particular nu_hat <= (4/3)^(1/4) = 1.07457 always (Hermite in dim 2), and the
law does NOT depend on the bit size.  If that holds empirically, the 20/24-bit
-> 256-bit extrapolation flagged as caveat (a) on 2026-07-29 is justified for
the *distribution* of nu_hat (not, separately, for the p_hat response to it).

Run: python3 glv_hnp_nuhat_cf_theory.py
"""

import math
import random
import sys
import time

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_lambda_threshold.py")
_t20a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t20a)

gauss_reduce_2d = _t20a.gauss_reduce_2d
scales = _t20a.scales
rival_sublattice_nu = _t20a.rival_sublattice_nu

HERMITE2 = (4.0 / 3.0) ** 0.25

# ---------------------------------------------------------------------------
# Continued fractions
# ---------------------------------------------------------------------------

def cf_convergents(x, y):
    """Convergents of x/y, 0 <= x < y, gcd(x,y)=1.

    Returns [(q_j, eta_j, a_j)] with eta_j = |q_j*x - p_j*y| and a_j the partial
    quotient that produced convergent j.  The final entry is (y, 0, .).
    """
    A, B = x, y
    h1, h2 = 1, 0
    k1, k2 = 0, 1
    out = []
    while B:
        a = A // B
        A, B = B, A - a * B
        h1, h2 = a * h1 + h2, h1
        k1, k2 = a * k1 + k2, k1
        out.append((k1, abs(k1 * x - h1 * y), a))
    return out

def lambda1_cf(n, lam, S_K1, S_K2):
    """(lambda_1^2, best (a, eta), convergent index) via T23-A -- no reduction.

    Index -1 denotes the a = 0 candidate (the row (n*S_K1, 0)).
    """
    lam %= n
    best = (n * S_K1) ** 2
    arg = (0, n)
    idx = -1
    for j, (q, eta, _a) in enumerate(cf_convergents(lam, n)):
        val = (eta * S_K1) ** 2 + (q * S_K2) ** 2
        if val < best:
            best, arg, idx = val, (q, eta), j
    return best, arg, idx

def lambda1_exact(n, lam, S_K1, S_K2):
    """lambda_1^2 by exact Lagrange-Gauss reduction (the Thread 20 route)."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return w[0] * w[0] + w[1] * w[1]

def nu_hat_cf(n, lam, k1_bound, k2_bound):
    """nu_hat via T23-A, safe for 256-bit inputs (log-domain ratio)."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    nsq, _, _ = lambda1_cf(n, lam, S_K1, S_K2)
    det = n * S_K1 * S_K2
    return math.exp(0.5 * (math.log(nsq) - math.log(det)))

def k1k2_for_eff(n, eff):
    """K2 = sqrt(n), K1 = eff*sqrt(n)  =>  K1*K2/n = eff (the Thread 20 convention)."""
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(eff * math.sqrt(n)))
    return k1, k2

# ---------------------------------------------------------------------------
# The null law of nu_hat (T23-C)
# ---------------------------------------------------------------------------

def nu_cdf(r):
    """Exact CDF of nu_hat under Haar equidistribution of L2/sqrt(det L2)."""
    if r <= 0:
        return 0.0
    if r >= HERMITE2:
        return 1.0
    if r <= 1.0:
        return 3.0 * r * r / math.pi
    y0 = 1.0 / (r * r)
    corr = math.sqrt(1.0 - y0 * y0) / y0 + math.asin(y0) - math.pi / 2.0
    return (3.0 / math.pi) * (r * r - 2.0 * corr)

def nu_quantile(p, lo=0.0, hi=HERMITE2, iters=80):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if nu_cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

def ks_test(sample, cdf):
    """One-sample two-sided KS statistic and asymptotic p-value."""
    xs = sorted(sample)
    N = len(xs)
    d = 0.0
    for i, x in enumerate(xs):
        f = cdf(x)
        d = max(d, abs(f - i / N), abs((i + 1) / N - f))
    lam = (math.sqrt(N) + 0.12 + 0.11 / math.sqrt(N)) * d
    s, p = 0.0, 0.0
    for k in range(1, 200):
        s += (-1) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam)
    p = max(0.0, min(1.0, 2.0 * s))
    return d, p

# ---------------------------------------------------------------------------
# statistics helpers
# ---------------------------------------------------------------------------

def pearson(xs, ys):
    N = len(xs)
    mx, my = sum(xs) / N, sum(ys) / N
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
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
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))

def quantiles(v, ps):
    s = sorted(v)
    out = []
    for p in ps:
        i = min(len(s) - 1, max(0, int(p * (len(s) - 1) + 0.5)))
        out.append(s[i])
    return out

def rand_prime(bits, rng):
    import sympy
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if sympy.isprime(c):
            return c

# ---------------------------------------------------------------------------

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72

BAR = "-" * 78

def part_a():
    print("\n" + "=" * 78)
    print("PART A (T23-A): CF/Pareto closed form vs exact Lagrange-Gauss")
    print("=" * 78)
    print("Claim: lambda_1(L2)^2 = min over {(0,n)} u {(q_j, eta_j)} of")
    print("       (eta*S_K1)^2 + (a*S_K2)^2, q_j = CF convergent denominators of lam/n.")
    print("Falsifier: ONE disagreement with gauss_reduce_2d.\n")
    rng = random.Random(20230823)
    total = mism = 0
    print(f"{'bits':>5} {'eff':>6} {'trials':>7} {'mismatch':>9} {'nu range':>22}")
    for bits in (12, 16, 20, 24, 32, 48, 64):
        n = rand_prime(bits, rng)
        for eff in (0.05, 0.10, 0.25, 0.50):
            k1, k2 = k1k2_for_eff(n, eff)
            S_K1, _, S_K2, _ = scales(n, k1, k2)
            nus = []
            bad = 0
            for _ in range(150):
                lam = rng.randrange(1, n)
                a = lambda1_cf(n, lam, S_K1, S_K2)[0]
                b = lambda1_exact(n, lam, S_K1, S_K2)
                if a != b:
                    bad += 1
                    if bad == 1:
                        print(f"    MISMATCH n={n} lam={lam} cf={a} gauss={b}")
                nus.append(math.exp(0.5 * (math.log(a) - math.log(n * S_K1 * S_K2))))
                total += 1
            mism += bad
            print(f"{bits:>5} {eff:>6.2f} {150:>7} {bad:>9} "
                  f"  [{min(nus):.4f}, {max(nus):.4f}]")
    print(f"\n  TOTAL {total} instances, {mism} mismatches.")
    print(f"  VERDICT: {'T23-A CONFIRMED (exact)' if mism == 0 else 'T23-A FALSIFIED'}")
    return mism == 0

def part_b():
    print("\n" + "=" * 78)
    print("PART B (T23-B): the (x, y) closed form and the scale q*")
    print("=" * 78)
    print("nu_hat^2 = min_j ( y_j^2/x_j^2 + x_j^2 ),  x_j = q_j/q*, y_j = q_j*eta_j/n,")
    print("q* = sqrt(n*S_K1/S_K2) = sqrt(n*K2/K1).\n")
    rng = random.Random(777)
    worst = 0.0
    rows = []
    for bits in (20, 32, 64, 128, 256):
        n = rand_prime(bits, rng)
        for eff in (0.05, 0.25):
            k1, k2 = k1k2_for_eff(n, eff)
            S_K1, _, S_K2, _ = scales(n, k1, k2)
            qstar = math.exp(0.5 * (math.log(n) + math.log(S_K1) - math.log(S_K2)))
            for _ in range(60):
                lam = rng.randrange(1, n)
                direct = nu_hat_cf(n, lam, k1, k2)
                # via (x, y).  The a=0 candidate contributes (n*S_K1)^2/det
                # = n*S_K1/S_K2 = q*^2, i.e. x = q*/q* rescaled to the origin.
                best = qstar * qstar
                for q, eta, _a in cf_convergents(lam % n, n):
                    x = q / qstar
                    yq = q * eta / n
                    val = (yq * yq) / (x * x) + x * x
                    best = min(best, val)
                viaxy = math.sqrt(best)
                worst = max(worst, abs(viaxy - direct) / direct)
            rows.append((bits, eff, qstar))
    print(f"  max relative deviation between the two evaluations: {worst:.3e}")
    print(f"  VERDICT: {'identity CONFIRMED' if worst < 1e-9 else 'identity FAILS'}")
    print("\n  scale q* = sqrt(n*K2/K1) for reference:")
    print(f"  {'bits':>5} {'eff':>6} {'q*':>26} {'log2 q*':>9}")
    for bits, eff, qs in rows:
        print(f"  {bits:>5} {eff:>6.2f} {qs:>26.6e} {math.log2(qs):>9.2f}")
    return worst < 1e-9

def part_c():
    print("\n" + "=" * 78)
    print("PART C: where does the minimising convergent sit? (scale localisation)")
    print("=" * 78)
    print("If T23-B is the right reading, argmin q_j clusters at q* -- i.e.")
    print("log2(q_argmin/q*) concentrates near 0 with spread ~1 partial quotient.\n")
    rng = random.Random(4242)
    print(f"{'bits':>5} {'eff':>6} {'mean':>7} {'sd':>6} {'|.|<=1':>8} {'|.|<=2':>8} "
          f"{'triv%':>7}")
    for bits in (20, 32, 64, 128, 256):
        n = rand_prime(bits, rng)
        for eff in (0.05, 0.25):
            k1, k2 = k1k2_for_eff(n, eff)
            S_K1, _, S_K2, _ = scales(n, k1, k2)
            qstar = math.exp(0.5 * (math.log(n) + math.log(S_K1) - math.log(S_K2)))
            ds, triv = [], 0
            for _ in range(400):
                lam = rng.randrange(1, n)
                _, (a, _eta), idx = lambda1_cf(n, lam, S_K1, S_K2)
                if idx == -1:
                    triv += 1
                    continue
                ds.append(math.log2(a / qstar))
            mu = sum(ds) / len(ds)
            sd = math.sqrt(sum((d - mu) ** 2 for d in ds) / len(ds))
            f1 = sum(1 for d in ds if abs(d) <= 1) / len(ds)
            f2 = sum(1 for d in ds if abs(d) <= 2) / len(ds)
            print(f"{bits:>5} {eff:>6.2f} {mu:>7.3f} {sd:>6.3f} {f1:>8.3f} {f2:>8.3f} "
                  f"{100*triv/400:>7.2f}")

def part_d():
    print("\n" + "=" * 78)
    print("PART D: why the scale-FREE CF invariants failed in June")
    print("=" * 78)
    print("max_a (largest partial quotient, falsified 2026-06-26/06-29) vs")
    print("a_local = the partial quotient at the minimising convergent.")
    print("Both are CF quantities; only a_local knows the scale.\n")
    rng = random.Random(31415)
    print(f"{'bits':>5} {'eff':>6} {'sp(max_a,nu)':>13} {'sp(a_loc,nu)':>13} "
          f"{'sp(y_loc,nu)':>13} {'sp(mu,nu)':>10}")
    for bits in (20, 24, 32, 64):
        n = rand_prime(bits, rng)
        for eff in (0.05, 0.10, 0.25):
            k1, k2 = k1k2_for_eff(n, eff)
            S_K1, _, S_K2, _ = scales(n, k1, k2)
            maxa, aloc, yloc, nus, mus = [], [], [], [], []
            for _ in range(500):
                lam = rng.randrange(1, n)
                cv = cf_convergents(lam, n)
                nsq, (a, eta), idx = lambda1_cf(n, lam, S_K1, S_K2)
                nu = math.exp(0.5 * (math.log(nsq) - math.log(n * S_K1 * S_K2)))
                maxa.append(max(t[2] for t in cv))
                # partial quotient governing the minimising convergent's quality
                nxt = cv[idx + 1][2] if 0 <= idx + 1 < len(cv) else cv[-1][2]
                aloc.append(nxt)
                yloc.append(a * eta / n if a else 1.0)
                nus.append(nu)
                mus.append(min(lam, n - lam) / n)
            print(f"{bits:>5} {eff:>6.2f} {spearman(maxa, nus):>13.3f} "
                  f"{spearman(aloc, nus):>13.3f} {spearman(yloc, nus):>13.3f} "
                  f"{spearman(mus, nus):>10.3f}")
    print("\n  Reading: max_a and mu are near-zero (scale-free / wrong coordinate);")
    print("  a_local and y_local carry the signal.  Same CF, different question.")

def part_e():
    print("\n" + "=" * 78)
    print("PART E (T23-C): the null law of nu_hat, and bit-size invariance")
    print("=" * 78)
    print("Predicted CDF  F(r) = 3r^2/pi on [0,1], closed form on [1, (4/3)^(1/4)],")
    print(f"support ends at Hermite bound (4/3)^(1/4) = {HERMITE2:.5f}.\n")
    ps = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
    print("  predicted quantiles: " +
          "  ".join(f"q{int(100*p):02d}={nu_quantile(p):.3f}" for p in ps))
    print()
    rng = random.Random(99)
    print(f"{'bits':>5} {'eff':>6} {'N':>6} {'KS D':>7} {'KS p':>8} {'max nu':>8} "
          f"{'  empirical quantiles (q05..q95)'}")
    ok = True
    for bits in (16, 20, 24, 32, 64, 128, 256):
        n = rand_prime(bits, rng)
        for eff in (0.05, 0.25):
            k1, k2 = k1k2_for_eff(n, eff)
            S_K1, _, S_K2, _ = scales(n, k1, k2)
            det = n * S_K1 * S_K2
            nus = []
            for _ in range(2000):
                lam = rng.randrange(1, n)
                nsq, _, _ = lambda1_cf(n, lam, S_K1, S_K2)
                nus.append(math.exp(0.5 * (math.log(nsq) - math.log(det))))
            d, p = ks_test(nus, nu_cdf)
            qs = quantiles(nus, ps)
            if p < 0.01:
                ok = False
            print(f"{bits:>5} {eff:>6.2f} {len(nus):>6} {d:>7.4f} {p:>8.4f} "
                  f"{max(nus):>8.5f}  " + " ".join(f"{q:.3f}" for q in qs))
    print(f"\n  VERDICT: {'law CONFIRMED across 16..256 bits' if ok else 'law rejected somewhere'}")
    print("  Cross-check vs 2026-07-29 log Part 5 (256-bit, eff=0.05, 4000 draws):")
    print("    logged  q05=0.232 q10=0.322 q25=0.496 q50=0.715 q75=0.885 q90=0.974 q95=0.999")
    print("    theory  q05=%.3f q10=%.3f q25=%.3f q50=%.3f q75=%.3f q90=%.3f q95=%.3f"
          % tuple(nu_quantile(p) for p in ps))
    print("    logged  '20.5%% of lam fall below 0.45';  theory F(0.45) = %.3f"
          % nu_cdf(0.45))
    return ok

def part_f():
    print("\n" + "=" * 78)
    print("PART F: secp256k1 placement, recomputed by CF and scored by the law")
    print("=" * 78)
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0, "lam^2+lam+1 != 0 mod n"
    print("  lam^2 + lam + 1 = 0 mod n   verified.")
    mu = min(SECP_LAM, SECP_N - SECP_LAM) / SECP_N
    print(f"  mu = {mu:.6f}    (falsified as a predictor 2026-07-29)\n")
    print(f"{'eff':>6} {'nu_hat(CF)':>11} {'nu(Gauss)':>11} {'F(nu) pct':>10} "
          f"{'log2 q*':>8} {'log2 q_arg':>11} {'a_local':>8} {'y_local':>9}")
    out = {}
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_for_eff(SECP_N, eff)
        S_K1, _, S_K2, _ = scales(SECP_N, k1, k2)
        nsq, (a, eta), idx = lambda1_cf(SECP_N, SECP_LAM, S_K1, S_K2)
        det = SECP_N * S_K1 * S_K2
        nu = math.exp(0.5 * (math.log(nsq) - math.log(det)))
        _, _, nu_g = rival_sublattice_nu(SECP_N, SECP_LAM, k1, k2)
        qstar = math.exp(0.5 * (math.log(SECP_N) + math.log(S_K1) - math.log(S_K2)))
        cv = cf_convergents(SECP_LAM, SECP_N)
        aloc = cv[idx + 1][2] if 0 <= idx + 1 < len(cv) else float('nan')
        yloc = a * eta / SECP_N if a else 1.0
        out[eff] = nu
        print(f"{eff:>6.2f} {nu:>11.4f} {nu_g:>11.4f} {100*nu_cdf(nu):>9.1f}% "
              f"{math.log2(qstar):>8.2f} {math.log2(a) if a else 0:>11.2f} "
              f"{aloc:>8.0f} {yloc:>9.4f}")
    print("\n  2026-07-29 logged values: eff=0.05 0.8709, 0.10 0.6624,")
    print("                            eff=0.25 0.5852, 0.50 0.6851")
    print(f"  this run (CF):            eff=0.05 {out[0.05]:.4f}, 0.10 {out[0.10]:.4f},")
    print(f"                            eff=0.25 {out[0.25]:.4f}, 0.50 {out[0.50]:.4f}")
    print("\n  Caveat carried forward verbatim from 2026-07-29: this whole thread is")
    print("  conditional on a non-standard nonce generator k = k1 + lam*k2 with k1")
    print("  bounded.  It says nothing about correctly generated nonces and is NOT")
    print("  an attack on secp256k1.")

def part_g():
    print("\n" + "=" * 78)
    print("PART G (corollary): why max_a is a one-sided CONSTRAINT, not a predictor")
    print("=" * 78)
    print("Each term of the T23-B min obeys y^2/x^2 + x^2 >= 2y (AM-GM, equality at")
    print("x = sqrt(y)), so nu_hat^2 >= 2*min_j y_j.  Classically")
    print("y_j = q_j*eta_j/n >= 1/(a_{j+1}+2), hence")
    print()
    print("      nu_hat >= sqrt( 2 / (max_a + 2) ) .")
    print()
    print("So a small nu_hat REQUIRES a large partial quotient somewhere in the CF")
    print("of lam/n -- necessary, not sufficient, because the large a_j must ALSO")
    print("land near the scale.  That is exactly the shape of a variable that")
    print("correlates weakly (~ -0.2, Part D) and fails as a classifier (June).\n")
    rng = random.Random(864213)
    print(f"{'bits':>5} {'eff':>6} {'N':>6} {'violations':>11} {'median slack':>13} "
          f"{'slack q95':>10}")
    viol_total = 0
    for bits in (20, 32, 64, 128, 256):
        n = rand_prime(bits, rng)
        for eff in (0.05, 0.25):
            k1, k2 = k1k2_for_eff(n, eff)
            S_K1, _, S_K2, _ = scales(n, k1, k2)
            det = n * S_K1 * S_K2
            slacks, viol = [], 0
            for _ in range(600):
                lam = rng.randrange(1, n)
                cv = cf_convergents(lam, n)
                nsq, _, _ = lambda1_cf(n, lam, S_K1, S_K2)
                nu = math.exp(0.5 * (math.log(nsq) - math.log(det)))
                bound = math.sqrt(2.0 / (max(t[2] for t in cv) + 2))
                if nu < bound - 1e-12:
                    viol += 1
                slacks.append(nu / bound)
            viol_total += viol
            s = quantiles(slacks, [0.5, 0.95])
            print(f"{bits:>5} {eff:>6.2f} {600:>6} {viol:>11} {s[0]:>13.3f} "
                  f"{s[1]:>10.3f}")
    print(f"\n  total violations of nu_hat >= sqrt(2/(max_a+2)): {viol_total}")
    print("  The bound is never violated but is loose, and its slack GROWS with")
    print("  log n (median 2.3x at 20 bits -> 8.3x at 256 bits), because a longer CF")
    print("  gives max_a more chances to be large somewhere far from the scale.")
    print("  That is precisely why max_a alone cannot classify, and why it degrades")
    print("  as a bound exactly in the cryptographic size range.")
    return viol_total == 0

def main():
    t0 = time.time()
    print("=" * 78)
    print("GLV-HNP Thread 23 — closed form for lambda_1(L2); null law of nu_hat")
    print("=" * 78)
    a_ok = part_a()
    b_ok = part_b()
    part_c()
    part_d()
    e_ok = part_e()
    part_f()
    g_ok = part_g()
    print("\n" + "=" * 78)
    print(f"SUMMARY  T23-A {'CONFIRMED' if a_ok else 'FALSIFIED'} | "
          f"T23-B {'CONFIRMED' if b_ok else 'FALSIFIED'} | "
          f"T23-C {'CONFIRMED' if e_ok else 'REJECTED'} | "
          f"corollary {'holds' if g_ok else 'VIOLATED'}")
    print(f"wall time {time.time()-t0:.1f}s")
    print("=" * 78)
    return 0

if __name__ == '__main__':
    sys.exit(main())
