"""
GLV-HNP Phase 2 — Thread 23: closed form for lambda_1(L2) via continued fractions.

Background
----------
Thread 20b/20c (2026-07-29 run #2) found that

    nu_hat = lambda_1(L2) / sqrt(det L2),
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2

separates the June C1/C2 classes with AUC 0.935 -- the first curve-level
predictor to beat baseline after eight falsified invariants.  But nu_hat was
*computed* (one Lagrange-Gauss reduction), not *understood*.

Thread 23 hypothesis (H-CF)
---------------------------
A general vector of L2 is  a*u + b*v = ( (a*n - b*lam)*S_K1 , b*S_K2 ).
For fixed b the first coordinate is minimised by the nearest multiple, giving
r_b = |b*lam mod± n|.  Hence

    lambda_1(L2)^2 = min over b >= 0 of  F(b) = S_K1^2 * r_b^2 + S_K2^2 * b^2

with the b = 0 term contributing (n*S_K1)^2.  F is increasing in both b and
r_b, so any minimiser must be Pareto-optimal in (b, r_b): there is no b' < b
with r_{b'} <= r_b.  That is exactly the definition of a *best approximation of
the second kind* to lam/n, and those are exactly the continued-fraction
convergent denominators q_j.  Therefore

    H-CF:  lambda_1(L2) = min over convergents q_j of sqrt(F(q_j)),

i.e. nu_hat is a continued-fraction quantity evaluated at the balance scale
    b* = sqrt(n * S_K1 / S_K2)  ~  sqrt(n * K2 / K1),
which is where S_K1*r_b and S_K2*b cross (using r_b ~ n/b).

This is the predicted resolution of the June puzzle: the raw CF invariants
(q_cf, max_q_cf, max_a) all failed because they are scale-free, whereas the
relevant approximation quality is scale-dependent -- it is the quality of the
*single* convergent nearest b*, not of the CF as a whole.

Falsifier (as stated 2026-07-29): if the CF-predicted nu_hat correlates < 0.9
with the computed nu_hat, the convergent-scale story is wrong.  We test the
much sharper claim of EXACT equality, in exact integer arithmetic.

Experiments
-----------
E1  Exactness: CF closed form vs exact-integer Lagrange-Gauss, random (n,lam)
    over a range of bit sizes.  Also checks whether semiconvergents are ever
    needed (H-CF says no).
E2  Balance scale: is the winning convergent q_{j*} the one nearest
    b* = sqrt(n*S_K1/S_K2)?  Reports |log2(q_{j*}/b*)|.
E3  Float-precision audit of the shipped gauss_reduce_2d(), which divides two
    Python ints via float `/`.  At 256 bits <u,v> and |v|^2 exceed 2^1024, so
    the division can overflow to inf/NaN.  Determines the bit size at which the
    shipped routine stops being trustworthy -- this bears directly on the
    secp256k1 nu_hat numbers reported on 2026-07-29.
E4  secp256k1 nu_hat recomputed exactly via the CF closed form, at the same
    eff values reported on 2026-07-29, to confirm or correct those figures.

Run: python3 glv_hnp_nuhat_cf_closed_form.py
"""

import math
import random
import sys

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "_nuhat", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_nuhat = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nuhat)

scales = _nuhat.scales
gauss_reduce_2d_float = _nuhat.gauss_reduce_2d

# secp256k1
P_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N_SECP = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


# ---------------------------------------------------------------------------
# Exact integer primitives
# ---------------------------------------------------------------------------

def idiv_round(a, b):
    """Round-half-away nearest integer of a/b, exact for arbitrary ints."""
    assert b > 0
    q, r = divmod(a, b)
    if 2 * r >= b:
        q += 1
    return q


def gauss_reduce_2d_exact(u, v):
    """Lagrange-Gauss reduction in exact integer arithmetic.

    Returns lambda_1(L)^2 as an exact int (squared norm, no float anywhere).
    """
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]

    if n2(u) < n2(v):
        u, v = v, u
    while True:
        nv = n2(v)
        if nv == 0:
            return n2(u)
        q = idiv_round(u[0] * v[0] + u[1] * v[1], nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if n2(u) <= n2(v):
            return n2(u)


def convergent_denominators(lam, n):
    """Denominators k_j of the continued-fraction convergents of lam/n,
    together with the partial quotients a_j.  Exact integer arithmetic.

    Recurrence with k_{-2} = 1, k_{-1} = 0 and k_j = a_j*k_{j-1} + k_{j-2},
    so that k_0 = 1 (convergent a_0/1) and k_1 = a_1.
    """
    a, b = lam, n
    qs, quots = [], []
    k_pp, k_p = 1, 0          # k_{-2}, k_{-1}
    while b:
        c = a // b
        a, b = b, a - c * b
        quots.append(c)
        k_pp, k_p = k_p, c * k_p + k_pp
        qs.append(k_p)
    return qs, quots


def r_of(b, lam, n):
    """|b*lam mod± n| -- distance from b*lam to the nearest multiple of n."""
    t = (b * lam) % n
    return min(t, n - t)


def F(b, lam, n, s_k1, s_k2):
    """Squared norm of the shortest L2 vector with second coordinate b*S_K2."""
    return s_k1 * s_k1 * r_of(b, lam, n) ** 2 + s_k2 * s_k2 * b * b


def nu_hat_cf(n, lam, k1_bound, k2_bound, with_semiconvergents=False):
    """CF closed form for nu_hat.  Returns (nu_hat, best_b, cand_count)."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    det = n * s_k1 * s_k2

    cands = [0]
    qs, quots = convergent_denominators(lam, n)
    cands.extend(q for q in qs if 0 < q <= n)
    if with_semiconvergents:
        k_pp, k_p = 1, 0
        for c in quots:
            for t in range(1, c):
                sc = t * k_p + k_pp
                if 0 < sc <= n:
                    cands.append(sc)
            k_pp, k_p = k_p, c * k_p + k_pp

    best_f, best_b = None, None
    for b in set(cands):
        f = (n * s_k1) ** 2 if b == 0 else F(b, lam, n, s_k1, s_k2)
        if best_f is None or f < best_f:
            best_f, best_b = f, b
    return math.sqrt(best_f / det), best_b, len(set(cands))


def nu_hat_exact(n, lam, k1_bound, k2_bound):
    """Ground truth: exact-integer Lagrange-Gauss on L2."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    det = n * s_k1 * s_k2
    lam1_sq = gauss_reduce_2d_exact((n * s_k1, 0), (-lam * s_k1, s_k2))
    return math.sqrt(lam1_sq / det), lam1_sq


def nu_hat_bruteforce(n, lam, k1_bound, k2_bound, cap=200000):
    """Exhaustive min over b < cap.  Only for small n, as an independent check."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    det = n * s_k1 * s_k2
    best = (n * s_k1) ** 2
    best_b = 0
    for b in range(1, min(n, cap)):
        f = F(b, lam, n, s_k1, s_k2)
        if f < best:
            best, best_b = f, b
    return math.sqrt(best / det), best_b


# ---------------------------------------------------------------------------
# Curve generation: j=0 GLV curves (lam^2 + lam + 1 = 0 mod n)
# ---------------------------------------------------------------------------

def glv_lambdas(n):
    """Both roots of x^2+x+1 mod n, or None."""
    # need sqrt(-3) mod n
    d = n - 3
    s = pow(d, (n + 1) // 4, n) if n % 4 == 3 else _tonelli(d, n)
    if s is None or (s * s - d) % n != 0:
        return None
    inv2 = pow(2, n - 2, n)
    r1 = (-1 + s) * inv2 % n
    r2 = (-1 - s) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return r1, r2


def _tonelli(a, p):
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


def random_glv_pairs(bits, count, rng):
    """(n, lam) pairs with n prime of the given bit size and lam a GLV root."""
    import sympy
    out = []
    tries = 0
    while len(out) < count and tries < 40000:
        tries += 1
        n = sympy.nextprime(rng.getrandbits(bits) | (1 << (bits - 1)))
        if n % 3 != 1:
            continue
        ls = glv_lambdas(n)
        if ls is None:
            continue
        out.append((n, min(ls)))
    return out


def k1k2_for_eff(n, eff):
    """BALANCED split: K1 ~ K2 ~ sqrt(eff*n), so K1*K2/n ~ eff."""
    prod = eff * n
    k1 = max(2, int(math.isqrt(int(prod))))
    k2 = max(2, int(prod) // k1)
    return k1, k2


def k1k2_thread20(n, eff):
    """The Thread 20b/20c split (glv_hnp_phase2_nuhat_control.py:88):
    K2 = sqrt(n) fixed, K1 = eff*sqrt(n).  Also gives K1*K2/n ~ eff, but a
    very different S_K1/S_K2 ratio -- and nu_hat depends on that ratio, not on
    eff alone.  Kept so the 2026-07-29 figures can be reproduced exactly."""
    return max(2, int(eff * math.isqrt(n))), math.isqrt(n) + 1


# ---------------------------------------------------------------------------
# E1 -- exactness of the CF closed form
# ---------------------------------------------------------------------------

def e1(rng):
    print("=" * 74)
    print("E1  CF closed form vs exact-integer Lagrange-Gauss")
    print("=" * 74)
    print("H-CF: lambda_1(L2) = min over CF convergent denominators q_j of F(q_j)")
    print()
    print(f"{'bits':>5} {'curves':>7} {'eff':>6} {'exact match':>12} {'max rel err':>12} "
          f"{'semiconv needed':>16}")
    total, matched, semi_needed = 0, 0, 0
    worst = 0.0
    for bits in (16, 20, 24, 32, 48, 64, 96, 128, 192, 256):
        pairs = random_glv_pairs(bits, 12, rng)
        for eff in (0.05, 0.10, 0.25):
            ok, bad = 0, 0.0
            for (n, lam) in pairs:
                k1, k2 = k1k2_for_eff(n, eff)
                nu_e, _ = nu_hat_exact(n, lam, k1, k2)
                nu_c, b_c, _ = nu_hat_cf(n, lam, k1, k2)
                nu_s, b_s, _ = nu_hat_cf(n, lam, k1, k2, with_semiconvergents=True)
                total += 1
                rel = abs(nu_c - nu_e) / max(nu_e, 1e-300)
                if rel < 1e-12:
                    ok += 1
                    matched += 1
                bad = max(bad, rel)
                if abs(nu_s - nu_c) > 1e-12 * max(nu_c, 1e-300):
                    semi_needed += 1
            worst = max(worst, bad)
            print(f"{bits:>5} {len(pairs):>7} {eff:>6.2f} {ok:>7}/{len(pairs):<4} "
                  f"{bad:>12.2e} {'-':>16}")
    print()
    print(f"  TOTAL exact matches: {matched}/{total}   max rel err {worst:.3e}")
    print(f"  cases where semiconvergents changed the answer: {semi_needed}/{total}")
    print(f"  VERDICT: {'H-CF CONFIRMED' if matched == total else 'H-CF FALSIFIED'}")
    return matched == total, total


def e1b(rng):
    """Independent brute-force check on small n (no CF, no Gauss)."""
    print()
    print("-" * 74)
    print("E1b  brute force over all b < n  (small n, independent of both methods)")
    print("-" * 74)
    pairs = random_glv_pairs(16, 8, rng)
    agree = 0
    for (n, lam) in pairs:
        k1, k2 = k1k2_for_eff(n, 0.10)
        nu_b, b_b = nu_hat_bruteforce(n, lam, k1, k2)
        nu_c, b_c, _ = nu_hat_cf(n, lam, k1, k2)
        nu_e, _ = nu_hat_exact(n, lam, k1, k2)
        same = abs(nu_b - nu_c) < 1e-12 and abs(nu_b - nu_e) < 1e-12
        agree += same
        print(f"  n={n:<7} lam={lam:<7} brute nu={nu_b:.6f}(b={b_b:<6}) "
              f"cf nu={nu_c:.6f}(b={b_c:<6}) gauss nu={nu_e:.6f}  {'OK' if same else 'MISMATCH'}")
    print(f"  {agree}/{len(pairs)} agree")
    return agree == len(pairs)


# ---------------------------------------------------------------------------
# E2 -- the balance scale
# ---------------------------------------------------------------------------

def e2(rng):
    print()
    print("=" * 74)
    print("E2  is the winning convergent the one nearest b* = sqrt(n*S_K1/S_K2)?")
    print("=" * 74)
    print("  Using the Thread 20 split (K2=sqrt(n), K1=eff*sqrt(n)) so that")
    print("  b* = sqrt(n/eff) actually MOVES with eff.  Under the balanced split")
    print("  S_K1/S_K2 = 1 for every eff, so b* = sqrt(n) is constant and the")
    print("  question is degenerate.")
    print()
    print(f"{'bits':>5} {'eff':>6} {'n':>5} {'|log2(q*/b*)| mean':>20} {'max':>8} "
          f"{'rank of q* by |log2|':>22}")
    for bits in (24, 64, 128, 256):
        for eff in (0.05, 0.25):
            pairs = random_glv_pairs(bits, 15, rng)
            devs, ranks = [], []
            for (n, lam) in pairs:
                k1, k2 = k1k2_thread20(n, eff)
                s_k1, _, s_k2, _ = scales(n, k1, k2)
                b_star = math.sqrt(n * s_k1 / s_k2)
                _, b_win, _ = nu_hat_cf(n, lam, k1, k2)
                if b_win == 0:
                    continue
                devs.append(abs(math.log2(b_win / b_star)))
                qs, _ = convergent_denominators(lam, n)
                qs = sorted(set(q for q in qs if 0 < q < n),
                            key=lambda q: abs(math.log2(q / b_star)))
                ranks.append(qs.index(b_win) if b_win in qs else -1)
            if not devs:
                continue
            r0 = sum(1 for r in ranks if r == 0)
            r01 = sum(1 for r in ranks if r in (0, 1))
            print(f"{bits:>5} {eff:>6.2f} {len(devs):>5} {sum(devs)/len(devs):>20.3f} "
                  f"{max(devs):>8.3f} "
                  f"{'rank0 ' + str(r0) + '/' + str(len(ranks)) + ', rank<=1 ' + str(r01):>22}")
    print()
    print("  b* is where S_K1*r_b and S_K2*b cross (r_b ~ n/b).  A small |log2|")
    print("  deviation means nu_hat is set by ONE convergent, the one at that scale.")


# ---------------------------------------------------------------------------
# E3 -- float precision audit of the shipped gauss_reduce_2d
# ---------------------------------------------------------------------------

def e3(rng):
    print()
    print("=" * 74)
    print("E3  float-precision audit of shipped gauss_reduce_2d()")
    print("=" * 74)
    print("  The shipped routine computes  q = round(<u,v> / |v|^2)  with float")
    print("  division of Python ints.  Both operands grow like (n*S_K1)^2 ~ n^4.")
    print()
    print(f"{'bits':>5} {'n^4 bits':>9} {'trials':>7} {'float==exact':>13} "
          f"{'overflow/NaN':>13} {'max rel err':>12}")
    for bits in (16, 24, 32, 64, 128, 192, 256, 384, 521):
        pairs = random_glv_pairs(bits, 8, rng)
        ok, blew, worst = 0, 0, 0.0
        for (n, lam) in pairs:
            k1, k2 = k1k2_for_eff(n, 0.10)
            nu_e, _ = nu_hat_exact(n, lam, k1, k2)
            try:
                s_k1, _, s_k2, _ = scales(n, k1, k2)
                det = n * s_k1 * s_k2
                nu_f = gauss_reduce_2d_float((n * s_k1, 0), (-lam * s_k1, s_k2)) \
                    / math.sqrt(det)
                if not math.isfinite(nu_f):
                    blew += 1
                    continue
                rel = abs(nu_f - nu_e) / max(nu_e, 1e-300)
                worst = max(worst, rel)
                if rel < 1e-9:
                    ok += 1
            except (OverflowError, ZeroDivisionError, ValueError):
                blew += 1
        nb = 4 * bits
        print(f"{bits:>5} {nb:>9} {len(pairs):>7} {ok:>8}/{len(pairs):<4} "
              f"{blew:>13} {worst:>12.2e}")
    print()
    print("  Diagnosis of the failure mode:")
    n, lam = random_glv_pairs(521, 1, rng)[0]
    k1, k2 = k1k2_for_eff(n, 0.10)
    s_k1, _, s_k2, _ = scales(n, k1, k2)
    u = (n * s_k1, 0)
    big = u[0] * u[0]
    print(f"    |u|^2 has {big.bit_length()} bits (float64 overflows above 1024)")
    for label, fn in (("int/int true division  <u,v>/|v|^2",
                       lambda: (u[0] * u[0] + u[1] * 0) / big),
                      ("math.sqrt(|u|^2)", lambda: math.sqrt(big))):
        try:
            v = fn()
            print(f"    {label:<38} -> {v!r}")
        except (OverflowError, ValueError) as e:
            print(f"    {label:<38} -> {type(e).__name__}: {e}")
    print("    => Python's int/int division is correctly rounded and does NOT")
    print("       overflow; the failure is the final math.sqrt() on a >1024-bit int.")
    print("    RELEVANCE TO PRIORITY THREAD 1 (P-521 LLL NaN): the same pattern --")
    print("    math.sqrt / float() applied to an exact big int -- is a candidate")
    print("    cause of NaN at 521 bits that is NOT a Gram-Schmidt precision problem")
    print("    and needs no bigfloat library to fix (use math.isqrt, or rescale by")
    print("    bit_length, instead).  Worth checking before adding a rug dependency.")


# ---------------------------------------------------------------------------
# E4 -- secp256k1 recomputed exactly
# ---------------------------------------------------------------------------

def e4():
    print()
    print("=" * 74)
    print("E4  secp256k1 nu_hat, exact CF closed form vs the 2026-07-29 figures")
    print("=" * 74)
    n = N_SECP
    ls = glv_lambdas(n)
    assert ls is not None
    lam = min(ls)
    assert (lam * lam + lam + 1) % n == 0
    print(f"  n   = {n}")
    print(f"  lam = {lam}   (lam^2+lam+1 = 0 mod n verified)")
    print(f"  lam* = min(lam,n-lam)/n = {min(lam, n-lam)/n:.6f}")
    print()
    reported = {0.05: 0.8709, 0.10: 0.6624, 0.25: 0.5852, 0.50: 0.6851}
    print("  (a) Thread 20b/20c split  K2=sqrt(n), K1=eff*sqrt(n)  -- reproduces 2026-07-29")
    print(f"{'eff':>6} {'S_K1/S_K2':>11} {'nu_hat (exact CF)':>19} {'2026-07-29':>12} "
          f"{'b_win/b*':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_thread20(n, eff)
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        nu_c, b_win, _ = nu_hat_cf(n, lam, k1, k2)
        nu_e, _ = nu_hat_exact(n, lam, k1, k2)
        assert abs(nu_c - nu_e) < 1e-12, (nu_c, nu_e)
        b_star = math.sqrt(n * s_k1 / s_k2)
        print(f"{eff:>6.2f} {s_k1/s_k2:>11.2f} {nu_c:>19.4f} {reported[eff]:>12.4f} "
              f"{b_win/b_star:>10.4f}")
    print()
    print("  (b) BALANCED split  K1 ~ K2 ~ sqrt(eff*n)  -- same eff, different nu_hat")
    print(f"{'eff':>6} {'S_K1/S_K2':>11} {'nu_hat (exact CF)':>19} {'b_win/b*':>10}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_for_eff(n, eff)
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        nu_c, b_win, _ = nu_hat_cf(n, lam, k1, k2)
        b_star = math.sqrt(n * s_k1 / s_k2)
        print(f"{eff:>6.2f} {s_k1/s_k2:>11.2f} {nu_c:>19.4f} {b_win/b_star:>10.4f}")
    print()
    print("  nu_hat is a function of the RATIO S_K1/S_K2 (i.e. of b*), not of eff.")
    print("  Two (K1,K2) splits with identical eff give different nu_hat.  Any future")
    print("  use of nu_hat must therefore quote the split, not just eff.")
    print()
    print("  NOTE: this is a lattice invariant of (n, lam), not an attack.  It is")
    print("  conditional on the non-standard nonce generator k = k1 + lam*k2 with")
    print("  k1 bounded; it says nothing about correctly generated nonces.")
    return lam


def e4b(lam, rng):
    """Null distribution of nu_hat over random lam at 256 bits, exact."""
    print()
    print("-" * 74)
    print("E4b  null distribution of nu_hat over random lam (256-bit n_secp), exact")
    print("-" * 74)
    print("  Thread 20 split, for comparability with the 2026-07-29 null table.")
    n = N_SECP
    for eff in (0.05, 0.10):
        k1, k2 = k1k2_thread20(n, eff)
        vals = []
        for _ in range(4000):
            l = rng.randrange(1, n)
            nu, _, _ = nu_hat_cf(n, l, k1, k2)
            vals.append(nu)
        vals.sort()
        def q(p):
            return vals[int(p * (len(vals) - 1))]
        nu_s, _, _ = nu_hat_cf(n, lam, k1, k2)
        pct = 100.0 * sum(1 for v in vals if v < nu_s) / len(vals)
        print(f"  eff={eff:.2f}  q05={q(.05):.3f} q25={q(.25):.3f} q50={q(.50):.3f} "
              f"q75={q(.75):.3f} q95={q(.95):.3f}")
        print(f"            secp256k1 nu_hat={nu_s:.4f} -> percentile {pct:.1f}%")


def main():
    rng = random.Random(20260801)
    print("GLV-HNP Phase 2 — Thread 23: closed form for lambda_1(L2)")
    print("2026-08-01 autolab run")
    print()
    ok1, ntot = e1(rng)
    ok1b = e1b(rng)
    e2(rng)
    e3(rng)
    lam = e4()
    e4b(lam, rng)
    print()
    print("=" * 74)
    print(f"SUMMARY: H-CF {'CONFIRMED' if ok1 and ok1b else 'FALSIFIED'} "
          f"over {ntot} (n,lam,eff) cases + brute-force cross-check")
    print("=" * 74)


if __name__ == '__main__':
    main()
