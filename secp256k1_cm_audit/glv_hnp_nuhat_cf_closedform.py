#!/usr/bin/env python3
"""
GLV-HNP Phase 2 — Thread 23: closed form for lambda_1(L2), and why the June
continued-fraction invariants failed.

Background (2026-07-29 run #2, commit e845207).  The Phase 2 lattice contains m
copies of a non-planted 2-D sublattice

    L2 = < u = (n*S_K1, 0),  v = (-lam*S_K1, S_K2) >,     det L2 = n*S_K1*S_K2

and nu_hat = lambda_1(L2)/sqrt(det L2) separates the June C1/C2 classes at
AUC 0.935 where six curve-level algebraic invariants had all failed.  nu_hat was
*computed* (one Lagrange-Gauss reduction) but never *understood*.

This script tests the conjecture recorded as the Thread 23 next-step: that
lambda_1(L2) is a best-rational-approximation quantity for lam/n taken at the
scale set by (K1, K2), and that the relevant convergent is the one with
denominator near

    q* = sqrt(n*S_K1/S_K2)  ~  sqrt(n*K2/K1).

Parametrisation.  A general lattice vector is

    a*v + b*u = ( (b*n - a*lam)*S_K1 , a*S_K2 )

so with r(a) = |a*lam mod n| (centred residue, i.e. min over b),

    lambda_1(L2)^2 = min_{a=0..n-1, a or r(a) nonzero} [ S_K1^2*r(a)^2 + S_K2^2*a^2 ]

with a = 0 contributing the trivial vector of norm n*S_K1.  Minimising a
weighted sum of |a*lam mod n| and a is exactly best rational approximation of
lam/n *at a scale*, which is the conjecture to be tested.

Experiments
  E0  exact-integer Lagrange-Gauss vs the float version used on 2026-07-29.
      The float version rounds q = round(<u,v>/|v|^2) in double precision; at
      256 bits the Gram entries are ~2^1012, so this is an audit of the
      previously reported secp256k1 nu_hat, not a formality.
  E1  identify the exact minimiser a* and test membership in
      {convergent denominators} and {semiconvergent denominators} of lam/n.
  E2  single-convergent predictor: the convergent nearest q*.  Falsifier from
      the Thread 23 proposal: corr(predicted nu_hat, computed nu_hat) < 0.9
      kills the convergent-scale story.
  E3  min over ALL convergents (separates "right set" from "right selector").
  E4  min over convergents + semiconvergents.
  E5  scale-dependence: fix (n, lam), sweep K1/K2, watch which convergent wins.
      This is the retro-explanation of the six June failures.
  E6  secp256k1: which convergent index realises nu_hat, at each eff.
  E7  selector hit-rate + dim-2 Hermite validation of every predictor.

OUTCOME (2026-07-29 run #3).  The conjecture splits cleanly.

  CONFIRMED -- the SET is right, and provably so.  The minimiser a* is always a
  convergent denominator of lam/n: 450/450 cells over 20/24/32/64/128/256 bits,
  cross-checked against exhaustive search at 20 bits (160/160).  This is a
  theorem, not a regularity.  Exchange argument: if a minimises
  S_K1^2*r(a)^2 + S_K2^2*a^2 but is not a best-approximation denominator of the
  second kind, there is a' < a with r(a') <= r(a), and then both terms are
  strictly smaller -- contradiction.  Best-approximation denominators of the
  second kind are exactly the CF convergent denominators.  Hence

      lambda_1(L2)^2 = min over convergent denominators k_j of lam/n of
                       [ S_K1^2*|k_j*lam mod n|_c^2 + S_K2^2*k_j^2 ]

  and semiconvergents are never needed (E4 == E3 in every cell).

  FALSIFIED -- the SELECTOR is wrong.  "The convergent nearest
  q* = sqrt(n*S_K1/S_K2)" does NOT pick the winner.  It picks the correct
  convergent index in only 83/150 (20b), 89/150 (24b), 36/75 (64b), 42/75
  (256b) of cells, and corr(E2, exact) is +0.036 / +0.202 / -0.220 / -0.168 --
  far under the 0.9 the Thread 23 proposal set as its own falsifier.  Its
  predicted nu_hat exceeds the dim-2 Hermite bound (2/sqrt3)^(1/2) = 1.07457 in
  63/150, 57/150, 37/75, 32/75 cells, i.e. it is provably not a lattice minimum
  there.  The exact values never violate the bound (0/450), which validates the
  reduction independently.

  This costs nothing: the min over all convergents is still one Euclid run,
  O(log n), the same complexity as the Lagrange-Gauss reduction it replaces.
  What is gained is the algebraic form, and with it E5's retro-explanation --
  the winning convergent index MOVES with (K1, K2) (e.g. [8, 10, 12] across eff
  on one fixed curve), so no scale-free function of the CF of lam/n can predict
  nu_hat.  That is the mechanism behind the six invariants falsified in June
  (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) and mu on 2026-07-29:
  q_cf, max_q_cf and max_a are precisely scale-free CF statistics.

Run:  python3 secp256k1_cm_audit/glv_hnp_nuhat_cf_closedform.py
No dependencies beyond the standard library.  Runtime ~7 min.
"""

import math
import random
import sys

# --- reuse the exact scale convention from Thread 20 -------------------------
# scales(n, K1, K2) = (S_K1, S_D, S_K2, S_KANNAN) = (n//K1, 1, max(1,n//K2), n)
# Only S_K1 and S_K2 enter L2.


def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)


def k1k2_project(n, eff):
    """
    The (K1, K2) convention used throughout Thread 20
    (glv_hnp_phase2_nuhat_control.py:88 and ..._threshold_r2.py:408):
        K2 = isqrt(n) + 1,   K1 = eff * sqrt(n),   so K1*K2 ~ eff*n.
    Consequently S_K1/S_K2 ~ 1/eff and q* = sqrt(n*S_K1/S_K2) ~ sqrt(n/eff),
    i.e. eff moves the scale.  A balanced split K1 ~ K2 ~ sqrt(eff*n) instead
    makes S_K1 ~ S_K2 and nu_hat becomes eff-independent -- do not use it.
    """
    return max(2, int(eff * math.sqrt(n))), math.isqrt(n) + 1


# ---------------------------------------------------------------------------
# Lagrange-Gauss reduction
# ---------------------------------------------------------------------------

def _nrm2(w):
    return w[0] * w[0] + w[1] * w[1]


def _iround_div(a, b):
    """round(a/b) for integers, exact, b > 0, ties away from zero."""
    if a >= 0:
        return (2 * a + b) // (2 * b)
    return -((-2 * a + b) // (2 * b))


def gauss_reduce_exact(u, v):
    """Exact-integer Lagrange-Gauss. Returns (shortest_vector, its norm^2)."""
    u, v = tuple(u), tuple(v)
    if _nrm2(u) < _nrm2(v):
        u, v = v, u
    while True:
        nv = _nrm2(v)
        if nv == 0:
            return u, _nrm2(u)
        q = _iround_div(u[0] * v[0] + u[1] * v[1], nv)
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if _nrm2(u) <= _nrm2(v):
            return u, _nrm2(u)


def gauss_reduce_float(u, v):
    """Verbatim copy of glv_hnp_phase2_lambda_threshold_r2.py:216 (float rounding)."""
    if _nrm2(u) < _nrm2(v):
        u, v = v, u
    while True:
        q = round((u[0] * v[0] + u[1] * v[1]) / _nrm2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if _nrm2(u) <= _nrm2(v):
            return math.sqrt(_nrm2(u))


def L2_basis(n, lam, k1_bound, k2_bound):
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    return (n * s_k1, 0), (-lam * s_k1, s_k2), s_k1, s_k2


def nu_hat_exact(n, lam, k1_bound, k2_bound):
    """(nu_hat, minimiser a*, lambda_1^2, det) with exact integer arithmetic."""
    u, v, s_k1, s_k2 = L2_basis(n, lam, k1_bound, k2_bound)
    w, nrm2 = gauss_reduce_exact(u, v)
    det = n * s_k1 * s_k2
    a_star = abs(w[1]) // s_k2
    nu = math.exp(0.5 * (math.log(nrm2) - math.log(det)))
    return nu, a_star, nrm2, det


# ---------------------------------------------------------------------------
# Continued fractions of lam/n
# ---------------------------------------------------------------------------

def cf_convergents(lam, n):
    """
    Convergents p_j/q_j of lam/n.  Returns list of (j, a_j, p_j, q_j) where a_j
    is the j-th partial quotient.  Denominators q_j are the classical best
    approximation denominators.
    """
    out = []
    x, y = lam, n
    # h = numerators, k = denominators;  h_{-2},h_{-1} = 0,1   k_{-2},k_{-1} = 1,0
    h_prev, h_cur = 0, 1
    k_prev, k_cur = 1, 0
    j = 0
    while y != 0 and j < 512:
        a = x // y
        x, y = y, x - a * y
        h_prev, h_cur = h_cur, a * h_cur + h_prev
        k_prev, k_cur = k_cur, a * k_cur + k_prev
        if k_cur > 0:
            out.append((j, a, k_cur, h_cur))
        j += 1
    # (j, a_j, denominator k_j, numerator h_j);  k_j is what multiplies lam
    return out


def semiconvergent_denoms(cf):
    """Intermediate fractions q = q_{j-1} + t*q_j, 1 <= t <= a_{j+1}."""
    denoms = set()
    qs = [c[2] for c in cf]
    for j in range(len(cf)):
        q_j = qs[j]
        q_jm1 = qs[j - 1] if j >= 1 else 1
        a_next = cf[j + 1][1] if j + 1 < len(cf) else 1
        t = 1
        while t <= a_next and t <= 4096:
            denoms.add(q_jm1 + t * q_j)
            t += 1
    return denoms


def centred_residue(a, lam, n):
    """|a*lam mod n| centred: min(a*lam mod n, n - a*lam mod n)."""
    r = (a * lam) % n
    return min(r, n - r)


def norm2_at(a, lam, n, s_k1, s_k2):
    return s_k1 * s_k1 * centred_residue(a, lam, n) ** 2 + s_k2 * s_k2 * a * a


def q_star(n, s_k1, s_k2):
    """Balance point: S_K1*(n/q) == S_K2*q  =>  q = sqrt(n*S_K1/S_K2)."""
    return math.sqrt(n * s_k1 / s_k2)


def nearest_convergent_to_qstar(cf, qs_target):
    """Convergent denominator closest to qs_target in log scale."""
    best, best_d = None, None
    for (j, a_j, q_j, p_j) in cf:
        if q_j <= 0:
            continue
        d = abs(math.log(q_j) - math.log(qs_target)) if qs_target > 0 else 0.0
        if best_d is None or d < best_d:
            best, best_d = (j, a_j, q_j, p_j), d
    return best


def brute_force_min(n, lam, s_k1, s_k2, cap):
    """Exhaustive min over a in [0, cap). Ground truth for small n."""
    best_a, best_n2 = 0, (n * s_k1) ** 2
    for a in range(1, cap):
        n2 = norm2_at(a, lam, n, s_k1, s_k2)
        if n2 < best_n2:
            best_a, best_n2 = a, n2
    return best_a, best_n2


# ---------------------------------------------------------------------------
# Curve generation: j=0 GLV curves (lam^2 + lam + 1 = 0 mod n)
# ---------------------------------------------------------------------------

def is_probable_prime(n_):
    if n_ < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n_ % p == 0:
            return n_ == p
    d, s = n_ - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n_)
        if x in (1, n_ - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n_
            if x == n_ - 1:
                break
        else:
            return False
    return True


def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, or None."""
    if n % 3 != 1:
        return None
    # need sqrt(-3) mod n
    t = pow(n - 3, (n - 1) // 2, n)
    if t != 1:
        return None
    s = tonelli(n - 3, n)
    if s is None:
        return None
    inv2 = pow(2, n - 2, n)
    r1 = (-1 + s) * inv2 % n
    r2 = (-1 - s) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (r1, r2)


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


def sample_glv_orders(bits, count, rng):
    """Sample primes n of the given bit size admitting a GLV eigenvalue."""
    out = []
    lo, hi = 1 << (bits - 1), (1 << bits) - 1
    tries = 0
    while len(out) < count and tries < 200000:
        tries += 1
        n = rng.randrange(lo, hi) | 1
        if not is_probable_prime(n):
            continue
        rr = glv_roots(n)
        if rr is None:
            continue
        out.append((n, min(rr)))
    return out


# ---------------------------------------------------------------------------
# secp256k1
# ---------------------------------------------------------------------------

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def corr_pearson(xs, ys):
    k = len(xs)
    mx, my = sum(xs) / k, sum(ys) / k
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float('nan')
    return sxy / math.sqrt(sxx * syy)


def corr_spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0
            for t in range(i, j + 1):
                r[order[t]] = avg
            i = j + 1
        return r
    return corr_pearson(rank(xs), rank(ys))


def pct(xs, q):
    s = sorted(xs)
    if not s:
        return float('nan')
    i = min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))
    return s[i]


# ===========================================================================
# E0 — exact vs float Lagrange-Gauss
# ===========================================================================

def e0_exact_vs_float(rng):
    print("=" * 74)
    print("E0  exact-integer vs float Lagrange-Gauss (audit of 2026-07-29 nu_hat)")
    print("=" * 74)
    print("The float reduction rounds q = <u,v>/|v|^2 in double precision.")
    print("Gram entries scale as (n*S_K1)^2 ~ 2^(4*bits), so this is checked at")
    print("each bit size actually used by the project.\n")
    print(f"  {'bits':>5} {'curves':>7} {'disagree':>9} {'max |rel diff|':>15}")
    for bits in (20, 24, 64, 128, 256):
        curves = sample_glv_orders(bits, 12, rng) if bits <= 24 else \
                 sample_glv_orders(bits, 6, rng)
        bad, maxrel = 0, 0.0
        for (n, lam) in curves:
            k1, k2 = k1k2_project(n, 0.05)
            u, v, s_k1, s_k2 = L2_basis(n, lam, k1, k2)
            _, nrm2 = gauss_reduce_exact(u, v)
            exact = math.exp(0.5 * math.log(nrm2))
            try:
                flt = gauss_reduce_float(u, v)
            except (OverflowError, ZeroDivisionError):
                flt = float('inf')
            rel = abs(flt - exact) / exact if exact > 0 and math.isfinite(flt) else 1.0
            if rel > 1e-9:
                bad += 1
            maxrel = max(maxrel, rel)
        print(f"  {bits:>5} {len(curves):>7} {bad:>9} {maxrel:>15.3e}")

    print("\n  secp256k1 (n, lam) at the eff values reported 2026-07-29:")
    print(f"    {'eff':>6} {'nu_hat exact':>14} {'nu_hat float':>14} {'rel diff':>11}")
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_project(SECP_N, eff)
        u, v, s_k1, s_k2 = L2_basis(SECP_N, SECP_LAM, k1, k2)
        _, nrm2 = gauss_reduce_exact(u, v)
        det = SECP_N * s_k1 * s_k2
        nu_e = math.exp(0.5 * (math.log(nrm2) - math.log(det)))
        try:
            nu_f = gauss_reduce_float(u, v) / math.sqrt(det)
        except (OverflowError, ZeroDivisionError):
            nu_f = float('nan')
        rel = abs(nu_f - nu_e) / nu_e if math.isfinite(nu_f) else float('nan')
        print(f"    {eff:>6.2f} {nu_e:>14.4f} {nu_f:>14.4f} {rel:>11.2e}")
    print()


# ===========================================================================
# E1 — is the minimiser a convergent denominator?
# ===========================================================================

def e1_minimiser_membership(rng):
    print("=" * 74)
    print("E1  is the exact minimiser a* a convergent denominator of lam/n ?")
    print("=" * 74)
    print("  a* read off the reduced vector as |w[1]|/S_K2.  Also brute-forced")
    print("  exhaustively at 20 bits as an independent check of the reduction.\n")
    print(f"  {'bits':>5} {'trials':>7} {'a* is conv':>11} {'a* is conv/semi':>16} "
          f"{'gauss==brute':>13}")
    for bits in (20, 24, 32, 64, 128, 256):
        curves = sample_glv_orders(bits, 40 if bits <= 32 else 20, rng)
        n_conv = n_semi = n_bf = n_bf_tot = 0
        for (n, lam) in curves:
            for eff in (0.02, 0.05, 0.15, 0.40):
                k1, k2 = k1k2_project(n, eff)
                _, _, s_k1, s_k2 = L2_basis(n, lam, k1, k2)
                _, a_star, nrm2, _ = nu_hat_exact(n, lam, k1, k2)
                cf = cf_convergents(lam, n)
                convs = {c[2] for c in cf}
                semis = convs | semiconvergent_denoms(cf)
                if a_star in convs:
                    n_conv += 1
                if a_star in semis:
                    n_semi += 1
                if bits <= 20:
                    n_bf_tot += 1
                    bf_a, bf_n2 = brute_force_min(n, lam, s_k1, s_k2,
                                                  min(n, 1 << 20))
                    if bf_n2 == nrm2:
                        n_bf += 1
        tot = len(curves) * 4
        bf_s = f"{n_bf}/{n_bf_tot}" if n_bf_tot else "n/a"
        print(f"  {bits:>5} {tot:>7} {n_conv:>7}/{tot:<3} {n_semi:>12}/{tot:<3} "
              f"{bf_s:>13}")
    print()


# ===========================================================================
# E2/E3/E4 — predictors
# ===========================================================================

def predictors(n, lam, k1, k2):
    """Return (nu_exact, nu_nearest_conv, nu_all_conv, nu_all_semi, idx_used)."""
    _, _, s_k1, s_k2 = L2_basis(n, lam, k1, k2)
    det = n * s_k1 * s_k2
    _, a_star, nrm2_exact, _ = nu_hat_exact(n, lam, k1, k2)

    cf = cf_convergents(lam, n)
    qs = q_star(n, s_k1, s_k2)

    triv = (n * s_k1) ** 2

    near = nearest_convergent_to_qstar(cf, qs)
    n2_near = min(triv, norm2_at(near[2], lam, n, s_k1, s_k2)) if near else triv

    n2_all, idx_all = triv, -1
    for (j, a_j, q_j, p_j) in cf:
        v = norm2_at(q_j, lam, n, s_k1, s_k2)
        if v < n2_all:
            n2_all, idx_all = v, j

    n2_semi = n2_all
    for q in semiconvergent_denoms(cf):
        v = norm2_at(q, lam, n, s_k1, s_k2)
        if v < n2_semi:
            n2_semi = v

    def to_nu(x):
        return math.exp(0.5 * (math.log(x) - math.log(det)))

    return (to_nu(nrm2_exact), to_nu(n2_near), to_nu(n2_all), to_nu(n2_semi),
            idx_all, a_star, qs)


def e234_predictors(rng):
    print("=" * 74)
    print("E2/E3/E4  closed-form predictors vs computed nu_hat")
    print("=" * 74)
    print("  E2 = single convergent nearest q* = sqrt(n*S_K1/S_K2)")
    print("  E3 = min over all convergent denominators")
    print("  E4 = E3 plus semiconvergent denominators")
    print("  Thread 23 falsifier: corr(E2, exact) < 0.9 kills the scale story.\n")

    for bits in (20, 24, 64, 256):
        curves = sample_glv_orders(bits, 30 if bits <= 24 else 15, rng)
        ex, p2, p3, p4 = [], [], [], []
        tight3 = tight4 = tot = 0
        for (n, lam) in curves:
            for eff in (0.02, 0.05, 0.10, 0.25, 0.50):
                k1, k2 = k1k2_project(n, eff)
                e, a2, a3, a4, _, _, _ = predictors(n, lam, k1, k2)
                ex.append(e); p2.append(a2); p3.append(a3); p4.append(a4)
                tot += 1
                if abs(a3 - e) / e < 1e-12:
                    tight3 += 1
                if abs(a4 - e) / e < 1e-12:
                    tight4 += 1
        print(f"  {bits} bits, {tot} (curve,eff) cells")
        print(f"    corr(E2 nearest-conv, exact) : pearson {corr_pearson(p2, ex):+.4f}"
              f"   spearman {corr_spearman(p2, ex):+.4f}")
        print(f"    corr(E3 all-conv,      exact) : pearson {corr_pearson(p3, ex):+.4f}"
              f"   spearman {corr_spearman(p3, ex):+.4f}")
        print(f"    corr(E4 conv+semi,     exact) : pearson {corr_pearson(p4, ex):+.4f}"
              f"   spearman {corr_spearman(p4, ex):+.4f}")
        r2 = [a / e for a, e in zip(p2, ex)]
        print(f"    E2/exact ratio  median {pct(r2,0.5):.4f}  p90 {pct(r2,0.90):.4f}"
              f"  max {max(r2):.4f}")
        print(f"    E3 exact-hit {tight3}/{tot}    E4 exact-hit {tight4}/{tot}")
        print()


# ===========================================================================
# E5 — scale dependence: the retro-explanation
# ===========================================================================

def e5_scale_dependence(rng):
    print("=" * 74)
    print("E5  scale dependence: one FIXED (n, lam), sweep the scale")
    print("=" * 74)
    print("  If the winning convergent index changes with (K1,K2) then no")
    print("  scale-free function of the CF of lam/n can predict nu_hat.  This is")
    print("  the retro-explanation for the six invariants falsified in June")
    print("  (delta/n, kappa(M), q_cf, max_q_cf, max_a, a_corn/n) and for mu.\n")

    curves = sample_glv_orders(24, 3, rng)
    for (n, lam) in curves:
        cf = cf_convergents(lam, n)
        print(f"  n={n}  lam={lam}   CF length {len(cf)}  "
              f"max partial quotient {max(c[1] for c in cf)}")
        print(f"    {'eff':>6} {'q*':>10} {'a*':>10} {'conv idx':>9} "
              f"{'nu_hat':>8}")
        seen = set()
        for eff in (0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75):
            k1, k2 = k1k2_project(n, eff)
            e, _, _, _, idx, a_star, qs = predictors(n, lam, k1, k2)
            seen.add(idx)
            print(f"    {eff:>6.2f} {qs:>10.1f} {a_star:>10} {idx:>9} {e:>8.4f}")
        print(f"    -> distinct winning convergent indices across eff: "
              f"{sorted(seen)}   (scale-free invariant cannot see this)\n")


# ===========================================================================
# E6 — secp256k1 placement, with the realising convergent
# ===========================================================================

def e6_secp256k1():
    print("=" * 74)
    print("E6  secp256k1: which convergent of lam/n realises nu_hat")
    print("=" * 74)
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0, "lam check failed"
    print("  lam^2 + lam + 1 = 0 mod n   VERIFIED")
    cf = cf_convergents(SECP_LAM, SECP_N)
    print(f"  CF(lam/n): length {len(cf)}, max partial quotient "
          f"{max(c[1] for c in cf)}\n")
    print(f"  {'eff':>6} {'q*':>12} {'a*':>26} {'idx':>5} {'nu_hat':>8} "
          f"{'E2':>8} {'E3':>8}")
    for eff in (0.05, 0.0993, 0.10, 0.25, 0.50):
        k1, k2 = k1k2_project(SECP_N, eff)
        e, a2, a3, a4, idx, a_star, qs = predictors(SECP_N, SECP_LAM, k1, k2)
        print(f"  {eff:>6.4f} {qs:>12.4e} {a_star:>26} {idx:>5} {e:>8.4f} "
              f"{a2:>8.4f} {a3:>8.4f}")
    print("\n  Reminder: this is a 20/24-bit -> 256-bit extrapolation of an")
    print("  empirical regularity, and the whole thread is conditional on a")
    print("  NON-STANDARD nonce generator k = k1 + lam*k2 with k1 bounded.")
    print("  It is NOT an attack on secp256k1.\n")


def e7_selector_and_hermite(rng):
    """
    Sharpen the E2 negative result: how often does 'nearest convergent to q*'
    pick the RIGHT convergent index?  And validate every computed nu_hat
    against the dimension-2 Hermite bound nu_hat <= (2/sqrt3)^(1/2) = 1.07457,
    which any true lattice minimum must satisfy.
    """
    HERMITE = math.sqrt(2.0 / math.sqrt(3.0))
    print("=" * 74)
    print("E7  selector diagnostics + Hermite validation")
    print("=" * 74)
    print(f"  dim-2 Hermite bound: nu_hat <= (2/sqrt3)^(1/2) = {HERMITE:.5f}")
    print("  Any predictor exceeding it is provably not a lattice minimum.\n")
    print(f"  {'bits':>5} {'cells':>6} {'q*-rule idx hit':>16} {'|idx-idx*| med':>15} "
          f"{'exact>H':>8} {'E2>H':>6}")
    for bits in (20, 24, 64, 256):
        curves = sample_glv_orders(bits, 30 if bits <= 24 else 15, rng)
        hit = tot = viol_e = viol_2 = 0
        gaps = []
        for (n, lam) in curves:
            for eff in (0.02, 0.05, 0.10, 0.25, 0.50):
                k1, k2 = k1k2_project(n, eff)
                _, _, s_k1, s_k2 = L2_basis(n, lam, k1, k2)
                e, a2, a3, a4, idx, a_star, qs = predictors(n, lam, k1, k2)
                cf = cf_convergents(lam, n)
                near = nearest_convergent_to_qstar(cf, qs)
                tot += 1
                if near and near[0] == idx:
                    hit += 1
                if near:
                    gaps.append(abs(near[0] - idx))
                if e > HERMITE + 1e-9:
                    viol_e += 1
                if a2 > HERMITE + 1e-9:
                    viol_2 += 1
        med = sorted(gaps)[len(gaps) // 2] if gaps else float('nan')
        print(f"  {bits:>5} {tot:>6} {hit:>10}/{tot:<4} {med:>15} "
              f"{viol_e:>8} {viol_2:>6}")
    print("\n  exact>H must be 0 (validates the reduction);  E2>H > 0 is an")
    print("  independent proof that the nearest-q* rule is not the minimum.\n")


def main():
    seed = 20260729
    if len(sys.argv) > 1:
        seed = int(sys.argv[1])
    rng = random.Random(seed)
    print("GLV-HNP Phase 2 — Thread 23: closed form for lambda_1(L2)")
    print(f"seed = {seed}\n")
    e0_exact_vs_float(rng)
    e1_minimiser_membership(rng)
    e234_predictors(rng)
    e5_scale_dependence(rng)
    e6_secp256k1()
    e7_selector_and_hermite(rng)


if __name__ == "__main__":
    main()
