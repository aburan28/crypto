"""
GLV-HNP Phase 2, Thread 23: is the "K1 wall" a lattice-geometry ceiling?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, T4/T4b/T5):
  T4  found a hard "K1 wall": the 12-bit/2677 curve recovers d for K1 <= 4
      but never for K1 >= 8.  T4b showed more data does not rescue it
      (m = 8/12/16/24/32 -> 0,0,1,0,1 of 5).
  T5  found that the planted vector is NEVER lambda_1 of the Phase-2 lattice:
      the trivial vector n*S_D*e_m (norm n) is always shorter.
  The 2026-07-29 run proposed Thread 23: "reformulate the lattice so the
  target is lambda_1 -- project out the trivial n*e_m direction.  Falsifier:
  if sv/pv rises above 1 and the K1 wall moves outward, the reformulation is
  a real improvement; if the wall stays at K1 ~ 4-6, the wall is
  information-theoretic and Phase 2 is at its ceiling."

This script tests a sharper hypothesis that subsumes both.

  H23:  Recovery is governed by the Gaussian-heuristic ratio

            R  =  ||v_planted|| / GH(L),      GH(L) = sqrt(dim/(2*pi*e)) * det(L)^(1/dim)

        and NOT by any curve-level invariant, not by lambda*, and not by the
        presence of the trivial vector.

  Closed form.  With the standard balanced scaling
  S_K1 = n/K1, S_D = 1, S_K2 = n/K2, S_KANNAN = n, and eff = K1*K2/n:

        dim      = 2m + 2
        det(L)   = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN
                 = n^(3m+1) / (K1*K2)^m  =  n^(2m+1) / eff^m
        E||v||^2 = n^2 * (2m + 4) / 3            (2m coords ~ U[0,n), d ~ U[0,n), Kannan = n)

        R(m, n, eff) = sqrt((2m+4)/(2m+2)) * sqrt(2*pi*e/3) * n^(1/(2m+2)) * eff^(m/(2m+2))

  Two predictions fall straight out, and they are the point of this script:

   (P1) ASYMPTOTIC eff CEILING.  As m -> infinity,
            R -> sqrt(2*pi*e/3) * sqrt(eff) = 2.3860 * sqrt(eff),
        so no amount of data can make R < 1 once
            eff > 1 / (2*pi*e/3) = 0.17565.
        This is a *lattice-formulation* ceiling, strictly stronger than the
        information-theoretic condition eff < 1 (m > log n / log(1/eff)).

   (P2) FINITE-m TRADEOFF.  Below the ceiling, R < 1 needs
            m > (ln n + ln(2*pi*e/3) * 2) / (ln(1/eff) - 2*ln(2*pi*e/3)/2)
        i.e. the m at which the n^(1/(2m+2)) penalty is paid off.  For the
        12-bit/2677 curve at K1=8 (eff = 0.1572, just *under* the ceiling)
        this predicts m ~ 87 -- far above the m <= 32 that T4b tested.
        So T4b's "more data does not rescue it" was measured on the wrong
        side of the predicted transition.

Experiments:
  E1  Validate R as a predictor.  Fresh 17-bit j=0 GLV curves x 4 eff levels,
      m = 12.  Report AUC of R (and of lambda*, and of eff) against success.
  E2  DECISIVE.  12-bit/2677 at K1=8 (eff = 0.157, R->0.946 asymptotically),
      m = 12..160.  If H23 is right, recovery switches ON somewhere near the
      m where R crosses 1.  If recovery never turns on, H23 is falsified and
      the wall is real.
  E2b CONTROL.  Same curve at K1=16 (eff = 0.314, ABOVE the ceiling,
      R -> 1.337).  H23 predicts recovery NEVER turns on, at any m.
  E3  Thread 23's own proposal: drop the d column entirely (d is only defined
      mod n, so the d coordinate carries no information and is exactly what
      puts the trivial vector n*S_D*e_m into the lattice).  This gives a
      (2m+1)-dimensional lattice with det = n^(2m)/eff^m and no trivial
      vector.  Measure (a) whether the trivial vector is really gone,
      (b) the change in R, (c) whether the K1 wall moves.

RESULTS (2026-08-05 run; full output in glv_hnp_thread23_gh_ceiling_output.txt)

  E1  R tracks success (AUC 0.918) but is CONFOUNDED: at fixed m and fixed
      n-size, R is a monotone function of eff, and AUC(eff) = 0.9184 to four
      decimals -- exactly the same number.  R buys nothing over eff here.
      AUC(lam*) = 0.5114, i.e. chance; this independently reconfirms T3.

  E2  P2 FALSIFIED.  12-bit/2677 at K1=8 (eff = 0.1572, under the ceiling):
      R crosses 1 at m ~ 104 as predicted, and recovery does NOT switch on.
          m      12    24    40    56    72    88   104   128   160
          dim    26    50    82   114   146   178   210   258   322
          R    1.428 1.172 1.078 1.039 1.018 1.005 0.996 0.986 0.978
          win   0/3   0/3   1/3   0/3   1/3   1/3   0/3   0/3   0/3
      So R < 1 is necessary-ish but nowhere near sufficient.  More data does
      not rescue the instance, confirming T4b out to dim 322 (T4b stopped at
      m = 32).

  E2b CONTROL passes: K1=16 (eff = 0.314, above the ceiling) is 0/3 at every m.

  E3  THREAD 23's PROPOSAL IS DEAD.  Deleting the d column removes the trivial
      vector (sv/pv rises from ~0.45 to 1.0000 on the instances that succeed --
      the planted vector really does become lambda_1) and changes R by at most
      3.3%.  It changes the outcome on ZERO of the 16 configurations tested:
          12-bit/2677 m=12, K1 = 2,3,4,6,8,12,16
            orig 5/5 5/5 5/5 2/5 0/5 0/5 0/5
            nod  5/5 5/5 5/5 2/5 0/5 0/5 0/5
          12-bit/2557 m=12, K1 = 4,8,12,16,24
            orig 5/5 5/5 4/5 1/5 0/5   nod  5/5 5/5 4/5 1/5 0/5
          12-bit/2677 K1=8, m = 12,40,88,160
            orig 0/3 1/3 1/3 0/3       nod  0/3 1/3 1/3 0/3
      Bit-identical, seed for seed.  The trivial vector n*S_D*e_m is not an
      obstruction: recover_d() reads row[m] mod n, and the trivial vector
      generates exactly the kernel of that map, so the d column is a redundant
      coordinate and the two lattices are the same attack.

  => T5's framing ("the planted vector is never lambda_1, so recovery is a BDD
     condition") is geometrically true but operationally irrelevant.  See
     glv_hnp_thread23_mechanism.py for what the recovery event actually is.

Run: python3 glv_hnp_thread23_gh_ceiling.py
"""

import math
import random
import sys
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

TWO_PI_E = 2.0 * math.pi * math.e
GH_CONST = math.sqrt(TWO_PI_E / 3.0)          # 2.38601...
EFF_CEILING = 1.0 / (TWO_PI_E / 3.0)          # 0.175650...

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def identify_b_for_n(p, n):
    for bb in range(1, min(p, 60)):
        rhs_ok = None
        rng = random.Random(7)
        for _ in range(200):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + bb) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                rhs_ok = (x, y)
                break
        if rhs_ok is None:
            continue
        if ec_mul(rhs_ok, n, p) is None:
            return bb
    return None

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim construction; see glv_hnp_phase2_20bit.py:262)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 400000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Original Phase-2 lattice, dim = 2m+2, with the d column."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def build_glv_lattice_nod(sigs, n, lam, k1_bound, k2_bound):
    """Thread-23 reformulation: the d column is deleted.

    Columns: 0..m-1 congruence (scale S_K1), m..2m-1 the k2 unknowns
    (scale S_K2), 2m the Kannan column.  dim = 2m+1.
    Generators (2m+2 of them, rank 2m+1 -- LLL will emit one zero row):
        R_i  = n*S_K1 * e_i                         (i < m)
        Rd   = (B_i*S_K1)_i                         (coefficient = d)
        Rk_i = -lam*S_K1 * e_i + S_K2 * e_{m+i}     (coefficient = k2_i)
        RK   = (A_i*S_K1)_i + S_KANNAN * e_{2m}
    Planted vector = (k1_i*S_K1, k2_i*S_K2, S_KANNAN); d never appears, so
    the trivial vector n*S_D*e_m of the original lattice does not exist here.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KANNAN
    rows.append(r)
    return rows

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def planted_vector_nod(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Gaussian-heuristic ratio
# ---------------------------------------------------------------------------

def log_det(n, m, k1_bound, k2_bound, drop_d):
    """ln det(L).  Original: (n*S_K1)^m * S_D * S_K2^m * S_KANNAN.
    Dropping the d column removes the S_D pivot AND divides the covolume by
    n (the Rd generator has order n in (Z/n)^m), so ln det loses ln n."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ld = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2) + math.log(S_KAN)
    if drop_d:
        ld -= math.log(n)
    return ld

def gh_norm(n, m, k1_bound, k2_bound, drop_d=False):
    dim = (2 * m + 1) if drop_d else (2 * m + 2)
    ld = log_det(n, m, k1_bound, k2_bound, drop_d)
    return math.sqrt(dim / TWO_PI_E) * math.exp(ld / dim)

def planted_norm_expected(n, m, k1_bound, k2_bound, drop_d=False):
    """E||v||, k1~U[0,K1), k2~U[0,K2), d~U[0,n); E[x^2] = X^2/3."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    tot = m * (k1_bound * S_K1) ** 2 / 3.0 + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2
    if not drop_d:
        tot += (n * S_D) ** 2 / 3.0
    return math.sqrt(tot)

def R_ratio(n, m, k1_bound, k2_bound, drop_d=False):
    return (planted_norm_expected(n, m, k1_bound, k2_bound, drop_d)
            / gh_norm(n, m, k1_bound, k2_bound, drop_d))

def R_closed_form(n, m, eff):
    """The asymptotic-form predictor stated in the docstring (sanity check)."""
    return (math.sqrt((2 * m + 4) / (2 * m + 2)) * GH_CONST
            * n ** (1.0 / (2 * m + 2)) * eff ** (m / (2.0 * m + 2)))

def m_star(n, eff):
    """Smallest m with R_closed_form < 1, or None if eff is above the ceiling."""
    if eff >= EFF_CEILING:
        return None
    for m in range(2, 100000):
        if R_closed_form(n, m, eff) < 1.0:
            return m
    return None

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_d(reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return True
    return False

def recover_d_nod(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """No d column: read (k1_0, k2_0) off the row and solve for d."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        for i in range(m):
            x, y = sign * row[i], sign * row[m + i]
            if x % S_K1 or y % S_K2:
                continue
            k1, k2 = x // S_K1, y // S_K2
            if not (0 <= k1 < k1_bound and 0 <= k2 < k2_bound):
                continue
            # k1 = A_i + B_i*d - lam*k2 (mod n)
            d_cand = (k1 - sigs[i]['A'] + lam * k2) % n * modinv(sigs[i]['B'], n) % n
            if d_cand == d_secret:
                return True
            break
    return False

def trivial_vector_share(reduced, m, n):
    """|sv[m]| / n for the shortest row -- 1.0000 means it IS n*S_D*e_m."""
    best, bn = None, None
    for row in reduced:
        nn = norm(row)
        if nn == 0:
            continue
        if bn is None or nn < bn:
            bn, best = nn, row
    if best is None:
        return None, None
    return bn, abs(best[m]) / float(n)

# ---------------------------------------------------------------------------
# One experiment
# ---------------------------------------------------------------------------

def run_once(curve, m, d_secret, k1_bound, seed, drop_d=False,
             use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if drop_d:
        rows = build_glv_lattice_nod(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1
        A = IntegerMatrix.from_matrix(rows)
    else:
        rows = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 2
        A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    if drop_d:
        ok = recover_d_nod(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret)
        pv = norm(planted_vector_nod(sigs, n, k1_bound, k2_bound))
        sv = min((norm(r) for r in reduced if norm(r) > 0), default=0.0)
        share = None
    else:
        ok = recover_d(reduced, m, n, S_KAN, d_secret)
        pv = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
        sv, share = trivial_vector_share(reduced, m, n)
    return {'ok': ok, 'pv': pv, 'sv': sv, 'share': share}

def sweep(curve, m, k1_bound, seeds, drop_d=False, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    wins, svpv, shares = 0, [], []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = run_once(curve, m, d_trial, k1_bound, seed, drop_d=drop_d,
                     use_bkz=use_bkz, bkz_beta=bkz_beta)
        if r is None:
            continue
        wins += r['ok']
        if r['pv'] > 0:
            svpv.append(r['sv'] / r['pv'])
        if r['share'] is not None:
            shares.append(r['share'])
    return (wins, len(seeds),
            sum(svpv) / len(svpv) if svpv else float('nan'),
            sum(shares) / len(shares) if shares else float('nan'))

def auc(scores, labels):
    """Area under ROC.  scores: higher = predicted FAILURE."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    tot = 0.0
    for a in pos:
        for bb in neg:
            tot += 1.0 if a < bb else (0.5 if a == bb else 0.0)
    return tot / (len(pos) * len(neg))

# ---------------------------------------------------------------------------

def banner(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)
    sys.stdout.flush()

print("=" * 78)
print("Thread 23 -- Gaussian-heuristic ceiling of the GLV-HNP Phase-2 lattice")
print("=" * 78)
print(f"  sqrt(2*pi*e/3) = {GH_CONST:.5f}")
print(f"  asymptotic eff ceiling 1/(2*pi*e/3) = {EFF_CEILING:.6f}")

# --- historical curves -----------------------------------------------------
HIST = [
    # label,             p,    b, n,    lam
    ("8-bit/199",        211,  2, 199,  106),
    ("12-bit/2557",      2557, 2, 2659, 1755),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185),
]
hist = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G)))

banner("Sanity: closed form vs exact numeric R (must agree to ~1e-3)")
print(f"{'curve':<20}{'K1':>4}{'m':>5}{'eff':>9}{'R_num':>9}{'R_cf':>9}{'m*':>7}")
for label, cv in hist:
    p, b, n, lam, G = cv
    k2 = math.isqrt(n) + 1
    for k1 in (2, 8):
        for m in (8, 12):
            eff = k1 * k2 / n
            ms = m_star(n, eff)
            print(f"{label:<20}{k1:>4}{m:>5}{eff:>9.4f}"
                  f"{R_ratio(n, m, k1, k2):>9.4f}{R_closed_form(n, m, eff):>9.4f}"
                  f"{(str(ms) if ms else 'never'):>7}")

# ===========================================================================
# E1 -- validate R as a predictor on fresh 17-bit curves
# ===========================================================================
banner("E1: is R a predictor?  Fresh 17-bit j=0 GLV curves, m=12, 5 seeds")

def find_curves(bits, want, seed=1):
    out = []
    p = sympy.nextprime(2 ** (bits - 1))
    while p < 2 ** bits and len(out) < want:
        if p % 3 == 1:
            dec = eisenstein_decompose(p)
            if dec:
                a, bb = dec
                for t in j0_traces(a, bb):
                    n = p + 1 - t
                    if n > 4 and sympy.isprime(n) and n % 3 == 1:
                        lam = glv_eigenvalue(n)
                        if lam is None:
                            continue
                        bcurve = identify_b_for_n(p, n)
                        if bcurve is None:
                            continue
                        G = find_generator(p, bcurve, n)
                        if G is None:
                            continue
                        out.append((p, bcurve, n, lam, G))
                        break
        p = sympy.nextprime(p + 1)
    return out

CURVES17 = find_curves(17, 12)
print(f"found {len(CURVES17)} 17-bit curves")
SEEDS5 = [42, 1234, 9999, 555, 31337]
M_E1 = 12

rows_e1 = []
print(f"\n{'p':>8}{'n':>8}{'lam*':>8}{'eff':>8}{'R':>8}{'m*':>7}{'win':>7}")
for cv in CURVES17:
    p, b, n, lam, G = cv
    k2 = math.isqrt(n) + 1
    lam_star = min(lam, n - lam) / n
    for target_eff in (0.05, 0.12, 0.18, 0.28):
        k1 = max(2, int(round(target_eff * n / k2)))
        eff = k1 * k2 / n
        R = R_ratio(n, M_E1, k1, k2)
        w, t, svpv, share = sweep(cv, M_E1, k1, SEEDS5)
        ms = m_star(n, eff)
        rows_e1.append({'p': p, 'n': n, 'lam_star': lam_star, 'eff': eff,
                        'R': R, 'win': w, 'tot': t})
        print(f"{p:>8}{n:>8}{lam_star:>8.4f}{eff:>8.4f}{R:>8.4f}"
              f"{(str(ms) if ms else 'never'):>7}{f'{w}/{t}':>7}")
        sys.stdout.flush()

lab = [r['win'] >= 3 for r in rows_e1]
print(f"\nn = {len(rows_e1)} configurations, {sum(lab)} majority-successes")
print(f"  AUC(R,      higher=failure) = {auc([r['R'] for r in rows_e1], lab):.4f}")
print(f"  AUC(eff,    higher=failure) = {auc([r['eff'] for r in rows_e1], lab):.4f}")
print(f"  AUC(-lam*,  higher=failure) = {auc([-r['lam_star'] for r in rows_e1], lab):.4f}")
suc_R = [r['R'] for r in rows_e1 if r['win'] >= 3]
fail_R = [r['R'] for r in rows_e1 if r['win'] < 3]
if suc_R and fail_R:
    print(f"  R range | success : [{min(suc_R):.4f}, {max(suc_R):.4f}]")
    print(f"  R range | failure : [{min(fail_R):.4f}, {max(fail_R):.4f}]")

# ===========================================================================
# E2 -- decisive: does large m rescue the K1=8 wall on 12-bit/2677?
# ===========================================================================
banner("E2: 12-bit/2677, K1=8 (eff just UNDER the ceiling) -- sweep m to 160")

cv2677 = hist[2][1]
p, b, n2677, lam2677, G = cv2677
k2_2677 = math.isqrt(n2677) + 1
for K1 in (8, 16):
    eff = K1 * k2_2677 / n2677
    ms = m_star(n2677, eff)
    print(f"  K1={K1:<3} eff={eff:.4f}  asymptotic R={GH_CONST*math.sqrt(eff):.4f}"
          f"  predicted m* = {ms if ms else 'never (above ceiling)'}")

M_GRID = [12, 24, 40, 56, 72, 88, 104, 128, 160]
SEEDS3 = [42, 1234, 9999]
print(f"\n{'m':>5}{'dim':>6}{'R':>8}{'win':>7}{'sv/pv':>9}{'|sv[m]|/n':>11}")
e2 = []
for m in M_GRID:
    R = R_ratio(n2677, m, 8, k2_2677)
    w, t, svpv, share = sweep(cv2677, m, 8, SEEDS3)
    e2.append((m, R, w, t))
    print(f"{m:>5}{2*m+2:>6}{R:>8.4f}{f'{w}/{t}':>7}{svpv:>9.4f}{share:>11.4f}")
    sys.stdout.flush()

banner("E2b CONTROL: same curve, K1=16 (eff ABOVE the ceiling) -- must never work")
print(f"\n{'m':>5}{'dim':>6}{'R':>8}{'win':>7}{'sv/pv':>9}")
for m in [12, 40, 88, 160]:
    R = R_ratio(n2677, m, 16, k2_2677)
    w, t, svpv, share = sweep(cv2677, m, 16, SEEDS3)
    print(f"{m:>5}{2*m+2:>6}{R:>8.4f}{f'{w}/{t}':>7}{svpv:>9.4f}")
    sys.stdout.flush()

# ===========================================================================
# E3 -- Thread 23's own proposal: delete the d column
# ===========================================================================
banner("E3: Thread-23 reformulation (drop the d column, kill the trivial vector)")
print("Predicted gain in R from dropping the d column:")
print(f"{'curve':<20}{'K1':>4}{'m':>5}{'R(2m+2)':>10}{'R(2m+1)':>10}{'gain':>8}")
for label, cv in hist:
    p, b, n, lam, G = cv
    k2 = math.isqrt(n) + 1
    for k1 in (4, 8):
        for m in (12, 40):
            Ra = R_ratio(n, m, k1, k2, drop_d=False)
            Rb = R_ratio(n, m, k1, k2, drop_d=True)
            print(f"{label:<20}{k1:>4}{m:>5}{Ra:>10.4f}{Rb:>10.4f}{Rb/Ra:>8.4f}")

print("\nEmpirical K1 wall, original lattice vs d-dropped lattice (12-bit/2677, m=12):")
print(f"{'K1':>4}{'eff':>8}{'orig':>8}{'nod':>8}{'R_orig':>9}{'R_nod':>8}")
for k1 in (2, 3, 4, 6, 8, 12, 16):
    eff = k1 * k2_2677 / n2677
    wo, to, _, _ = sweep(cv2677, 12, k1, SEEDS5, drop_d=False)
    wn, tn, _, _ = sweep(cv2677, 12, k1, SEEDS5, drop_d=True)
    print(f"{k1:>4}{eff:>8.4f}{f'{wo}/{to}':>8}{f'{wn}/{tn}':>8}"
          f"{R_ratio(n2677,12,k1,k2_2677):>9.4f}{R_ratio(n2677,12,k1,k2_2677,True):>8.4f}")
    sys.stdout.flush()

print("\nSame for 12-bit/2557 (lam* = 0.34), m=12:")
cv2557 = hist[1][1]
n2557 = cv2557[2]
k2_2557 = math.isqrt(n2557) + 1
print(f"{'K1':>4}{'eff':>8}{'orig':>8}{'nod':>8}")
for k1 in (4, 8, 12, 16, 24):
    wo, to, _, _ = sweep(cv2557, 12, k1, SEEDS5, drop_d=False)
    wn, tn, _, _ = sweep(cv2557, 12, k1, SEEDS5, drop_d=True)
    print(f"{k1:>4}{k1*k2_2557/n2557:>8.4f}{f'{wo}/{to}':>8}{f'{wn}/{tn}':>8}")
    sys.stdout.flush()

print("\nDoes the d-dropped lattice help at large m (K1=8, 12-bit/2677)?")
print(f"{'m':>5}{'R_nod':>8}{'nod win':>9}")
for m in [12, 40, 88, 160]:
    w, t, _, _ = sweep(cv2677, m, 8, SEEDS3, drop_d=True)
    print(f"{m:>5}{R_ratio(n2677,m,8,k2_2677,True):>8.4f}{f'{w}/{t}':>9}")
    sys.stdout.flush()

banner("DONE")
