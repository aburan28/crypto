"""
GLV-HNP Phase 2, Thread 23: make the planted vector lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 lattice L (dim 2m+2) always contains the *trivial* vector
      t = n * S_D * e_m        (= n*row_m - sum_i B_i*row_i)
  of norm exactly n*S_D, while the planted vector has norm
      ||v_planted|| ~ n * sqrt(2m/3 + 4/3).
  So t is shorter by a factor ~0.8*sqrt(m) on every instance, success and
  failure alike, and the planted vector is NEVER lambda_1.  t carries no
  information (d is only defined mod n) and no choice of S_D removes it,
  because both vectors scale linearly in S_D.  Recovery in L is therefore a
  BDD/coset condition, not an SVP condition.

  Proposed fix (2026-07-29 "next step"): quotient out the trivial direction.
  Since t is a multiple of the coordinate vector e_m, the quotient
  L / (Z*t) is realised concretely by DELETING COLUMN m (the d-column):
      pi : Z^{2m+2} -> Z^{2m+1},  forget coordinate m.
  Then ker(pi|_L) = Z*t exactly, so pi(L) has rank 2m+1 and
      det pi(L) = det L / (n*S_D) = (n*S_K1)^m * S_K2^m * S_KANNAN / n.
  The image of the planted vector,
      v' = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN),
  survives with norm ~ n*sqrt(2m/3 + 1).  The secret d is no longer a lattice
  coordinate, but it is recovered *afterwards* in closed form from any single
  signature:  d = (k1_i + lam*k2_i - A_i) * B_i^{-1} mod n.

FALSIFIER (stated by the 2026-07-29 run, verbatim):
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  U1  structural verification of pi: rank, determinant formula, trivial vector
      killed, sv/pv in L vs pi(L), uSVP gap gamma = GH(pi(L))/||v'||.
  U2  head-to-head replication of the 2026-07-29 T4 grid (the two historical
      12-bit curves, K1 in {2,3,4,6,8,12,16,24}, 5 seeds): L vs pi(L) vs
      an explicit Babai-CVP arm that drops the Kannan embedding as well.
  U3  T4b replication: m-sweep at K1=8 on the lam*=0.07 curve, all three arms.
  U4  is gamma the wall?  eff = K1*K2/n sweep on fresh 17-bit curves, testing
      whether gamma ~ 1 predicts the success boundary.
  U5  BKZ(20)/BKZ(30) on pi(L) at the wall - is the residual failure a
      reduction-strength problem or an information-theoretic one?

Arms:
  ORIG  L (dim 2m+2), LLL, read d off the d-column   [the historical attack]
  PROJ  pi(L) (dim 2m+1), LLL, read (k1,k2) off a row with |last| = S_KANNAN,
        then solve for d
  CVP   L'' (dim 2m; no d-column, no Kannan row), Babai nearest-plane against
        the target a = (A_i*S_K1 | 0), then solve for d

RESULTS (2026-08-02 run):

  U1  The projection is correct.  rank pi(L) = 2m+1 on all three historical
      curves; det pi(L) = det L/(n*S_D) verified EXACTLY (sympy integer
      determinant vs the closed form); pi(v_planted) verified to be an integer
      combination of the pi(L) basis.  sv/pv rises 0.6035 -> 0.8428,
      0.5170 -> 0.5324, 0.4221 -> 0.8127.  In the success regime of the
      lam*=0.07 curve it reaches exactly 1.0000 at K1 = 2, 3, 4 - i.e. the
      planted vector IS lambda_1 of pi(L) there, which is what Thread 23 set
      out to achieve.

  U2  THE REFORMULATION CHANGES NOTHING.  PROJ and ORIG agree on all 16 cells
      of the T4 grid, exactly, on both curves.  The K1 wall does not move:
      2557 stays at K1 in (8,12], 2677 stays at K1 in (4,6].  Under the
      falsifier stated on 2026-07-29 this settles it: the wall is
      information-theoretic, not an artifact of the trivial vector.  The CVP
      arm is strictly WORSE (5/5 -> 4/5 at K1=3 and 5/5 -> 2/5 at K1=4 on the
      2677 curve); Babai nearest-plane loses to LLL-on-Kannan.

  U3  T4b replicates exactly: m = 8/12/16/24/32 -> 0,0,1,0,1 for both ORIG and
      PROJ (0,0,0,0,0 for CVP), reproducing the 2026-07-29 numbers on an
      independently written implementation.

  U4  gamma = GH(pi(L))/||v'|| is NOT a universal threshold.  Aggregate over
      6 fresh 17-bit curves: eff=0.05 (gamma 1.181) -> 30/30; 0.10 (0.849) ->
      12/30; 0.15 (0.703) -> 11/30; 0.20 (0.612) -> 9/30; 0.25 (0.549) -> 6/30.
      But the per-curve split is bimodal, and 2557 succeeds at gamma=0.671
      while 2677 fails at gamma=0.841 - so gamma alone does not locate the
      wall.  The residual lam-dependence is chased in
      glv_hnp_phase2_lamstar_third.py (U6), which falsifies |lam* - 1/3| and
      finds that the lambda-block quantity rho separates with AUC 0.89-0.98.

  U5  BKZ does not help.  BKZ(20) and BKZ(30) on pi(L) match LLL on 5 of 6
      cells; the single difference is 0/5 -> 1/5 on 2677 at K1=6.  Consistent
      with U2's information-theoretic reading.

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
from fractions import Fraction

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + CM helpers
# (verbatim from glv_hnp_phase2_lambda_threshold.py so the comparison is exact)
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
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Scales, signatures, and the three lattices
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- verbatim from the 2026-06-15 construction."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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

def build_orig(sigs, n, lam, k1_bound, k2_bound):
    """ORIG arm: the historical Phase-2 lattice L, dim 2m+2.
    columns: [0,m)=k1 block | m = d | [m+1,2m+1)=k2 block | 2m+1 = Kannan."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def build_proj(sigs, n, lam, k1_bound, k2_bound):
    """PROJ arm: pi(L) = L with column m (the d-column) deleted.
    2m+2 generators spanning a rank-(2m+1) lattice; the single relation is
    n*row_m - sum_i B_i*row_i = pi(t) = 0.
    columns: [0,m)=k1 block | [m,2m)=k2 block | 2m = Kannan."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(2 * m + 2)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def build_cvp(sigs, n, lam, k1_bound, k2_bound):
    """CVP arm: L'' = pi(L) with the Kannan row AND the Kannan column dropped.
    Returns (basis rows spanning rank 2m in dim 2m, target a).
    The planted vector is  a + u  for the closest u in L'' to -a."""
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    a = [sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    return rows, a

def planted_orig(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def planted_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def d_from_k(k1, k2, sigs, n, lam, d_secret=None):
    """Closed-form d from a candidate (k1_i, k2_i) vector.  Verified against
    every signature, so false positives are ruled out independently of
    knowing d_secret."""
    m = len(sigs)
    for i in range(m):
        B = sigs[i]['B']
        if B % n == 0:
            continue
        d = (k1[i] + lam * k2[i] - sigs[i]['A']) * modinv(B, n) % n
        if d == 0:
            continue
        if all((sigs[j]['A'] + sigs[j]['B'] * d - k1[j] - lam * k2[j]) % n == 0
               for j in range(m)):
            return d
    return None

def recover_orig(reduced, m, n, S_KAN, d_secret):
    """Historical criterion: a row with |last| = S_KANNAN, d read off column m."""
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sgn = 1 if last > 0 else -1
        d_cand = (sgn * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def recover_proj(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Projected criterion: a row with |last| = S_KANNAN; read (k1,k2) off the
    two blocks and solve for d."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sgn = 1 if last > 0 else -1
        r = [sgn * x for x in row]
        if any(r[i] % S_K1 for i in range(m)): continue
        if any(r[m + i] % S_K2 for i in range(m)): continue
        k1 = [r[i] // S_K1 for i in range(m)]
        k2 = [r[m + i] // S_K2 for i in range(m)]
        d = d_from_k(k1, k2, sigs, n, lam, d_secret)
        if d is not None and d == d_secret:
            return d
    return None

# ---------------------------------------------------------------------------
# Babai nearest-plane (exact rational Gram-Schmidt; dims here are <= 26)
# ---------------------------------------------------------------------------

def gram_schmidt_exact(basis):
    bstar, mu = [], []
    for v in basis:
        w = [Fraction(x) for x in v]
        coeffs = []
        for j, bs in enumerate(bstar):
            den = sum(y * y for y in bs)
            if den == 0:
                coeffs.append(Fraction(0))
                continue
            c = sum(Fraction(v[k]) * bs[k] for k in range(len(v))) / den
            coeffs.append(c)
            w = [w[k] - c * bs[k] for k in range(len(w))]
        bstar.append(w)
        mu.append(coeffs)
    return bstar

def babai_nearest_plane(basis, target):
    """Closest-vector approximation to `target` in the lattice spanned by the
    (already LLL-reduced) rows of `basis`."""
    bstar = gram_schmidt_exact(basis)
    b = [Fraction(x) for x in target]
    coeffs = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        den = sum(y * y for y in bstar[i])
        if den == 0:
            continue
        c = sum(b[k] * bstar[i][k] for k in range(len(b))) / den
        ci = int((c.numerator * 2 + c.denominator) // (c.denominator * 2)) \
            if c >= 0 else -int((-c.numerator * 2 + c.denominator) // (c.denominator * 2))
        coeffs[i] = ci
        b = [b[k] - ci * Fraction(basis[i][k]) for k in range(len(b))]
    out = [0] * len(target)
    for i, ci in enumerate(coeffs):
        if ci:
            for k in range(len(out)):
                out[k] += ci * basis[i][k]
    return out

# ---------------------------------------------------------------------------
# Gaussian heuristic
# ---------------------------------------------------------------------------

def gh_from_logdet(dim, logdet):
    """Gaussian heuristic  GH = det^(1/dim) * Gamma(1+dim/2)^(1/dim) / sqrt(pi)."""
    return math.exp(logdet / dim + math.lgamma(1 + dim / 2) / dim
                    - 0.5 * math.log(math.pi))

def proj_logdet(m, n, k1_bound, k2_bound):
    """log det pi(L) = m*log(n*S_K1) + m*log(S_K2) + log(S_KANNAN) - log(n*S_D)."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_KAN)
            - math.log(n * S_D))

def planted_norm_proj(m, n, k1_bound, k2_bound):
    """E||v'|| with k1~U[0,K1), k2~U[0,K2):  E[x^2] = X^2/3."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                     + m * (k2_bound * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)

def gamma_gap(m, n, k1_bound, k2_bound):
    """uSVP gap of the PROJECTED instance: GH(pi(L)) / ||v'||."""
    return gh_from_logdet(2 * m + 1, proj_logdet(m, n, k1_bound, k2_bound)) \
        / planted_norm_proj(m, n, k1_bound, k2_bound)

# ---------------------------------------------------------------------------
# One trial, all three arms
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, arms=('ORIG', 'PROJ', 'CVP'),
          reduction='LLL', beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {}

    def reduce_(M, rows, cols):
        A = IntegerMatrix.from_matrix(M)
        if reduction == 'BKZ':
            BKZ.reduction(A, BKZ.Param(beta))
        else:
            LLL.reduction(A)
        return [[A[i][j] for j in range(cols)] for i in range(rows)]

    if 'ORIG' in arms:
        M = build_orig(sigs, n, lam, k1_bound, k2_bound)
        red = reduce_(M, 2 * m + 2, 2 * m + 2)
        pv = norm(planted_orig(sigs, d_secret, n, k1_bound, k2_bound))
        nz = [norm(r) for r in red if any(r)]
        out['ORIG'] = {
            'ok': recover_orig(red, m, n, S_KAN, d_secret) is not None,
            'pv': pv, 'sv': min(nz) if nz else float('nan'),
        }

    if 'PROJ' in arms:
        M = build_proj(sigs, n, lam, k1_bound, k2_bound)
        red = reduce_(M, 2 * m + 2, 2 * m + 1)
        pv = norm(planted_proj(sigs, n, k1_bound, k2_bound))
        nz = [norm(r) for r in red if any(r)]
        out['PROJ'] = {
            'ok': recover_proj(red, sigs, n, lam, k1_bound, k2_bound,
                               d_secret) is not None,
            'pv': pv, 'sv': min(nz) if nz else float('nan'),
        }

    if 'CVP' in arms:
        rows, a = build_cvp(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(rows)
        LLL.reduction(A)
        basis = [[A[i][j] for j in range(2 * m)] for i in range(len(rows))]
        basis = [r for r in basis if any(r)]
        u = babai_nearest_plane(basis, [-x for x in a])
        w = [a[k] + u[k] for k in range(2 * m)]
        ok = False
        if not any(w[i] % S_K1 for i in range(m)) and \
           not any(w[m + i] % S_K2 for i in range(m)):
            k1 = [w[i] // S_K1 for i in range(m)]
            k2 = [w[m + i] // S_K2 for i in range(m)]
            ok = d_from_k(k1, k2, sigs, n, lam, d_secret) == d_secret
        out['CVP'] = {'ok': ok, 'pv': float('nan'), 'sv': norm(w)}

    return out

def rate(curve, m, k1_bound, seeds, arms=('ORIG', 'PROJ', 'CVP'),
         reduction='LLL', beta=20):
    acc = {a: {'wins': 0, 'ratios': []} for a in arms}
    p, b, n, lam, G = curve
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, arms, reduction, beta)
        if res is None:
            continue
        for a in arms:
            acc[a]['wins'] += bool(res[a]['ok'])
            if res[a]['pv'] == res[a]['pv']:      # not NaN
                acc[a]['ratios'].append(res[a]['sv'] / res[a]['pv'])
    for a in arms:
        r = acc[a]['ratios']
        acc[a]['svpv'] = sum(r) / len(r) if r else float('nan')
    return acc

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

def search_curves(lo, hi, want=6):
    """Fresh j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, with lam*
    spread across (0, 0.5)."""
    nbins = want
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
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
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if bins[idx]:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(bins.values()):
            break
        p = int(sympy.nextprime(p))
    return [bins[i][0] for i in range(nbins) if bins[i]]


if __name__ == "__main__":
    # ===========================================================================
    print("=" * 78)
    print("Thread 23 - projected Phase-2 lattice: is the planted vector lambda_1?")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,             p,    b, n,    lam,  K1, m
        ("8-bit/199",        211,  2, 199,  106,  2,  6),
        ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
    ]

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), k1, m))


    # ===========================================================================
    # U1 - structural verification of the projection
    # ===========================================================================
    print("\n" + "-" * 78)
    print("U1: does pi kill the trivial vector, and what is det pi(L)?")
    print("-" * 78)
    print("pi deletes column m.  Claims: (a) rank pi(L) = 2m+1;")
    print("(b) det pi(L) = det L / (n*S_D);  (c) pi(v_planted) is in pi(L);")
    print("(d) sv/pv rises from ~0.4-0.6 (in L) to >= 1 (in pi(L)).")
    print()
    print(f"{'curve':<18} {'m':>2} {'K1':>3} {'rank':>5} {'det ok':>7} "
          f"{'sv/pv ORIG':>11} {'sv/pv PROJ':>11} {'gamma':>7}")

    for label, curve, k1_bound, m in hist:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        d_secret = random.Random(7777 + 42).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, 42)
        S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

        # (a) rank of pi(L)
        Mp = build_proj(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(Mp)
        LLL.reduction(A)
        red = [[A[i][j] for j in range(2 * m + 1)] for i in range(2 * m + 2)]
        nzrows = [r for r in red if any(r)]
        rank = len(nzrows)

        # (b) determinant, exactly, from the reduced square basis
        detp = abs(sympy.Matrix(nzrows).det()) if rank == 2 * m + 1 else None
        det_pred = (n * S_K1) ** m * S_K2 ** m * S_KAN // (n * S_D)
        det_ok = (detp == det_pred)

        # (c) membership of pi(v_planted): solve pi(L) x = v' over the rationals
        vp = planted_proj(sigs, n, k1_bound, k2_bound)
        sol = sympy.Matrix(nzrows).T.solve(sympy.Matrix(vp)) if rank == 2 * m + 1 else None
        member = sol is not None and all(x.q == 1 for x in sol)

        # (d) sv/pv, both lattices
        Mo = build_orig(sigs, n, lam, k1_bound, k2_bound)
        Ao = IntegerMatrix.from_matrix(Mo)
        LLL.reduction(Ao)
        redo = [[Ao[i][j] for j in range(2 * m + 2)] for i in range(2 * m + 2)]
        pvo = norm(planted_orig(sigs, d_secret, n, k1_bound, k2_bound))
        svo = min(norm(r) for r in redo if any(r))
        pvp = norm(vp)
        svp = min(norm(r) for r in nzrows)

        print(f"{label:<18} {m:>2} {k1_bound:>3} {rank:>5} {str(det_ok):>7} "
              f"{svo/pvo:>11.4f} {svp/pvp:>11.4f} "
              f"{gamma_gap(m, n, k1_bound, k2_bound):>7.3f}")
        assert member, f"pi(v_planted) not in pi(L) for {label}"

    print("\npi(v_planted) verified to be an integer combination of the pi(L) basis")
    print("on all three curves (assert passed).")

    # Where does the trivial vector go?
    print("\nTrivial vector t = n*S_D*e_m:  ||t|| = n*S_D in L, and pi(t) = 0.")
    for label, curve, k1_bound, m in hist:
        p, b, n, lam, G = curve
        S_K1, S_D, S_K2, S_KAN = scales(n, math.isqrt(n) + 1, math.isqrt(n) + 1)
        pn = planted_norm_proj(m, n, k1_bound, math.isqrt(n) + 1)
        print(f"  {label:<18} ||t||={n:>8}   E||v'||={pn:>12.1f}   "
              f"||t||/E||v'|| = {n/pn:.3f}")


    # ===========================================================================
    # U2 - head-to-head on the 2026-07-29 T4 grid
    # ===========================================================================
    print("\n" + "-" * 78)
    print("U2: T4 grid replication.  K1 wall, three arms, 5 seeds each.")
    print("-" * 78)
    print("2026-07-29 T4 (ORIG arm) found: 2557 wall at K1 ~ 12-16, 2677 at K1 ~ 4-6.")
    print("Falsifier: PROJ/CVP move the 2677 wall outward => real improvement.")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    u2_rows = []
    for label, curve, _k1, m in hist[1:]:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        print(f"\n{label}   n={n}  lam={lam}  lam*={lam_star(lam,n):.4f}  m={m}  "
              f"K2={k2b}")
        print(f"  {'K1':>3} {'eff':>6} {'gamma':>6} | {'ORIG':>6} {'PROJ':>6} "
              f"{'CVP':>6} | {'sv/pv O':>8} {'sv/pv P':>8}")
        for k1b in K1_GRID:
            acc = rate(curve, m, k1b, SEEDS)
            eff = k1b * k2b / n
            g = gamma_gap(m, n, k1b, k2b)
            print(f"  {k1b:>3} {eff:>6.3f} {g:>6.3f} | "
                  f"{acc['ORIG']['wins']}/{len(SEEDS):<4} "
                  f"{acc['PROJ']['wins']}/{len(SEEDS):<4} "
                  f"{acc['CVP']['wins']}/{len(SEEDS):<4} | "
                  f"{acc['ORIG']['svpv']:>8.4f} {acc['PROJ']['svpv']:>8.4f}")
            u2_rows.append((label, k1b, eff, g, acc['ORIG']['wins'],
                            acc['PROJ']['wins'], acc['CVP']['wins']))


    # ===========================================================================
    # U3 - T4b replication: more data at the wall
    # ===========================================================================
    print("\n" + "-" * 78)
    print("U3: T4b replication - does more data rescue the lam*=0.07 curve at K1=8?")
    print("-" * 78)
    print("2026-07-29 T4b (ORIG): m = 8/12/16/24/32 -> 0,0,1,0,1 out of 5.")

    label, curve, _k1, _m = hist[2]
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"\n{label}  K1=8  eff={8*k2b/n:.3f}")
    print(f"  {'m':>3} {'dim':>4} {'gamma':>6} | {'ORIG':>6} {'PROJ':>6} {'CVP':>6}")
    for m_try in [8, 12, 16, 24, 32]:
        acc = rate(curve, m_try, 8, SEEDS)
        print(f"  {m_try:>3} {2*m_try+1:>4} {gamma_gap(m_try, n, 8, k2b):>6.3f} | "
              f"{acc['ORIG']['wins']}/{len(SEEDS):<4} "
              f"{acc['PROJ']['wins']}/{len(SEEDS):<4} "
              f"{acc['CVP']['wins']}/{len(SEEDS):<4}")


    # ===========================================================================
    # U4 - is gamma the wall?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("U4: is the uSVP gap gamma = GH(pi(L))/||v'|| the wall?")
    print("-" * 78)
    print("Fresh 17-bit curves, m=12, eff swept.  gamma depends only on (n,m,K1,K2),")
    print("not on lam, so if gamma predicts success then lam is irrelevant (as T3")
    print("already concluded) and the wall is information-theoretic.")

    curves17 = search_curves(2 ** 16, 2 ** 17, want=6)
    print(f"\n{len(curves17)} curves found:")
    for (p, b, n, lam, G) in curves17:
        print(f"  p={p:<7} n={n:<7} lam={lam:<7} lam*={lam_star(lam,n):.4f}")

    M17 = 12
    print(f"\nm={M17}, 5 seeds/cell.  Cell entry = wins/5 (ORIG | PROJ).")
    hdr = "  " + f"{'curve (lam*)':<20}"
    EFFS = [0.05, 0.10, 0.15, 0.20, 0.25]
    for e in EFFS:
        hdr += f"{'eff=' + str(e):>13}"
    print(hdr)
    u4 = []
    for (p, b, n, lam, G) in curves17:
        k2b = math.isqrt(n) + 1
        line = "  " + f"{str(p) + ' (' + format(lam_star(lam,n), '.3f') + ')':<20}"
        for e in EFFS:
            k1b = max(2, int(round(e * n / k2b)))
            acc = rate((p, b, n, lam, G), M17, k1b, SEEDS, arms=('ORIG', 'PROJ'))
            line += f"{acc['ORIG']['wins']}/5 | {acc['PROJ']['wins']}/5".rjust(13)
            u4.append((p, lam_star(lam, n), e, gamma_gap(M17, n, k1b, k2b),
                       acc['ORIG']['wins'], acc['PROJ']['wins']))
        print(line)

    print(f"\n{'eff':>6} {'gamma (mean)':>13} {'ORIG wins':>10} {'PROJ wins':>10} "
          f"{'of':>4}")
    for e in EFFS:
        sub = [r for r in u4 if r[2] == e]
        print(f"{e:>6.2f} {sum(r[3] for r in sub)/len(sub):>13.3f} "
              f"{sum(r[4] for r in sub):>10} {sum(r[5] for r in sub):>10} "
              f"{5*len(sub):>4}")


    # ===========================================================================
    # U5 - BKZ at the wall
    # ===========================================================================
    print("\n" + "-" * 78)
    print("U5: BKZ(20) and BKZ(30) on pi(L) at the wall")
    print("-" * 78)
    print("If stronger reduction rescues cells that LLL misses, the wall is a")
    print("reduction-strength wall; if not, it is information-theoretic.")

    print(f"\n{'curve':<18} {'K1':>3} {'gamma':>6} | {'LLL':>6} {'BKZ20':>7} "
          f"{'BKZ30':>7}")
    for label, curve, _k1, m in hist[1:]:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1b in [6, 8, 12]:
            a1 = rate(curve, m, k1b, SEEDS, arms=('PROJ',), reduction='LLL')
            a2 = rate(curve, m, k1b, SEEDS, arms=('PROJ',), reduction='BKZ', beta=20)
            a3 = rate(curve, m, k1b, SEEDS, arms=('PROJ',), reduction='BKZ', beta=30)
            print(f"{label:<18} {k1b:>3} {gamma_gap(m,n,k1b,k2b):>6.3f} | "
                  f"{a1['PROJ']['wins']}/{len(SEEDS):<4} "
                  f"{a2['PROJ']['wins']}/{len(SEEDS):<5} "
                  f"{a3['PROJ']['wins']}/{len(SEEDS):<5}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
