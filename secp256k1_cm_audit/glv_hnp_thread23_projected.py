"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is
lambda_1, and test whether that moves the K1 wall.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 lattice L (dim 2m+2, built in glv_hnp_phase2_20bit.py:262)
  always contains the vector  n*S_D*e_m  of norm n, while the planted vector
  has norm ~ n*sqrt(2m/3 + 4/3).  The planted vector is therefore NEVER
  lambda_1, for any m >= 1 and any choice of S_D (both scale linearly in S_D).
  The 2026-07-29 run proposed removing that direction.

This script implements the two reformulations proposed there and adds a
quantitative success predictor.

  Variant A  (baseline)   L, dim 2m+2, Kannan embedding.        [as-is]
  Variant B  (projected)  pi_{e_m}(L): delete the d-column.     [proposal (a)]
                          rank 2m+1; d is recovered arithmetically from
                          (k1_0, k2_0) instead of read off a coordinate.
  Variant D  (CVP)        Babai nearest-plane on the projected congruence
                          lattice, no Kannan row.               [proposal (b)]

Experiments
  U1  correctness + geometry: is the planted vector in L_B, and is it now
      lambda_1?  Reports ||sv||/||pv|| for A and B on the historical curves.
  U2  the T4 K1-grid (2026-07-29) re-run for A / B / D.  Falsifier: does the
      K1 wall move outward?
  U3  Gaussian-heuristic predictor alpha = ||v_planted|| / GH(L_B) and its
      closed-form consequence -- an eff = K1*K2/n ceiling that does not
      depend on m.  Tested against a fresh 17-bit eff sweep.
  U4  alpha vs the 2026-07-29 separator nu_hat, and the combined predictor,
      scored by AUC on the pooled U2+U3 trials.

Run: python3 glv_hnp_thread23_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py:87)
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

def gauss_reduce_2d(u, v):
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def lambda_block_mu(n, lam, S_K1, S_K2):
    """|lambda_1| of the 2D block < (n*S_K1,0), (-lam*S_K1,S_K2) >."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1])

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

def search_curves(lo, hi, want, seed=12345):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    cur = build_curve(p, n_cand, seed=seed)
                    if cur is None:
                        continue
                    out.append((p, cur[1], n_cand, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Phase-2 lattice, variant A (verbatim from glv_hnp_phase2_20bit.py:262)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
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

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
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

def planted_vector_A(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_A(rows, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Variant B: projected lattice pi_{e_m}(L) -- the d-column is deleted.
# ---------------------------------------------------------------------------
#
# Column layout (2m+1 columns):
#   0 .. m-1        k1 block, scale S_K1
#   m .. 2m-1       k2 block, scale S_K2
#   2m              Kannan,   scale S_KANNAN
#
# The generating set is the image of L's basis, so it has 2m+2 generators of
# rank 2m+1 (n*row_m - sum_i B_i*row_i vanishes).  fplll's LLL handles the
# dependency and emits a zero row, which is dropped.
#
# d is no longer a coordinate; it is recovered from signature 0 by
#   d = B_0^{-1} (k1_0 + lam*k2_0 - A_0)  mod n.

def build_glv_lattice_B(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    cols = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                                   # n*S_K1*e_i
        r = [0] * cols; r[i] = n * S_K1; rows.append(r)
    r = [0] * cols                                       # pi(row_m): d-generator
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                   # k2 rows
        r = [0] * cols
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * cols                                       # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[cols - 1] = S_KANNAN
    rows.append(r)
    return rows

def planted_vector_B(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_B(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_0, k2_0) off any row with |Kannan coord| = S_KANNAN and solve
    for d.  Both signs are tried."""
    m = len(sigs)
    cols = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in rows:
        if abs(row[cols - 1]) != S_KAN:
            continue
        sign = 1 if row[cols - 1] > 0 else -1
        if sign * row[0] % S_K1 or sign * row[m] % S_K2:
            continue
        k1_0 = sign * row[0] // S_K1
        k2_0 = sign * row[m] // S_K2
        d_cand = B0inv * (k1_0 + lam * k2_0 - sigs[0]['A']) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Variant D: explicit CVP (Babai nearest plane) on the projected congruence
# lattice, no Kannan row.  rank 2m in 2m columns.
# ---------------------------------------------------------------------------

def build_congruence_lattice_D(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    cols = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * cols; r[i] = n * S_K1; rows.append(r)
    r = [0] * cols
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * cols
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    return rows

def babai_nearest_plane(basis, target):
    """Exact-rational Gram-Schmidt Babai on an LLL-reduced basis (list of
    int rows, full rank in its own span).  Returns the lattice vector."""
    from fractions import Fraction
    k = len(basis)
    bs = [[Fraction(x) for x in row] for row in basis]
    gs, mu = [], []
    for i in range(k):
        v = bs[i][:]
        mrow = []
        for j in range(i):
            den = sum(g * g for g in gs[j])
            c = Fraction(0) if den == 0 else sum(bs[i][t] * gs[j][t] for t in range(len(v))) / den
            mrow.append(c)
            v = [v[t] - c * gs[j][t] for t in range(len(v))]
        gs.append(v); mu.append(mrow)
    w = [Fraction(x) for x in target]
    coeffs = [0] * k
    for i in range(k - 1, -1, -1):
        den = sum(g * g for g in gs[i])
        if den == 0:
            continue
        c = sum(w[t] * gs[i][t] for t in range(len(w))) / den
        ci = int((c.numerator * 2 + c.denominator) // (2 * c.denominator)) if c >= 0 else \
             -int(((-c).numerator * 2 + c.denominator) // (2 * c.denominator))
        coeffs[i] = ci
        w = [w[t] - ci * bs[i][t] for t in range(len(w))]
    out = [0] * len(target)
    for i in range(k):
        if coeffs[i]:
            for t in range(len(target)):
                out[t] += coeffs[i] * basis[i][t]
    return out

# ---------------------------------------------------------------------------
# Reduction driver
# ---------------------------------------------------------------------------

def lll_rows(rows, use_bkz=False, beta=20):
    ncols = len(rows[0])
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    out = [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]
    return [r for r in out if any(r)]

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def run_trial(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False, beta=20):
    """variant in {'A','B','D'}.  Returns (ok, extra_dict)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, {}
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

    if variant == 'A':
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        red = lll_rows(M, use_bkz, beta)
        pv = planted_vector_A(sigs, d_secret, n, k1_bound, k2_bound)
        ok = recover_d_A(red, m, n, S_KAN, d_secret) is not None
    elif variant == 'B':
        M = build_glv_lattice_B(sigs, n, lam, k1_bound, k2_bound)
        red = lll_rows(M, use_bkz, beta)
        pv = planted_vector_B(sigs, n, k1_bound, k2_bound)
        ok = recover_d_B(red, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
    elif variant == 'D':
        M = build_congruence_lattice_D(sigs, n, lam, k1_bound, k2_bound)
        red = lll_rows(M, use_bkz, beta)
        tgt = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
        cv = babai_nearest_plane(red, tgt)
        # For the exact CVP answer, tgt - cv = -(k1_i*S_K1, k2_i*S_K2), so the
        # readout is negated.  Both signs are tried anyway.
        diff = [tgt[t] - cv[t] for t in range(2 * m)]
        ok = False
        B0inv = modinv(sigs[0]['B'], n)
        for sgn in (-1, 1):
            if (sgn * diff[0]) % S_K1 or (sgn * diff[m]) % S_K2:
                continue
            k1_0 = sgn * diff[0] // S_K1
            k2_0 = sgn * diff[m] // S_K2
            if B0inv * (k1_0 + lam * k2_0 - sigs[0]['A']) % n == d_secret:
                ok = True
                break
        pv = [sigs[i]['k1'] * S_K1 for i in range(m)] + \
             [sigs[i]['k2'] * S_K2 for i in range(m)]
        red = red  # for sv reporting
    else:
        raise ValueError(variant)

    sv = min((norm(r) for r in red), default=float('inf'))
    return ok, {'sv': sv, 'pv': norm(pv), 'ratio': sv / norm(pv)}

def success_rate(curve, m, k1_bound, seeds, variant, use_bkz=False, beta=20):
    wins, tot, ratios = 0, 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        ok, ex = run_trial(curve, m, d_trial, k1_bound, seed, variant, use_bkz, beta)
        if ok is None:
            continue
        tot += 1; wins += ok
        if ex: ratios.append(ex['ratio'])
    return wins, tot, (sum(ratios) / len(ratios) if ratios else float('nan'))

# ---------------------------------------------------------------------------
# Gaussian-heuristic predictor for variant B
# ---------------------------------------------------------------------------

def gh_alpha(n, m, k1_bound, k2_bound):
    """alpha = E||v_planted|| / GH(L_B).  Success predicted iff alpha < 1.

    rank r = 2m+1
    det L_B = (n*S_K1)^m * S_K2^m * S_KANNAN / n     (L_B = pi_{e_m}(L),
              det pi(L) = det(L) / (n*S_D), S_D = 1)
    """
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    r = 2 * m + 1
    log_det = (m * (math.log(n) + math.log(S_K1)) + m * math.log(S_K2)
               + math.log(S_KAN) - math.log(n))
    gh = math.sqrt(r / (2 * math.pi * math.e)) * math.exp(log_det / r)
    pv = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                   + m * (k2_bound * S_K2) ** 2 / 3.0
                   + S_KAN ** 2)
    return pv / gh, pv, gh

def gh_C(m):
    """C(m) = sqrt( (2m+1) / (2*pi*e*(2m/3 + 1)) ).  C(m) -> sqrt(3/(2*pi*e))."""
    return math.sqrt((2 * m + 1) / (2 * math.pi * math.e * (2 * m / 3.0 + 1)))

def eff_ceiling(m, n):
    """Largest eff = K1*K2/n with alpha < 1 at this (m, n), in the idealised
    regime S_K1*K1 = S_K2*K2 = S_KANNAN = n:

      alpha < 1  <=>  n * eff^m  <  C(m)^(2m+1)
                 <=>  eff  <  C(m)^((2m+1)/m) * n^(-1/m).

    The n^(-1/m) factor -> 1 as m grows, so the ceiling increases with m but is
    capped: sup_m C(m)^((2m+1)/m) = C(inf)^2 = 3/(2*pi*e) = 0.1756.  No number
    of signatures lifts eff past that.
    """
    return gh_C(m) ** ((2 * m + 1) / m) * n ** (-1.0 / m)

def m_for_eff(eff, n, m_max=4000):
    """Smallest m with alpha < 1 at this eff; None if eff is above the cap."""
    for m in range(2, m_max):
        if eff < eff_ceiling(m, n):
            return m
    return None

C_INF = 3.0 / (2 * math.pi * math.e)

def auc(scores, labels):
    """AUC with score oriented so that LOWER score => success (label 1)."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    wins = 0.0
    for a in pos:
        for b in neg:
            wins += 1.0 if a < b else (0.5 if a == b else 0.0)
    return wins / (len(pos) * len(neg))


# ===========================================================================
print("=" * 78)
print("Thread 23 — projected GLV-HNP Phase-2 lattice (planted vector as lambda_1)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,          p,    b, n,    lam,  K1, m
    ("8-bit/199",     211,  2, 199,  106,  2,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",   2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U1: is the planted vector in L_B, and is it now lambda_1?")
print("-" * 78)
print("sv = shortest row after LLL, pv = planted vector.  sv/pv >= 1 means the")
print("planted vector is (up to LLL approximation) the shortest vector.")
print()
print(f"{'curve':<14} {'K1':>3} {'m':>3} | {'A: sv/pv':>9} {'A ok':>5} | "
      f"{'B: sv/pv':>9} {'B ok':>5} | {'B in L?':>7}")

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2, 42)
    # membership check for variant B: solve pv = x * basis over the generators
    rowsB = build_glv_lattice_B(sigs, n, lam, k1, k2)
    pvB = planted_vector_B(sigs, n, k1, k2)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    # explicit combination: kannan + d*pi(row_m) + sum k2_i*row_{k2,i} + c_i*n e_i
    comb = [0] * (2 * m + 1)
    for t in range(2 * m + 1):
        comb[t] += rowsB[m][t] * d0 + rowsB[2 * m + 1][t]
        for i in range(m):
            comb[t] += rowsB[m + 1 + i][t] * sigs[i]['k2']
    for i in range(m):
        delta = pvB[i] - comb[i]
        assert delta % (n * S_K1) == 0, (label, i, delta)
        comb[i] += delta
    inL = (comb == pvB)

    # variant-D sanity: with an EXACT CVP answer, does the readout recover d?
    rowsD = build_congruence_lattice_D(sigs, n, lam, k1, k2)
    tgt = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    cv = [0] * (2 * m)
    for t in range(2 * m):
        cv[t] += rowsD[m][t] * d0
        for i in range(m):
            cv[t] += rowsD[m + 1 + i][t] * sigs[i]['k2']
    for i in range(m):
        delta = (sigs[i]['k1'] - sigs[i]['A']) * S_K1 - cv[i]
        assert delta % (n * S_K1) == 0
        cv[i] += delta
    err = [tgt[t] - cv[t] for t in range(2 * m)]
    B0inv = modinv(sigs[0]['B'], n)
    d_ex = B0inv * (-err[0] // S_K1 + lam * (-err[m] // S_K2) - sigs[0]['A']) % n
    d_ok = (d_ex == d0)

    rA = run_trial(curve, m, d0, k1, 42, 'A')
    rB = run_trial(curve, m, d0, k1, 42, 'B')
    print(f"{label:<14} {k1:>3} {m:>3} | {rA[1]['ratio']:>9.3f} {str(rA[0]):>5} | "
          f"{rB[1]['ratio']:>9.3f} {str(rB[0]):>5} | {str(inL):>7} | "
          f"D-exact-CVP ok={d_ok}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U2: the 2026-07-29 T4 K1-grid, variants A / B / D")
print("-" * 78)
print("Falsifier (2026-07-29): if the K1 wall moves outward under B or D the")
print("reformulation is a real improvement; if it stays, the wall is")
print("information-theoretic and Phase 2 is at its ceiling.")

K1_GRID = [1, 2, 3, 4, 6, 8, 12, 16, 24]
u2_rows = []
for label, curve, _k1, m in hist[1:]:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    print(f"\n{label}: n={n}, lam={lam}, lam*={lam_star(lam,n):.3f}, K2={k2}, m={m}")
    print(f"{'K1':>4} {'eff':>7} {'alpha':>7} | {'A':>5} {'B':>5} {'D':>5} | "
          f"{'B sv/pv':>8}")
    for k1 in K1_GRID:
        al, pv, gh = gh_alpha(n, m, k1, k2)
        wa, ta, _ = success_rate(curve, m, k1, SEEDS, 'A')
        wb, tb, rb = success_rate(curve, m, k1, SEEDS, 'B')
        wd, td, _ = success_rate(curve, m, k1, SEEDS, 'D')
        print(f"{k1:>4} {k1*k2/n:>7.4f} {al:>7.3f} | {wa}/{ta:<3} {wb}/{tb:<3} "
              f"{wd}/{td:<3} | {rb:>8.3f}")
        u2_rows.append({'label': label, 'n': n, 'lam': lam, 'm': m, 'K1': k1,
                        'K2': k2, 'eff': k1 * k2 / n, 'alpha': al,
                        'A': wa / max(ta, 1), 'B': wb / max(tb, 1),
                        'D': wd / max(td, 1)})

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U3: the Gaussian-heuristic eff ceiling")
print("-" * 78)
print("alpha = E||v_planted|| / GH(L_B) < 1  <=>  n*eff^m < C(m)^(2m+1),")
print(f"C(m) -> sqrt(3/(2*pi*e)); eff ceiling -> 3/(2*pi*e) = {C_INF:.4f}.")
print()
print(f"{'m':>4} {'C(m)':>7} {'C^((2m+1)/m)':>13} {'ceiling@n=2^17':>15} "
      f"{'ceiling@n=2^256':>16}")
for m in (4, 6, 8, 10, 12, 16, 24, 48, 64, 256, 1024):
    print(f"{m:>4} {gh_C(m):>7.4f} {gh_C(m)**((2*m+1)/m):>13.4f} "
          f"{eff_ceiling(m, 1 << 17):>15.4f} {eff_ceiling(m, 1 << 256):>16.4f}")
print(f"\nsup over m = C(inf)^2 = 3/(2*pi*e) = {C_INF:.4f}  (hard cap, any n)")

print("\n17-bit sweep: 8 fresh curves x eff in {0.05,0.10,0.15,0.20,0.25}, m=12")
curves17 = search_curves(1 << 16, 1 << 17, 8)
print(f"  found {len(curves17)} curves: " +
      ", ".join(f"n={c[2]}(lam*={lam_star(c[3],c[2]):.3f})" for c in curves17))
M17 = 12
u3_rows = []
print(f"\n{'eff':>6} {'alpha':>7} | {'B: curves 5/5':>14} {'B: total':>10} "
      f"{'A: total':>10}")
for eff_t in (0.05, 0.10, 0.15, 0.20, 0.25):
    tot_b = tot_a = trials = full_b = 0
    alphas = []
    for c in curves17:
        n = c[2]
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(round(eff_t * n / k2)))
        al, _, _ = gh_alpha(n, M17, k1, k2)
        alphas.append(al)
        wb, tb, _ = success_rate(c, M17, k1, SEEDS, 'B')
        wa, ta, _ = success_rate(c, M17, k1, SEEDS, 'A')
        tot_b += wb; tot_a += wa; trials += tb
        full_b += (wb == tb)
        u3_rows.append({'n': n, 'lam': c[3], 'm': M17, 'K1': k1, 'K2': k2,
                        'eff': k1 * k2 / n, 'alpha': al,
                        'A': wa / max(ta, 1), 'B': wb / max(tb, 1)})
    print(f"{eff_t:>6.2f} {sum(alphas)/len(alphas):>7.3f} | "
          f"{full_b}/{len(curves17):<12} {tot_b}/{trials:<8} {tot_a}/{trials:<8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U5: does MORE data lift eff past the ceiling?")
print("-" * 78)
print("The ceiling claim is falsifiable: for eff > 3/(2*pi*e) = 0.1756 no m")
print("should ever work.  Sweep m at fixed eff on 17-bit curves, variant B.")
print("Predicted smallest m with alpha < 1 (n ~ 2^17) shown as m*.")

M_GRID = [12, 24, 36, 48]
EFF_GRID = [0.10, 0.125, 0.15, 0.20]
c5 = curves17[:5]
seeds5 = SEEDS[:3]
u5_rows = []
print(f"\n{'eff':>6} {'m*':>5} | " + " ".join(f"{'m='+str(mm):>10}" for mm in M_GRID))
for eff_t in EFF_GRID:
    mstar = m_for_eff(eff_t, 1 << 17)
    cells = []
    for mm in M_GRID:
        wins = tot = 0
        als = []
        for c in c5:
            n = c[2]; k2 = math.isqrt(n) + 1
            k1 = max(2, int(round(eff_t * n / k2)))
            al, _, _ = gh_alpha(n, mm, k1, k2)
            als.append(al)
            w, t, _ = success_rate(c, mm, k1, seeds5, 'B')
            wins += w; tot += t
            u5_rows.append({'n': n, 'lam': c[3], 'm': mm, 'K1': k1, 'K2': k2,
                            'eff': k1 * k2 / n, 'alpha': al,
                            'B': w / max(t, 1)})
        cells.append(f"{wins}/{tot} a={sum(als)/len(als):.2f}")
    ms = str(mstar) if mstar else ">cap"
    print(f"{eff_t:>6.3f} {ms:>5} | " + " ".join(f"{c:>10}" for c in cells))

print("\nper-curve breakdown (wins over 3 seeds, summed over the 4 m values,")
print("max 12), against lam* and nu_hat = lambda_1(L2)/sqrt(det L2):")
print(f"\n{'n':>7} {'lam*':>6} {'nu_hat':>7} | " +
      " ".join(f"{'eff='+format(e,'.3f'):>10}" for e in EFF_GRID))
for c in c5:
    n = c[2]; k2 = math.isqrt(n) + 1
    k1_ref = max(2, int(round(0.10 * n / k2)))
    S_K1, _sd, S_K2, _sk = scales(n, k1_ref, k2)
    nh = lambda_block_mu(n, c[3], S_K1, S_K2) / math.sqrt(n * S_K1 * S_K2)
    cells = []
    for e in EFF_GRID:
        tot = sum(r['B'] * len(seeds5) for r in u5_rows
                  if r['n'] == n and abs(r['eff'] - e) < 0.02)
        cells.append(f"{int(round(tot))}/{len(M_GRID)*len(seeds5)}")
    print(f"{n:>7} {lam_star(c[3], n):>6.3f} {nh:>7.3f} | " +
          " ".join(f"{x:>10}" for x in cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U4: predictor comparison (AUC, lower score = predicted success)")
print("-" * 78)
allrows = u2_rows + u3_rows + u5_rows
lab = [r['B'] >= 0.8 for r in allrows]
sc_alpha = [r['alpha'] for r in allrows]
sc_eff = [r['eff'] for r in allrows]
nuhat, comb = [], []
for r in allrows:
    S_K1, _sd, S_K2, _sk = scales(r['n'], r['K1'], r['K2'])
    mu = lambda_block_mu(r['n'], r['lam'], S_K1, S_K2)
    det2 = math.sqrt(r['n'] * S_K1 * S_K2)
    # 2026-07-29 orientation: LOW nu_hat predicted success, so score = nu_hat.
    nuhat.append(mu / det2)
    _al, pv, gh = gh_alpha(r['n'], r['m'], r['K1'], r['K2'])
    comb.append(math.log(_al) + 0.5 * math.log(mu / det2))
print(f"  n trials      : {len(allrows)}  ({sum(lab)} success, "
      f"{len(allrows)-sum(lab)} fail)")
print(f"  AUC(alpha)    : {auc(sc_alpha, lab):.3f}")
print(f"  AUC(eff)      : {auc(sc_eff, lab):.3f}")
print(f"  AUC(nu_hat)   : {auc(nuhat, lab):.3f}   (low nu_hat = success, "
      f"2026-07-29 orientation)")
print(f"  AUC(alpha*sqrt(nu_hat)) : {auc(comb, lab):.3f}")
n_al = sum(1 for r, l in zip(allrows, lab) if (r['alpha'] < 1) == l)
print(f"  alpha<1 as a hard rule: {n_al}/{len(allrows)} correct "
      f"({100.0*n_al/len(allrows):.1f}%)")
fp = [(r['label'] if 'label' in r else f"n={r['n']}", r['K1'], r['alpha'], r['B'])
      for r, l in zip(allrows, lab) if r['alpha'] < 1 and not l]
fn = [(r['label'] if 'label' in r else f"n={r['n']}", r['K1'], r['alpha'], r['B'])
      for r, l in zip(allrows, lab) if r['alpha'] >= 1 and l]
print(f"  alpha<1 but failed : {len(fp)}  {fp[:6]}")
print(f"  alpha>=1 but worked: {len(fn)}  {fn[:6]}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
