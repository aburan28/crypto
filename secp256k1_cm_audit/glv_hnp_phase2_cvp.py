"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
the target of the search.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5):
  In the Kannan/SVP formulation of `build_glv_lattice` the shortest vector is
  ALWAYS the trivial vector n*S_D*e_m (100% of its energy in the d-column,
  |sv[m]|/n = 1.0000 on every curve tested).  Algebraically, with
  S_K1*K1 ~ S_K2*K2 ~ S_KANNAN ~ n and S_D = 1,
      ||v_planted||^2 ~ n^2 (2m/3 + 4/3)   vs   ||n*S_D*e_m||^2 = n^2,
  so the planted vector is never lambda_1 for any m >= 1.  Recovery is
  therefore not an SVP condition but a BDD/coset condition.

This script replaces the SVP formulation with a genuine CVP:

  Lattice L (dim 2m+1, no Kannan row/column), rows
    q rows   i=0..m-1 :  n*S_K1*F * e_i
    d row              :  (B_i*S_K1*F)_i  (+) 0_m  (+) S_D
    k2 rows  i=0..m-1 : (-lam*S_K1*F)*e_i (+) (S_K2*F)*e_i (+) 0
  Target
    t = (-A_i*S_K1*F)_i (+) 0_m (+) 0
  Then for v = d*(d row) + sum k2_i*(k2 row i) - sum q_i*(q row i),
    v - t = (k1_i*S_K1*F)_i (+) (k2_i*S_K2*F)_i (+) (d*S_D),
  i.e. the planted vector is exactly the CVP error vector, and d is read off
  the last coordinate mod n.  The trivial vector n*S_D*e_{2m} still lies in L
  but is harmless here: it only shifts d by n, which is the problem's own
  redundancy (d is defined mod n).

Three solvers are compared on the same instances:
  V0  Kannan + LLL      (the 2026-06-15 .. 2026-07-29 baseline)
  V1  CVP + Babai nearest plane (exact rational GS on the LLL-reduced basis)
  V2  CVP + fpylll enumeration (CVP.closest_vector) -- the information-
      theoretic reference: if V2 misses, no lattice algorithm can win.

The V2/planted comparison is the sharp diagnostic the 2026-07-29 entry asked
for: if the enumerated closest vector is strictly closer to t than the planted
vector is, the instance is information-theoretically lost (not enough data);
if the planted vector IS the closest vector but Babai/LLL miss it, the wall is
a reduction-quality wall and the reformulation can move it.

Run: python3 glv_hnp_phase2_cvp.py
"""

import math
import random
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_20bit.py so the instances are bit-identical to prior runs.
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
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def lam_star(lam, n):
    return min(lam % n, (-lam) % n) / n

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p."""
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
    return (r1, r2)

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

def search_fresh_curves(lo, hi, want, seed=1):
    """Fresh j=0 GLV curves with n prime, spread over lam*."""
    import sympy
    out, seen_ls = [], []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    ls = lam_star(lam, n_cand)
                    if any(abs(ls - s) < 0.06 for s in seen_ls):
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    seen_ls.append(ls)
                    out.append((f"{n_cand.bit_length()}b/{p}",
                                (p, cur[1], n_cand, lam, cur[4])))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Instance generation -- verbatim from glv_hnp_phase2_lambda_threshold.py
# so that (curve, m, K1, seed) reproduces the same signatures as prior runs.
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# V0: the historical Kannan/SVP lattice
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, S_D=1):
    return (n // k1_bound, S_D, max(1, n // k2_bound), n)

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound, S_D)
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

def build_glv_lattice_centered(sigs, n, lam, k1_bound, k2_bound, S_D=1,
                               center_d=False):
    """
    The SAME Kannan/SVP lattice, but with the last (Kannan) row shifted to the
    midpoints of the k1 and k2 boxes, so the planted vector carries
    k1_i - K1/2 and k2_i - K2/2 instead of k1_i and k2_i.

    This is the ablation that separates 'CVP reformulation' from 'box
    centering': if centering alone recovers the CVP arm's gain, the
    reformulation is not the operative ingredient.
    """
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound, S_D)
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
        M[2 * m + 1][i] = (sigs[i]['A'] - k1_bound // 2) * S_K1
        M[2 * m + 1][m + 1 + i] = -(k2_bound // 2) * S_K2
    if center_d:
        M[2 * m + 1][m] = -(n // 2) * S_D
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector_kannan_centered(sigs, d_secret, n, k1_bound, k2_bound,
                                   S_D=1, center_d=False):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - k1_bound // 2) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - k2_bound // 2) * S_K2
    v[m] = (d_secret - (n // 2 if center_d else 0)) * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_kannan_centered(M_reduced, m, n, S_KANNAN, S_D, d_secret,
                              center_d=False):
    dim = 2 * m + 2
    off = (n // 2) if center_d else 0
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % S_D: continue
        d_cand = (val // S_D + off) % n
        if d_cand == d_secret:
            return d_cand
    return None

def run_v0_centered(curve, m, d_secret, k1_bound, seed=42, S_D=1,
                    center_d=False, use_bkz=False, beta=20):
    from fpylll import BKZ
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice_centered(sigs, n, lam, k1_bound, k2_bound, S_D,
                                   center_d)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_K1, S_Dv, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    ok = recover_d_kannan_centered(reduced, m, n, S_KAN, S_Dv, d_secret,
                                   center_d) is not None
    pv = planted_vector_kannan_centered(sigs, d_secret, n, k1_bound, k2_bound,
                                        S_D, center_d)
    pn = math.sqrt(sum(x * x for x in pv))
    sv = min(reduced, key=lambda r: sum(x * x for x in r))
    sn = math.sqrt(sum(x * x for x in sv))
    return (ok, pn, sn, 0.0)

def planted_vector_kannan(sigs, d_secret, n, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_kannan(M_reduced, m, n, S_KANNAN, S_D, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        if (sign * row[m]) % S_D != 0: continue
        d_cand = ((sign * row[m]) // S_D) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def run_v0(curve, m, d_secret, k1_bound, seed=42, S_D=1, use_bkz=False, beta=20):
    """Returns (ok, planted_norm, shortest_norm, frac_energy_in_d_col)."""
    from fpylll import BKZ
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, S_D)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_K1, S_Dv, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    ok = recover_d_kannan(reduced, m, n, S_KAN, S_Dv, d_secret) is not None
    pv = planted_vector_kannan(sigs, d_secret, n, k1_bound, k2_bound, S_D)
    pn = math.sqrt(sum(x * x for x in pv))
    sv = min(reduced, key=lambda r: sum(x * x for x in r))
    sn = math.sqrt(sum(x * x for x in sv))
    e_d = (sv[m] * sv[m]) / sum(x * x for x in sv) if sn else 0.0
    return (ok, pn, sn, e_d)

# ---------------------------------------------------------------------------
# V1/V2: the CVP reformulation
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound, F=1, S_D=1):
    """dim 2m+1.  Columns: [0,m) congruence, [m,2m) k2, {2m} d."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1 = (n // k1_bound) * F
    S_K2 = max(1, n // k2_bound) * F
    B = [[0] * dim for _ in range(dim)]
    for i in range(m):                        # q rows
        B[i][i] = n * S_K1
    for i in range(m):                        # d row
        B[m][i] = sigs[i]['B'] * S_K1
    B[m][2 * m] = S_D
    for i in range(m):                        # k2 rows
        B[m + 1 + i][i] = -lam * S_K1
        B[m + 1 + i][m + i] = S_K2
    return B, S_K1, S_K2

def cvp_target(sigs, S_K1, S_K2=None, k1_bound=None, k2_bound=None,
               centered=False):
    """
    Uncentered: t_i = -A_i*S_K1, so the error carries k1_i in [0,K1) and
    k2_i in [0,K2) -- one-sided boxes, E[x^2] = K^2/3.

    Centered: shift the target to the box midpoints so the error carries
    k1_i - K1/2 and k2_i - K2/2, i.e. E[x^2] = K^2/12.  That is a factor-2
    reduction in the expected error norm for free; the prior Phase-2
    construction never did it.
    """
    m = len(sigs)
    t = [0] * (2 * m + 1)
    for i in range(m):
        t[i] = -(sigs[i]['A'] * S_K1)
    if centered:
        for i in range(m):
            t[i] += (k1_bound // 2) * S_K1
            t[m + i] = (k2_bound // 2) * S_K2
    return t

def cvp_planted_error(sigs, d_secret, S_K1, S_K2, n, S_D=1):
    """
    The error vector v_planted - t that the CVP must find.

    The d-column entry is taken at the BALANCED representative of d mod n:
    n*S_D*e_{2m} is a lattice vector, so d and d-n label the same solution and
    the shorter of the two is the one a CVP solver will actually return.
    Comparing against the unbalanced d would report 'a strictly closer vector
    exists' for every instance with d > n/2, which is a redundancy of the
    encoding and not a failure.
    """
    m = len(sigs)
    e = [0] * (2 * m + 1)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1
        e[m + i] = sigs[i]['k2'] * S_K2
    d_bal = d_secret - n if 2 * d_secret > n else d_secret
    e[2 * m] = d_bal * S_D
    return e

def cvp_planted_error_centered(sigs, d_secret, S_K1, S_K2, n,
                               k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    e = [0] * (2 * m + 1)
    for i in range(m):
        e[i] = (sigs[i]['k1'] - k1_bound // 2) * S_K1
        e[m + i] = (sigs[i]['k2'] - k2_bound // 2) * S_K2
    d_bal = d_secret - n if 2 * d_secret > n else d_secret
    e[2 * m] = d_bal * S_D
    return e

def decode_error(err, m, S_K1, S_K2, S_D, n, k1_bound, k2_bound, centered):
    """Decode a CVP error vector back into (k1_i, k2_i, d).  Returns None if
    the vector is not on the expected sublattice grid."""
    k1s, k2s = [], []
    for i in range(m):
        if err[i] % S_K1 or err[m + i] % S_K2:
            return None
        a = err[i] // S_K1
        b = err[m + i] // S_K2
        if centered:
            a += k1_bound // 2
            b += k2_bound // 2
        k1s.append(a)
        k2s.append(b)
    if err[2 * m] % S_D:
        return None
    d = (err[2 * m] // S_D) % n
    return k1s, k2s, d

def out_of_box(k1s, k2s, k1_bound, k2_bound):
    """How many decoded coordinates fall outside the true search boxes."""
    bad1 = sum(1 for x in k1s if not (0 <= x < k1_bound))
    bad2 = sum(1 for x in k2s if not (0 <= x < k2_bound))
    return bad1, bad2

# ---- exact rational Gram-Schmidt + Babai nearest plane --------------------

def gso_exact(B):
    rows = len(B)
    Bs, norms = [], []
    for i in range(rows):
        v = [Fraction(x) for x in B[i]]
        for j in range(i):
            if norms[j] == 0:
                continue
            num = sum(Fraction(B[i][k]) * Bs[j][k] for k in range(len(B[i])))
            mu = num / norms[j]
            if mu:
                v = [v[k] - mu * Bs[j][k] for k in range(len(v))]
        Bs.append(v)
        norms.append(sum(x * x for x in v))
    return Bs, norms

def babai_nearest_plane(B, Bs, norms, t):
    """Exact Babai nearest plane; B should already be LLL-reduced."""
    b = [Fraction(x) for x in t]
    coeffs = [0] * len(B)
    for i in reversed(range(len(B))):
        if norms[i] == 0:
            continue
        c = sum(b[k] * Bs[i][k] for k in range(len(b))) / norms[i]
        ci = int((c.numerator * 2 + c.denominator) // (2 * c.denominator))
        coeffs[i] = ci
        if ci:
            b = [b[k] - ci * Fraction(B[i][k]) for k in range(len(b))]
    return coeffs

def mat_vec(coeffs, B):
    dim = len(B[0])
    out = [0] * dim
    for i, c in enumerate(coeffs):
        if not c:
            continue
        Bi = B[i]
        for k in range(dim):
            if Bi[k]:
                out[k] += c * Bi[k]
    return out

def nrm(v):
    return math.sqrt(float(sum(x * x for x in v)))

def run_cvp(curve, m, d_secret, k1_bound, seed=42, F=1, S_D=1,
            do_enum=True, enum_dim_cap=31, bkz_beta=None, centered=False,
            want_decode=False):
    """
    Returns dict with:
      ok_babai   d recovered by Babai nearest plane
      ok_enum    d recovered by fpylll CVP enumeration
      d_planted  ||planted error|| (balanced d)
      d_babai    ||Babai error||
      d_enum     ||enumerated error||
      planted_is_closest  bool -- enum found nothing strictly closer than the
                          planted error.  False => information-theoretic loss.

    F inflates S_K1 and S_K2 relative to S_D, which makes the d-column's
    contribution to the CVP distance negligible.  F=1 penalises large |d|
    (the target has 0 in the d-column), which is not the objective we want;
    F >> 1 is the faithful BDD.
    """
    from fpylll import BKZ
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    B, S_K1, S_K2 = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound, F, S_D)
    dim = 2 * m + 1
    t = cvp_target(sigs, S_K1, S_K2, k1_bound, k2_bound, centered)
    if centered:
        e_planted = cvp_planted_error_centered(sigs, d_secret, S_K1, S_K2, n,
                                               k1_bound, k2_bound, S_D)
    else:
        e_planted = cvp_planted_error(sigs, d_secret, S_K1, S_K2, n, S_D)

    A = IntegerMatrix.from_matrix(B)
    if bkz_beta:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    Bred = [[A[i][j] for j in range(dim)] for i in range(dim)]

    Bs, norms = gso_exact(Bred)
    coeffs = babai_nearest_plane(Bred, Bs, norms, t)
    v_bab = mat_vec(coeffs, Bred)
    err_bab = [v_bab[k] - t[k] for k in range(dim)]
    d_bab = (err_bab[2 * m] // S_D) % n if err_bab[2 * m] % S_D == 0 else None
    ok_bab = (d_bab == d_secret)

    res = {
        'ok_babai': ok_bab,
        'd_planted': nrm(e_planted),
        'd_babai': nrm(err_bab),
        'ok_enum': None, 'd_enum': None, 'planted_is_closest': None,
    }

    if do_enum and dim <= enum_dim_cap:
        try:
            w = CVP.closest_vector(A, tuple(int(x) for x in t))
            v_enum = list(w)
            err_enum = [v_enum[k] - t[k] for k in range(dim)]
            d_en = (err_enum[2 * m] // S_D) % n if err_enum[2 * m] % S_D == 0 else None
            res['ok_enum'] = (d_en == d_secret)
            res['d_enum'] = nrm(err_enum)
            res['planted_is_closest'] = (
                sum(x * x for x in err_enum) >= sum(x * x for x in e_planted))
            if want_decode:
                res['decode_enum'] = decode_error(
                    err_enum, m, S_K1, S_K2, S_D, n, k1_bound, k2_bound,
                    centered)
        except Exception as exc:            # enumeration failure is informative
            res['enum_error'] = repr(exc)
    return res

# ---------------------------------------------------------------------------
# Curves: the two historical 12-bit curves that define the K1 wall (T4)
# ---------------------------------------------------------------------------

HIST = [
    # label,            p,    b, n,    lam
    ("8-bit/199",       211,  2, 199,  106),
    ("12-bit/2557",     2557, 2, 2659, 1755),
    ("12-bit/2677",     2677, 2, 2647, 185),
]

curves = {}
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves[label] = (p, b, n, lam, G)

SEEDS = [42, 1234, 9999, 555, 31337]

def d_for_seed(seed, n):
    return random.Random(seed + 7777).randint(1, n - 1)

print("=" * 78)
print("Thread 23 -- CVP reformulation of the GLV-HNP Phase-2 lattice")
print("=" * 78)
for label in curves:
    p, b, n, lam, G = curves[label]
    print(f"  {label:<14} p={p:<6} n={n:<6} ({n.bit_length()}b) "
          f"lam={lam:<6} lam*={lam_star(lam, n):.4f} K2={math.isqrt(n)+1}")

# ===========================================================================
# EXP C1 -- structural sanity: planted vector is in L, and the trivial
#           vector of the Kannan form is confirmed as lambda_1
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C1: structural sanity of the two formulations")
print("-" * 78)
print(f"{'curve':<14} {'K1':>3} {'m':>3} {'||pv||_kan':>12} {'||sv||_kan':>12} "
      f"{'sv/pv':>7} {'E_d':>6} {'||e||_cvp':>12} {'ratio':>7}")
for label, k1, m in [("8-bit/199", 2, 6), ("12-bit/2557", 8, 8),
                     ("12-bit/2677", 8, 10)]:
    curve = curves[label]
    n = curve[2]
    d0 = d_for_seed(42, n)
    v0 = run_v0(curve, m, d0, k1, seed=42)
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(curve[4], d0, m, n, curve[3], curve[0], k1, k2b, 42)
    _, S_K1, S_K2 = build_cvp_lattice(sigs, n, curve[3], k1, k2b)
    e = cvp_planted_error(sigs, d0, S_K1, S_K2, n)
    ok, pn, sn, e_d = v0
    print(f"{label:<14} {k1:>3} {m:>3} {pn:>12.1f} {sn:>12.1f} "
          f"{sn/pn:>7.3f} {e_d:>6.3f} {nrm(e):>12.1f} {nrm(e)/pn:>7.3f}")
print("\nE_d = fraction of the Kannan shortest vector's energy in the d-column.")
print("E_d = 1.000 reproduces the 2026-07-29 T5 result (sv is always n*S_D*e_m).")
print("||e||_cvp / ||pv||_kan < 1: dropping the Kannan column and the d-scale")
print("shortens the object being searched for.")

# ===========================================================================
# EXP C2 -- the K1 wall: V0 vs V1(Babai) vs V2(enumeration)
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C2: K1 wall -- Kannan+LLL vs CVP at F=1 and F=2^12   (m=10, 5 seeds)")
print("-" * 78)
print("F=1 penalises large |d| (target has 0 in the d-column); F=2^12 makes")
print("the d-column free, which is the faithful BDD objective.")
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_C2 = 10
F_BIG = 1 << 12
c2_table = {}
for label in ["12-bit/2557", "12-bit/2677"]:
    curve = curves[label]
    n = curve[2]
    print(f"\n{label}  (lam*={lam_star(curve[3], n):.4f})")
    print(f"{'K1':>4} {'V0 Kan+LLL':>11} {'Bab F=1':>8} {'enum F=1':>9} "
          f"{'Bab F=2^12':>11} {'enum F=2^12':>12} {'pl=closest':>11} "
          f"{'d_enum/d_pl':>12}")
    for k1 in K1_GRID:
        w0 = w1 = w2 = w1f = w2f = pic = 0
        ratios = []
        for seed in SEEDS:
            d0 = d_for_seed(seed, n)
            r0 = run_v0(curve, M_C2, d0, k1, seed=seed)
            if r0 and r0[0]:
                w0 += 1
            rc = run_cvp(curve, M_C2, d0, k1, seed=seed, F=1)
            if rc is not None:
                w1 += bool(rc['ok_babai'])
                if rc['ok_enum'] is not None:
                    w2 += bool(rc['ok_enum'])
            rf = run_cvp(curve, M_C2, d0, k1, seed=seed, F=F_BIG)
            if rf is not None:
                w1f += bool(rf['ok_babai'])
                if rf['ok_enum'] is not None:
                    w2f += bool(rf['ok_enum'])
                    pic += bool(rf['planted_is_closest'])
                    if rf['d_planted']:
                        ratios.append(rf['d_enum'] / rf['d_planted'])
        mr = sum(ratios) / len(ratios) if ratios else float('nan')
        c2_table[(label, k1)] = (w0, w1, w2, w1f, w2f, pic, mr)
        S = len(SEEDS)
        print(f"{k1:>4} {str(w0)+'/'+str(S):>11} {str(w1)+'/'+str(S):>8} "
              f"{str(w2)+'/'+str(S):>9} {str(w1f)+'/'+str(S):>11} "
              f"{str(w2f)+'/'+str(S):>12} {str(pic)+'/'+str(S):>11} "
              f"{mr:>12.4f}")

print("\n'pl=closest' (F=2^12): instances where the enumerated closest vector is")
print("NOT strictly closer than the planted error (d taken at its balanced")
print("representative).  0/5 => a closer lattice point exists => the instance")
print("is information-theoretically lost and no reduction can win.")

# ===========================================================================
# EXP C3 -- does inflating S_D rescue the Kannan/SVP form AT THE WALL?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C3: S_D sweep in the Kannan/SVP form (12-bit/2677, m=10)")
print("-" * 78)
print("2026-07-29 claimed 'no choice of S_D removes the trivial vector --")
print("both vectors scale linearly in S_D'.  Only the trivial vector does:")
print("||trivial|| = n*S_D but ||planted||^2 = C^2 + (d*S_D)^2 with C fixed,")
print("so the ratio tends to n/d > 1 and there must be a crossover.  Swept at")
print("K1 = 4 (baseline already wins) and K1 = 6, 8 (baseline fails).")
curve = curves["12-bit/2677"]
n = curve[2]
for k1 in [4, 6, 8]:
    print(f"\n  K1={k1}")
    print(f"  {'S_D':>10} {'wins':>7} {'mean sv/pv':>11} {'mean E_d':>9}")
    for S_D in [1, 2, 3, 4, 6, 8, 16, 32, 64, math.isqrt(n), n]:
        wins = 0
        rr, ee = [], []
        for seed in SEEDS:
            d0 = d_for_seed(seed, n)
            r = run_v0(curve, 10, d0, k1, seed=seed, S_D=S_D)
            if r is None:
                continue
            ok, pn, sn, e_d = r
            wins += bool(ok)
            rr.append(sn / pn)
            ee.append(e_d)
        print(f"  {S_D:>10} {str(wins)+'/'+str(len(SEEDS)):>7} "
              f"{sum(rr)/len(rr):>11.4f} {sum(ee)/len(ee):>9.3f}")

# ===========================================================================
# EXP C4 -- more data at the wall: does m rescue CVP where it did not rescue LLL?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C4: m sweep at the wall (12-bit/2677, K1=8), F=2^12")
print("-" * 78)
print("2026-07-29 T4b: at K1=8 the Kannan form got 0,0,1,0,1 of 5 for")
print("m = 8,12,16,24,32.  Does the CVP form do better?")
print(f"{'m':>4} {'dim':>5} {'V0 Kannan':>11} {'Babai':>7} {'BKZ20+Bab':>11} "
      f"{'enum':>6} {'pl=closest':>11}")
for m in [8, 10, 12, 14]:
    w0 = w1 = wb = w2 = pic = 0
    for seed in SEEDS:
        d0 = d_for_seed(seed, n)
        r0 = run_v0(curve, m, d0, 8, seed=seed)
        if r0 and r0[0]:
            w0 += 1
        rc = run_cvp(curve, m, d0, 8, seed=seed, F=F_BIG)
        if rc is not None:
            w1 += bool(rc['ok_babai'])
            if rc['ok_enum'] is not None:
                w2 += bool(rc['ok_enum'])
                pic += bool(rc['planted_is_closest'])
        rb = run_cvp(curve, m, d0, 8, seed=seed, F=F_BIG, do_enum=False,
                     bkz_beta=20)
        if rb is not None:
            wb += bool(rb['ok_babai'])
    S = len(SEEDS)
    print(f"{m:>4} {2*m+1:>5} {str(w0)+'/'+str(S):>11} {str(w1)+'/'+str(S):>7} "
          f"{str(wb)+'/'+str(S):>11} {str(w2)+'/'+str(S):>6} "
          f"{str(pic)+'/'+str(S):>11}")


# ===========================================================================
# EXP C5 -- what IS the closer vector?  relaxation artefact or real solution?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C5: decode the enumerated closest vector at the wall")
print("-" * 78)
print("If the closer-than-planted vector decodes to k1_i/k2_i OUTSIDE the")
print("true boxes [0,K1) x [0,K2), the loss is a RELAXATION artefact -- the")
print("lattice admits points the real problem forbids -- and tightening the")
print("encoding can recover them.  If it decodes INSIDE the boxes, a genuine")
print("second solution exists and the loss is information-theoretic.")
print(f"{'curve':<14} {'K1':>3} {'m':>3} {'seed':>6} {'d ok':>5} "
      f"{'d_en/d_pl':>10} {'k1 out':>7} {'k2 out':>7} {'k1 range':>16}")
for label, k1, m in [("12-bit/2677", 6, 10), ("12-bit/2677", 8, 10),
                     ("12-bit/2557", 12, 10)]:
    curve = curves[label]
    nn = curve[2]
    k2b = math.isqrt(nn) + 1
    for seed in SEEDS:
        d0 = d_for_seed(seed, nn)
        rc = run_cvp(curve, m, d0, k1, seed=seed, F=F_BIG, want_decode=True)
        if rc is None or rc.get('decode_enum') is None:
            continue
        k1s, k2s, d_dec = rc['decode_enum']
        b1, b2 = out_of_box(k1s, k2s, k1, k2b)
        print(f"{label:<14} {k1:>3} {m:>3} {seed:>6} "
              f"{str(bool(rc['ok_enum'])):>5} "
              f"{rc['d_enum']/rc['d_planted']:>10.4f} {b1:>7} {b2:>7} "
              f"{('['+str(min(k1s))+','+str(max(k1s))+']'):>16}")

# ===========================================================================
# EXP C6 -- box centering: the free factor-2 the prior construction missed
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C6: centered vs uncentered CVP target (m=10, 5 seeds, F=2^12)")
print("-" * 78)
print("Uncentered error carries k1_i in [0,K1); centered carries k1_i - K1/2.")
print("E[x^2] drops from K^2/3 to K^2/12, i.e. ||error|| halves.")
for label in ["12-bit/2557", "12-bit/2677"]:
    curve = curves[label]
    nn = curve[2]
    print(f"\n{label}  (lam*={lam_star(curve[3], nn):.4f})")
    print(f"{'K1':>4} {'V0 Kan+LLL':>11} {'unc Bab':>8} {'unc enum':>9} "
          f"{'ctr Bab':>8} {'ctr enum':>9} {'ctr pl=cl':>10} "
          f"{'||e||ctr/unc':>13}")
    for k1 in K1_GRID:
        w0 = wu = weu = wc = wec = pic = 0
        rr = []
        for seed in SEEDS:
            d0 = d_for_seed(seed, nn)
            r0 = run_v0(curve, 10, d0, k1, seed=seed)
            if r0 and r0[0]:
                w0 += 1
            ru = run_cvp(curve, 10, d0, k1, seed=seed, F=F_BIG)
            if ru is not None:
                wu += bool(ru['ok_babai'])
                if ru['ok_enum'] is not None:
                    weu += bool(ru['ok_enum'])
            rc = run_cvp(curve, 10, d0, k1, seed=seed, F=F_BIG, centered=True)
            if rc is not None:
                wc += bool(rc['ok_babai'])
                if rc['ok_enum'] is not None:
                    wec += bool(rc['ok_enum'])
                    pic += bool(rc['planted_is_closest'])
            if ru is not None and rc is not None and ru['d_planted']:
                rr.append(rc['d_planted'] / ru['d_planted'])
        S = len(SEEDS)
        mr = sum(rr) / len(rr) if rr else float('nan')
        print(f"{k1:>4} {str(w0)+'/'+str(S):>11} {str(wu)+'/'+str(S):>8} "
              f"{str(weu)+'/'+str(S):>9} {str(wc)+'/'+str(S):>8} "
              f"{str(wec)+'/'+str(S):>9} {str(pic)+'/'+str(S):>10} "
              f"{mr:>13.4f}")


# ===========================================================================
# EXP C7 -- ABLATION: is the gain from the CVP reformulation or from centering?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C7: ablation -- centering alone, inside the ORIGINAL Kannan/SVP form")
print("-" * 78)
print("2x2 design: {Kannan/SVP, CVP} x {uncentered, centered}, m=10, 5 seeds.")
print("If 'Kan ctr' matches 'CVP ctr', the operative ingredient is CENTERING")
print("and the CVP reformulation is cosmetic.")
for label in ["12-bit/2557", "12-bit/2677"]:
    curve = curves[label]
    nn = curve[2]
    print(f"\n{label}  (lam*={lam_star(curve[3], nn):.4f})")
    print(f"{'K1':>4} {'Kan unc':>8} {'Kan ctr':>8} {'Kan ctr+d':>10} "
          f"{'CVP unc Bab':>12} {'CVP ctr Bab':>12} {'CVP ctr enum':>13}")
    for k1 in K1_GRID:
        a = b_ = c = e = f = g = 0
        for seed in SEEDS:
            d0 = d_for_seed(seed, nn)
            r = run_v0(curve, 10, d0, k1, seed=seed)
            if r and r[0]: a += 1
            r = run_v0_centered(curve, 10, d0, k1, seed=seed)
            if r and r[0]: b_ += 1
            r = run_v0_centered(curve, 10, d0, k1, seed=seed, center_d=True)
            if r and r[0]: c += 1
            r = run_cvp(curve, 10, d0, k1, seed=seed, F=F_BIG, do_enum=False)
            if r: e += bool(r['ok_babai'])
            r = run_cvp(curve, 10, d0, k1, seed=seed, F=F_BIG, centered=True)
            if r:
                f += bool(r['ok_babai'])
                if r['ok_enum'] is not None:
                    g += bool(r['ok_enum'])
        S = len(SEEDS)
        print(f"{k1:>4} {str(a)+'/'+str(S):>8} {str(b_)+'/'+str(S):>8} "
              f"{str(c)+'/'+str(S):>10} {str(e)+'/'+str(S):>12} "
              f"{str(f)+'/'+str(S):>12} {str(g)+'/'+str(S):>13}")


# ===========================================================================
# EXP C8 -- generality: fresh curves, 10 seeds, wall location in eff = K1*K2/n
# ===========================================================================
print("\n" + "-" * 78)
print("EXP C8: generality of the centering fix -- fresh 17-bit curves, 10 seeds")
print("-" * 78)
print("Wall = largest K1 with >=8/10 recovery.  eff = K1*K2/n at that K1.")
SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 2718, 161803, 99991, 4242]
K1_WIDE = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

fresh = search_fresh_curves(2 ** 16, 2 ** 17, want=6)
allc = [("8-bit/199", curves["8-bit/199"]),
        ("12-bit/2557", curves["12-bit/2557"]),
        ("12-bit/2677", curves["12-bit/2677"])] + fresh

print(f"\n{'curve':<14} {'lam*':>7} {'n':>8} {'K2':>5} "
      f"{'K1 wall unc':>12} {'eff unc':>8} {'K1 wall ctr':>12} {'eff ctr':>8} "
      f"{'gain':>6}")
gains = []
for label, curve in allc:
    p, b, nn, lam, G = curve
    k2b = math.isqrt(nn) + 1
    wall_u = wall_c = 0
    for k1 in K1_WIDE:
        if k1 * k2b >= nn:
            break
        wu = wc = 0
        for seed in SEEDS10:
            d0 = d_for_seed(seed, nn)
            r = run_v0(curve, 10, d0, k1, seed=seed)
            if r and r[0]: wu += 1
            r = run_v0_centered(curve, 10, d0, k1, seed=seed)
            if r and r[0]: wc += 1
        if wu >= 8: wall_u = k1
        if wc >= 8: wall_c = k1
    eu = wall_u * k2b / nn if wall_u else float('nan')
    ec = wall_c * k2b / nn if wall_c else float('nan')
    g = (wall_c / wall_u) if wall_u else float('nan')
    if wall_u:
        gains.append(g)
    print(f"{label:<14} {lam_star(lam, nn):>7.4f} {nn:>8} {k2b:>5} "
          f"{wall_u:>12} {eu:>8.4f} {wall_c:>12} {ec:>8.4f} {g:>6.2f}")

if gains:
    print(f"\nmedian K1-wall gain from centering: "
          f"{sorted(gains)[len(gains)//2]:.2f}x over {len(gains)} curves")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
