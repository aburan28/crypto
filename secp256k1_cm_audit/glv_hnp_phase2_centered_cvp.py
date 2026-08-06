"""
GLV-HNP Phase 2, Thread 23: reformulate the Phase-2 lattice.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, T5):
  The Phase-2 Kannan lattice's shortest vector is always the trivial vector
  n*S_D*e_m, so the planted vector is never lambda_1.  That run proposed
  Thread 23: "project the lattice along e_m (quotient out the trivial n*e_m
  direction) and solve BDD in the projection, or replace the Kannan embedding
  with an explicit CVP call."  Falsifier stated there: if the K1 wall on the
  lam*=0.07 curve (currently K1 ~ 4-6) moves outward, the reformulation is a
  real improvement; if the wall stays, the wall is information-theoretic.

This script answers that, and finds a much larger defect along the way.

  DEFECT (H23a).  The lattice built by build_glv_lattice() (unchanged since
  2026-06-15) encodes the unknowns UNCENTERED: k1_i ~ U[0,K1), k2_i ~ U[0,K2),
  d ~ U[0,n).  Lattice reduction looks for vectors short around the ORIGIN, so
  the natural encoding is the centered one, k1_i - K1/2 etc.  Uncentered
  E[x^2] = X^2/3; centered E[x^2] = X^2/12.  So

      ||e_centered|| / ||e_uncentered||  ->  1/2

  exactly, for every coordinate simultaneously.  A factor 2 in BDD distance is
  a factor 4 in the bias strength eff = K1*K2/n at which the attack survives.
  Centering is free: it is a translation of the Kannan/CVP target, not a
  change of lattice.

  Centering algebra.  With c1 = K1//2, c2 = K2//2, cd = (n-1)//2 and
      Atil_i = A_i + B_i*cd - c1 - lam*c2   (mod n)
  the congruence  k1_i = A_i + B_i*d - lam*k2_i (mod n)  becomes
      k1_i - c1 = Atil_i + B_i*(d - cd) - lam*(k2_i - c2)   (mod n)
  so the SAME lattice with A replaced by Atil has planted vector
      ((k1_i - c1)*S_K1, (d - cd)*S_D, (k2_i - c2)*S_K2, S_KANNAN).

  REFORMULATION (H23b).  Deleting column m (orthogonal projection along e_m)
  gives the rank-2m integral lattice pi(L) = <rows 0..2m, column m deleted>,
  with det(pi(L)) = det(L)/(n*S_D).  The d-ambiguity d -> d+n lives exactly in
  the killed direction, so nothing is lost: after CVP in pi(L) recovers
  (k1'_i, k2'_i), d follows algebraically from
      d' = (k1'_i - Atil_i + lam*k2'_i) * B_i^{-1}  (mod n),   d = d' + cd.
  This removes the trivial vector AND the (d*S_D)^2 term from the error norm.

Formulations compared (identical signature instances throughout):
  F0  Kannan / SVP, uncentered   -- the 2026-06-15 baseline, verbatim
  F1  Kannan / SVP, centered
  F2  projected CVP (pi along e_m) + Babai nearest plane, centered
  F3  F2 with BKZ-20 preprocessing

Run: python3 glv_hnp_phase2_centered_cvp.py
"""

import math
import random
from fractions import Fraction

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic  -- verbatim from glv_hnp_phase2_20bit.py
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

# ---------------------------------------------------------------------------
# CM theory for j=0 curves -- verbatim from glv_hnp_phase2_lambda_threshold.py
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

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

# ---------------------------------------------------------------------------
# Signatures + scales  -- verbatim from glv_hnp_phase2_20bit.py
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

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# F0: the 2026-06-15 Kannan lattice, uncentered (verbatim)
# ---------------------------------------------------------------------------

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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret, shift=0):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m] + shift) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# F1: the SAME lattice with a centered target  (H23a)
# ---------------------------------------------------------------------------

def centering(n, lam, k1_bound, k2_bound):
    """(c1, c2, cd): the centers subtracted from k1, k2, d."""
    return (k1_bound // 2, k2_bound // 2, (n - 1) // 2)

def centered_A(sigs, n, lam, k1_bound, k2_bound):
    """Atil_i = A_i + B_i*cd - c1 - lam*c2  (mod n)."""
    c1, c2, cd = centering(n, lam, k1_bound, k2_bound)
    return [(s['A'] + s['B'] * cd - c1 - lam * c2) % n for s in sigs]

def build_glv_lattice_centered(sigs, n, lam, k1_bound, k2_bound):
    """Identical row space to build_glv_lattice(); only the Kannan row moves."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    At = centered_A(sigs, n, lam, k1_bound, k2_bound)
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
        M[2 * m + 1][i] = At[i] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector_centered(sigs, d_secret, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1, c2, cd = centering(n, lam, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - c2) * S_K2
    v[m] = (d_secret - cd) * S_D
    v[2 * m + 1] = S_KAN
    return v

# ---------------------------------------------------------------------------
# F2/F3: projected lattice pi(L) (column m deleted) + Babai nearest plane
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Rows 0..2m of the Phase-2 lattice with column m deleted.
    (2m+1) generators of a rank-2m lattice in Z^{2m}."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                                   # n*S_K1*e_i
        r = [0] * (2 * m); r[i] = n * S_K1; rows.append(r)
    r = [0] * (2 * m)                                    # pi(d-row)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                   # k2 rows
        r = [0] * (2 * m); r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    return rows

def projected_target(sigs, n, lam, k1_bound, k2_bound, centered=True):
    """t = (-Atil_i * S_K1, 0, ..., 0); lattice point - t = the error."""
    m = len(sigs)
    S_K1, _, _, _ = scales(n, k1_bound, k2_bound)
    At = (centered_A(sigs, n, lam, k1_bound, k2_bound) if centered
          else [s['A'] % n for s in sigs])
    return [-At[i] * S_K1 for i in range(m)] + [0] * m

def gram_schmidt(B):
    Bs = []
    for b in B:
        v = [Fraction(x) for x in b]
        for bs in Bs:
            nj = sum(x * x for x in bs)
            if nj == 0:
                continue
            mu = sum(Fraction(a) * c for a, c in zip(b, bs)) / nj
            v = [x - mu * y for x, y in zip(v, bs)]
        Bs.append(v)
    return Bs

def babai_nearest_plane(B, Bs, t):
    """Exact nearest-plane. Returns the lattice vector w closest-ish to t."""
    b = [Fraction(x) for x in t]
    w = [0] * len(t)
    for i in reversed(range(len(B))):
        nj = sum(x * x for x in Bs[i])
        if nj == 0:
            continue
        c = sum(x * y for x, y in zip(b, Bs[i])) / nj
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        if ci:
            b = [x - ci * y for x, y in zip(b, B[i])]
            w = [x + ci * y for x, y in zip(w, B[i])]
    return w

def reduce_projected(rows, dim_cols, use_bkz=False, bkz_beta=20):
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    out = [[A[i][j] for j in range(dim_cols)] for i in range(A.nrows)]
    return [r for r in out if any(r)]

def cvp_recover_d(sigs, n, lam, k1_bound, k2_bound, use_bkz=False, bkz_beta=20):
    """F2/F3: solve BDD in pi(L), then read d off algebraically."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    c1, c2, cd = centering(n, lam, k1_bound, k2_bound)
    At = centered_A(sigs, n, lam, k1_bound, k2_bound)
    rows = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    t = projected_target(sigs, n, lam, k1_bound, k2_bound, centered=True)
    B = reduce_projected(rows, 2 * m, use_bkz=use_bkz, bkz_beta=bkz_beta)
    Bs = gram_schmidt(B)
    w = babai_nearest_plane(B, Bs, t)
    e = [w[i] - t[i] for i in range(2 * m)]
    if any(e[i] % S_K1 for i in range(m)):
        return None, e
    k1p = [e[i] // S_K1 for i in range(m)]
    k2p = [e[m + i] // S_K2 for i in range(m)]
    # d' = (k1'_i - Atil_i + lam*k2'_i) * B_i^{-1} mod n, majority over i
    votes = {}
    for i in range(m):
        try:
            dp = (k1p[i] - At[i] + lam * k2p[i]) * modinv(sigs[i]['B'], n) % n
        except ValueError:
            continue
        votes[dp] = votes.get(dp, 0) + 1
    if not votes:
        return None, e
    dp = max(votes, key=votes.get)
    return (dp + cd) % n, e

# ---------------------------------------------------------------------------
# Single-instance drivers
# ---------------------------------------------------------------------------

def make_instance(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    return (sigs if len(sigs) == m else None), k2_bound

def run_kannan(curve, m, d_secret, k1_bound, seed, centered, sigs=None,
               k2_bound=None, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    if sigs is None:
        sigs, k2_bound = make_instance(curve, m, d_secret, k1_bound, seed)
    if sigs is None:
        return None
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    if centered:
        M = build_glv_lattice_centered(sigs, n, lam, k1_bound, k2_bound)
        _, _, cd = centering(n, lam, k1_bound, k2_bound)
        shift = cd
        pv = planted_vector_centered(sigs, d_secret, n, lam, k1_bound, k2_bound)
    else:
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        shift = 0
        pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_d(reduced, m, n, S_KAN, d_secret, shift=shift) is not None
    pn = norm(pv)
    # norm without the constant Kannan coordinate: the informative part
    pn_core = norm(pv[:dim - 1])
    return ok, pn, pn_core, min(norm(r) for r in reduced) / pn

def run_cvp(curve, m, d_secret, k1_bound, seed, sigs=None, k2_bound=None,
            use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    if sigs is None:
        sigs, k2_bound = make_instance(curve, m, d_secret, k1_bound, seed)
    if sigs is None:
        return None
    d_cand, e = cvp_recover_d(sigs, n, lam, k1_bound, k2_bound,
                              use_bkz=use_bkz, bkz_beta=bkz_beta)
    return (d_cand == d_secret), norm(e)

# ---------------------------------------------------------------------------
# Predicted BDD ratio (Gaussian heuristic)
# ---------------------------------------------------------------------------

def gh_ratio(m, n, k1_bound, k2_bound, centered, projected):
    """r = E||e|| / lambda_1^GH.  r < 1 is the classical feasibility line."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    var = 1.0 / 12.0 if centered else 1.0 / 3.0
    if projected:                       # pi(L): dim 2m, det = det(L)/(n*S_D)
        N = 2 * m
        logdet = m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_D) \
            - math.log(n * S_D)
        e2 = m * var * (k1_bound * S_K1) ** 2 + m * var * (k2_bound * S_K2) ** 2
    else:                               # core lattice L: dim 2m+1
        N = 2 * m + 1
        logdet = m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_D)
        e2 = m * var * (k1_bound * S_K1) ** 2 + m * var * (k2_bound * S_K2) ** 2 \
            + var * (n * S_D) ** 2
    lam1 = math.sqrt(N / (2 * math.pi * math.e)) * math.exp(logdet / N)
    return math.sqrt(e2) / lam1

# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  m
    ("8-bit/199",        211,  2, 199,  106,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8),
    ("12-bit/2677",      2677, 2, 2647, 185,  10),
]

hist = []
for label, p, b, n, lam, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), m))

print("=" * 78)
print("Thread 23 - reformulating the GLV-HNP Phase-2 lattice")
print("=" * 78)
print("\nCurves (lam* = min(lam,n-lam)/n):")
for label, curve, m in hist:
    p, b, n, lam, G = curve
    print(f"  {label:<14} p={p:<6} n={n:<6} lam={lam:<6} "
          f"lam*={lam_star(lam, n):.4f}  K2={math.isqrt(n)+1:<4} m={m}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C0: correctness of the centered / projected formulations")
print("-" * 78)
print("Note: centering is a TRANSLATION OF THE TARGET, so the Kannan embedding")
print("lattice does change (the Kannan row moves by a non-lattice vector).")
print("What must hold: (a) the core sublattice (rows 0..2m) is untouched,")
print("(b) the centered planted vector is an explicit integer combination of")
print("the centered rows, (c) the projected error decodes back to the true")
print("(k1_i, k2_i, d).\n")

def lattice_combination_ok(rows, target, coeffs):
    acc = [0] * len(target)
    for c, r in zip(coeffs, rows):
        if c:
            acc = [a + c * x for a, x in zip(acc, r)]
    return acc == list(target)

for label, curve, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    K1t = 4
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sg = gen_signatures(G, d0, m, n, lam, p, K1t, k2b, 42)
    M0 = build_glv_lattice(sg, n, lam, K1t, k2b)
    M1 = build_glv_lattice_centered(sg, n, lam, K1t, k2b)
    core_same = M0[:2 * m + 1] == M1[:2 * m + 1]

    # (b) explicit membership of the centered planted vector
    c1, c2, cd = centering(n, lam, K1t, k2b)
    S_K1, S_D, S_K2, _ = scales(n, K1t, k2b)
    At = centered_A(sg, n, lam, K1t, k2b)
    dp = d0 - cd
    coeffs = [0] * (2 * m + 2)
    exact = True
    for i in range(m):
        k2p = sg[i]['k2'] - c2
        coeffs[m + 1 + i] = k2p
        # q_i so that column i lands on (k1_i - c1)*S_K1
        num = (sg[i]['k1'] - c1) - (At[i] + dp * sg[i]['B'] - lam * k2p)
        exact &= (num % n == 0)
        coeffs[i] = num // n
    coeffs[m] = dp
    coeffs[2 * m + 1] = 1
    pv = planted_vector_centered(sg, d0, n, lam, K1t, k2b)
    memb = exact and lattice_combination_ok(M1, pv, coeffs)

    # (c) projected CVP decodes the planted point when handed to it exactly
    d_cand, _ = cvp_recover_d(sg, n, lam, K1t, k2b)
    print(f"  {label:<14} core sublattice unchanged: {str(core_same):<6} "
          f"planted vector in centered lattice: {str(memb):<6} "
          f"F2 d-recovery: {d_cand == d0}")

# no-false-positive check: recovery must fail on a wrong secret
label, curve, m = hist[1]
p, b, n, lam, G = curve
k2b = math.isqrt(n) + 1
d0 = random.Random(1)  # placeholder
d0 = 12345 % n or 7
sg = gen_signatures(G, d0, m, n, lam, p, 4, k2b, 42)
fp = recover_d([[A[j] for j in range(2 * m + 2)]
                for A in build_glv_lattice_centered(sg, n, lam, 4, k2b)],
               m, n, n, (d0 + 1) % n, shift=(n - 1) // 2)
print(f"  no-false-positive on a wrong secret: {fp is None}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C1: the K1 wall - F0 (uncentered) vs F1 (centered) vs F2 (proj CVP)")
print("-" * 78)
print("Prediction H23a: centering halves ||e||, so the K1 wall should move")
print("outward by ~4x in eff = K1*K2/n.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]

c1_rows = []
for label, curve, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"{label}  (m={m}, K2={k2b}, n={n})")
    print(f"  {'K1':>5} {'eff':>7} {'F0':>6} {'F1':>6} {'F2':>6} {'F3':>6} "
          f"{'|e0|/n':>8} {'|e1|/n':>8} {'|e2|/n':>8} {'sv/pv0':>7} "
          f"{'sv/pv1':>7} {'r_GH1':>7} {'r_GHp':>7}")
    for K1 in K1_GRID:
        if K1 >= n:
            continue
        eff = K1 * k2b / n
        w0 = w1 = w2 = w3 = 0
        e0s, e1s, e2s, sv0s, sv1s = [], [], [], [], []
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            sigs, _ = make_instance(curve, m, d_trial, K1, seed)
            if sigs is None:
                continue
            r0 = run_kannan(curve, m, d_trial, K1, seed, False, sigs, k2b)
            r1 = run_kannan(curve, m, d_trial, K1, seed, True, sigs, k2b)
            r2 = run_cvp(curve, m, d_trial, K1, seed, sigs, k2b)
            r3 = run_cvp(curve, m, d_trial, K1, seed, sigs, k2b, use_bkz=True)
            w0 += bool(r0[0]); e0s.append(r0[2]); sv0s.append(r0[3])
            w1 += bool(r1[0]); e1s.append(r1[2]); sv1s.append(r1[3])
            w2 += bool(r2[0]); e2s.append(r2[1])
            w3 += bool(r3[0])
        ns = len(e0s)
        if ns == 0:
            continue
        e0 = sum(e0s) / ns / n; e1 = sum(e1s) / ns / n; e2 = sum(e2s) / ns / n
        sv0 = sum(sv0s) / ns; sv1 = sum(sv1s) / ns
        rg1 = gh_ratio(m, n, K1, k2b, centered=True, projected=False)
        rgp = gh_ratio(m, n, K1, k2b, centered=True, projected=True)
        print(f"  {K1:>5} {eff:>7.3f} {str(w0)+'/'+str(ns):>6} "
              f"{str(w1)+'/'+str(ns):>6} {str(w2)+'/'+str(ns):>6} "
              f"{str(w3)+'/'+str(ns):>6} "
              f"{e0:>8.2f} {e1:>8.2f} {e2:>8.2f} {sv0:>7.3f} {sv1:>7.3f} "
              f"{rg1:>7.3f} {rgp:>7.3f}")
        c1_rows.append(dict(label=label, m=m, n=n, K1=K1, K2=k2b, eff=eff,
                            ns=ns, w0=w0, w1=w1, w2=w2, w3=w3,
                            e0=e0, e1=e1, e2=e2, sv0=sv0, sv1=sv1,
                            rg1=rg1, rgp=rgp))
    print()
print("|e*|/n columns exclude the constant Kannan coordinate (the informative")
print("part of the BDD distance).  sv/pv = shortest reduced row / planted norm.")

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP C2: measured ||e_uncentered|| / ||e_centered||")
print("-" * 78)
print("Continuous prediction: E[x^2] = X^2/3 uncentered vs X^2/12 centered, so")
print("the ratio -> 2.00 as K1 grows (for small K1 the discrete distribution")
print("over {0..K1-1} departs from the uniform-interval approximation).")
print("|e1|/|e2| predicted sqrt((2m+1)/(2m)): the projection drops the d-term.\n")
print(f"  {'curve':<14} {'K1':>5} {'|e0|/|e1|':>10} {'|e1|/|e2|':>10} "
      f"{'pred e1/e2':>11}")
for r in c1_rows:
    if r['K1'] in (4, 16, 64, 128):
        pred = math.sqrt((2 * r['m'] + 1) / (2 * r['m']))
        print(f"  {r['label']:<14} {r['K1']:>5} {r['e0']/r['e1']:>10.3f} "
              f"{r['e1']/r['e2']:>10.3f} {pred:>11.3f}")

print("\n  C2b: the |e0|/|e1| shortfall above is small-sample noise (5 seeds).")
print("  Monte Carlo over the same (m, n, K1, K2) with 4000 draws:")
print(f"  {'curve':<14} {'K1':>5} {'MC |e0|/n':>10} {'MC |e1|/n':>10} {'ratio':>7}")
for label, curve, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for K1 in (16, 128):
        S_K1, S_D, S_K2, _ = scales(n, K1, k2b)
        c1, c2, cd = centering(n, lam, K1, k2b)
        rng = random.Random(20260806)
        su = sc = 0.0
        T = 4000
        for _ in range(T):
            kk1 = [rng.randrange(K1) for _ in range(m)]
            kk2 = [rng.randrange(k2b) for _ in range(m)]
            dd = rng.randrange(1, n)
            su += math.sqrt(sum((x * S_K1) ** 2 for x in kk1)
                            + sum((x * S_K2) ** 2 for x in kk2) + dd * dd)
            sc += math.sqrt(sum(((x - c1) * S_K1) ** 2 for x in kk1)
                            + sum(((x - c2) * S_K2) ** 2 for x in kk2)
                            + (dd - cd) ** 2)
        su /= T * n; sc /= T * n
        print(f"  {label:<14} {K1:>5} {su:>10.3f} {sc:>10.3f} {su/sc:>7.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C3: where is the wall?  (last K1 with >=4/5, per formulation)")
print("-" * 78)
print(f"  {'curve':<14} {'F0 K1*':>7} {'F1 K1*':>7} {'F2 K1*':>7} {'F3 K1*':>7} "
      f"{'F1/F0':>7} {'eff*(F0)':>9} {'eff*(F1)':>9}")
for label, curve, m in hist:
    rows = [r for r in c1_rows if r['label'] == label]
    def wall(key):
        best = 0
        for r in rows:
            if r[key] >= 4:
                best = max(best, r['K1'])
        return best
    w0, w1, w2, w3 = wall('w0'), wall('w1'), wall('w2'), wall('w3')
    k2b = rows[0]['K2']; n = rows[0]['n']
    print(f"  {label:<14} {w0:>7} {w1:>7} {w2:>7} {w3:>7} "
          f"{(w1/w0 if w0 else float('nan')):>7.2f} "
          f"{w0*k2b/n:>9.3f} {w1*k2b/n:>9.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C4: is r_GH a predictor?  (pooled over all curves and K1)")
print("-" * 78)

def threshold_report(rows, key_r, key_w, name):
    pts = [(r[key_r], r[key_w] >= 4) for r in rows]
    pts.sort()
    best, bt = -1, None
    cand = sorted(set(p[0] for p in pts))
    for t in cand:
        acc = sum(1 for r, s in pts if (r < t) == s) / len(pts)
        if acc > best:
            best, bt = acc, t
    pos = [r for r, s in pts if s]
    neg = [r for r, s in pts if not s]
    auc = float('nan')
    if pos and neg:
        auc = sum((1.0 if a < b else 0.5 if a == b else 0.0)
                  for a in pos for b in neg) / (len(pos) * len(neg))
    base = max(sum(1 for _, s in pts if s), sum(1 for _, s in pts if not s)) / len(pts)
    print(f"  {name:<28} best-threshold r*={bt:.3f}  acc={best:.3f}  "
          f"(baseline {base:.3f})  AUC={auc:.3f}")

threshold_report(c1_rows, 'rg1', 'w1', 'r_GH (centered Kannan) -> F1')
threshold_report(c1_rows, 'rgp', 'w2', 'r_GH (projected)       -> F2')
threshold_report(c1_rows, 'rg1', 'w0', 'r_GH (centered Kannan) -> F0')

print("\n  Closed form.  With eff = K1*K2/n and N = 2m+1,")
print("      r_GH = sqrt(2*pi*e/12) * eff^(m/N) * n^(1/N)   ~= 1.194 * ...")
print("  m -> infinity gives r_GH -> 1.194*sqrt(eff), i.e. a ceiling of")
print("  eff <= 0.70; the n^(1/N) factor is the finite-m penalty.")
print(f"  {'curve':<14} {'m':>3} {'K1':>5} {'closed form':>12} {'gh_ratio()':>11}")
for r in c1_rows:
    if r['K1'] in (4, 16, 64):
        N = 2 * r['m'] + 1
        cf = (math.sqrt(2 * math.pi * math.e / 12)
              * r['eff'] ** (r['m'] / N) * r['n'] ** (1.0 / N))
        print(f"  {r['label']:<14} {r['m']:>3} {r['K1']:>5} {cf:>12.3f} "
              f"{r['rg1']:>11.3f}")
print("\n  Asymptotic eff ceilings (r_GH = 1):")
print(f"    centered   eff* = {1/(2*math.pi*math.e/12):.3f}")
print(f"    uncentered eff* = {1/(2*math.pi*math.e/3):.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C5: fresh 17-bit curves - does the centering gain hold out of sample?")
print("-" * 78)

def find_17bit(count=6, lo=2**16, hi=2**17):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < count:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    out.append((p, cur[1], nc, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

curves17 = find_17bit()
print(f"  found {len(curves17)} 17-bit j=0 GLV curves\n")
M17 = 12
tot = {'F0': 0, 'F1': 0, 'F2': 0, 'N': 0}
per_eff = {}
print(f"  {'p':>7} {'n':>7} {'lam*':>7} {'eff':>6} {'r_GH1':>6} "
      f"{'F0':>5} {'F1':>5} {'F2':>5}")
for (p, b, n, lam, G) in curves17:
    k2b = math.isqrt(n) + 1
    for eff in (0.15, 0.25, 0.35, 0.50):
        K1 = max(2, int(eff * n / k2b))
        curve = (p, b, n, lam, G)
        w0 = w1 = w2 = 0; ns = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            sigs, _ = make_instance(curve, M17, d_trial, K1, seed)
            if sigs is None:
                continue
            ns += 1
            w0 += bool(run_kannan(curve, M17, d_trial, K1, seed, False, sigs, k2b)[0])
            w1 += bool(run_kannan(curve, M17, d_trial, K1, seed, True, sigs, k2b)[0])
            w2 += bool(run_cvp(curve, M17, d_trial, K1, seed, sigs, k2b)[0])
        if ns:
            rg = gh_ratio(M17, n, K1, k2b, centered=True, projected=False)
            print(f"  {p:>7} {n:>7} {lam_star(lam,n):>7.4f} {K1*k2b/n:>6.2f} "
                  f"{rg:>6.3f} "
                  f"{str(w0)+'/'+str(ns):>5} {str(w1)+'/'+str(ns):>5} "
                  f"{str(w2)+'/'+str(ns):>5}")
            tot['F0'] += w0; tot['F1'] += w1; tot['F2'] += w2; tot['N'] += ns
            a = per_eff.setdefault(eff, [0, 0, 0, 0])
            a[0] += w0; a[1] += w1; a[2] += w2; a[3] += ns
    print()

print(f"  pooled 17-bit: F0 {tot['F0']}/{tot['N']}  F1 {tot['F1']}/{tot['N']}  "
      f"F2 {tot['F2']}/{tot['N']}")
print(f"  {'eff':>6} {'F0':>8} {'F1':>8} {'F2':>8}")
for eff in sorted(per_eff):
    a = per_eff[eff]
    print(f"  {eff:>6.2f} {str(a[0])+'/'+str(a[3]):>8} "
          f"{str(a[1])+'/'+str(a[3]):>8} {str(a[2])+'/'+str(a[3]):>8}")

print("=" * 78)
print("done")
print("=" * 78)
