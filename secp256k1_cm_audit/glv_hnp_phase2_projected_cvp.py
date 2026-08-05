"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is
the target of a *unique-BDD / CVP* instance rather than an SVP instance.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, T5):
  In the Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` the vector
  `n * S_D * e_m` is ALWAYS shorter than the planted vector (norm n vs
  n*sqrt(2m/3 + 4/3)), so lambda_1 is never the answer.  It is the trivial
  "d is only defined mod n" vector and carries no information.

Reformulation:
  Drop the d-coordinate (S_D -> 0, i.e. project along e_m) and drop the Kannan
  row; solve CVP directly.  d is recovered *afterwards* from the k1/k2 values.
  In the projected lattice the trivial direction is quotiented out, so the
  n*e_m competitor disappears entirely (in CVP form it would anyway only be a
  benign d -> d+n shift of the same answer).

  L_proj  subset Z^{2m},  coords = (k1-cols 0..m-1 | k2-cols m..2m-1)
    row 0        : (B'_i * S_K1)_i        (+) 0          , B'_i = B_i/B_0 mod n
    rows 1..m-1  : n * S_K1 * e_i
    rows m..2m-1 : (-lam * S_K1 * e_i)    (+) S_K2 * e_i
  target  t = (A_i * S_K1)_i (+) 0 ;   t + u = (k1_i*S_K1 | k2_i*S_K2).
  vol(L_proj) = S_K1^m * S_K2^m * n^(m-1).

Three methods are compared on identical signature sets:
  KAN-ORIG : original 2m+2 Kannan/SVP lattice (baseline, verbatim)
  PROJ-CVP : L_proj + exact CVP (fpylll) / Babai fallback     <- the proposal
  PROJ-KAN : Kannan embedding of L_proj (dim 2m+1, SVP form)  <- to measure sv/pv

Falsifier stated on 2026-07-29:
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_projected_cvp.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ, GSO

try:
    from fpylll import CVP as _CVP
except ImportError:  # pragma: no cover
    _CVP = None

# ---------------------------------------------------------------------------
# EC + CM helpers (verbatim from glv_hnp_phase2_lambda_threshold.py)
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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

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
# Signatures (verbatim)
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

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) — verbatim from the Phase-2 construction."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def gauss_reduce_2d(u, v):
    """Exact shortest nonzero vector of a 2-D integer lattice (Lagrange).
    Verbatim from glv_hnp_phase2_lambda_threshold.py:188."""
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
    """|shortest vector| of the 2-D block <(n*S_K1, 0), (-lam*S_K1, S_K2)>.

    This block sits inside L_proj once per signature, so its shortest vector is
    the residual SVP competitor after the trivial n*e_m direction is quotiented
    out.  mu ~ sqrt(n*S_K1*S_K2) = n/sqrt(eff).
    """
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w

# ---------------------------------------------------------------------------
# METHOD 1 — KAN-ORIG: original 2m+2 Kannan/SVP lattice (baseline)
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

def planted_vector_orig(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_orig(rows, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def run_kan_orig(sigs, n, lam, k1_bound, k2_bound, d_secret, beta=0):
    m = len(sigs)
    dim = 2 * m + 2
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ok = recover_d_orig(rows, m, n, S_KAN, d_secret) is not None
    pv = norm(planted_vector_orig(sigs, d_secret, n, k1_bound, k2_bound))
    sv = min((norm(r) for r in rows if any(r)), default=float('inf'))
    return ok, sv, pv

# ---------------------------------------------------------------------------
# METHOD 2/3 — the projected lattice L_proj
# ---------------------------------------------------------------------------

def build_proj_basis(sigs, n, lam, k1_bound, k2_bound):
    """Exact basis of L_proj (rank 2m, no dependent generators) + target t.

    Uses B'_i = B_i * B_0^{-1} mod n so that B'_0 = 1; the free parameter of
    the d-row then *is* the first k1-column residue, which lets us drop the
    now-redundant n*S_K1*e_0 row and keep a genuine basis.
    """
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    Bp = [sigs[i]['B'] * B0inv % n for i in range(m)]
    assert Bp[0] == 1

    dim = 2 * m
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[0][i] = Bp[i] * S_K1
    for i in range(1, m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m + i][i] = -(lam % n) * S_K1
        M[m + i][m + i] = S_K2
    t = [sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    return M, t, S_K1, S_K2

def planted_proj(sigs, n, k1_bound, k2_bound):
    """(k1_i*S_K1 | k2_i*S_K2) — the vector CVP must land on."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return [sigs[i]['k1'] * S_K1 for i in range(m)] + \
           [sigs[i]['k2'] * S_K2 for i in range(m)]

def check_proj_membership(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Sanity: t + u_planted must equal planted_proj for some u in L_proj."""
    M, t, S_K1, S_K2 = build_proj_basis(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    target = planted_proj(sigs, n, k1_bound, k2_bound)
    diff = [target[j] - t[j] for j in range(2 * m)]
    # solve for integer coefficients: k2 rows fix coords m..2m-1, then row 0
    # coefficient is read off coord 0, and rows 1..m-1 absorb the rest.
    coeff = [0] * (2 * m)
    for i in range(m):
        assert diff[m + i] % S_K2 == 0
        coeff[m + i] = diff[m + i] // S_K2
    rem = [diff[i] + (lam % n) * S_K1 * coeff[m + i] for i in range(m)]
    assert rem[0] % S_K1 == 0
    c0 = rem[0] // S_K1
    coeff[0] = c0
    B0inv = modinv(sigs[0]['B'] % n, n)
    Bp = [sigs[i]['B'] * B0inv % n for i in range(m)]
    for i in range(1, m):
        r = rem[i] - c0 * Bp[i] * S_K1
        assert r % (n * S_K1) == 0, "planted vector is NOT in L_proj"
        coeff[i] = r // (n * S_K1)
    return True

def babai_nearest_plane_exact(Ared, target):
    """Proper Babai nearest-plane using exact-ish float GSO computed here."""
    d = Ared.nrows
    nc = Ared.ncols
    B = [[float(Ared[i][k]) for k in range(nc)] for i in range(d)]
    # Gram-Schmidt
    Bs = []
    mu = [[0.0] * d for _ in range(d)]
    for i in range(d):
        v = B[i][:]
        for j in range(len(Bs)):
            dj = sum(Bs[j][k] * Bs[j][k] for k in range(nc))
            if dj == 0:
                mu[i][j] = 0.0
                continue
            mu[i][j] = sum(B[i][k] * Bs[j][k] for k in range(nc)) / dj
            for k in range(nc):
                v[k] -= mu[i][j] * Bs[j][k]
        Bs.append(v)
    w = [float(x) for x in target]
    coeffs = [0] * d
    for i in range(d - 1, -1, -1):
        dj = sum(Bs[i][k] * Bs[i][k] for k in range(nc))
        if dj == 0:
            continue
        c = round(sum(w[k] * Bs[i][k] for k in range(nc)) / dj)
        coeffs[i] = int(c)
        for k in range(nc):
            w[k] -= c * B[i][k]
    out = [0] * nc
    for i in range(d):
        if coeffs[i]:
            for k in range(nc):
                out[k] += coeffs[i] * int(Ared[i][k])
    return out

def _decode_proj(vec, sigs, n, lam, k1_bound, k2_bound):
    """From a candidate (k1_i*S_K1 | k2_i*S_K2), recover d WITHOUT the oracle.

    Success criterion is self-checking: all k1_i in [0,K1), all k2_i in
    [0,K2), and the derived d must be consistent across every signature.
    """
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    k1s, k2s = [], []
    for i in range(m):
        a, b = vec[i], vec[m + i]
        if a % S_K1 or b % S_K2:
            return None
        k1, k2 = a // S_K1, b // S_K2
        if not (0 <= k1 < k1_bound) or not (0 <= k2 < k2_bound):
            return None
        k1s.append(k1); k2s.append(k2)
    d = None
    for i in range(m):
        k_full = (k1s[i] + lam * k2s[i]) % n
        Bi = sigs[i]['B'] % n
        if Bi == 0:
            continue
        di = (k_full - sigs[i]['A']) * modinv(Bi, n) % n
        if d is None:
            d = di
        elif di != d:
            return None
    return d

def centering_shift(m, k1_bound, k2_bound, S_K1, S_K2):
    """Offset that recentres the solution box [0,K)^2m on the origin.

    CVP/SVP minimise a norm centred at 0, but the planted k1_i, k2_i are
    uniform on [0,K), i.e. ALL-POSITIVE.  Uncentred,
    E||pv||^2 = m(K1*S_K1)^2/3 + m(K2*S_K2)^2/3; centred it is /12 — a factor
    2 in norm, which is exactly the handicap that makes the uncentred CVP
    degrade as m grows.
    """
    return [(k1_bound // 2) * S_K1] * m + [(k2_bound // 2) * S_K2] * m

def run_proj_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret, beta=0,
                 exact_cvp=True, center=False):
    """METHOD 2 — projected lattice, solved as CVP."""
    m = len(sigs)
    M, t, S_K1, S_K2 = build_proj_basis(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    sh = centering_shift(m, k1_bound, k2_bound, S_K1, S_K2) if center \
        else [0] * (2 * m)
    tc = [t[j] - sh[j] for j in range(2 * m)]
    neg_t = [-x for x in tc]
    u = None
    if exact_cvp and _CVP is not None:
        try:
            u = list(_CVP.closest_vector(A, tuple(neg_t)))
        except Exception:
            u = None
    if u is None:
        u = babai_nearest_plane_exact(A, neg_t)
    cand_c = [tc[j] + u[j] for j in range(2 * m)]
    cand = [cand_c[j] + sh[j] for j in range(2 * m)]
    d = _decode_proj(cand, sigs, n, lam, k1_bound, k2_bound)
    pv = planted_proj(sigs, n, k1_bound, k2_bound)
    pvc = [pv[j] - sh[j] for j in range(2 * m)]
    # ||t_c + u_cvp||, the distance from the (centred) target to the returned
    # lattice point.  Exact CVP => this is <= ||t_c + u_planted|| = ||pv_c||;
    # a value strictly below 1.0 together with a decode failure means the
    # planted vector genuinely is NOT the closest vector (a real BDD failure,
    # not a reduction failure).
    return (d is not None and d == d_secret), norm(cand), norm(pvc), norm(cand_c)

def run_proj_kan(sigs, n, lam, k1_bound, k2_bound, d_secret, beta=0,
                 center=False):
    """METHOD 3 — Kannan embedding of L_proj (SVP form), to measure sv/pv."""
    m = len(sigs)
    M, t, S_K1, S_K2 = build_proj_basis(sigs, n, lam, k1_bound, k2_bound)
    S_KAN = n
    dim = 2 * m + 1
    sh = centering_shift(m, k1_bound, k2_bound, S_K1, S_K2) if center \
        else [0] * (2 * m)
    tc = [t[j] - sh[j] for j in range(2 * m)]
    MM = [row + [0] for row in M]
    MM.append(tc + [S_KAN])
    A = IntegerMatrix.from_matrix(MM)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    pv = [planted_proj(sigs, n, k1_bound, k2_bound)[j] - sh[j]
          for j in range(2 * m)] + [S_KAN]
    ok = False
    for row in rows:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        cand = [sign * row[j] + sh[j] for j in range(2 * m)]
        d = _decode_proj(cand, sigs, n, lam, k1_bound, k2_bound)
        if d is not None and d == d_secret:
            ok = True
            break
    sv = min((norm(r) for r in rows if any(r)), default=float('inf'))
    return ok, sv, norm(pv)

# ---------------------------------------------------------------------------
# Experiment driver
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

def load_hist():
    out = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, (p, b, n, lam, G), k1, m))
    return out

def trial(curve, m, k1_bound, seed, d_secret=None, beta=0):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    if d_secret is None:
        d_secret = random.Random(seed ^ 0xABCDEF).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    res = {}
    res['orig'] = run_kan_orig(sigs, n, lam, k1_bound, k2_bound, d_secret, beta)
    res['proj_cvp'] = run_proj_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret, beta)
    res['proj_kan'] = run_proj_kan(sigs, n, lam, k1_bound, k2_bound, d_secret, beta)
    res['proj_cvp_c'] = run_proj_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret,
                                     beta, center=True)
    res['proj_kan_c'] = run_proj_kan(sigs, n, lam, k1_bound, k2_bound, d_secret,
                                     beta, center=True)
    res['sigs'] = sigs
    res['d'] = d_secret
    res['k2b'] = k2_bound
    return res


def main():
    print("=" * 78)
    print("Thread 23 — projected-CVP reformulation of the GLV-HNP Phase-2 lattice")
    print("=" * 78)
    print(f"exact CVP available: {_CVP is not None}")

    hist = load_hist()

    # -- sanity: planted vector really lies in L_proj -----------------------
    print("\n" + "-" * 78)
    print("S0 — sanity: planted vector membership in L_proj")
    print("-" * 78)
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d_secret = random.Random(7).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2b, 42)
        check_proj_membership(sigs, n, lam, k1, k2b, d_secret)
        print(f"  {label:<18} m={m:<3} OK  (t + u_planted = (k1*S_K1 | k2*S_K2))")

    # -- T1: sv/pv on the three historical curves ---------------------------
    print("\n" + "-" * 78)
    print("T1 — is the planted vector lambda_1 now?  (sv/pv, seed=42)")
    print("     2026-07-29 baseline sv/pv: 0.603 / 0.517 / 0.422 (all < 1)")
    print("-" * 78)
    print(f"  {'curve':<18} {'K1':>3} {'m':>3} | {'ORIG':>7} {'PKAN':>7} "
          f"{'PKAN-C':>7} | {'PCVP':>7} {'PCVP-C':>7} | ok: O/PK/PKc/PC/PCc")
    print(f"  {'':<18} {'':>3} {'':>3} | {'--- sv/pv ---':>23} | "
          f"{'-- dist/pv --':>15} |")
    for label, curve, k1, m in hist:
        r = trial(curve, m, k1, 42)
        o_ok, o_sv, o_pv = r['orig']
        pk_ok, pk_sv, pk_pv = r['proj_kan']
        pkc_ok, pkc_sv, pkc_pv = r['proj_kan_c']
        pc_ok, _pc_nrm, pc_pv, pc_dist = r['proj_cvp']
        pcc_ok, _pcc_nrm, pcc_pv, pcc_dist = r['proj_cvp_c']
        flags = "/".join("T" if x else "F" for x in
                         (o_ok, pk_ok, pkc_ok, pc_ok, pcc_ok))
        print(f"  {label:<18} {k1:>3} {m:>3} | {o_sv/o_pv:>7.3f} "
              f"{pk_sv/pk_pv:>7.3f} {pkc_sv/pkc_pv:>7.3f} | "
              f"{pc_dist/pc_pv:>7.3f} {pcc_dist/pcc_pv:>7.3f} | {flags}")

    # -- T2: the K1 wall (reproduces the 2026-07-29 T4 grid) ----------------
    print("\n" + "-" * 78)
    print("T2 — K1 wall, 5 seeds, all three methods")
    print("     2026-07-29 T4 baseline (ORIG): 2557 -> 5,5,5,5,5,4,1,0")
    print("                                    2677 -> 5,5,5,2,0,0,0,0")
    print("-" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    for label, curve, _k1, m in hist[1:]:
        p, b, n, lam, G = curve
        print(f"\n  {label}  (n={n}, lam*={lam_star(lam, n):.4f}, m={m})")
        print(f"    {'method':<10} " + " ".join(f"K1={k:<3}" for k in K1_GRID))
        tallies = {k: [] for k in ('orig', 'proj_cvp', 'proj_kan',
                                   'proj_cvp_c', 'proj_kan_c')}
        for k1 in K1_GRID:
            cnt = {k: 0 for k in tallies}
            for s in SEEDS:
                r = trial(curve, m, k1, s)
                if r is None:
                    continue
                for k in cnt:
                    if r[k][0]: cnt[k] += 1
            for k in tallies:
                tallies[k].append(cnt[k])
        for k, nm in (('orig', 'KAN-ORIG'), ('proj_cvp', 'PROJ-CVP'),
                      ('proj_kan', 'PROJ-KAN'), ('proj_cvp_c', 'PROJ-CVP-C'),
                      ('proj_kan_c', 'PROJ-KAN-C')):
            print(f"    {nm:<10} " + " ".join(f"{v}/5  " for v in tallies[k]))

    # -- T3: how far past the old wall does PROJ-CVP go? --------------------
    print("\n" + "-" * 78)
    print("T3 — PROJ-CVP past the old wall: more signatures at large K1")
    print("-" * 78)
    label, curve, _k1, _m = hist[2]           # the lam*=0.07 failure curve
    p, b, n, lam, G = curve
    print(f"  {label}  n={n}  lam*={lam_star(lam, n):.4f}")
    print(f"    {'K1':>4} {'m':>4} {'eff':>7} {'m_thr':>6} | {'ORIG':>5} "
          f"{'PCVP':>5} {'PCVP-C':>6} {'PKAN':>5} {'PKAN-C':>6}")
    k2b = math.isqrt(n) + 1
    for k1 in (8, 12, 16, 24):
        eff = k1 * k2b / n
        m_thr = (math.ceil(math.log(n) / math.log(1.0 / eff))
                 if eff < 1.0 else float('inf'))
        for m in (10, 16, 24, 32):
            c = {k: 0 for k in ('orig', 'proj_cvp', 'proj_kan',
                                'proj_cvp_c', 'proj_kan_c')}
            for s in SEEDS:
                r = trial(curve, m, k1, s)
                if r is None:
                    continue
                for k in c:
                    if r[k][0]: c[k] += 1
            print(f"    {k1:>4} {m:>4} {eff:>7.3f} {str(m_thr):>6} | "
                  f"{c['orig']}/5  {c['proj_cvp']}/5  {c['proj_cvp_c']}/5  "
                  f"{c['proj_kan']}/5  {c['proj_kan_c']}/5")

    # -- T4: is the residual K1>=12 wall reduction-quality or structural? ---
    print("\n" + "-" * 78)
    print("T4 — BKZ preprocessing for PROJ-CVP-C past the K1=8 wall")
    print("-" * 78)
    label, curve, _k1, _m = hist[2]
    p, b, n, lam, G = curve
    print(f"  {label}  n={n}  lam*={lam_star(lam, n):.4f}")
    print(f"    {'K1':>4} {'m':>4} {'dim':>4} | {'LLL':>6} {'BKZ20':>7} {'BKZ30':>7}")
    for k1 in (8, 12, 16):
        for m in (16, 24, 32):
            row = []
            for beta in (0, 20, 30):
                c = 0
                for s in SEEDS:
                    k2b = math.isqrt(n) + 1
                    d_secret = random.Random(s ^ 0xABCDEF).randint(1, n - 1)
                    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2b, s)
                    if len(sigs) < m:
                        continue
                    if run_proj_cvp(sigs, n, lam, k1, k2b, d_secret, beta,
                                    center=True)[0]:
                        c += 1
                row.append(c)
            print(f"    {k1:>4} {m:>4} {2*m:>4} | {row[0]}/5{'':>3} "
                  f"{row[1]}/5{'':>4} {row[2]}/5")

    # -- T5: the lambda-block geometric window ------------------------------
    print("\n" + "-" * 78)
    print("T5 — the residual competitor: the 2-D lambda-block shortest vector")
    print("     mu = |shortest of <(n*S_K1, 0), (-lam*S_K1, S_K2)>|")
    print("     SVP-form recovery needs  m_thr < m < m_geo,  where m_geo is the")
    print("     largest m with E||pv_proj||(m) < mu.")
    print("-" * 78)
    for label, curve, _k1, _m in hist[1:]:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        print(f"\n  {label}  n={n}  lam*={lam_star(lam, n):.4f}  K2={k2b}")
        print(f"    {'K1':>4} {'eff':>7} {'mu/n':>7} {'m_thr':>6} {'m_geo':>6} "
              f"{'window':>8} | {'best PROJ-KAN-C over m':>22}")
        for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
            S_K1, _S_D, S_K2, S_KAN = scales(n, k1, k2b)
            mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
            eff = k1 * k2b / n
            m_thr = (math.ceil(math.log(n) / math.log(1.0 / eff))
                     if eff < 1.0 else float('inf'))
            m_geo = 0
            for mm in range(1, 200):
                e = math.sqrt(mm * (k1 * S_K1) ** 2 / 3.0
                              + mm * (k2b * S_K2) ** 2 / 3.0 + S_KAN ** 2)
                if e < mu:
                    m_geo = mm
                else:
                    break
            window = "EMPTY" if m_geo <= m_thr else f"{int(m_thr)}..{m_geo}"
            best, best_m = 0, None
            for mm in sorted({int(m_thr), int(m_thr) + 1, max(1, m_geo),
                              m_geo + 2, 10, 16}):
                if mm < 2 or mm > 40:
                    continue
                c = 0
                for s in SEEDS:
                    d_secret = random.Random(s ^ 0xABCDEF).randint(1, n - 1)
                    sigs = gen_signatures(G, d_secret, mm, n, lam, p, k1, k2b, s)
                    if len(sigs) < mm:
                        continue
                    if run_proj_kan(sigs, n, lam, k1, k2b, d_secret,
                                    center=True)[0]:
                        c += 1
                if c > best:
                    best, best_m = c, mm
            print(f"    {k1:>4} {eff:>7.3f} {mu/n:>7.2f} {str(m_thr):>6} "
                  f"{m_geo:>6} {window:>8} | {best}/5 at m={best_m}")

    # -- T6: same m-sweep on the OTHER historical curve (symmetry check) ----
    print("\n" + "-" * 78)
    print("T6 — m-sweep on 12-bit/2557 (lam*=0.34), all three methods")
    print("     T2 fixed m=8 at the historical value; does growing m rescue it?")
    print("-" * 78)
    label, curve, _k1, _m = hist[1]
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"  {label}  n={n}  lam*={lam_star(lam, n):.4f}")
    print(f"    {'K1':>4} {'m':>4} {'eff':>7} {'mu/n':>6} | {'ORIG':>5} "
          f"{'PCVP':>5} {'PCVP-C':>6} {'PKAN':>5} {'PKAN-C':>6}")
    for k1 in (4, 6, 8, 12):
        S_K1, _S_D, S_K2, _S_KAN = scales(n, k1, k2b)
        mu, _w = lambda_block_mu(n, lam, S_K1, S_K2)
        for m in (8, 16, 24, 32):
            c = {k: 0 for k in ('orig', 'proj_cvp', 'proj_kan',
                                'proj_cvp_c', 'proj_kan_c')}
            for s in SEEDS:
                r = trial(curve, m, k1, s)
                if r is None:
                    continue
                for k in c:
                    if r[k][0]: c[k] += 1
            print(f"    {k1:>4} {m:>4} {k1*k2b/n:>7.3f} {mu/n:>6.2f} | "
                  f"{c['orig']}/5  {c['proj_cvp']}/5  {c['proj_cvp_c']}/5  "
                  f"{c['proj_kan']}/5  {c['proj_kan_c']}/5")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == '__main__':
    main()
