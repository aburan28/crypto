"""
GLV-HNP Phase 2 — Thread 23: reformulate the lattice so the planted vector is
the target of a BDD/CVP instance rather than an SVP instance.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, finding T5):
  In the Phase-2 Kannan lattice (dim 2m+2, `build_glv_lattice` in
  glv_hnp_phase2_20bit.py:262) the shortest vector is ALWAYS the trivial
  vector n*S_D*e_m -- 100% of its energy in the d-column, |sv[m]|/n = 1.0000
  on every curve tested, and sv/pv in [0.34, 0.60].  The planted vector is
  therefore never lambda_1, and recovery is a coset/BDD condition, not SVP.
  No curve-level invariant can predict it.

Thread 23 reformulation.  The trivial vector spans the e_m (d-) direction, and
L contains n*e_m, so quotient it out: work in the projection along e_m.  That
removes the d-column entirely and turns the problem into a CVP in dimension 2m.

  Congruence:  A_i + B_i*d - lam*k2_i = k1_i  (mod n),  0<=k1_i<K1, 0<=k2_i<K2

  Put f_i := x_i + lam*y_i (mod n) and C_i := B_i * B_1^{-1} (mod n).  Then
      L' = { (x,y) in Z^m x Z^m : f_i = C_i * f_1 (mod n) for all i }
  is a full-rank lattice in Z^{2m} of determinant n^{m-1} (before scaling), and
  (x,y) = (k1 - A, k2) is a lattice point.  Explicit basis:

      g0     : x = (1, C_2, ..., C_m),               y = 0
      g1     : x = (0, C_2*lam, ..., C_m*lam) mod n, y = e_1
      g_j    : x = -lam*e_j,                         y = e_j      (j=2..m)
      h_i    : x = n*e_i,                            y = 0        (i=2..m)

  Columns are scaled by S_K1 (x-block) and S_K2 (y-block) exactly as before.
  CVP target (centred):
      t = (-A_i*S_K1 + S_K1*floor((K1-1)/2) | S_K2*floor((K2-1)/2))
  so that the error v - t has ZERO-MEAN coordinates of size ~n/sqrt(12) instead
  of the uncentred ~n/sqrt(3).  d is read straight off the recovered lattice
  point: d = B_1^{-1} * (x_1 + lam*y_1) mod n.

Experiments
  E1  sanity: planted point really is in L'; d read-out is exact; det = n^{m-1}.
  E2  geometry: ||err|| vs lambda_1(L') and vs the Gaussian heuristic, next to
      the old lattice's sv/pv.  Falsifier from the 2026-07-29 log: "if sv/pv
      rises above 1 after the reformulation ... the reformulation is a real
      improvement".
  E3  the T4 K1 grid (the 2026-07-29 wall measurement) re-run with four
      methods: OLD Kannan-SVP, NEW Babai-CVP, NEW Kannan-embedded-SVP, and NEW
      Kannan uncentred (to attribute any gain to centring vs projection).
  E4  is the wall predicted by a curve-independent lattice-geometric quantity
      nu_GH = ||err|| / GH(L')?  If one threshold in nu_GH reproduces the K1
      wall on curves with lam* = 0.070, 0.340 and 0.0068 alike, the wall is
      information-theoretic, not curve-dependent.

Run: python3 glv_hnp_phase2_projected_cvp.py
"""

import math
import random
from fractions import Fraction

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (verbatim from glv_hnp_phase2_20bit.py so the comparison to the
# 2026-07-29 T4 grid is exact)
# ---------------------------------------------------------------------------


def modinv(a, m):
    return pow(a, -1, m)


def ec_add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_mul(P, k, p):
    if k == 0:
        return None
    R, Q = None, P
    while k > 0:
        if k & 1:
            R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R


def tonelli_shanks(a, p):
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
    m2, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while True:
        if t == 0:
            return 0
        if t == 1:
            return r
        i, tmp = 0, t
        while tmp != 1:
            tmp = tmp * tmp % p
            i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p


def find_generator(p, b_param, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_param) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None


def find_b_for_n(p, n):
    for b_try in range(1, min(p, 500)):
        rng = random.Random(99)
        for _ in range(30):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                if ec_mul((x, y), n, p) is None:
                    return b_try
                break
    return None


def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    return (r1, r2)


def glv_eigenvalue(n):
    rr = glv_roots(n)
    if rr is None:
        return None
    r1, r2 = rr
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)


def lam_star(lam, n):
    return min(lam, n - lam) / n


# ---------------------------------------------------------------------------
# Signature generation (verbatim scheme from glv_hnp_phase2_20bit.py:233)
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
        if k_full == 0:
            continue
        R = ec_mul(G, k_full, p)
        if R is None:
            continue
        r = R[0] % n
        if r == 0:
            continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0:
            continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        if B % n == 0:
            continue
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs


# ---------------------------------------------------------------------------
# OLD lattice (baseline) -- verbatim from glv_hnp_phase2_20bit.py:262
# ---------------------------------------------------------------------------


def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
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
    return M, S_K1, S_D, S_K2, S_KANNAN


def build_glv_lattice_centred(sigs, n, lam, k1_bound, k2_bound):
    """OLD lattice with the Kannan row shifted so the planted k1/k2 blocks are
    zero-mean.  Isolates the effect of centring from the effect of projection."""
    M, S_K1, S_D, S_K2, S_KAN = build_glv_lattice(sigs, n, lam, k1_bound,
                                                  k2_bound)
    m = len(sigs)
    h1 = (k1_bound - 1) // 2
    h2 = (k2_bound - 1) // 2
    for i in range(m):
        M[2 * m + 1][i] = (sigs[i]['A'] - h1) * S_K1
        M[2 * m + 1][m + 1 + i] = -h2 * S_K2
    return M, S_K1, S_D, S_K2, S_KAN


def old_recover(sigs, n, lam, K1, K2, d_secret, use_bkz=False, beta=20,
                centred=False):
    m = len(sigs)
    builder = build_glv_lattice_centred if centred else build_glv_lattice
    M, S_K1, S_D, S_K2, S_KAN = builder(sigs, n, lam, K1, K2)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        if (sign * A[i][m]) % n == d_secret:
            return True
    return False


def old_geometry(sigs, n, lam, K1, K2, d_secret):
    """||shortest vector after LLL|| / ||planted vector|| in the OLD lattice."""
    m = len(sigs)
    M, S_K1, S_D, S_K2, S_KAN = build_glv_lattice(sigs, n, lam, K1, K2)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    sv2 = min(sum(A[i][j] ** 2 for j in range(dim)) for i in range(dim)
              if any(A[i][j] for j in range(dim)))
    pv2 = (sum((s['k1'] * S_K1) ** 2 for s in sigs)
           + (d_secret * S_D) ** 2
           + sum((s['k2'] * S_K2) ** 2 for s in sigs)
           + S_KAN ** 2)
    return math.sqrt(sv2) / math.sqrt(pv2)


# ---------------------------------------------------------------------------
# NEW projected lattice L' (dim 2m) + CVP target
# ---------------------------------------------------------------------------


def build_projected_lattice(sigs, n, lam, K1, K2, centred=True):
    """Basis of L' (rows, dim 2m) plus the CVP target t.

    Coordinates: 0..m-1 are the congruence columns (scale S_K1),
                 m..2m-1 are the k2 columns (scale S_K2).
    """
    m = len(sigs)
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    B1inv = modinv(sigs[0]['B'], n)
    C = [sigs[i]['B'] * B1inv % n for i in range(m)]  # C[0] == 1

    rows = []
    # g0: x = (1, C_2, ..., C_m), y = 0
    g0 = [0] * (2 * m)
    for i in range(m):
        g0[i] = C[i] * S_K1
    rows.append(g0)
    # g1: x = (0, C_2*lam, ..., C_m*lam) mod n, y = e_1
    g1 = [0] * (2 * m)
    for i in range(1, m):
        g1[i] = (C[i] * lam % n) * S_K1
    g1[m] = S_K2
    rows.append(g1)
    # g_j (j = 2..m): x = -lam*e_j, y = e_j
    for j in range(1, m):
        g = [0] * (2 * m)
        g[j] = -lam * S_K1
        g[m + j] = S_K2
        rows.append(g)
    # h_i (i = 2..m): x = n*e_i
    for i in range(1, m):
        h = [0] * (2 * m)
        h[i] = n * S_K1
        rows.append(h)

    c1 = S_K1 * ((K1 - 1) // 2) if centred else 0
    c2 = S_K2 * ((K2 - 1) // 2) if centred else 0
    t = [(-sigs[i]['A']) * S_K1 + c1 for i in range(m)] + [c2] * m
    return rows, t, S_K1, S_K2, B1inv


def planted_point(sigs, n, S_K1, S_K2):
    """The lattice point (k1 - A, k2), scaled."""
    m = len(sigs)
    return ([(s['k1'] - s['A']) * S_K1 for s in sigs]
            + [s['k2'] * S_K2 for s in sigs])


def read_d(v, m, n, lam, S_K1, S_K2, B1inv):
    """d = B_1^{-1} (x_1 + lam*y_1) mod n from a lattice point v."""
    if v[0] % S_K1 or v[m] % S_K2:
        return None
    x1 = v[0] // S_K1
    y1 = v[m] // S_K2
    return B1inv * (x1 + lam * y1) % n


# ---------------------------------------------------------------------------
# Exact Babai nearest-plane on an LLL-reduced basis
# ---------------------------------------------------------------------------


def gram_schmidt_exact(basis):
    n = len(basis)
    bstar = []
    mu = [[Fraction(0)] * n for _ in range(n)]
    norms = []
    for i in range(n):
        v = [Fraction(x) for x in basis[i]]
        for j in range(i):
            if norms[j] == 0:
                continue
            num = sum(Fraction(basis[i][k]) * bstar[j][k] for k in range(len(v)))
            mu[i][j] = num / norms[j]
            for k in range(len(v)):
                v[k] -= mu[i][j] * bstar[j][k]
        bstar.append(v)
        norms.append(sum(x * x for x in v))
    return bstar, norms


def babai_nearest_plane(basis, target):
    """Return the lattice vector produced by Babai's nearest-plane on `basis`
    (which should already be LLL-reduced)."""
    n = len(basis)
    dim = len(target)
    bstar, norms = gram_schmidt_exact(basis)
    b = [Fraction(x) for x in target]
    coeffs = [0] * n
    for i in range(n - 1, -1, -1):
        if norms[i] == 0:
            continue
        c = sum(b[k] * bstar[i][k] for k in range(dim)) / norms[i]
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = ci
        for k in range(dim):
            b[k] -= ci * Fraction(basis[i][k])
    v = [sum(coeffs[i] * basis[i][k] for i in range(n)) for k in range(dim)]
    return v


def lll_rows(rows):
    d = len(rows)
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    return [[A[i][j] for j in range(len(rows[0]))] for i in range(d)]


def norm2(v):
    return sum(x * x for x in v)


# ---------------------------------------------------------------------------
# NEW recovery methods
# ---------------------------------------------------------------------------


def new_babai_recover(sigs, n, lam, K1, K2, d_secret, centred=True):
    m = len(sigs)
    rows, t, S_K1, S_K2, B1inv = build_projected_lattice(sigs, n, lam, K1, K2,
                                                         centred)
    red = lll_rows(rows)
    v = babai_nearest_plane(red, t)
    return read_d(v, m, n, lam, S_K1, S_K2, B1inv) == d_secret


def new_kannan_recover(sigs, n, lam, K1, K2, d_secret, centred=True,
                       tau=None, use_bkz=False, beta=20):
    """Kannan-embed the CVP into an SVP of dimension 2m+1 and scan."""
    m = len(sigs)
    rows, t, S_K1, S_K2, B1inv = build_projected_lattice(sigs, n, lam, K1, K2,
                                                         centred)
    if tau is None:
        tau = max(1, (n // K1) * ((K1 - 1) // 2) // 2 or 1)
    dim = 2 * m + 1
    emb = [r + [0] for r in rows] + [t + [tau]]
    A = IntegerMatrix.from_matrix(emb)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != tau:
            continue
        sign = 1 if last > 0 else -1
        # row = sign*(t|tau) - v  =>  v = sign*t - sign*row
        v = [sign * t[k] - sign * A[i][k] for k in range(2 * m)]
        dc = read_d(v, m, n, lam, S_K1, S_K2, B1inv)
        if dc == d_secret:
            return True
    return False


def new_geometry(sigs, n, lam, K1, K2, centred=True):
    """||err|| vs lambda_1(L') (LLL estimate) and vs the Gaussian heuristic."""
    m = len(sigs)
    rows, t, S_K1, S_K2, _ = build_projected_lattice(sigs, n, lam, K1, K2,
                                                      centred)
    red = lll_rows(rows)
    lam1 = math.sqrt(min(norm2(r) for r in red if any(r)))
    pv = planted_point(sigs, n, S_K1, S_K2)
    err = math.sqrt(norm2([pv[k] - t[k] for k in range(2 * m)]))
    dim = 2 * m
    log2vol = (m - 1) * math.log2(n) + m * math.log2(S_K1) + m * math.log2(S_K2)
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * 2 ** (log2vol / dim)
    return err, lam1, gh


def gh_predictor(n, K1, K2, m):
    """nu_GH = E||err|| / GH(L'), a curve-independent closed form."""
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    dim = 2 * m
    err = math.sqrt(2 * m / 12.0) * n  # centred: E[coord^2] ~ n^2/12
    log2vol = (m - 1) * math.log2(n) + m * math.log2(S_K1) + m * math.log2(S_K2)
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * 2 ** (log2vol / dim)
    return err / gh


# ---------------------------------------------------------------------------
# Curve set
# ---------------------------------------------------------------------------

# label, p, b, n, lam  -- the historical Phase-2 curves (log 2026-06-15 /
# 2026-07-26 / 2026-07-29) plus the 17-bit lam*=0.0068 curve from 2026-07-29 T3.
HIST = [
    ("8-bit/199", 211, 2, 199, 106),
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]

SEEDS = [42, 1234, 9999, 555, 31337]


def load_curves():
    out = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, p, b, n, lam, G))
    # 17-bit lam*=0.0068 curve (2026-07-29 T3): p=65713, n=65269, lam=442
    p, n = 65713, 65269
    lam = glv_eigenvalue(n)
    if lam is not None and (lam * lam + lam + 1) % n == 0:
        b = find_b_for_n(p, n)
        if b is not None:
            G = find_generator(p, b, n)
            if G is not None:
                out.append(("17-bit/65713", p, b, n, lam, G))
    return out


def make_instance(curve, m, K1, seed):
    label, p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed ^ 0x5EED)
    d = rng.randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, K1, K2, seed)
    if len(sigs) < m:
        return None
    return sigs, d, K2


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("=" * 78)
    print("Thread 23 — projected (d-quotiented) GLV-HNP lattice + CVP")
    print("=" * 78)

    curves = load_curves()
    print("\nCurves:")
    for (label, p, b, n, lam, G) in curves:
        print(f"  {label:<14} p={p:<7} n={n:<7} ({n.bit_length()}b) "
              f"lam={lam:<7} lam*={lam_star(lam, n):.4f}")

    # ---------------------------------------------------------------- E1
    print("\n" + "-" * 78)
    print("E1: sanity — planted point in L', exact d read-out, det(L') = n^(m-1)")
    print("-" * 78)
    print(f"{'curve':<14} {'m':>3} {'K1':>4} {'planted in L?':>14} "
          f"{'d readout':>10} {'log2 det':>10} {'expected':>10}")
    for curve in curves:
        label, p, b, n, lam, G = curve
        m, K1 = 6, 4
        inst = make_instance(curve, m, K1, 42)
        if inst is None:
            print(f"{label:<14} instance generation failed")
            continue
        sigs, d, K2 = inst
        rows, t, S_K1, S_K2, B1inv = build_projected_lattice(sigs, n, lam,
                                                             K1, K2)
        pv = planted_point(sigs, n, S_K1, S_K2)
        # membership test: solve pv = z * Basis over Q, check z integral
        inlat = solve_membership(rows, pv)
        dread = read_d(pv, m, n, lam, S_K1, S_K2, B1inv)
        ld = log2_det(rows)
        exp_ld = (m - 1) * math.log2(n) + m * math.log2(S_K1) + m * math.log2(S_K2)
        print(f"{label:<14} {m:>3} {K1:>4} {str(inlat):>14} "
              f"{('OK' if dread == d else 'WRONG'):>10} {ld:>10.2f} "
              f"{exp_ld:>10.2f}")

    # ---------------------------------------------------------------- E2
    print("\n" + "-" * 78)
    print("E2: geometry — OLD sv/pv (planted never lambda_1) vs NEW err/lambda_1")
    print("-" * 78)
    print("NEW: recovery needs ||err|| below the decoding radius ~ lambda_1/2.")
    print(f"\n{'curve':<14} {'K1':>4} {'OLD sv/pv':>10} {'NEW err':>12} "
          f"{'NEW lam1':>12} {'lam1/err':>9} {'err/GH':>8}")
    m = 12
    for curve in curves:
        label, p, b, n, lam, G = curve
        for K1 in (2, 4, 8, 16):
            inst = make_instance(curve, m, K1, 42)
            if inst is None:
                continue
            sigs, d, K2 = inst
            old = old_geometry(sigs, n, lam, K1, K2, d)
            err, lam1, gh = new_geometry(sigs, n, lam, K1, K2)
            print(f"{label:<14} {K1:>4} {old:>10.3f} {err:>12.3e} "
                  f"{lam1:>12.3e} {lam1 / err:>9.3f} {err / gh:>8.3f}")

    # ---------------------------------------------------------------- E3
    print("\n" + "-" * 78)
    print("E3: the 2026-07-29 T4 K1 grid, four methods, m=12, 5 seeds")
    print("-" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
    METHODS = [
        ("OLD  Kannan-SVP  ", lambda s, n, l, k1, k2, d:
            old_recover(s, n, l, k1, k2, d)),
        ("NEW  Babai-CVP   ", lambda s, n, l, k1, k2, d:
            new_babai_recover(s, n, l, k1, k2, d, centred=True)),
        ("NEW  Kannan-embed", lambda s, n, l, k1, k2, d:
            new_kannan_recover(s, n, l, k1, k2, d, centred=True)),
        ("NEW  Kannan uncen", lambda s, n, l, k1, k2, d:
            new_kannan_recover(s, n, l, k1, k2, d, centred=False)),
        ("OLD  centred     ", lambda s, n, l, k1, k2, d:
            old_recover(s, n, l, k1, k2, d, centred=True)),
    ]
    m = 12
    grid_results = {}
    for curve in curves:
        label, p, b, n, lam, G = curve
        print(f"\n  {label}  (lam*={lam_star(lam, n):.4f}, n={n}, m={m})")
        hdr = "  " + f"{'method':<18}" + "".join(f"{'K1=' + str(k):>8}"
                                                 for k in K1_GRID)
        print(hdr)
        for mname, fn in METHODS:
            cells = []
            for K1 in K1_GRID:
                wins = 0
                tot = 0
                for seed in SEEDS:
                    inst = make_instance(curve, m, K1, seed)
                    if inst is None:
                        continue
                    sigs, d, K2 = inst
                    tot += 1
                    try:
                        if fn(sigs, n, lam, K1, K2, d):
                            wins += 1
                    except Exception as exc:      # noqa: BLE001
                        print(f"      [{mname} K1={K1} seed={seed}] {exc}")
                cells.append((wins, tot))
                grid_results[(label, mname, K1)] = (wins, tot)
            print("  " + f"{mname:<18}"
                  + "".join(f"{str(w) + '/' + str(t):>8}" for w, t in cells))

    # ---------------------------------------------------------------- E4
    print("\n" + "-" * 78)
    print("E4: is the wall predicted by nu_GH = E||err|| / GH(L')?")
    print("-" * 78)
    print("nu_GH is a closed form in (n, K1, K2, m) only — no curve data at all.")
    print(f"\n{'curve':<14} {'K1':>4} {'eff=K1K2/n':>11} {'nu_GH':>8} "
          f"{'OLD':>7} {'NEW-Bab':>8} {'NEW-Kan':>8}")
    for curve in curves:
        label, p, b, n, lam, G = curve
        K2 = math.isqrt(n) + 1
        for K1 in K1_GRID:
            key_o = grid_results.get((label, "OLD  Kannan-SVP  ", K1))
            key_b = grid_results.get((label, "NEW  Babai-CVP   ", K1))
            key_k = grid_results.get((label, "NEW  Kannan-embed", K1))
            if key_o is None:
                continue
            nu = gh_predictor(n, K1, K2, 12)
            print(f"{label:<14} {K1:>4} {K1 * K2 / n:>11.4f} {nu:>8.3f} "
                  f"{str(key_o[0]) + '/' + str(key_o[1]):>7} "
                  f"{str(key_b[0]) + '/' + str(key_b[1]):>8} "
                  f"{str(key_k[0]) + '/' + str(key_k[1]):>8}")

    # ---------------------------------------------------------------- E5
    print("\n" + "-" * 78)
    print("E5: how well does the curve-independent nu_GH predict each method?")
    print("-" * 78)
    print("A cell (curve,K1) counts as SUCCESS if >=3/5 seeds recover.")
    print("For each method: best single nu_GH threshold, its error count, and")
    print("the ambiguity band [max nu_GH of a success, min nu_GH of a failure].")
    print(f"\n{'method':<18} {'best thr':>9} {'errors':>7} {'cells':>6} "
          f"{'max nu(succ)':>13} {'min nu(fail)':>13} {'monotone':>9}")
    for mname, _ in METHODS:
        pts = []
        for curve in curves:
            label, p, b, n, lam, G = curve
            K2 = math.isqrt(n) + 1
            for K1 in K1_GRID:
                cell = grid_results.get((label, mname, K1))
                if cell is None or cell[1] == 0:
                    continue
                pts.append((gh_predictor(n, K1, K2, 12),
                            1 if cell[0] * 2 >= cell[1] else 0))
        pts.sort()
        best_err, best_thr = len(pts) + 1, None
        for idx in range(len(pts) + 1):
            thr = pts[idx][0] if idx < len(pts) else float('inf')
            err = (sum(1 for nu, ok in pts if nu < thr and not ok)
                   + sum(1 for nu, ok in pts if nu >= thr and ok))
            if err < best_err:
                best_err, best_thr = err, thr
        succ = [nu for nu, ok in pts if ok]
        fail = [nu for nu, ok in pts if not ok]
        mono = (max(succ) <= min(fail)) if succ and fail else True
        print(f"{mname:<18} {best_thr:>9.3f} {best_err:>7} {len(pts):>6} "
              f"{(max(succ) if succ else float('nan')):>13.3f} "
              f"{(min(fail) if fail else float('nan')):>13.3f} "
              f"{str(mono):>9}")

    # ---------------------------------------------------------------- E6
    print("\n" + "-" * 78)
    print("E6: OUT-OF-SAMPLE — the closed form predicts how the wall moves in m")
    print("-" * 78)
    print("nu_GH ~ sqrt(2*pi*e/12) * sqrt(K1*K2/n) * n^(1/(2m)), so the wall")
    print("K1* solving nu_GH = THR must grow with m like n^(1 - 1/m).")
    print("E3/E5 fixed m=12; this sweeps m and compares predicted vs observed.")
    THR = 1.25
    for curve in curves[1:3]:            # the two 12-bit historical curves
        label, p, b, n, lam, G = curve
        K2 = math.isqrt(n) + 1
        print(f"\n  {label}  (lam*={lam_star(lam, n):.4f}, n={n}, K2={K2})")
        print(f"  {'m':>3} {'K1* pred':>9} {'observed wall (K1 at which <3/5)':>36}"
              f" {'nu_GH at wall':>14}")
        for mm in (6, 8, 12, 18, 24):
            pred = predicted_wall(n, K2, mm, THR)
            grid = sorted(set([2, 4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 64]))
            obs = None
            obs_nu = float('nan')
            for K1 in grid:
                wins = tot = 0
                for seed in SEEDS:
                    inst = make_instance(curve, mm, K1, seed)
                    if inst is None:
                        continue
                    sigs, d, K2b = inst
                    tot += 1
                    if new_kannan_recover(sigs, n, lam, K1, K2b, d,
                                          centred=True):
                        wins += 1
                if tot and wins * 2 < tot:
                    obs = K1
                    obs_nu = gh_predictor(n, K1, K2, mm)
                    break
            obs_s = str(obs) if obs is not None else ">64"
            print(f"  {mm:>3} {pred:>9.1f} {obs_s:>36} {obs_nu:>14.3f}")

    print("\nDone.")


def predicted_wall(n, K2, m, thr=1.25):
    """K1 at which nu_GH crosses `thr` (closed form, no curve data)."""
    eff = (thr / math.sqrt(2 * math.pi * math.e / 12.0)) ** 2 * n ** (-1.0 / m)
    return eff * n / K2


# ---------------------------------------------------------------------------
# small exact-linear-algebra helpers (used by E1)
# ---------------------------------------------------------------------------


def solve_membership(rows, v):
    """True iff v is an integer combination of `rows` (rows square & full rank)."""
    d = len(rows)
    aug = [[Fraction(rows[i][j]) for i in range(d)] + [Fraction(v[j])]
           for j in range(d)]          # transpose: columns are basis vectors
    # Gaussian elimination
    r = 0
    for c in range(d):
        piv = None
        for i in range(r, d):
            if aug[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        aug[r], aug[piv] = aug[piv], aug[r]
        pv = aug[r][c]
        aug[r] = [x / pv for x in aug[r]]
        for i in range(d):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [aug[i][k] - f * aug[r][k] for k in range(d + 1)]
        r += 1
    if r < d:
        return False
    sol = [aug[i][d] for i in range(d)]
    return all(x.denominator == 1 for x in sol)


def log2_det(rows):
    d = len(rows)
    mat = [[Fraction(x) for x in row] for row in rows]
    det = Fraction(1)
    for c in range(d):
        piv = None
        for i in range(c, d):
            if mat[i][c] != 0:
                piv = i
                break
        if piv is None:
            return float('-inf')
        if piv != c:
            mat[c], mat[piv] = mat[piv], mat[c]
            det = -det
        det *= mat[c][c]
        inv = mat[c][c]
        for i in range(c + 1, d):
            if mat[i][c] != 0:
                f = mat[i][c] / inv
                mat[i] = [mat[i][k] - f * mat[c][k] for k in range(d)]
    return math.log2(abs(det))


if __name__ == '__main__':
    main()
