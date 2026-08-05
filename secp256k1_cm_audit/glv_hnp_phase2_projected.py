"""
GLV-HNP Phase 2, Thread 23 — reformulate the lattice so the planted vector is lambda_1.

Context (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, experiment T5):
  In the Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` (`build_glv_lattice`) the
  shortest vector after LLL is *always* the trivial vector n*S_D*e_m: 100% of its
  energy sits in the d-readout column, |sv[m]|/n = 1.0000 exactly.  It carries no
  information (d is only defined mod n) and no choice of S_D removes it, because
  both it and the planted vector scale linearly in S_D.  Recovery is therefore not
  an SVP condition but a BDD/coset condition.

  Proposed falsifier (2026-07-29, "Next step proposal"):
    project the lattice along e_m (quotient out the trivial n*e_m direction), or
    replace the Kannan embedding by an explicit CVP call.
    * if sv/pv rises above 1 AND the K1 wall on the lam*=0.07 curve moves outward
      (currently K1 ~ 4-6)  =>  the reformulation is a real improvement;
    * if the wall stays at K1 ~ 4-6  =>  the wall is information-theoretic and
      Phase 2 is at its ceiling.

This script builds three variants and runs that falsifier.

  V0  baseline           dim 2m+2, columns [k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan]
  V1  projected          delete the d column  ->  rank 2m+1 in 2m+1 columns.
                         The trivial vector maps to 0 and leaves the lattice.
                         d is re-derived from (k1_0, k2_0) instead of being read off.
  V2  explicit CVP       drop the d column AND the Kannan row/column -> rank 2m in
                         2m columns; solve BDD directly against t = (A_i*S_K1, 0).

Plus an analytic BDD model (E3/E4): rho = ||e|| / GH(L') where L' is the V2 lattice.
rho is not a curve-level invariant -- it depends on (m, K1, K2) -- so it is not
excluded by the 2026-07-29 finding that no curve-level invariant separates C1 from C2.

Run: python3 glv_hnp_phase2_projected.py
Deps: fpylll (+ cysignals), sympy.
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py so the comparison is exact.
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
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    """Symmetric size of lam: min(lam, n-lam)/n in [0, 0.5]."""
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Phase-2 instance: scales, signatures, lattices
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- verbatim from glv_hnp_phase2_20bit.py:262."""
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


def build_v0(sigs, n, lam, k1_bound, k2_bound):
    """Baseline Phase-2 lattice, dim 2m+2.  Columns:
         0..m-1   k1_i * S_K1
         m        d * S_D           <-- the readout column that hosts n*e_m
         m+1..2m  k2_i * S_K2
         2m+1     kannan
    """
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

def build_v1(sigs, n, lam, k1_bound, k2_bound):
    """V0 with column m (the d readout) deleted.  2m+2 generators, rank 2m+1.
    Columns:  0..m-1 k1 | m..2m-1 k2 | 2m kannan.

    Rank drops by exactly 1 because  n*row_m  =  sum_i B_i * row_i  once the
    d column is gone; that identity is precisely the trivial vector n*e_m,
    so the projection quotients it out.  fpylll's LLL accepts the dependent
    generating set and emits one zero row.
    """
    M0 = build_v0(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    return [row[:m] + row[m + 1:] for row in M0]

def build_v2(sigs, n, lam, k1_bound, k2_bound):
    """CVP formulation.  Lattice L' = rows 0..2m of V0, with the d column and the
    kannan column both deleted:  2m+1 generators, rank 2m, in 2m columns
    (0..m-1 k1 | m..2m-1 k2).  Target t = (A_i*S_K1, 0, ..., 0).

    The planted vector satisfies  t + u = (k1_i*S_K1, k2_i*S_K2)  for the lattice
    vector u = d*row_m + sum_i k2_i*row_{m+1+i} - sum_i c_i*row_i in L', so the
    error is exactly what BDD must find.  No trivial vector exists in L'.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B = []
    for i in range(m):                                    # rows 0..m-1
        r = [0] * (2 * m); r[i] = n * S_K1; B.append(r)
    r = [0] * (2 * m)                                     # row m (the d row)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    B.append(r)
    for i in range(m):                                    # rows m+1..2m
        r = [0] * (2 * m); r[i] = -lam * S_K1; r[m + i] = S_K2; B.append(r)
    t = [sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    return B, t


def planted_vector_v0(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def planted_vector_v1(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def error_vector_v2(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return ([sigs[i]['k1'] * S_K1 for i in range(m)] +
            [sigs[i]['k2'] * S_K2 for i in range(m)])

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def d_from_k(sigs, n, lam, k1_0, k2_0):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    k_full = (k1_0 + lam * k2_0) % n
    try:
        return (k_full - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
    except ValueError:
        return None

def recover_v0(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return True
    return False

def recover_v1(rows, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_0, k2_0) off a row whose kannan coordinate is +-S_KAN, then
    reconstruct d.  The d column no longer exists, so this is the only route."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        if (sign * row[0]) % S_K1 or (sign * row[m]) % S_K2:
            continue
        k1_0 = (sign * row[0]) // S_K1
        k2_0 = (sign * row[m]) // S_K2
        d_cand = d_from_k(sigs, n, lam, k1_0, k2_0)
        if d_cand is not None and d_cand == d_secret:
            return True
    return False

def recover_v2(err, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for sign in (1, -1):
        e = [sign * x for x in err]
        if e[0] % S_K1 or e[m] % S_K2:
            continue
        d_cand = d_from_k(sigs, n, lam, e[0] // S_K1, e[m] // S_K2)
        if d_cand is not None and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Reduction drivers
# ---------------------------------------------------------------------------

def lll_rows(M):
    ncols = len(M[0])
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]

def nonzero(rows):
    return [r for r in rows if any(r)]

def run_v0(sigs, n, lam, K1, K2, d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    rows = lll_rows(build_v0(sigs, n, lam, K1, K2))
    nz = nonzero(rows)
    pv = norm(planted_vector_v0(sigs, d_secret, n, K1, K2))
    sv = min(norm(r) for r in nz)
    ok = recover_v0(rows, m, n, S_KAN, d_secret)
    # energy fraction of the shortest vector in the d column
    sh = min(nz, key=norm)
    d_frac = (float(sh[m]) ** 2) / (norm(sh) ** 2) if norm(sh) else float('nan')
    return ok, pv, sv, d_frac

def run_v1(sigs, n, lam, K1, K2, d_secret):
    m = len(sigs)
    rows = lll_rows(build_v1(sigs, n, lam, K1, K2))
    nz = nonzero(rows)
    pv = norm(planted_vector_v1(sigs, n, K1, K2))
    sv = min(norm(r) for r in nz)
    ok = recover_v1(rows, sigs, m, n, lam, K1, K2, d_secret)
    nzeros = len(rows) - len(nz)
    return ok, pv, sv, nzeros

def run_v2(sigs, n, lam, K1, K2, d_secret):
    """Babai nearest-plane (fpylll CVP.babai) and exact CVP on L'."""
    m = len(sigs)
    B, t = build_v2(sigs, n, lam, K1, K2)
    rows = nonzero(lll_rows(B))
    A = IntegerMatrix.from_matrix(rows)
    neg_t = [-x for x in t]
    e_true = error_vector_v2(sigs, n, K1, K2)
    out = {}
    for tag, fn in (('babai', CVP.babai), ('exact', CVP.closest_vector)):
        try:
            u = list(fn(A, neg_t))
        except Exception as exc:                       # pragma: no cover
            out[tag] = (False, float('nan'), str(exc))
            continue
        e = [t[i] + u[i] for i in range(2 * m)]
        out[tag] = (recover_v2(e, sigs, m, n, lam, K1, K2, d_secret), norm(e), None)
    return out, norm(e_true), rows

# ---------------------------------------------------------------------------
# BDD model:  rho = ||e|| / GH(L')
# ---------------------------------------------------------------------------

def det_lprime(n, m, K1, K2):
    """det(L') = (n*S_K1)^m * S_K2^m / n.

    Derivation: det(L_V0) = (n*S_K1)^m * S_D * S_K2^m * S_KAN (triangular).
    ker(pi) cap L = < n*S_D*e_m > (a*B_i = 0 mod n for all i forces a = 0 mod n),
    so det(pi(L)) = det(L)/(n*S_D); intersecting with {kannan = 0} divides by S_KAN.
    """
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    return (n * S_K1) ** m * S_K2 ** m // n

def gh(det, dim):
    """Gaussian heuristic for lambda_1."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * det ** (1.0 / dim)

def rho_bdd(n, m, K1, K2, e_norm=None):
    """rho = ||e|| / GH(L').  BDD is expected to be solvable iff rho < 1."""
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    if e_norm is None:                       # expectation, k1~U[0,K1), k2~U[0,K2)
        e_norm = math.sqrt(m * (K1 * S_K1) ** 2 / 3.0 + m * (K2 * S_K2) ** 2 / 3.0)
    return e_norm / gh(det_lprime(n, m, K1, K2), 2 * m)

def rho_closed_form(n, m, K1, K2):
    """sqrt(2*pi*e/3) * sqrt(K1*K2/n) * n^(1/(2m)), the S_* -> exact limit."""
    return math.sqrt(2 * math.pi * math.e / 3) * math.sqrt(K1 * K2 / n) * n ** (1.0 / (2 * m))

def predicted_k1_wall(n, m, K2):
    """Largest K1 with rho_closed_form < 1."""
    return (3 * n) / (2 * math.pi * math.e * K2 * n ** (1.0 / m))

# ---------------------------------------------------------------------------
# Curve construction
# ---------------------------------------------------------------------------

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

def search_curves(lo, hi, per_bin=1, nbins=10, max_primes=100000):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by lam*."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
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
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,              p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",         211,  2, 199,  106,  2,  6),
    ("12-bit/2557",       2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  8,  10),
]

print("=" * 78)
print("Thread 23 - projected / CVP reformulation of the GLV-HNP Phase-2 lattice")
print("=" * 78)

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))
    print(f"  {label:<18} n={n:<6} lam={lam:<6} lam*={lam_star(lam, n):.4f} "
          f"K1_hist={k1} m={m}")


# ---------------------------------------------------------------------------
# E0 -- sanity: does the projection actually remove the trivial vector?
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E0: does deleting the d column remove the trivial vector n*e_m?")
print("-" * 78)
print(f"{'curve':<18} {'K1':>3} {'V0 sv/pv':>9} {'V0 d-energy':>12} "
      f"{'V1 sv/pv':>9} {'V1 zeros':>9} {'rank drop ok':>13}")
for label, curve, K1, m in hist:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    d_secret = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, 42)
    ok0, pv0, sv0, dfrac = run_v0(sigs, n, lam, K1, K2, d_secret)
    ok1, pv1, sv1, nz = run_v1(sigs, n, lam, K1, K2, d_secret)
    print(f"{label:<18} {K1:>3} {sv0/pv0:>9.3f} {dfrac:>12.4f} "
          f"{sv1/pv1:>9.3f} {nz:>9} {str(nz == 1):>13}")

# determinant identity check
print("\n  det(L') identity check (formula vs exact Gram determinant):")
label, curve, K1, m = hist[0]
p, b, n, lam, G = curve
K2 = math.isqrt(n) + 1
d_secret = random.Random(42 + 7777).randint(1, n - 1)
sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, 42)
B2, _t = build_v2(sigs, n, lam, K1, K2)
basis = nonzero(lll_rows(B2))
gram = sympy.Matrix(basis) * sympy.Matrix(basis).T
det_exact = sympy.sqrt(gram.det())
print(f"    formula  (n*S_K1)^m * S_K2^m / n = {det_lprime(n, m, K1, K2)}")
print(f"    exact    sqrt(det(B B^T))        = {det_exact}")
print(f"    match: {sympy.simplify(det_exact - det_lprime(n, m, K1, K2)) == 0}")


# ---------------------------------------------------------------------------
# E1 -- the falsifier: K1 wall, V0 vs V1 vs V2
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E1: K1 wall for each variant (wins out of 5 seeds)")
print("-" * 78)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

wall_rows = []
for label, curve, _K1h, m in hist:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"\n  {label}  (n={n}, m={m}, K2={K2}, lam*={lam_star(lam, n):.3f})")
    print(f"    {'K1':>3} {'eff':>6} {'rho':>6} {'V0':>5} {'V1':>5} "
          f"{'V2-bab':>7} {'V2-cvp':>7} {'V0 sv/pv':>9} {'V1 sv/pv':>9}")
    for K1 in K1_GRID:
        if K1 >= n:
            continue
        w0 = w1 = w2b = w2c = 0
        r0 = []
        r1 = []
        for seed in SEEDS:
            d_secret = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, seed)
            if len(sigs) < m:
                continue
            ok0, pv0, sv0, _ = run_v0(sigs, n, lam, K1, K2, d_secret)
            ok1, pv1, sv1, _ = run_v1(sigs, n, lam, K1, K2, d_secret)
            out2, _en, _bs = run_v2(sigs, n, lam, K1, K2, d_secret)
            w0 += bool(ok0); w1 += bool(ok1)
            w2b += bool(out2['babai'][0]); w2c += bool(out2['exact'][0])
            r0.append(sv0 / pv0); r1.append(sv1 / pv1)
        eff = K1 * K2 / n
        rho = rho_bdd(n, m, K1, K2)
        print(f"    {K1:>3} {eff:>6.3f} {rho:>6.3f} {w0:>5} {w1:>5} "
              f"{w2b:>7} {w2c:>7} {sum(r0)/len(r0):>9.3f} {sum(r1)/len(r1):>9.3f}")
        wall_rows.append((label, n, m, K1, K2, eff, rho, w0, w1, w2b, w2c,
                          sum(r0)/len(r0), sum(r1)/len(r1)))


# ---------------------------------------------------------------------------
# E2 -- predicted vs observed K1 wall
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E2: predicted K1 wall (rho < 1) vs observed")
print("-" * 78)

def observed_wall(rows, label, col):
    """Largest K1 in the grid with >= 3/5 wins (majority)."""
    best = None
    for r in rows:
        if r[0] != label:
            continue
        if r[col] >= 3:
            best = r[3]
    return best

print(f"{'curve':<18} {'lam*':>6} {'m':>3} {'K1* pred':>9} "
      f"{'obs V0':>7} {'obs V1':>7} {'obs V2':>7}")
for label, curve, _K1h, m in hist:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    pred = predicted_k1_wall(n, m, K2)
    print(f"{label:<18} {lam_star(lam, n):>6.3f} {m:>3} {pred:>9.2f} "
          f"{str(observed_wall(wall_rows, label, 7)):>7} "
          f"{str(observed_wall(wall_rows, label, 8)):>7} "
          f"{str(observed_wall(wall_rows, label, 10)):>7}")


# ---------------------------------------------------------------------------
# E3 -- is rho a separator?  In-sample over the E1 grid.
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E3: rho as a success separator over the E1 grid")
print("-" * 78)

def auc(points):
    """points = [(score, label)]; AUC for 'low score => positive'."""
    pos = [s for s, y in points if y]
    neg = [s for s, y in points if not y]
    if not pos or not neg:
        return float('nan')
    wins = sum(1 for a in pos for b in neg if a < b) + \
           0.5 * sum(1 for a in pos for b in neg if a == b)
    return wins / (len(pos) * len(neg))

for name, col in (('V0', 7), ('V1', 8), ('V2-cvp', 10)):
    pts = [(r[6], r[col] >= 3) for r in wall_rows]
    hi_rho_win = [r for r in wall_rows if r[6] >= 1 and r[col] >= 3]
    lo_rho_fail = [r for r in wall_rows if r[6] < 1 and r[col] < 3]
    print(f"  {name:<7} AUC(rho) = {auc(pts):.3f}   "
          f"rho>=1 but wins: {len(hi_rho_win)}   rho<1 but fails: {len(lo_rho_fail)}")
    for r in hi_rho_win:
        print(f"      violation(win at rho>=1): {r[0]} K1={r[3]} rho={r[6]:.3f} w={r[col]}")
    for r in lo_rho_fail:
        print(f"      violation(fail at rho<1): {r[0]} K1={r[3]} rho={r[6]:.3f} w={r[col]}")


# ---------------------------------------------------------------------------
# E4 -- out-of-sample: fresh 17-bit curves across the lam* range
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E4: out-of-sample -- fresh 17-bit j=0 GLV curves, m=12")
print("-" * 78)

fresh = search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=10)
print(f"  found {len(fresh)} curves spanning lam* in "
      f"[{min(lam_star(c[3], c[2]) for c in fresh):.4f}, "
      f"{max(lam_star(c[3], c[2]) for c in fresh):.4f}]")

M_OOS = 12
OOS_SEEDS = SEEDS[:3]
oos_rows = []
print(f"\n  {'p':>7} {'n':>7} {'lam*':>6} {'K1':>3} {'eff':>6} {'rho':>6} "
      f"{'V0':>3} {'V1':>3} {'V2':>3}")
for (p, b, n, lam, G) in fresh:
    K2 = math.isqrt(n) + 1
    for K1 in (8, 16, 24, 32, 48, 64):     # spans rho ~ 0.67 .. 1.9 at m=12
        w0 = w1 = w2 = 0
        for seed in OOS_SEEDS:
            d_secret = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d_secret, M_OOS, n, lam, p, K1, K2, seed)
            if len(sigs) < M_OOS:
                continue
            ok0, _, _, _ = run_v0(sigs, n, lam, K1, K2, d_secret)
            ok1, _, _, _ = run_v1(sigs, n, lam, K1, K2, d_secret)
            out2, _, _ = run_v2(sigs, n, lam, K1, K2, d_secret)
            w0 += bool(ok0); w1 += bool(ok1); w2 += bool(out2['exact'][0])
        rho = rho_bdd(n, M_OOS, K1, K2)
        oos_rows.append((p, n, lam_star(lam, n), K1, K1 * K2 / n, rho, w0, w1, w2))
        print(f"  {p:>7} {n:>7} {lam_star(lam, n):>6.4f} {K1:>3} "
              f"{K1*K2/n:>6.3f} {rho:>6.3f} {w0:>3} {w1:>3} {w2:>3}")

print("\n  out-of-sample separator quality (majority = >=2/3):")
for name, col in (('V0', 6), ('V1', 7), ('V2', 8)):
    pts = [(r[5], r[col] >= 2) for r in oos_rows]
    lam_pts = [(r[2], r[col] >= 2) for r in oos_rows]
    eff_pts = [(r[4], r[col] >= 2) for r in oos_rows]
    print(f"    {name:<3} AUC(rho)={auc(pts):.3f}  AUC(eff)={auc(eff_pts):.3f}  "
          f"AUC(lam*)={auc(lam_pts):.3f}")


# ---------------------------------------------------------------------------
# E5 -- why do a few curves win far past rho = 1?
#       Replace the Gaussian heuristic by the *measured* lambda_1 of L'.
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E5: GH-based rho vs empirical rho = ||e|| / lambda_1(L')   (per-seed)")
print("-" * 78)

E5_SEEDS = [42, 1234, 9999, 555, 31337, 271828]
E5_K1 = (24, 32, 48)
e5_pts = []
print(f"  {'p':>7} {'lam*':>6} {'K1':>3} {'rho_GH':>7} {'rho_emp':>8} "
      f"{'l1/GH':>7} {'wins':>5}")
for (p, b, n, lam, G) in fresh:
    K2 = math.isqrt(n) + 1
    for K1 in E5_K1:
        wins = 0
        remp = []
        lratio = []
        for seed in E5_SEEDS:
            d_secret = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d_secret, M_OOS, n, lam, p, K1, K2, seed)
            if len(sigs) < M_OOS:
                continue
            ok0, _, _, _ = run_v0(sigs, n, lam, K1, K2, d_secret)
            Bp, _tp = build_v2(sigs, n, lam, K1, K2)
            basis = nonzero(lll_rows(Bp))
            l1 = min(norm(r) for r in basis)
            en = norm(error_vector_v2(sigs, n, K1, K2))
            g = gh(det_lprime(n, M_OOS, K1, K2), 2 * M_OOS)
            wins += bool(ok0)
            remp.append(en / l1)
            lratio.append(l1 / g)
            e5_pts.append((rho_bdd(n, M_OOS, K1, K2, en), en / l1, bool(ok0)))
        print(f"  {p:>7} {lam_star(lam, n):>6.4f} {K1:>3} "
              f"{rho_bdd(n, M_OOS, K1, K2):>7.3f} "
              f"{sum(remp)/len(remp):>8.3f} {sum(lratio)/len(lratio):>7.3f} "
              f"{wins:>2}/{len(E5_SEEDS)}")

print("\n  per-seed separator quality over %d points:" % len(e5_pts))
print(f"    AUC(rho_GH  = ||e||/GH(L'))       = {auc([(a, y) for a, _b, y in e5_pts]):.3f}")
print(f"    AUC(rho_emp = ||e||/lambda_1(L')) = {auc([(b, y) for _a, b, y in e5_pts]):.3f}")


# ---------------------------------------------------------------------------
# E6 -- the right BDD criterion is the *last* GS norm, not lambda_1.
#
# E5 shows lambda_1(L') << GH  =>  SUCCESS (anti-correlated).  Mechanism: det(L')
# is fixed, so an anomalously short b*_1 forces the remaining GS norms up, and
# nearest-plane succeeds whenever ||e|| < (1/2) min_i ||b*_i||.  Test that, and
# check where the short b*_1 comes from (the 2D lambda-block).
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E6: last-GS-norm BDD criterion, and the origin of the short lambda_1")
print("-" * 78)

from fpylll import GSO, FPLLL

def gs_norms(basis):
    """||b*_i|| for an LLL-reduced basis, computed at high precision."""
    A = IntegerMatrix.from_matrix(basis)
    for ft in ('mpfr', 'dd', 'd'):
        try:
            if ft == 'mpfr':
                FPLLL.set_precision(212)
            M = GSO.Mat(A, float_type=ft)
            M.update_gso()
            return [math.sqrt(M.get_r(i, i)) for i in range(A.nrows)]
        except Exception:
            continue
    return None

def gauss_reduce_2d(u, v):
    def n2(w): return w[0] * w[0] + w[1] * w[1]
    def dt(w, z): return w[0] * z[0] + w[1] * z[1]
    if n2(u) > n2(v): u, v = v, u
    while True:
        num, den = dt(v, u), n2(u)
        if den == 0: break
        q = (2*num + den) // (2*den) if num >= 0 else -((-2*num + den) // (2*den))
        v = (v[0] - q*u[0], v[1] - q*u[1])
        if n2(v) >= n2(u): break
        u, v = v, u
    return math.sqrt(n2(u))

def block_mu_ratio(n, lam, K1, K2):
    """|shortest vector| of the 2D lambda-block < (n*S_K1,0), (-lam*S_K1,S_K2) >,
    normalised by that block's own Gaussian heuristic sqrt(2/(2 pi e) * det)."""
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    mu = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return mu / (math.sqrt(2 / (2 * math.pi * math.e)) * math.sqrt(n * S_K1 * S_K2))

e6_pts = []
print(f"  {'p':>7} {'lam*':>6} {'mu/GHblk':>9} {'maxCF':>6} {'K1':>3} "
      f"{'l1/GH':>7} {'e/b*last':>9} {'e/babaiR':>9} {'wins':>5}")
for (p, b, n, lam, G) in fresh:
    K2 = math.isqrt(n) + 1
    a_cf, b_cf, cf = lam, n, []
    while b_cf:
        cf.append(a_cf // b_cf); a_cf, b_cf = b_cf, a_cf % b_cf
    max_cf = max(cf[1:]) if len(cf) > 1 else 0
    for K1 in E5_K1:
        wins = 0
        r_last, r_bab, r_l1 = [], [], []
        for seed in E5_SEEDS:
            d_secret = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d_secret, M_OOS, n, lam, p, K1, K2, seed)
            if len(sigs) < M_OOS:
                continue
            ok0, _, _, _ = run_v0(sigs, n, lam, K1, K2, d_secret)
            Bp, _tp = build_v2(sigs, n, lam, K1, K2)
            basis = nonzero(lll_rows(Bp))
            gsn = gs_norms(basis)
            if gsn is None:
                continue
            en = norm(error_vector_v2(sigs, n, K1, K2))
            g = gh(det_lprime(n, M_OOS, K1, K2), 2 * M_OOS)
            wins += bool(ok0)
            r_l1.append(min(norm(r) for r in basis) / g)
            r_last.append(en / gsn[-1])
            r_bab.append(en / (0.5 * min(gsn)))
            e6_pts.append((en / gsn[-1], en / (0.5 * min(gsn)), bool(ok0)))
        if not r_last:
            continue
        print(f"  {p:>7} {lam_star(lam, n):>6.4f} "
              f"{block_mu_ratio(n, lam, K1, K2):>9.3f} {max_cf:>6} {K1:>3} "
              f"{sum(r_l1)/len(r_l1):>7.3f} {sum(r_last)/len(r_last):>9.3f} "
              f"{sum(r_bab)/len(r_bab):>9.3f} {wins:>2}/{len(E5_SEEDS)}")

print("\n  per-seed separator quality over %d points:" % len(e6_pts))
print(f"    AUC(||e|| / ||b*_last||)          = "
      f"{auc([(a, y) for a, _b, y in e6_pts]):.3f}")
print(f"    AUC(||e|| / babai radius)         = "
      f"{auc([(b, y) for _a, b, y in e6_pts]):.3f}")
succ_last = [a for a, _b, y in e6_pts if y]
fail_last = [a for a, _b, y in e6_pts if not y]
if succ_last and fail_last:
    print(f"    ||e||/||b*_last||  successes: [{min(succ_last):.3f}, {max(succ_last):.3f}]"
          f"   failures: [{min(fail_last):.3f}, {max(fail_last):.3f}]")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
