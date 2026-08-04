"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, T5):
  In the Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` (dim 2m+2, columns
  [k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan]) the shortest vector after LLL
  is, on EVERY curve tested, the trivial vector  n*S_D*e_m  --- 100% of its
  energy in the d-column, |sv[m]|/n = 1.0000 exactly.  It is a lattice vector
  because  n*(d-row) - sum_i n*B_i*(mod-row i)  =  (0,...,0, n*S_D, 0,...,0),
  and it is always shorter than the planted vector:

      ||n*S_D*e_m||^2 = n^2 * S_D^2       (S_D = 1)
      ||v_planted||^2 ~ n^2 * (2m/3 + 4/3)

  It carries no information (d is only defined mod n; the vector IS the
  relation d -> d+n) and no choice of S_D removes it, since both norms scale
  linearly in S_D.  So recovery in that lattice is never an SVP condition.

THIS RUN'S FIX (the "projected" / d-free lattice).
  Quotient out the trivial direction by deleting the d-column outright.  Then
  the d-row (B_i*S_K1 | 0 | 0) becomes linearly DEPENDENT on the mod-rows
  (n*(d-row) is exactly 0), so the 2m+2 generators have rank only 2m+1 and
  fpylll's LLL would silently work on a rank-deficient basis.  Replace the
  m mod-rows + d-row by an explicit HNF basis of the rank-m sublattice they
  generate in the k1-columns:

      Lambda = { u in Z^m : exists x, u = x*B  (mod n) }
             = { u : u_i = u_0 * c_i (mod n) },   c_i = B_i * B_0^{-1} mod n

      HNF basis:  (1, c_1, ..., c_{m-1}) and n*e_i for i = 1..m-1
      det(Lambda) = n^(m-1)                       [index n inside n*Z^m]

  New lattice L' (dim 2m+1, columns [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan]):

      row 0      :  S_K1*(1, c_1, ..., c_{m-1})   |  0        |  0
      rows 1..m-1:  n*S_K1*e_i                    |  0        |  0
      rows m..2m-1: -lam*S_K1*e_j                 |  S_K2*e_j |  0
      row 2m     :  A_i*S_K1                      |  0        |  S_KANNAN

      det(L') = S_K1^m * S_K2^m * S_KANNAN * n^(m-1)  =  det(L_old) / n
      planted  = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN)          [membership: x = d]

  d is no longer read off a column; it is reconstructed from ONE recovered
  signature:   d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  mod n.

FALSIFIER (stated by the 2026-07-29 entry).
  If sv/pv rises above 1 after the reformulation AND the K1 wall of T4 moves
  outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
  real improvement.  If the wall stays at K1 ~ 4-6, the wall is
  information-theoretic and Phase 2 is at its ceiling.

Experiments:
  E1  construction check: rank, determinant vs the closed form, planted-vector
      membership, sv/pv for OLD vs NEW on the three historical curves.
  E2  the T4 K1-grid, OLD vs NEW lattice, identical signature sets, 5 seeds.
      This is the falsifier.
  E3  Gaussian-heuristic uniqueness ratio tau = ||v_planted|| / GH(L) for both
      lattices, and whether tau < 1 predicts success across a fresh sweep.
      Tests whether the K1 wall is the GH/uniqueness bound.

Run: python3 glv_hnp_phase2_projected.py
Deps: fpylll, sympy   (pip install fpylll cysignals sympy)
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + CM helpers (verbatim from
# glv_hnp_phase2_lambda_threshold.py so results are directly comparable)
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
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (norm form of Z[omega])."""
    if p % 3 != 1: return None
    # 4p = L^2 + 27 M^2  ->  Cornacchia on x^2 + 3 y^2 = 4p
    for L in range(1, int(math.isqrt(4 * p)) + 1):
        rem = 4 * p - L * L
        if rem <= 0: break
        if rem % 3: continue
        M2 = rem // 3
        M = math.isqrt(M2)
        if M * M != M2: continue
        # a + b*omega with a = (L + M)/2 ... enumerate the consistent pair
        for (a, b) in ((( L + 3 * M) // 2, M), ((L - 3 * M) // 2, -M),
                       ((L + M) // 2, M), ((L - M) // 2, -M)):
            if a * a - a * b + b * b == p:
                return (a, b)
    return None

def j0_traces(a, b):
    """Traces of the six sextic twists of a j=0 curve with pi = a + b*omega."""
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, or None."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
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
        if y is None or y == 0: continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

# ---------------------------------------------------------------------------
# Signatures
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
# OLD lattice (dim 2m+2, with the d-column) -- verbatim from prior runs
# ---------------------------------------------------------------------------

def build_old(sigs, n, lam, k1_bound, k2_bound):
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

def planted_old(sigs, d, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_old(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sgn = 1 if last > 0 else -1
        d_cand = (sgn * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# NEW lattice (dim 2m+1, d-column projected out)
# ---------------------------------------------------------------------------

def build_new(sigs, n, lam, k1_bound, k2_bound):
    """Columns: [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan].  Full rank 2m+1."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    c = [sigs[i]['B'] * B0inv % n for i in range(m)]   # c[0] == 1
    assert c[0] == 1

    M = [[0] * dim for _ in range(dim)]
    # HNF basis of Lambda = {u : u_i = u_0*c_i mod n}, scaled by S_K1
    for i in range(m):
        M[0][i] = c[i] * S_K1
    for i in range(1, m):
        M[i][i] = n * S_K1
    # k2 rows
    for j in range(m):
        M[m + j][j] = -lam * S_K1
        M[m + j][m + j] = S_K2
    # Kannan row
    for i in range(m):
        M[2 * m][i] = sigs[i]['A'] * S_K1
    M[2 * m][dim - 1] = S_KAN
    return M

def planted_new(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_new(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read k1_0, k2_0 off a Kannan row and reconstruct d algebraically."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sgn = 1 if last > 0 else -1
        a0 = sgn * row[0]
        b0 = sgn * row[m]
        if a0 % S_K1 or b0 % S_K2: continue
        k1_0, k2_0 = a0 // S_K1, b0 // S_K2
        d_cand = (k1_0 + lam * k2_0 - sigs[0]['A']) * B0inv % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Membership / determinant checks
# ---------------------------------------------------------------------------

def in_lattice(M, v):
    """Is v in the row-span (over Z) of M?  Solve v = x*M by HNF."""
    from sympy import Matrix
    A = Matrix(M).T                      # columns are the basis vectors
    b = Matrix(len(v), 1, list(v))
    sol = A.gauss_jordan_solve(b) if A.rank() == A.cols else None
    if sol is None: return False
    x = sol[0]
    return all(val.is_integer for val in x)

def det_lattice(M):
    from sympy import Matrix
    return abs(Matrix(M).det())

def gaussian_heuristic(det, dim):
    return math.sqrt(dim / (2 * math.pi * math.e)) * (float(det) ** (1.0 / dim))

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_once(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if variant == 'old':
        M = build_old(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_old(sigs, d_secret, n, k1_bound, k2_bound)
    else:
        M = build_new(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_new(sigs, n, k1_bound, k2_bound)
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    if variant == 'old':
        ok = recover_old(rows, m, n, S_KAN, d_secret) is not None
    else:
        ok = recover_new(rows, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
    nz = [r for r in rows if any(r)]
    sv = min(norm(r) for r in nz) if nz else float('inf')
    return {'ok': ok, 'sv': sv, 'pv': norm(pv), 'rows': rows, 'sigs': sigs,
            'M': M, 'dim': dim}

def rate(curve, m, d_secret, k1_bound, seeds, variant, use_bkz=False, beta=20):
    good = 0
    svpv = []
    for s in seeds:
        r = run_once(curve, m, d_secret, k1_bound, s, variant, use_bkz, beta)
        if r is None: continue
        good += 1 if r['ok'] else 0
        svpv.append(r['sv'] / r['pv'])
    return good, (sum(svpv) / len(svpv) if svpv else float('nan'))


SEEDS = [42, 1234, 9999, 555, 31337]
D_SECRET_FRAC = 0.37     # d = floor(0.37*n); fixed across runs for comparability

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26 / 2026-07-29)
HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 23 - project out the d-column so the planted vector is lambda_1")
    print("=" * 78)

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), k1, m))
    print(f"[setup] {len(hist)} historical curves rebuilt, generators verified.")

    # ===========================================================================
    # E1 - construction check
    # ===========================================================================
    print()
    print("=" * 78)
    print("E1  construction check: rank, det, planted membership, sv/pv")
    print("=" * 78)

    print(f"{'curve':<18} {'m':>2} {'K1':>3} | {'dim':>4} {'det==closed':>12} "
          f"{'planted in L':>13} | {'sv/pv':>7} {'|sv[d]|/n':>10}")
    print("-" * 78)

    e1_rows = []
    for label, curve, k1_bound, m in hist:
        p, b, n, lam, G = curve
        d_secret = max(1, int(D_SECRET_FRAC * n))
        k2_bound = math.isqrt(n) + 1
        S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, 42)

        for variant in ('old', 'new'):
            if variant == 'old':
                M = build_old(sigs, n, lam, k1_bound, k2_bound)
                pv = planted_old(sigs, d_secret, n, k1_bound, k2_bound)
                closed = (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN
            else:
                M = build_new(sigs, n, lam, k1_bound, k2_bound)
                pv = planted_new(sigs, n, k1_bound, k2_bound)
                closed = S_K1 ** m * S_K2 ** m * S_KAN * n ** (m - 1)
            dim = len(M)
            det = det_lattice(M)
            det_ok = (det == closed)
            member = in_lattice(M, pv)

            A = IntegerMatrix.from_matrix(M)
            LLL.reduction(A)
            rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
            nz = [r for r in rows if any(r)]
            sv_row = min(nz, key=norm)
            sv = norm(sv_row)
            # energy of the shortest vector in the d-column (old lattice only)
            dcol = abs(sv_row[m]) / n if variant == 'old' else float('nan')

            print(f"{label:<18} {m:>2} {k1_bound:>3} | {dim:>4} {str(det_ok):>12} "
                  f"{str(member):>13} | {sv/norm(pv):>7.3f} "
                  f"{(f'{dcol:.4f}' if variant == 'old' else '   n/a'):>10}   [{variant}]")
            e1_rows.append((label, variant, dim, det_ok, member, sv / norm(pv)))

    # ===========================================================================
    # E2 - the falsifier: does the K1 wall move?
    # ===========================================================================
    print()
    print("=" * 78)
    print("E2  FALSIFIER - T4 K1 grid, OLD vs NEW lattice, identical signature sets")
    print("=" * 78)
    print("(cell = #seeds recovering d, out of 5)")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
    E2_CURVES = [(label, curve, m) for (label, curve, k1, m) in hist if '12-bit' in label]

    for label, curve, m in E2_CURVES:
        p, b, n, lam, G = curve
        d_secret = max(1, int(D_SECRET_FRAC * n))
        k2_bound = math.isqrt(n) + 1
        ls = lam_star(lam, n)
        print()
        print(f"  {label}   n={n}  lam*={ls:.4f}  m={m}  K2={k2_bound}")
        hdr = "    " + "".join(f"{k:>6}" for k in K1_GRID)
        print(f"    {'K1':<0}" + "".join(f"{k:>6}" for k in K1_GRID)[4:])
        print("    " + "-" * (6 * len(K1_GRID)))
        for variant in ('old', 'new'):
            cells = []
            for k1 in K1_GRID:
                g, _ = rate(curve, m, d_secret, k1, SEEDS, variant)
                cells.append(g)
            print(f"    {variant:<4}" + "".join(f"{c:>6}" for c in cells))
        effs = [k1 * k2_bound / n for k1 in K1_GRID]
        print(f"    eff " + "".join(f"{e:>6.2f}" for e in effs))
        print(f"    (1/m = {1.0/m:.3f})")

    # ===========================================================================
    # E3 - is the wall the Gaussian-heuristic uniqueness bound?
    # ===========================================================================
    print()
    print("=" * 78)
    print("E3  uniqueness ratio  tau = ||v_planted|| / GH(L)   vs observed success")
    print("=" * 78)
    print("If tau < 1 the planted vector is (heuristically) the unique shortest")
    print("vector; tau > 1 means it is not expected to be findable by ANY lattice")
    print("algorithm, i.e. the wall is information-theoretic, not algorithmic.")
    print()
    print(f"{'curve':<16} {'K1':>3} {'eff':>6} | {'tau_old':>8} {'tau_new':>8} | "
          f"{'old':>5} {'new':>5} {'newBKZ':>7}")
    print("-" * 78)

    e3 = []
    for label, curve, m in E2_CURVES:
        p, b, n, lam, G = curve
        d_secret = max(1, int(D_SECRET_FRAC * n))
        k2_bound = math.isqrt(n) + 1
        S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound=1, k2_bound=k2_bound)  # placeholder
        for k1 in K1_GRID:
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2_bound)
            det_old = (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN
            det_new = S_K1 ** m * S_K2 ** m * S_KAN * n ** (m - 1)
            gh_old = gaussian_heuristic(det_old, 2 * m + 2)
            gh_new = gaussian_heuristic(det_new, 2 * m + 1)
            # expected planted norm: k1_i*S_K1 ~ U[0,n), k2_i*S_K2 ~ U[0,n)
            pv_exp_old = math.sqrt(2 * m * n * n / 3.0 + (D_SECRET_FRAC * n) ** 2 + n * n)
            pv_exp_new = math.sqrt(2 * m * n * n / 3.0 + n * n)
            g_old, _ = rate(curve, m, d_secret, k1, SEEDS, 'old')
            g_new, _ = rate(curve, m, d_secret, k1, SEEDS, 'new')
            g_bkz, _ = rate(curve, m, d_secret, k1, SEEDS, 'new', use_bkz=True, beta=20)
            eff = k1 * k2_bound / n
            t_old, t_new = pv_exp_old / gh_old, pv_exp_new / gh_new
            print(f"{label:<16} {k1:>3} {eff:>6.2f} | {t_old:>8.3f} {t_new:>8.3f} | "
                  f"{g_old:>5} {g_new:>5} {g_bkz:>7}")
            e3.append((label, k1, eff, t_old, t_new, g_old, g_new, g_bkz))
        print("-" * 78)

    # separation quality of tau_new
    succ = [r for r in e3 if r[6] >= 3]
    fail = [r for r in e3 if r[6] < 3]
    if succ and fail:
        print(f"tau_new over successes (>=3/5): [{min(r[4] for r in succ):.3f}, "
              f"{max(r[4] for r in succ):.3f}]")
        print(f"tau_new over failures  (<3/5) : [{min(r[4] for r in fail):.3f}, "
              f"{max(r[4] for r in fail):.3f}]")
        sep = max(r[4] for r in succ) < min(r[4] for r in fail)
        print(f"tau_new separates success from failure cleanly: {sep}")
        print(f"eff over successes: [{min(r[2] for r in succ):.3f}, "
              f"{max(r[2] for r in succ):.3f}]  "
              f"failures: [{min(r[2] for r in fail):.3f}, {max(r[2] for r in fail):.3f}]")

    print()
    print("=" * 78)
    print("done")
    print("=" * 78)
