"""
GLV-HNP Phase 2 / Thread 23: is the K1 wall a reduction wall or an
information-theoretic wall?

Background (RESEARCH_AUTOLAB_LOG.md 2026-07-29, T5): in the (2m+2)-dimensional
Kannan lattice of `glv_hnp_phase2_20bit.py:262`, the planted vector is NEVER the
shortest vector -- the trivial vector n*S_D*e_m (norm exactly n) is shorter for
every m >= 1, and no choice of S_D removes it.  Recovery is therefore a BDD /
coset condition, not an SVP condition.

This script tests the reformulation proposed there, plus the decisive oracle
experiment:

  A0  Kannan baseline   dim 2m+2, verbatim construction, LLL/BKZ.
  A1  Projected         delete column m (the d-column).  pi is a lattice
                        homomorphism, pi(n*e_m) = 0, so the trivial vector is
                        killed exactly.  d is recovered afterwards from a single
                        signature: d = B_0^{-1} (k1_0 + lam*k2_0 - A_0) mod n.
  A2  Explicit CVP      drop the d-column AND the Kannan row/column; Babai
                        nearest-plane (or exact CVP) against the target
                        t = (-A_i*S_K1, 0^m).

  A3  List-decoded CVP  the A2 lattice, but enumerating EVERY lattice point
                        within a public-parameter radius of t and accepting if
                        any decodes to the true d.  Each candidate costs one
                        scalar multiplication to test against the public key,
                        so this is a real attack, not an oracle.

Headline result (EXP 8): the "K1 wall" reported on 2026-07-26 and re-measured
on 2026-07-29 is an artifact of NEAREST-POINT decoding, not a property of the
lattice.  On the canonical failure curve (12-bit/2677, m=10) LLL and BKZ-40
both give 0/10 at K1=8 and exact nearest-point CVP gives 1/10, but list
decoding gives 10/10 -- and still 5/10 at K1=12, where every nearest-point arm
is at 0/10.

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random

from fpylll import (IntegerMatrix, LLL, BKZ, GSO, CVP,
                    Enumeration, EvaluatorStrategy)

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

# ---------------------------------------------------------------------------
# Phase-2 lattice, verbatim
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, kan_mult=1.0):
    """Column scales.  kan_mult=1.0 reproduces the historical S_KANNAN = n.

    The textbook Kannan embedding wants the embedding factor to track the
    target distance, which here is ~ n*sqrt(2m/3), not n; kan_mult exposes it.
    """
    return (n // k1_bound, 1, max(1, n // k2_bound),
            max(1, int(round(n * kan_mult))))

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
    """A0: the (2m+2)-dim Kannan lattice, verbatim."""
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

# ---------------------------------------------------------------------------
# A1: projected lattice (d-column deleted)
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound, kan_mult=1.0):
    """Delete column m from the A0 basis.

    Column layout (dim 2m+1):
        0 .. m-1     k1 columns   (scale S_K1)
        m .. 2m-1    k2 columns   (scale S_K2)
        2m           Kannan column (scale S_KANNAN)

    The basis has 2m+2 rows spanning a rank-(2m+1) lattice: the rank drop is
    exactly the kernel n*e_m of the projection.  fplll's LLL handles the
    dependency and emits a zero row.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound, kan_mult)
    M = []
    for i in range(m):                       # modular rows
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    row = [0] * dim                          # B-row (the old d-row, d-col gone)
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for i in range(m):                       # lambda rows
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    row = [0] * dim                          # Kannan / A-row
    for i in range(m):
        row[i] = sigs[i]['A'] * S_K1
    row[dim - 1] = S_KANNAN
    M.append(row)
    return M

def projected_planted_vector(sigs, n, k1_bound, k2_bound, kan_mult=1.0):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, kan_mult)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def d_from_k(sigs, k1_0, k2_0, n, lam):
    """d = B_0^{-1} (k1_0 + lam*k2_0 - A_0) mod n."""
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    if B0 % n == 0:
        return None
    return (k1_0 + lam * k2_0 - A0) * modinv(B0, n) % n

def recover_d_projected(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret,
                        kan_mult=1.0):
    """Scan for a row with |last| = S_KANNAN, read (k1, k2), rebuild d."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, kan_mult)
    for row in reduced:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sgn = 1 if row[dim - 1] > 0 else -1
        if sgn * row[0] % S_K1 != 0 or sgn * row[m] % S_K2 != 0:
            continue
        k1_0 = sgn * row[0] // S_K1
        k2_0 = sgn * row[m] // S_K2
        d_cand = d_from_k(sigs, k1_0, k2_0, n, lam)
        if d_cand is not None and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# A2: explicit CVP (no Kannan embedding, no d-column)
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """dim 2m.  Rows: m modular, 1 B-row, m lambda rows (rank 2m)."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    row = [0] * dim
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    return M

def cvp_target(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, _S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -(sigs[i]['A'] % n) * S_K1
    return t

def cvp_planted_point(sigs, n, k1_bound, k2_bound):
    """The lattice point w with t - w = (-k1_i*S_K1, -k2_i*S_K2).

    w = sum_i c_i (n*S_K1*e_i) + d*(B-row) + sum_i k2_i*(lambda-row), so the
    k2-block of w carries +k2_i*S_K2 (the lambda-row has +S_K2 there).
    """
    m = len(sigs)
    t = cvp_target(sigs, n, k1_bound, k2_bound)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    w = [0] * (2 * m)
    for i in range(m):
        w[i] = t[i] + sigs[i]['k1'] * S_K1
        w[m + i] = sigs[i]['k2'] * S_K2
    return w

# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def dist(u, v):
    return math.sqrt(sum((float(a) - float(b)) ** 2 for a, b in zip(u, v)))

def to_rows(A):
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def nonzero_rows(rows):
    return [r for r in rows if any(x != 0 for x in r)]

def reduce_rows(M, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return to_rows(A)

def in_lattice(rows, v):
    """Membership test for v in the lattice spanned by `rows`, via HNF."""
    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form
    B = Matrix([r for r in rows if any(x != 0 for x in r)])
    H = hermite_normal_form(B.T).T            # row-style HNF
    # solve x*H = v over the rationals, require integrality
    sol = H.T.solve_least_squares(Matrix(v))
    if any(x != int(x) for x in sol):
        return False
    return list(H.T * sol) == list(v)

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_A0(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    rows = reduce_rows(M, use_bkz, beta)
    dim = 2 * m + 2
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    ok = False
    for row in rows:
        if abs(row[dim - 1]) != S_KAN: continue
        sgn = 1 if row[dim - 1] > 0 else -1
        if sgn * row[m] % n == d_secret:
            ok = True
    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    sv = min(norm(r) for r in nonzero_rows(rows))
    return {'ok': ok, 'pn': norm(pv), 'sn': sv, 'sigs': sigs}

def run_A1(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20,
           kan_mult=1.0):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound, kan_mult)
    rows = reduce_rows(M, use_bkz, beta)
    ok = recover_d_projected(rows, sigs, n, lam, k1_bound, k2_bound,
                             d_secret, kan_mult) is not None
    pv = projected_planted_vector(sigs, n, k1_bound, k2_bound, kan_mult)
    sv = min(norm(r) for r in nonzero_rows(rows))
    return {'ok': ok, 'pn': norm(pv), 'sn': sv, 'sigs': sigs}

def run_A1_svp(curve, m, d_secret, k1_bound, seed, kan_mult=1.0):
    """A1 lattice, but solved with EXACT SVP instead of LLL/BKZ-40.

    EXP 1 shows pi(v_p) is lambda_1 in the small-K1 regime; this arm asks
    whether it stays lambda_1 past the LLL wall, i.e. whether the wall is a
    strength-of-reduction issue in the projected formulation.

    Exact SVP is obtained as BKZ with block_size = dim (HKZ reduction), whose
    first basis vector is provably a shortest vector.  SVP.shortest_vector()
    is not usable here: it loads fplll's pruning-strategy JSON, which a
    pip-installed fpylll does not ship.
    """
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound, kan_mult)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    A = IntegerMatrix.from_matrix(nonzero_rows(to_rows(A)))
    LLL.reduction(A)
    BKZ.reduction(A, BKZ.Param(A.nrows))          # block_size = dim  =>  HKZ
    rows = nonzero_rows(to_rows(A))
    v = min(rows, key=norm)
    ok = recover_d_projected([v, [-x for x in v]], sigs, n, lam, k1_bound,
                             k2_bound, d_secret, kan_mult) is not None
    pv = projected_planted_vector(sigs, n, k1_bound, k2_bound, kan_mult)
    return {'ok': ok, 'pn': norm(pv), 'sn': norm(v)}

def run_A2(curve, m, d_secret, k1_bound, seed, exact=False):
    """Babai nearest-plane (exact=False) or exact CVP (exact=True)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = nonzero_rows(to_rows(A))
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    t = cvp_target(sigs, n, k1_bound, k2_bound)
    if exact:
        w = list(CVP.closest_vector(A, tuple(int(x) for x in t)))
    else:
        g = GSO.Mat(A, float_type='mpfr' if n.bit_length() > 200 else 'd')
        g.update_gso()
        coeffs = g.babai([int(x) for x in t])
        w = [0] * A.ncols
        for i, c in enumerate(coeffs):
            if c:
                for j in range(A.ncols):
                    w[j] += c * A[i][j]
    S_K1, _S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    ok = False
    if (w[0] - t[0]) % S_K1 == 0 and w[m] % S_K2 == 0:
        k1_0 = (w[0] - t[0]) // S_K1
        k2_0 = w[m] // S_K2
        d_cand = d_from_k(sigs, k1_0, k2_0, n, lam)
        ok = (d_cand == d_secret)
    wp = cvp_planted_point(sigs, n, k1_bound, k2_bound)
    return {'ok': ok, 'd_planted': dist(t, wp), 'd_found': dist(t, w),
            'same_point': list(w) == list(wp),
            'sigs': sigs, 'w': w, 'wp': wp}

def run_A3_list(curve, m, d_secret, k1_bound, seed, radius_mult=1.10,
                nr_solutions=2000):
    """A2 lattice, LIST-decoded: enumerate every lattice point within radius R
    of t and accept if ANY of them decodes to the true d.

    The radius uses PUBLIC parameters only -- E[dist(t, w_planted)]^2 =
    m*(K1*S_K1)^2/3 + m*(K2*S_K2)^2/3 -- so this is a real attack, not an
    oracle.  Nearest-point CVP is the special case R = dist(t, CVP(t)); the
    Kannan arms of EXP 7 already beat nearest-point decoding, which is why the
    A2-EXACT column is NOT the information-theoretic ceiling.
    """
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    A = IntegerMatrix.from_matrix(nonzero_rows(to_rows(A)))
    LLL.reduction(A)
    t = cvp_target(sigs, n, k1_bound, k2_bound)
    S_K1, _S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    r_pub = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                      + m * (k2_bound * S_K2) ** 2 / 3.0)
    R = radius_mult * r_pub
    gso = GSO.Mat(A, update=True)
    tc = gso.from_canonical([float(x) for x in t])
    enum = Enumeration(gso, nr_solutions=nr_solutions,
                       strategy=EvaluatorStrategy.BEST_N_SOLUTIONS)
    try:
        sols = enum.enumerate(0, A.nrows, R * R, 0, target=tc)
    except Exception:                                # noqa: BLE001  (no point in radius)
        sols = []
    ok = False
    for _length, coeffs in sols:
        w = [0] * A.ncols
        for i, c in enumerate(coeffs):
            c = int(round(c))
            if not c: continue
            for j in range(A.ncols):
                w[j] += c * A[i][j]
        if (w[0] - t[0]) % S_K1 or w[m] % S_K2:
            continue
        if d_from_k(sigs, (w[0] - t[0]) // S_K1, w[m] // S_K2, n,
                    lam) == d_secret:
            ok = True
            break
    return {'ok': ok, 'n_cand': len(sols)}

def rate(fn, curve, m, k1_bound, seeds, **kw):
    p, b, n, lam, G = curve
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = fn(curve, m, d_trial, k1_bound, seed, **kw)
        if r is None: continue
        wins += bool(r['ok'])
        if 'pn' in r and r['pn']:
            ratios.append(r['sn'] / r['pn'])
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr

# ---------------------------------------------------------------------------
# Curves (identical to the 2026-07-29 run so the grids are comparable)
# ---------------------------------------------------------------------------

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
SEEDS = [42, 1234, 9999, 555, 31337]

print("=" * 78)
print("Thread 23 -- projected lattice / explicit CVP  (GLV-HNP Phase 2)")
print("=" * 78)

curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP 1: does the projection kill the trivial vector, and is pi(v_p) now lambda_1?")
print("-" * 78)
print("A0 = Kannan (2m+2)-dim; A1 = projected (2m+1)-dim.  sv/pv < 1 means the")
print("planted vector is NOT the shortest vector after reduction.\n")
print(f"{'curve':<20}{'K1':>4}{'m':>4}  {'A0 sv/pv':>9}{'A1 sv/pv':>10}"
      f"{'A0 |sv|/n':>11}{'A1 |sv|/n':>11}  {'A0 ok':>6}{'A1 ok':>6}")

exp1_rows = []
for label, curve, k1, m in curves:
    p, b, n, lam, G = curve
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    r0 = run_A0(curve, m, d_trial, k1, 42)
    r1 = run_A1(curve, m, d_trial, k1, 42)
    exp1_rows.append((label, k1, m, r0, r1))
    print(f"{label:<20}{k1:>4}{m:>4}  {r0['sn']/r0['pn']:>9.3f}"
          f"{r1['sn']/r1['pn']:>10.3f}{r0['sn']/n:>11.4f}{r1['sn']/n:>11.4f}"
          f"  {str(r0['ok']):>6}{str(r1['ok']):>6}")

print("\nMembership check: is pi(v_planted) actually in the projected lattice?")
for label, curve, k1, m in curves:
    p, b, n, lam, G = curve
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_trial, m, n, lam, p, k1, k2b, 42)
    M = build_projected_lattice(sigs, n, lam, k1, k2b)
    v = projected_planted_vector(sigs, n, k1, k2b)
    try:
        member = in_lattice(M, v)
    except Exception as e:                       # noqa: BLE001
        member = f"error: {e}"
    print(f"  {label:<20} pi(v_p) in L~ : {member}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP 2: the K1 wall under all three formulations (T4 grid replication)")
print("-" * 78)
print("2026-07-29 T4 measured the A0 column only.  If A1/A2 push the wall")
print("outward, the reformulation is a real improvement.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
exp2 = {}
for label, curve, _k1, m in curves[1:]:          # the two 12-bit curves
    p, b, n, lam, G = curve
    print(f"{label}  (lam*={min(lam, n-lam)/n:.3f}, n={n}, m={m}, "
          f"{len(SEEDS)} seeds)")
    print(f"  {'K1':>4}  {'A0':>6}{'A1':>6}{'A2':>6}   {'A0 sv/pv':>9}{'A1 sv/pv':>9}")
    exp2[label] = {}
    for k1 in K1_GRID:
        w0, t0, r0 = rate(run_A0, curve, m, k1, SEEDS)
        w1, t1, r1 = rate(run_A1, curve, m, k1, SEEDS)
        w2, t2, _ = rate(run_A2, curve, m, k1, SEEDS)
        exp2[label][k1] = (w0, w1, w2, r0, r1)
        print(f"  {k1:>4}  {w0:>3}/{t0}{w1:>4}/{t1}{w2:>4}/{t2}   "
              f"{r0:>9.3f}{r1:>9.3f}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP 3: exact CVP oracle -- is the planted point even the closest one?")
print("-" * 78)
print("If dist(t, CVP(t)) < dist(t, w_planted) then a closer lattice point")
print("exists and NO algorithm reading only this lattice can return the true d.")
print("That makes the wall information-theoretic, not a reduction failure.\n")

# keep dimensions tractable for exact enumeration: dim = 2m
E3 = [
    ("12-bit/2557", 6, [2, 4, 8, 12, 16]),
    ("12-bit/2677 FAIL", 6, [2, 4, 6, 8, 12]),
    ("12-bit/2677 FAIL", 8, [4, 6, 8]),
]
lookup = {label: (curve, k1, m) for label, curve, k1, m in curves}
print(f"{'curve':<20}{'m':>3}{'K1':>4}  {'seed':>6}  {'d(t,w_p)':>12}"
      f"{'d(t,CVP)':>12}  {'ratio':>7}  {'CVP=planted':>12}")
exp3_rows = []
for label, m, k1s in E3:
    curve = lookup[label][0]
    p, b, n, lam, G = curve
    for k1 in k1s:
        agree = 0
        closer = 0
        for seed in SEEDS[:3]:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = run_A2(curve, m, d_trial, k1, seed, exact=True)
            agree += bool(r['ok'])
            closer += (r['d_found'] < r['d_planted'] - 1e-6)
            exp3_rows.append((label, m, k1, seed, r['d_planted'], r['d_found'],
                              r['ok']))
            print(f"{label:<20}{m:>3}{k1:>4}  {seed:>6}  {r['d_planted']:>12.1f}"
                  f"{r['d_found']:>12.1f}  {r['d_found']/r['d_planted']:>7.3f}"
                  f"  {str(r['same_point']):>12}")
        print(f"{'':<20}{m:>3}{k1:>4}  -> exact-CVP recovers d: {agree}/3; "
              f"strictly closer non-planted point: {closer}/3")
print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP 4: exact CVP at the *same* (m, K1) where EXP 2's LLL wall sits")
print("-" * 78)
print("EXP 3 used m=6/8 for enumeration speed; the EXP 2 walls were measured at")
print("m=8 (2557) and m=10 (2677).  This closes the gap: exact CVP (dim 2m, up")
print("to 20) run exactly where LLL fails.\n")

E4 = [
    ("12-bit/2557", 8, K1_GRID),
    ("12-bit/2677 FAIL", 10, K1_GRID),
]
print("Columns: LLL on A0 / LLL on A1 / BKZ-20 on A1 / BKZ-40 on A1 /")
print("Babai on A2 / EXACT CVP on A2.  The exact-CVP column is the")
print("information-theoretic boundary: it is what an unbounded-time algorithm")
print("reading this lattice can do.\n")
walls4 = {}
for label, m, k1s in E4:
    curve = lookup[label][0]
    p, b, n, lam, G = curve
    print(f"{label}  (m={m}, dim(A0)={2*m+2}, dim(A2)={2*m}, {len(SEEDS)} seeds)")
    print(f"  {'K1':>4}  {'A0-LLL':>7}{'A1-LLL':>7}{'A1-BKZ20':>9}{'A1-BKZ40':>9}"
          f"{'A2-Babai':>9}{'A2-EXACT':>9}   {'closer pt':>10}{'mean ratio':>11}")
    walls4[(label, m)] = {}
    for k1 in k1s:
        w_a0, _, _ = rate(run_A0, curve, m, k1, SEEDS)
        w_a1, _, _ = rate(run_A1, curve, m, k1, SEEDS)
        w_b20, _, _ = rate(run_A1, curve, m, k1, SEEDS, use_bkz=True, beta=20)
        w_b40, _, _ = rate(run_A1, curve, m, k1, SEEDS, use_bkz=True, beta=40)
        w_bab, _, _ = rate(run_A2, curve, m, k1, SEEDS)
        exact = closer = 0
        ratios = []
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = run_A2(curve, m, d_trial, k1, seed, exact=True)
            exact += bool(r['ok'])
            closer += (r['d_found'] < r['d_planted'] - 1e-6)
            ratios.append(r['d_found'] / r['d_planted'])
        walls4[(label, m)][k1] = (w_a0, w_a1, w_b20, w_b40, w_bab, exact)
        ns = len(SEEDS)
        print(f"  {k1:>4}  {w_a0:>5}/{ns}{w_a1:>5}/{ns}{w_b20:>7}/{ns}"
              f"{w_b40:>7}/{ns}{w_bab:>7}/{ns}{exact:>7}/{ns}   "
              f"{closer:>8}/{ns}{sum(ratios)/len(ratios):>11.3f}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP 5: when a closer point exists, does it decode to a DIFFERENT d?")
print("-" * 78)
print("d fixes k1_i + lam*k2_i mod n but not the split, so a closer lattice")
print("point can still carry the true d.  The attack only fails when the closer")
print("point implies a different d.\n")
print(f"{'curve':<20}{'m':>3}{'K1':>4}  {'seed':>6}  {'closer':>7}"
      f"  {'d ok':>5}  {'implied d':>12}  {'true d':>12}")
for label, m, k1s in [("12-bit/2677 FAIL", 6, [4, 6, 8, 12])]:
    curve = lookup[label][0]
    p, b, n, lam, G = curve
    for k1 in k1s:
        for seed in SEEDS[:3]:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = run_A2(curve, m, d_trial, k1, seed, exact=True)
            S_K1, _S, S_K2, _K = scales(n, k1, math.isqrt(n) + 1)
            t = cvp_target(r['sigs'], n, k1, math.isqrt(n) + 1)
            imp = None
            if (r['w'][0] - t[0]) % S_K1 == 0 and r['w'][m] % S_K2 == 0:
                imp = d_from_k(r['sigs'], (r['w'][0] - t[0]) // S_K1,
                               r['w'][m] // S_K2, n, lam)
            print(f"{label:<20}{m:>3}{k1:>4}  {seed:>6}"
                  f"  {str(r['d_found'] < r['d_planted'] - 1e-6):>7}"
                  f"  {str(r['ok']):>5}  {str(imp):>12}  {d_trial:>12}")
print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP 6: the gap point -- 12-bit/2677, m=10, 10 seeds")
print("-" * 78)
print("EXP 4 shows LLL/BKZ-40 at 0-1/5 where exact CVP is at 4/5.  Widen the")
print("seed set and add EXACT SVP on the projected lattice to decide whether")
print("the gap is real and whether the projection makes it reachable.\n")

SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 2024, 60613, 8191, 271828]
curve2677 = lookup["12-bit/2677 FAIL"][0]
n2677 = curve2677[2]
print(f"  {'K1':>4}  {'A0-LLL':>7}{'A1-LLL':>7}{'A1-BKZ40':>9}{'A1-SVP*':>9}"
      f"{'A2-Babai':>9}{'A2-CVP*':>9}      (* = exact)")
for k1 in [3, 4, 6, 8]:
    r_a0, _, _ = rate(run_A0, curve2677, 10, k1, SEEDS10)
    r_a1, _, _ = rate(run_A1, curve2677, 10, k1, SEEDS10)
    r_b4, _, _ = rate(run_A1, curve2677, 10, k1, SEEDS10, use_bkz=True, beta=40)
    r_sv, _, _ = rate(run_A1_svp, curve2677, 10, k1, SEEDS10)
    r_bb, _, _ = rate(run_A2, curve2677, 10, k1, SEEDS10)
    r_cv, _, _ = rate(run_A2, curve2677, 10, k1, SEEDS10, exact=True)
    ns = len(SEEDS10)
    print(f"  {k1:>4}  {r_a0:>5}/{ns}{r_a1:>5}/{ns}{r_b4:>7}/{ns}"
          f"{r_sv:>7}/{ns}{r_bb:>7}/{ns}{r_cv:>7}/{ns}")
print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP 7: is S_KANNAN = n mistuned?")
print("-" * 78)
print("EXP 6 shows exact SVP on the embedded lattice (1/10) losing to exact CVP")
print("on the un-embedded one (5/10) at K1=6 -- the embedding itself is leaking.")
print("Standard Kannan wants the embedding factor ~ the target distance, which")
print("here is ||(k1*S_K1, k2*S_K2)|| ~ n*sqrt(2m/3), not n.\n")

for m_e in (10,):
    ideal = math.sqrt(2 * m_e / 3.0)
    print(f"  m={m_e}: n*sqrt(2m/3) => kan_mult_ideal = {ideal:.2f}")
print(f"  {'kan_mult':>9}", end="")
for k1 in (4, 6, 8):
    print(f"{'K1=' + str(k1) + ' LLL':>12}{'BKZ40':>8}", end="")
print()
for km in (0.25, 0.5, 1.0, 2.0, 2.58, 4.0, 8.0):
    print(f"  {km:>9.2f}", end="")
    for k1 in (4, 6, 8):
        w_l, _, _ = rate(run_A1, curve2677, 10, k1, SEEDS10, kan_mult=km)
        w_b, _, _ = rate(run_A1, curve2677, 10, k1, SEEDS10, kan_mult=km,
                         use_bkz=True, beta=40)
        print(f"{w_l:>9}/10{w_b:>6}/10", end="")
    print()
print(f"\n  (baseline for reference: exact CVP = 10/10, 5/10, 1/10 at "
      f"K1 = 4, 6, 8)")
print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP 8: list-decoding -- the actual information-theoretic ceiling")
print("-" * 78)
print("EXP 7 shows a tuned Kannan arm (3/10 at K1=8) beating exact nearest-point")
print("CVP (1/10), so nearest-point decoding is NOT the ceiling.  A3 enumerates")
print("every lattice point within R = 1.1 * sqrt(m*(K1*S_K1)^2/3 +")
print("m*(K2*S_K2)^2/3) of t -- a radius computed from public parameters only --")
print("and accepts if any candidate decodes to the true d.\n")

print("#cand is the enumeration list size = the number of extra scalar")
print("multiplications the attack pays; 2000.0 means the nr_solutions cap was")
print("hit, so that row's success rate is a LOWER bound.\n")
print(f"  {'K1':>4}  {'A1-LLL':>8}{'A1-BKZ40':>10}{'A2-CVP*':>9}{'A3-LIST':>9}"
      f"   {'mean #cand':>11}")
for k1 in (4, 6, 8, 12, 16):
    w_l, _, _ = rate(run_A1, curve2677, 10, k1, SEEDS10)
    w_b, _, _ = rate(run_A1, curve2677, 10, k1, SEEDS10, use_bkz=True, beta=40)
    w_c, _, _ = rate(run_A2, curve2677, 10, k1, SEEDS10, exact=True)
    wins, cands = 0, []
    for seed in SEEDS10:
        d_trial = random.Random(seed + 7777).randint(1, n2677 - 1)
        r = run_A3_list(curve2677, 10, d_trial, k1, seed)
        wins += bool(r['ok'])
        cands.append(r['n_cand'])
    print(f"  {k1:>4}  {w_l:>6}/10{w_b:>8}/10{w_c:>7}/10{wins:>7}/10"
          f"   {sum(cands)/len(cands):>11.1f}")
print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("SUMMARY")
print("-" * 78)
ARMS = ('A0-LLL', 'A1-LLL', 'A1-BKZ20', 'A1-BKZ40', 'A2-Babai', 'A2-EXACT')
print("Wall = first K1 in the grid where the arm drops to a MINORITY of seeds")
print("(<3/5).  NOTE: A2-EXACT is the ceiling only for NEAREST-POINT decoding;")
print("EXP 7 and EXP 8 both beat it, so it is not an information-theoretic")
print("bound.  The list-decoding arm (EXP 8) is the best ceiling measured here,")
print("and its K1=12/16 rows are LOWER BOUNDS -- enumeration was capped at")
print("nr_solutions=2000, which both rows hit.\n")
maj = (len(SEEDS) + 1) // 2
for (label, m), grid in walls4.items():
    walls = {}
    for idx, arm in enumerate(ARMS):
        walls[arm] = next((k1 for k1 in K1_GRID if grid[k1][idx] < maj), None)
    print(f"  {label} (m={m}):")
    for arm in ARMS:
        print(f"      {arm:<10} wall K1 = {walls[arm]}")
    red = walls['A1-BKZ40']
    info = walls['A2-EXACT']
    if red is not None and info is not None and red < info:
        print(f"      -> REDUCTION GAP: BKZ-40 fails at K1={red} but the")
        print(f"         information is still present up to K1={info}.")
    elif red is not None and info is not None and red >= info:
        print(f"      -> NO GAP vs nearest-point CVP (both wall at K1={info}); "
              f"see EXP 8 for list decoding.")
print("\nDone.")
