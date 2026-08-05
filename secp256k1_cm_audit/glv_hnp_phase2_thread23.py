"""
GLV-HNP Thread 23 — reformulate the Phase-2 lattice so the target is lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, finding T5):
  The baseline Phase-2 lattice L0 (dim 2m+2, `build_glv_lattice` in
  glv_hnp_phase2_20bit.py:262) always contains the trivial vector n*S_D*e_m,
  of norm n*S_D, which is shorter than the planted vector for every m >= 1.
  So the planted vector is never lambda_1, and recovery is a BDD/coset
  condition rather than an SVP condition.

This script builds and compares three formulations on identical signature sets:

  V0  baseline           dim 2m+2, Kannan embedding, d carried in its own column
                         (verbatim from glv_hnp_phase2_lambda_threshold.py:254)
  V1  projected          dim 2m+1, d-column quotiented out.  The x-part lattice
                         becomes Lambda = Z*B + n*Z^m (index n^(m-1) in Z^m),
                         which removes the trivial vector entirely.
  V2  CVP                dim 2m, no Kannan row/column.  Solve CVP directly with
                         target w = (-A_i*S_K1, 0^m).  Run both Babai
                         nearest-plane and *exact* CVP (fpylll enumeration).

Exact CVP is the decisive diagnostic: if exact CVP cannot recover d, then no
amount of lattice reduction can, and the K1 wall measured on 2026-07-29 (T4) is
information-theoretic rather than a reduction-quality artefact.

Run: python3 glv_hnp_phase2_thread23.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + helpers (verbatim from glv_hnp_phase2_20bit.py:22-90)
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
    return min(lam % n, n - (lam % n)) / n

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Signatures + scales (verbatim from glv_hnp_phase2_lambda_threshold.py:227-252)
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

# ---------------------------------------------------------------------------
# V0 — baseline lattice (verbatim, glv_hnp_phase2_lambda_threshold.py:254)
# ---------------------------------------------------------------------------

def build_V0(sigs, n, lam, k1_bound, k2_bound):
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

def planted_V0(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_V0(rows, m, n, S_KANNAN, d_secret):
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

# ---------------------------------------------------------------------------
# V1 — projected lattice, d-column quotiented out.  dim 2m+1.
#
#   L1 = { (x*S_K1, y*S_K2, z*S_KANNAN) : x = z*A + d*B - lam*y  (mod n),
#                                          some d in Z }
#
# Basis (2m+1 rows):
#   (a) i=0..m-1 : -lam*S_K1*e_i + S_K2*f_i          [the k2 rows]
#   (b)           (A_i*S_K1)_i , 0^m , S_KANNAN      [the Kannan row]
#   (c)           (B'_i*S_K1)_i, 0^m , 0             [B' = B_0^{-1} B mod n]
#   (d) i=1..m-1 : n*S_K1*e_i , 0^m , 0
#
# (c)+(d) is a basis of Lambda = Z*B + n*Z^m  (B'_0 = 1, so no n*e_0 is needed);
# index of Lambda in Z^m is n^(m-1), hence covol(L1) = covol(L0)/n^2.
# ---------------------------------------------------------------------------

def bprime(sigs, n):
    """B' = B_0^{-1} * B mod n, reduced into [0,n).  B'_0 == 1."""
    B0inv = modinv(sigs[0]['B'], n)
    return [B0inv * s['B'] % n for s in sigs]

def build_V1(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    Bp = bprime(sigs, n)
    M = []
    for i in range(m):                                    # (a)
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    row = [0] * dim                                       # (b)
    for i in range(m):
        row[i] = sigs[i]['A'] * S_K1
    row[2 * m] = S_KANNAN
    M.append(row)
    row = [0] * dim                                       # (c)
    for i in range(m):
        row[i] = Bp[i] * S_K1
    M.append(row)
    for i in range(1, m):                                 # (d)
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    assert len(M) == dim
    return M

def planted_V1(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def assert_planted_in_V1(sigs, n, lam, d_secret):
    """Constructive membership proof: exhibit the integer coefficients.

    v1 = 1*(b) + sum_i k2_i*(a_i) + c*(c) + sum_{i>=1} t_i*(d_i)
    with c = k1_0 - A_0 + lam*k2_0  (exact integer, forces coord 0),
    t_i = (k1_i - A_i + lam*k2_i - c*B'_i)/n  must be integral.
    """
    m = len(sigs)
    Bp = bprime(sigs, n)
    c = sigs[0]['k1'] - sigs[0]['A'] + lam * sigs[0]['k2']
    assert (c - d_secret * sigs[0]['B']) % n == 0, "c != d*B_0 mod n"
    for i in range(1, m):
        num = sigs[i]['k1'] - sigs[i]['A'] + lam * sigs[i]['k2'] - c * Bp[i]
        assert num % n == 0, f"V1 membership failed at i={i}"
    return True

def recover_V1(rows, sigs, m, n, lam, S_K1, S_K2, S_KANNAN, d_secret):
    """Scan reduced rows for one with |last| == S_KANNAN; read (k1_0,k2_0) and
    solve d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n.

    Returns (d_or_None, self_consistent) where self_consistent is True when the
    same d falls out of index 1 as well (i.e. recovery does not need the
    d_secret oracle)."""
    dim = 2 * m + 1
    B0inv = modinv(sigs[0]['B'], n)
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        k1 = [sign * row[i] for i in range(m)]
        k2 = [sign * row[m + i] for i in range(m)]
        if any(x % S_K1 for x in k1) or any(x % S_K2 for x in k2):
            continue
        k1 = [x // S_K1 for x in k1]
        k2 = [x // S_K2 for x in k2]
        d_cand = (k1[0] + lam * k2[0] - sigs[0]['A']) * B0inv % n
        if d_cand == 0: continue
        consistent = True
        if m > 1:
            Biinv = modinv(sigs[1]['B'], n)
            d2 = (k1[1] + lam * k2[1] - sigs[1]['A']) * Biinv % n
            consistent = (d2 == d_cand)
        if d_cand == d_secret:
            return d_cand, consistent
    return None, False

# ---------------------------------------------------------------------------
# V2 — CVP formulation.  dim 2m, no Kannan row/column.
#
#   L2 = { (x*S_K1, y*S_K2) : x = d*B - lam*y (mod n), some d in Z }
#   target w = (-A_i*S_K1, 0^m);  closest u has u - w = (k1_i*S_K1, k2_i*S_K2).
# ---------------------------------------------------------------------------

def build_V2(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    Bp = bprime(sigs, n)
    M = []
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    row = [0] * dim
    for i in range(m):
        row[i] = Bp[i] * S_K1
    M.append(row)
    for i in range(1, m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    assert len(M) == dim
    return M

def target_V2(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m

def true_error_V2(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return [sigs[i]['k1'] * S_K1 for i in range(m)] + \
           [sigs[i]['k2'] * S_K2 for i in range(m)]

def d_from_error(err, sigs, m, n, lam, S_K1, S_K2):
    """Read d out of a CVP error vector (k1_i*S_K1, k2_i*S_K2)."""
    if err[0] % S_K1 or err[m] % S_K2:
        return None
    k1_0 = err[0] // S_K1
    k2_0 = err[m] // S_K2
    B0inv = modinv(sigs[0]['B'], n)
    return (k1_0 + lam * k2_0 - sigs[0]['A']) * B0inv % n

# ---------------------------------------------------------------------------
# Lattice geometry
# ---------------------------------------------------------------------------

def log2_covol(M):
    """log2 |det| of a square integer basis, via exact fraction-free Gauss."""
    import fractions
    A = [[fractions.Fraction(x) for x in row] for row in M]
    nn = len(A)
    det = fractions.Fraction(1)
    for c in range(nn):
        piv = None
        for r in range(c, nn):
            if A[r][c] != 0:
                piv = r; break
        if piv is None:
            return None
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, nn):
            if A[r][c] != 0:
                f = A[r][c] * inv
                for cc in range(c, nn):
                    A[r][cc] -= f * A[c][cc]
    return math.log2(abs(det))

def gaussian_heuristic(log2det, dim):
    return math.sqrt(dim / (2 * math.pi * math.e)) * 2 ** (log2det / dim)

# ---------------------------------------------------------------------------
# One trial across all three variants
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, beta=0, exact_cvp=True,
          want_geometry=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    out = {'m': m, 'K1': k1_bound, 'K2': k2_bound, 'seed': seed}

    def reduce_rows(M):
        A = IntegerMatrix.from_matrix(M)
        if beta:
            BKZ.reduction(A, BKZ.Param(beta))
        else:
            LLL.reduction(A)
        return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)], A

    # ---- V0 ----------------------------------------------------------------
    M0 = build_V0(sigs, n, lam, k1_bound, k2_bound)
    v0 = planted_V0(sigs, d_secret, n, k1_bound, k2_bound)
    rows0, _ = reduce_rows(M0)
    out['V0'] = recover_V0(rows0, m, n, S_KAN, d_secret) is not None
    nz0 = [r for r in rows0 if any(r)]
    out['V0_sv'] = min(norm(r) for r in nz0)
    out['V0_pv'] = norm(v0)

    # ---- V1 ----------------------------------------------------------------
    assert_planted_in_V1(sigs, n, lam, d_secret)
    M1 = build_V1(sigs, n, lam, k1_bound, k2_bound)
    v1 = planted_V1(sigs, n, k1_bound, k2_bound)
    rows1, _ = reduce_rows(M1)
    d1, cons1 = recover_V1(rows1, sigs, m, n, lam, S_K1, S_K2, S_KAN, d_secret)
    out['V1'] = d1 is not None
    out['V1_consistent'] = cons1
    nz1 = [r for r in rows1 if any(r)]
    out['V1_sv'] = min(norm(r) for r in nz1)
    out['V1_pv'] = norm(v1)

    # ---- V2 ----------------------------------------------------------------
    M2 = build_V2(sigs, n, lam, k1_bound, k2_bound)
    w = target_V2(sigs, n, k1_bound, k2_bound)
    e_true = true_error_V2(sigs, n, k1_bound, k2_bound)
    out['V2_etrue'] = norm(e_true)
    A2 = IntegerMatrix.from_matrix(M2)
    LLL.reduction(A2)
    if beta:
        BKZ.reduction(A2, BKZ.Param(beta))
    # Babai nearest plane
    try:
        u_bab = CVP.babai(A2, tuple(w))
        e_bab = [u_bab[i] - w[i] for i in range(2 * m)]
        out['V2_babai_err'] = norm(e_bab)
        d_b = d_from_error(e_bab, sigs, m, n, lam, S_K1, S_K2)
        out['V2_babai'] = (d_b == d_secret)
    except Exception as exc:                       # pragma: no cover
        out['V2_babai'] = None
        out['V2_babai_err'] = float('nan')
        out['V2_babai_exc'] = repr(exc)
    # Exact CVP (enumeration)
    if exact_cvp:
        try:
            u_ex = CVP.closest_vector(A2, tuple(w))
            e_ex = [u_ex[i] - w[i] for i in range(2 * m)]
            out['V2_exact_err'] = norm(e_ex)
            d_e = d_from_error(e_ex, sigs, m, n, lam, S_K1, S_K2)
            out['V2_exact'] = (d_e == d_secret)
            # Is the true solution actually the closest point?
            out['V2_true_is_closest'] = out['V2_exact_err'] >= out['V2_etrue'] - 1e-9
        except Exception as exc:
            out['V2_exact'] = None
            out['V2_exact_err'] = float('nan')
            out['V2_true_is_closest'] = None
            out['V2_exact_exc'] = repr(exc)

    if want_geometry:
        for tag, M, v in (('V0', M0, v0), ('V1', M1, v1)):
            ld = log2_covol(M)
            dim = len(M)
            out[tag + '_dim'] = dim
            out[tag + '_log2det'] = ld
            out[tag + '_gh'] = gaussian_heuristic(ld, dim)
        ld2 = log2_covol(M2)
        out['V2_dim'] = 2 * m
        out['V2_log2det'] = ld2
        out['V2_gh'] = gaussian_heuristic(ld2, 2 * m)
    return out


def rate(curve, m, k1_bound, seeds, d_secret=None, beta=0, exact_cvp=True):
    p, b, n, lam, G = curve
    keys = ['V0', 'V1', 'V2_babai', 'V2_exact']
    acc = {k: 0 for k in keys}
    tot = 0
    extra = {'true_is_closest': 0, 'V1_consistent': 0,
             'exact_err': [], 'etrue': []}
    for s in seeds:
        ds = d_secret if d_secret is not None else (random.Random(s).randrange(1, n))
        r = trial(curve, m, ds, k1_bound, s, beta=beta, exact_cvp=exact_cvp)
        if r is None:
            continue
        tot += 1
        for k in keys:
            if r.get(k):
                acc[k] += 1
        if r.get('V2_true_is_closest'):
            extra['true_is_closest'] += 1
        if r.get('V1_consistent'):
            extra['V1_consistent'] += 1
        if not math.isnan(r.get('V2_exact_err', float('nan'))):
            extra['exact_err'].append(r['V2_exact_err'])
            extra['etrue'].append(r['V2_etrue'])
    return acc, tot, extra


# ===========================================================================
if __name__ == '__main__':
    print("=" * 78)
    print("Thread 23 — reformulating the GLV-HNP Phase-2 lattice")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,           p,    b, n,    lam,  K1, m
        ("8-bit/199",      211,  2, 199,  106,  2,  6),
        ("12-bit/2557",    2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677",    2677, 2, 2647, 185,  8,  10),
    ]
    curves = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        curves.append((label, (p, b, n, lam, G), k1, m))

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U1: lattice geometry — does removing the d column make the")
    print("        planted vector lambda_1?")
    print("-" * 78)
    print(f"{'curve':<14} {'var':<4} {'dim':>4} {'log2det':>9} {'GH':>11} "
          f"{'||v||':>11} {'GH/||v||':>9} {'LLL sv':>11} {'sv/pv':>7}")
    for label, curve, k1, m in curves:
        p, b, n, lam, G = curve
        d_secret = random.Random(7).randrange(1, n)
        r = trial(curve, m, d_secret, k1, 42, want_geometry=True)
        for tag in ('V0', 'V1'):
            print(f"{label:<14} {tag:<4} {r[tag+'_dim']:>4} "
                  f"{r[tag+'_log2det']:>9.2f} {r[tag+'_gh']:>11.3e} "
                  f"{r[tag+'_pv']:>11.3e} "
                  f"{r[tag+'_gh']/r[tag+'_pv']:>9.3f} "
                  f"{r[tag+'_sv']:>11.3e} {r[tag+'_sv']/r[tag+'_pv']:>7.3f}")
        print(f"{'':<14} {'V2':<4} {r['V2_dim']:>4} {r['V2_log2det']:>9.2f} "
              f"{r['V2_gh']:>11.3e} {r['V2_etrue']:>11.3e} "
              f"{r['V2_gh']/r['V2_etrue']:>9.3f} {'-':>11} {'-':>7}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U2: the K1 wall (T4 replication) across all four solvers")
    print("        m as in T4; 5 seeds; success = d recovered exactly")
    print("-" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    for label, curve, _k1, m in curves[1:]:
        p, b, n, lam, G = curve
        print(f"\n{label}: n={n}, lam={lam}, lam*={lam_star(lam,n):.4f}, m={m}")
        hdr = f"  {'solver':<12}" + "".join(f"{('K1=%d'%k):>8}" for k in K1_GRID)
        print(hdr)
        rows = {k: [] for k in ['V0', 'V1', 'V2_babai', 'V2_exact']}
        closest = []
        for k1b in K1_GRID:
            acc, tot, extra = rate(curve, m, k1b, SEEDS)
            for k in rows:
                rows[k].append(f"{acc[k]}/{tot}")
            closest.append(f"{extra['true_is_closest']}/{tot}")
        for k in ['V0', 'V1', 'V2_babai', 'V2_exact']:
            print(f"  {k:<12}" + "".join(f"{v:>8}" for v in rows[k]))
        print(f"  {'true=closest':<12}" + "".join(f"{v:>8}" for v in closest))

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U3: is the wall information-theoretic?")
    print("        compare ||e_true|| to the exact CVP distance dist(w, L2)")
    print("-" * 78)
    print(f"  {'curve':<14}{'K1':>4}{'||e_true||':>13}{'dist(w,L2)':>13}"
          f"{'ratio':>8}{'true=closest':>14}{'d recovered':>13}")
    for label, curve, _k1, m in curves[1:]:
        p, b, n, lam, G = curve
        for k1b in [2, 4, 6, 8, 12, 16]:
            acc, tot, extra = rate(curve, m, k1b, SEEDS)
            if not extra['exact_err']:
                continue
            me = sum(extra['etrue']) / len(extra['etrue'])
            mc = sum(extra['exact_err']) / len(extra['exact_err'])
            print(f"  {label:<14}{k1b:>4}{me:>13.4e}{mc:>13.4e}"
                  f"{mc/me:>8.3f}{('%d/%d'%(extra['true_is_closest'],tot)):>14}"
                  f"{('%d/%d'%(acc['V2_exact'],tot)):>13}")

    print("\n" + "=" * 78)
    print("done")
