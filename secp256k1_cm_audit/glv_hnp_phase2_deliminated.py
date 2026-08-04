"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md 2026-07-29, finding T5):
  The (2m+2)-dim Phase-2 lattice of `glv_hnp_phase2_20bit.py:263` always contains the
  trivial vector  t = n*S_D*e_m  (norm n*S_D = n), which is shorter than the planted
  vector (norm ~ n*sqrt(2m/3+4/3)) for every m >= 1.  So the recovery event was never
  an SVP event; sv/pv sat in [0.34, 0.37] for successes and failures alike.

Why t exists, in general:  for any column c whose pivot row R_c has all its OTHER
entries in columns that DO carry an n*e_j reduction row, the vector n*R_c reduces to
n*S_c*e_c in L.  Its norm is n*S_c.  So every non-congruence column contributes a
trivial vector of norm n*S_c, and it is competitive with the planted vector
(~n*sqrt(dim/3)) exactly when S_c = O(sqrt(dim)).  The d-column is the ONLY column
with S_c = O(1), because d is the only full-size (~n) unknown.  Every other unknown
is small, so its scale S_c ~ n/bound is huge and its trivial vector is irrelevant.

  => The fix is not to rescale (2026-07-29 showed S_D cancels).  The fix is to
     ELIMINATE the full-size unknown d algebraically.

Construction L23 (d-eliminated, dim 2m+1).  Signatures give
    k1_i + lam*k2_i == A_i + B_i*d   (mod n),      k1_i < K1,  k2_i < K2.
Pivot on i=0:  d == B_0^{-1} (k1_0 + lam*k2_0 - A_0)  (mod n).  Substituting, for
i = 1..m-1, with C_i = B_i*B_0^{-1} mod n and D_i = A_i - C_i*A_0 mod n:

    k1_i == D_i + C_i*k1_0 + C_i*lam*k2_0 - lam*k2_i   (mod n)

Columns (2m+1):  [0] k1_0 (scale S_K1) | [1..m-1] congruence cols, value k1_i (S_K1)
                 | [m..2m-1] k2_0..k2_{m-1} (S_K2) | [2m] Kannan (S_KANNAN)
Rows (2m+1):     n*S_K1*e_i for i=1..m-1                              (reduction)
                 R_k10 :  S_K1 @ col0,   +C_i*S_K1 @ col i
                 R_k20 :  S_K2 @ col m,  +C_i*lam*S_K1 @ col i
                 R_k2j :  S_K2 @ col m+j, -lam*S_K1 @ col j           (j=1..m-1)
                 R_tgt :  D_i*S_K1 @ col i, S_KANNAN @ col 2m

Planted vector = (k1_0*S_K1, k1_1*S_K1, ..., k2_j*S_K2, ..., S_KANNAN), norm
~ n*sqrt(2m/3 + 1).  d is read back from (k1_0, k2_0) via the pivot identity.
No column has scale O(1), so no competitive trivial vector exists.

Falsifier (stated 2026-07-29): if sv/pv rises above ~1 AND the K1 wall of T4 moves
outward, the reformulation is a real improvement; if the wall stays at K1 ~ 4-6, the
wall is information-theoretic and Phase 2 is at its ceiling.

Run: python3 glv_hnp_phase2_deliminated.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
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
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def glv_roots(n):
    """The two roots of x^2+x+1 = 0 mod n (n prime, n = 1 mod 3)."""
    if n % 3 != 1: return None
    for g in range(2, min(n, 10000)):
        c = pow(g, (n - 1) // 3, n)
        if c != 1 and (c * c + c + 1) % n == 0:
            return tuple(sorted((c, (n - 1 - c) % n)))
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Shared: scales and signature generation (verbatim)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- identical to the Phase-2 convention."""
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
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# L_old : the 2026-06-15 Phase-2 lattice (dim 2m+2), verbatim
# ---------------------------------------------------------------------------

def build_lattice_old(sigs, n, lam, k1_bound, k2_bound):
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

def planted_old(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_old(reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def det_old(m, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN

# ---------------------------------------------------------------------------
# L_23 : d-eliminated lattice (dim 2m+1)
# ---------------------------------------------------------------------------

def elim_coeffs(sigs, n, lam):
    """C_i = B_i/B_0 mod n and D_i = A_i - C_i*A_0 mod n, for i = 1..m-1."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'] % n, n)
    C = [0] * m
    D = [0] * m
    for i in range(1, m):
        C[i] = sigs[i]['B'] % n * B0inv % n
        D[i] = (sigs[i]['A'] - C[i] * sigs[0]['A']) % n
    return C, D

def build_lattice_23(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    C, D = elim_coeffs(sigs, n, lam)
    M = [[0] * dim for _ in range(dim)]
    # rows 0..m-2 : reduction rows for congruence columns 1..m-1
    for i in range(1, m):
        M[i - 1][i] = n * S_K1
    # row m-1 : k1_0
    r = m - 1
    M[r][0] = S_K1
    for i in range(1, m):
        M[r][i] = C[i] * S_K1
    # row m : k2_0
    r = m
    M[r][m] = S_K2
    for i in range(1, m):
        M[r][i] = C[i] * lam * S_K1
    # rows m+1 .. 2m-1 : k2_j for j = 1..m-1
    for j in range(1, m):
        r = m + j
        M[r][m + j] = S_K2
        M[r][j] = -lam * S_K1
    # row 2m : target / Kannan
    r = 2 * m
    for i in range(1, m):
        M[r][i] = D[i] * S_K1
    M[r][2 * m] = S_KANNAN
    return M

def planted_23(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    v[0] = sigs[0]['k1'] * S_K1
    for i in range(1, m):
        v[i] = sigs[i]['k1'] * S_K1
    for j in range(m):
        v[m + j] = sigs[j]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def planted_23_combination(sigs, n, lam, k1_bound, k2_bound):
    """Rebuild the planted vector as an explicit integer combination of the rows of
    L23.  Returns (vector, coeffs) -- proves membership rather than asserting it."""
    m = len(sigs)
    M = build_lattice_23(sigs, n, lam, k1_bound, k2_bound)
    C, D = elim_coeffs(sigs, n, lam)
    coeffs = [0] * (2 * m + 1)
    coeffs[2 * m] = 1                       # target row
    coeffs[m - 1] = sigs[0]['k1']           # k1_0 row
    coeffs[m] = sigs[0]['k2']               # k2_0 row
    for j in range(1, m):
        coeffs[m + j] = sigs[j]['k2']
    for i in range(1, m):
        # congruence: k1_i == D_i + C_i*k1_0 + C_i*lam*k2_0 - lam*k2_i  (mod n)
        num = (D[i] + C[i] * sigs[0]['k1'] + C[i] * lam * sigs[0]['k2']
               - lam * sigs[i]['k2'] - sigs[i]['k1'])
        assert num % n == 0, "congruence identity failed -- construction bug"
        coeffs[i - 1] = -num // n
    v = [sum(coeffs[r] * M[r][c] for r in range(2 * m + 1)) for c in range(2 * m + 1)]
    return v, coeffs

def recover_23(reduced, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_0, k2_0) off a Kannan row and invert the pivot identity for d."""
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    dim = 2 * m + 1
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        r0, r1 = sign * row[0], sign * row[m]
        if r0 % S_K1 or r1 % S_K2: continue
        k1_0, k2_0 = r0 // S_K1, r1 // S_K2
        d_cand = B0inv * (k1_0 + lam * k2_0 - sigs[0]['A']) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def build_lattice_cvp(sigs, n, lam, k1_bound, k2_bound):
    """L23 with the Kannan row/column stripped (dim 2m), plus the CVP target.

    planted[0:2m] = t + w  with  w = closest lattice vector to -t, so recovering
    the planted vector becomes an exact CVP instead of a uSVP.  This is the second
    route proposed 2026-07-29 ("replace the Kannan embedding with an explicit CVP
    call"); route one is build_lattice_23 above.
    """
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    C, D = elim_coeffs(sigs, n, lam)
    M = [[0] * dim for _ in range(dim)]
    for i in range(1, m):
        M[i - 1][i] = n * S_K1
    r = m - 1
    M[r][0] = S_K1
    for i in range(1, m):
        M[r][i] = C[i] * S_K1
    r = m
    M[r][m] = S_K2
    for i in range(1, m):
        M[r][i] = C[i] * lam * S_K1
    for j in range(1, m):
        M[m + j][m + j] = S_K2
        M[m + j][j] = -lam * S_K1
    t = [0] * dim
    for i in range(1, m):
        t[i] = D[i] * S_K1
    return M, t

def run_cvp(curve, m, d_secret, k1_bound, seed=42):
    """Exact-CVP variant.  Returns True iff d is recovered."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M, t = build_lattice_cvp(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    target = tuple(-x for x in t)
    try:
        w = CVP.closest_vector(A, target)
    except Exception:
        return None
    v = [t[i] + w[i] for i in range(2 * m)]
    if v[0] % S_K1 or v[m] % S_K2:
        return False
    k1_0, k2_0 = v[0] // S_K1, v[m] // S_K2
    B0inv = modinv(sigs[0]['B'] % n, n)
    d_cand = B0inv * (k1_0 + lam * k2_0 - sigs[0]['A']) % n
    return d_cand == d_secret

def rate_cvp(curve, m, k1_bound, seeds):
    wins, tot = 0, 0
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = run_cvp(curve, m, d_trial, k1_bound, seed)
        if res is None:
            continue
        tot += 1
        wins += bool(res)
    return wins, tot

def det_23(m, n, k1_bound, k2_bound):
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return n ** (m - 1) * S_K1 ** m * S_K2 ** m * S_KAN

# ---------------------------------------------------------------------------
# Gaussian heuristic / uSVP gap model
# ---------------------------------------------------------------------------

def gaussian_heuristic(dim, det):
    """GH(L) = sqrt(dim/(2*pi*e)) * det^(1/dim), computed in logs to avoid overflow."""
    log_det = math.log(det)
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def pv_expected(which, m, n, k1_bound, k2_bound):
    """E[||v_planted||] with k1~U[0,K1), k2~U[0,K2), d~U[0,n): E[x^2] = X^2/3."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    t = (m * (k1_bound * S_K1) ** 2 / 3.0
         + m * (k2_bound * S_K2) ** 2 / 3.0
         + S_KAN ** 2)
    if which == 'old':
        t += (n * S_D) ** 2 / 3.0
    return math.sqrt(t)

# ---------------------------------------------------------------------------
# Experiment driver
# ---------------------------------------------------------------------------

def run_one(which, curve, m, d_secret, k1_bound, seed=42, use_bkz=False, bkz_beta=20):
    """Returns (recovered, ||planted||, ||shortest reduced||, trivial_norm)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    if which == 'old':
        M = build_lattice_old(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_old(sigs, d_secret, n, k1_bound, k2_bound)
        triv = float(n * S_D)
    else:
        M = build_lattice_23(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_23(sigs, n, k1_bound, k2_bound)
        triv = float(n * S_K1)           # smallest column-trivial vector in L23
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    if which == 'old':
        ok = recover_old(reduced, m, n, S_KANNAN, d_secret) is not None
    else:
        ok = recover_23(reduced, sigs, m, n, lam, k1_bound, k2_bound,
                        d_secret) is not None
    nz = [norm(r) for r in reduced if any(r)]
    sn = min(nz) if nz else float('nan')
    return (ok, norm(pv), sn, triv)

def rate(which, curve, m, k1_bound, seeds, use_bkz=False, bkz_beta=20):
    """(wins, total, mean sv/pv)."""
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = run_one(which, curve, m, d_trial, k1_bound, seed,
                      use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        ok, pn, sn, _ = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================
print("=" * 78)
print("Thread 23 — d-eliminated GLV-HNP Phase-2 lattice (target = lambda_1?)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m       (as used 2026-07-29)
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",      2677, 2, 2647, 185,  8,  10),
]

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))

print("\ncurve            n     lam*    K1   m    K2   eff=K1*K2/n")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    print(f"{label:<14} {n:>5}  {lam_star(lam,n):.4f}  {k1:>3}  {m:>2}  {k2:>3}   "
          f"{k1*k2/n:.4f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-1: construction check — planted vector is an explicit integer")
print("        combination of the L23 rows, and d is recoverable from it")
print("-" * 78)
print(f"{'curve':<14} {'dim_old':>7} {'dim_23':>6} {'in L23':>7} {'d recovered':>12} "
      f"{'det check':>10}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sg = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    v_comb, _ = planted_23_combination(sg, n, lam, k1, k2b)
    v_dir = planted_23(sg, n, k1, k2b)
    in_lat = (v_comb == v_dir)
    # recover d from the planted vector alone
    got = recover_23([v_dir], sg, m, n, lam, k1, k2b, d0)
    M23 = build_lattice_23(sg, n, lam, k1, k2b)
    det_exact = abs(sympy.Matrix(M23).det())
    det_ok = (det_exact == det_23(m, n, k1, k2b))
    print(f"{label:<14} {2*m+2:>7} {2*m+1:>6} {str(in_lat):>7} "
          f"{str(got == d0):>12} {str(det_ok):>10}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-2: does the trivial vector go away?  ||shortest||/||planted|| after LLL")
print("        old: trivial n*S_D = n dominates (2026-07-29 T5, sv/pv ~ 0.34-0.60)")
print("-" * 78)
print(f"{'curve':<14} {'K1':>3} {'sv/pv old':>10} {'sv/pv 23':>9} "
      f"{'triv/pv old':>12} {'triv/pv 23':>11} {'win old':>8} {'win 23':>7}")
t23_2 = []
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    wo, to, ro = rate('old', curve, m, k1, SEEDS)
    wn, tn, rn = rate('23', curve, m, k1, SEEDS)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    pv_o = pv_expected('old', m, n, k1, k2b)
    pv_n = pv_expected('23', m, n, k1, k2b)
    to_r, tn_r = (n * S_D) / pv_o, (n * S_K1) / pv_n
    t23_2.append((label, ro, rn, to_r, tn_r, wo, wn, to))
    print(f"{label:<14} {k1:>3} {ro:>10.3f} {rn:>9.3f} {to_r:>12.3f} "
          f"{tn_r:>11.1f} {str(wo)+'/'+str(to):>8} {str(wn)+'/'+str(tn):>7}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-3: the K1 wall — replication of the 2026-07-29 T4 grid, old vs L23")
print("        falsifier: does the wall move outward under the reformulation?")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
# both curves at BOTH m values, so the lam* comparison is not confounded by m
# (2026-07-29 T4 ran 2557 at m=8 and 2677 at m=10)
for label, curve, _k1, _m0 in hist[1:]:
  for m in (8, 10):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"\n  {label}  (lam* = {lam_star(lam,n):.3f}, n = {n}, m = {m}, K2 = {k2b})")
    print("    " + "".join(f"{'K1='+str(k):>9}" for k in K1_GRID))
    print("    eff " + "".join(f"{k*k2b/n:>9.3f}" for k in K1_GRID))
    for which in ('old', '23'):
        cells = []
        for k1 in K1_GRID:
            w, t, _ = rate(which, curve, m, k1, SEEDS)
            cells.append(f"{str(w)+'/'+str(t):>9}")
        print(f"    {which:<3} " + "".join(cells))
    # BKZ-40 on the reformulated lattice at the wall
    cells = []
    for k1 in K1_GRID:
        w, t, _ = rate('23', curve, m, k1, SEEDS, use_bkz=True, bkz_beta=40)
        cells.append(f"{str(w)+'/'+str(t):>9}")
    print("    b40 " + "".join(cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-4: more signatures at the wall (K1 = 8) — old vs L23")
print("-" * 78)
M_GRID = [8, 12, 16, 24, 32]
for label, curve, _k1, _m in hist[1:]:
    p, b, n, lam, G = curve
    print(f"\n  {label}   K1 = 8")
    print("    " + "".join(f"{'m='+str(mm):>9}" for mm in M_GRID))
    for which in ('old', '23'):
        cells = []
        for mm in M_GRID:
            w, t, _ = rate(which, curve, mm, 8, SEEDS)
            cells.append(f"{str(w)+'/'+str(t):>9}")
        print(f"    {which:<3} " + "".join(cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-5: uSVP model — is the K1 wall just  ||planted|| > GH(L23) ?")
print("        GH = sqrt(dim/(2*pi*e)) * det^(1/dim);  recovery expected iff pv/GH < 1")
print("-" * 78)
print(f"{'curve':<14} {'K1':>3} {'eff':>6} {'pv/GH old':>10} {'pv/GH 23':>9} "
      f"{'pred':>5} {'obs old':>8} {'obs 23':>7}")
model_rows = []
for label, curve, _k1, m in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        gh_o = gaussian_heuristic(2 * m + 2, det_old(m, n, k1, k2b))
        gh_n = gaussian_heuristic(2 * m + 1, det_23(m, n, k1, k2b))
        r_o = pv_expected('old', m, n, k1, k2b) / gh_o
        r_n = pv_expected('23', m, n, k1, k2b) / gh_n
        wo, to, _ = rate('old', curve, m, k1, SEEDS)
        wn, tn, _ = rate('23', curve, m, k1, SEEDS)
        pred = 'YES' if r_n < 1.0 else 'no'
        model_rows.append((label, k1, k1 * k2b / n, r_o, r_n, r_n < 1.0, wo, wn, tn))
        print(f"{label:<14} {k1:>3} {k1*k2b/n:>6.3f} {r_o:>10.3f} {r_n:>9.3f} "
              f"{pred:>5} {str(wo)+'/'+str(to):>8} {str(wn)+'/'+str(tn):>7}")

hits = sum(1 for r in model_rows if r[5] == (r[7] > r[8] / 2))
print(f"\n  model agreement (pv/GH < 1  <->  majority of seeds recover): "
      f"{hits}/{len(model_rows)}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-6: is the small BKZ-40 gain on 12-bit/2677 due to L23 or just to BKZ?")
print("        matched old-vs-23 under LLL and BKZ-40, 10 seeds")
print("-" * 78)
SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 20260804, 99, 123456, 31415]
c2677 = [c for c in hist if c[0] == "12-bit/2677"][0][1]
print(f"{'m':>3} {'K1':>3} {'LLL old':>8} {'LLL 23':>7} {'BKZ40 old':>10} "
      f"{'BKZ40 23':>9}")
for mm in (8, 10):
    for k1 in (6, 8):
        a = rate('old', c2677, mm, k1, SEEDS10)
        b_ = rate('23', c2677, mm, k1, SEEDS10)
        c_ = rate('old', c2677, mm, k1, SEEDS10, use_bkz=True, bkz_beta=40)
        d_ = rate('23', c2677, mm, k1, SEEDS10, use_bkz=True, bkz_beta=40)
        f = lambda r: f"{r[0]}/{r[1]}"
        print(f"{mm:>3} {k1:>3} {f(a):>8} {f(b_):>7} {f(c_):>10} {f(d_):>9}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-7: seed variance — the 'wall' cells with 60 seeds instead of 5")
print("        the 5 historical seeds gave 0/5 where 10 seeds gave 4/10, so the")
print("        wall may be a soft probabilistic boundary, not a hard 0")
print("-" * 78)
SEEDS60 = [42, 1234, 9999, 555, 31337] + [1000 * i + 3 for i in range(55)]
print(f"{'curve':<14} {'m':>3} {'K1':>3} {'old 60':>8} {'23 60':>7} {'rate':>7} "
      f"{'hist 5':>7} {'rest 55':>8}")
CELLS = [("12-bit/2677", 8, 4), ("12-bit/2677", 8, 6), ("12-bit/2677", 8, 8),
         ("12-bit/2677", 10, 6), ("12-bit/2557", 8, 8), ("12-bit/2557", 8, 12),
         ("12-bit/2557", 8, 16)]
for lbl, mm, k1 in CELLS:
    cur = [c for c in hist if c[0] == lbl][0][1]
    wo, to, _ = rate('old', cur, mm, k1, SEEDS60)
    wn, tn, _ = rate('23', cur, mm, k1, SEEDS60)
    h5, _, _ = rate('old', cur, mm, k1, SEEDS60[:5])
    r55, _, _ = rate('old', cur, mm, k1, SEEDS60[5:])
    print(f"{lbl:<14} {mm:>3} {k1:>3} {str(wo)+'/'+str(to):>8} "
          f"{str(wn)+'/'+str(tn):>7} {wo/to:>7.2f} {str(h5)+'/5':>7} "
          f"{str(r55)+'/55':>8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-8: route two — exact CVP on the d-eliminated lattice (no Kannan row)")
print("        if the wall is an embedding artifact, exact CVP should break it")
print("-" * 78)
SEEDS20 = SEEDS60[:20]
print(f"{'curve':<14} {'m':>3} {'K1':>3} {'LLL old':>8} {'LLL 23':>7} {'CVP':>6}")
for lbl, mm, k1 in CELLS:
    cur = [c for c in hist if c[0] == lbl][0][1]
    wo, to, _ = rate('old', cur, mm, k1, SEEDS20)
    wn, tn, _ = rate('23', cur, mm, k1, SEEDS20)
    wc, tc = rate_cvp(cur, mm, k1, SEEDS20)
    print(f"{lbl:<14} {mm:>3} {k1:>3} {str(wo)+'/'+str(to):>8} "
          f"{str(wn)+'/'+str(tn):>7} {str(wc)+'/'+str(tc):>6}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("T23-9: confirmation — CVP vs LLL at 60 seeds on the cells where they differ")
print("        LLL gives 0/60 at 12-bit/2677 K1=8; is the CVP hit real?")
print("-" * 78)
print(f"{'curve':<14} {'m':>3} {'K1':>3} {'LLL 23 /60':>11} {'CVP /60':>8} "
      f"{'delta':>7}")
for lbl, mm, k1 in CELLS:
    cur = [c for c in hist if c[0] == lbl][0][1]
    wn, tn, _ = rate('23', cur, mm, k1, SEEDS60)
    wc, tc = rate_cvp(cur, mm, k1, SEEDS60)
    print(f"{lbl:<14} {mm:>3} {k1:>3} {str(wn)+'/'+str(tn):>11} "
          f"{str(wc)+'/'+str(tc):>8} {wc/tc - wn/tn:>+7.3f}")
print("\n  note: exact CVP is a SINGLE-candidate decoder (it returns the closest")
print("  lattice point, right or wrong), whereas the Kannan scan tests EVERY basis")
print("  row with |last| = S_KANNAN.  That is why CVP can win on one curve and")
print("  lose on another rather than dominating.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
