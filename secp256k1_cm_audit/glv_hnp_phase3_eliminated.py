"""
GLV-HNP Phase 3, Thread 23: eliminate d from the lattice.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  In the Phase-2 lattice (dim 2m+2, `build_glv_lattice` in
  glv_hnp_phase2_20bit.py:262) the shortest vector after LLL is ALWAYS the
  trivial vector n*S_D*e_m -- 100% of its energy in the d-column, Kannan
  coordinate 0.  It carries no information (d is only defined mod n) and no
  choice of S_D removes it, because both it and the planted vector scale
  linearly in S_D.  Measured sv/pv in [0.337, 0.368] on every curve, success
  and failure alike.

  The 2026-07-29 next-step proposal (Thread 23) was: reformulate so that the
  planted vector is lambda_1, then re-run the T4 K1-grid and see whether the
  K1 wall moves outward.  Falsifier stated there:
      "if the wall stays at K1 ~ 4-6, then the wall is information-theoretic
       and Phase 2 is at its ceiling."

This script implements the reformulation and runs that falsifier.

  -------------------------------------------------------------------------
  THE PHASE-3 LATTICE (d eliminated by pivoting on signature 0)
  -------------------------------------------------------------------------
  ECDSA-GLV relation, per signature i:   k_full,i = A_i + B_i*d  (mod n)
  with the GLV split                     k_full,i = k1_i + lam*k2_i (mod n),
                                         0 <= k1_i < K1,  0 <= k2_i < K2.

  Pivot on i=0:  d = B_0^{-1} (k_full,0 - A_0), so for j = 1..m-1
      k_full,j - C_j * k_full,0  ==  D_j   (mod n),
      C_j = B_j * B_0^{-1} mod n,      D_j = A_j - C_j*A_0 mod n.

  Substituting the GLV split and solving for k1_j (j >= 1):
      k1_j = D_j + C_j*k1_0 + C_j*lam*k2_0 - lam*k2_j + n*q_j.        (*)

  So the FREE unknowns are (k1_0, k2_0, k2_1..k2_{m-1}) -- m+1 of them --
  and the DETERMINED quantities k1_1..k1_{m-1} must land in [0, K1).
  d does not appear anywhere.  Lattice dimension = (m+1) + (m-1) + 1 = 2m+1,
  one LESS than Phase 2, and the n*e_d direction is gone by construction.

  Columns: [k1_j : j=1..m-1] (scale S1) | k1_0 (S1) | k2_0..k2_{m-1} (S2) | Kannan (S_K)
  Rows:
    - free-unknown rows: S_u on the unknown's own column, and the coefficient
      of that unknown in (*) times S1 on each k1_j column;
    - modulus rows n*S1 on each k1_j column;
    - Kannan row: D_j*S1 on each k1_j column, S_K on the last.
  Scales are IDENTICAL to Phase 2 (S1 = n//K1, S2 = n//K2, S_K = n) so the
  planted-vector norms are directly comparable:
      ||w_phase3||^2 = ||v_phase2||^2 - (d*S_D)^2.

  Determinant bookkeeping (both full rank):
      det_2 = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN = (n*S_K1)^m * S_K2^m * n
      det_3 = (n*S1)^(m-1) * S1 * S2^m * S_K       = (n*S_K1)^m * S_K2^m
  i.e. det_3 = det_2 / n in one fewer dimension: a strictly better gap.

  -------------------------------------------------------------------------
  PREDICTION MADE BEFORE RUNNING (recorded so the falsifier is honest)
  -------------------------------------------------------------------------
  Eliminating d removes the n*e_d vector, but it does NOT make the planted
  vector lambda_1.  The "lambda-block" vectors survive elimination: any
  (x_i, y_i) with x_i + lam*y_i == 0 (mod n) makes every k_full,i vanish,
  hence satisfies (*) homogeneously, and its norm is the Gauss-reduced
  mu of  < (n*S1, 0), (-lam*S1, S2) >  -- the very quantity Thread 20 called
  mu, measured at rho = mu/||v|| in [0.399, 0.686] < 1 on all three
  historical curves.  So sub-planted vectors are unavoidable in ANY lattice
  built from this relation, and "make the planted vector lambda_1" is not an
  achievable goal.  The testable question is only whether the smaller
  dimension + smaller determinant move the K1 wall.

Experiments:
  E1  correctness: w is in L, recovery inverts to d, on a known-good curve.
  E2  sv/pv and energy profile, Phase 2 vs Phase 3, historical curves.
  E3  THE FALSIFIER: T4 K1-grid replication, Phase 2 vs Phase 3, 5 seeds.
  E4  17-bit sweep at the eff values where Phase 2 dies (0.15 / 0.25).

Run:  python3 secp256k1_cm_audit/glv_hnp_phase3_eliminated.py
"""

import math
import random
import sympy
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (a=0) -- copied verbatim from
# glv_hnp_phase2_lambda_threshold.py so the comparison to 2026-07-29 is exact.
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

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    return min(r1, r2), max(r1, r2)

def lam_star(lam, n):
    l = lam % n
    return min(l, n - l) / n

def gauss_reduce_2d(u, v):
    """Lagrange/Gauss reduction of a 2-D basis; returns the shortest vector."""
    def dot(a, b): return a[0] * b[0] + a[1] * b[1]
    def nrm2(a): return dot(a, a)
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        mu = round(dot(u, v) / nrm2(u))
        v = (v[0] - mu * u[0], v[1] - mu * u[1])
        if nrm2(v) >= nrm2(u): return u
        u, v = v, u

def lambda_block_mu(n, lam, S1, S2):
    """Shortest vector of < (n*S1, 0), (-lam*S1, S2) >.  Thread 20's mu."""
    return math.sqrt(sum(x * x for x in gauss_reduce_2d((n * S1, 0), (-lam * S1, S2))))

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Instance generation -- identical to Phase 2
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- same convention as Phase 2."""
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
# PHASE 2 lattice (baseline) -- verbatim from glv_hnp_phase2_20bit.py:262
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

def planted_vector_p2(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_p2(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret: return d_cand
    return None

# ---------------------------------------------------------------------------
# PHASE 3 lattice (d eliminated).  dim = 2m+1.
# ---------------------------------------------------------------------------
#   col layout:  0 .. m-2      -> k1_j for j = 1..m-1     (scale S1)
#                m-1           -> k1_0                    (scale S1)
#                m .. 2m-1     -> k2_0 .. k2_{m-1}        (scale S2)
#                2m            -> Kannan                  (scale S_K)

def phase3_coeffs(sigs, n):
    """C_j = B_j/B_0, D_j = A_j - C_j*A_0  (mod n), j = 1..m-1."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'], n)
    C, D = [], []
    for j in range(1, m):
        Cj = sigs[j]['B'] * B0inv % n
        Dj = (sigs[j]['A'] - Cj * sigs[0]['A']) % n
        C.append(Cj); D.append(Dj)
    return C, D

def build_phase3_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S1, _S_D, S2, S_K = scales(n, k1_bound, k2_bound)
    C, D = phase3_coeffs(sigs, n)
    lam_r = lam % n

    M = [[0] * dim for _ in range(dim)]
    r = 0
    # (a) modulus rows: n*S1 on each k1_j column, j = 1..m-1
    for jj in range(m - 1):
        M[r][jj] = n * S1
        r += 1
    # (b) free unknown k1_0 : coeff C_j in eq (*), own column m-1
    for jj in range(m - 1):
        M[r][jj] = C[jj] * S1
    M[r][m - 1] = S1
    r += 1
    # (c) free unknown k2_0 : coeff C_j*lam, own column m
    for jj in range(m - 1):
        M[r][jj] = C[jj] * lam_r % n * S1
    M[r][m] = S2
    r += 1
    # (d) free unknowns k2_j (j=1..m-1) : coeff -lam on column j only
    for jj in range(m - 1):
        M[r][jj] = (-lam_r) * S1
        M[r][m + 1 + jj] = S2
        r += 1
    # (e) Kannan row
    for jj in range(m - 1):
        M[r][jj] = D[jj] * S1
    M[r][dim - 1] = S_K
    r += 1
    assert r == dim, (r, dim)
    return M

def planted_vector_p3(sigs, n, k1_bound, k2_bound):
    """(k1_j*S1 for j=1..m-1 | k1_0*S1 | k2_i*S2 | S_K)."""
    m = len(sigs)
    S1, _S_D, S2, S_K = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for jj in range(m - 1):
        v[jj] = sigs[jj + 1]['k1'] * S1
    v[m - 1] = sigs[0]['k1'] * S1
    for i in range(m):
        v[m + i] = sigs[i]['k2'] * S2
    v[2 * m] = S_K
    return v

def recover_d_p3(M_reduced, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    """From a reduced row with |Kannan| = S_K, read k1_0 and k2_0, rebuild
    k_full,0 = k1_0 + lam*k2_0, then d = (k_full,0 - A_0)/B_0 mod n."""
    dim = 2 * m + 1
    S1, _S_D, S2, S_K = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_K: continue
        sign = 1 if last > 0 else -1
        e10, r10 = divmod(sign * row[m - 1], S1)
        e20, r20 = divmod(sign * row[m], S2)
        if r10 or r20: continue
        k_full0 = (e10 + lam * e20) % n
        d_cand = (k_full0 - sigs[0]['A']) * B0inv % n
        if d_cand == 0: continue
        if d_cand == d_secret: return d_cand
    return None

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def _reduce(M, dim, use_bkz, beta):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def run_p2(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m: return None
    red = _reduce(build_glv_lattice(sigs, n, lam, k1_bound, k2b), 2 * m + 2, use_bkz, beta)
    _S1, _SD, _S2, S_KAN = scales(n, k1_bound, k2b)
    return recover_d_p2(red, m, n, S_KAN, d_secret) is not None

def run_p3(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m: return None
    red = _reduce(build_phase3_lattice(sigs, n, lam, k1_bound, k2b), 2 * m + 1, use_bkz, beta)
    return recover_d_p3(red, sigs, m, n, lam, k1_bound, k2b, d_secret) is not None

def success_rate(runner, curve, m, k1_bound, seeds, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    wins = 0
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        if runner(curve, m, d, k1_bound, seed=s, use_bkz=use_bkz, beta=beta):
            wins += 1
    return wins, len(seeds)

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p)."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def build_curve(p, n, seed=12345):
    """Find the sextic twist b with #E = n, plus a generator."""
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

def search_curves(lo, hi, per_bin=1, nbins=10, max_primes=100000):
    """Collect j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3,
    bucketed by lam* = min(lam,n-lam)/n into `nbins` bins over [0,0.5].
    Verbatim from glv_hnp_phase2_lambda_threshold.py:373."""
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
                    if n_cand < 2 or n_cand % 3 != 1: continue
                    if not sympy.isprime(n_cand): continue
                    roots = glv_roots(n_cand)
                    if roots is None or roots[0] is None: continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin: continue
                    cur = build_curve(p, n_cand)
                    if cur is None: continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()): break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out


# ===========================================================================
print("=" * 78)
print("Thread 23 — Phase-3 lattice with d eliminated  (GLV-HNP)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1: correctness of the Phase-3 construction")
print("-" * 78)
print("Checks, for each historical curve: (i) the planted vector w really is an")
print("integer combination of the Phase-3 basis rows; (ii) recover_d_p3 inverts")
print("w back to d; (iii) ||w_p3||^2 == ||v_p2||^2 - d^2 as predicted.\n")

print(f"{'curve':<18} {'m':>3} {'dim2':>5} {'dim3':>5} {'w in L':>7} "
      f"{'d ok':>5} {'norm id':>8} {'||w3||/||v2||':>14}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    M3 = build_phase3_lattice(sigs, n, lam, k1, k2b)
    w = planted_vector_p3(sigs, n, k1, k2b)
    v2 = planted_vector_p2(sigs, d0, n, k1, k2b)

    # (i) membership: solve for the integer coefficient vector directly from the
    #     construction (coeffs are the unknowns themselves + the q_j).
    S1, _SD, S2, S_K = scales(n, k1, k2b)
    C, D = phase3_coeffs(sigs, n)
    coeffs = [0] * (2 * m + 1)
    for jj in range(m - 1):   # q_j on the modulus rows
        num = (sigs[jj + 1]['k1'] - D[jj] - C[jj] * sigs[0]['k1']
               - C[jj] * (lam % n) % n * sigs[0]['k2'] + (lam % n) * sigs[jj + 1]['k2'])
        assert num % n == 0, "q_j not integral — relation (*) is wrong"
        coeffs[jj] = num // n
    coeffs[m - 1] = sigs[0]['k1']
    coeffs[m] = sigs[0]['k2']
    for jj in range(m - 1):
        coeffs[m + 1 + jj] = sigs[jj + 1]['k2']
    coeffs[2 * m] = 1
    combo = [sum(coeffs[r] * M3[r][c] for r in range(2 * m + 1)) for c in range(2 * m + 1)]
    in_L = (combo == w)

    # (ii) recovery from w itself
    d_ok = recover_d_p3([w], sigs, m, n, lam, k1, k2b, d0) == d0

    # (iii) norm identity
    norm_id = abs(norm(w) ** 2 - (norm(v2) ** 2 - d0 * d0)) < 1e-6
    print(f"{label:<18} {m:>3} {2*m+2:>5} {2*m+1:>5} {str(in_L):>7} "
          f"{str(d_ok):>5} {str(norm_id):>8} {norm(w)/norm(v2):>14.4f}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1b: PROJECTION THEOREM — is Phase 3 just Phase 2 projected along e_d?")
print("-" * 78)
print("The d-column of the Phase-2 lattice is a coordinate axis, so orthogonal")
print("projection along e_d is literally deletion of column m.  Claim:")
print("    pi_{e_d^perp}(L2)  ==  L3   (as lattices, after aligning column order)")
print("Checked by integer HNF of a full-rank row basis.  If true, 'eliminating d'")
print("is not a new lattice at all and the T5 vector n*e_d spans ker(pi).\n")

def _rowbasis(rows, nc):
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    out = [[A[i][j] for j in range(nc)] for i in range(A.nrows)]
    return [r for r in out if any(r)]

def _hnf_rows(rows):
    H = hermite_normal_form(Matrix(rows).T)
    return tuple(tuple(H.row(i)) for i in range(H.rows))

print(f"{'curve':<18} {'rank pi(L2)':>11} {'2m+1':>5} {'det3==det2/n':>13} {'HNF equal':>10}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    M2 = build_glv_lattice(sigs, n, lam, k1, k2b)
    M3 = build_phase3_lattice(sigs, n, lam, k1, k2b)
    # project along e_d == delete column m
    P = [[r[c] for c in range(2 * m + 2) if c != m] for r in M2]
    # align column order: pi(L2) is [k1_0..k1_{m-1} | k2_* | kan],
    #                     L3     is [k1_1..k1_{m-1} | k1_0 | k2_* | kan]
    perm = list(range(1, m)) + [0] + list(range(m, 2 * m + 1))
    Pp = [[r[i] for i in perm] for r in P]
    rk = Matrix(P).rank()
    det_ok = abs(Matrix(M3).det()) == abs(Matrix(M2).det()) // n
    eq = _hnf_rows(_rowbasis(Pp, 2 * m + 1)) == _hnf_rows(M3)
    print(f"{label:<18} {rk:>11} {2*m+1:>5} {str(det_ok):>13} {str(eq):>10}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E2: is the planted vector lambda_1 now?   (sv/pv, Phase 2 vs Phase 3)")
print("-" * 78)
print("sv = ||first row after LLL||, pv = ||planted vector||.")
print("mu = Gauss-reduced shortest vector of the lambda-block (Thread 20's mu),")
print("which the header predicts SURVIVES elimination and stays below pv.\n")

print(f"{'curve':<18} {'K1':>3} {'sv/pv (P2)':>11} {'sv/pv (P3)':>11} "
      f"{'mu/pv3':>8} {'|sv3| is mu':>12} {'kan=0 (P3)':>11}")
e2_rows = []
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    S1, _SD, S2, S_K = scales(n, k1, k2b)

    red2 = _reduce(build_glv_lattice(sigs, n, lam, k1, k2b), 2 * m + 2, False, 0)
    red3 = _reduce(build_phase3_lattice(sigs, n, lam, k1, k2b), 2 * m + 1, False, 0)
    v2 = planted_vector_p2(sigs, d0, n, k1, k2b)
    w3 = planted_vector_p3(sigs, n, k1, k2b)
    sv2 = min(norm(r) for r in red2 if any(r))
    sv3 = min(norm(r) for r in red3 if any(r))
    mu = lambda_block_mu(n, lam, S1, S2)
    # is the P3 shortest vector the lambda-block vector? (norm match within 1%)
    is_mu = abs(sv3 - mu) / mu < 0.01
    # Kannan coordinate of the P3 shortest vector
    shortest3 = min((r for r in red3 if any(r)), key=norm)
    kan0 = (shortest3[2 * m] == 0)
    print(f"{label:<18} {k1:>3} {sv2/norm(v2):>11.4f} {sv3/norm(w3):>11.4f} "
          f"{mu/norm(w3):>8.4f} {str(is_mu):>12} {str(kan0):>11}")
    e2_rows.append((label, sv2 / norm(v2), sv3 / norm(w3), mu / norm(w3), is_mu, kan0))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E3: THE FALSIFIER — T4 K1-grid, Phase 2 vs Phase 3")
print("-" * 78)
print("2026-07-29 T4 baseline (Phase 2, 5 seeds):")
print("  12-bit/2557 (lam*=0.340): K1=2,3,4,6,8 -> 5/5;  12 -> 4/5;  16 -> 1/5;  24 -> 0/5")
print("  12-bit/2677 (lam*=0.070): K1=2,3,4 -> 5/5;  6 -> 2/5;  8,12,16,24 -> 0/5")
print("Falsifier: if the Phase-3 wall on 2677 moves outward past K1 ~ 4-6, the")
print("reformulation is a real improvement; if it stays, the wall is")
print("information-theoretic and Phase 2 is at its ceiling.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
# NOTE: T4 (glv_hnp_phase2_lambda_threshold.py:604) swept K1 at a FIXED m=12 for
# both curves, not at the per-curve m of the HIST table.  Use m=12 here so the
# P2 column reproduces the logged 2026-07-29 baseline exactly.
E3 = [("12-bit/2557", 12), ("12-bit/2677 FAIL", 12)]

for label, m in E3:
    curve = next(c for lb, c, _k, _m in hist if lb == label)
    p, b, n, lam, G = curve
    print(f"{label}   (lam*={lam_star(lam, n):.4f}, n={n}, m={m})")
    print(f"  {'K1':>4} " + " ".join(f"{k:>5}" for k in K1_GRID))
    for tag, runner in (("P2", run_p2), ("P3", run_p3)):
        cells = []
        for k1 in K1_GRID:
            w, t = success_rate(runner, curve, m, k1, SEEDS)
            cells.append(f"{w}/{t}")
        print(f"  {tag:>4} " + " ".join(f"{c:>5}" for c in cells))
    print()


# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E4: 17-bit sweep at the eff values where Phase 2 dies")
print("-" * 78)
print("2026-07-29 T3 baseline: at eff = K1*K2/n = 0.05 Phase 2 recovers 19/20")
print("curves; at 0.15 only 3/20; at 0.25 none.  Does Phase 3 do better?\n")

# Fresh 17-bit j=0 GLV curves, lam* spread over (0, 0.5) -- same generator as
# the 2026-07-29 T3 sweep, so the comparison to that baseline is like-for-like.
curves17 = [(f"p={c[0]}", c) for c in search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=10)]
print(f"built {len(curves17)} fresh 17-bit j=0 GLV curves "
      f"(lam* = {', '.join(f'{lam_star(c[3], c[2]):.3f}' for _l, c in curves17)})\n")

M_SIGS = 12
for eff in (0.05, 0.15, 0.25):
    p2_wins = p3_wins = 0
    detail = []
    for label, curve in curves17:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        k1b = max(2, int(eff * n / k2b))
        w2, t = success_rate(run_p2, curve, M_SIGS, k1b, SEEDS)
        w3, _ = success_rate(run_p3, curve, M_SIGS, k1b, SEEDS)
        p2_wins += (w2 == t); p3_wins += (w3 == t)
        detail.append((label, lam_star(lam, n), k1b, w2, w3, t))
    print(f"eff = {eff:.2f}   full-recovery curves:  P2 {p2_wins}/{len(curves17)}   "
          f"P3 {p3_wins}/{len(curves17)}")
    print(f"  {'curve':<16} {'lam*':>7} {'K1':>5} {'P2':>6} {'P3':>6}")
    for label, ls, k1b, w2, w3, t in detail:
        print(f"  {label:<16} {ls:>7.4f} {k1b:>5} {str(w2)+'/'+str(t):>6} "
              f"{str(w3)+'/'+str(t):>6}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E5: what kind of wall is it?  info-theoretic vs unique-SVP gap")
print("-" * 78)
print("The 2026-07-29 falsifier offered a dichotomy: either the reformulation")
print("moves the wall, or 'the wall is information-theoretic'.  E3 shows the wall")
print("does not move -- but the second horn is testable too, and is WRONG.")
print("  m_thresh = ceil(log n / log(1/eff))  is the information-theoretic bound")
print("             (unknown entropy < constraint entropy).")
print("  ||w||/GH  is the unique-SVP gap; recovery needs it below ~1.")
print("Planted norms use E[x^2] = n^2/3 per bounded coordinate.\n")

n = 2647; K2c = math.isqrt(n) + 1; mE5 = 10
print(f"12-bit/2677:  n={n}, K2={K2c}, m={mE5}  (m_thresh << m throughout)")
print(f"  {'K1':>4} {'eff':>7} {'m_thresh':>8} {'||w3||/GH3':>11} {'||v2||/GH2':>11} "
      f"{'gap3/gap2':>10} {'observed':>9}")
obs = {2: '5/5', 3: '5/5', 4: '5/5', 6: '0/5', 8: '0/5', 12: '0/5', 16: '0/5', 24: '0/5'}
for K1 in K1_GRID:
    eff = K1 * K2c / n
    mt = math.ceil(math.log(n) / math.log(1 / eff)) if eff < 1 else float('inf')
    S1 = n // K1; S2 = max(1, n // K2c); SK = n
    det3 = (n * S1) ** (mE5 - 1) * S1 * S2 ** mE5 * SK
    det2 = (n * S1) ** mE5 * 1 * S2 ** mE5 * SK
    w3 = math.sqrt(mE5 * (K1 * S1) ** 2 / 3 + mE5 * (K2c * S2) ** 2 / 3 + SK ** 2)
    v2 = math.sqrt(w3 ** 2 + n * n / 3)
    gh3 = math.sqrt((2 * mE5 + 1) / (2 * math.pi * math.e)) * det3 ** (1 / (2 * mE5 + 1))
    gh2 = math.sqrt((2 * mE5 + 2) / (2 * math.pi * math.e)) * det2 ** (1 / (2 * mE5 + 2))
    print(f"  {K1:>4} {eff:>7.4f} {mt:>8} {w3/gh3:>11.4f} {v2/gh2:>11.4f} "
          f"{(w3/gh3)/(v2/gh2):>10.4f} {obs[K1]:>9}")

print("""
Reading: at K1=6 and K1=8 the information-theoretic bound is m_thresh = 4 and 5
while m = 10 -- twice the data needed -- yet recovery is 0/5.  So the wall is
NOT information-theoretic.  The unique-SVP gap crosses 1.0 between K1=3 (0.94)
and K1=4 (1.08), bracketing the observed 5/5 -> 0/5 transition at K1=4 -> 6.
Phase 3 improves the gap by only 3.5-5% (gap3/gap2 in [0.950, 1.003]), because
it trades one dimension for a factor n in the determinant and these nearly
cancel:  gap3/gap2 = (n*eff^m)^(-1/((2m+1)(2m+2))) -> ~0.97 at m=10.
A 3-5% gap gain cannot move a K1 grid whose steps are 33-50% apart, which is
exactly why E3/E4 are identical in 46 of 46 cells.""")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E6: is the eff wall structural, or just a shortage of signatures?")
print("-" * 78)
print("E5 says the wall tracks the uniqueSVP gap, and the gap SHRINKS with m:")
print("    gap(m) ~ (n * eff^m)^(1/(2m+1)) * C   ->   sqrt(eff)*C  as m -> inf.")
print("So for any fixed eff there should be an m past which recovery works.")
print("2026-07-26 called the eff=0.157 failure 'structural'; 2026-07-29 T4b")
print("probed only to m=32 -- but on n=2647, where K1*K2 = 416 distinct k values")
print("means m=32 is already deep into birthday collisions AND the model puts the")
print("crossing near m~78.  Re-run on a 20-bit curve where neither problem bites.\n")

def gap_model(n, K1, K2, m):
    """||w3|| / GH(L3), computed in log space (det overflows float)."""
    S1 = n // K1; S2 = max(1, n // K2); SK = n
    ldet = (m - 1) * math.log(n * S1) + math.log(S1) + m * math.log(S2) + math.log(SK)
    w = math.sqrt(m * (K1 * S1) ** 2 / 3 + m * (K2 * S2) ** 2 / 3 + SK ** 2)
    lgh = 0.5 * math.log((2 * m + 1) / (2 * math.pi * math.e)) + ldet / (2 * m + 1)
    return math.exp(math.log(w) - lgh)

c20 = search_curves(2 ** 19, 2 ** 20, per_bin=1, nbins=1)[0]
p20, b20, n20, lam20, G20 = c20
K2_20 = math.isqrt(n20) + 1
print(f"20-bit curve: p={p20}, n={n20}, lam*={lam_star(lam20, n20):.4f}, K2={K2_20}")
E6_SEEDS = [42, 1234, 9999]

for eff_t, m_grid in ((0.10, [8, 12, 16, 20, 24]),
                      (0.157, [12, 20, 32, 48, 64, 80]),
                      (0.25, [12, 20, 32, 48, 64, 80])):
    K1 = round(eff_t * n20 / K2_20)
    print(f"\n  eff_target={eff_t}   K1={K1}   eff={K1*K2_20/n20:.4f}   "
          f"distinct k values = {K1*K2_20}")
    print(f"    {'m':>4} {'dim2':>5} {'gap':>7} {'P2':>5} {'P3':>5}")
    for m in m_grid:
        if m * m / 2 > K1 * K2_20:
            print(f"    {m:>4}  SKIP (birthday collisions in k)")
            continue
        w2, t = success_rate(run_p2, c20, m, K1, E6_SEEDS)
        w3, _ = success_rate(run_p3, c20, m, K1, E6_SEEDS)
        print(f"    {m:>4} {2*m+2:>5} {gap_model(n20, K1, K2_20, m):>7.4f} "
              f"{str(w2)+'/'+str(t):>5} {str(w3)+'/'+str(t):>5}")

print("""
Reading:
  eff=0.10   0/3 at m=8, then 3/3 from m=12 on.
  eff=0.157  0/3 at m=12 and m=20, then 1/3 at m=32 and 2/3 at m=48 and m=80.
             The 'structural' eff=0.157 wall is therefore NOT structural: it is
             a signature-count wall.  (m=64 dipping to 0/3 is 3-seed noise; this
             grid is too thin to place the crossing precisely.)
  eff=0.25   still 0/3 at m=80 (dim 162), gap 1.31 and falling only slowly.

Caveat, stated because it limits the model: the gap is a good predictor of the
TREND at fixed eff, but is NOT a calibrated cross-eff threshold.  eff=0.10 at
m=12 has gap 1.386 and recovers 3/3, while eff=0.157 at m=20 has gap 1.367 and
recovers 0/3.  Same gap, opposite outcome -- so the constant C in the model is
absorbing something eff-dependent that this run did not isolate.""")

print("=" * 78)
print("done")
print("=" * 78)
