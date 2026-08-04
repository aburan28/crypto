"""
GLV-HNP Phase 2, Thread 23: eliminate the secret-key column so the planted
vector can be lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, result T5):
  The Phase-2 lattice of `glv_hnp_phase2_20bit.py::build_glv_lattice` has
  dimension 2m+2 with columns

      [ k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan ]
        scale S_K1=n/K1  S_D=1   S_K2=n/K2     S_KAN=n

  and always contains the vector  n*S_D*e_m.  Proof: n*(row m) minus
  sum_i B_i*(row i) kills every k1 column and leaves n*S_D in the d column.
  Its norm is n*S_D, while

      ||v_planted|| ~ n * sqrt(2m/3 + 4/3)  >  n*S_D

  for every m >= 1, and BOTH scale linearly in S_D, so no re-weighting can
  fix it.  The planted vector is therefore never lambda_1: recovery is a
  BDD/coset problem, not SVP.  Measured sv/pv on 2026-07-29 was in
  [0.34, 0.61] for successes and failures alike.

Root cause: d is a FULL-SIZE unknown (uniform mod n) that was nevertheless
given its own lattice coordinate.  A coordinate whose true value is ~n
contributes ~n to the target norm while the modular row n*e_m contributes
exactly n.  The standard HNP remedy is to not carry d at all.

This script implements that remedy.

  ELIMINATED LATTICE (dimension 2m+1, no d column)

  Signature relations:  k_i = A_i + B_i*d (mod n), k_i = k1_i + lam*k2_i.
  Take i=0 as pivot:    d = (k_0 - A_0)*B_0^{-1}, so for i = 1..m-1

      k_i  =  C_i + t_i*k_0   (mod n),
      t_i  =  B_i * B_0^{-1} mod n,      C_i = A_i - t_i*A_0 mod n,

  i.e.  k1_i  =  C_i + t_i*k1_0 + t_i*lam*k2_0 - lam*k2_i   (mod n).

  Free variables: k1_0, k2_0, k2_1..k2_{m-1}   (m+1 of them)
  Determined:     k1_1..k1_{m-1}               (m-1, each mod n)

  Columns (2m+1):  k1_0 | k2_0 | k1_1..k1_{m-1} | k2_1..k2_{m-1} | kannan
  Rows    (2m+1):  free-k1_0, free-k2_0, free-k2_j (j>=1),
                   n*e_{k1_i} (i>=1), kannan row carrying C_i.

  The m-1 lost equations are exactly compensated by the m-1 lost unknowns
  plus d's own log n bits, so the information balance
      m*log n >= m*(log K1 + log n/2) + log n
  is unchanged:  (m-1)*log n >= m*(log K1 + log n/2).

  Determinant:  n^{m-1} * S_K1^m * S_K2^m * S_KAN,  the same lattice volume
  the old construction had (the old one had an extra S_D=1 column and an
  extra modular row, contributing a factor 1).  What changes is that the
  short parasitic vector n*S_D*e_m is gone and the target loses its ~n-sized
  d coordinate.

Predicted target/Gaussian-heuristic ratio for the eliminated lattice:
      ||v_planted|| / GH  ~  sqrt(2*pi*e/3) * sqrt(eff)  ~  2.39*sqrt(eff)
so the planted vector should be lambda_1 whenever eff < ~0.175.

FALSIFIER (stated in the 2026-07-29 log entry):
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is
   a real improvement; if the wall stays at K1 ~ 4-6, the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  E0  construction sanity: planted vector has exact integer coordinates in
      the eliminated basis; recovery round-trips d.
  E1  sv/pv, old vs eliminated, on the three historical curves.
  E2  T4 K1-grid replication, old vs eliminated (THE FALSIFIER).
  E3  T4b m-sweep at the wall on 12-bit/2677, old vs eliminated.
  E4  17-bit sweep at eff = 0.05 / 0.15 / 0.25, old vs eliminated.

Run: python3 glv_hnp_phase2_eliminated.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic -- verbatim from glv_hnp_phase2_lambda_threshold.py
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

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
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
                    if n_cand < 2 or n_cand % 3 != 1:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    ls = lam_star(lam, n_cand)
                    idx = min(nbins - 1, int(ls / (0.5 / nbins)))
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

# ---------------------------------------------------------------------------
# Signatures -- verbatim
# ---------------------------------------------------------------------------

# Column re-weighting knob for E6.  SCALE_C[0] == 1.0 reproduces the historical
# scaling byte-for-byte; c != 1 stretches the k1 columns by c and shrinks the k2
# columns by c, which leaves the kernel-block determinant n*S_K1*S_K2 unchanged
# but changes its SHAPE (and hence mu).
SCALE_C = [1.0]

def scales(n, k1_bound, k2_bound):
    s1 = n // k1_bound
    s2 = max(1, n // k2_bound)
    c = SCALE_C[0]
    if c != 1.0:
        s1 = max(1, int(round(s1 * c)))
        s2 = max(1, int(round(s2 / c)))
    return (s1, 1, s2, n)

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
# OLD lattice (2m+2, with the d column) -- verbatim from the 2026-07-29 script
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

def planted_vector_old(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_old(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# NEW: eliminated lattice (2m+1, no d column)
# ---------------------------------------------------------------------------
#
# Column layout for m signatures:
#     col 0            -> k1_0
#     col 1            -> k2_0
#     col 1+i, i=1..m-1-> k1_i        (so cols 2 .. m)
#     col m+j, j=1..m-1-> k2_j        (so cols m+1 .. 2m-1)
#     col 2m           -> kannan
#
# Row layout:
#     row 0            -> free k1_0
#     row 1            -> free k2_0
#     row 1+j, j=1..m-1-> free k2_j
#     row m+i, i=1..m-1-> n * e_{k1_i}
#     row 2m           -> kannan (carries C_i)

def _k1col(i, m):
    return 0 if i == 0 else 1 + i

def _k2col(j, m):
    return 1 if j == 0 else m + j

def elim_coeffs(sigs, n):
    """(t_i, C_i) for i = 1..m-1 from the pivot signature 0."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'], n)
    t = [0] * m
    C = [0] * m
    for i in range(1, m):
        t[i] = sigs[i]['B'] * B0inv % n
        C[i] = (sigs[i]['A'] - t[i] * sigs[0]['A']) % n
    return t, C

def build_elim_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    lam_r = lam % n
    t, C = elim_coeffs(sigs, n)

    M = [[0] * dim for _ in range(dim)]

    # free k1_0
    M[0][_k1col(0, m)] = S_K1
    for i in range(1, m):
        M[0][_k1col(i, m)] = t[i] * S_K1
    # free k2_0
    M[1][_k2col(0, m)] = S_K2
    for i in range(1, m):
        M[1][_k1col(i, m)] = t[i] * lam_r % n * S_K1
    # free k2_j, j >= 1
    for j in range(1, m):
        M[1 + j][_k2col(j, m)] = S_K2
        M[1 + j][_k1col(j, m)] = -lam_r * S_K1
    # modular rows for the determined k1_i
    for i in range(1, m):
        M[m + i][_k1col(i, m)] = n * S_K1
    # kannan row
    for i in range(1, m):
        M[2 * m][_k1col(i, m)] = C[i] * S_K1
    M[2 * m][dim - 1] = S_KAN
    return M

def planted_vector_elim(sigs, n, k1_bound, k2_bound):
    """(k1_i*S_K1, k2_i*S_K2, S_KAN) laid out in the eliminated columns.
    Note: no d coordinate."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[_k1col(i, m)] = sigs[i]['k1'] * S_K1
        v[_k2col(i, m)] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def planted_coeffs_elim(sigs, n, lam):
    """Exact integer coordinates of the planted vector in the eliminated
    basis.  Returns None if any q_i fails to be an integer (would mean the
    construction is wrong)."""
    m = len(sigs)
    lam_r = lam % n
    t, C = elim_coeffs(sigs, n)
    coeffs = [0] * (2 * m + 1)
    coeffs[0] = sigs[0]['k1']
    coeffs[1] = sigs[0]['k2']
    for j in range(1, m):
        coeffs[1 + j] = sigs[j]['k2']
    for i in range(1, m):
        pred = (C[i] + t[i] * sigs[0]['k1']
                + (t[i] * lam_r % n) * sigs[0]['k2']
                - lam_r * sigs[i]['k2'])
        diff = sigs[i]['k1'] - pred
        if diff % n != 0:
            return None
        coeffs[m + i] = diff // n
    coeffs[2 * m] = 1
    return coeffs

def recover_d_elim(M_reduced, sigs, n, k1_bound, k2_bound, d_secret):
    """From a row with |kannan| == S_KAN read k1_0, k2_0, rebuild k_0, then
    d = (k_0 - A_0)/B_0."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    B0inv = modinv(sigs[0]['B'], n)
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        a = sign * row[_k1col(0, m)]
        b = sign * row[_k2col(0, m)]
        if a % S_K1 or b % S_K2: continue
        k1_0, k2_0 = a // S_K1, b // S_K2
        k0 = (k1_0 + LAM_HOLDER[0] * k2_0) % n
        d_cand = (k0 - sigs[0]['A']) * B0inv % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# lam is needed inside recover_d_elim; pass it through a 1-slot holder to keep
# the signature identical in shape to recover_d_old.
LAM_HOLDER = [0]

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# The GLV kernel block -- the parasite that survives the elimination
# ---------------------------------------------------------------------------
#
# For every index i the eliminated lattice contains the 2-dimensional block
#
#     Lambda'_i  =  < (n*S_K1, 0), (-lam*S_K1, S_K2) >   in coords (k1_i, k2_i)
#
# i.e. the SCALED kernel of the GLV splitting map (k1,k2) |-> k1 + lam*k2 mod n.
# (For i >= 1 it is free-k2_i row + modular row; for i = 0 it is
#  a*row0 + b*row1 with a + lam*b = 0 mod n, whose k1_i entries t_i*(a+lam*b)
#  are then cleared by the modular rows.)
#
# These vectors are NOT artifacts of carrying d: they express the genuine
# non-uniqueness of the (k1,k2) splitting of a nonce, so no re-parametrisation
# in these coordinates can remove them.  They live in m mutually orthogonal
# coordinate planes and every copy has the same shortest vector mu.
#
# Consequence:   lambda_1(L)  <=  mu  <=  sqrt(2/sqrt(3) * n * S_K1 * S_K2)
#                              =  (2/sqrt 3)^{1/2} * n / sqrt(eff)      (Minkowski)
# while          E||v_planted|| = n * sqrt(2m/3 + 1).
# So the planted vector is provably not lambda_1 as soon as
#                2m/3 + 1  >  (2/sqrt 3) / eff.

def gauss_reduce_2d(u, v):
    """Exact shortest nonzero vector of a 2D integer lattice (Lagrange)."""
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

def kernel_block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w

def minkowski_mu_bound(n, k1_bound, k2_bound):
    """(2/sqrt3)^{1/2} * sqrt(det) for the 2D block; det = n*S_K1*S_K2."""
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return math.sqrt(2.0 / math.sqrt(3.0) * n * S_K1 * S_K2)

def planted_norm_expected_elim(m, n, k1_bound, k2_bound):
    """E||v_planted|| for the eliminated lattice: no d coordinate."""
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                     + m * (k2_bound * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)

C_GRID = [2.0 ** (j / 4.0) for j in range(-12, 13)]

def mu_over_pv(n, lam, k1_bound, k2_bound, m, c):
    """mu(c)/||v_planted||(c) for the eliminated lattice.  Pure arithmetic:
    a 2x2 Gauss reduction, no lattice reduction, O(log n)."""
    old = SCALE_C[0]
    SCALE_C[0] = c
    try:
        S1, _sd, S2, _sk = scales(n, k1_bound, k2_bound)
        mu, _w = kernel_block_mu(n, lam, S1, S2)
        pvn = math.sqrt(m * (k1_bound * S1) ** 2 / 3.0
                        + m * (k2_bound * S2) ** 2 / 3.0 + float(n) ** 2)
    finally:
        SCALE_C[0] = old
    return mu / pvn

def best_reweight(n, lam, k1_bound, k2_bound, m, grid=None):
    """argmax_c mu(c)/pv(c).  Returns (c, ratio)."""
    best = (None, -1.0)
    for c in (grid or C_GRID):
        r = mu_over_pv(n, lam, k1_bound, k2_bound, m, c)
        if r > best[1]:
            best = (c, r)
    return best

def shortest_nonzero(rows):
    best, bn = None, None
    for r in rows:
        nn = sum(x * x for x in r)
        if nn == 0: continue
        if bn is None or nn < bn:
            bn, best = nn, r
    return best, math.sqrt(bn) if bn is not None else None

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_one(curve, m, d_secret, k1_bound, seed, mode, use_bkz=False, beta=20):
    """mode in {'old','elim'}.  Returns (ok, sv_over_pv)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if mode == 'old':
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector_old(sigs, d_secret, n, k1_bound, k2_bound)
        dim = 2 * m + 2
    else:
        M = build_elim_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector_elim(sigs, n, k1_bound, k2_bound)
        dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _sv, svn = shortest_nonzero(red)
    svpv = svn / norm(pv) if svn else float('nan')
    if mode == 'old':
        ok = recover_d_old(red, m, n, S_KAN, d_secret) is not None
    else:
        LAM_HOLDER[0] = lam % n
        ok = recover_d_elim(red, sigs, n, k1_bound, k2_bound, d_secret) is not None
    return ok, svpv

def success_rate(curve, m, k1_bound, seeds, mode, use_bkz=False, beta=20):
    n = curve[2]
    wins, svpvs = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        ok, svpv = run_one(curve, m, d_trial, k1_bound, seed, mode,
                           use_bkz=use_bkz, beta=beta)
        if ok is None: continue
        wins += ok
        svpvs.append(svpv)
    return wins, len(seeds), (sum(svpvs) / len(svpvs) if svpvs else float('nan'))


def sv_structure(curve, m, k1_bound, seed):
    """Diagnose the LLL shortest vector of the ELIMINATED lattice.
    Returns (svpv, kannan_zero, in_kernel, frac_k1, frac_k2, frac_kan)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_elim_lattice(sigs, n, lam, k1_bound, k2_bound)
    pv = planted_vector_elim(sigs, n, k1_bound, k2_bound)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    w, wn = shortest_nonzero(red)
    e_tot = float(sum(x * x for x in w))
    e_k1 = sum(w[_k1col(i, m)] ** 2 for i in range(m))
    e_k2 = sum(w[_k2col(i, m)] ** 2 for i in range(m))
    e_kan = w[dim - 1] ** 2
    # kernel membership: every (a_i,b_i) satisfies a_i + lam*b_i = 0 mod n
    in_kernel = (w[dim - 1] == 0)
    if in_kernel:
        lam_r = lam % n
        for i in range(m):
            a, b_ = w[_k1col(i, m)], w[_k2col(i, m)]
            if a % S_K1 or b_ % S_K2:
                in_kernel = False; break
            if (a // S_K1 + lam_r * (b_ // S_K2)) % n != 0:
                in_kernel = False; break
    return (wn / norm(pv), w[dim - 1] == 0, in_kernel,
            e_k1 / e_tot, e_k2 / e_tot, e_kan / e_tot)


print("=" * 78)
print("Thread 23 — eliminate the d column so the planted vector can be lambda_1")
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


# ===========================================================================
# E0 — construction sanity
# ===========================================================================
print("\n" + "=" * 78)
print("E0 — eliminated-lattice construction sanity")
print("=" * 78)
print(f"{'curve':>18} {'m':>3} {'dim':>4} {'coeffs int':>11} {'v in L':>7} "
      f"{'d round-trip':>13}")

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    d_secret = random.Random(4242).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2, seed=42)
    M = build_elim_lattice(sigs, n, lam, k1, k2)
    pv = planted_vector_elim(sigs, n, k1, k2)
    co = planted_coeffs_elim(sigs, n, lam)
    dim = 2 * m + 1
    coeff_ok = co is not None
    # explicit check that sum_r co[r]*M[r] == pv
    in_L = False
    if coeff_ok:
        acc = [0] * dim
        for r in range(dim):
            if co[r] == 0: continue
            for c in range(dim):
                if M[r][c]:
                    acc[c] += co[r] * M[r][c]
        in_L = (acc == pv)
    # d round-trip straight from the planted vector (no reduction)
    LAM_HOLDER[0] = lam % n
    rt = recover_d_elim([pv], sigs, n, k1, k2, d_secret) == d_secret
    print(f"{label:>18} {m:>3} {dim:>4} {str(coeff_ok):>11} {str(in_L):>7} "
          f"{str(rt):>13}")


# ===========================================================================
# E1 — sv/pv, old vs eliminated
# ===========================================================================
print("\n" + "=" * 78)
print("E1 — shortest-vector / planted-vector ratio after LLL")
print("=" * 78)
print("sv/pv < 1  =>  the planted vector is NOT lambda_1 (BDD regime)")
print("sv/pv = 1  =>  LLL's shortest vector IS the planted vector\n")
print(f"{'curve':>18} {'K1':>3} {'m':>3} {'eff':>6} "
      f"{'old dim':>8} {'old sv/pv':>10} {'old':>6} "
      f"{'new dim':>8} {'new sv/pv':>10} {'new':>6}")

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    eff = k1 * k2 / n
    wo, to, so = success_rate(curve, m, k1, SEEDS, 'old')
    wn, tn, sn = success_rate(curve, m, k1, SEEDS, 'elim')
    print(f"{label:>18} {k1:>3} {m:>3} {eff:>6.3f} "
          f"{2*m+2:>8} {so:>10.3f} {str(wo)+'/'+str(to):>6} "
          f"{2*m+1:>8} {sn:>10.3f} {str(wn)+'/'+str(tn):>6}")

# --- E1b: what IS the surviving shortest vector? -----------------------------
print("\n" + "-" * 78)
print("E1b — structure of the surviving LLL shortest vector (eliminated lattice)")
print("-" * 78)
print(f"{'curve':>18} {'seed':>6} {'sv/pv':>7} {'kan=0':>6} {'kernel':>7} "
      f"{'E_k1':>6} {'E_k2':>6} {'E_kan':>6} | {'mu/pv':>7} {'mink/pv':>8}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1, k2)
    mu, _w = kernel_block_mu(n, lam, S_K1, S_K2)
    pvn = planted_norm_expected_elim(m, n, k1, k2)
    mink = minkowski_mu_bound(n, k1, k2)
    for seed in SEEDS[:3]:
        r = sv_structure(curve, m, k1, seed)
        if r is None: continue
        svpv, kz, ker, f1, f2, fk = r
        print(f"{label:>18} {seed:>6} {svpv:>7.3f} {str(kz):>6} {str(ker):>7} "
              f"{f1:>6.3f} {f2:>6.3f} {fk:>6.3f} | {mu/pvn:>7.3f} "
              f"{mink/pvn:>8.3f}")


# ===========================================================================
# E2 — T4 K1-grid replication (THE FALSIFIER)
# ===========================================================================
print("\n" + "=" * 78)
print("E2 — T4 K1-grid replication, old vs eliminated  (THE FALSIFIER)")
print("=" * 78)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

for label, curve, _k1, m in hist[1:]:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    print(f"\n{label}:  n={n}, lam={lam}, lam*={lam_star(lam,n):.3f}, "
          f"K2={k2}, m={m}")
    print(f"{'K1':>5} {'eff':>7} {'old':>7} {'new':>7} {'new sv/pv':>10}")
    for k1 in K1_GRID:
        eff = k1 * k2 / n
        if eff >= 1.0: break
        wo, to, _ = success_rate(curve, m, k1, SEEDS, 'old')
        wn, tn, sn = success_rate(curve, m, k1, SEEDS, 'elim')
        print(f"{k1:>5} {eff:>7.3f} {str(wo)+'/'+str(to):>7} "
              f"{str(wn)+'/'+str(tn):>7} {sn:>10.3f}")


# ===========================================================================
# E3 — T4b m-sweep at the wall
# ===========================================================================
print("\n" + "=" * 78)
print("E3 — T4b m-sweep at the wall (12-bit/2677, K1=8), old vs eliminated")
print("=" * 78)

label, curve, _k1, _m = hist[2]
p, b, n, lam, G = curve
k2 = math.isqrt(n) + 1
S_K1, _S_D, S_K2, _S_KAN = scales(n, 8, k2)
mu, _w = kernel_block_mu(n, lam, S_K1, S_K2)
m_thr = math.log(n) / math.log(n / (8 * k2))
print(f"{label}: n={n}, K1=8, K2={k2}, eff={8*k2/n:.3f}, "
      f"m_thresh={m_thr:.2f}, mu/n={mu/n:.3f}")
print("pv grows as sqrt(m) while mu is m-independent, so mu/pv is largest at")
print("the smallest information-theoretically admissible m.  Sweep DOWN too.\n")
print(f"{'m':>4} {'mu/pv':>7} {'old':>7} {'new':>7} {'new BKZ20':>10} {'new sv/pv':>10}")
for m in (5, 6, 7, 8, 10, 12, 16, 24, 32):
    pvn = planted_norm_expected_elim(m, n, 8, k2)
    wo, to, _ = success_rate(curve, m, 8, SEEDS, 'old')
    wn, tn, sn = success_rate(curve, m, 8, SEEDS, 'elim')
    if wn == tn:
        bk = "-"
    else:
        wb, tb, _ = success_rate(curve, m, 8, SEEDS, 'elim',
                                 use_bkz=True, beta=20)
        bk = f"{wb}/{tb}"
    print(f"{m:>4} {mu/pvn:>7.3f} {str(wo)+'/'+str(to):>7} {str(wn)+'/'+str(tn):>7} "
          f"{bk:>10} {sn:>10.3f}")


# ===========================================================================
# E4 — 17-bit sweep, old vs eliminated
# ===========================================================================
print("\n" + "=" * 78)
print("E4 — 17-bit sweep at eff = 0.05 / 0.15 / 0.25, old vs eliminated")
print("=" * 78)

LO, HI = 2**16, 2**17
M_SIGS = 12
print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
curves = search_curves(LO, HI, per_bin=2, nbins=10)
print(f"Found {len(curves)} curves.")

E4_ROWS = []
for eff_t in (0.05, 0.15, 0.25):
    print(f"\n--- eff target = {eff_t}  (m = {M_SIGS}, {len(SEEDS)} seeds) ---")
    print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} "
          f"{'old':>6} {'new':>6} {'new sv/pv':>10} {'mu/pv':>7}")
    n_old_full = n_new_full = 0
    rows_eff = []
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        eff = k1 * k2 / n
        cur = (p, b, n, lam, G)
        wo, to, _ = success_rate(cur, M_SIGS, k1, SEEDS, 'old')
        wn, tn, sn = success_rate(cur, M_SIGS, k1, SEEDS, 'elim')
        mpv = mu_over_pv(n, lam, k1, k2, M_SIGS, 1.0)
        n_old_full += (wo == to)
        n_new_full += (wn == tn)
        rows_eff.append({'mu_pv': mpv, 'full': wn == tn, 'wins': wn})
        print(f"{p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {eff:>6.3f} "
              f"{str(wo)+'/'+str(to):>6} {str(wn)+'/'+str(tn):>6} {sn:>10.3f} "
              f"{mpv:>7.3f}")
    print(f"  full recovery (5/5):  old {n_old_full}/{len(curves)}   "
          f"new {n_new_full}/{len(curves)}")
    hi = [r for r in rows_eff if r['mu_pv'] > 1.0]
    lo = [r for r in rows_eff if r['mu_pv'] <= 1.0]
    print(f"  mu/pv >  1: {len(hi):>2} curves, {sum(r['full'] for r in hi)} "
          f"of them 5/5")
    print(f"  mu/pv <= 1: {len(lo):>2} curves, {sum(r['full'] for r in lo)} "
          f"of them 5/5")
    E4_ROWS.extend(rows_eff)


# ===========================================================================
# E5 — feasibility region implied by the kernel block
# ===========================================================================
print("\n" + "=" * 78)
print("E5 — feasibility region: when CAN the planted vector be lambda_1?")
print("=" * 78)
print("""Necessary condition (Minkowski on the 2D kernel block):
    ||v_planted||  <  mu  <=  (2/sqrt3)^{1/2} * n / sqrt(eff)
    n*sqrt(2m/3 + 1)  <  (2/sqrt3)^{1/2} * n / sqrt(eff)
    =>   2m/3 + 1  <  (2/sqrt3) / eff                                  (C)
The attacker minimises m, but information theory forces
    m  >=  m_min(eff) = ln n / ln(1/eff).
Substituting m_min into (C) gives a bound on eff that depends only on n.
Note mu <= the Minkowski value, so failing (C) is a PROOF that the planted
vector is not lambda_1; passing (C) is necessary, not sufficient.""")

def feasible(n_bits, eff):
    ln_n = n_bits * math.log(2.0)
    if eff >= 1.0:
        return False
    m_min = ln_n / math.log(1.0 / eff)
    return (2.0 * m_min / 3.0 + 1.0) < (2.0 / math.sqrt(3.0)) / eff

def eff_star(n_bits):
    lo, hi = 1e-9, 0.999
    if not feasible(n_bits, lo):
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if feasible(n_bits, mid):
            lo = mid
        else:
            hi = mid
    return lo

print(f"\n{'n bits':>8} {'eff*':>10} {'m_min@eff*':>11} {'dim':>6} "
      f"{'bias bits':>10} {'K1@K2=sqrt n':>14}")
for nb in (16, 17, 32, 64, 128, 160, 192, 224, 256, 384, 521):
    es = eff_star(nb)
    if es is None:
        print(f"{nb:>8} {'none':>10}")
        continue
    ln_n = nb * math.log(2.0)
    mm = ln_n / math.log(1.0 / es)
    k1_bits = math.log2(es) + nb / 2.0
    print(f"{nb:>8} {es:>10.5f} {mm:>11.1f} {2*int(math.ceil(mm))+1:>6} "
          f"{-math.log2(es):>10.2f} {'2^%.1f' % k1_bits:>14}")

print("""
Reading: for a 256-bit group the eliminated lattice can have the planted
vector as lambda_1 only for eff below eff*, i.e. only when the GLV nonce set
is a small enough fraction of the group.  Above eff* the planted vector is
provably not lambda_1 for ANY number of signatures, because increasing m to
satisfy information theory increases ||v_planted|| faster than it helps.""")


# ===========================================================================
# E6 — can column re-weighting lengthen mu past the planted vector?
# ===========================================================================
print("\n" + "=" * 78)
print("E6 — re-weighting the k1/k2 columns (kernel-block determinant fixed)")
print("=" * 78)
print("""Stretch the k1 columns by c and shrink the k2 columns by c.  The kernel
block det = n*S_K1*S_K2 is invariant, so Minkowski's bound on mu is invariant,
but the block's SHAPE changes and mu(c) moves.  ||v_planted|| is minimised at
c = 1, so the question is whether mu(c)/pv(c) can be pushed above 1 anywhere.
If it cannot, the kernel parasite is not a scaling artifact.""")

C_GRID = [2.0 ** (j / 4.0) for j in range(-12, 13)]

for label, curve, k1_list, m in [
        ("12-bit/2557", hist[1][1], [8, 12, 16], 8),
        ("12-bit/2677", hist[2][1], [4, 6, 8], 10)]:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    print(f"\n{label}: n={n}, lam={lam}, K2={k2}, m={m}")
    print(f"{'K1':>4} {'eff':>6} {'best c':>8} {'mu/pv @c':>9} "
          f"{'mu/pv @1':>9} {'rec @c':>7} {'rec @1':>7}")
    for K1 in k1_list:
        eff = K1 * k2 / n
        best = (-1.0, None)
        base = None
        for c in C_GRID:
            SCALE_C[0] = c
            S1, _sd, S2, _sk = scales(n, K1, k2)
            mu, _w = kernel_block_mu(n, lam, S1, S2)
            pvn = math.sqrt(m * (K1 * S1) ** 2 / 3.0
                            + m * (k2 * S2) ** 2 / 3.0 + float(n) ** 2)
            r = mu / pvn
            if abs(c - 1.0) < 1e-12:
                base = r
            if r > best[0]:
                best = (r, c)
        SCALE_C[0] = best[1]
        wc, tc, _ = success_rate(curve, m, K1, SEEDS, 'elim')
        SCALE_C[0] = 1.0
        w1, t1, _ = success_rate(curve, m, K1, SEEDS, 'elim')
        print(f"{K1:>4} {eff:>6.3f} {best[1]:>8.3f} {best[0]:>9.3f} "
              f"{base:>9.3f} {str(wc)+'/'+str(tc):>7} {str(w1)+'/'+str(t1):>7}")

SCALE_C[0] = 1.0


# ===========================================================================
# E7 — adaptive re-weighting over the full 17-bit curve set
# ===========================================================================
print("\n" + "=" * 78)
print("E7 — adaptive re-weighting over all 20 curves (eff = 0.15 and 0.25)")
print("=" * 78)
print("""For each curve pick c maximising mu(c)/pv(c) -- a closed-form O(log n)
choice made BEFORE any lattice reduction -- and compare recovery at that c
against c = 1.  Rescue = a curve that gains ground; regression = one that
loses it.""")

for eff_t in (0.15, 0.25):
    print(f"\n--- eff target = {eff_t}  (m = {M_SIGS}, {len(SEEDS)} seeds) ---")
    print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'mu/pv @1':>9} "
          f"{'best c':>7} {'mu/pv @c':>9} {'rec @1':>7} {'rec @c':>7}")
    tot1 = totc = 0
    full1 = fullc = 0
    rescued = regressed = 0
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        cur = (p, b, n, lam, G)
        r1 = mu_over_pv(n, lam, k1, k2, M_SIGS, 1.0)
        cbest, rbest = best_reweight(n, lam, k1, k2, M_SIGS)
        SCALE_C[0] = 1.0
        w1, t1, _ = success_rate(cur, M_SIGS, k1, SEEDS, 'elim')
        SCALE_C[0] = cbest
        wc, tc, _ = success_rate(cur, M_SIGS, k1, SEEDS, 'elim')
        SCALE_C[0] = 1.0
        tot1 += w1; totc += wc
        full1 += (w1 == t1); fullc += (wc == tc)
        rescued += (wc > w1); regressed += (wc < w1)
        print(f"{p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {r1:>9.3f} "
              f"{cbest:>7.3f} {rbest:>9.3f} {str(w1)+'/'+str(t1):>7} "
              f"{str(wc)+'/'+str(tc):>7}")
    ns = len(curves) * len(SEEDS)
    print(f"  total recoveries: c=1 {tot1}/{ns}   best-c {totc}/{ns}")
    print(f"  full (5/5) curves: c=1 {full1}/{len(curves)}   "
          f"best-c {fullc}/{len(curves)}")
    print(f"  curves improved: {rescued}   curves regressed: {regressed}")

SCALE_C[0] = 1.0


# ===========================================================================
# Predictor scorecard: is mu/pv > 1 a sufficient condition for recovery?
# ===========================================================================
print("\n" + "=" * 78)
print("Predictor scorecard — mu/pv computed at c = 1 over all E4 rows")
print("=" * 78)
hi = [r for r in E4_ROWS if r['mu_pv'] > 1.0]
lo = [r for r in E4_ROWS if r['mu_pv'] <= 1.0]
print(f"rows total: {len(E4_ROWS)}")
print(f"  mu/pv >  1 : {len(hi):>3} rows, {sum(r['full'] for r in hi):>3} "
      f"full 5/5  ({(sum(r['full'] for r in hi)/len(hi)*100 if hi else 0):.1f}%)")
print(f"  mu/pv <= 1 : {len(lo):>3} rows, {sum(r['full'] for r in lo):>3} "
      f"full 5/5  ({(sum(r['full'] for r in lo)/len(lo)*100 if lo else 0):.1f}%)")
print("""
mu/pv > 1 says the planted vector is shorter than the kernel parasite, hence
is a plausible lambda_1.  It is a SUFFICIENT signal, not a necessary one:
recovery only needs a BDD/coset success, which happens well below 1 too.""")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
