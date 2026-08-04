"""
Shared library for GLV-HNP Thread 23 (CVP reformulation of the Phase-2 lattice).

Split out of glv_hnp_phase23_cvp.py so that glv_hnp_phase23_cvp_scaling.py can
reuse the exact same lattice construction and attack code.  See the docstring
of glv_hnp_phase23_cvp.py for the derivation of L23 and the target vector.
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ, GSO, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_20bit.py, a=0 curves)
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

def find_point(p, b, seed=42):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

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

# ---------------------------------------------------------------------------
# Signature generation (verbatim from glv_hnp_phase2_20bit.py:235)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Baseline: Kannan embedding, dim 2m+2 (verbatim, glv_hnp_phase2_20bit.py:263)
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

def recover_d_kannan(M_reduced, m, n, S_KANNAN, d_secret):
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

def attack_kannan(sigs, n, lam, k1_bound, k2_bound, d_secret,
                  use_bkz=False, bkz_beta=20):
    m = len(sigs)
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam,
                                                     k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d_kannan(reduced, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Thread 23: the CVP reformulation
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound, centered=True):
    """
    Returns (basis rows, target, scales, dim).  dim = 2m+1: the Kannan
    coordinate is dropped and d is read straight out of column m.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)

    B = [[0] * dim for _ in range(dim)]
    for i in range(m):
        B[i][i] = n * S_K1
    for i in range(m):
        B[m][i] = sigs[i]['B'] * S_K1
    B[m][m] = S_D
    for i in range(m):
        B[m + 1 + i][i] = -lam * S_K1
        B[m + 1 + i][m + 1 + i] = S_K2

    t = [0] * dim
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    if centered:
        for i in range(m):
            t[i] += (k1_bound // 2) * S_K1
            t[m + 1 + i] = (k2_bound // 2) * S_K2
        t[m] = (n // 2) * S_D

    return B, t, (S_K1, S_D, S_K2), dim

def planted_vector(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """The lattice vector corresponding to the true (k1, d, k2).  Also
    asserts it really is in L23 (the row-m multiplier must be integral)."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    v = [0] * dim
    for i in range(m):
        s = sigs[i]
        v[i] = S_K1 * (s['k1'] - s['A'])
        # membership check: x_i must be an integer
        num = (s['k1'] - s['A']) - (d_secret * s['B'] - lam * s['k2'])
        assert num % n == 0, "planted vector is not in L23"
        v[m + 1 + i] = s['k2'] * S_K2
    v[m] = d_secret * S_D
    return v

def norm2(v):
    return sum(x * x for x in v)

def dist2(v, t):
    return sum((a - b) * (a - b) for a, b in zip(v, t))

def babai_nearest_plane(A_reduced, t):
    """Babai nearest-plane via fpylll's GSO; returns the lattice vector."""
    M = GSO.Mat(A_reduced, float_type="d")
    M.update_gso()
    coeffs = M.babai(list(t))
    dim = A_reduced.nrows
    v = [0] * A_reduced.ncols
    for j in range(dim):
        c = coeffs[j]
        if c:
            for k in range(A_reduced.ncols):
                v[k] += c * A_reduced[j][k]
    return v

def attack_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret,
               centered=True, exact=False):
    """Returns (recovered_bool, d_dist2, planted_dist2, is_planted_closest)."""
    m = len(sigs)
    B, t, (S_K1, S_D, S_K2), dim = build_cvp_lattice(
        sigs, n, lam, k1_bound, k2_bound, centered)
    A = IntegerMatrix.from_matrix(B)
    LLL.reduction(A)

    if exact:
        v = list(CVP.closest_vector(A, tuple(t)))
    else:
        v = babai_nearest_plane(A, t)

    d_cand = (v[m] // S_D) % n
    vp = planted_vector(sigs, n, lam, k1_bound, k2_bound, d_secret)
    return (d_cand == d_secret, dist2(v, t), dist2(vp, t), v == vp)


def lam_star(lam, n):
    return min(lam, n - lam) / n
