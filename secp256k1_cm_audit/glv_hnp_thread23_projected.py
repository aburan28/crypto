"""
GLV-HNP Phase 2, Thread 23: make the planted vector the SHORTEST vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, Exp T5):
  The Phase-2 lattice (glv_hnp_phase2_20bit.py:262, dim 2m+2) always contains
  the vector  n*S_D*e_m .  Proof: n*(d-row) - sum_i B_i*(row_i) = (0,...,0,
  n*S_D, 0,...,0).  Its norm is n*S_D, while

      ||v_planted||^2 ~ n^2 * (2m/3 + 4/3)

  so for every m >= 1 the trivial vector is strictly shorter.  Measured on
  every curve tested on 2026-07-29: |sv[m]|/n = 1.0000 exactly, sv/pv in
  [0.337, 0.368] for successes and failures alike.  The trivial vector
  carries no information (d is only defined mod n), and no choice of S_D
  removes it because both vectors scale linearly in S_D.

  Consequence: Phase-2 recovery was never an SVP condition.  It is a
  BDD/coset condition in the projection along e_m.

This script executes the Thread 23 proposal: quotient out the e_m direction
and re-measure.  Two reformulations, both mathematically equivalent to the
original system  k1_i = A_i + B_i*d - lam*k2_i  (mod n):

  P  ("projected")  dim 2m+1.  The d-column is deleted outright.  d is not a
     bounded unknown -- it only ever contributed a term ~n^2/3 to the norm
     and created the trivial vector -- so it is carried implicitly by the
     d-row and back-solved after k1,k2 are read off:
         d = B_i^{-1} * (k1_i + lam*k2_i - A_i)  (mod n).
     Full-rank square basis (see build_projected_lattice for the derivation
     of why dropping row_0 keeps the same lattice).

  C  ("CVP")  dim 2m.  Same lattice with the Kannan row removed, solved as an
     explicit closest-vector problem with Babai nearest-plane against the
     target t = (-A_i*S_K1 ; 0).  This is the second variant named in the
     Thread 23 proposal and is a control on P: it tests whether the Kannan
     embedding itself (rather than the d-column) was the bottleneck.

FALSIFIER (stated 2026-07-29, verbatim):
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  E1  correctness -- planted vector is in P and in C (integral coordinates,
      checked exactly over Q), and P is the same lattice as the rank-deficient
      generator set (HNF fingerprint).
  E2  geometry -- sv/pv for O (original), P, C on the three historical curves.
  E3  the K1 wall, head to head O vs P vs C, on the lam*=0.340 and lam*=0.070
      curves.  This is the falsifier.
  E4  the eff sweep -- 17-bit curves at eff = 0.05 / 0.15 / 0.25, where O
      scored 19/20, 3/20, 0/20 on 2026-07-29.

Run:  python3 secp256k1_cm_audit/glv_hnp_thread23_projected.py
Deps: fpylll, sympy
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ, GSO
import sympy

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
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

# ---------------------------------------------------------------------------
# Signatures and scales (verbatim)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- unchanged from the 2026-06-15 design."""
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
# LATTICE O -- the original Phase-2 lattice, dim 2m+2  (verbatim)
# ---------------------------------------------------------------------------

def build_original_lattice(sigs, n, lam, k1_bound, k2_bound):
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

def planted_vector_original(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_original(reduced, m, n, S_KANNAN, d_secret):
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

# ---------------------------------------------------------------------------
# LATTICE P -- projected along e_m (d-column deleted), dim 2m+1
# ---------------------------------------------------------------------------
#
# Column layout:   0 .. m-1     k1 block, scale S_K1
#                  m .. 2m-1    k2 block, scale S_K2
#                  2m           Kannan column, scale S_KANNAN
#
# The generating set inherited from O (with column m deleted) is
#     R_i    = n*S_K1*e_i                       i = 0..m-1
#     R_d    = (B_i*S_K1)_i                     (the old d-row, d-entry gone)
#     R_k2,i = -lam*S_K1*e_i + S_K2*e_{m+i}     i = 0..m-1
#     R_kan  = (A_i*S_K1)_i + S_KANNAN*e_{2m}
# -- 2m+2 generators of a rank-(2m+1) lattice: the relation
#     n*R_d - sum_i B_i*R_i = 0
# is exactly the trivial vector of O, now collapsed to zero.  That is the
# whole point of the reformulation.
#
# To get a square full-rank basis, normalise R_d by a unit: pick j with
# B_j invertible mod n (n prime, B_j != 0 always), and replace R_d by
#     R_d' = (B_i * B_j^{-1} mod n)_i * S_K1        [entry j equals S_K1]
# which is d-coefficient rescaling d -> d*B_j, a unimodular relabelling.
# Then R_j is redundant:  n*R_d' has entry i equal to n*(B_i/B_j)*S_K1, a
# multiple of n*S_K1, so n*R_d' - sum_{i != j} (B_i/B_j mod n)*R_i = n*S_K1*e_j
# = R_j.  Dropping R_j therefore leaves the lattice unchanged, giving the
# square basis {R_d'} u {R_i : i != j} u {R_k2,i} u {R_kan} of rank 2m+1.
# Determinant check: S_K1 * (n*S_K1)^(m-1) * S_K2^m * S_KANNAN, i.e. the
# k1-block has index n in (S_K1*Z)^m, as it must -- one coset per value of d.

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Returns (basis, pivot_index j).  Square, rank 2m+1."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)

    j = next(i for i in range(m) if sigs[i]['B'] % n != 0)
    Bj_inv = modinv(sigs[j]['B'] % n, n)

    M = []
    row_d = [0] * dim
    for i in range(m):
        row_d[i] = (sigs[i]['B'] * Bj_inv % n) * S_K1
    M.append(row_d)
    for i in range(m):
        if i == j:
            continue
        r = [0] * dim
        r[i] = n * S_K1
        M.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        M.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KANNAN
    M.append(r)
    assert len(M) == dim
    return M, j

def projected_generators(sigs, n, lam, k1_bound, k2_bound):
    """The un-normalised (rank-deficient) 2m+2 generator set, for the HNF check."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; M.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    M.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        M.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KANNAN
    M.append(r)
    return M

def planted_vector_projected(sigs, n, k1_bound, k2_bound):
    """(k1_i*S_K1 ; k2_i*S_K2 ; S_KANNAN).  Independent of d -- that is the point."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_projected(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read k1,k2 off a Kannan row, then back-solve d from any signature."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        k1s, k2s, ok = [], [], True
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 or b % S_K2:
                ok = False; break
            k1s.append(a // S_K1); k2s.append(b // S_K2)
        if not ok:
            continue
        for i in range(m):
            Bi = sigs[i]['B'] % n
            if Bi == 0:
                continue
            d_cand = modinv(Bi, n) * (k1s[i] + lam * k2s[i] - sigs[i]['A']) % n
            if d_cand == 0:
                continue
            if d_cand == d_secret:
                return d_cand
    return None

# ---------------------------------------------------------------------------
# LATTICE C -- explicit CVP (Babai nearest plane), dim 2m
# ---------------------------------------------------------------------------
#
# Same lattice as P with the Kannan row removed and the affine part moved
# into the target:   t = (-A_i*S_K1 ; 0).
# A lattice point w = d'*R_d' + sum k2_i*R_k2,i - sum c_i*R_i has
#     w_i = (d*B_i - lam*k2_i - c_i*n)*S_K1 = (k1_i - A_i)*S_K1
# so  w - t = (k1_i*S_K1 ; k2_i*S_K2), the planted vector without its
# Kannan coordinate.

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    M, j = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    B = [row[:2 * len(sigs)] for row in M[:-1]]     # drop Kannan row and column
    return B, j

def cvp_target(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, _S_K2, _S = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -(sigs[i]['A'] % n) * S_K1
    return t

def babai_cvp(B, t):
    """Babai nearest-plane on an LLL-reduced copy of B.  Returns w in L."""
    dim = len(B)
    A = IntegerMatrix.from_matrix([list(r) for r in B])
    G = GSO.Mat(A, float_type='mpfr')
    red = LLL.Reduction(G)
    red()
    coeffs = G.babai(list(t))
    w = [0] * len(t)
    for c, i in zip(coeffs, range(dim)):
        if c:
            for k in range(len(t)):
                w[k] += c * A[i][k]
    return w

def recover_d_cvp(w, t, sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S = scales(n, k1_bound, k2_bound)
    e = [w[i] - t[i] for i in range(2 * m)]
    k1s, k2s = [], []
    for i in range(m):
        if e[i] % S_K1 or e[m + i] % S_K2:
            return None
        k1s.append(e[i] // S_K1); k2s.append(e[m + i] // S_K2)
    for i in range(m):
        Bi = sigs[i]['B'] % n
        if Bi == 0:
            continue
        d_cand = modinv(Bi, n) * (k1s[i] + lam * k2s[i] - sigs[i]['A']) % n
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Trial drivers
# ---------------------------------------------------------------------------

def reduce_rows(M, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix([list(r) for r in M])
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def trial(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False, beta=20):
    """variant in {'O','P','C'}.  Returns (ok, pv_norm, sv_norm)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if variant == 'O':
        M = build_original_lattice(sigs, n, lam, k1_bound, k2_bound)
        red = reduce_rows(M, use_bkz, beta)
        S_KAN = scales(n, k1_bound, k2_bound)[3]
        ok = recover_d_original(red, m, n, S_KAN, d_secret) is not None
        pv = norm(planted_vector_original(sigs, d_secret, n, k1_bound, k2_bound))
        sv = min(norm(r) for r in red if any(r))
    elif variant == 'P':
        M, _ = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
        red = reduce_rows(M, use_bkz, beta)
        ok = recover_d_projected(red, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
        pv = norm(planted_vector_projected(sigs, n, k1_bound, k2_bound))
        sv = min(norm(r) for r in red if any(r))
    elif variant == 'C':
        B, _ = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
        t = cvp_target(sigs, n, k1_bound, k2_bound)
        w = babai_cvp(B, t)
        ok = recover_d_cvp(w, t, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
        pv = norm(planted_vector_projected(sigs, n, k1_bound, k2_bound))
        sv = norm([w[i] - t[i] for i in range(len(t))])   # achieved distance
    else:
        raise ValueError(variant)
    return (ok, pv, sv)

def success_rate(curve, m, k1_bound, seeds, variant, use_bkz=False, beta=20):
    n = curve[2]
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, variant, use_bkz, beta)
        if res is None:
            continue
        ok, pv, sv = res
        wins += bool(ok)
        ratios.append(sv / pv if pv else float('nan'))
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr


def run_experiments():
    # ===========================================================================
    print("=" * 78)
    print("Thread 23 -- reformulate the Phase-2 lattice so the target is lambda_1")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,              p,    b, n,    lam,  K1, m
        ("8-bit/199",         211,  2, 199,  106,  2,  6),
        ("12-bit/2557",       2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  8,  10),
    ]

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), k1, m))


    # ===========================================================================
    # E1 -- correctness of the reformulation
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E1: correctness -- planted vector integral in P and C; P == HNF(generators)")
    print("-" * 78)

    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form

    def in_lattice(basis, v):
        """Exact test: solve x*B = v over Q, return x if integral else None."""
        B = Matrix(basis).T          # columns are basis vectors
        sol = B.solve(Matrix(v))     # square, full rank
        if all(c.q == 1 for c in sol):
            return [int(c) for c in sol]
        return None

    print(f"{'curve':<20} {'m':>3} {'K1':>3} {'pv in P':>8} {'pv in C':>8} "
          f"{'HNF(P)=HNF(gens)':>18}")
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, 42)
        MP, _ = build_projected_lattice(sigs, n, lam, k1, k2b)
        vP = planted_vector_projected(sigs, n, k1, k2b)
        xP = in_lattice(MP, vP)

        BC, _ = build_cvp_lattice(sigs, n, lam, k1, k2b)
        t = cvp_target(sigs, n, k1, k2b)
        vC = [vP[i] for i in range(2 * m)]
        wC = [vC[i] + t[i] for i in range(2 * m)]
        xC = in_lattice(BC, wC)

        gens = projected_generators(sigs, n, lam, k1, k2b)
        h1 = hermite_normal_form(Matrix(MP).T)
        h2 = hermite_normal_form(Matrix(gens).T)
        same = (h1 == h2)
        print(f"{label:<20} {m:>3} {k1:>3} {str(xP is not None):>8} "
              f"{str(xC is not None):>8} {str(same):>18}")

    print("\n(pv in P / pv in C = the planted vector has integral coordinates in the")
    print(" reformulated basis;  HNF equality = deleting row_j lost nothing.)")


    # ===========================================================================
    # E2 -- geometry: is the planted vector now lambda_1?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E2: sv/pv after LLL.  O = original (2m+2), P = projected (2m+1),")
    print("    C = CVP/Babai (2m, sv = achieved distance).")
    print("-" * 78)
    print(f"{'curve':<20} {'K1':>3} {'m':>3} | {'O sv/pv':>8} {'P sv/pv':>8} "
          f"{'C d/pv':>8} | {'O':>5} {'P':>5} {'C':>5}")
    e2 = []
    for label, curve, k1, m in hist:
        row = [label, k1, m]
        for var in ('O', 'P', 'C'):
            w, tot, mr = success_rate(curve, m, k1, SEEDS, var)
            row.append((w, tot, mr))
        e2.append(row)
        print(f"{label:<20} {k1:>3} {m:>3} | "
              f"{row[3][2]:>8.3f} {row[4][2]:>8.3f} {row[5][2]:>8.3f} | "
              f"{row[3][0]}/{row[3][1]:<3} {row[4][0]}/{row[4][1]:<3} {row[5][0]}/{row[5][1]:<3}")

    print("\nFALSIFIER CHECK (E2 half): does P's sv/pv rise above 1?")
    for row in e2:
        verdict = "YES" if row[4][2] >= 1.0 else "no"
        print(f"  {row[0]:<20} O={row[3][2]:.3f} -> P={row[4][2]:.3f}   {verdict}")


    # ===========================================================================
    # E3 -- the K1 wall, head to head.  THE FALSIFIER.
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E3: K1 wall.  2026-07-29 T4 measured (lattice O):")
    print("      12-bit/2557 (lam*=0.340): 5,5,5,5,5,4,1,0  for K1 = 2,3,4,6,8,12,16,24")
    print("      12-bit/2677 (lam*=0.070): 5,5,5,2,0,0,0,0")
    print("    Falsifier: does the wall move OUTWARD under P or C?")
    print("-" * 78)

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    E3_CURVES = [(lbl, cur, m) for lbl, cur, k1, m in hist if lbl.startswith("12-bit")]

    def wall_position(results, grid):
        """Largest K1 with a full 5/5, and the first K1 that scores 0."""
        last_full = None
        first_zero = None
        for k1, (w, tot) in zip(grid, results):
            if w == tot and first_zero is None:
                last_full = k1
            if w == 0 and first_zero is None:
                first_zero = k1
        return last_full, first_zero

    e3 = {}
    for label, curve, m in E3_CURVES:
        p, b, n, lam, G = curve
        ls = lam_star(lam, n)
        print(f"\n  {label}  n={n}  lam={lam}  lam*={ls:.4f}  m={m}")
        print(f"  {'variant':<8} " + " ".join(f"{k:>5}" for k in K1_GRID) +
              "   last_5/5  first_0")
        for var in ('O', 'P', 'C'):
            res = []
            for k1 in K1_GRID:
                w, tot, _ = success_rate(curve, m, k1, SEEDS, var)
                res.append((w, tot))
            lf, fz = wall_position(res, K1_GRID)
            e3[(label, var)] = (res, lf, fz)
            print(f"  {var:<8} " + " ".join(f"{w}/{t:<3}" for w, t in res) +
                  f"   {str(lf):>8}  {str(fz):>7}")
            sys.stdout.flush()

    print("\nFALSIFIER VERDICT (E3):")
    for label, curve, m in E3_CURVES:
        o_lf = e3[(label, 'O')][1]
        p_lf = e3[(label, 'P')][1]
        c_lf = e3[(label, 'C')][1]
        moved = (p_lf is not None and o_lf is not None and p_lf > o_lf) or \
                (c_lf is not None and o_lf is not None and c_lf > o_lf)
        print(f"  {label:<20} last 5/5:  O={o_lf}  P={p_lf}  C={c_lf}   "
              f"wall moved outward: {'YES' if moved else 'NO'}")


    # ===========================================================================
    # E4 -- 17-bit sweep at three bias strengths
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E4: fresh 17-bit curves at eff = K1*K2/n in {0.05, 0.15, 0.25}, m=12.")
    print("    2026-07-29 baseline (lattice O): 19/20, 3/20, 0/20 curves at 5/5.")
    print("-" * 78)

    def search_curves(lo, hi, want=12, nbins=6):
        bins = {i: [] for i in range(nbins)}
        p = int(sympy.nextprime(lo))
        while p < hi:
            if p % 3 == 1:
                eis = eisenstein_decompose(p)
                if eis is not None:
                    a_e, b_e = eis
                    for t in j0_traces(a_e, b_e):
                        n_c = p + 1 - t
                        if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                            continue
                        r = glv_roots(n_c)
                        if r is None:
                            continue
                        lam = r[0]
                        idx = min(nbins - 1, int(lam_star(lam, n_c) / (0.5 / nbins)))
                        if len(bins[idx]) >= want // nbins:
                            continue
                        cur = build_curve(p, n_c)
                        if cur is None:
                            continue
                        bins[idx].append((p, cur[1], n_c, lam, cur[4]))
            if all(len(v) >= want // nbins for v in bins.values()):
                break
            p = int(sympy.nextprime(p))
        out = []
        for i in range(nbins):
            out.extend(bins[i])
        return out

    curves17 = search_curves(2 ** 16, 2 ** 17, want=12, nbins=6)
    print(f"  found {len(curves17)} 17-bit j=0 GLV curves, "
          f"lam* in [{min(lam_star(c[3], c[2]) for c in curves17):.4f}, "
          f"{max(lam_star(c[3], c[2]) for c in curves17):.4f}]")

    M_SIGS = 12
    for eff in (0.05, 0.15, 0.25):
        print(f"\n  eff = {eff}:")
        print(f"  {'p':>7} {'n':>7} {'lam*':>7} {'K1':>5} | "
              f"{'O':>6} {'P':>6} {'C':>6}")
        tally = {'O': 0, 'P': 0, 'C': 0}
        for cur in curves17:
            p, b, n, lam, G = cur
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            line = f"  {p:>7} {n:>7} {lam_star(lam, n):>7.4f} {k1b:>5} |"
            for var in ('O', 'P', 'C'):
                w, tot, _ = success_rate(cur, M_SIGS, k1b, SEEDS, var)
                tally[var] += (w == tot)
                line += f" {w}/{tot:<4}"
            print(line)
            sys.stdout.flush()
        print(f"  {'TOTAL curves at 5/5':<32} "
              f" O={tally['O']}/{len(curves17)}  P={tally['P']}/{len(curves17)}"
              f"  C={tally['C']}/{len(curves17)}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)



if __name__ == "__main__":
    run_experiments()
