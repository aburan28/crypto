"""
GLV-HNP Phase 2, Thread 23: reformulate the Phase-2 lattice so the planted
vector IS lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, EXP T5):
  The Phase-2 lattice L (dim 2m+2, built by build_glv_lattice in
  glv_hnp_phase2_20bit.py:263) always contains the *trivial* vector

        t = n * S_D * e_m           ||t|| = n*S_D = n

  because  n*row_m - sum_i B_i*row_i = n*S_D*e_m.  Since

        ||v_planted||^2 ~ n^2 * (2m/3 + 4/3)   >  n^2   for every m >= 1,

  t is shorter than the planted vector on EVERY instance.  So the planted
  vector is never lambda_1; recovery is a BDD/coset condition (shortest
  vector among those with last coordinate +-S_KANNAN), not SVP.  t carries
  no information: d is only defined mod n, and t scales with S_D exactly as
  the planted vector does, so no rescaling removes it.

  The 2026-07-29 next-step proposal: quotient L by <t> (equivalently,
  project along e_m), and check whether the planted vector becomes
  lambda_1 of the quotient, and whether the K1 wall of EXP T4 moves.

This script does exactly that.

  Lbar  :=  L / <n*S_D*e_m>   =  image of L under "delete coordinate m".

  Kernel of the deletion map on L is exactly Z*(n*S_D*e_m) (a lattice
  vector supported only on column m must be a multiple of it), so
  rank(Lbar) = 2m+1 and det(Lbar) = det(L)/(n*S_D).

  Explicit basis of Lbar (columns: k1_0..k1_{m-1} | k2_0..k2_{m-1} | Kannan):

      P-row       : (B'_i * S_K1)_i , 0.., 0        B'_i = B_0^{-1} B_i mod n
      n-rows      : n*S_K1*e_i        for i = 1..m-1   (NOT i=0)
      lambda-rows : -lam*S_K1 at col i, S_K2 at col m+i
      Kannan row  : (A_i*S_K1)_i , 0.., S_KANNAN

  B'_0 = 1, which is what makes {P-row} u {n*S_K1*e_i, i != 0} a basis of
  the column-0..m-1 part (the naive generating set {n*S_K1*e_i} u {B_i*S_K1}
  is NOT a basis -- it has index B_0).

  Planted vector in Lbar:  (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN), obtained with
  P-row coefficient d' = k1_0 - A_0 + lam*k2_0  (an integer, = B_0*d mod n).
      ||v_planted||^2 ~ n^2 * (2m/3 + 1)      (the d column is gone)

  Recovery: from a reduced row with last coordinate +-S_KANNAN read off
  k1_0 and k2_0 and set  d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n.
  No information is lost by the projection.

Experiments:
  P1  sanity + structure: on the three historical curves, is the planted
      vector lambda_1 of Lbar?  (sv/pv ratio, where the shortest vector's
      energy sits, det/GH bookkeeping)
  P2  the EXP-T4 K1 grid, original lattice vs projected lattice -- does the
      K1 wall move outward?  (this is the falsifier stated on 2026-07-29)
  P3  T4b m-sweep at the wall (K1=8, lam*=0.07 curve), both lattices
  P4  is  pv/GH(Lbar)  a usable *a-priori* predictor of success?  Sweep
      fresh 17-bit curves x K1 grid and score the threshold.
      (T5 argued no *curve-level* invariant can work; pv/GH is an
      instance-level quantity computable before running LLL.)

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

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
# Signatures + scaling (verbatim from the Phase-2 scripts)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN)"""
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
# ORIGINAL Phase-2 lattice (dim 2m+2)
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

def recover_d_orig(reduced, m, n, S_KANNAN, d_secret):
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
# PROJECTED lattice Lbar = L / <n*S_D*e_m>   (dim 2m+1)
# columns:  0..m-1 = k1 | m..2m-1 = k2 | 2m = Kannan
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    Bp = [B0inv * sigs[i]['B'] % n for i in range(m)]     # Bp[0] == 1
    assert Bp[0] == 1

    M = []
    # P-row (the image of row_m, normalised so that its column-0 entry is 1)
    row = [0] * dim
    for i in range(m):
        row[i] = Bp[i] * S_K1
    M.append(row)
    # n-rows, skipping i = 0 (its pivot is taken by the P-row)
    for i in range(1, m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    # lambda-rows
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    # Kannan row
    row = [0] * dim
    for i in range(m):
        row[i] = sigs[i]['A'] * S_K1
    row[dim - 1] = S_KANNAN
    M.append(row)
    assert len(M) == dim
    return M

def planted_vector_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_proj(reduced, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n, read off a Kannan row."""
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        a0, b0 = sign * row[0], sign * row[m]
        if a0 % S_K1 or b0 % S_K2: continue
        k1_0, k2_0 = a0 // S_K1, b0 // S_K2
        d_cand = (k1_0 + lam * k2_0 - sigs[0]['A']) * B0inv % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Determinant / Gaussian-heuristic bookkeeping (closed form, no big linalg)
# ---------------------------------------------------------------------------

def log_det_orig(m, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))

def log_det_proj(m, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (math.log(S_K1) + (m - 1) * math.log(n * S_K1)
            + m * math.log(S_K2) + math.log(S_KAN))

def gaussian_heuristic(log_det, dim):
    return math.exp(log_det / dim) * math.sqrt(dim / (2 * math.pi * math.e))

def planted_norm_expected(m, n, k1_bound, k2_bound, projected):
    """E||v_planted|| with k1~U[0,K1), k2~U[0,K2), d~U[0,n); E[x^2]=X^2/3."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    tot = (m * (k1_bound * S_K1) ** 2 / 3.0
           + m * (k2_bound * S_K2) ** 2 / 3.0
           + S_KAN ** 2)
    if not projected:
        tot += (n * S_D) ** 2 / 3.0
    return math.sqrt(tot)

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_one(curve, m, d_secret, k1_bound, seed, projected, use_bkz=False,
            bkz_beta=20):
    """Returns dict with ok / planted norm / shortest reduced norm / sv row."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if projected:
        M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1
        pv = planted_vector_proj(sigs, n, k1_bound, k2_bound)
    else:
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 2
        pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if projected:
        ok = recover_d_proj(reduced, sigs, m, n, lam, k1_bound, k2_bound,
                            d_secret) is not None
    else:
        ok = recover_d_orig(reduced, m, n, S_KAN, d_secret) is not None
    nz = [r for r in reduced if any(r)]
    sv = min(nz, key=norm)
    return {'ok': ok, 'pn': norm(pv), 'sn': norm(sv), 'sv': sv,
            'dim': dim, 'sigs': sigs}

def success_rate(curve, m, k1_bound, seeds, projected, use_bkz=False,
                 bkz_beta=20):
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        r = run_one(curve, m, d_trial, k1_bound, seed, projected,
                    use_bkz=use_bkz, bkz_beta=bkz_beta)
        if r is None:
            continue
        wins += bool(r['ok'])
        ratios.append(r['sn'] / r['pn'])
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================
print("=" * 78)
print("Thread 23 — projecting out the trivial vector n*S_D*e_m  (GLV-HNP Phase 2)")
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
# EXPERIMENT P0 — correctness of the projected construction
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P0: correctness of Lbar")
print("-" * 78)
print("Checks, per curve: (a) the planted vector really lies in Lbar (solve for")
print("the integer coefficient vector and verify), (b) det(L)/det(Lbar) = n*S_D,")
print("(c) the recovery map from a Kannan row returns d exactly.\n")

print(f"{'curve':<18} {'m':>3} {'planted in Lbar':>16} {'det ratio / (n*S_D)':>21} "
      f"{'d-map ok':>9}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    M = build_projected_lattice(sigs, n, lam, k1, k2b)
    pv = planted_vector_proj(sigs, n, k1, k2b)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    # explicit coefficient vector: Kannan(1) + sum k2_i*lambda-row_i
    #                              + d'*P-row + c_i*n-row_i
    dprime = sigs[0]['k1'] - sigs[0]['A'] + lam * sigs[0]['k2']
    coeff = [0] * (2 * m + 1)
    coeff[0] = dprime                                  # P-row
    for i in range(m):
        coeff[m + i] = sigs[i]['k2']                   # lambda-rows
    coeff[2 * m] = 1                                   # Kannan row
    # solve for the n-row coefficients from columns 1..m-1
    B0inv = modinv(sigs[0]['B'] % n, n)
    for i in range(1, m):
        Bp_i = B0inv * sigs[i]['B'] % n
        want = sigs[i]['k1']
        have = dprime * Bp_i - lam * sigs[i]['k2'] + sigs[i]['A']
        diff = want - have
        assert diff % n == 0, "column not congruent mod n"
        coeff[i] = diff // n                           # n-row i
    combo = [sum(coeff[r] * M[r][c] for r in range(2 * m + 1))
             for c in range(2 * m + 1)]
    in_lattice = (combo == pv)
    ratio = math.exp(log_det_orig(m, n, k1, k2b) - log_det_proj(m, n, k1, k2b))
    det_ok = abs(ratio / (n * S_D) - 1.0) < 1e-9
    # d-map on the exact planted vector
    dback = recover_d_proj([pv], sigs, m, n, lam, k1, k2b, d0)
    print(f"{label:<18} {m:>3} {str(in_lattice):>16} "
          f"{ratio/(n*S_D):>21.9f} {str(dback == d0):>9}")


# ===========================================================================
# EXPERIMENT P1 — is the planted vector lambda_1 of Lbar?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P1: sv/pv in L (original) vs Lbar (projected)")
print("-" * 78)
print("T5 (2026-07-29) measured sv/pv in [0.34, 0.61] in L, with 100% of the")
print("shortest vector's energy in the d column.  If the projection works,")
print("sv/pv should jump to ~1.00 in Lbar.  'kan' = |sv_last|/S_KANNAN.\n")

print(f"{'curve':<18} {'K1':>3} {'m':>3} | {'sv/pv L':>8} {'kan L':>6} "
      f"{'ok L':>5} | {'sv/pv Lbar':>10} {'kan Lbar':>9} {'ok Lbar':>8} "
      f"{'pv=sv?':>7}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    ro = run_one(curve, m, d0, k1, 42, projected=False)
    rp = run_one(curve, m, d0, k1, 42, projected=True)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, math.isqrt(n) + 1)
    kan_o = abs(ro['sv'][-1]) / S_KAN
    kan_p = abs(rp['sv'][-1]) / S_KAN
    same = abs(rp['sn'] / rp['pn'] - 1.0) < 1e-12
    print(f"{label:<18} {k1:>3} {m:>3} | {ro['sn']/ro['pn']:>8.3f} "
          f"{kan_o:>6.2f} {str(ro['ok']):>5} | {rp['sn']/rp['pn']:>10.3f} "
          f"{kan_p:>9.2f} {str(rp['ok']):>8} {str(same):>7}")


# ===========================================================================
# EXPERIMENT P2 — the falsifier: does the K1 wall of T4 move?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P2: EXP-T4 K1 grid, original L vs projected Lbar  (m=12, 5 seeds)")
print("-" * 78)
print("2026-07-29 falsifier: 'if the K1 wall on the lam*=0.07 curve (currently")
print("K1 ~ 4-6) moves outward, the reformulation is a real improvement; if it")
print("stays, the wall is information-theoretic and Phase 2 is at its ceiling.'\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
print(f"{'curve':<18} {'lam*':>7} {'lattice':>8} "
      + " ".join(f"K1={k:<4}" for k in K1_GRID))
p2_rows = []
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue
    for projected in (False, True):
        cells = []
        for k1 in K1_GRID:
            w, t, _ = success_rate(curve, 12, k1, SEEDS, projected)
            cells.append(f"{w}/{t}   ")
            p2_rows.append({'label': label, 'proj': projected, 'k1': k1,
                            'wins': w})
        tag = "Lbar" if projected else "L"
        print(f"{label:<18} {lam_star(lam,n):>7.4f} {tag:>8} " + " ".join(cells))

k2ref = math.isqrt(2647) + 1
print("\n(eff = K1*K2/n with K2=%d, n=2647: " % k2ref
      + ", ".join(f"K1={k}:{k*k2ref/2647:.3f}" for k in K1_GRID) + ")")

# wall summary
print("\nWall (largest K1 with >= 3/5 recoveries):")
for label in sorted({r['label'] for r in p2_rows}):
    for projected in (False, True):
        good = [r['k1'] for r in p2_rows
                if r['label'] == label and r['proj'] == projected
                and r['wins'] >= 3]
        tag = "Lbar" if projected else "L   "
        print(f"  {label:<18} {tag}  K1_max = {max(good) if good else '-'}")


# ===========================================================================
# EXPERIMENT P3 — T4b: does more data help in Lbar? (K1=8, lam*=0.07 curve)
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P3: m-sweep at the wall (K1=8, 12-bit/2677, lam*=0.070)")
print("-" * 78)
print("T4b in L: m = 8/12/16/24/32 -> 0,0,1,0,1 out of 5.  In Lbar the planted")
print("norm grows as sqrt(2m/3+1) while GH(Lbar) grows too -- so a genuine")
print("optimum in m is predicted, not monotone improvement.\n")

fail_curve = [c for lbl, c, _, _ in hist if c[2] == 2647][0]
M_GRID = [4, 6, 8, 12, 16, 24, 32]
print(f"{'lattice':>8} " + " ".join(f"m={m:<5}" for m in M_GRID)
      + "   | pv/GH(Lbar) per m")
for projected in (False, True):
    cells, ratios = [], []
    for m_try in M_GRID:
        w, t, _ = success_rate(fail_curve, m_try, 8, SEEDS, projected)
        cells.append(f"{w}/{t}  ")
    tag = "Lbar" if projected else "L"
    print(f"{tag:>8} " + " ".join(cells))
ratios = []
for m_try in M_GRID:
    n = fail_curve[2]
    k2b = math.isqrt(n) + 1
    gh = gaussian_heuristic(log_det_proj(m_try, n, 8, k2b), 2 * m_try + 1)
    pv = planted_norm_expected(m_try, n, 8, k2b, projected=True)
    ratios.append(pv / gh)
print("  pv/GH  " + " ".join(f"{r:<7.3f}" for r in ratios))


# ===========================================================================
# EXPERIMENT P4 — is pv/GH(Lbar) an a-priori predictor?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P4: pv/GH(Lbar) as an a-priori success predictor (17-bit curves)")
print("-" * 78)
print("pv/GH is computable BEFORE running LLL: it needs only (n, K1, K2, m).")
print("T5 argued no *curve-level* invariant can predict success; pv/GH is an")
print("*instance-level* (parameter-level) quantity, so it is not excluded.\n")

def search_curves(lo, hi, want=12):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    out.append((p, cur[1], nc, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

curves17 = search_curves(2 ** 16, 2 ** 17, want=10)
print(f"collected {len(curves17)} 17-bit j=0 GLV curves\n")

rows = []
for (p, b, n, lam, G) in curves17:
    k2b = math.isqrt(n) + 1
    for k1 in (2, 4, 8, 16, 32):
        for m_try in (8, 16):
            w, t, _ = success_rate((p, b, n, lam, G), m_try, k1, SEEDS[:3],
                                   projected=True)
            gh = gaussian_heuristic(log_det_proj(m_try, n, k1, k2b),
                                    2 * m_try + 1)
            pv = planted_norm_expected(m_try, n, k1, k2b, projected=True)
            rows.append({'n': n, 'lam_star': lam_star(lam, n), 'k1': k1,
                         'm': m_try, 'wins': w, 'tot': t, 'ok': w >= 2,
                         'pvgh': pv / gh, 'eff': k1 * k2b / n})

succ = [r['pvgh'] for r in rows if r['ok']]
fail = [r['pvgh'] for r in rows if not r['ok']]
print(f"{len(rows)} (curve, K1, m) cells: {len(succ)} success, {len(fail)} failure")
if succ and fail:
    print(f"  pv/GH  success range [{min(succ):.3f}, {max(succ):.3f}]")
    print(f"  pv/GH  failure range [{min(fail):.3f}, {max(fail):.3f}]")
    print(f"  separable (no overlap): {max(succ) < min(fail)}")
    # best threshold accuracy + AUC
    vals = sorted({r['pvgh'] for r in rows})
    best = (0, None)
    for v in vals:
        acc = sum(1 for r in rows if (r['pvgh'] < v) == r['ok'])
        if acc > best[0]:
            best = (acc, v)
    base = max(len(succ), len(fail))
    auc = sum(1.0 if s < f else 0.5 if s == f else 0.0
              for s in succ for f in fail) / (len(succ) * len(fail))
    print(f"  best threshold  pv/GH < {best[1]:.4f} -> {best[0]}/{len(rows)} "
          f"= {best[0]/len(rows):.1%}  (majority baseline "
          f"{base}/{len(rows)} = {base/len(rows):.1%})")
    print(f"  AUC(pv/GH, lower = success) = {auc:.3f}")
    # comparison predictors
    for key, direction in (('eff', 'lower'), ('lam_star', 'higher')):
        s = [r[key] for r in rows if r['ok']]
        f = [r[key] for r in rows if not r['ok']]
        a = sum(1.0 if (x < y if direction == 'lower' else x > y)
                else 0.5 if x == y else 0.0
                for x in s for y in f) / (len(s) * len(f))
        print(f"  AUC({key}, {direction} = success) = {a:.3f}")

print("\nper-cell detail (first 40):")
print(f"{'n':>7} {'lam*':>6} {'K1':>3} {'m':>3} {'eff':>6} {'pv/GH':>7} {'wins':>5}")
for r in rows[:40]:
    print(f"{r['n']:>7} {r['lam_star']:>6.3f} {r['k1']:>3} {r['m']:>3} "
          f"{r['eff']:>6.3f} {r['pvgh']:>7.3f} {str(r['wins'])+'/'+str(r['tot']):>5}")


# ===========================================================================
# EXPERIMENT P5 — the homogeneous sublattice L0 is the real obstruction
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P5: BDD formulation — L0 = {v in Lbar : Kannan coord 0}, rank 2m")
print("-" * 78)
print("P1 shows the shortest vector of Lbar still has Kannan coordinate 0, so")
print("recovery is BDD with target distance pv in the rank-2m lattice L0.")
print("L0 splits into m identical 2-D lambda blocks <(n*S_K1,0),(-lam*S_K1,S_K2)>")
print("plus the (long) P-row, so lambda_1(L0) should equal the exact 2-D block")
print("shortest vector mu.  Unique BDD decoding needs pv < lambda_1(L0)/2,")
print("i.e. rho = mu/pv > 2.  This is the controlled test T2 could not make:")
print("same n-size, same K2, same m, ONLY lam and K1 vary.\n")

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
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

def block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1])

def lambda1_L0(sigs, n, lam, k1_bound, k2_bound):
    """LLL on the rank-2m homogeneous sublattice (drop the Kannan row)."""
    m = len(sigs)
    full = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    sub = [r[:-1] for r in full[:-1]]
    A = IntegerMatrix.from_matrix(sub)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(2 * m)] for i in range(2 * m)]
    return min(norm(r) for r in rows if any(r))

print(f"{'curve':<18} {'lam*':>7} {'K1':>3} {'mu':>10} {'l1(L0)':>10} "
      f"{'mu=l1':>6} {'pv':>10} {'rho=mu/pv':>10} {'wins':>6}")
p5_rows = []
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
        mu = block_mu(n, lam, S_K1, S_K2)
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sg = gen_signatures(G, d0, 12, n, lam, p, k1, k2b, 42)
        l1 = lambda1_L0(sg, n, lam, k1, k2b)
        pv = planted_norm_expected(12, n, k1, k2b, projected=True)
        w = [r['wins'] for r in p2_rows
             if r['label'] == label and r['proj'] and r['k1'] == k1][0]
        p5_rows.append({'label': label, 'lam_star': lam_star(lam, n), 'k1': k1,
                        'rho': mu / pv, 'ok': w >= 3, 'wins': w})
        print(f"{label:<18} {lam_star(lam,n):>7.4f} {k1:>3} {mu:>10.1f} "
              f"{l1:>10.1f} {str(abs(mu-l1)<1.0):>6} {pv:>10.1f} "
              f"{mu/pv:>10.3f} {str(w)+'/5':>6}")

s = [r['rho'] for r in p5_rows if r['ok']]
f = [r['rho'] for r in p5_rows if not r['ok']]
if s and f:
    print(f"\n  rho  success range [{min(s):.3f}, {max(s):.3f}]")
    print(f"  rho  failure range [{min(f):.3f}, {max(f):.3f}]")
    print(f"  separable by a single rho threshold: {min(s) > max(f)}")
    vals = sorted({r['rho'] for r in p5_rows})
    best = (0, None)
    for v in vals:
        acc = sum(1 for r in p5_rows if (r['rho'] >= v) == r['ok'])
        if acc > best[0]:
            best = (acc, v)
    base = max(len(s), len(f))
    print(f"  best threshold  rho >= {best[1]:.4f} -> {best[0]}/{len(p5_rows)}"
          f"  (majority baseline {base}/{len(p5_rows)})")
    print(f"  unique-BDD prediction rho > 2 would give: "
          f"{sum(1 for r in p5_rows if (r['rho'] > 2.0) == r['ok'])}"
          f"/{len(p5_rows)}")


# ===========================================================================
# EXPERIMENT P6 — fully controlled slice: fix (n-size, K1, K2, m), vary lam
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P6: the K1=32 slice of P4 — 10 curves, identical eff and pv/GH")
print("-" * 78)
print("In P4 the curves n=65287 (recovers 3/3) and n=65053, 65119, 66109 (0/3)")
print("have eff and pv/GH agreeing to 3 decimals.  Neither eff nor pv/GH can")
print("explain the split; lam* is at chance.  Here mu = lambda_1 of the 2-D")
print("lambda block is the only remaining candidate, tested on a slice where")
print("EVERYTHING else is held fixed.\n")

print(f"{'n':>7} {'lam':>7} {'lam*':>6} {'mu':>10} {'pv':>10} {'rho':>7} "
      f"{'m=8':>5} {'m=16':>5}")
p6 = []
for (p, b, n, lam, G) in curves17:
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, 32, k2b)
    mu = block_mu(n, lam, S_K1, S_K2)
    cells = {}
    for m_try in (8, 16):
        w = [r['wins'] for r in rows
             if r['n'] == n and r['k1'] == 32 and r['m'] == m_try][0]
        cells[m_try] = w
    pv = planted_norm_expected(8, n, 32, k2b, projected=True)
    p6.append({'n': n, 'rho': mu / pv, 'ok': cells[8] >= 2,
               'lam_star': lam_star(lam, n)})
    print(f"{n:>7} {lam:>7} {lam_star(lam,n):>6.3f} {mu:>10.1f} {pv:>10.1f} "
          f"{mu/pv:>7.3f} {str(cells[8])+'/3':>5} {str(cells[16])+'/3':>5}")

s = [r['rho'] for r in p6 if r['ok']]
f = [r['rho'] for r in p6 if not r['ok']]
if s and f:
    auc = sum(1.0 if x > y else 0.5 if x == y else 0.0
              for x in s for y in f) / (len(s) * len(f))
    print(f"\n  rho  success range [{min(s):.3f}, {max(s):.3f}]  (n={len(s)})")
    print(f"  rho  failure range [{min(f):.3f}, {max(f):.3f}]  (n={len(f)})")
    print(f"  separable by a single rho threshold: {min(s) > max(f)}")
    print(f"  AUC(rho, higher = success) = {auc:.3f}")
    sl = [r['lam_star'] for r in p6 if r['ok']]
    fl = [r['lam_star'] for r in p6 if not r['ok']]
    aucl = sum(1.0 if x > y else 0.5 if x == y else 0.0
               for x in sl for y in fl) / (len(sl) * len(fl))
    print(f"  AUC(lam*, higher = success) = {aucl:.3f}")
else:
    print("\n  degenerate slice (all successes or all failures) — uninformative")


# ===========================================================================
# EXPERIMENT P7 — confirm P6 on a larger slice, and pin down the sign flip
# ===========================================================================
print("\n" + "-" * 78)
print("EXP P7: (a) enlarge the controlled slice; (b) document the sign flip")
print("-" * 78)
print("P6 (10 curves, K1=32, m=16) separated perfectly with SMALL mu = success.")
print("That is the OPPOSITE direction to H20 (2026-07-29), which predicted")
print("large rho = mu/pv = success and was falsified.  10 points can be luck,")
print("so (a) re-runs the slice on ~30 fresh curves.  (b) contrasts it with the")
print("within-curve K1 sweep, where the direction is reversed.")
print("nu = mu/n is the scale-free form: mu = n * min_(a,b) sqrt((a/K1)^2+(b/K2)^2)")
print("over a = -lam*b mod n, so nu is directly comparable across K1 and n.\n")

curves17b = search_curves(2 ** 17, 2 ** 18, want=30)
print(f"(a) collected {len(curves17b)} fresh curves in [2^17, 2^18)\n")

K1_FIX, M_FIX = 32, 16
print(f"{'n':>8} {'lam*':>6} {'nu=mu/n':>8} {'wins':>6}")
p7 = []
for (p, b, n, lam, G) in curves17b:
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, K1_FIX, k2b)
    nu = block_mu(n, lam, S_K1, S_K2) / n
    w, t, _ = success_rate((p, b, n, lam, G), M_FIX, K1_FIX, SEEDS[:3],
                           projected=True)
    p7.append({'n': n, 'nu': nu, 'ok': w >= 2, 'wins': w,
               'lam_star': lam_star(lam, n)})
    print(f"{n:>8} {lam_star(lam,n):>6.3f} {nu:>8.3f} {str(w)+'/'+str(t):>6}")

s = [r['nu'] for r in p7 if r['ok']]
f = [r['nu'] for r in p7 if not r['ok']]
if s and f:
    auc = sum(1.0 if x < y else 0.5 if x == y else 0.0
              for x in s for y in f) / (len(s) * len(f))
    sl = [r['lam_star'] for r in p7 if r['ok']]
    fl = [r['lam_star'] for r in p7 if not r['ok']]
    aucl = sum(1.0 if x > y else 0.5 if x == y else 0.0
               for x in sl for y in fl) / (len(sl) * len(fl))
    print(f"\n  nu   success range [{min(s):.3f}, {max(s):.3f}]  (n={len(s)})")
    print(f"  nu   failure range [{min(f):.3f}, {max(f):.3f}]  (n={len(f)})")
    print(f"  separable by a single nu threshold: {max(s) < min(f)}")
    print(f"  AUC(nu, LOWER = success)   = {auc:.3f}")
    print(f"  AUC(lam*, higher=success)  = {aucl:.3f}")
    vals = sorted({r['nu'] for r in p7})
    best = (0, None)
    for v in vals:
        acc = sum(1 for r in p7 if (r['nu'] < v) == r['ok'])
        if acc > best[0]:
            best = (acc, v)
    base = max(len(s), len(f))
    print(f"  best threshold nu < {best[1]:.3f} -> {best[0]}/{len(p7)}"
          f"  (majority baseline {base}/{len(p7)})")
else:
    print("\n  degenerate slice — uninformative")

print("\n(b) within-curve K1 sweep on the historical pair (nu vs success, m=12):")
print(f"{'curve':<18} " + " ".join(f"K1={k:<5}" for k in K1_GRID))
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue
    k2b = math.isqrt(n) + 1
    nus, wins = [], []
    for k1 in K1_GRID:
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
        nus.append(block_mu(n, lam, S_K1, S_K2) / n)
        wins.append([r['wins'] for r in p2_rows
                     if r['label'] == label and r['proj'] and r['k1'] == k1][0])
    print(f"{label+' nu':<18} " + " ".join(f"{v:<7.3f}" for v in nus))
    print(f"{label+' win':<18} " + " ".join(f"{w}/5    " for w in wins))
print("\n  Within a curve, success DECREASES as nu decreases (large nu = success).")
print("  Across curves at fixed K1, success INCREASES as nu decreases.")
print("  => nu is not a global predictor with a fixed sign; the two sweeps move")
print("     along different directions of the (K1, lam) parameter plane.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
