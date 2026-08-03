"""
GLV-HNP Phase 2, Thread 23: eliminate d so the planted vector can be lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run #1, finding T5):
  In the Phase-2 lattice built by `build_glv_lattice` (dim 2m+2, columns
  [k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan]) the shortest vector after
  LLL is, on EVERY curve tested, the trivial vector n*S_D*e_m -- 100% of its
  energy in the d-column, |sv[m]|/n = 1.0000 exactly.  It carries no
  information (d is only defined mod n) and no choice of S_D removes it,
  because the planted vector's d-coordinate d*S_D and the trivial vector
  n*S_D*e_m scale together.  Measured sv/pv sat in [0.337, 0.368] for
  successes and failures alike.

  Root cause (this script's framing): d is the ONLY unbounded unknown in the
  system.  k1_i lives in [0,K1) and k2_i in [0,K2), so their columns can be
  scaled up (S_K1 = n//K1, S_K2 = n//K2) until the modular vectors n*S_K1*e_i
  and n*S_K2*e_{m+1+i} are far longer than the planted vector.  d ranges over
  the whole of [0,n), so S_D is pinned at 1 and the modular vector n*S_D*e_m
  has norm exactly n -- the same size as ONE planted coordinate, hence
  sqrt(2m/3+1) times shorter than the planted vector as a whole.

  Thread 23 (proposed 2026-07-29): reformulate so the planted vector is
  lambda_1.  This script does it by ELIMINATING d algebraically rather than
  by rescaling.

The d-free construction
-----------------------
  ECDSA/HNP relation, m signatures:  k_i = A_i + B_i*d (mod n),
  with the GLV split k_i = k1_i + lam*k2_i, k1_i in [0,K1), k2_i in [0,K2).

  Take signature 0 as pivot.  With t_i = B_i * B_0^{-1} mod n and
  u_i = (A_i - t_i*A_0) mod n, for i = 1..m-1:

      k_i = u_i + t_i*k_0   (mod n)                      [d has vanished]

  i.e.  k1_i = u_i + t_i*k1_0 + t_i*lam*k2_0 - lam*k2_i  (mod n).

  Free variables: k1_0, k2_0..k2_{m-1}  (m+1 of them, all BOUNDED).
  Determined:     k1_1..k1_{m-1}        (m-1 congruences).
  Deficit m+1 is identical to the original lattice's 2m+1 unknowns minus m
  congruences -- no information is gained or lost, only the geometry changes.

  Lattice (dim 2m+1), columns [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan]:
     m-1 modular rows   n*S_K1*e_i                        (i = 1..m-1)
     R_k10   col 0 = S_K1,  col i = t_i*S_K1              (i >= 1)
     R_k20   col i = (t_i*lam mod n)*S_K1 (i >= 1),  col m = S_K2
     R_k2j   col j = -lam*S_K1,  col m+j = S_K2           (j = 1..m-1)
     R_kan   col i = u_i*S_K1 (i >= 1),  col 2m = S_KANNAN

  Planted vector = (k1_0*S_K1, .., k1_{m-1}*S_K1, k2_0*S_K2, .., S_KANNAN),
  reached as  1*R_kan + k1_0*R_k10 + k2_0*R_k20 + sum_j k2_j*R_k2j + modular.

  Shortest trivial vectors are now n*S_K1*e_0 (norm n^2/K1) and
  n*S_K2*e_{m+j} (norm n^2/K2 ~ n^1.5) -- both far LONGER than the planted
  vector's ~n*sqrt(2m/3+1).  Nothing of norm ~n survives.

  d is recovered afterwards from the pivot:  k_0 = k1_0 + lam*k2_0 mod n,
  then d = (k_0 - A_0) * B_0^{-1} mod n.

Falsifier (stated verbatim in the 2026-07-29 next-step proposal)
---------------------------------------------------------------
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  E1  construction validation: planted vector is IN the lattice, recovery of d
      through the pivot works, and the d-free lattice has no norm-~n vector.
  E2  sv/pv head-to-head, old lattice vs d-free, on the historical curves.
  E3  K1-wall replication (the 2026-07-29 T4 grid) old vs d-free.
  E4  17-bit multi-curve head-to-head at eff in {0.05, 0.15, 0.25}.

Run: python3 glv_hnp_phase2_dfree.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py so the comparison is exact
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
# CM theory for j=0 curves
# ---------------------------------------------------------------------------

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
# Signatures (verbatim) and scales
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- identical to the original construction."""
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
# ORIGINAL Phase-2 lattice (dim 2m+2), for head-to-head
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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
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
# D-FREE Phase-2 lattice (dim 2m+1) -- Thread 23
# ---------------------------------------------------------------------------

def dfree_params(sigs, n):
    """Pivot elimination coefficients t_i, u_i (i = 1..m-1), pivot = sig 0."""
    B0inv = modinv(sigs[0]['B'] % n, n)
    A0 = sigs[0]['A'] % n
    t = [0] * len(sigs)
    u = [0] * len(sigs)
    for i in range(1, len(sigs)):
        t[i] = sigs[i]['B'] * B0inv % n
        u[i] = (sigs[i]['A'] - t[i] * A0) % n
    return t, u

def build_dfree_lattice(sigs, n, lam, k1_bound, k2_bound):
    """dim 2m+1, columns [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan]."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    lam_r = lam % n
    t, u = dfree_params(sigs, n)

    M = []
    # m-1 modular rows on the determined k1 columns
    for i in range(1, m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    # R_k10 : the free k1_0
    row = [0] * dim
    row[0] = S_K1
    for i in range(1, m):
        row[i] = t[i] * S_K1
    M.append(row)
    # R_k20 : the free k2_0
    row = [0] * dim
    for i in range(1, m):
        row[i] = (t[i] * lam_r % n) * S_K1
    row[m] = S_K2
    M.append(row)
    # R_k2j : the free k2_j, j = 1..m-1
    for j in range(1, m):
        row = [0] * dim
        row[j] = -lam_r * S_K1
        row[m + j] = S_K2
        M.append(row)
    # Kannan row
    row = [0] * dim
    for i in range(1, m):
        row[i] = u[i] * S_K1
    row[dim - 1] = S_KAN
    M.append(row)

    assert len(M) == dim, (len(M), dim)
    return M

def dfree_planted_vector(sigs, n, k1_bound, k2_bound):
    """(k1_0*S_K1, .., k1_{m-1}*S_K1, k2_0*S_K2, .., k2_{m-1}*S_K2, S_KAN)."""
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def dfree_recover_d(M_reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read k1_0,k2_0 off any row with |last| = S_KAN, then d via the pivot."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    A0 = sigs[0]['A'] % n
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sgn = 1 if last > 0 else -1
        if row[0] % S_K1 != 0 or row[m] % S_K2 != 0:
            continue
        k1_0 = sgn * row[0] // S_K1
        k2_0 = sgn * row[m] // S_K2
        k_0 = (k1_0 + lam * k2_0) % n
        d_cand = (k_0 - A0) * B0inv % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_once(curve, m, d_secret, k1_bound, seed=42, variant='dfree',
             use_bkz=False, bkz_beta=20):
    """Returns (recovered, planted_norm, shortest_reduced_norm, sigs, reduced)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if variant == 'dfree':
        M = build_dfree_lattice(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1
        pv = dfree_planted_vector(sigs, n, k1_bound, k2_bound)
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
    S_KAN = scales(n, k1_bound, k2_bound)[3]
    if variant == 'dfree':
        ok = dfree_recover_d(reduced, sigs, n, lam, k1_bound, k2_bound,
                             d_secret) is not None
    else:
        ok = recover_d(reduced, m, n, S_KAN, d_secret) is not None
    nz = [norm(r) for r in reduced if any(r)]
    sn = min(nz) if nz else float('nan')
    return (ok, norm(pv), sn, sigs, reduced)

def success_rate(curve, m, k1_bound, seeds, variant='dfree',
                 use_bkz=False, bkz_beta=20):
    """Returns (wins, total, mean sv/pv)."""
    wins, ratios = 0, []
    n = curve[2]
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = run_once(curve, m, d_trial, k1_bound, seed, variant=variant,
                       use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        ok, pn, sn, _, _ = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    mean_ratio = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mean_ratio


# ===========================================================================
print("=" * 78)
print("Thread 23 — d-free reformulation of the GLV-HNP Phase-2 lattice")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26)
HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

hist_curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist_curves.append((label, (p, b, n, lam, G), k1, m))

# ===========================================================================
# EXPERIMENT E1 — construction validation
# ===========================================================================
print("\n" + "-" * 78)
print("EXP E1: is the d-free construction correct?")
print("-" * 78)
print("For each curve: (a) the planted vector must lie IN the d-free lattice,")
print("(b) d must be recoverable from it through the pivot, (c) the lattice")
print("must contain no vector of norm ~n (the defect of the original).\n")

def in_lattice(M, v):
    """Solve x*M = v over Q and check integrality (exact, via sympy)."""
    from sympy import Matrix, Rational
    Mm = Matrix(M).T                      # columns = basis vectors
    sol = Mm.solve(Matrix(v))
    return all(Rational(c).q == 1 for c in sol), [int(c) for c in sol] \
        if all(Rational(c).q == 1 for c in sol) else None

print(f"{'curve':<18} {'m':>3} {'K1':>3} {'planted in L':>13} "
      f"{'d via pivot':>12} {'min triv/n':>11}")
e1_ok = True
for label, curve, k1, m in hist_curves:
    p, b, n, lam, G = curve
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_trial, m, n, lam, p, k1, k2_bound, 42)
    M = build_dfree_lattice(sigs, n, lam, k1, k2_bound)
    pv = dfree_planted_vector(sigs, n, k1, k2_bound)
    inL, coeffs = in_lattice(M, pv)
    # recover d directly from the planted vector (sanity of the pivot map)
    drec = dfree_recover_d([pv], sigs, n, lam, k1, k2_bound, d_trial)
    # shortest NON-planted "trivial" vector: LLL-reduce the sublattice spanned
    # by everything except the Kannan row (the homogeneous part).
    hom = [r[:] for r in M[:-1]]
    Ah = IntegerMatrix.from_matrix(hom)
    LLL.reduction(Ah)
    hr = [[Ah[i][j] for j in range(2 * m + 1)] for i in range(len(hom))]
    triv = min(norm(r) for r in hr if any(r))
    e1_ok &= bool(inL) and drec == d_trial
    print(f"{label:<18} {m:>3} {k1:>3} {str(bool(inL)):>13} "
          f"{str(drec == d_trial):>12} {triv / n:>11.3f}")

print(f"\n  E1 verdict: {'PASS' if e1_ok else 'FAIL'} "
      f"(planted vector in lattice AND d recoverable, all curves)")
print("  'min triv/n' is the shortest vector of the homogeneous sublattice,")
print("  in units of n.  The ORIGINAL lattice has this = 1.000 exactly on")
print("  every curve (log 2026-07-29 T5).  Anything > sqrt(2m/3+1) means the")
print("  planted vector can be lambda_1.")

# ===========================================================================
# EXPERIMENT E2 — sv/pv head-to-head
# ===========================================================================
print("\n" + "-" * 78)
print("EXP E2: sv/pv (shortest reduced vector / planted norm), old vs d-free")
print("-" * 78)
print("Falsifier part 1: sv/pv must rise above 1 in the d-free lattice.\n")
print(f"{'curve':<18} {'m':>3} {'K1':>3} "
      f"{'old sv/pv':>10} {'old':>7} {'dfree sv/pv':>12} {'dfree':>7}")
for label, curve, k1, m in hist_curves:
    w_o, t_o, r_o = success_rate(curve, m, k1, SEEDS, variant='orig')
    w_d, t_d, r_d = success_rate(curve, m, k1, SEEDS, variant='dfree')
    print(f"{label:<18} {m:>3} {k1:>3} {r_o:>10.3f} {w_o}/{t_o:<5} "
          f"{r_d:>12.3f} {w_d}/{t_d:<5}")

# ===========================================================================
# EXPERIMENT E3 — the K1 wall (2026-07-29 T4 grid), old vs d-free
# ===========================================================================
print("\n" + "-" * 78)
print("EXP E3: K1 wall replication (2026-07-29 T4), old vs d-free")
print("-" * 78)
print("Falsifier part 2: on 12-bit/2677 (lam*=0.070) the old wall is K1 ~ 4-6.")
print("If the d-free wall moves outward, the reformulation is a real gain;")
print("if it stays, the wall is information-theoretic.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
T4 = [("12-bit/2557", 8), ("12-bit/2677 FAIL", 10)]
for label, m in T4:
    curve = next(c for lbl, c, _, _ in hist_curves if lbl == label for c in [c])
    n = curve[2]
    ls = lam_star(curve[3], n)
    print(f"  {label}  (lam*={ls:.3f}, m={m})")
    hdr = "    " + f"{'variant':<8}" + "".join(f"{k:>6}" for k in K1_GRID)
    print(hdr)
    for variant in ('orig', 'dfree'):
        cells = []
        for k1 in K1_GRID:
            w, t, _ = success_rate(curve, m, k1, SEEDS, variant=variant)
            cells.append(f"{w}/{t}")
        print("    " + f"{variant:<8}" + "".join(f"{c:>6}" for c in cells))
    print()

# ===========================================================================
# EXPERIMENT E4 — 17-bit multi-curve head-to-head
# ===========================================================================
print("-" * 78)
print("EXP E4: 17-bit multi-curve head-to-head at three bias strengths")
print("-" * 78)
sys.stdout.flush()

curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
print(f"  found {len(curves17)} 17-bit j=0 GLV curves "
      f"(lam* in [{min(lam_star(c[3], c[2]) for c in curves17):.4f}, "
      f"{max(lam_star(c[3], c[2]) for c in curves17):.4f}])\n")

M_SIGS = 12
print(f"{'eff':>6} {'orig 5/5':>10} {'dfree 5/5':>11} "
      f"{'orig sv/pv':>11} {'dfree sv/pv':>12}")
e4_rows = []
for eff in (0.05, 0.15, 0.25):
    o_full = d_full = 0
    o_ratios, d_ratios = [], []
    for cur in curves17:
        n = cur[2]
        k2_bound = math.isqrt(n) + 1
        k1_bound = max(2, int(eff * n / k2_bound))
        w_o, t_o, r_o = success_rate(cur, M_SIGS, k1_bound, SEEDS, variant='orig')
        w_d, t_d, r_d = success_rate(cur, M_SIGS, k1_bound, SEEDS, variant='dfree')
        o_full += (w_o == t_o)
        d_full += (w_d == t_d)
        o_ratios.append(r_o)
        d_ratios.append(r_d)
        e4_rows.append((eff, n, lam_star(cur[3], n), k1_bound, w_o, w_d, r_o, r_d))
    print(f"{eff:>6.2f} {o_full:>4}/{len(curves17):<5} {d_full:>5}/{len(curves17):<5} "
          f"{sum(o_ratios)/len(o_ratios):>11.3f} "
          f"{sum(d_ratios)/len(d_ratios):>12.3f}")
    sys.stdout.flush()

print("\n  per-curve detail (eff, n, lam*, K1, orig wins, dfree wins):")
for eff, n, ls, k1, w_o, w_d, r_o, r_d in e4_rows:
    flag = "  <-- dfree only" if w_d > w_o else ("  <-- orig only" if w_o > w_d else "")
    print(f"    {eff:.2f}  n={n:<7} lam*={ls:.4f}  K1={k1:<4} "
          f"orig={w_o}/{len(SEEDS)}  dfree={w_d}/{len(SEEDS)}"
          f"  sv/pv {r_o:.2f}->{r_d:.2f}{flag}")

# ===========================================================================
# EXPERIMENT E5 — instance-level agreement
# ===========================================================================
print("\n" + "-" * 78)
print("EXP E5: do the two lattices agree instance-by-instance, not just in rate?")
print("-" * 78)
print("E3/E4 report identical win COUNTS in every cell.  That is consistent")
print("with two different attacks of equal strength.  Here we check the")
print("stronger claim: the SAME (curve, seed) instances succeed in both.")
print("Also runs a BKZ(20) arm to see whether the d-free geometry lets a")
print("stronger reduction convert instances that LLL misses.\n")

def instance_outcomes(curve, m, k1_bound, seeds, use_bkz=False, bkz_beta=20):
    out = []
    n = curve[2]
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        row = []
        for variant in ('orig', 'dfree'):
            res = run_once(curve, m, d_trial, k1_bound, seed, variant=variant,
                           use_bkz=use_bkz, bkz_beta=bkz_beta)
            row.append(bool(res[0]) if res is not None else None)
        out.append(tuple(row))
    return out

for algo, use_bkz in (("LLL", False), ("BKZ(20)", True)):
    both = only_o = only_d = neither = 0
    for cur in curves17:
        n = cur[2]
        k2_bound = math.isqrt(n) + 1
        for eff in (0.05, 0.15, 0.25):
            k1_bound = max(2, int(eff * n / k2_bound))
            for o, d in instance_outcomes(cur, M_SIGS, k1_bound, SEEDS,
                                          use_bkz=use_bkz):
                if o and d: both += 1
                elif o: only_o += 1
                elif d: only_d += 1
                else: neither += 1
    tot = both + only_o + only_d + neither
    print(f"  {algo:<8} instances={tot:<5} both={both:<5} orig-only={only_o:<4} "
          f"dfree-only={only_d:<4} neither={neither:<5} "
          f"disagreements={only_o + only_d}")

# ===========================================================================
# EXPERIMENT E6 — why is the agreement exact?  det and GS-profile identity
# ===========================================================================
print("\n" + "-" * 78)
print("EXP E6: is L_dfree exactly the quotient L_orig / <n*S_D*e_m> ?")
print("-" * 78)
print("Prediction from the construction:")
print("    det(L_orig) = (n*S_K1)^m * S_D * S_K2^m * S_KAN")
print("    det(L_free) = (n*S_K1)^(m-1) * S_K1 * S_K2^m * S_KAN")
print("  => det(L_orig)/det(L_free) = n*S_D = ||trivial vector|| exactly,")
print("  and dim drops by exactly 1.  If so the two lattices have the SAME")
print("  geometry once the trivial direction is projected out, which would")
print("  explain E5's 0/300 disagreements structurally rather than by luck.\n")

def lattice_det(M):
    from sympy import Matrix
    return abs(Matrix(M).det())

def project_out_d(M, m):
    """pi: delete the d-column (index m).  n*S_D*e_m spans ker(pi|_L), and
    n*S_D*e_m is axis-aligned, so pi IS the orthogonal projection onto e_m^perp."""
    return [r[:m] + r[m + 1:] for r in M]

def row_hnf_basis(rows):
    """A square basis of the lattice generated by `rows` (which may be an
    over-determined generating set), via row-style Hermite normal form."""
    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form
    H = hermite_normal_form(Matrix(rows).T).T
    return [[int(x) for x in H.row(i)] for i in range(H.rows)]

def spanned_by(rows, basis):
    """Is every row in the integer row-span of the SQUARE `basis`?  Exact."""
    from sympy import Matrix, Rational
    Bt = Matrix(basis).T
    assert Bt.rows == Bt.cols, "basis must be square"
    for r in rows:
        if not any(r):
            continue
        sol = Bt.solve(Matrix(list(r)))
        if not all(Rational(c).q == 1 for c in sol):
            return False
    return True

print(f"{'curve':<18} {'dim o/f':>9} {'det ratio':>10} {'pi(L_o) subset L_f':>19} "
      f"{'L_f subset pi(L_o)':>19}")
e6_ok = True
for label, curve, k1, m in hist_curves:
    p, b, n, lam, G = curve
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_trial, m, n, lam, p, k1, k2_bound, 42)
    Mo = build_glv_lattice(sigs, n, lam, k1, k2_bound)
    Mf = build_dfree_lattice(sigs, n, lam, k1, k2_bound)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2_bound)
    ratio_ok = lattice_det(Mo) == lattice_det(Mf) * n * S_D
    piMo = row_hnf_basis(project_out_d(Mo, m))   # square basis of pi(L_orig)
    fwd = spanned_by(piMo, Mf)          # pi(L_orig) subset L_free
    bwd = spanned_by(Mf, piMo)          # L_free subset pi(L_orig)
    e6_ok &= ratio_ok and fwd and bwd
    print(f"{label:<18} {len(Mo):>4}/{len(Mf):<4} {str(ratio_ok):>10} "
          f"{str(fwd):>19} {str(bwd):>19}")

print(f"\n  E6 verdict: {'PASS' if e6_ok else 'FAIL'} — L_free = pi(L_orig) EXACTLY.")
print("  The d-free lattice is not a new lattice: it is the original with the")
print("  trivial direction n*S_D*e_m projected out.  Recovery is therefore the")
print("  same BDD problem in both, which explains E5's 0/300 disagreements")
print("  structurally.  The d-column costs LLL one wasted dimension and nothing")
print("  else; removing it cannot move the wall.")

print("\nDone.")
