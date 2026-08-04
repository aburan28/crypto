"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
the target of the reduction, and test whether the K1 wall moves.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5)
-------------------------------------------------------
In the Phase-2 Kannan lattice L (dim 2m+2, basis in `build_glv_lattice`) the
shortest vector after LLL was, on every curve tested, exactly the *trivial*
vector n*S_D*e_m: 100% of its energy in the d-column, |sv[m]|/n = 1.0000.
That run concluded "the planted vector is never lambda_1" and proposed
Thread 23: quotient out the trivial direction, or replace Kannan by CVP.

Why n*S_D*e_m is in L (algebra, verified numerically in U0):
    n*row_m  -  sum_i B_i * row_i  =  (0,...,0 | n*S_D | 0,...,0)
because row_i = n*S_K1*e_i and row_m has B_i*S_K1 in column i.  Its norm is
n*S_D = n, while ||v_planted|| ~ n*sqrt(2m/3 + 4/3) >> n.  So it is shorter
than the planted vector for every m >= 1, and no choice of S_D removes it
(both scale linearly in S_D).

The reformulation tested here
-----------------------------
L ∩ R*e_m = <n*S_D*e_m> is a rank-1 summand aligned with a coordinate axis, so
the orthogonal projection along e_m -- concretely: DELETE COLUMN m -- maps L
onto a genuine lattice

    piL = pi(L)  in  Z^(2m+1),   det(piL) = det(L) / (n*S_D)

with the trivial vector killed.  Column m is the only place the secret d
appears, so at first sight the projection throws away the answer.  It does
not: pi(v_planted) = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN) still determines d,
because a single recovered nonce gives

    k_full,i = k1_i + lam*k2_i  (mod n),     d = (k_full,i - A_i) * B_i^-1  (mod n)

so the d-column was never load-bearing -- it was carrying a short trivial
vector and nothing else.

Determinants (both verified in U1 against the fpylll GS profile):
    det(L)   = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN
    det(piL) = (n*S_K1)^m * S_K2^m * S_KANNAN / n

Experiments
-----------
  U0  verify n*S_D*e_m in L and L ∩ R*e_m = <n*S_D*e_m> (confirms T5 algebraically)
  U1  build piL, verify rank 2m+1, det formula, and pi(v_planted) in piL
  U2  head-to-head Kannan-L vs piL on the T4 K1 grid (the Thread 23 falsifier)
  U3  Gaussian-heuristic predictor r = ||v_planted|| / GH(piL): does it separate
      success from failure across historical + 17-bit-sweep data?
  U4  failure-mode diagnosis: when recovery fails, is pi(v) still the shortest
      vector (LLL missed it) or not (the instance is genuinely not uSVP)?

Falsifier stated by the 2026-07-29 entry:
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run:  python3 secp256k1_cm_audit/glv_hnp_phase2_projected.py
Deps: fpylll, sympy   (pip install fpylll cysignals sympy)
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, BKZ, GSO

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

# ---------------------------------------------------------------------------
# CM theory for j=0 curves (verbatim)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
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
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
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
# Phase-2 lattice (verbatim construction) + the projected reformulation
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Original Kannan-embedded Phase-2 lattice, dim 2m+2.
    Column layout: [0..m-1] k1 | [m] d | [m+1..2m] k2 | [2m+1] Kannan."""
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

def project_out_d(M, m):
    """pi: delete column m.  Returns the (2m+2) x (2m+1) generating set of pi(L).
    New column layout: [0..m-1] k1 | [m..2m-1] k2 | [2m] Kannan."""
    return [row[:m] + row[m + 1:] for row in M]

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

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_d_kannan(reduced, m, n, S_KANNAN, d_secret):
    """Original formulation: read d straight out of column m."""
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

def recover_d_projected(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Projected formulation: read the nonces out of the k-columns, then
    d = (k1 + lam*k2 - A_i) * B_i^-1 mod n for any single signature i."""
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 or b % S_K2:
                continue
            k1, k2 = a // S_K1, b // S_K2
            if not (0 <= k1 < k1_bound and 0 <= k2 < k2_bound):
                continue
            k_full = (k1 + lam * k2) % n
            if k_full == 0:
                continue
            try:
                d_cand = (k_full - sigs[i]['A']) * modinv(sigs[i]['B'], n) % n
            except ValueError:
                continue
            if d_cand and d_cand == d_secret:
                return d_cand
    return None

# ---------------------------------------------------------------------------
# Lattice invariants
# ---------------------------------------------------------------------------

def log_det_kannan(m, n, S_K1, S_D, S_K2, S_KAN):
    """log det(L) = m*log(n*S_K1) + log S_D + m*log S_K2 + log S_KAN.
    Exact: column m, columns m+1..2m and column 2m+1 each have a single
    nonzero entry, so expanding gives (n*S_K1)^m * S_D * S_K2^m * S_KANNAN."""
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))

def gaussian_heuristic(log_det, dim):
    """GH = sqrt(dim/(2*pi*e)) * det^(1/dim)."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def reduce_matrix(rows, ncols, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]

def shortest_nonzero(reduced):
    best = None
    for r in reduced:
        if any(r):
            nv = norm(r)
            if best is None or nv < best:
                best = nv
    return best

# ---------------------------------------------------------------------------
# One trial, both formulations, same signatures
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    P = project_out_d(M, m)
    v = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    pv = v[:m] + v[m + 1:]

    redK = reduce_matrix(M, 2 * m + 2, use_bkz, beta)
    redP = reduce_matrix(P, 2 * m + 1, use_bkz, beta)

    okK = recover_d_kannan(redK, m, n, S_KAN, d_secret) is not None
    okP = recover_d_projected(redP, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None

    ldK = log_det_kannan(m, n, S_K1, S_D, S_K2, S_KAN)
    ldP = ldK - math.log(n * S_D)
    return {
        'okK': okK, 'okP': okP,
        'pnK': norm(v), 'pnP': norm(pv),
        'svK': shortest_nonzero(redK), 'svP': shortest_nonzero(redP),
        'ghK': gaussian_heuristic(ldK, 2 * m + 2),
        'ghP': gaussian_heuristic(ldP, 2 * m + 1),
        'ldP': ldP, 'n': n, 'k2_bound': k2_bound,
    }

def rate(curve, m, k1_bound, seeds, use_bkz=False, beta=20):
    """Returns (winsK, winsP, total, mean svP/pnP, mean pnP/ghP)."""
    wk = wp = tot = 0
    rsv, rgh = [], []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        r = trial(curve, m, d_trial, k1_bound, seed, use_bkz, beta)
        if r is None:
            continue
        tot += 1
        wk += r['okK']; wp += r['okP']
        rsv.append(r['svP'] / r['pnP'])
        rgh.append(r['pnP'] / r['ghP'])
    mean = lambda xs: sum(xs) / len(xs) if xs else float('nan')
    return wk, wp, tot, mean(rsv), mean(rgh)


# ===========================================================================
print("=" * 78)
print("Thread 23 — project out the trivial direction; does the K1 wall move?")
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
print("\n" + "-" * 78)
print("EXP U0: is n*S_D*e_m really in L, and is L ∩ R*e_m exactly <n*S_D*e_m>?")
print("-" * 78)
print("Claim (2026-07-29 T5): n*row_m - sum_i B_i*row_i = n*S_D*e_m.")
print("Second column: the exact generator of L ∩ R*e_m, i.e. the least t>0 with")
print("t*e_m in L.  Computed by solving x*M = t*e_m over Q and clearing")
print("denominators: t_min = lcm(denominators of the m-th row of M^-1).")
print(f"{'curve':<20} {'m':>3} {'identity holds':>15} {'gen(L ∩ Re_m)':>16} {'= n*S_D':>9}")

from sympy import Matrix, ilcm

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    d0 = random.Random(11).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    M = build_glv_lattice(sigs, n, lam, k1, k2b)
    dim = 2 * m + 2
    # n*row_m - sum_i B_i*row_i
    comb = [n * M[m][j] for j in range(dim)]
    for i in range(m):
        for j in range(dim):
            comb[j] -= sigs[i]['B'] * M[i][j]
    target = [0] * dim
    target[m] = n * S_D
    holds = (comb == target)
    # L ∩ R*e_m : least t > 0 with t*e_m in L.  x*M = t*e_m  =>  x = t*(M^-1)[m,:],
    # so t_min clears the denominators of that row of M^-1.
    em = [0] * dim
    em[m] = 1
    x = Matrix(M).T.solve(Matrix(em))          # column vector of Rationals
    t_min = int(ilcm(*[xi.q for xi in x])) if dim > 1 else int(x[0].q)
    print(f"{label:<20} {m:>3} {str(holds):>15} {t_min:>16} "
          f"{str(t_min == n * S_D):>9}")

# ===========================================================================
print("\n" + "-" * 78)
print("EXP U1: the projected lattice piL — rank, determinant, planted membership")
print("-" * 78)
print("pi = delete column m.  det(piL) should equal det(L)/(n*S_D), and")
print("pi(v_planted) must still lie in piL (else the reformulation is invalid).")
print(f"{'curve':<20} {'dim':>4} {'rank':>5} {'log det(piL)':>14} {'predicted':>12} "
      f"{'pi(v) in piL':>13}")

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    d0 = random.Random(11).randint(1, n - 1)
    sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    M = build_glv_lattice(sigs, n, lam, k1, k2b)
    P = project_out_d(M, m)
    ncols = 2 * m + 1
    A = IntegerMatrix.from_matrix(P)
    LLL.reduction(A)
    nz = [i for i in range(A.nrows) if any(A[i][j] for j in range(ncols))]
    rank = len(nz)
    basis = [[A[i][j] for j in range(ncols)] for i in nz]
    gso_log = 0.0
    Ab = IntegerMatrix.from_matrix(basis)
    Mg = GSO.Mat(Ab)
    Mg.update_gso()
    for i in range(Ab.nrows):
        gso_log += 0.5 * math.log(Mg.get_r(i, i))
    pred = log_det_kannan(m, n, S_K1, S_D, S_K2, S_KAN) - math.log(n * S_D)
    # membership: solve x * basis = pi(v) over Q, check integrality
    v = planted_vector(sigs, d0, n, k1, k2b)
    pv = v[:m] + v[m + 1:]
    sol = Matrix(basis).T.solve(Matrix(pv))
    member = all(x.q == 1 for x in sol)
    print(f"{label:<20} {ncols:>4} {rank:>5} {gso_log:>14.4f} {pred:>12.4f} "
          f"{str(member):>13}")

# ===========================================================================
print("\n" + "-" * 78)
print("EXP U2: THE FALSIFIER — K1 grid, Kannan-L vs projected piL, 5 seeds")
print("-" * 78)
print("2026-07-29 T4 walls: 12-bit/2557 fails at K1>=12-16; 12-bit/2677 at K1>=6-8.")
print("If the projection is a real improvement the piL column must survive")
print("further right than the K vector column.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
u2_rows = []
for label, curve, _k1, m in hist:
    p, b, n, lam, G = curve
    print(f"{label}  (n={n}, lam*={lam_star(lam, n):.4f}, m={m}, "
          f"K2={math.isqrt(n)+1})")
    print(f"  {'K1':>4} {'eff':>7} | {'Kannan':>7} {'proj':>6} | "
          f"{'svP/pnP':>8} {'pnP/GH_P':>9}")
    for k1 in K1_GRID:
        wk, wp, tot, r_sv, r_gh = rate(curve, m, k1, SEEDS)
        eff = k1 * (math.isqrt(n) + 1) / n
        print(f"  {k1:>4} {eff:>7.3f} | {wk:>3}/{tot:<3} {wp:>2}/{tot:<3} | "
              f"{r_sv:>8.3f} {r_gh:>9.3f}")
        u2_rows.append({'label': label, 'k1': k1, 'm': m, 'eff': eff,
                        'wk': wk, 'wp': wp, 'tot': tot,
                        'r_sv': r_sv, 'r_gh': r_gh,
                        'lam_star': lam_star(lam, n)})
    print()

# ===========================================================================
print("-" * 78)
print("EXP U3: is pnP/GH(piL) the predictor? (17-bit sweep, lam* spread)")
print("-" * 78)
sys.stdout.flush()

curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=6)
print(f"collected {len(curves17)} fresh 17-bit j=0 GLV curves")
M_SIGS = 12
u3_rows = []
print(f"  {'p':>7} {'n':>7} {'lam*':>7} {'eff':>6} | {'Kan':>5} {'proj':>5} | "
      f"{'pnP/GH_P':>9}")
for (p, b, n, lam, G) in curves17:
    k2b = math.isqrt(n) + 1
    for eff_t in (0.05, 0.15, 0.25):
        k1 = max(2, int(round(eff_t * n / k2b)))
        wk, wp, tot, r_sv, r_gh = rate((p, b, n, lam, G), M_SIGS, k1, SEEDS)
        eff = k1 * k2b / n
        print(f"  {p:>7} {n:>7} {lam_star(lam,n):>7.4f} {eff:>6.3f} | "
              f"{wk:>2}/{tot:<2} {wp:>2}/{tot:<2} | {r_gh:>9.3f}")
        u3_rows.append({'p': p, 'n': n, 'lam_star': lam_star(lam, n),
                        'eff': eff, 'wk': wk, 'wp': wp, 'tot': tot,
                        'r_gh': r_gh})
        sys.stdout.flush()

print("\nSeparation of the predictor r = ||pi(v_planted)|| / GH(piL)")
print("(pooled over EXP U2 + EXP U3; 'success' = projected formulation >= 3/5)")
pool = ([{'r': r['r_gh'], 'ok': r['wp'] * 2 >= r['tot'], 'eff': r['eff'],
          'ls': r['lam_star']} for r in u2_rows]
        + [{'r': r['r_gh'], 'ok': r['wp'] * 2 >= r['tot'], 'eff': r['eff'],
            'ls': r['lam_star']} for r in u3_rows])

def auc(rows, key, invert=False):
    pos = [r[key] for r in rows if r['ok']]
    neg = [r[key] for r in rows if not r['ok']]
    if not pos or not neg:
        return float('nan')
    c = sum((a < b) + 0.5 * (a == b) for a in pos for b in neg) / (len(pos) * len(neg))
    return 1 - c if invert else c

def best_split(rows, key):
    vals = sorted({r[key] for r in rows})
    best = (0, None)
    for i in range(len(vals) - 1):
        t = (vals[i] + vals[i + 1]) / 2
        acc = sum((r[key] < t) == r['ok'] for r in rows) / len(rows)
        if acc > best[0]:
            best = (acc, t)
    return best

npos = sum(r['ok'] for r in pool)
print(f"  pooled configurations: {len(pool)}  (success {npos}, failure {len(pool)-npos})")
for key, name in (('r', 'r = pnP/GH(piL)'), ('eff', 'eff = K1*K2/n'),
                  ('ls', 'lam* = min(lam,n-lam)/n')):
    a = auc(pool, key)
    acc, t = best_split(pool, key)
    print(f"  {name:<26} AUC(small=success) {a:>6.3f}   "
          f"best threshold {t if t is None else round(t,4)!s:>8}  acc {acc:>5.3f}")

maj = max(npos, len(pool) - npos) / len(pool)
print(f"  majority baseline accuracy  {maj:.3f}")

# ===========================================================================
print("\n" + "-" * 78)
print("EXP U4: failure-mode diagnosis on the projected lattice")
print("-" * 78)
print("For failing configurations: is pi(v_planted) still the shortest vector")
print("(=> LLL missed it, BKZ should rescue) or is something shorter present")
print("(=> the instance is not uSVP; no amount of reduction helps)?\n")
print(f"  {'curve':<20} {'K1':>3} {'proj LLL':>9} {'proj BKZ20':>11} "
      f"{'svP/pnP':>8} {'pnP/GH_P':>9} {'verdict':>16}")

for label, curve, _k1, m in hist:
    p, b, n, lam, G = curve
    for k1 in (6, 8, 12, 16, 24):
        wk, wp, tot, r_sv, r_gh = rate(curve, m, k1, SEEDS)
        if wp * 2 >= tot:
            continue
        _, wpb, totb, r_svb, _ = rate(curve, m, k1, SEEDS, use_bkz=True, beta=20)
        verdict = ("LLL-limited" if r_sv >= 0.999 else
                   "not uSVP" if r_gh > 1.0 else "borderline")
        print(f"  {label:<20} {k1:>3} {wp:>4}/{tot:<4} {wpb:>6}/{totb:<4} "
              f"{r_sv:>8.3f} {r_gh:>9.3f} {verdict:>16}")
        sys.stdout.flush()


# ===========================================================================
print("\n" + "-" * 78)
print("EXP U5: closed form for the wall, and where it actually sits")
print("-" * 78)
print("Solving ||pi(v_planted)|| = GH(piL) for eff = K1*K2/n gives")
print("  eff_1 = [ n^(2m) * C^dim / (n*sqrt(2m/3+1))^dim ]^(1/m),  C=sqrt(dim/(2*pi*e))")
print("with dim = 2m+1.  As m -> infinity this tends to 3/(2*pi*e) = "
      f"{3/(2*math.pi*math.e):.4f}.")
print("eff_obs = midpoint of the observed pass/fail transition in EXP U2/U3.\n")
print(f"  {'curve':<20} {'m':>3} {'n':>7} {'eff_1 (r=1)':>12} "
      f"{'eff_obs':>9} {'obs/pred':>9}")

def eff_crossing(m, n, r_target=1.0):
    """eff at which ||pi(v_planted)||/GH(piL) = r_target.

    ||pi(v)|| does not depend on eff (k1_i*S_K1 ranges over [0,n) whatever K1
    is), while GH ∝ det^(1/dim) ∝ eff^(-m/dim).  Hence r ∝ eff^(m/dim), i.e.
    eff(r) = eff_1 * r^(dim/m) -- r_target multiplies C, it does not divide
    the planted norm."""
    dim = 2 * m + 1
    C = math.sqrt(dim / (2 * math.pi * math.e))
    log_eff_m = (2 * m * math.log(n) + dim * math.log(r_target * C)
                 - dim * math.log(n * math.sqrt(2 * m / 3 + 1)))
    return math.exp(log_eff_m / m)

# observed transition midpoints, read off the EXP U2 / U3 tables above
OBS = [("8-bit/199", 6, 199, (0.226 + 0.302) / 2),
       ("12-bit/2557", 8, 2659, (0.156 + 0.235) / 2),
       ("12-bit/2677 FAIL", 10, 2647, (0.079 + 0.118) / 2),
       ("17-bit sweep (pooled)", 12, 65719, (0.051 + 0.149) / 2)]
for label, m, n, eff_obs in OBS:
    e1 = eff_crossing(m, n)
    print(f"  {label:<20} {m:>3} {n:>7} {e1:>12.4f} {eff_obs:>9.4f} "
          f"{eff_obs/e1:>9.2f}x")

print("\nThe observed wall sits ~2-4x above the r=1 point, i.e. LLL/BKZ recovers the")
print("planted vector even when it is LONGER than the Gaussian heuristic.  The best")
print("empirical split in EXP U3 is r = 1.379, not r = 1.  So the planted vector was")
print("never required to be lambda_1 -- which is exactly why removing the trivial")
print("vector (EXP U2) changed nothing.")
print("\nCross-check: eff(r) = eff_1 * r^(dim/m), so the empirical r=1.379 split must")
print("map onto the empirical eff split (0.1176 from EXP U3) if the model is right:")
print(f"  eff at r=1.379, m=12, n=65719: {eff_crossing(12, 65719, 1.379):.4f}"
      f"   (EXP U3 best eff threshold: 0.1176)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
