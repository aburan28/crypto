"""
GLV-HNP Thread 23 — reformulate the Phase-2 lattice so the planted vector is lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, finding T5): in the Phase-2
Kannan lattice the shortest vector is ALWAYS the trivial vector n*S_D*e_m
(100% of its energy in the d-column, |sv[m]|/n = 1.0000 exactly), never the
planted vector.  With S_D = 1 its norm is n, while
    ||v_planted||^2 ~ n^2 * (2m/3 + 4/3),
so the planted vector loses to it for every m >= 1.  The trivial vector carries
no information (d is only defined mod n) and no choice of S_D removes it (both
scale linearly in S_D).

This script tests the three reformulations proposed as Thread 23:

  L_A  "orig"  the 2026-07-29 lattice, verbatim         dim 2m+2 (baseline)
  L_B  "proj"  L_A with the d-column deleted            dim 2m+1
               (= quotient by the trivial direction; d is re-derived from
                k1_i = A_i + B_i*d - lam*k2_i mod n afterwards)
  L_C  "cvp"   no Kannan row, no d-column; Babai        dim 2m
               nearest-plane against the centred target

Falsifier stated on 2026-07-29:
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, the wall is
   information-theoretic and Phase 2 is at its ceiling."

Everything downstream of the lattice (curve set, signature generation, scales,
seeds) is copied verbatim from glv_hnp_phase2_lambda_threshold.py so the
comparison against the 2026-07-29 numbers is exact.

Run: python3 glv_hnp_phase23_projected.py
"""

import math
import random

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
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p)."""
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
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n.  Requires n = 1 mod 3."""
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
    """Symmetric size of lam: min(lam, n-lam)/n in [0, 0.5]."""
    return min(lam % n, n - (lam % n)) / n

def build_curve(p, n, seed=12345):
    """Find the sextic twist b with #E = n, plus a generator."""
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

def search_curves(lo, hi, per_bin=1, nbins=10, max_primes=100000):
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
# Signature generation + scales (verbatim from 2026-07-29)
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

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ===========================================================================
# L_A — original Kannan lattice, dim 2m+2   (verbatim)
# ===========================================================================

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound):
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

def planted_A(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_A(rows, m, n, S_KANNAN, d_secret):
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

# ===========================================================================
# L_B — d-column deleted, dim 2m+1
# ===========================================================================
# Column layout after deletion:  0..m-1 = k1 cols, m..2m-1 = k2 cols, 2m = Kannan.
# The generating set stays the same 2m+2 rows, so L_B has rank 2m+1 with one
# dependency (n*row_B - sum_i B_i*row_i = 0); fpylll returns a zero row for it.

def build_lattice_B(sigs, n, lam, k1_bound, k2_bound):
    A = build_lattice_A(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    return [row[:m] + row[m + 1:] for row in A]

def planted_B(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    v = planted_A(sigs, d_secret, n, k1_bound, k2_bound)
    return v[:m] + v[m + 1:]

def recover_B(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """d is re-derived from k1_i = A_i + B_i*d - lam*k2_i (mod n)."""
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        for i in range(m):
            xi = sign * row[i]
            yi = sign * row[m + i]
            if xi % S_K1 or yi % S_K2:
                continue
            k1i, k2i = xi // S_K1, yi // S_K2
            B = sigs[i]['B'] % n
            if math.gcd(B, n) != 1:
                continue
            d_cand = (k1i - sigs[i]['A'] + lam * k2i) * modinv(B, n) % n
            if d_cand == 0: continue
            if d_cand == d_secret:
                return d_cand
    return None

# ===========================================================================
# L_C — no Kannan row, no d-column; Babai nearest-plane, dim 2m
# ===========================================================================
# Rows:  n*S_K1*e_i (i<m);  row_B: col i = B_i*S_K1;
#        row_{2,i}: col i = -lam*S_K1, col m+i = S_K2.
# A lattice point is  x_i = (B_i*d - lam*y_i - c_i*n)*S_K1,  x_{m+i} = y_i*S_K2.
# Target (centred):  t_i = (-A_i + K1/2)*S_K1,  t_{m+i} = (K2/2)*S_K2,
# so that x - t_uncentred has the planted offset (k1_i*S_K1, k2_i*S_K2).

def build_lattice_C(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(m + 1 + m)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    return M

def target_C(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = (-sigs[i]['A'] + k1_bound // 2) * S_K1
        t[m + i] = (k2_bound // 2) * S_K2
    return t

def recover_C(x, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """x = the Babai closest lattice point; read k1_i, k2_i and solve for d."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    votes = {}
    for i in range(m):
        if x[i] % S_K1 or x[m + i] % S_K2:
            continue
        k2i = x[m + i] // S_K2
        B = sigs[i]['B'] % n
        if math.gcd(B, n) != 1:
            continue
        d_cand = (x[i] // S_K1 + lam * k2i) * modinv(B, n) % n
        votes[d_cand] = votes.get(d_cand, 0) + 1
    if not votes:
        return None
    best = max(votes, key=votes.get)
    return best if best == d_secret else None

# ---------------------------------------------------------------------------
# Reduction drivers
# ---------------------------------------------------------------------------

def reduce_rows(M, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def shortest_nonzero(rows):
    best, bn = None, None
    for r in rows:
        nn = sum(x * x for x in r)
        if nn == 0:
            continue
        if bn is None or nn < bn:
            best, bn = r, nn
    return best, (math.sqrt(bn) if bn else 0.0)

def babai_closest(M, t, use_bkz=False, beta=20):
    """Babai nearest-plane against target t.  Returns the closest lattice point.

    L_C is given by 2m+1 generators of rank 2m (row_B lies in the Q-span of the
    modular rows), so LLL emits a leading zero row; GSO.Mat needs it dropped or
    the nearest-plane recursion divides by a zero Gram-Schmidt norm.
    """
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    basis = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)
             if any(A[i][j] for j in range(A.ncols))]
    B = IntegerMatrix.from_matrix(basis)
    G = GSO.Mat(B)
    G.update_gso()
    coeffs = G.babai([int(v) for v in t])
    dim = B.ncols
    x = [0] * dim
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        for j in range(dim):
            x[j] += c * B[i][j]
    return x

# ---------------------------------------------------------------------------
# One trial: run all three formulations on the same signature set
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, variants=('A', 'B', 'C'),
          use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {}
    if 'A' in variants:
        MA = build_lattice_A(sigs, n, lam, k1_bound, k2_bound)
        rA = reduce_rows(MA, use_bkz, beta)
        svA, svnA = shortest_nonzero(rA)
        pvA = norm(planted_A(sigs, d_secret, n, k1_bound, k2_bound))
        S_KAN = scales(n, k1_bound, k2_bound)[3]
        out['A'] = {
            'ok': recover_A(rA, m, n, S_KAN, d_secret) is not None,
            'svpv': svnA / pvA,
            'sv': svA,
            'pv': pvA,
        }
    if 'B' in variants:
        MB = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
        rB = reduce_rows(MB, use_bkz, beta)
        svB, svnB = shortest_nonzero(rB)
        pvB = norm(planted_B(sigs, d_secret, n, k1_bound, k2_bound))
        out['B'] = {
            'ok': recover_B(rB, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None,
            'svpv': svnB / pvB,
            'sv': svB,
            'pv': pvB,
        }
    if 'C' in variants:
        MC = build_lattice_C(sigs, n, lam, k1_bound, k2_bound)
        t = target_C(sigs, n, k1_bound, k2_bound)
        x = babai_closest(MC, t)
        out['C'] = {'ok': recover_C(x, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None}
    return out

def rate(curve, m, k1_bound, seeds, variants=('A', 'B', 'C'), use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    cnt = {v: 0 for v in variants}
    svpv = {v: [] for v in variants}
    for s in seeds:
        d = random.Random(s * 7919 + 13).randrange(1, n)
        r = trial(curve, m, d, k1_bound, s, variants, use_bkz, beta)
        if r is None:
            continue
        for v in variants:
            if r[v]['ok']:
                cnt[v] += 1
            if 'svpv' in r[v]:
                svpv[v].append(r[v]['svpv'])
    return cnt, {v: (sum(x) / len(x) if x else float('nan')) for v, x in svpv.items()}


SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26 / 2026-07-29)
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

print("=" * 78)
print("Thread 23 — reformulating the Phase-2 lattice so the planted vector is lambda_1")
print("=" * 78)

# ===========================================================================
# E1 — does the trivial vector actually disappear?  (sv/pv, energy split)
# ===========================================================================
print()
print("-" * 78)
print("E1  structure: shortest vector after LLL, per formulation")
print("-" * 78)
print(f"{'curve':<18}{'lam*':>7}{'K1':>4}{'m':>4} | "
      f"{'A sv/pv':>9}{'A |sv_d|/n':>12} | {'B sv/pv':>9}{'B is pv':>9}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    d = random.Random(42 * 7919 + 13).randrange(1, n)
    r = trial(curve, m, d, k1, 42, variants=('A', 'B'))
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, math.isqrt(n) + 1)
    svA = r['A']['sv']
    dcol = abs(svA[m]) / n
    pvB = planted_B(gen_signatures(G, d, m, n, lam, p, k1, math.isqrt(n) + 1, 42),
                    d, n, k1, math.isqrt(n) + 1)
    svB = r['B']['sv']
    is_pv = (svB == pvB) or (svB == [-x for x in pvB])
    print(f"{label:<18}{lam_star(lam, n):>7.3f}{k1:>4}{m:>4} | "
          f"{r['A']['svpv']:>9.3f}{dcol:>12.4f} | {r['B']['svpv']:>9.3f}{str(is_pv):>9}")

# ===========================================================================
# E2 — the T4 K1-grid, replicated for all three formulations (THE FALSIFIER)
# ===========================================================================
print()
print("-" * 78)
print("E2  K1 wall: successes out of 5 seeds, per formulation  (2026-07-29 T4 grid)")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
for label, curve, _k1, m in hist[1:]:
    p, b, n, lam, G = curve
    print(f"\n{label}   lam*={lam_star(lam, n):.3f}  m={m}  n={n}")
    hdr = "  ".join(f"K1={k:<3}" for k in K1_GRID)
    print(f"  {'variant':<10}{hdr}")
    for v in ('A', 'B', 'C'):
        cells = []
        for k1 in K1_GRID:
            cnt, _ = rate(curve, m, k1, SEEDS, variants=(v,))
            cells.append(f"{cnt[v]}/5 ")
        print(f"  {v:<10}" + "  ".join(f"{c:<5}" for c in cells))

# ===========================================================================
# E3 — fresh 17-bit curves at the discriminating bias strength eff = 0.15
# ===========================================================================
print()
print("-" * 78)
print("E3  fresh 17-bit sweep at eff = K1*K2/n = 0.15  (2026-07-29: A scored 3/20)")
print("-" * 78)
curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
print(f"  found {len(curves17)} curves")
M_SIGS = 12
tot = {'A': 0, 'B': 0, 'C': 0}
full = {'A': 0, 'B': 0, 'C': 0}
print(f"  {'p':>7}{'n':>8}{'lam*':>7}{'K1':>5} | {'A':>5}{'B':>5}{'C':>5}")
for (p, b, n, lam, G) in curves17:
    k2 = math.isqrt(n) + 1
    k1 = max(2, round(0.15 * n / k2))
    cnt, _ = rate((p, b, n, lam, G), M_SIGS, k1, SEEDS)
    for v in ('A', 'B', 'C'):
        tot[v] += cnt[v]
        full[v] += 1 if cnt[v] == len(SEEDS) else 0
    print(f"  {p:>7}{n:>8}{lam_star(lam, n):>7.4f}{k1:>5} | "
          f"{cnt['A']:>5}{cnt['B']:>5}{cnt['C']:>5}")
den = len(curves17) * len(SEEDS)
print(f"\n  totals (seed-level): A {tot['A']}/{den}   B {tot['B']}/{den}   C {tot['C']}/{den}")
print(f"  totals (5/5 curves): A {full['A']}/{len(curves17)}   "
      f"B {full['B']}/{len(curves17)}   C {full['C']}/{len(curves17)}")

# ===========================================================================
# E4 — where is the NEW bias-strength ceiling?  sweep eff for A vs C
# ===========================================================================
print()
print("-" * 78)
print("E4  bias-strength ceiling: eff sweep on the same 20 curves (m=12, 5 seeds)")
print("-" * 78)
print(f"  {'eff':>6} | {'A seed-lvl':>11}{'A 5/5':>8} | {'C seed-lvl':>11}{'C 5/5':>8}")
for eff in (0.15, 0.25, 0.35, 0.50, 0.70):
    t = {'A': 0, 'C': 0}
    f = {'A': 0, 'C': 0}
    for (p, b, n, lam, G) in curves17:
        k2 = math.isqrt(n) + 1
        k1 = max(2, round(eff * n / k2))
        cnt, _ = rate((p, b, n, lam, G), M_SIGS, k1, SEEDS, variants=('A', 'C'))
        for v in ('A', 'C'):
            t[v] += cnt[v]
            f[v] += 1 if cnt[v] == len(SEEDS) else 0
    print(f"  {eff:>6.2f} | {t['A']:>7}/{den:<3}{f['A']:>5}/{len(curves17):<3} | "
          f"{t['C']:>7}/{den:<3}{f['C']:>5}/{len(curves17):<3}")

print()
print("=" * 78)
print("done")
print("=" * 78)
