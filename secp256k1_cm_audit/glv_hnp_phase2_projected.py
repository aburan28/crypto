"""
GLV-HNP Phase 2, Thread 23: make the planted vector lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, finding T5)
-----------------------------------------------------------
The (2m+2)-dimensional Phase-2 lattice L contains the vector n*S_D*e_m
(take n*row_m and subtract B_i*row_i for each i).  With S_D = 1 its norm
is exactly n, while

    ||v_planted||^2 ~ n^2 * (2m/3 + 4/3),

so for every m >= 1 the trivial vector is strictly shorter and the planted
vector is NEVER lambda_1.  The trivial vector carries no information: d is
only defined mod n, and both vectors scale linearly in S_D, so no choice of
S_D removes it.

Thread 23 asks: quotient that direction out and see whether recovery improves.

The fix implemented here
------------------------
L contains n*S_D*e_m, so L / <n*S_D*e_m> is a rank-(2m+1) lattice, and the
quotient map is simply "delete column m".  Concretely:

    L_proj = pi(L),  pi(x_0..x_{2m+1}) = (x_0..x_{m-1}, x_{m+1}..x_{2m+1})

    columns  0    .. m-1   k1 block, scale S_K1
             m    .. 2m-1  k2 block, scale S_K2
             2m            Kannan,   scale S_KANNAN

pi(v_planted) = (k1_i*S_K1, k2_i*S_K2, S_KANNAN), norm^2 ~ n^2*(2m/3 + 1).
No information is lost: from k1_0 and k2_0 the secret comes back as

    d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n)

because A_0 + B_0*d - lam*k2_0 = k1_0 (mod n) is the defining congruence.

Determinants:  det(L)      = n^(2m+1) / eff^m       (dim 2m+2)
               det(L_proj) = det(L)/(n*S_D)
                           = n^(2m)   / eff^m       (dim 2m+1)
with eff = K1*K2/n, so the Gaussian-heuristic length is essentially unchanged
-- the projection removes a masking vector, it does not change the density.

A third arm drops the Kannan row from L_proj as well and calls Babai
nearest-plane directly on the target (-A_i*S_K1, 0, ..., 0), i.e. solves the
CVP the Kannan embedding was standing in for.

Falsifier (stated 2026-07-29)
-----------------------------
  * if sv/pv rises above 1 (planted vector becomes lambda_1) AND the K1 wall
    moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
    reformulation is a real improvement;
  * if the wall stays at K1 ~ 4-6, the wall is information-theoretic and
    Phase 2 is at its ceiling.

Experiments
-----------
  T23-A  geometry: pv, sv, sv/pv and energy split for L and L_proj on the
         three historical anchor curves.
  T23-B  the falsifier: T4's K1 grid on the two 12-bit anchors, L vs L_proj
         vs Babai-CVP.
  T23-C  the T3 17-bit sweep (20 curves, m=12, 5 seeds) at eff in
         {0.05, 0.15, 0.25}, L vs L_proj.  Baseline to beat: 19/20, 3/20, 0/20.
  T23-D  BKZ-20 on the survivors, to separate "structure" from
         "strength of reduction".

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sys
from fractions import Fraction

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
    mm, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (Eisenstein norm form)."""
    if p % 3 != 1:
        return None
    lim = int(2 * math.isqrt(p) / math.sqrt(3)) + 2
    for a in range(-lim, lim + 1):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            continue
        r = math.isqrt(disc)
        if r * r != disc:
            continue
        for bb in ((a + r) // 2, (a - r) // 2):
            if a * a - a * bb + bb * bb == p:
                return (a, bb)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, or None."""
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

# ---------------------------------------------------------------------------
# Phase-2 lattice (baseline construction, verbatim) and its projection
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

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Baseline (2m+2)-dim Kannan lattice.  Columns:
       0..m-1 k1 | m d | m+1..2m k2 | 2m+1 Kannan."""
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

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """L / <n*S_D*e_m>, realised as "delete column m".  2m+2 rows spanning a
    rank-(2m+1) lattice in Z^(2m+1).  Columns:
       0..m-1 k1 | m..2m-1 k2 | 2m Kannan."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):                       # n*S_K1*e_i
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    row = [0] * dim                          # d-row, d-column deleted
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for i in range(m):                       # k2 rows
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    row = [0] * dim                          # Kannan row
    for i in range(m):
        row[i] = sigs[i]['A'] * S_K1
    row[2 * m] = S_KANNAN
    M.append(row)
    return M

def build_cvp_basis_and_target(sigs, n, lam, k1_bound, k2_bound):
    """L_proj with the Kannan row removed: a rank-2m lattice in Z^2m, plus the
    CVP target  t = (-A_i*S_K1, 0, ..., 0).  The closest vector v satisfies
    t - v = (k1_i*S_K1, -k2_i*S_K2) up to sign."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    row = [0] * dim
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    t = [0] * dim
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return M, t

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

def planted_vector_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_d_kannan(rows, m, n, S_KANNAN, d_secret):
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

def _d_from_k(k1_0, k2_0, sigs, n, lam):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    if B0 % n == 0:
        return None
    return (k1_0 + lam * k2_0 - A0) * modinv(B0 % n, n) % n

def recover_d_proj(rows, m, n, lam, sigs, S_K1, S_K2, S_KANNAN, d_secret):
    """Rows of the projected lattice.  d is reconstructed from (k1_0, k2_0)."""
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        a, b = sign * row[0], sign * row[m]
        if a % S_K1 or b % S_K2: continue
        d_cand = _d_from_k(a // S_K1, b // S_K2, sigs, n, lam)
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

def babai_nearest_plane(basis_rows, target):
    """Exact nearest-plane on an LLL-reduced basis (Fraction GS).  Small dims
    only; used for the CVP arm."""
    B = [[Fraction(x) for x in r] for r in basis_rows]
    k = len(B)
    gs, mu, nrm = [], [], []
    for i in range(k):
        v = list(B[i])
        mu_i = []
        for j in range(len(gs)):
            if nrm[j] == 0:
                mu_i.append(Fraction(0)); continue
            c = sum(B[i][t] * gs[j][t] for t in range(len(v))) / nrm[j]
            mu_i.append(c)
            v = [v[t] - c * gs[j][t] for t in range(len(v))]
        gs.append(v); mu.append(mu_i)
        nrm.append(sum(x * x for x in v))
    b = [Fraction(x) for x in target]
    coeffs = [0] * k
    for i in range(k - 1, -1, -1):
        if nrm[i] == 0:
            continue
        c = sum(b[t] * gs[i][t] for t in range(len(b))) / nrm[i]
        c_int = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = c_int
        b = [b[t] - c_int * Fraction(B[i][t]) for t in range(len(b))]
    v = [sum(coeffs[i] * basis_rows[i][t] for i in range(k))
         for t in range(len(target))]
    return v

# ---------------------------------------------------------------------------
# Single trials
# ---------------------------------------------------------------------------

def lll_rows(M, beta=None):
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def shortest_nonzero(rows):
    best, bv = None, None
    for r in rows:
        nn = norm(r)
        if nn > 0 and (best is None or nn < best):
            best, bv = nn, r
    return best, bv

def trial(curve, m, d_secret, k1_bound, seed, variant, beta=None):
    """variant in {'kannan', 'proj', 'cvp'}.  Returns dict with ok/pv/sv."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None

    if variant == 'kannan':
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        rows = lll_rows(M, beta)
        pv = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
        sv, svv = shortest_nonzero(rows)
        ok = recover_d_kannan(rows, m, n, S_KAN, d_secret) is not None
        return {'ok': ok, 'pv': pv, 'sv': sv, 'svv': svv, 'm': m}

    if variant == 'proj':
        M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
        rows = lll_rows(M, beta)
        pv = norm(planted_vector_proj(sigs, n, k1_bound, k2_bound))
        sv, svv = shortest_nonzero(rows)
        ok = recover_d_proj(rows, m, n, lam, sigs, S_K1, S_K2, S_KAN,
                            d_secret) is not None
        return {'ok': ok, 'pv': pv, 'sv': sv, 'svv': svv, 'm': m}

    if variant == 'cvp':
        M, t = build_cvp_basis_and_target(sigs, n, lam, k1_bound, k2_bound)
        rows = lll_rows(M, beta)
        rows = [r for r in rows if any(r)]
        v = babai_nearest_plane(rows, t)
        e = [t[i] - v[i] for i in range(len(t))]
        ok = False
        for sign in (1, -1):
            a, bb = sign * e[0], -sign * e[m]
            if a % S_K1 or bb % S_K2: continue
            d_cand = _d_from_k(a // S_K1, bb // S_K2, sigs, n, lam)
            if d_cand and d_cand == d_secret:
                ok = True
        pv = norm([s['k1'] * S_K1 for s in sigs] + [s['k2'] * S_K2 for s in sigs])
        return {'ok': ok, 'pv': pv, 'sv': norm(e), 'svv': e, 'm': m}

    raise ValueError(variant)

def rate(curve, m, k1_bound, seeds, variant, beta=None, d_secret=None):
    p, b, n, lam, G = curve
    d = d_secret if d_secret is not None else n // 3 + 7
    wins, ratios = 0, []
    for s in seeds:
        r = trial(curve, m, d % n, k1_bound, s, variant, beta)
        if r is None:
            continue
        wins += bool(r['ok'])
        if r['pv']:
            ratios.append(r['sv'] / r['pv'])
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr

# ---------------------------------------------------------------------------
# Curve construction / search
# ---------------------------------------------------------------------------

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

def search_curves(lo, hi, per_bin=2, nbins=10):
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
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

# ===========================================================================

if __name__ == "__main__":

    print("=" * 78)
    print("Thread 23 -- project out the trivial n*e_m direction of the Phase-2 lattice")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,            p,    b, n,    lam,   K1, m
        ("8-bit/199",       211,  2, 199,  106,   2,  6),
        ("12-bit/2557",     2557, 2, 2659, 1755,  8,  8),
        ("12-bit/2677",     2677, 2, 2647, 185,   8,  10),
    ]

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), k1, m))

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("T23-A: geometry of L (baseline) vs L_proj (d-column deleted)")
    print("-" * 78)
    print("sv = shortest nonzero row after LLL, pv = ||planted||.")
    print("Blocks of sv as a fraction of its energy: k1 | d | k2 | kan.\n")

    print(f"{'curve':<14}{'K1':>3}{'m':>3} {'lat':<7}"
          f"{'pv/n':>8}{'sv/n':>8}{'sv/pv':>8}  {'k1':>6}{'d':>6}{'k2':>6}{'kan':>6}  ok")
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        d_secret = n // 3 + 7
        for variant in ('kannan', 'proj'):
            r = trial(curve, m, d_secret, k1, 42, variant)
            sv_v = r['svv']
            if variant == 'kannan':
                e_k1 = sum(float(x) ** 2 for x in sv_v[:m])
                e_d = float(sv_v[m]) ** 2
                e_k2 = sum(float(x) ** 2 for x in sv_v[m + 1:2 * m + 1])
                e_kan = float(sv_v[2 * m + 1]) ** 2
            else:
                e_k1 = sum(float(x) ** 2 for x in sv_v[:m])
                e_d = 0.0
                e_k2 = sum(float(x) ** 2 for x in sv_v[m:2 * m])
                e_kan = float(sv_v[2 * m]) ** 2
            tot = e_k1 + e_d + e_k2 + e_kan
            print(f"{label:<14}{k1:>3}{m:>3} {variant:<7}"
                  f"{r['pv']/n:>8.3f}{r['sv']/n:>8.3f}{r['sv']/r['pv']:>8.3f}  "
                  f"{e_k1/tot:>6.3f}{e_d/tot:>6.3f}{e_k2/tot:>6.3f}{e_kan/tot:>6.3f}"
                  f"  {'Y' if r['ok'] else 'n'}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("T23-B: the falsifier -- K1 wall on the two 12-bit anchors")
    print("-" * 78)
    print("2026-07-29 T4 baseline (LLL, lattice L):")
    print("  12-bit/2557 (lam*=0.340): 5/5 up to K1=6, 4/5 at 12, 1/5 at 16, 0/5 at 24")
    print("  12-bit/2677 (lam*=0.070): 5/5 up to K1=4, 2/5 at 6, 0/5 from K1=8\n")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    b_rows = []
    for label, curve, _k1, m in hist[1:]:
        p, b, n, lam, G = curve
        print(f"{label}  n={n}  lam*={lam_star(lam,n):.3f}  m={m}  "
              f"K2={math.isqrt(n)+1}")
        hdr = f"  {'lattice':<8}" + "".join(f"{k:>6}" for k in K1_GRID)
        print(hdr)
        for variant in ('kannan', 'proj', 'cvp'):
            cells = []
            for k1 in K1_GRID:
                w, t, mr = rate(curve, m, k1, SEEDS, variant)
                cells.append(f"{w}/{t}")
                b_rows.append((label, variant, k1, w, t, mr))
            print(f"  {variant:<8}" + "".join(f"{c:>6}" for c in cells))
        # sv/pv summary at the wall
        line = []
        for variant in ('kannan', 'proj'):
            _w, _t, mr = rate(curve, m, 8, SEEDS, variant)
            line.append(f"{variant} {mr:.3f}")
        print(f"  mean sv/pv at K1=8: " + ",  ".join(line))
        print()

    # ---------------------------------------------------------------------------
    print("-" * 78)
    print("T23-C: 17-bit sweep, 20 curves, m=12, 5 seeds  (T3 protocol)")
    print("-" * 78)
    print("2026-07-29 T3 baseline (lattice L): eff=0.05 -> 19/20, 0.15 -> 3/20, 0.25 -> 0/20\n")
    sys.stdout.flush()

    curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
    print(f"found {len(curves17)} 17-bit j=0 GLV curves; "
          f"lam* in [{min(lam_star(c[3], c[2]) for c in curves17):.4f}, "
          f"{max(lam_star(c[3], c[2]) for c in curves17):.4f}]")
    sys.stdout.flush()

    M_SIGS = 12
    c_summary = {}
    for eff in (0.05, 0.15, 0.25):
        tally = {'kannan': 0, 'proj': 0}
        full = {'kannan': 0, 'proj': 0}
        trials = {'kannan': 0, 'proj': 0}
        lam_ok = {'kannan': [], 'proj': []}
        for cur in curves17:
            p, b, n, lam, G = cur
            k2_bound = math.isqrt(n) + 1
            k1_bound = max(2, round(eff * n / k2_bound))
            for variant in ('kannan', 'proj'):
                w, t, _mr = rate(cur, M_SIGS, k1_bound, SEEDS, variant)
                tally[variant] += w
                trials[variant] += t
                if w == t:
                    full[variant] += 1
                    lam_ok[variant].append(lam_star(lam, n))
        c_summary[eff] = (dict(full), dict(tally), dict(trials))
        print(f"eff={eff:<5}  "
              f"L:      {full['kannan']:>2}/{len(curves17)} curves 5/5, "
              f"{tally['kannan']:>3}/{trials['kannan']} trials   |   "
              f"L_proj: {full['proj']:>2}/{len(curves17)} curves 5/5, "
              f"{tally['proj']:>3}/{trials['proj']} trials")
        sys.stdout.flush()

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("T23-D: BKZ-20 on the hard regime (eff=0.15), L vs L_proj")
    print("-" * 78)
    sys.stdout.flush()

    eff = 0.15
    bk = {'kannan': [0, 0], 'proj': [0, 0]}
    for cur in curves17[:10]:
        p, b, n, lam, G = cur
        k2_bound = math.isqrt(n) + 1
        k1_bound = max(2, round(eff * n / k2_bound))
        for variant in ('kannan', 'proj'):
            w, t, _mr = rate(cur, M_SIGS, k1_bound, SEEDS[:3], variant, beta=20)
            bk[variant][0] += w
            bk[variant][1] += t
    print(f"BKZ-20, 10 curves x 3 seeds:  L {bk['kannan'][0]}/{bk['kannan'][1]}   "
          f"L_proj {bk['proj'][0]}/{bk['proj'][1]}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
