"""
GLV-HNP Phase 2 — shared helpers.

Extracted VERBATIM from `glv_hnp_phase2_lambda_threshold.py` (Thread 20,
2026-07-29) lines 83-411 so that Thread 23 reuses the exact same EC
arithmetic, curve search, signature generation and Phase-2 lattice
construction.  Nothing here is new; do not edit without re-running
`glv_hnp_phase2_lambda_threshold.py` and diffing its output.
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)  -- same as glv_hnp_phase2_20bit.py
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

# ---------------------------------------------------------------------------
# The lambda-block 2D sublattice and its exact shortest vector
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Lagrange/Gauss reduction of a 2D integer lattice.  Returns the exact
    shortest nonzero vector (u, v are tuples of ints)."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        # mu = round(<v,u>/<u,u>), exact integer rounding (round-half-away)
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def lambda_block_mu(n, lam, S_K1, S_K2):
    """Exact |shortest vector| of the 2D block < (n*S_K1,0), (-lam*S_K1,S_K2) >."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w

def planted_norm_expected(m, n, k1_bound, k2_bound, S_K1, S_D, S_K2, S_KANNAN):
    """E[||v_planted||] with k1~U[0,K1), k2~U[0,K2), d~U[0,n).  E[x^2]=X^2/3."""
    return math.sqrt(
        m * (k1_bound * S_K1) ** 2 / 3.0
        + (n * S_D) ** 2 / 3.0
        + m * (k2_bound * S_K2) ** 2 / 3.0
        + S_KANNAN ** 2
    )

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim construction from glv_hnp_phase2_20bit.py)
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
    """lam is used verbatim (NOT reduced mod n) so that representation
    invariance can be tested."""
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
    """The exact planted lattice vector (k1_i*S_K1, d*S_D, k2_i*S_K2, S_KANNAN)."""
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
    return math.sqrt(sum(x * x for x in v))


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

def run_experiment(curve, m, d_secret, k1_bound, seed=42, lam_override=None,
                   use_bkz=False, bkz_beta=20):
    """Returns (recovered: bool, planted_norm, shortest_reduced_norm)."""
    p, b, n, lam, G = curve
    lam_lat = lam if lam_override is None else lam_override
    k2_bound = math.isqrt(n) + 1
    # signatures always use the canonical lam (the true GLV eigenvalue);
    # only the lattice representation varies.
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam_lat, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KANNAN = scales(n, k1_bound, k2_bound)
    ok = recover_d(reduced, m, n, S_KANNAN, d_secret) is not None
    pn = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
    sn = min(norm(r) for r in reduced)
    return (ok, pn, sn)

def success_rate(curve, m, k1_bound, seeds, lam_override=None,
                 use_bkz=False, bkz_beta=20):
    """Returns (wins, total, mean(shortest/planted) norm ratio)."""
    p, b, n, lam, G = curve
    wins = 0
    ratios = []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = run_experiment(curve, m, d_trial, k1_bound, seed,
                             lam_override=lam_override,
                             use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        ok, pn, sn = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    mean_ratio = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mean_ratio

def hnf_fingerprint(M):
    """Row-style Hermite normal form via sympy, as a lattice fingerprint."""
    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form
    return hermite_normal_form(Matrix(M).T).T

# ---------------------------------------------------------------------------
# Curve construction helpers
# ---------------------------------------------------------------------------

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

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    """Collect j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3,
    bucketed by lam* = min(lam,n-lam)/n into `nbins` bins over [0,0.5]."""
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


