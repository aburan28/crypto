"""
GLV-HNP — Thread 23: is the Phase-2 K1 wall algorithmic or information-theoretic?

Context (2026-07-29 log entry, RESEARCH_AUTOLAB_LOG.md ~line 5990)
------------------------------------------------------------------
T5 of the 2026-07-29 run showed that in the Phase-2 Kannan lattice the shortest
vector after LLL is ALWAYS the trivial vector  n*S_D*e_m  (100% of its energy in
the d-column, |sv[m]|/n = 1.0000).  The planted vector is therefore never lam_1,
and recovery is a BDD/coset condition, not an SVP condition.  The proposed next
step was:

    "Thread 23 - reformulate the Phase-2 lattice so the target is lam_1.
     Concretely: project the lattice along e_m (quotient out the trivial n*e_m
     direction) and solve BDD in the projection, or replace the Kannan embedding
     with an explicit CVP call (Babai nearest-plane on the reduced basis)
     targeting (A_i*S_K1, 0, ...).
     Falsifier: if sv/pv rises above 1 after the reformulation and the K1 wall in
     T4 moves outward on the lam*=0.07 curve (currently K1~4-6), the
     reformulation is a real improvement; if the wall stays at K1~4-6, then the
     wall is information-theoretic and Phase 2 is at its ceiling."

This script executes that, and adds the arm that actually settles the question:
EXACT CVP.  If exact CVP also fails at the same K1, no reduction algorithm and no
reformulation can help - the planted vector simply is not the closest lattice
point, i.e. the wall is information-theoretic.

Experiments
-----------
E1  Analytic: Gaussian-heuristic ratios for the Kannan lattice L, its e_m
    projection, and the CVP lattice L0.  Predicts the eff = K1*K2/n wall.
E2  Three-arm K1 grid on the two 12-bit reference curves (lam*=0.340 and
    lam*=0.070) - Kannan+LLL (T4 baseline), L0+Babai, L0+exact CVP.
E3  Geometry: sv/pv in the e_m projection, and whether a lattice point strictly
    closer to the CVP target than the planted one exists (decisive test).

Lattice conventions are taken verbatim from
`glv_hnp_phase2_lambda_threshold.py:254 build_glv_lattice` so numbers are
directly comparable with the 2026-07-29 run.

Run: python3 glv_hnp_phase2_projected_cvp.py
"""

import math
import random
import time

from fpylll import IntegerMatrix, LLL, BKZ, CVP, GSO

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) - verbatim from prior scripts
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = 3 * x1 * x1 % p * modinv(2 * y1 % p, p) % p
    else:
        lam = (y2 - y1) % p * modinv((x2 - x1) % p, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def find_point(p, b, seed=42):
    rng = random.Random(seed)
    for _ in range(5000):
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is not None and y != 0:
            return (x, y)
    return None

# ---------------------------------------------------------------------------
# CM machinery for harvesting fresh j=0 GLV curves (E5 out-of-sample set)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) >= 0 with a^2 - a*b + b^2 = p."""
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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_eigenvalues(n):
    """Both roots of x^2+x+1 = 0 mod n (needs n = 1 mod 3)."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 - sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0 or (r2 * r2 + r2 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def identify_twist(p, n):
    """Smallest b with #(y^2 = x^3 + b / F_p) = n."""
    for b in range(1, 400):
        P = find_point(p, b)
        if P is None:
            continue
        if ec_mul(P, n, p) is None:
            return b
    return None

def harvest_fresh_curves(lo, hi, per_bin=2, nbins=10, max_primes=60000):
    """j=0 GLV curves with n prime, bucketed by lam* into nbins over [0,0.5].

    Bucketing on lam* (not on nu_hat) keeps the out-of-sample set independent of
    the predictor being tested.
    """
    import sympy
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 2 or n % 3 != 1 or not sympy.isprime(n):
                        continue
                    roots = glv_eigenvalues(n)
                    if roots is None:
                        continue
                    lam = roots[0]
                    ls = lam_star(lam, n)
                    idx = min(nbins - 1, int(ls / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    b = identify_twist(p, n)
                    if b is None:
                        continue
                    G = find_generator(p, b, n)
                    if G is None:
                        continue
                    bins[idx].append({'tag': f'p{p}', 'p': p, 'b': b, 'n': n,
                                      'lam': lam, 'lam*': ls, 'G': G})
                    break
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    return [c for i in range(nbins) for c in bins[i]]

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim: glv_hnp_phase2_lambda_threshold.py:227-284)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

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

def lambda_block_mu(n, lam, S_K1, S_K2):
    """Exact |shortest vector| of the 2D rival block L2 (Thread 20c nu numerator)."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w

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
    """Kannan-embedded Phase-2 lattice, dim 2m+2."""
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

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

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
# NEW - Thread 23: the CVP lattice L0 (Kannan row removed) and its target
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """L0, dim 2m+1.  Columns: [0..m-1] k1, [m] d, [m+1..2m] k2.

    Rows:  n*S_K1*e_i              (i < m)
           (B_i*S_K1, ..., S_D, 0...)
           (-lam*S_K1*e_i, 0, S_K2*e_i)

    Contains  u = w_planted - t  where t = (A_i*S_K1, 0, ..., 0), because
    k1_i*S_K1 = (A_i + d*B_i - lam*k2_i + c_i*n)*S_K1.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    return M

def cvp_target(sigs, n, k1_bound, k2_bound):
    """-t = (-A_i*S_K1, 0, ..., 0).  CVP solution u satisfies u-(-t) = w_planted."""
    m = len(sigs)
    S_K1, _, _, _ = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m + 1)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return t

def planted_vector_cvp(sigs, d_secret, n, k1_bound, k2_bound):
    """w_planted = (k1_i*S_K1, d*S_D, k2_i*S_K2), dim 2m+1 (no Kannan coord)."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    return v

def mat_rows(A, dim):
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def combine(rows, coeffs):
    dim = len(rows[0])
    out = [0] * dim
    for c, r in zip(coeffs, rows):
        if c:
            for j in range(dim):
                out[j] += c * r[j]
    return out

# ---------------------------------------------------------------------------
# The three arms
# ---------------------------------------------------------------------------

def arm_kannan_lll(sigs, curve_n, lam, k1_bound, k2_bound, d_secret,
                   use_bkz=False, beta=20):
    """T4 baseline: Kannan embedding + LLL (or BKZ)."""
    m = len(sigs)
    dim = 2 * m + 2
    M = build_glv_lattice(sigs, curve_n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    red = mat_rows(A, dim)
    _, _, _, S_KAN = scales(curve_n, k1_bound, k2_bound)
    ok = recover_d(red, m, curve_n, S_KAN, d_secret) is not None
    pv = norm(planted_vector(sigs, d_secret, curve_n, k1_bound, k2_bound))
    sv = min(norm(r) for r in red)
    return ok, pv, sv

def arm_babai(sigs, curve_n, lam, k1_bound, k2_bound, d_secret,
              use_bkz=False, beta=20):
    """L0 + Babai nearest-plane on the LLL/BKZ-reduced basis."""
    m = len(sigs)
    dim = 2 * m + 1
    M = build_cvp_lattice(sigs, curve_n, lam, k1_bound, k2_bound)
    t = cvp_target(sigs, curve_n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = mat_rows(A, dim)
    G = GSO.Mat(A)
    G.update_gso()
    coeffs = G.babai([int(x) for x in t])
    u = combine(rows, coeffs)
    w = [u[j] - t[j] for j in range(dim)]     # = u + (A_i*S_K1, 0, ...)
    _, S_D, _, _ = scales(curve_n, k1_bound, k2_bound)
    d_cand = (u[m] // S_D) % curve_n
    ok = (d_cand == d_secret)
    return ok, norm(w)

def arm_exact_cvp(sigs, curve_n, lam, k1_bound, k2_bound, d_secret, beta=None):
    """L0 + exact CVP (fpylll enumeration).  Decisive arm."""
    m = len(sigs)
    dim = 2 * m + 1
    M = build_cvp_lattice(sigs, curve_n, lam, k1_bound, k2_bound)
    t = cvp_target(sigs, curve_n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
        LLL.reduction(A)
    u = list(CVP.closest_vector(A, tuple(int(x) for x in t)))
    w = [u[j] - t[j] for j in range(dim)]
    _, S_D, _, _ = scales(curve_n, k1_bound, k2_bound)
    d_cand = (u[m] // S_D) % curve_n
    ok = (d_cand == d_secret)
    return ok, norm(w)

# ---------------------------------------------------------------------------
# Analytic predictions (E1)
# ---------------------------------------------------------------------------

def gh(dim, log_det):
    """Gaussian heuristic: sqrt(dim/(2*pi*e)) * det^(1/dim), via log det."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def analytic(n, m, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    # Kannan lattice L: dim 2m+2
    dimL = 2 * m + 2
    logdetL = (m * math.log(n * S_K1) + math.log(S_D)
               + m * math.log(S_K2) + math.log(S_KAN))
    pv2 = (m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
           + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    pv = math.sqrt(pv2)
    ghL = gh(dimL, logdetL)
    # e_m projection of L: dim 2m+1, det/(n*S_D)
    dimP = 2 * m + 1
    logdetP = logdetL - math.log(n * S_D)
    pvP = math.sqrt(pv2 - (n * S_D) ** 2 / 3.0)
    ghP = gh(dimP, logdetP)
    # CVP lattice L0: dim 2m+1, no Kannan row/col
    dim0 = 2 * m + 1
    logdet0 = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2)
    dist0 = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
                      + m * (k2_bound * S_K2) ** 2 / 3.0)
    gh0 = gh(dim0, logdet0)
    return {
        'eff': k1_bound * k2_bound / n,
        'pv/gh_L': pv / ghL,
        'pvP/gh_P': pvP / ghP,
        'dist0/gh0': dist0 / gh0,
        'dist0/(gh0/2)': dist0 / (gh0 / 2.0),
        'trivial/gh_L': (n * S_D) / ghL,
    }

# ---------------------------------------------------------------------------
# Reference curves (2026-06-15 / 2026-07-29)
# ---------------------------------------------------------------------------

def ref_curves():
    out = []
    for (p, b, n, lam, tag) in [
        (211, 2, 199, 106, '8-bit/199'),
        (2557, 2, 2659, 1755, '12-bit/2557'),
        (2677, 2, 2647, 185, '12-bit/2677'),
    ]:
        G = find_generator(p, b, n)
        assert G is not None, f'no generator for p={p}'
        out.append({'tag': tag, 'p': p, 'b': b, 'n': n, 'lam': lam, 'G': G,
                    'lam*': lam_star(lam, n)})
    return out

SEEDS = [1, 2, 3, 4, 5]


def main():
    t_start = time.time()
    print('=' * 78)
    print('Thread 23 - projected / CVP reformulation of the Phase-2 GLV-HNP lattice')
    print('=' * 78)

    curves = ref_curves()
    for c in curves:
        print(f"  {c['tag']:<12} p={c['p']:<5} n={c['n']:<5} lam={c['lam']:<5} "
              f"lam*={c['lam*']:.4f}  K2=isqrt(n)+1={math.isqrt(c['n']) + 1}")

    # ---------------- E1 ----------------
    print('\n' + '=' * 78)
    print('E1 - analytic uniqueness / BDD ratios  (m=12, curve 12-bit/2677)')
    print('=' * 78)
    print('  pv/gh_L      : planted norm / Gaussian heuristic of the Kannan lattice')
    print('  pvP/gh_P     : same, in the e_m projection (Thread 23 reformulation)')
    print('  dist0/gh0    : CVP target distance / GH of L0  (BDD: need < ~1)')
    print('  triv/gh_L    : trivial vector n*S_D / GH of L  (< 1 => L is not random)')
    print()
    n_ref = 2647
    k2_ref = math.isqrt(n_ref) + 1
    print(f"  {'K1':>4} {'eff':>7} {'pv/gh_L':>9} {'pvP/gh_P':>9} "
          f"{'dist0/gh0':>10} {'dist0/(gh0/2)':>14} {'triv/gh_L':>10}")
    for k1 in [2, 3, 4, 6, 8, 12, 16, 24]:
        a = analytic(n_ref, 12, k1, k2_ref)
        print(f"  {k1:>4} {a['eff']:>7.3f} {a['pv/gh_L']:>9.3f} "
              f"{a['pvP/gh_P']:>9.3f} {a['dist0/gh0']:>10.3f} "
              f"{a['dist0/(gh0/2)']:>14.3f} {a['trivial/gh_L']:>10.3f}")

    print('\n  Cross-check against T3 (2026-07-29, 17-bit sweep, m=12):')
    for eff, obs in [(0.05, '19/20'), (0.15, '3/20'), (0.25, '0/20')]:
        n17 = 65269
        k2_17 = math.isqrt(n17) + 1
        k1_17 = max(2, int(round(eff * n17 / k2_17)))
        a = analytic(n17, 12, k1_17, k2_17)
        print(f"    eff~{eff:.2f} (K1={k1_17:>3}, actual eff={a['eff']:.3f}): "
              f"pv/gh_L={a['pv/gh_L']:.3f}  pvP/gh_P={a['pvP/gh_P']:.3f}  "
              f"dist0/gh0={a['dist0/gh0']:.3f}   observed {obs}")

    # ---------------- E2 ----------------
    print('\n' + '=' * 78)
    print('E2 - three-arm K1 grid, m=12, 5 seeds')
    print('=' * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    D_SECRET_OFFSET = 7   # d = n//3 + 7, same style as prior runs

    results = {}
    for c in curves[1:]:                       # the two 12-bit curves
        n = c['n']
        k2b = math.isqrt(n) + 1
        d_secret = n // 3 + D_SECRET_OFFSET
        print(f"\n--- {c['tag']}  (lam*={c['lam*']:.3f}, n={n}, K2={k2b}, "
              f"d={d_secret}) ---")
        print(f"  {'K1':>4} {'eff':>6} | {'Kannan+LLL':>11} {'L0+Babai':>10} "
              f"{'L0+exactCVP':>12} | {'sv/pv':>7} {'w_bab/w_pl':>11} "
              f"{'w_cvp/w_pl':>11} {'closer?':>8}")
        for k1 in K1_GRID:
            wins = [0, 0, 0]
            svpv, rb, rc, closer = [], [], [], 0
            for seed in SEEDS:
                sigs = gen_signatures(c['G'], d_secret, 12, n, c['lam'], c['p'],
                                      k1, k2b, seed)
                if len(sigs) < 12:
                    continue
                ok_a, pv, sv = arm_kannan_lll(sigs, n, c['lam'], k1, k2b, d_secret)
                ok_b, wb = arm_babai(sigs, n, c['lam'], k1, k2b, d_secret)
                ok_c, wc = arm_exact_cvp(sigs, n, c['lam'], k1, k2b, d_secret)
                w_pl = norm(planted_vector_cvp(sigs, d_secret, n, k1, k2b))
                wins[0] += ok_a; wins[1] += ok_b; wins[2] += ok_c
                svpv.append(sv / pv); rb.append(wb / w_pl); rc.append(wc / w_pl)
                if wc < w_pl - 1e-9:
                    closer += 1
            ns = len(svpv)
            eff = k1 * k2b / n
            print(f"  {k1:>4} {eff:>6.3f} | {wins[0]:>6}/{ns:<4} "
                  f"{wins[1]:>5}/{ns:<4} {wins[2]:>7}/{ns:<4} | "
                  f"{sum(svpv)/ns:>7.3f} {sum(rb)/ns:>11.3f} "
                  f"{sum(rc)/ns:>11.3f} {closer:>5}/{ns}")
            results[(c['tag'], k1)] = (wins, ns)

    # ---------------- E3 ----------------
    print('\n' + '=' * 78)
    print('E3 - is the planted vector the closest lattice point?')
    print('=' * 78)
    print('  "closer" above counts seeds where exact CVP found a lattice point')
    print('  strictly closer to the target than the planted one.  Any such seed')
    print('  is an information-theoretic failure: no algorithm can recover d.')
    print()
    print('  BKZ-40 control on the two arms at the wall (K1=8, 12-bit/2677):')
    c = curves[2]
    n = c['n']; k2b = math.isqrt(n) + 1; d_secret = n // 3 + D_SECRET_OFFSET
    for k1 in [4, 6, 8]:
        wa = wb = wc = ns = 0
        for seed in SEEDS:
            sigs = gen_signatures(c['G'], d_secret, 12, n, c['lam'], c['p'],
                                  k1, k2b, seed)
            if len(sigs) < 12:
                continue
            ns += 1
            ok_a, _, _ = arm_kannan_lll(sigs, n, c['lam'], k1, k2b, d_secret,
                                        use_bkz=True, beta=40)
            ok_b, _ = arm_babai(sigs, n, c['lam'], k1, k2b, d_secret,
                                use_bkz=True, beta=40)
            ok_c, _ = arm_exact_cvp(sigs, n, c['lam'], k1, k2b, d_secret, beta=40)
            wa += ok_a; wb += ok_b; wc += ok_c
        print(f"    K1={k1:>3}: Kannan+BKZ40 {wa}/{ns}   L0+Babai(BKZ40) {wb}/{ns}"
              f"   L0+exactCVP(BKZ40) {wc}/{ns}")

    # ---------------- E4 ----------------
    print('\n' + '=' * 78)
    print('E4 - which predictor locates the wall: pv/gh_L or nu_hat?')
    print('=' * 78)
    print('  nu_hat = lambda_1(L2)/sqrt(det L2) for the 2D rival block L2 =')
    print('  <(n*S_K1,0), (-lam*S_K1,S_K2)>   (Thread 20c, AUC 0.935 on C1/C2).')
    print('  Restored dependency: glv_hnp_phase2_nuhat_base.py')
    print()
    print(f"  {'curve':<12} {'lam*':>6} {'K1':>4} {'eff':>6} {'pv/gh_L':>8} "
          f"{'nu_hat':>7} | {'LLL':>5} {'CVP':>5}")
    for c in curves[1:]:
        n = c['n']; k2b = math.isqrt(n) + 1
        for k1 in K1_GRID:
            S_K1, _, S_K2, _ = scales(n, k1, k2b)
            mu_norm, _ = lambda_block_mu(n, c['lam'], S_K1, S_K2)
            nu_hat = mu_norm / math.sqrt(n * S_K1 * S_K2)
            a = analytic(n, 12, k1, k2b)
            wins, ns = results[(c['tag'], k1)]
            print(f"  {c['tag']:<12} {c['lam*']:>6.3f} {k1:>4} {a['eff']:>6.3f} "
                  f"{a['pv/gh_L']:>8.3f} {nu_hat:>7.3f} | "
                  f"{wins[0]:>2}/{ns:<2} {wins[2]:>2}/{ns:<2}")

    # ---------------- E5 ----------------
    print('\n' + '=' * 78)
    print('E5 - OUT-OF-SAMPLE test of the two-factor predictor  P = (pv/gh_L)*nu_hat')
    print('=' * 78)
    print('  E4 suggests P < 1 <=> recovery, on 16 in-sample rows / 2 curves.')
    print('  E4 fitted the exponent (1) and threshold (1) post hoc, so this')
    print('  section tests P on FRESH curves that played no part in choosing it.')
    print('  Baselines: pv/gh_L alone (lam-independent) and nu_hat alone.')

    fresh = harvest_fresh_curves(lo=1 << 16, hi=1 << 17, per_bin=2, nbins=10)
    print(f"\n  harvested {len(fresh)} fresh 16-17-bit j=0 GLV curves, "
          f"lam* in [{min(c['lam*'] for c in fresh):.3f}, "
          f"{max(c['lam*'] for c in fresh):.3f}]")

    rows = []
    M_SIGS = 12
    OOS_SEEDS = [11, 12, 13]
    for c in fresh:
        n = c['n']; k2b = math.isqrt(n) + 1
        d_secret = n // 3 + D_SECRET_OFFSET
        for k1 in [4, 8, 16, 32, 48, 64]:
            if k1 * k2b >= n:
                continue
            S_K1, _, S_K2, _ = scales(n, k1, k2b)
            mu_norm, _ = lambda_block_mu(n, c['lam'], S_K1, S_K2)
            nu_hat = mu_norm / math.sqrt(n * S_K1 * S_K2)
            a = analytic(n, M_SIGS, k1, k2b)
            wins = tot = 0
            for seed in OOS_SEEDS:
                sigs = gen_signatures(c['G'], d_secret, M_SIGS, n, c['lam'],
                                      c['p'], k1, k2b, seed)
                if len(sigs) < M_SIGS:
                    continue
                ok, _, _ = arm_kannan_lll(sigs, n, c['lam'], k1, k2b, d_secret)
                wins += ok; tot += 1
            if tot == 0:
                continue
            rows.append({'tag': c['tag'], 'lam*': c['lam*'], 'k1': k1,
                         'eff': a['eff'], 'pvgh': a['pv/gh_L'], 'nu': nu_hat,
                         'P': a['pv/gh_L'] * nu_hat,
                         'p_hat': wins / tot, 'win': wins == tot})

    print(f"  {len(rows)} out-of-sample (curve, K1) cells, {OOS_SEEDS} seeds each\n")

    def auc(rows, key, positive_is_small=True):
        pos = [r[key] for r in rows if r['win']]
        neg = [r[key] for r in rows if not r['win']]
        if not pos or not neg:
            return float('nan')
        c = sum((1.0 if ((a < b) == positive_is_small) else
                 0.5 if a == b else 0.0) for a in pos for b in neg)
        return c / (len(pos) * len(neg))

    print(f"  {'predictor':<12} {'AUC':>7}   (1.0 = perfect ordering of "
          f"5/5-success vs not)")
    for key in ['P', 'pvgh', 'nu', 'eff']:
        print(f"  {key:<12} {auc(rows, key):>7.3f}")

    tp = sum(1 for r in rows if r['P'] < 1 and r['win'])
    fp = sum(1 for r in rows if r['P'] < 1 and not r['win'])
    fn = sum(1 for r in rows if r['P'] >= 1 and r['win'])
    tn = sum(1 for r in rows if r['P'] >= 1 and not r['win'])
    n_all = tp + fp + fn + tn
    print(f"\n  Confusion matrix for the PRE-REGISTERED rule 'P < 1 => 3/3 recovery':")
    print(f"    P<1  & recovered : {tp:>4}      P<1  & failed : {fp:>4}")
    print(f"    P>=1 & recovered : {fn:>4}      P>=1 & failed : {tn:>4}")
    print(f"    accuracy = {(tp + tn) / n_all:.3f}   "
          f"(majority baseline = {max(tp + fn, fp + tn) / n_all:.3f})")

    print(f"\n  Fairness check - each predictor at its BEST-FIT threshold")
    print(f"  (favours the baselines: P is judged at the pre-registered 1.0):")
    for key in ['P', 'pvgh', 'nu', 'eff']:
        vals = sorted(set(r[key] for r in rows))
        best, best_t = 0.0, None
        for t in vals:
            acc = sum(1 for r in rows if (r[key] < t) == r['win']) / len(rows)
            if acc > best:
                best, best_t = acc, t
        pre = (sum(1 for r in rows if (r[key] < 1.0) == r['win']) / len(rows)
               if key in ('P', 'pvgh') else float('nan'))
        print(f"    {key:<6} best acc={best:.3f} at threshold {best_t:.3f}"
              + (f"   (acc at 1.0 = {pre:.3f})" if key in ('P', 'pvgh') else ''))

    print(f"\n  p_hat by P decile:")
    srt = sorted(rows, key=lambda r: r['P'])
    nb = 10
    print(f"    {'decile':>7} {'P mid':>7} {'mean p_hat':>11} {'n':>4}")
    for i in range(nb):
        chunk = srt[i * len(srt) // nb:(i + 1) * len(srt) // nb]
        if not chunk:
            continue
        mid = chunk[len(chunk) // 2]['P']
        print(f"    {i + 1:>7} {mid:>7.3f} "
              f"{sum(r['p_hat'] for r in chunk) / len(chunk):>11.3f} "
              f"{len(chunk):>4}")

    print(f"\n[total runtime {time.time() - t_start:.1f}s]")


if __name__ == '__main__':
    main()
