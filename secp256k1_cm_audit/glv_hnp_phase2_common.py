"""
Shared helpers for the GLV-aware HNP Phase 2 experiments (Thread 5 / Thread 20).

Extracted from glv_hnp_phase2_lambda_threshold.py so that the experiment scripts
can import the machinery without executing each other's experiment bodies.
Contains: minimal j=0 EC arithmetic, Eisenstein CM curve search, GLV eigenvalue
and centered-rho computation, biased-nonce signature generation, the
column-scaled attack lattice, and oracle-free candidate verification.
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
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
    """Find a generator of E: y^2=x^3+b over F_p, assuming #E = n prime."""
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
# CM theory: Eisenstein decomposition for j=0 curves over F_p
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) >= 0 with a^2 - a*b + b^2 = p.  O(sqrt(p))."""
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

def glv_eigenvalues(n):
    """Both roots of x^2 + x + 1 = 0 mod n (requires n = 1 mod 3)."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 - sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0 or (r2 * r2 + r2 + 1) % n != 0:
        return None
    return (r1, r2)

def centered_rho(lam, n):
    """rho = min(lam, n-lam)/n in (0, 1/2]."""
    return min(lam % n, n - (lam % n)) / n

def identify_twist(p, n, seed=777):
    """Find b with #(y^2 = x^3 + b over F_p) = n, by testing point orders."""
    rng = random.Random(seed)
    for b in range(1, min(p, 400)):
        for _ in range(60):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                if ec_mul((x, y), n, p) is None:
                    return b
                break
    return None

def find_glv_curves(p_lo, p_hi, rho_lo=0.0, rho_hi=0.51, want=6, seed=777):
    """
    Collect j=0 GLV curves with prime p in [p_lo,p_hi), prime order n = 1 mod 3,
    and centered rho in [rho_lo, rho_hi).  Returns list of dicts.
    """
    out = []
    p = sympy.nextprime(p_lo - 1)
    while p < p_hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n = p + 1 - t
                    if n < 8 or not sympy.isprime(n) or n % 3 != 1:
                        continue
                    lams = glv_eigenvalues(n)
                    if lams is None:
                        continue
                    lam = min(lams)
                    rho = centered_rho(lam, n)
                    if not (rho_lo <= rho < rho_hi):
                        continue
                    b = identify_twist(p, n, seed)
                    if b is None:
                        continue
                    G = find_generator(p, b, n, seed)
                    if G is None:
                        continue
                    out.append(dict(p=p, b=b, n=n, lam=lam, lam2=max(lams),
                                    rho=rho, G=G))
                    break
        p = sympy.nextprime(p)
    return out

# ---------------------------------------------------------------------------
# Signature generation with a biased GLV nonce
# ---------------------------------------------------------------------------

def gen_signatures(curve, d_secret, m, k1_bound, k2_bound, seed=42):
    p, b, n, lam, G = curve['p'], curve['b'], curve['n'], curve['lam'], curve['G']
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 400000:
        attempts += 1
        k1 = rng.randrange(k1_bound)
        k2 = rng.randrange(k2_bound)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randrange(n)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append(dict(A=A, B=B, k1=k1, k2=k2))
    return sigs

# ---------------------------------------------------------------------------
# Lattice construction (column-diagonal scaling, validated 2026-07-26)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (max(1, n // k1_bound), 1, max(1, n // k2_bound), n)

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def planted_vector(sigs, n, d, k1_bound, k2_bound):
    """The vector the attack is trying to find, in scaled coordinates."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [s['k1'] * S_K1 for s in sigs]
    v.append(d * S_D)
    v += [s['k2'] * S_K2 for s in sigs]
    v.append(S_KAN)
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def log2_det(m, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (m * math.log2(n * S_K1) + math.log2(S_D)
            + m * math.log2(S_K2) + math.log2(S_KAN))

# ---------------------------------------------------------------------------
# Oracle-free verification
# ---------------------------------------------------------------------------

def verify_d(d_cand, sigs, n, lam, k1_bound, k2_bound):
    """
    Accept d_cand only if EVERY signature's implied nonce k = A + B*d admits a
    GLV decomposition k = k1 + lam*k2 with k1 < K1 and k2 < K2.  No secret used.
    """
    if not (0 < d_cand < n):
        return False
    for s in sigs:
        k_full = (s['A'] + s['B'] * d_cand) % n
        ok = False
        for k2 in range(k2_bound):
            if (k_full - lam * k2) % n < k1_bound:
                ok = True
                break
        if not ok:
            return False
    return True


# ---------------------------------------------------------------------------
# Rank-2 spurious-vector machinery and the attack-instance runner (Thread 20)
# ---------------------------------------------------------------------------

from fractions import Fraction

M_MIN, M_MAX = 3, 16
SEEDS = [42, 1234, 9999]

def gauss_reduce(b1, b2):
    """Lagrange-Gauss reduction of a rank-2 integer lattice. Exact arithmetic."""
    dot = lambda u, w: u[0] * w[0] + u[1] * w[1]
    while True:
        if dot(b1, b1) > dot(b2, b2):
            b1, b2 = b2, b1
        q = Fraction(dot(b1, b2), dot(b1, b1))
        mu = math.floor(q + Fraction(1, 2))
        if mu == 0:
            return b1, b2
        b2 = (b2[0] - mu * b1[0], b2[1] - mu * b1[1])

def spurious_mu(n, lam, k1_bound, k2_bound):
    """lambda_1(L2) and the (u,v) attaining it."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    b1, b2 = gauss_reduce((n * S_K1, 0), (lam * S_K1, S_K2))
    short = b1 if (b1[0]**2 + b1[1]**2) <= (b2[0]**2 + b2[1]**2) else b2
    return math.sqrt(float(short[0])**2 + float(short[1])**2), \
           (short[0] // S_K1, short[1] // S_K2)

def m_info(n, k1_bound, k2_bound):
    eff = k1_bound * k2_bound / n
    return (math.ceil(math.log(n) / math.log(1.0 / eff))
            if eff < 1.0 else float('inf'))

# ---------------------------------------------------------------------------
# Instance runner.  Signatures are generated once at M_MAX per (curve,K1,seed)
# and sliced, so sigs[:m] is exactly what gen_signatures(...,m) would return.
# ---------------------------------------------------------------------------

def reduce_and_test(curve, sigs, m, d, k1_bound, algo='LLL', beta=20):
    n, lam = curve['n'], curve['lam']
    k2_bound = math.isqrt(n) + 1
    sigs = sigs[:m]
    dim = 2 * m + 2
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    v = planted_vector(sigs, n, d, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if algo == 'LLL':
        LLL.reduction(A)
    else:
        BKZ.reduction(A, BKZ.Param(beta))
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    ok = False
    for row in rows:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        if verify_d((sign * row[m]) % n, sigs, n, lam, k1_bound, k2_bound):
            ok = True
            break
    n_spur = 0
    for row in rows:
        if row[m] != 0 or row[dim - 1] != 0:
            continue
        i1 = [i for i in range(m) if row[i] != 0]
        i2 = [i for i in range(m) if row[m + 1 + i] != 0]
        if len(i1) <= 1 and len(i2) <= 1 and (i1 or i2):
            if not i1 or not i2 or i1[0] == i2[0]:
                n_spur += 1
    return ok, norm(v), n_spur, dim

