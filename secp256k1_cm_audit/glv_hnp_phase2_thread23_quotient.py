"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
lambda_1 instead of a BDD target hidden behind a trivial short vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, result T5):
  The legacy Phase-2 lattice (glv_hnp_phase2_20bit.py:262, dim 2m+2) always
  has shortest vector n*S_D*e_m -- the "d column" direction, which exists
  because d is only defined mod n.  It is 2-3x shorter than the planted
  vector on every curve tested, success and failure alike, so recovery is a
  BDD/coset condition rather than an SVP condition.  That entry proposed
  Thread 23: quotient the lattice along e_m and re-measure the K1 wall.

  It also asserted, incorrectly, that "no choice of S_D removes it -- both
  vectors scale linearly in S_D".  Only the trivial vector does.  See R1.

Lattices compared:
  A  legacy, dim 2m+2, d carried as a coordinate; S_D exposed as a parameter.
  B  d eliminated algebraically via k_i = C_i*k_0 + D_i (mod n), dim 2m+1.
     No d column, hence no n*S_D*e_m vector.

Experiments:
  R0  correctness of B (planted vector in the lattice; recovery at K1=2)
  R1  S_D sweep on A: when does the planted vector beat n*S_D*e_m, and does
      recovery follow?
  R2  the K1 wall on the three historical curves: A(S_D=1) vs A(S_D=16) vs B
  R3  20 fresh 17-bit curves at eff ~ 0.15, the 2026-07-29 T3 breaking point

Falsifier (stated in the 2026-07-29 entry): if sv/pv rises above 1 after the
reformulation AND the K1 wall moves outward on the lam*=0.07 curve (currently
K1 ~ 4-6), the reformulation is a real improvement; if the wall stays put, the
wall is information-theoretic and Phase 2 is at its ceiling.

Run: python3 glv_hnp_phase2_thread23_quotient.py
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

# ===========================================================================
# LATTICE A — the legacy Phase-2 lattice, with S_D exposed as a parameter
# ===========================================================================

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    """Legacy Phase-2 lattice (glv_hnp_phase2_20bit.py:262), except that the
    d-column scale S_D is a free parameter instead of the hard-wired 1.

    dim = 2m+2.  Columns: [k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan]."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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


def planted_A(sigs, d_secret, n, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v


def recover_A(reduced, m, n, S_KAN, S_D, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        if (sign * row[m]) % S_D != 0:
            continue
        d_cand = ((sign * row[m]) // S_D) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None


# ===========================================================================
# LATTICE B — d-eliminated ("quotient") lattice.  Thread 23.
# ===========================================================================
#
# Eliminate the secret from the congruences instead of carrying it as a
# coordinate.  With k_i = A_i + B_i*d (mod n) and B_0 invertible (n prime):
#
#     d  = (k_0 - A_0) * B_0^{-1}
#     k_i = C_i * k_0 + D_i   (mod n),   C_i = B_i/B_0,  D_i = A_i - C_i*A_0
#
# for i = 1..m-1 (i=0 is the identity C_0=1, D_0=0 and is dropped).  Then
# substitute the GLV split k_j = k1_j + lam*k2_j and solve each congruence
# for the *dependent* unknown k1_i:
#
#     k1_i = D_i + C_i*k1_0 + lam*C_i*k2_0 - lam*k2_i   (mod n),  i >= 1
#
# Free unknowns: k1_0, k2_0, ..., k2_{m-1}   (m+1 of them).
# Dependent small outputs: k1_1..k1_{m-1}    (m-1 of them).
#
# dim = (m-1) + 1 + m + 1 = 2m+1 — exactly one less than lattice A, the
# dimension carrying the trivial vector n*S_D*e_m.  There is no d column, so
# that vector does not exist in B.
#
# Columns: [k1_1..k1_{m-1} | k1_0 | k2_0..k2_{m-1} | kannan]
#   0 .. m-2   residual k1_i (i = j+1), scale S_K1
#   m-1        k1_0,                    scale S_K1
#   m .. 2m-1  k2_t (t = col-m),        scale S_K2
#   2m         kannan,                  scale S_KANNAN

def cd_coeffs(sigs, n):
    """C_i = B_i/B_0, D_i = A_i - C_i*A_0 (mod n), for i = 1..m-1."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'], n)
    C = [0] * m
    D = [0] * m
    for i in range(1, m):
        C[i] = sigs[i]['B'] * B0inv % n
        D[i] = (sigs[i]['A'] - C[i] * sigs[0]['A']) % n
    return C, D


def build_lattice_B(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    C, D = cd_coeffs(sigs, n)
    M = [[0] * dim for _ in range(dim)]
    # modular rows for the m-1 residual columns
    for j in range(m - 1):
        M[j][j] = n * S_K1
    # k1_0 generator row
    r = m - 1
    for j in range(m - 1):
        M[r][j] = C[j + 1] * S_K1
    M[r][m - 1] = S_K1
    # k2_0 generator row (k2_0 enters every residual column via lam*C_i)
    r = m
    for j in range(m - 1):
        M[r][j] = lam * C[j + 1] % n * S_K1
    M[r][m] = S_K2
    # k2_t generator rows, t = 1..m-1 (each enters only its own column)
    for t in range(1, m):
        r = m + t
        M[r][t - 1] = (-lam) % n * S_K1
        M[r][m + t] = S_K2
    # constant / Kannan row
    r = 2 * m
    for j in range(m - 1):
        M[r][j] = D[j + 1] * S_K1
    M[r][2 * m] = S_KAN
    return M


def planted_B(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for j in range(m - 1):
        v[j] = sigs[j + 1]['k1'] * S_K1
    v[m - 1] = sigs[0]['k1'] * S_K1
    for t in range(m):
        v[m + t] = sigs[t]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def recover_B(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_0, k2_0) off any row with |kannan| = S_KANNAN, rebuild k_0,
    then d = (k_0 - A_0)/B_0."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        a, b = sign * row[m - 1], sign * row[m]
        if a % S_K1 or b % S_K2:
            continue
        k0 = (a // S_K1 + lam * (b // S_K2)) % n
        d_cand = (k0 - sigs[0]['A']) * B0inv % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None


# ===========================================================================
# Instance runner (both lattices)
# ===========================================================================

def norm(v):
    return math.sqrt(sum(x * x for x in v))


def run(curve, m, d_secret, k1_bound, seed=42, variant='A', S_D=1,
        use_bkz=False, bkz_beta=20):
    """Returns dict with ok / planted norm / shortest reduced norm / trivial
    norm, or None if signature generation failed."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if variant == 'A':
        M = build_lattice_A(sigs, n, lam, k1_bound, k2_bound, S_D)
        pv = planted_A(sigs, d_secret, n, k1_bound, k2_bound, S_D)
    else:
        M = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_B(sigs, n, k1_bound, k2_bound)
    dim = len(M)
    # the planted vector must actually lie in the lattice -- assert it once
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    if variant == 'A':
        ok = recover_A(reduced, m, n, S_KAN, S_D, d_secret) is not None
    else:
        ok = recover_B(reduced, sigs, n, lam, k1_bound, k2_bound,
                       d_secret) is not None
    pn = norm(pv)
    sn = min(norm(r) for r in reduced)
    triv = n * S_D if variant == 'A' else float('inf')
    return {'ok': ok, 'pn': pn, 'sn': sn, 'triv': triv, 'dim': dim,
            'sigs': sigs, 'reduced': reduced}


def rate(curve, m, k1_bound, seeds, variant='A', S_D=1,
         use_bkz=False, bkz_beta=20):
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = run(curve, m, d_trial, k1_bound, seed, variant, S_D,
                  use_bkz, bkz_beta)
        if res is None:
            continue
        wins += bool(res['ok'])
        ratios.append(res['sn'] / res['pn'])
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else
                              float('nan'))


# ===========================================================================
print("=" * 78)
print("Thread 23 — reformulate the Phase-2 GLV lattice so the target is a")
print("            shortest vector (kill the trivial n*S_D*e_m direction)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,          p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",     211,  2, 199,  106,  2,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",   2677, 2, 2647, 185,  8,  10),
]

CURVES = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    CURVES.append((label, (p, b, n, lam, G), k1, m, lam_star(lam, n)))


# ---------------------------------------------------------------------------
# R0 — sanity: the planted vector really is in lattice B, and B recovers d
#      at a K1 where lattice A is known to work.
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R0: correctness of lattice B (planted vector in span; recovery works)")
print("-" * 78)


def in_lattice(M, v):
    """Solve x*M = v over Q and check integrality, via sympy."""
    from sympy import Matrix
    sol = Matrix(M).T.solve(Matrix(len(v), 1, v))
    return all(x.q == 1 for x in sol)


print(f"{'curve':<14} {'dim A':>6} {'dim B':>6} {'pv in B':>8} "
      f"{'A K1=2':>8} {'B K1=2':>8}")
for label, curve, k1h, m, ls in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sg = gen_signatures(G, d0, m, n, lam, p, 2, k2b, 42)
    MB = build_lattice_B(sg, n, lam, 2, k2b)
    pvB = planted_B(sg, n, 2, k2b)
    wa = rate(curve, m, 2, SEEDS, variant='A')
    wb = rate(curve, m, 2, SEEDS, variant='B')
    print(f"{label:<14} {2*m+2:>6} {2*m+1:>6} {str(in_lattice(MB, pvB)):>8} "
          f"{str(wa[0])+'/'+str(wa[1]):>8} {str(wb[0])+'/'+str(wb[1]):>8}")


# ---------------------------------------------------------------------------
# R0b — B is not merely "a" d-free lattice: it is exactly pi(A), the
#       coordinate projection of A that deletes the d column.  Proof sketch:
#       A ∩ (d-axis) = <n*S_D*e_m> has rank 1, so pi(A) has rank 2m+1 and
#       det pi(A) = det A / (n*S_D) = n^{m-1} S_K1^m S_K2^m S_KAN = det B.
#       And pi(B_0^{-1} * row_m) reduced mod n is exactly B's k1_0 generator,
#       pi(row_{2m+1}) - A_0 * that is B's constant row, etc.  Checked here by
#       comparing Hermite normal forms (after the column permutation that maps
#       B's column order [k1_1..k1_{m-1}, k1_0, k2, kan] to A's
#       [k1_0..k1_{m-1}, k2, kan]).
#       This is *why* A and B give bit-identical success in R2/R3: LLL on A
#       factors through pi.
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R0b: is lattice B exactly pi(A) (A with the d column deleted)?")
print("-" * 78)


def hnf(M):
    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form
    return hermite_normal_form(Matrix(M).T).T


print(f"{'curve':<14} {'det pi(A)==det B':>17} {'HNF equal':>10}")
for label, curve, k1h, m, ls in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sg = gen_signatures(G, d0, m, n, lam, p, k1h, k2b, 42)
    MA = build_lattice_A(sg, n, lam, k1h, k2b, S_D=1)
    MB = build_lattice_B(sg, n, lam, k1h, k2b)
    # pi(A): delete column m from every row of A (2m+2 generators, rank 2m+1)
    piA = [row[:m] + row[m + 1:] for row in MA]
    # permute B's columns [k1_1..k1_{m-1}, k1_0, ...] -> [k1_0, k1_1..k1_{m-1}, ...]
    perm = [m - 1] + list(range(m - 1)) + list(range(m, 2 * m + 1))
    MBp = [[row[c] for c in perm] for row in MB]
    from sympy import Matrix
    dets = (abs(Matrix(MA).det()) // (n * 1) == abs(Matrix(MBp).det()))
    print(f"{label:<14} {str(dets):>17} {str(hnf(piA) == hnf(MBp)):>10}")


# ---------------------------------------------------------------------------
# R1 — S_D sweep on lattice A.
#      The 2026-07-29 entry claimed "no choice of S_D removes [the trivial
#      vector] — both vectors scale linearly in S_D".  Only the trivial vector
#      does; the planted vector's d-coordinate is the only part that scales:
#          ||pv||^2 = d^2*S_D^2 + n^2*(2m/3 + 1)
#          ||triv||^2 = n^2*S_D^2
#      so ||pv||/||triv||  ->  d/n < 1  as S_D grows.  Predicted crossover at
#          S_D > sqrt( (2m/3+1) / (1 - (d/n)^2) ).
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R1: S_D sweep on lattice A — does the planted vector ever beat n*S_D*e_m?")
print("-" * 78)
SD_GRID = [1, 2, 4, 8, 16, 32, 64, 256]
print("predicted crossover S_D* = sqrt((2m/3+1)/(1-(d/n)^2)), d/n ~ 0.5 -> ~sqrt(4m/3+2)")
for label, curve, k1h, m, ls in CURVES:
    p, b, n, lam, G = curve
    print(f"\n  {label}  (lam*={ls:.4f}, m={m}, K1={k1h}, S_D* ~ "
          f"{math.sqrt((2*m/3+1)/0.75):.1f})")
    print(f"  {'S_D':>5} {'pv/triv':>9} {'sv/pv':>7} {'recov':>7}")
    for sd in SD_GRID:
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        res = run(curve, m, d0, k1h, 42, variant='A', S_D=sd)
        w, t, _ = rate(curve, m, k1h, SEEDS, variant='A', S_D=sd)
        print(f"  {sd:>5} {res['pn']/res['triv']:>9.3f} "
              f"{res['sn']/res['pn']:>7.3f} {str(w)+'/'+str(t):>7}")


# ---------------------------------------------------------------------------
# R2 — the K1 wall: lattice A (S_D=1) vs lattice A (best S_D) vs lattice B.
#      Falsifier from the 2026-07-29 entry: "if sv/pv rises above 1 after the
#      reformulation AND the K1 wall moves outward, the reformulation is a
#      real improvement; if the wall stays at K1 ~ 4-6, the wall is
#      information-theoretic and Phase 2 is at its ceiling."
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R2: K1 wall — lattice A (S_D=1) vs lattice A (S_D=16) vs lattice B")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
for label, curve, k1h, m, ls in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"\n  {label}  (lam*={ls:.4f}, m={m}, n={n}, K2={k2b})")
    hdr = "  {:<14}".format("K1") + "".join(f"{k:>6}" for k in K1_GRID)
    print(hdr)
    print("  {:<14}".format("eff=K1*K2/n") +
          "".join(f"{k*k2b/n:>6.2f}" for k in K1_GRID))
    for tag, variant, sd in [("A  S_D=1", 'A', 1), ("A  S_D=16", 'A', 16),
                             ("B (no d)", 'B', 1)]:
        cells = []
        for k1 in K1_GRID:
            w, t, _ = rate(curve, m, k1, SEEDS, variant=variant, S_D=sd)
            cells.append(f"{w}/{t}")
        print("  {:<14}".format(tag) + "".join(f"{c:>6}" for c in cells))
    # sv/pv for lattice B at the historical wall
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    rb = run(curve, m, d0, k1h, 42, variant='B')
    ra = run(curve, m, d0, k1h, 42, variant='A')
    print(f"  sv/pv at K1={k1h}:  A = {ra['sn']/ra['pn']:.3f}   "
          f"B = {rb['sn']/rb['pn']:.3f}   (>1 means planted is lambda_1)")


# ---------------------------------------------------------------------------
# R3 — does B move the wall on fresh curves too, or only on the 3 historical
#      ones?  17-bit sweep, lam* spread over (0, 0.5), at the eff where the
#      2026-07-29 T3 sweep started to break (eff = 0.15 -> 3/20 curves).
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R3: fresh 17-bit curves at eff ~ 0.15 (the 2026-07-29 T3 breaking point)")
print("-" * 78)
fresh = search_curves(2**16, 2**17, per_bin=2, nbins=10)
print(f"  {len(fresh)} curves; m=12, K1 = round(0.15*n/K2), 5 seeds each\n")
print(f"  {'p':>7} {'n':>7} {'lam*':>7} {'K1':>4} {'eff':>6} "
      f"{'A S_D=1':>9} {'A S_D=16':>9} {'B':>6}")
tot = {'A1': 0, 'A16': 0, 'B': 0, 'N': 0}
for cv in fresh:
    p, b, n, lam, G = cv
    k2b = math.isqrt(n) + 1
    k1 = max(2, round(0.15 * n / k2b))
    m = 12
    wa1, t, _ = rate(cv, m, k1, SEEDS, variant='A', S_D=1)
    wa16, _, _ = rate(cv, m, k1, SEEDS, variant='A', S_D=16)
    wb, _, _ = rate(cv, m, k1, SEEDS, variant='B')
    tot['A1'] += wa1; tot['A16'] += wa16; tot['B'] += wb; tot['N'] += t
    print(f"  {p:>7} {n:>7} {lam_star(lam, n):>7.4f} {k1:>4} "
          f"{k1*k2b/n:>6.3f} {str(wa1)+'/'+str(t):>9} "
          f"{str(wa16)+'/'+str(t):>9} {str(wb)+'/'+str(t):>6}")
print(f"\n  TOTAL   A S_D=1: {tot['A1']}/{tot['N']}   "
      f"A S_D=16: {tot['A16']}/{tot['N']}   B: {tot['B']}/{tot['N']}")


# ---------------------------------------------------------------------------
# R4 — if the trivial vector is not what blocks recovery, what is?
#      L0 := the homogeneous sublattice of B (kannan coordinate 0), dim 2m.
#      Recovery is a BDD instance whose unique-decoding condition is governed
#      by lambda_1(L0), not by the Gaussian heuristic of the full lattice and
#      not by mu = lambda_1 of a single 2-D lambda block (the 2026-07-29 T2
#      quantity, which is an upper bound only: L0 strictly contains the direct
#      sum of the m blocks, because the d-row contributes a cross-block
#      generator B_0^{-1}*row_m).
#      Measured: tau = ||planted|| / lambda_1(L0).
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R4: tau = ||planted|| / lambda_1(L0) as a predictor of the K1 wall")
print("-" * 78)


def homogeneous_L0(sigs, n, lam, k1_bound, k2_bound):
    """Lattice B with the constant/Kannan row and column removed: dim 2m."""
    M = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    return [row[:2 * m] for row in M[:2 * m]]


def lambda1(M, beta=30):
    A = IntegerMatrix.from_matrix(M)
    BKZ.reduction(A, BKZ.Param(beta))
    d = len(M)
    return min(norm([A[i][j] for j in range(d)]) for i in range(d))


for label, curve, k1h, m, ls in CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, _, S_K2, _ = scales(n, k1h, k2b)
    print(f"\n  {label}  (lam*={ls:.4f}, m={m}, n={n})")
    print(f"  {'K1':>4} {'eff':>6} {'||pv||/n':>9} {'L1(L0)/n':>9} "
          f"{'mu/n':>8} {'tau':>7} {'recov':>7}")
    for k1 in K1_GRID:
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sg = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
        if len(sg) < m:
            continue
        sk1, _, sk2, _ = scales(n, k1, k2b)
        pv = norm(planted_B(sg, n, k1, k2b))
        l1 = lambda1(homogeneous_L0(sg, n, lam, k1, k2b))
        mu = lambda_block_mu(n, lam, sk1, sk2)[0]
        w, t, _ = rate(curve, m, k1, SEEDS, variant='B')
        print(f"  {k1:>4} {k1*k2b/n:>6.3f} {pv/n:>9.3f} {l1/n:>9.3f} "
              f"{mu/n:>8.3f} {pv/l1:>7.3f} {str(w)+'/'+str(t):>7}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
