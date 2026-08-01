"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
the target of a BDD instance rather than a hopeless SVP instance.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5):
  In the Phase-2 lattice of `glv_hnp_phase2_20bit.py:263` the shortest vector
  after LLL is ALWAYS the trivial vector n*S_D*e_m (|sv[m]|/n = 1.0000, all
  other coordinates 0).  It is 2-3x shorter than the planted vector on every
  curve, success and failure alike, and it carries no information (d is only
  defined mod n).  So the planted vector is never lambda_1 and no curve-level
  invariant can predict success.

  Proof that n*S_D*e_m is in the lattice: with row_i = n*S_K1*e_i (i<m) and
  row_m = (B_i*S_K1)_i + S_D*e_m,
      n*row_m - sum_i B_i*row_i = n*S_D*e_m.
  Scaling S_D cannot help: both n*S_D*e_m and the planted vector's d-coordinate
  scale linearly in S_D.

The fix tested here: QUOTIENT OUT that direction, i.e. project along e_m.

  L'  = pi(L)          dim 2m+1, drop the d column.
        d is then recovered algebraically from any one signature:
            d = B_0^{-1} (k1_0 + lam*k2_0 - A_0)  mod n.
  L_0 = kan-coordinate also dropped, dim 2m.  The attack becomes an explicit
        CVP: find u in L_0 closest to -t, t = (A_i*S_K1)_i (+) 0_{k2};
        then (k1_i*S_K1, k2_i*S_K2) = u + t.  Solved by Babai nearest-plane
        on an LLL/BKZ-reduced basis (no Kannan embedding at all).

Three attacks are compared on identical signature sets:
  KAN  = original 2m+2 Kannan lattice          (glv_hnp_phase2_20bit.py:263)
  PRJ  = projected 2m+1 Kannan lattice         (this file, build_projected)
  BAB  = Babai nearest-plane CVP on L_0, 2m    (this file, babai_attack)

Experiments:
  E1  correctness / reproduction on the three historical curves
  E2  is the planted vector lambda_1 of L' ?   (sv/pv, and mu = shortest vector
      of the 2D lambda-block, which upper-bounds lambda_1)
  E3  the K1 wall: 2026-07-29 EXP T4 grid re-run under all three attacks.
      FALSIFIER (stated in the 2026-07-29 log): if the wall on the lam*=0.07
      curve moves outward from K1 ~ 4-6, the reformulation is a real
      improvement; if it stays put, the wall is information-theoretic and
      Phase 2 is at its ceiling.
  E4  a parameter-level (not curve-level) predictor:
          gamma = dist / gh(L_0),  gh = sqrt(2m/(2*pi*e)) * covol^(1/2m)
      dist = ||(k1_i*S_K1, k2_i*S_K2)|| ~ n*sqrt(2m/3) is the BDD distance.
      T5 showed no curve-level invariant can work; gamma is not curve-level.
  E5  17-bit sweep: does gamma separate success from failure across curves?

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
# -- verbatim from glv_hnp_phase2_lambda_threshold.py so the comparison to the
#    2026-07-29 numbers is exact.
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

# ---------------------------------------------------------------------------
# Signatures and scales (identical to the 2026-07-29 script)
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

# ---------------------------------------------------------------------------
# ATTACK 1 (baseline) — the original 2m+2 Kannan lattice
# ---------------------------------------------------------------------------

def build_kannan(sigs, n, lam, k1_bound, k2_bound):
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

def planted_kannan(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def attack_kannan(sigs, n, lam, k1_bound, k2_bound, d_secret,
                  use_bkz=False, bkz_beta=20):
    m = len(sigs)
    dim = 2 * m + 2
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(build_kannan(sigs, n, lam, k1_bound, k2_bound))
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = False
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        if (sign * row[m]) % n == d_secret:
            ok = True
            break
    nz = [r for r in rows if any(r)]
    sv = min(norm(r) for r in nz) if nz else float('nan')
    pv = norm(planted_kannan(sigs, d_secret, n, k1_bound, k2_bound))
    return ok, sv, pv

# ---------------------------------------------------------------------------
# ATTACK 2 — projected Kannan lattice L' (2m+1), d recovered algebraically
# ---------------------------------------------------------------------------
#   coords:  0..m-1 = k1 block,  m..2m-1 = k2 block,  2m = Kannan
#   generators (2m+2 of them, rank 2m+1, one relation n*row_d = sum B_i row_i)

def build_projected(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # n*S_K1*e_i
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                            # the (projected) d row
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                       # lambda block
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    r = [0] * dim                            # Kannan / target row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KAN
    rows.append(r)
    return rows

def planted_projected(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def d_from_k(sigs, n, lam, k1_0, k2_0):
    """d = B_0^{-1}(k1_0 + lam*k2_0 - A_0) mod n."""
    B0 = sigs[0]['B'] % n
    if math.gcd(B0, n) != 1:
        return None
    return modinv(B0, n) * ((k1_0 + lam * k2_0 - sigs[0]['A']) % n) % n

def check_d(sigs, n, lam, d, k1_bound, k2_bound):
    """Independent verification: does d make every k_i GLV-small?"""
    for s in sigs:
        k = (s['A'] + s['B'] * d) % n
        # brute-force-free check: we know the true k1,k2, but the attacker
        # would test smallness; here we only need a soundness check.
        if (s['k1'] + lam * s['k2']) % n != k:
            return False
    return True

def attack_projected(sigs, n, lam, k1_bound, k2_bound, d_secret,
                     use_bkz=False, bkz_beta=20):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(build_projected(sigs, n, lam, k1_bound, k2_bound))
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    ok = False
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        k1_0 = sign * row[0]
        k2_0 = sign * row[m]
        if k1_0 % S_K1 or k2_0 % S_K2:
            continue
        d = d_from_k(sigs, n, lam, k1_0 // S_K1, k2_0 // S_K2)
        if d is not None and d == d_secret:
            ok = True
            break
    nz = [r for r in rows if any(r)]
    sv = min(norm(r) for r in nz) if nz else float('nan')
    pv = norm(planted_projected(sigs, n, k1_bound, k2_bound))
    return ok, sv, pv

# ---------------------------------------------------------------------------
# ATTACK 3 — explicit CVP on L_0 (2m) by Babai nearest-plane, no Kannan
# ---------------------------------------------------------------------------

def build_L0(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    return rows

def reduce_basis(rows, dim, use_bkz=False, bkz_beta=20):
    """LLL/BKZ-reduce a (possibly rank-deficient) generating set; return the
    non-zero rows, which form a reduced basis of the lattice they generate."""
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    out = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    return [r for r in out if any(r)]

def gso(basis):
    """Classical Gram-Schmidt in float; returns (bstar, mu)."""
    bstar, mu = [], []
    for v in basis:
        w = [float(x) for x in v]
        coeffs = []
        for j, u in enumerate(bstar):
            nu = sum(x * x for x in u)
            c = (sum(a * b for a, b in zip(v, u)) / nu) if nu > 0 else 0.0
            coeffs.append(c)
            for t in range(len(w)):
                w[t] -= c * u[t]
        bstar.append(w)
        mu.append(coeffs)
    return bstar, mu

def babai_nearest_plane(basis, bstar, target):
    """Return the lattice vector produced by Babai's nearest-plane on target."""
    b = [float(x) for x in target]
    coeff = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        nu = sum(x * x for x in bstar[i])
        if nu == 0:
            continue
        c = round(sum(x * y for x, y in zip(b, bstar[i])) / nu)
        coeff[i] = c
        for t in range(len(b)):
            b[t] -= c * basis[i][t]
    v = [0] * len(target)
    for i, c in enumerate(coeff):
        if c:
            for t in range(len(v)):
                v[t] += c * basis[i][t]
    return v

def attack_babai(sigs, n, lam, k1_bound, k2_bound, d_secret,
                 use_bkz=False, bkz_beta=20):
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    basis = reduce_basis(build_L0(sigs, n, lam, k1_bound, k2_bound), dim,
                         use_bkz, bkz_beta)
    bstar, _ = gso(basis)
    t = [0] * dim
    for i in range(m):
        t[i] = sigs[i]['A'] * S_K1
    neg_t = [-x for x in t]
    u = babai_nearest_plane(basis, bstar, neg_t)
    cand = [u[i] + t[i] for i in range(dim)]      # should equal (k1*S_K1, k2*S_K2)
    if cand[0] % S_K1 or cand[m] % S_K2:
        return False, basis
    d = d_from_k(sigs, n, lam, cand[0] // S_K1, cand[m] // S_K2)
    return (d is not None and d == d_secret), basis

# ---------------------------------------------------------------------------
# Geometry: covolume, Gaussian heuristic, BDD distance, lambda-block mu
# ---------------------------------------------------------------------------

def covolume(basis):
    bstar, _ = gso(basis)
    lg = 0.0
    for w in bstar:
        nw = math.sqrt(sum(x * x for x in w))
        if nw <= 0:
            return 0.0
        lg += math.log(nw)
    return lg                                     # returns log(covol)

def gaussian_heuristic_log(log_covol, rank):
    return 0.5 * math.log(rank / (2 * math.pi * math.e)) + log_covol / rank

def bdd_distance(sigs, n, k1_bound, k2_bound):
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    tot = 0
    for s in sigs:
        tot += (s['k1'] * S_K1) ** 2 + (s['k2'] * S_K2) ** 2
    return math.sqrt(tot)

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

def lambda_block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1])

# ---------------------------------------------------------------------------
# Trial driver
# ---------------------------------------------------------------------------

def trial(curve, m, k1_bound, seed, use_bkz=False, bkz_beta=20, want_geom=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    res = {}
    ok_k, sv_k, pv_k = attack_kannan(sigs, n, lam, k1_bound, k2_bound,
                                     d_secret, use_bkz, bkz_beta)
    ok_p, sv_p, pv_p = attack_projected(sigs, n, lam, k1_bound, k2_bound,
                                        d_secret, use_bkz, bkz_beta)
    ok_b, basis0 = attack_babai(sigs, n, lam, k1_bound, k2_bound,
                                d_secret, use_bkz, bkz_beta)
    res.update(kan=ok_k, prj=ok_p, bab=ok_b,
               svpv_kan=sv_k / pv_k, svpv_prj=sv_p / pv_p)
    if want_geom:
        dist = bdd_distance(sigs, n, k1_bound, k2_bound)
        lc = covolume(basis0)
        gh = math.exp(gaussian_heuristic_log(lc, len(basis0)))
        S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
        res.update(dist=dist, gh0=gh, gamma=dist / gh,
                   mu=lambda_block_mu(n, lam, S_K1, S_K2),
                   pv_prj=pv_p, rank0=len(basis0))
    return res

def rates(curve, m, k1_bound, seeds, use_bkz=False, bkz_beta=20, want_geom=False):
    agg = {'kan': 0, 'prj': 0, 'bab': 0, 'tot': 0}
    geo = []
    svk, svp = [], []
    for s in seeds:
        r = trial(curve, m, k1_bound, s, use_bkz, bkz_beta, want_geom)
        if r is None:
            continue
        agg['tot'] += 1
        for k in ('kan', 'prj', 'bab'):
            agg[k] += bool(r[k])
        svk.append(r['svpv_kan']); svp.append(r['svpv_prj'])
        if want_geom:
            geo.append(r)
    agg['svpv_kan'] = sum(svk) / len(svk) if svk else float('nan')
    agg['svpv_prj'] = sum(svp) / len(svp) if svp else float('nan')
    if geo:
        for k in ('gamma', 'dist', 'gh0', 'mu', 'pv_prj'):
            agg[k] = sum(g[k] for g in geo) / len(geo)
    return agg


# ===========================================================================

print("=" * 82)
print("Thread 23 — projected / BDD reformulation of the GLV-HNP Phase-2 lattice")
print("=" * 82)

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


# ---------------------------------------------------------------------------
print("\n" + "-" * 82)
print("EXP E1: reproduction — do PRJ and BAB recover d wherever KAN does?")
print("-" * 82)
print("KAN = original 2m+2 Kannan lattice (baseline, 2026-07-29 numbers)")
print("PRJ = 2m+1 projected Kannan lattice, d recovered from signature 0")
print("BAB = Babai nearest-plane CVP on the 2m lattice L_0, no Kannan row\n")
print(f"{'curve':<18} {'lam*':>7} {'K1':>4} {'m':>3} {'dim K/P/B':>12} "
      f"{'KAN':>6} {'PRJ':>6} {'BAB':>6}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    a = rates(curve, m, k1, SEEDS)
    print(f"{label:<18} {lam_star(lam,n):>7.4f} {k1:>4} {m:>3} "
      f"{str(2*m+2)+'/'+str(2*m+1)+'/'+str(2*m):>12} "
      f"{str(a['kan'])+'/'+str(a['tot']):>6} "
      f"{str(a['prj'])+'/'+str(a['tot']):>6} "
      f"{str(a['bab'])+'/'+str(a['tot']):>6}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 82)
print("EXP E2: is the planted vector lambda_1 after the projection?")
print("-" * 82)
print("sv/pv = ||shortest reduced row|| / ||planted vector||, in each lattice.")
print("mu    = exact shortest vector of the 2D lambda block <(n S_K1,0),")
print("        (-lam S_K1, S_K2)>; L' contains m copies of it, so it is an")
print("        upper bound on lambda_1(L') independent of the projection.\n")
print(f"{'curve':<18} {'K1':>4} {'sv/pv KAN':>10} {'sv/pv PRJ':>10} "
      f"{'mu/pv_prj':>10} {'gamma':>7}")
for label, curve, k1, m in hist:
    a = rates(curve, m, k1, SEEDS, want_geom=True)
    print(f"{label:<18} {k1:>4} {a['svpv_kan']:>10.3f} {a['svpv_prj']:>10.3f} "
          f"{a['mu']/a['pv_prj']:>10.3f} {a['gamma']:>7.3f}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 82)
print("EXP E3: the K1 wall — 2026-07-29 EXP T4 grid, all three attacks")
print("-" * 82)
print("2026-07-29 T4 found (KAN, m=12): lam*=0.34 curve walls at K1 ~ 12-16,")
print("lam*=0.07 curve walls at K1 ~ 4-6.  FALSIFIER: does the wall move?\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_E3 = 12
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue
    print(f"\n  {label}   lam* = {lam_star(lam,n):.4f}   "
          f"K2 = {math.isqrt(n)+1}   m = {M_E3}")
    print("    " + f"{'K1':>4} {'eff':>6} {'gamma':>7} " +
          f"{'KAN':>6} {'PRJ':>6} {'BAB':>6}")
    for k1 in K1_GRID:
        a = rates(curve, M_E3, k1, SEEDS, want_geom=True)
        eff = k1 * (math.isqrt(n) + 1) / n
        print("    " + f"{k1:>4} {eff:>6.3f} {a['gamma']:>7.3f} "
              f"{str(a['kan'])+'/'+str(a['tot']):>6} "
              f"{str(a['prj'])+'/'+str(a['tot']):>6} "
              f"{str(a['bab'])+'/'+str(a['tot']):>6}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 82)
print("EXP E4: does BKZ-20 change the picture on the projected lattice?")
print("-" * 82)
fail_curve = [c for lbl, c, _, _ in hist if c[2] == 2647][0]
print(f"{'K1':>4} {'gamma':>7} {'KAN(LLL)':>9} {'PRJ(LLL)':>9} "
      f"{'KAN(BKZ)':>9} {'PRJ(BKZ)':>9} {'BAB(BKZ)':>9}")
for k1 in (6, 8, 12, 16):
    a = rates(fail_curve, M_E3, k1, SEEDS, want_geom=True)
    bz = rates(fail_curve, M_E3, k1, SEEDS, use_bkz=True, bkz_beta=20)
    print(f"{k1:>4} {a['gamma']:>7.3f} "
          f"{str(a['kan'])+'/'+str(a['tot']):>9} "
          f"{str(a['prj'])+'/'+str(a['tot']):>9} "
          f"{str(bz['kan'])+'/'+str(bz['tot']):>9} "
          f"{str(bz['prj'])+'/'+str(bz['tot']):>9} "
          f"{str(bz['bab'])+'/'+str(bz['tot']):>9}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 82)
print("EXP E5: 17-bit cross-curve sweep — does gamma separate success/failure?")
print("-" * 82)

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

def search_curves(lo, hi, per_bin=1, nbins=10):
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
                    roots = glv_roots(n_c)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_c) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_c)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_c, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

LO, HI = 2**16, 2**17
curves = search_curves(LO, HI, per_bin=1, nbins=10)
print(f"Found {len(curves)} 17-bit j=0 GLV curves (one per lam* decile).\n")

M_E5 = 12
allrows = []
for eff_t in (0.05, 0.10, 0.15, 0.20):
    print(f"\n  === eff target {eff_t} (m={M_E5}, {len(SEEDS)} seeds) ===")
    print("    " + f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} "
          f"{'gamma':>7} {'mu/pv':>7} {'KAN':>6} {'PRJ':>6} {'BAB':>6}")
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        a = rates((p, b, n, lam, G), M_E5, k1, SEEDS, want_geom=True)
        row = {'p': p, 'n': n, 'lam_star': lam_star(lam, n), 'k1': k1,
               'eff': k1 * k2 / n, 'gamma': a['gamma'],
               'mupv': a['mu'] / a['pv_prj'],
               'kan': a['kan'] == a['tot'], 'prj': a['prj'] == a['tot'],
               'bab': a['bab'] == a['tot'],
               'kan_any': a['kan'] > 0, 'prj_any': a['prj'] > 0}
        allrows.append(row)
        print("    " + f"{p:>7} {n:>7} {row['lam_star']:>6.3f} {k1:>4} "
              f"{row['eff']:>6.3f} {a['gamma']:>7.3f} {row['mupv']:>7.3f} "
              f"{str(a['kan'])+'/'+str(a['tot']):>6} "
              f"{str(a['prj'])+'/'+str(a['tot']):>6} "
              f"{str(a['bab'])+'/'+str(a['tot']):>6}")

print("\n" + "-" * 82)
print("Predictor evaluation on the pooled E5 grid "
      f"({len(allrows)} curve x eff cells)")
print("-" * 82)

def best_threshold(rows, key, label, direction=-1):
    """direction=-1: predict success iff value <= thr."""
    vals = sorted({r[key] for r in rows})
    best = (0, None)
    for thr in vals:
        acc = sum(1 for r in rows if ((r[key] <= thr) if direction < 0
                                      else (r[key] >= thr)) == r[label])
        if acc > best[0]:
            best = (acc, thr)
    return best

for label in ('kan', 'prj', 'bab'):
    npos = sum(1 for r in allrows if r[label])
    ntot = len(allrows)
    if npos in (0, ntot):
        print(f"[{label}] degenerate {npos}/{ntot}")
        continue
    base = max(npos, ntot - npos)
    print(f"\n[{label}] positives {npos}/{ntot}, majority baseline "
          f"{base}/{ntot} = {base/ntot:.1%}")
    for key, dr in (('gamma', -1), ('lam_star', +1), ('eff', -1), ('mupv', +1)):
        acc, thr = best_threshold(allrows, key, label, dr)
        op = '<=' if dr < 0 else '>='
        print(f"    {key:<9} {op} {thr:>8.4f} -> {acc}/{ntot} = {acc/ntot:.1%}")
    s = [r['gamma'] for r in allrows if r[label]]
    f = [r['gamma'] for r in allrows if not r[label]]
    if s and f:
        print(f"    gamma range: success [{min(s):.3f}, {max(s):.3f}]  "
              f"failure [{min(f):.3f}, {max(f):.3f}]  "
              f"overlap={min(s) <= max(f) and min(f) <= max(s)}")

print("\nDone.")
