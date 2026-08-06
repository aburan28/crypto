"""
GLV-HNP Phase 2, Thread 23: d-eliminated lattice — does killing the trivial
vector n*S_D*e_m move the K1 wall?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, experiment T5):
  The Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` (build_glv_lattice)
  always has shortest vector = n*S_D*e_m, the "trivial" vector living purely
  in the d column.  It is 1.6-3x shorter than the planted vector on every
  curve tested, success and failure alike, and no choice of S_D removes it
  (both vectors scale linearly in S_D).  The planted vector is therefore
  never lambda_1, and recovery is a BDD/coset condition rather than SVP.

  That run proposed Thread 23: reformulate so the target IS lambda_1, then
  check whether the K1 wall measured in T4 (K1 ~ 12-16 for lam*=0.34,
  K1 ~ 4-6 for lam*=0.07) moves outward.

This script implements the reformulation by ALGEBRAIC ELIMINATION of d.

  ECDSA/GLV relations:   A_i + B_i*d = k1_i + lam*k2_i  (mod n),  i = 0..m-1
  with 0 <= k1_i < K1, 0 <= k2_i < K2.

  Take i=0 as reference and set
        t_i = B_i * B_0^{-1} mod n,      u_i = (A_i - t_i*A_0) mod n.
  Substituting d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} into equation i gives,
  for i = 1..m-1,
        k1_i  ==  u_i + t_i*k1_0 + t_i*lam*k2_0 - lam*k2_i   (mod n).      (*)

  Every unknown in (*) is SMALL (k1_0 < K1, k2_j < K2, and k1_i is the
  residual).  d, the only full-size unknown, is gone.  The lattice built from
  (*) therefore has no full-size coordinate and no trivial vector: the
  competing short vectors are only the lam-block/continued-fraction ones.

  Dimension bookkeeping:  original lattice 2m+2, d-eliminated 2m+1 (exactly
  one coordinate removed).  Determinants (both bases are lower-triangular
  w.r.t. column index, so det = product of diagonals):
        det_orig  = (n*S_K1)^m     * S_D * S_K2^m * S_KANNAN
        det_delim = (n*S_K1)^(m-1) * S_K1 * S_K2^m * S_KANNAN
  i.e. det_delim = det_orig / (n * S_D / S_K1); with S_D=1 the volume drops
  by n*K1/n = ... see gaussian_heuristic() below for the exact numbers.
  Information content is identical, so any gain must come from geometry.

Experiments:
  E1  correctness of the d-eliminated construction (recovers the same d)
  E2  is the planted vector lambda_1 now?  sv/pv for both lattices
  E3  THE FALSIFIER: K1-wall grid on the two historical 12-bit curves,
      both lattices, 5 seeds
  E4  m-dependence at the wall (replicates T4b in the new lattice)
  E5  17-bit head-to-head at eff = 0.15 (the discriminating operating point:
      the original lattice scored 3/20 curves there on 2026-07-29)
  E6  does ||v_planted|| / GH(L) predict the wall?

Run: python3 glv_hnp_phase2_delim.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py:87 so comparisons are exact
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
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    """j=0 GLV curves with p in [lo,hi), n prime = 1 mod 3, bucketed by lam*."""
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
# Signatures (verbatim scaling convention from glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN)."""
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

def det_orig(m, n, S_K1, S_D, S_K2, S_KAN):
    return (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN

# ---------------------------------------------------------------------------
# D-ELIMINATED lattice (dim 2m+1)  -- Thread 23
# ---------------------------------------------------------------------------
#
# Column layout (dim = 2m+1):
#   0 .. m-2   equation columns for signatures i = 1..m-1  (residual k1_i*S_K1)
#   m-1        k1_0 column                                  (k1_0 * S_K1)
#   m .. 2m-1  k2_i columns, i = 0..m-1                     (k2_i * S_K2)
#   2m         Kannan column                                (S_KANNAN)
#
# Every row's rightmost nonzero entry sits on its own index, so the basis is
# triangular and det = product of diagonals.

def delim_params(sigs, n):
    """(t_i, u_i) for i=1..m-1, eliminating d via signature 0."""
    B0inv = modinv(sigs[0]['B'], n)
    A0 = sigs[0]['A']
    t, u = [], []
    for i in range(1, len(sigs)):
        ti = sigs[i]['B'] * B0inv % n
        ui = (sigs[i]['A'] - ti * A0) % n
        t.append(ti)
        u.append(ui)
    return t, u

def build_delim_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    t, u = delim_params(sigs, n)
    lam_r = lam % n
    M = [[0] * dim for _ in range(dim)]
    # mod-n rows on the m-1 equation columns
    for j in range(m - 1):
        M[j][j] = n * S_K1
    # k1_0 row
    for j in range(m - 1):
        M[m - 1][j] = t[j] * S_K1
    M[m - 1][m - 1] = S_K1
    # k2_0 row (coefficient t_i*lam on equation column j = i-1)
    for j in range(m - 1):
        M[m][j] = t[j] * lam_r % n * S_K1
    M[m][m] = S_K2
    # k2_i rows, i = 1..m-1  (coefficient -lam on equation column i-1)
    for i in range(1, m):
        M[m + i][i - 1] = -lam_r * S_K1
        M[m + i][m + i] = S_K2
    # Kannan row
    for j in range(m - 1):
        M[2 * m][j] = u[j] * S_K1
    M[2 * m][2 * m] = S_KANNAN
    return M

def planted_vector_delim(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for j in range(m - 1):
        v[j] = sigs[j + 1]['k1'] * S_K1
    v[m - 1] = sigs[0]['k1'] * S_K1
    for i in range(m):
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_delim(M_reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read k1_0,k2_0 off any row with Kannan coordinate +-S_KANNAN, then
    d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sg = 1 if last > 0 else -1
        a = sg * row[m - 1]
        b = sg * row[m]
        if a % S_K1 or b % S_K2:
            continue
        k1_0, k2_0 = a // S_K1, b // S_K2
        k0 = (k1_0 + lam * k2_0) % n
        d_cand = (k0 - sigs[0]['A']) * B0inv % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def det_delim(m, n, S_K1, S_K2, S_KAN):
    return (n * S_K1) ** (m - 1) * S_K1 * S_K2 ** m * S_KAN

# ---------------------------------------------------------------------------
# Common measurement helpers
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def gaussian_heuristic(det, dim):
    """GH(L) = sqrt(dim / (2*pi*e)) * det^(1/dim), computed in logs."""
    return math.exp(math.log(det) / dim) * math.sqrt(dim / (2 * math.pi * math.e))

def run_one(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False, bkz_beta=20):
    """variant in {'orig','delim'}.  Returns dict of measurements or None."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if variant == 'orig':
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 2
        pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
        det = det_orig(m, n, S_K1, S_D, S_K2, S_KAN)
    else:
        M = build_delim_lattice(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1
        pv = planted_vector_delim(sigs, n, k1_bound, k2_bound)
        det = det_delim(m, n, S_K1, S_K2, S_KAN)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    if variant == 'orig':
        ok = recover_d(reduced, m, n, S_KAN, d_secret) is not None
    else:
        ok = recover_d_delim(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
    pn = norm(pv)
    norms = [norm(r) for r in reduced]
    sn = min(nn for nn in norms if nn > 0)
    return {'ok': ok, 'pn': pn, 'sn': sn, 'gh': gaussian_heuristic(det, dim),
            'dim': dim, 'sigs': sigs, 'reduced': reduced, 'pv': pv}

def success_rate(curve, m, k1_bound, seeds, variant, use_bkz=False, bkz_beta=20):
    wins, tot, ratios, ghr = 0, 0, [], []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        r = run_one(curve, m, d_trial, k1_bound, seed, variant, use_bkz, bkz_beta)
        if r is None:
            continue
        tot += 1
        wins += bool(r['ok'])
        ratios.append(r['sn'] / r['pn'])
        ghr.append(r['pn'] / r['gh'])
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    mg = sum(ghr) / len(ghr) if ghr else float('nan')
    return wins, tot, mr, mg

def auc(rows, key, label='ok'):
    """AUC of `key` as a score for predicting NOT-label (lower key = success)."""
    pos = [r[key] for r in rows if r[label]]
    neg = [r[key] for r in rows if not r[label]]
    if not pos or not neg:
        return float('nan')
    c = sum((1.0 if a < b else 0.5 if a == b else 0.0) for a in pos for b in neg)
    return c / (len(pos) * len(neg))

# ===========================================================================

print("=" * 78)
print("Thread 23 — d-eliminated GLV-HNP Phase-2 lattice")
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

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1: correctness of the d-eliminated construction")
print("-" * 78)
print("Checks (a) the planted vector really is in the lattice (integer coords in")
print("the basis), and (b) recovery returns the true d.  K1 kept small so both")
print("lattices are expected to succeed.\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} {'dim_o':>6} {'dim_d':>6} "
      f"{'orig':>6} {'delim':>6} {'pv in L':>8}")

def pv_in_lattice(M, pv):
    """Solve x*M = pv over Q and check x is integral."""
    sol = sympy.Matrix(M).T.solve(sympy.Matrix(pv))
    return all(sympy.Rational(s).q == 1 for s in sol)

for label, curve, k1, m in hist:
    n = curve[2]
    k1_small = 2
    d_t = random.Random(42 + 7777).randint(1, n - 1)
    ro = run_one(curve, m, d_t, k1_small, 42, 'orig')
    rd = run_one(curve, m, d_t, k1_small, 42, 'delim')
    k2b = math.isqrt(n) + 1
    Md = build_delim_lattice(rd['sigs'], n, curve[3], k1_small, k2b)
    pvd = planted_vector_delim(rd['sigs'], n, k1_small, k2b)
    inl = pv_in_lattice(Md, pvd)
    print(f"{label:<18} {k1_small:>3} {m:>3} {ro['dim']:>6} {rd['dim']:>6} "
          f"{str(ro['ok']):>6} {str(rd['ok']):>6} {str(inl):>8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E2: is the planted vector lambda_1 now?   sv/pv and pv/GH")
print("-" * 78)
print("sv/pv < 1  =>  something shorter than the planted vector exists.")
print("T5 (2026-07-29) measured sv/pv in [0.34, 0.61] for 'orig'.\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} | {'sv/pv orig':>10} {'pv/GH orig':>10} "
      f"| {'sv/pv del':>10} {'pv/GH del':>10} | {'o':>3} {'d':>3}")

e2rows = []
for label, curve, k1, m in hist:
    wo, to, mro, mgo = success_rate(curve, m, k1, SEEDS, 'orig')
    wd, td, mrd, mgd = success_rate(curve, m, k1, SEEDS, 'delim')
    print(f"{label:<18} {k1:>3} {m:>3} | {mro:>10.3f} {mgo:>10.3f} "
          f"| {mrd:>10.3f} {mgd:>10.3f} | {wo:>1}/{to:<1} {wd:>1}/{td:<1}")
    e2rows.append((label, mro, mrd))

# where does the shortest vector of the delim lattice live?
print("\nEnergy profile of the shortest reduced vector (delim lattice), "
      "fraction of ||.||^2:")
print(f"{'curve':<18} {'eq-cols':>8} {'k1_0':>8} {'k2-cols':>8} {'kannan':>8} "
      f"{'== pv?':>7}")
for label, curve, k1, m in hist:
    n = curve[2]
    d_t = random.Random(42 + 7777).randint(1, n - 1)
    r = run_one(curve, m, d_t, k1, 42, 'delim')
    red = r['reduced']
    sv = min((row for row in red if norm(row) > 0), key=norm)
    tot = sum(x * x for x in sv)
    eq = sum(sv[j] ** 2 for j in range(m - 1)) / tot
    a1 = sv[m - 1] ** 2 / tot
    k2 = sum(sv[m + i] ** 2 for i in range(m)) / tot
    kan = sv[2 * m] ** 2 / tot
    same = abs(norm(sv) - r['pn']) < 1e-6
    print(f"{label:<18} {eq:>8.3f} {a1:>8.3f} {k2:>8.3f} {kan:>8.3f} {str(same):>7}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E3 (FALSIFIER): K1 wall, original vs d-eliminated, 5 seeds")
print("-" * 78)
print("2026-07-29 T4 walls: 12-bit/2557 -> K1 ~ 12-16;  12-bit/2677 -> K1 ~ 4-6.")
print("If d-elimination is a real improvement the delim row moves right.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
e3rows = []
for label, curve, _k1, m in hist[1:]:
    n = curve[2]
    print(f"{label}  (n={n}, lam*={lam_star(curve[3], n):.3f}, m={m})")
    hdr = "  " + " ".join(f"{k:>5}" for k in K1_GRID)
    print(f"  {'K1':<8}" + hdr)
    for variant in ('orig', 'delim'):
        cells = []
        for k1 in K1_GRID:
            w, t, mr, mg = success_rate(curve, m, k1, SEEDS, variant)
            cells.append(f"{w}/{t}")
            e3rows.append({'curve': label, 'variant': variant, 'K1': k1,
                           'wins': w, 'tot': t, 'svpv': mr, 'ghr': mg,
                           'ok': w >= 3})
        print(f"  {variant:<8}  " + " ".join(f"{c:>5}" for c in cells))
    # BKZ-40 retry on the delim lattice at the wall
    cells = []
    for k1 in K1_GRID:
        w, t, mr, mg = success_rate(curve, m, k1, SEEDS, 'delim',
                                    use_bkz=True, bkz_beta=40)
        cells.append(f"{w}/{t}")
    print(f"  {'del+BKZ40':<8}  " + " ".join(f"{c:>5}" for c in cells))
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E4: more data at the wall (replicates T4b for the delim lattice)")
print("-" * 78)
label, curve, _k1, _m = hist[2]
n = curve[2]
print(f"{label}  n={n}, K1=8, varying m\n")
M_GRID = [8, 12, 16, 24, 32]
print(f"  {'m':<10}" + " ".join(f"{mm:>5}" for mm in M_GRID))
for variant in ('orig', 'delim'):
    cells = []
    for mm in M_GRID:
        w, t, mr, mg = success_rate(curve, mm, 8, SEEDS, variant)
        cells.append(f"{w}/{t}")
    print(f"  {variant:<10}" + " ".join(f"{c:>5}" for c in cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E5: 17-bit head-to-head at eff = K1*K2/n = 0.15")
print("-" * 78)
print("2026-07-29 T3 measured 3/20 curves recovering 5/5 for 'orig' at this")
print("operating point -- the most discriminating one available.\n")

curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
print(f"collected {len(curves17)} curves in [2^16, 2^17)\n")

M17 = 12
EFF = 0.15
print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} | {'orig':>6} {'delim':>6} "
      f"| {'sv/pv o':>8} {'sv/pv d':>8} {'pv/GH o':>8} {'pv/GH d':>8}")
tot_o = tot_d = 0
e5rows = []
for (p, b, n, lam, G) in curves17:
    k2b = math.isqrt(n) + 1
    k1b = max(2, int(EFF * n / k2b))
    cur = (p, b, n, lam, G)
    wo, to, mro, mgo = success_rate(cur, M17, k1b, SEEDS, 'orig')
    wd, td, mrd, mgd = success_rate(cur, M17, k1b, SEEDS, 'delim')
    tot_o += (wo == len(SEEDS))
    tot_d += (wd == len(SEEDS))
    e5rows.append({'p': p, 'n': n, 'ls': lam_star(lam, n), 'K1': k1b,
                   'wo': wo, 'wd': wd, 'ghr_o': mgo, 'ghr_d': mgd,
                   'ok': wd >= 3})
    print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1b:>4} | "
          f"{wo:>1}/{to:<4} {wd:>1}/{td:<4} | {mro:>8.3f} {mrd:>8.3f} "
          f"{mgo:>8.3f} {mgd:>8.3f}")

nc = len(curves17)
print(f"\nfull recovery (5/5): orig {tot_o}/{nc}   delim {tot_d}/{nc}")
sum_o = sum(r['wo'] for r in e5rows)
sum_d = sum(r['wd'] for r in e5rows)
print(f"total seed-wins:     orig {sum_o}/{nc * len(SEEDS)}   "
      f"delim {sum_d}/{nc * len(SEEDS)}")
better = sum(1 for r in e5rows if r['wd'] > r['wo'])
worse = sum(1 for r in e5rows if r['wd'] < r['wo'])
print(f"per-curve head-to-head: delim better on {better}, worse on {worse}, "
      f"tied on {nc - better - worse}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E6: does ||v_planted|| / GH(L) predict the wall?")
print("-" * 78)
print("Model: recovery iff the planted vector is below the Gaussian heuristic.")
print("Threshold-free check = AUC (1.0 = perfect, 0.5 = coin flip).\n")

for variant in ('orig', 'delim'):
    rows = [r for r in e3rows if r['variant'] == variant]
    a = auc(rows, 'ghr')
    lo = [r['ghr'] for r in rows if r['ok']]
    hi = [r['ghr'] for r in rows if not r['ok']]
    print(f"E3 grid, {variant:<6}: AUC(pv/GH) = {a:.3f}   "
          f"success pv/GH in [{min(lo):.2f},{max(lo):.2f}]   "
          f"failure in [{min(hi):.2f},{max(hi):.2f}]")

key = 'ghr_d'
a5 = auc(e5rows, key)
lo = [r[key] for r in e5rows if r['ok']]
hi = [r[key] for r in e5rows if not r['ok']]
if lo and hi:
    print(f"E5 17-bit, delim : AUC(pv/GH) = {a5:.3f}   "
          f"success in [{min(lo):.2f},{max(lo):.2f}]   "
          f"failure in [{min(hi):.2f},{max(hi):.2f}]")
else:
    print(f"E5 17-bit, delim : degenerate (all {'success' if lo else 'failure'})")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E7: is d-elimination EXACTLY the projection along the trivial vector?")
print("-" * 78)
print("ker(pi) & L_orig = Z*(n*S_D*e_m) is rank 1 (row m must be used n times to")
print("clear the B_i*S_K1 entries), so pi(L_orig) has dim 2m+1 like L_delim.")
print("Both live in the same coordinate space -- L_orig's equation column i and")
print("L_delim's (k1_0 col for i=0, equation col i-1 otherwise) both carry")
print("k1_i*S_K1 -- so the lattices are directly comparable.  Test: HNF equality.\n")

from sympy.matrices.normalforms import hermite_normal_form

def hnf(rows):
    return hermite_normal_form(sympy.Matrix(rows).T).T

print(f"{'curve':<18} {'m':>3} {'dim':>4} {'HNF(pi L_orig) == HNF(L_delim)':>32} "
      f"{'pv equal':>9}")
for label, curve, k1, m in hist:
    n = curve[2]
    k2b = math.isqrt(n) + 1
    d_t = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(curve[4], d_t, m, n, curve[3], curve[0], k1, k2b, 42)
    Mo = build_glv_lattice(sigs, n, curve[3], k1, k2b)
    Md = build_delim_lattice(sigs, n, curve[3], k1, k2b)
    # pi: delete the d column (index m) from every row of L_orig.
    # resulting column order: [k1_0, k1_1..k1_{m-1}, k2_0..k2_{m-1}, kannan]
    proj = [row[:m] + row[m + 1:] for row in Mo]
    # permute L_delim columns [eq(1..m-1), k1_0, k2..., kannan] into the same order
    perm = [m - 1] + list(range(m - 1)) + list(range(m, 2 * m + 1))
    Mdp = [[row[c] for c in perm] for row in Md]
    same = hnf(proj) == hnf(Mdp)
    pvo = planted_vector(sigs, d_t, n, k1, k2b)
    pvd = planted_vector_delim(sigs, n, k1, k2b)
    pvo_p = pvo[:m] + pvo[m + 1:]
    pvd_p = [pvd[c] for c in perm]
    print(f"{label:<18} {m:>3} {2*m+1:>4} {str(same):>32} "
          f"{str(pvo_p == pvd_p):>9}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E8: the projection makes lambda_1 INFORMATIVE (instance-level)")
print("-" * 78)
print("In L_orig, lambda_1 is always the trivial vector n*S_D*e_m, so sv/pv is")
print("pinned near a constant (E5: 0.347-0.368 across all 20 curves) and carries")
print("zero information -- which retro-explains the 2026-06-21..06-29 run of failed")
print("separators.  In pi(L_orig) the trivial vector is gone and lambda_1 becomes")
print("the true rival.  Compare three predictors at the INSTANCE level")
print("(one row per (curve, seed), 2 operating points, 5 seeds):")
print("  sv/E[pv]  : shortest reduced vector / expected planted norm  (post-reduction,")
print("              secret-free -- E[pv] needs only n, K1, K2, m)")
print("  nu_hat    : lambda_1(L2)/sqrt(det L2) of the 2-D lam block  (2026-07-29,")
print("              a-priori curve invariant, reported AUC 0.935)")
print("  pv/GH     : planted norm over Gaussian heuristic\n")

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

def nu_hat(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] ** 2 + w[1] ** 2) / math.sqrt(n * S_K1 * S_K2)

def expected_pv_norm(m, n, K1, K2, S_K1, S_K2, S_KAN):
    """E||v_planted|| for the projected/delim lattice: m k1-coords, m k2-coords,
    one Kannan coord.  Uses only public parameters."""
    return math.sqrt(m * (K1 * S_K1) ** 2 / 3.0
                     + m * (K2 * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)

inst = []
for eff in (0.15, 0.20):
    for (p, b, n, lam, G) in curves17:
        k2b = math.isqrt(n) + 1
        k1b = max(2, int(eff * n / k2b))
        S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
        nh = nu_hat(n, lam, S_K1, S_K2)
        epv = expected_pv_norm(M17, n, k1b, k2b, S_K1, S_K2, S_KAN)
        for seed in SEEDS:
            d_t = random.Random(seed + 7777).randint(1, n - 1)
            r = run_one((p, b, n, lam, G), M17, d_t, k1b, seed, 'delim')
            if r is None:
                continue
            inst.append({'ok': r['ok'], 'sv_epv': r['sn'] / epv,
                         'nuhat': nh, 'ghr': r['pn'] / r['gh'],
                         'sv_pv': r['sn'] / r['pn'], 'eff': eff})

nok = sum(1 for r in inst if r['ok'])
print(f"{len(inst)} instances, {nok} successes, {len(inst) - nok} failures\n")
print(f"{'predictor':<12} {'AUC':>7}   {'success range':>22} {'failure range':>22}")
for key in ('sv_epv', 'nuhat', 'ghr', 'sv_pv'):
    a = auc(inst, key)
    lo = [r[key] for r in inst if r['ok']]
    hi = [r[key] for r in inst if not r['ok']]
    print(f"{key:<12} {a:>7.3f}   [{min(lo):>9.3f},{max(lo):>9.3f}] "
          f"[{min(hi):>9.3f},{max(hi):>9.3f}]")

# the same predictors measured on the ORIGINAL lattice, for contrast
inst_o = []
for eff in (0.15, 0.20):
    for (p, b, n, lam, G) in curves17:
        k2b = math.isqrt(n) + 1
        k1b = max(2, int(eff * n / k2b))
        S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
        epv = expected_pv_norm(M17, n, k1b, k2b, S_K1, S_K2, S_KAN)
        for seed in SEEDS:
            d_t = random.Random(seed + 7777).randint(1, n - 1)
            r = run_one((p, b, n, lam, G), M17, d_t, k1b, seed, 'orig')
            if r is None:
                continue
            inst_o.append({'ok': r['ok'], 'sv_epv': r['sn'] / epv})
a_o = auc(inst_o, 'sv_epv')
lo = [r['sv_epv'] for r in inst_o if r['ok']]
hi = [r['sv_epv'] for r in inst_o if not r['ok']]
print(f"\nsame statistic on L_orig (lambda_1 = trivial vector):")
print(f"{'sv_epv orig':<12} {a_o:>7.3f}   [{min(lo):>9.3f},{max(lo):>9.3f}] "
      f"[{min(hi):>9.3f},{max(hi):>9.3f}]")

print("\nCONFOUND CHECK.  sv_epv on L_orig has total range [0.333,0.334] -- the")
print("shortest vector there is exactly n*S_D, so the statistic is constant per")
print("operating point and its 'AUC' can only be separating the two eff arms,")
print("which have different base success rates.  Recompute every AUC WITHIN each")
print("arm, where eff is held fixed:\n")
print(f"{'arm':<10} {'n_ok':>5} {'n_fail':>7} " +
      " ".join(f"{k:>9}" for k in ('sv_epv', 'nuhat', 'ghr', 'sv_pv')))
for eff in (0.15, 0.20):
    sub = [r for r in inst if r['eff'] == eff]
    k_ok = sum(1 for r in sub if r['ok'])
    cells = [auc(sub, k) for k in ('sv_epv', 'nuhat', 'ghr', 'sv_pv')]
    print(f"eff={eff:<6.2f} {k_ok:>5} {len(sub) - k_ok:>7} " +
          " ".join(f"{c:>9.3f}" for c in cells))
sub_o = [r for r in inst_o]
print(f"\n(L_orig sv_epv is constant within an arm, so its within-arm AUC is "
      f"undefined/0.5 by construction.)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
