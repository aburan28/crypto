"""
GLV-HNP Phase 2, Thread 23: does removing the trivial lattice vector move the wall?

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, Thread 20 / T5):
the Phase-2 Kannan lattice of `glv_hnp_phase2_20bit.py:263` always contains the
vector  n*S_D*e_m  (take n*row_m minus sum_i B_i*row_i).  Its norm is n*S_D while
the planted vector has norm ~ n*S_D*sqrt(2m/3 + 4/3), so the planted vector is
NEVER lambda_1 and no choice of S_D fixes it (both scale linearly in S_D).
Thread 23's proposal: reformulate so the target IS lambda_1, then check whether
the K1 wall of T4 moves outward.

This script implements the reformulation by ELIMINATING d.  Using signature 0,
    d = B_0^{-1} (k_full,0 - A_0)  (mod n),
every other signature becomes
    k_full,i = C_i + beta_i * k_full,0   (mod n),
    beta_i = B_i B_0^{-1},  C_i = A_i - beta_i A_0,   (beta_0 = 1, C_0 = 0).
The d column disappears, so the trivial vector disappears with it.

Experiments
  E0  correctness: planted vector is in the eliminated lattice; recovery works.
  E1  geometry: where does the planted vector sit in each lattice (sv/pv,
      Gaussian-heuristic ratio, is pv = lambda_1)?
  E2  the T4 grid (12-bit/2557 and 12-bit/2677, K1 in {2,...,24}, 5 seeds),
      baseline vs eliminated.  This is Thread 23's stated falsifier.
  E3  eff = K1*K2/n sweep on fresh 17-bit curves, baseline vs eliminated,
      against the Gaussian-heuristic prediction and the information bound.
  E4  does extra data help the eliminated lattice where T4b showed it does not
      help the baseline (m sweep at a fixed failing K1)?

Run: python3 glv_hnp_phase2_thread23_elim.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a = 0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py:87 so results are directly comparable.
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
# CM theory for j = 0 curves
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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, or None."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (r1, r2)

def lam_star(lam, n):
    return min(lam, n - lam) / n

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
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, want, seed=12345):
    """Collect `want` j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    cur = build_curve(p, n_cand, seed=seed)
                    if cur is None:
                        continue
                    out.append((p, cur[1], n_cand, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Signature generation (verbatim shape from glv_hnp_phase2_20bit.py:235)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# BASELINE lattice -- verbatim from glv_hnp_phase2_20bit.py:263
# columns: [0..m-1] k1 residual (S_K1) | m: d (S_D) | [m+1..2m] k2 (S_K2)
#          | 2m+1: Kannan (S_KANNAN);   dim = 2m+2
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

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
    return M, S_K1, S_D, S_K2, S_KANNAN

def baseline_planted(sigs, d_secret, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    v = [0] * dim
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[dim - 1] = n
    return v

def recover_d_baseline(rows, m, n, S_KANNAN, d_secret):
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

# ---------------------------------------------------------------------------
# ELIMINATED lattice (Thread 23)
# unknowns: u = k_full,0 ; k2_0..k2_{m-1} ; t_1..t_{m-1}
#   k1_i = C_i + beta_i*u - lam*k2_i - n*t_i        (t_0 folded into the
#                                                    n*e_0 vector, see below)
# columns: [0..m-1] k1 (w1) | [m..2m-1] k2 (w2) | 2m: Kannan (W); dim = 2m+1
# rows: R_u, R_t(i>=1) [m-1], R_k2(i) [m], R_c  -> 2m+1 rows, full rank
#   (n*e_0 = n*R_u + sum_{i>=1} beta_i*R_t_i, so dropping R_t_0 loses nothing)
# ---------------------------------------------------------------------------

def elim_params(sigs, n, k1_bound, k2_bound):
    B0inv = modinv(sigs[0]['B'], n)
    beta = [sigs[i]['B'] * B0inv % n for i in range(len(sigs))]
    C = [(sigs[i]['A'] - beta[i] * sigs[0]['A']) % n for i in range(len(sigs))]
    W = k1_bound * k2_bound
    w1 = W // k1_bound          # = k2_bound
    w2 = max(1, W // k2_bound)  # = k1_bound
    return beta, C, W, w1, w2

def build_elim_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    beta, C, W, w1, w2 = elim_params(sigs, n, k1_bound, k2_bound)

    M = [[0] * dim for _ in range(dim)]
    r = 0
    # R_u
    for i in range(m):
        M[r][i] = beta[i] * w1
    r += 1
    # R_t_i, i = 1..m-1
    for i in range(1, m):
        M[r][i] = -n * w1
        r += 1
    # R_k2_i, i = 0..m-1
    for i in range(m):
        M[r][i] = -lam * w1
        M[r][m + i] = w2
        r += 1
    # R_c
    for i in range(m):
        M[r][i] = C[i] * w1
    M[r][2 * m] = W
    r += 1
    assert r == dim
    return M, W, w1, w2

def elim_planted(sigs, n, lam, k1_bound, k2_bound):
    """The planted vector, built from its known coordinates (and separately
    verified to lie in the row span by elim_planted_coeffs)."""
    m = len(sigs)
    _, _, W, w1, w2 = elim_params(sigs, n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * w1
        v[m + i] = sigs[i]['k2'] * w2
    v[2 * m] = W
    return v

def elim_planted_in_lattice(M, sigs, n, lam, k1_bound, k2_bound):
    """Explicitly exhibit the integer coefficient vector, proving membership."""
    m = len(sigs)
    beta, C, W, w1, w2 = elim_params(sigs, n, k1_bound, k2_bound)
    u = sigs[0]['k_full']
    coeffs = [0] * (2 * m + 1)
    coeffs[0] = u                                   # R_u
    for i in range(1, m):                           # R_t_i
        t_i = (C[i] + beta[i] * u - lam * sigs[i]['k2'] - sigs[i]['k1'])
        assert t_i % n == 0
        coeffs[i] = t_i // n
    for i in range(m):                              # R_k2_i
        coeffs[m + i] = sigs[i]['k2']
    coeffs[2 * m] = 1                               # R_c
    # t_0 correction via the n*e_0 vector = n*R_u + sum_{i>=1} beta_i R_t_i
    t0 = (u - lam * sigs[0]['k2'] - sigs[0]['k1'])
    assert t0 % n == 0
    t0 //= n
    coeffs[0] -= t0 * n
    for i in range(1, m):
        coeffs[i] -= t0 * beta[i]
    out = [0] * (2 * m + 1)
    for r in range(2 * m + 1):
        if coeffs[r] == 0: continue
        for c in range(2 * m + 1):
            out[c] += coeffs[r] * M[r][c]
    return out

def recover_d_elim(rows, sigs, n, lam, W, w1, w2, d_secret):
    m = len(sigs)
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    B0inv = modinv(B0, n)
    for row in rows:
        last = row[2 * m]
        if abs(last) != W: continue
        sgn = 1 if last > 0 else -1
        c_k1_0 = sgn * row[0]
        c_k2_0 = sgn * row[m]
        if c_k1_0 % w1 or c_k2_0 % w2: continue
        k1_0 = c_k1_0 // w1
        k2_0 = c_k2_0 // w2
        u = (k1_0 + lam * k2_0) % n
        d_cand = (u - A0) * B0inv % n
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def gh_lambda1(log_det, dim):
    """Gaussian heuristic  sqrt(dim/(2 pi e)) * det^(1/dim), via logs."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def baseline_geometry(m, n, k1_bound, k2_bound):
    """(log det of the lattice PROJECTED along the trivial n*e_m direction,
    dim, expected ||pv||).  The projection is the fair comparison: pv has
    ~zero component along e_m, and quotienting removes the trivial vector."""
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    ld = m * math.log(n * S_K1) + math.log(1) + m * math.log(S_K2) + math.log(n)
    ld_proj = ld - math.log(n)          # quotient by Z*(n*e_m), covol n
    pv = n * math.sqrt(2 * m / 3 + 4 / 3)
    return ld, ld_proj, 2 * m + 2, pv

def elim_geometry(m, n, k1_bound, k2_bound):
    W = k1_bound * k2_bound
    w1, w2 = k2_bound, k1_bound
    ld = (m - 1) * math.log(n) + m * math.log(w1) + m * math.log(w2) + math.log(W)
    pv = W * math.sqrt(2 * m / 3 + 1)
    return ld, 2 * m + 1, pv

# Closed form for the eliminated lattice's GH ratio.  With
#   det = n^(m-1) (K1 K2)^(m+1),  dim = 2m+1,  ||pv|| = K1 K2 sqrt(2m/3+1),
#   gh = sqrt(dim/(2 pi e)) det^(1/dim) / ||pv||
#      ~ sqrt(3/(2 pi e)) * [ n^((m-1)/m) / (K1 K2) ]^(m/(2m+1)).
# Verified against elim_geometry to within 7% over m in [6,32] (E5 below).

def gh_closed(m, n, k1_bound, k2_bound):
    return (math.sqrt(3 / (2 * math.pi * math.e))
            * (n ** ((m - 1) / m) / (k1_bound * k2_bound)) ** (m / (2 * m + 1)))

def eff_star(m, n, c=0.80):
    """Solve gh_closed = c for eff = K1*K2/n.  c = 0.80 is the empirical
    50%-success crossing measured in E3b."""
    return n ** (-1.0 / m) * (math.sqrt(3 / (2 * math.pi * math.e)) / c) ** ((2 * m + 1) / m)

# LLL does not have an SVP oracle: its first basis vector has norm
# ~ delta^dim * det^(1/dim) with delta ~ 1.021, so the achievable ratio is
# gh * sqrt(2 pi e / dim) / delta^dim.  This is what makes extra signatures
# stop paying (E4).
LLL_DELTA = 1.021

def lll_ratio(gh, dim):
    return gh * math.sqrt(2 * math.pi * math.e / dim) / (LLL_DELTA ** dim)

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def reduce_rows(M, dim, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def run_baseline(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m: return None
    M, S_K1, S_D, S_K2, S_KAN = build_glv_lattice(sigs, n, lam, k1_bound, k2)
    rows = reduce_rows(M, 2 * m + 2, use_bkz, beta)
    ok = recover_d_baseline(rows, m, n, S_KAN, d_secret) is not None
    return ok, rows, sigs

def run_elim(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m: return None
    M, W, w1, w2 = build_elim_lattice(sigs, n, lam, k1_bound, k2)
    rows = reduce_rows(M, 2 * m + 1, use_bkz, beta)
    ok = recover_d_elim(rows, sigs, n, lam, W, w1, w2, d_secret) is not None
    return ok, rows, sigs

def rate(fn, curve, m, k1_bound, seeds, **kw):
    n = curve[2]
    hits = 0
    for s in seeds:
        d = random.Random(s * 7919 + 13).randrange(1, n)
        r = fn(curve, m, d, k1_bound, s, **kw)
        if r is not None and r[0]:
            hits += 1
    return hits

def per_seed(curve, m, k1_bound, seeds, **kw):
    """Instance-by-instance outcomes; both formulations see the SAME signatures
    (same seed -> same gen_signatures stream), so this is a paired comparison."""
    n = curve[2]
    out = []
    for s in seeds:
        d = random.Random(s * 7919 + 13).randrange(1, n)
        rb = run_baseline(curve, m, d, k1_bound, s, **kw)
        re_ = run_elim(curve, m, d, k1_bound, s, **kw)
        out.append((bool(rb and rb[0]), bool(re_ and re_[0])))
    return out


# ===========================================================================
print("=" * 78)
print("Thread 23 -- eliminate d: does removing the trivial vector move the wall?")
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
    assert G is not None and (lam * lam + lam + 1) % n == 0
    hist.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E0: correctness of the eliminated lattice")
print("-" * 78)
print("Checks: (a) the planted vector is an integer combination of the basis,")
print("        (b) its coordinates are exactly (k1_i*w1, k2_i*w2, W),")
print("        (c) end-to-end d recovery on a known-good configuration.")
print()
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d = random.Random(1).randrange(1, n)
    sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, 42)
    M, W, w1, w2 = build_elim_lattice(sigs, n, lam, k1, k2b)
    pv_direct = elim_planted(sigs, n, lam, k1, k2b)
    pv_span = elim_planted_in_lattice(M, sigs, n, lam, k1, k2b)
    membership = (pv_direct == pv_span)
    rows = reduce_rows(M, 2 * m + 1)
    rec = recover_d_elim(rows, sigs, n, lam, W, w1, w2, d)
    print(f"  {label:<18} m={m:<3} K1={k1:<3} in-span={membership!s:<5} "
          f"coords-ok={all(pv_direct[i] == sigs[i]['k1']*w1 for i in range(m))!s:<5} "
          f"recover={'YES' if rec else 'no'}")
    assert membership, "planted vector is NOT in the eliminated lattice"

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1: geometry -- where does the planted vector sit?")
print("-" * 78)
print("sv/pv = ||shortest reduced row|| / ||planted||.  sv/pv >= 1 means the")
print("planted vector is (at least) as short as anything LLL found.")
print("gh = Gaussian-heuristic lambda_1 / ||pv||; gh > 1 predicts success.")
print()
print(f"  {'curve':<18} {'K1':>3} {'m':>3} | {'base sv/pv':>10} {'base gh*':>8} "
      f"{'triv/pv':>8} | {'elim sv/pv':>10} {'elim gh':>8}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d = random.Random(1).randrange(1, n)

    rb = run_baseline(curve, m, d, k1, 42)
    re_ = run_elim(curve, m, d, k1, 42)
    sigs = rb[2]
    pvb = baseline_planted(sigs, d, n, lam, k1, k2b)
    svb = min((norm(r) for r in rb[1] if any(r)), default=float('inf'))
    ld, ld_proj, dimb, _ = baseline_geometry(m, n, k1, k2b)
    ghb = gh_lambda1(ld_proj, dimb - 1) / norm(pvb)
    trivial = float(n)                       # ||n * S_D * e_m||, S_D = 1

    sigs2 = re_[2]
    pve = elim_planted(sigs2, n, lam, k1, k2b)
    sve = min((norm(r) for r in re_[1] if any(r)), default=float('inf'))
    lde, dime, _ = elim_geometry(m, n, k1, k2b)
    ghe = gh_lambda1(lde, dime) / norm(pve)

    print(f"  {label:<18} {k1:>3} {m:>3} | {svb/norm(pvb):>10.3f} {ghb:>8.3f} "
          f"{trivial/norm(pvb):>8.3f} | {sve/norm(pve):>10.3f} {ghe:>8.3f}")
print("\n  base gh* uses the lattice PROJECTED along the trivial n*e_m direction")
print("  (the fair comparison; the unprojected lattice has lambda_1 = n exactly).")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E2: the T4 grid -- baseline vs eliminated (Thread 23 falsifier)")
print("-" * 78)
print("Falsifier stated 2026-07-29: if the K1 wall moves outward the")
print("reformulation is a real improvement; if it stays, the wall is not the")
print("trivial vector's fault.  5 seeds per cell; entries are successes/5.")
print()
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
E2_CURVES = [(hist[1][0], hist[1][1], 8), (hist[2][0], hist[2][1], 10)]
hdr = "  {:<18} {:<5} ".format("curve", "form") + " ".join(f"{k:>4}" for k in K1_GRID)
print(hdr)
e2 = {}
for label, curve, m in E2_CURVES:
    for form, fn in (("base", run_baseline), ("elim", run_elim)):
        cells = [rate(fn, curve, m, k1, SEEDS) for k1 in K1_GRID]
        e2[(label, form)] = cells
        print("  {:<18} {:<5} ".format(label, form) +
              " ".join(f"{c:>4}" for c in cells))
    lam = curve[3]; n = curve[2]
    print(f"  {'':<18} (lam*={lam_star(lam,n):.3f}, m={m}, "
          f"K2={math.isqrt(n)+1}, n={n})")

print("\n  E2b: PAIRED per-instance agreement in the boundary cells.")
print("  Equal column totals could in principle hide disagreements that")
print("  cancel; this checks instance by instance (same signatures both ways).")
print(f"  {'curve':<18} {'K1':>3} | {'base':>5} {'elim':>5} {'agree':>6} "
      f"{'b-only':>7} {'e-only':>7}")
agree_tot = disagree_tot = 0
for label, curve, m in E2_CURVES:
    for k1 in (4, 6, 8, 12):
        rows_ps = per_seed(curve, m, k1, SEEDS)
        nb = sum(1 for a, _ in rows_ps if a)
        ne = sum(1 for _, b in rows_ps if b)
        ag = sum(1 for a, b in rows_ps if a == b)
        bo = sum(1 for a, b in rows_ps if a and not b)
        eo = sum(1 for a, b in rows_ps if b and not a)
        agree_tot += ag
        disagree_tot += len(rows_ps) - ag
        print(f"  {label:<18} {k1:>3} | {nb:>5} {ne:>5} {ag:>4}/{len(rows_ps)} "
              f"{bo:>7} {eo:>7}")
print(f"  TOTAL paired instances: {agree_tot + disagree_tot}, "
      f"agree {agree_tot}, disagree {disagree_tot}")

print("\n  E2c: does BKZ-40 rescue the eliminated lattice at the boundary?")
print(f"  {'curve':<18} {'K1':>3} | {'elim LLL':>9} {'elim BKZ40':>11} "
      f"{'base BKZ40':>11}")
for label, curve, m in E2_CURVES:
    for k1 in (8, 12, 16):
        e_lll = rate(run_elim, curve, m, k1, SEEDS)
        e_bkz = rate(run_elim, curve, m, k1, SEEDS, use_bkz=True, beta=40)
        b_bkz = rate(run_baseline, curve, m, k1, SEEDS, use_bkz=True, beta=40)
        print(f"  {label:<18} {k1:>3} | {e_lll:>9} {e_bkz:>11} {b_bkz:>11}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E3: eff = K1*K2/n sweep on fresh 17-bit curves")
print("-" * 78)
print("info bound: recovery is information-theoretically possible while")
print("eff < n^(-1/m).  GH bound: predicted lattice limit (gh > 1).")
print()
CURVES17 = search_curves(2**16, 2**17, 10)
print(f"  found {len(CURVES17)} fresh 17-bit j=0 GLV curves")
M17 = 12
n0 = CURVES17[0][2]
print(f"  m={M17}: info bound eff < {n0 ** (-1.0/M17):.3f}")
print()
print(f"  {'eff':>6} {'K1':>4} | {'base 5/5':>9} {'elim 5/5':>9} | "
      f"{'base gh*':>8} {'elim gh':>8}")
CALIB = []      # (gh_elim, base_ok, elim_ok) per instance
for eff in (0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.25, 0.35):
    nb = ne = 0
    ghb_s = ghe_s = 0.0
    k1s = []
    for curve in CURVES17:
        n = curve[2]
        k2b = math.isqrt(n) + 1
        k1 = max(2, int(round(eff * n / k2b)))
        k1s.append(k1)
        ps = per_seed(curve, M17, k1, SEEDS)
        lde, dime, pve = elim_geometry(M17, n, k1, k2b)
        gh_e = gh_lambda1(lde, dime) / pve
        for a, b in ps:
            CALIB.append((gh_e, a, b))
        if sum(1 for a, _ in ps if a) == len(SEEDS): nb += 1
        if sum(1 for _, b in ps if b) == len(SEEDS): ne += 1
        ld, ldp, dimb, pvb = baseline_geometry(M17, n, k1, k2b)
        ghb_s += gh_lambda1(ldp, dimb - 1) / pvb
        ghe_s += gh_e
    N = len(CURVES17)
    print(f"  {eff:>6.2f} {int(sum(k1s)/N):>4} | {nb:>4}/{N:<4} {ne:>4}/{N:<4} | "
          f"{ghb_s/N:>8.3f} {ghe_s/N:>8.3f}")

print("\n  E3b: calibration of the Gaussian-heuristic predictor.")
print("  All E3 instances pooled and binned by gh (eliminated lattice).")
print(f"  {'gh bin':>14} {'n':>5} {'base ok':>9} {'elim ok':>9}")
BINS = [(0.0, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 1.0),
        (1.0, 1.1), (1.1, 1.3), (1.3, 9.9)]
for lo_g, hi_g in BINS:
    sel = [r for r in CALIB if lo_g <= r[0] < hi_g]
    if not sel:
        continue
    ab = sum(1 for _, a, _ in sel if a)
    ae = sum(1 for _, _, b in sel if b)
    print(f"  [{lo_g:.2f},{hi_g:.2f}) {len(sel):>5} "
          f"{100.0*ab/len(sel):>8.0f}% {100.0*ae/len(sel):>8.0f}%")
n_dis = sum(1 for _, a, b in CALIB if a != b)
print(f"  paired instances {len(CALIB)}, base != elim on {n_dis}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E4: does extra data rescue the eliminated lattice?")
print("-" * 78)
print("T4b (2026-07-29) showed the baseline does NOT improve with m at the")
print("failing K1=8 on 12-bit/2677.  Same grid, both formulations.")
print()
label, curve, _ = (hist[2][0], hist[2][1], 10)
print(f"  curve {label},  K1=8,  successes/5")
print("  gh assumes an SVP oracle; r_LLL discounts it by LLL's Hermite factor")
print(f"  delta^dim with delta={LLL_DELTA}.")
print(f"  {'m':>4} {'dim':>5} {'base':>6} {'elim':>6} {'elim gh':>8} {'r_LLL':>8}")
for m in (8, 12, 16, 24, 32):
    n = curve[2]; k2b = math.isqrt(n) + 1
    rb = rate(run_baseline, curve, m, 8, SEEDS)
    re2 = rate(run_elim, curve, m, 8, SEEDS)
    lde, dime, pve = elim_geometry(m, n, 8, k2b)
    gh = gh_lambda1(lde, dime) / pve
    print(f"  {m:>4} {dime:>5} {rb:>6} {re2:>6} {gh:>8.3f} "
          f"{lll_ratio(gh, dime):>8.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E5: closed form for the wall")
print("-" * 78)
print("  gh ~ sqrt(3/(2 pi e)) * [ n^((m-1)/m) / (K1 K2) ]^(m/(2m+1))")
print("  =>  eff* = K1K2/n at the wall = n^(-1/m) * (sqrt(3/(2 pi e))/c)^((2m+1)/m)")
print()
print(f"  {'m':>4} {'n':>8} {'K1':>4} {'K2':>4} | {'gh exact':>9} {'gh closed':>10} {'ratio':>7}")
for (m, n, K1, K2) in [(12, 90000, 13, 301), (12, 90000, 25, 301),
                       (12, 90000, 89, 301), (8, 2659, 8, 52),
                       (10, 2647, 8, 52), (24, 2647, 8, 52), (32, 2647, 8, 52),
                       (6, 199, 2, 15)]:
    ld, dim, pv = elim_geometry(m, n, K1, K2)
    ex = gh_lambda1(ld, dim) / pv
    cl = gh_closed(m, n, K1, K2)
    print(f"  {m:>4} {n:>8} {K1:>4} {K2:>4} | {ex:>9.4f} {cl:>10.4f} {cl/ex:>7.4f}")

print()
print("  predicted wall (c = 0.80) vs measured 50%-success crossing:")
print(f"  {'setting':<28} {'eff* pred':>10} {'K1* pred':>9} {'K1* obs':>8}")
for lbl, m, n, K2, k1obs in [("12-bit/2557 (lam*=0.340) m=8", 8, 2659, 52, "~10"),
                             ("12-bit/2677 (lam*=0.070) m=10", 10, 2647, 52, "~5.5"),
                             ("17-bit sweep m=12", 12, 90000, 301, "~26")]:
    e = eff_star(m, n)
    print(f"  {lbl:<28} {e:>10.4f} {e*n/K2:>9.1f} {k1obs:>8}")
print()
print("  info-theoretic bound for comparison (eff < n^(-1/m)):")
for m, n in ((8, 2659), (10, 2647), (12, 90000)):
    print(f"    m={m:<3} n={n:<7} info {n**(-1.0/m):.3f}   GH(c=0.80) {eff_star(m,n):.3f}"
          f"   gap {n**(-1.0/m)/eff_star(m,n):.1f}x")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
