"""
GLV-HNP Phase 2, Thread 23: reformulate the Phase-2 lattice so the planted
vector becomes the BDD target (and remove the useless trivial vector).

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, "T5"):
  The Phase-2 Kannan lattice of glv_hnp_phase2_20bit.py:262 always has
  lambda_1 = n*S_D*e_m, the trivial vector supported entirely on the d
  column.  It is 2-3x shorter than the planted vector on every curve, so the
  planted vector is NEVER lambda_1 and recovery is a BDD/coset condition, not
  an SVP condition.  T4 measured a "K1 wall" at K1 ~ 4-6 (lam*=0.070) and
  K1 ~ 12-16 (lam*=0.340), and T4b showed more signatures do not move it.

  The 2026-07-29 next-step proposal (Thread 23) was:
    "project the lattice along e_m (quotient out the trivial n*e_m direction)
     and solve BDD in the projection, or replace the Kannan embedding with an
     explicit CVP call targeting (A_i*S_K1, 0, ...).  Falsifier: if sv/pv
     rises above 1 after the reformulation and the K1 wall in T4 moves
     outward on the lam*=0.070 curve (currently K1 ~ 4-6), the reformulation
     is a real improvement; if the wall stays at K1 ~ 4-6, then the wall is
     information-theoretic and Phase 2 is at its ceiling."

This script does exactly that, in three formulations.

-----------------------------------------------------------------------------
THE THREE LATTICES
-----------------------------------------------------------------------------
Unknowns per signature i:  A_i + B_i*d - lam*k2_i == k1_i  (mod n),
with k1_i in [0,K1), k2_i in [0,K2), d in [0,n).

L_K   ("kannan")  -- the historical lattice, verbatim.  dim 2m+2.
                    Columns: m k1-cols (scale S_K1), d-col (S_D=1),
                    m k2-cols (scale S_K2), Kannan col (S_KANNAN=n).
                    Planted coords are k1_i in [0,K1), k2_i in [0,K2),
                    d in [0,n):  ALL ONE-SIDED, so E[x^2] = X^2/3.

L_KC  ("kannanC")  -- L_K with the *target recentred*: substitute
                    k1_i = k1'_i + K1/2, k2_i = k2'_i + K2/2, d = d' + n/2,
                    which is achieved by the single substitution
                      A_i -> A_i - lam*(K2//2) - K1//2 + B_i*(n//2)   (mod n)
                    and reading d back as row[m] + n//2.  The lattice is
                    UNCHANGED (same basis, same determinant); only the coset
                    representative moves.  Now E[x'^2] = X^2/12, so the
                    expected planted norm drops by exactly 2x.  This is free:
                    K1, K2 are public in the bias model.

L_R   ("reduced")  -- d eliminated algebraically, so the trivial n*e_m
                    direction does not exist.  With t_i = B_i*B_0^{-1} mod n,
                    Lambda = {x in Z^m : x == d*B (mod n) for some d} has the
                    basis {(1,t_1,...,t_{m-1})} u {n*e_i : i>=1}, of
                    determinant n^{m-1} (not n^m).  So

                      dim 2m,  det = S_K1^m * n^(m-1) * S_K2^m

                    versus L_K's dim 2m+1 (ignoring the Kannan column) and
                    det = (n*S_K1)^m * S_D * S_K2^m = S_K1^m * n^m * S_K2^m.
                    One factor of n of determinant is gone AND the d
                    coordinate is gone from the distance.  d is recovered
                    afterwards from d = B_0^{-1}(k1_0 + lam*k2_0 - A_0) mod n.
                    Target is recentred as in L_KC.  Solved by Kannan
                    embedding (dim 2m+1) and, for small dim, cross-checked
                    against exact CVP enumeration.

-----------------------------------------------------------------------------
THE PREDICTOR
-----------------------------------------------------------------------------
For each formulation define the BDD ratio

    nu  =  E||residual||  /  lambda_1^GH(L),
    lambda_1^GH(L) = sqrt(dim/(2*pi*e)) * det(L)^(1/dim)

computed on the *underlying* CVP lattice (Kannan column excluded).  Closed
form, with eff = K1*K2/n and S_K1*K1 = S_K2*K2 = n:

    L_K :  nu = sqrt(2*pi*e/3)  * sqrt(eff) * n^( 1/(2m))  ~ 2.386*sqrt(eff)*...
    L_KC:  nu = sqrt(  pi*e/6)  * sqrt(eff) * n^( 1/(2m))  ~ 1.193*sqrt(eff)*...
    L_R :  nu = sqrt(  pi*e/6)  * sqrt(eff) * n^(-1/(2m))  ~ 1.193*sqrt(eff)*...

so L_KC buys a factor 2 in nu (= factor 4 in the eff budget = 2 bits of
nonce) over L_K, and L_R buys a further n^(1/m).  Note nu -> c*sqrt(eff) as
m -> infinity: the m-dependence is only n^(+-1/(2m)), which is why T4b saw no
rescue from m = 8..32 -- the required m is ~ log(n)/(2*log(1/(c*sqrt(eff)))).

Experiments:
  E1  diagnostics on the 3 historical curves: measured ||residual||/lambda_1,
      does the trivial vector still dominate, closed-form nu vs measured.
  E2  K1-wall sweep, all 3 formulations, 2 historical curves.  <- the falsifier
  E3  calibration: nu at the measured wall (is the wall at nu ~ 1?).
  E4  large-m test at fixed eff beyond the L_K wall: does m rescue, and at
      the m the predictor says?
  E5  exact-CVP oracle vs LLL on L_R: reduction weakness or information limit?

Run: python3 glv_hnp_phase2_bdd.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ, CVP

TWO_PI_E = 2.0 * math.pi * math.e

# ---------------------------------------------------------------------------
# Minimal EC / CM helpers (verbatim from glv_hnp_phase2_lambda_threshold.py)
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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_point(p, b, seed=1):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 with p == 1 mod 3 (norm form of Z[omega])."""
    for a in range(1, int(2 * math.isqrt(p)) + 2):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        r = math.isqrt(disc)
        if r * r == disc and (a + r) % 2 == 0:
            return (a, (a + r) // 2)
    return None

def j0_traces(a, b):
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0: continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

# ---------------------------------------------------------------------------
# Signatures (verbatim)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 400000:
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
        if B % n == 0: continue
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def lll_rows(M):
    dim_r, dim_c = len(M), len(M[0])
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(dim_c)] for i in range(dim_r)]

def bkz_rows(M, beta):
    dim_r, dim_c = len(M), len(M[0])
    A = IntegerMatrix.from_matrix(M)
    BKZ.reduction(A, BKZ.Param(beta))
    return [[A[i][j] for j in range(dim_c)] for i in range(dim_r)]

def shortest_nonzero(rows):
    best = None
    for r in rows:
        nr = norm(r)
        if nr > 0 and (best is None or nr < best):
            best = nr
    return best

def gh_lambda1(dim, log_det):
    """Gaussian-heuristic lambda_1 from log(det) (natural log), dim-dimensional."""
    return math.sqrt(dim / TWO_PI_E) * math.exp(log_det / dim)

# ===========================================================================
# L_K / L_KC  -- historical Kannan lattice, optionally recentred
# ===========================================================================

def build_kannan(sigs, n, lam, k1_bound, k2_bound, centred=False):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1 = k1_bound // 2 if centred else 0
    c2 = k2_bound // 2 if centred else 0
    cd = n // 2 if centred else 0
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -(lam % n) * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        At = (sigs[i]['A'] - lam * c2 - c1 + sigs[i]['B'] * cd) % n
        M[2 * m + 1][i] = At * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M, (c1, c2, cd)

def kannan_planted_residual(sigs, d_secret, n, k1_bound, k2_bound, centred):
    """Coordinates of the planted vector, Kannan column EXCLUDED (that is the
    CVP residual whose norm is compared against lambda_1 of the dim-(2m+1)
    underlying lattice)."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    c1 = k1_bound // 2 if centred else 0
    c2 = k2_bound // 2 if centred else 0
    cd = n // 2 if centred else 0
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - c2) * S_K2
    v[m] = (d_secret - cd) * S_D
    return v

def kannan_underlying_logdet(m, n, k1_bound, k2_bound):
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    return (m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2))

def kannan_recover(rows, m, n, S_KAN, cd, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m] + cd) % n
        if d_cand == d_secret:
            return d_cand
    return None

# ===========================================================================
# L_R  -- d eliminated, no trivial vector.  dim 2m.
# ===========================================================================

def build_reduced(sigs, n, lam, k1_bound, k2_bound):
    """Basis (2m x 2m) of
         L = {(S_K1*x, S_K2*y) : x == d*B - lam*y (mod n) for some d}
       plus the recentred target.  Returns (basis, target, aux)."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    t = [sigs[i]['B'] * B0inv % n for i in range(m)]   # t[0] == 1
    assert t[0] == 1
    lm = lam % n

    basis = []
    for i in range(m):                       # y = e_i, x = -lam*e_i
        row = [0] * (2 * m)
        row[i] = -lm * S_K1
        row[m + i] = S_K2
        basis.append(row)
    row = [0] * (2 * m)                      # d = B0inv  ->  x = t
    for i in range(m):
        row[i] = t[i] * S_K1
    basis.append(row)
    for i in range(1, m):                    # x = n*e_i
        row = [0] * (2 * m)
        row[i] = n * S_K1
        basis.append(row)
    assert len(basis) == 2 * m

    c1, c2 = k1_bound // 2, k2_bound // 2
    tgt = [0] * (2 * m)
    for i in range(m):
        tgt[i] = (-sigs[i]['A'] + c1) * S_K1
        tgt[m + i] = c2 * S_K2
    return basis, tgt, (S_K1, S_K2, c1, c2, B0inv)

def reduced_planted_residual(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    c1, c2 = k1_bound // 2, k2_bound // 2
    v = [0] * (2 * m)
    for i in range(m):
        v[i] = (c1 - sigs[i]['k1']) * S_K1
        v[m + i] = (c2 - sigs[i]['k2']) * S_K2
    return v

def reduced_logdet(m, n, k1_bound, k2_bound):
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    return m * math.log(S_K1) + (m - 1) * math.log(n) + m * math.log(S_K2)

def reduced_d_from_residual(res, sigs, n, lam, aux):
    """res = tgt - v.  Recover d, or None if res is not of the planted shape."""
    m = len(sigs)
    S_K1, S_K2, c1, c2, B0inv = aux
    for sign in (1, -1):
        k1s, k2s, ok = [], [], True
        for i in range(m):
            a, b = sign * res[i], sign * res[m + i]
            if a % S_K1 or b % S_K2:
                ok = False; break
            k1s.append(c1 - a // S_K1)
            k2s.append(c2 - b // S_K2)
        if not ok:
            continue
        d = B0inv * (k1s[0] + lam * k2s[0] - sigs[0]['A']) % n
        yield d

def solve_reduced(sigs, n, lam, k1_bound, k2_bound, d_secret,
                  algo="lll", beta=20, emb=None):
    """Kannan-embed the CVP instance for L_R and scan for the planted residual."""
    m = len(sigs)
    basis, tgt, aux = build_reduced(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    if emb is None:
        emb = max(1, n // 3)
    M = [row + [0] for row in basis]
    M.append(list(tgt) + [emb])
    rows = bkz_rows(M, beta) if algo == "bkz" else lll_rows(M)
    for row in rows:
        if abs(row[dim - 1]) != emb:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        res = [sign * row[j] for j in range(dim - 1)]
        for d in reduced_d_from_residual(res, sigs, n, lam, aux):
            if d == d_secret:
                return d
    return None

def solve_reduced_exact_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Exact CVP oracle (enumeration).  Small dim only."""
    m = len(sigs)
    basis, tgt, aux = build_reduced(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(basis)
    LLL.reduction(A)
    v = CVP.closest_vector(A, tuple(int(x) for x in tgt))
    res = [int(tgt[j]) - int(v[j]) for j in range(2 * m)]
    for d in reduced_d_from_residual(res, sigs, n, lam, aux):
        if d == d_secret:
            return d, norm(res)
    return None, norm(res)

# ===========================================================================
# Unified experiment driver
# ===========================================================================

def run_one(curve, m, d_secret, k1_bound, variant, seed=42, algo="lll", beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if variant in ("kannan", "kannanC"):
        centred = (variant == "kannanC")
        M, (c1, c2, cd) = build_kannan(sigs, n, lam, k1_bound, k2_bound, centred)
        rows = bkz_rows(M, beta) if algo == "bkz" else lll_rows(M)
        got = kannan_recover(rows, m, n, S_KAN, cd, d_secret)
        return got == d_secret
    elif variant == "reduced":
        got = solve_reduced(sigs, n, lam, k1_bound, k2_bound, d_secret,
                            algo=algo, beta=beta)
        return got == d_secret
    raise ValueError(variant)

def rate(curve, m, d_secret_fn, k1_bound, variant, seeds, algo="lll", beta=20):
    wins = 0
    for s in seeds:
        r = run_one(curve, m, d_secret_fn(curve), k1_bound, variant, seed=s,
                    algo=algo, beta=beta)
        if r: wins += 1
    return wins, len(seeds)

def nu_closed_form(variant, n, m, k1_bound, k2_bound):
    """Closed-form BDD ratio E||res|| / lambda_1^GH of the underlying lattice."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    if variant == "reduced":
        dim = 2 * m
        ld = reduced_logdet(m, n, k1_bound, k2_bound)
        var = (m * (k1_bound * S_K1) ** 2 + m * (k2_bound * S_K2) ** 2) / 12.0
    elif variant == "kannanC":
        dim = 2 * m + 1
        ld = kannan_underlying_logdet(m, n, k1_bound, k2_bound)
        var = (m * (k1_bound * S_K1) ** 2 + m * (k2_bound * S_K2) ** 2
               + (n * S_D) ** 2) / 12.0
    elif variant == "kannan":
        dim = 2 * m + 1
        ld = kannan_underlying_logdet(m, n, k1_bound, k2_bound)
        var = (m * (k1_bound * S_K1) ** 2 + m * (k2_bound * S_K2) ** 2
               + (n * S_D) ** 2) / 3.0
    else:
        raise ValueError(variant)
    return math.sqrt(var) / gh_lambda1(dim, ld)

# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]
VARIANTS = ["kannan", "kannanC", "reduced"]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,              p,    b, n,    lam,  K1, m
    ("8-bit/199",         211,  2, 199,  106,  2,  6),
    ("12-bit/2557",       2557, 2, 2659, 1755, 8,  12),
    ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  8,  12),
]

def hist_curves():
    out = []
    for label, p, b, n, lam, K1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        out.append((label, (p, b, n, lam, G), K1, m))
    return out

def dsec(curve):
    return curve[2] // 3 + 7          # deterministic, in [1, n)

print("=" * 78)
print("Thread 23 - reformulating the GLV-HNP Phase-2 lattice as BDD")
print("=" * 78)

CURVES = hist_curves()
for label, cur, K1, m in CURVES:
    p, b, n, lam, G = cur
    print(f"  {label:20s} p={p:6d} n={n:6d} ({n.bit_length()}b) lam={lam:6d} "
          f"lam*={lam_star(lam, n):.4f}  K2=isqrt(n)+1={math.isqrt(n)+1}")

# ===========================================================================
# E0 - sanity: is L_R really a basis, with the claimed determinant, and does
#      it really contain the planted vector?
# ===========================================================================
print()
print("=" * 78)
print("E0  sanity checks on L_R (basis / determinant n^(m-1) / planted in L)")
print("=" * 78)
for label, cur, K1, m in CURVES:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, dsec(cur), m, n, lam, p, K1, k2b, seed=42)
    basis, tgt, aux = build_reduced(sigs, n, lam, K1, k2b)
    # exact determinant via integer Hermite/Bareiss on the 2m x 2m basis
    Mat = sympy.Matrix(basis)
    det_meas = abs(Mat.det(method="berkowitz"))
    det_pred = (n // K1) ** m * n ** (m - 1) * max(1, n // k2b) ** m
    # planted residual must be tgt - (lattice point):  solve for coefficients
    res = reduced_planted_residual(sigs, n, K1, k2b)
    v = [tgt[j] - res[j] for j in range(2 * m)]
    coef = Mat.T.solve(sympy.Matrix(v))
    integral = all(c.q == 1 for c in coef)
    # and the recovered d must be d_secret
    ds = [d for d in reduced_d_from_residual(res, sigs, n, lam, aux)]
    print(f"  {label:20s} det==pred: {det_meas == det_pred}  "
          f"(det={det_meas}, pred={det_pred})")
    print(f"  {'':20s} planted in L (integral coeffs): {integral}   "
          f"d recovered == d_secret: {dsec(cur) in ds}")

# ===========================================================================
# E1 - diagnostics: is the planted vector the BDD target now?
# ===========================================================================
print()
print("=" * 78)
print("E1  measured BDD ratio  nu_meas = ||planted residual|| / lambda_1(LLL)")
print("    (lambda_1 measured on the UNDERLYING lattice, Kannan col dropped)")
print("=" * 78)
print(f"{'curve':20s} {'variant':9s} {'m':>3s} {'K1':>3s} {'||res||':>10s} "
      f"{'lam1_LLL':>10s} {'nu_meas':>8s} {'nu_cf':>8s} {'triv?':>6s}")

for label, cur, K1, m in CURVES:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, dsec(cur), m, n, lam, p, K1, k2b, seed=42)
    for variant in VARIANTS:
        if variant == "reduced":
            basis, tgt, aux = build_reduced(sigs, n, lam, K1, k2b)
            res = reduced_planted_residual(sigs, n, K1, k2b)
            rows = lll_rows(basis)
            l1 = shortest_nonzero(rows)
            triv = "n/a"
        else:
            centred = (variant == "kannanC")
            Mfull, (c1, c2, cd) = build_kannan(sigs, n, lam, K1, k2b, centred)
            # underlying lattice = drop Kannan row and Kannan column
            und = [row[:-1] for row in Mfull[:-1]]
            rows = lll_rows(und)
            l1 = shortest_nonzero(rows)
            res = kannan_planted_residual(sigs, dsec(cur), n, K1, k2b, centred)
            # is lambda_1 the trivial vector n*e_m ?
            sv = min(rows, key=lambda r: norm(r) if norm(r) > 0 else 1e300)
            triv = "YES" if (abs(sv[m]) == n and
                             all(x == 0 for j, x in enumerate(sv) if j != m)) else "no"
        nu_m = norm(res) / l1
        nu_c = nu_closed_form(variant, n, m, K1, k2b)
        print(f"{label:20s} {variant:9s} {m:3d} {K1:3d} {norm(res):10.1f} "
              f"{l1:10.1f} {nu_m:8.3f} {nu_c:8.3f} {triv:>6s}")

# ===========================================================================
# E2 - THE FALSIFIER: does the K1 wall move outward?
# ===========================================================================
print()
print("=" * 78)
print("E2  K1-wall sweep (m=12, 5 seeds, LLL).  Thread-23 falsifier.")
print("    eff = K1*K2/n;  wall = largest K1 with 5/5")
print("=" * 78)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
E2 = {}
for label, cur, _K1, _m in CURVES[1:]:          # the two 12-bit curves
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    print(f"\n  {label}  (lam*={lam_star(lam, n):.3f}, n={n}, K2={k2b})")
    hdr = "    " + f"{'variant':9s}" + "".join(f"{k:>7d}" for k in K1_GRID)
    print(hdr)
    print("    " + f"{'eff':9s}" + "".join(f"{k*k2b/n:7.3f}" for k in K1_GRID))
    for variant in VARIANTS:
        cells, wall = [], None
        for K1 in K1_GRID:
            w, t = rate(cur, 12, dsec, K1, variant, SEEDS)
            cells.append(f"{w}/{t}")
            if w == t:
                wall = K1
        E2[(label, variant)] = (dict(zip(K1_GRID, cells)), wall)
        print("    " + f"{variant:9s}" + "".join(f"{c:>7s}" for c in cells)
              + f"   wall={wall}")

print("\n  WALL SUMMARY (largest K1 with 5/5):")
for label, cur, _K1, _m in CURVES[1:]:
    row = [f"{v}={E2[(label, v)][1]}" for v in VARIANTS]
    base = E2[(label, 'kannan')][1]
    gains = []
    for v in VARIANTS[1:]:
        w = E2[(label, v)][1]
        gains.append(f"{v}/kannan = {w/base:.2f}x" if (w and base) else f"{v}: n/a")
    print(f"    {label:20s} " + "  ".join(row) + "   |  " + ", ".join(gains))

# ===========================================================================
# E3 - calibration: what nu does the wall sit at?
# ===========================================================================
print()
print("=" * 78)
print("E3  calibration - closed-form nu at the measured wall and one step past")
print("=" * 78)
print(f"{'curve':20s} {'variant':9s} {'K1_wall':>7s} {'nu@wall':>8s} "
      f"{'K1_next':>7s} {'nu@next':>8s}")
for label, cur, _K1, _m in CURVES[1:]:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    for variant in VARIANTS:
        cells, wall = E2[(label, variant)]
        if wall is None:
            print(f"{label:20s} {variant:9s} {'none':>7s}")
            continue
        idx = K1_GRID.index(wall)
        nxt = K1_GRID[idx + 1] if idx + 1 < len(K1_GRID) else None
        nu_w = nu_closed_form(variant, n, 12, wall, k2b)
        nu_n = nu_closed_form(variant, n, 12, nxt, k2b) if nxt else float('nan')
        print(f"{label:20s} {variant:9s} {wall:7d} {nu_w:8.3f} "
              f"{(nxt if nxt else 0):7d} {nu_n:8.3f}")

# ===========================================================================
# E4 - does m rescue beyond the wall?  (tests the n^(+-1/(2m)) term)
# ===========================================================================
print()
print("=" * 78)
print("E4  large-m at K1 fixed beyond the L_K wall (curve 2677, K1=8).")
print("    predictor: nu(m) = c*sqrt(eff)*n^(+-1/(2m)); success expected when nu<~1")
print("=" * 78)
label4, cur4, _, _ = CURVES[2]
p4, b4, n4, lam4, G4 = cur4
k2b4 = math.isqrt(n4) + 1
K1_4 = 8
eff4 = K1_4 * k2b4 / n4
print(f"  {label4}: n={n4}, K1={K1_4}, K2={k2b4}, eff={eff4:.4f}")
for variant in VARIANTS:
    # smallest m with nu < 1, if any
    m_pred = None
    for mm in range(2, 4001):
        if nu_closed_form(variant, n4, mm, K1_4, k2b4) < 1.0:
            m_pred = mm; break
    print(f"    {variant:9s} predicted m with nu<1: {m_pred}")
M_GRID = [12, 24, 48, 72, 96, 128]
print(f"\n    {'variant':9s}" + "".join(f"{f'm={mm}':>9s}" for mm in M_GRID)
      + f"{'':4s}nu(m) below")
for variant in VARIANTS:
    cells, nus = [], []
    for mm in M_GRID:
        w, t = rate(cur4, mm, dsec, K1_4, variant, SEEDS[:3])
        cells.append(f"{w}/{t}")
        nus.append(nu_closed_form(variant, n4, mm, K1_4, k2b4))
    print(f"    {variant:9s}" + "".join(f"{c:>9s}" for c in cells))
    print(f"    {'':9s}" + "".join(f"{x:9.3f}" for x in nus))

# ===========================================================================
# E5 - exact CVP oracle vs LLL on L_R:  reduction weakness or info limit?
# ===========================================================================
print()
print("=" * 78)
print("E5  L_R: exact-CVP oracle vs LLL-embedding vs BKZ(30), m=8 (dim 16)")
print("    'oracle=1' means the planted residual IS the closest vector.")
print("=" * 78)
print(f"{'curve':20s} {'K1':>3s} {'eff':>6s} {'nu_cf':>7s} {'LLL':>6s} "
      f"{'BKZ30':>6s} {'CVPoracle':>10s} {'d_res/d_pl':>11s}")
for label, cur, _K1, _m in CURVES[1:]:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    for K1 in [4, 8, 16, 24, 32, 48]:
        m5 = 8
        wl, _ = rate(cur, m5, dsec, K1, "reduced", SEEDS[:3])
        wb, _ = rate(cur, m5, dsec, K1, "reduced", SEEDS[:3], algo="bkz", beta=30)
        oracle, ratios = 0, []
        for s in SEEDS[:3]:
            sigs = gen_signatures(G, dsec(cur), m5, n, lam, p, K1, k2b, seed=s)
            if len(sigs) < m5: continue
            got, dres = solve_reduced_exact_cvp(sigs, n, lam, K1, k2b, dsec(cur))
            dpl = norm(reduced_planted_residual(sigs, n, K1, k2b))
            if got is not None: oracle += 1
            ratios.append(dres / dpl if dpl else float('nan'))
        nu = nu_closed_form("reduced", n, m5, K1, k2b)
        rr = sum(ratios) / len(ratios) if ratios else float('nan')
        print(f"{label:20s} {K1:3d} {K1*k2b/n:6.3f} {nu:7.3f} {wl:>4d}/3 "
              f"{wb:>4d}/3 {oracle:>8d}/3 {rr:11.4f}")

# ===========================================================================
# E6 - generalisation: fresh 17-bit curves spanning lam*.  Is the centring
#      gain a 12-bit artifact?  Is the nu* ~ 1 calibration stable?
# ===========================================================================
print()
print("=" * 78)
print("E6  17-bit generalisation: K1 wall per formulation, 6 fresh curves")
print("=" * 78)

def search_curves(lo, hi, want=6, nbins=6):
    """j=0 GLV curves with p in [lo,hi), n prime, n==1 mod 3, spread by lam*."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, nc) / (0.5 / nbins)))
                    if bins[idx]:
                        continue
                    cu = build_curve(p, nc)
                    if cu is None:
                        continue
                    bins[idx].append((p, cu[1], nc, lam, cu[4]))
        if sum(len(v) for v in bins.values()) >= want:
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

C17 = search_curves(2 ** 16, 2 ** 17, want=6)
print(f"  found {len(C17)} 17-bit curves\n")
K1_17 = [4, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256]
M17 = 12
print(f"  {'p':>7s} {'n':>7s} {'lam*':>6s} | "
      + " ".join(f"{v:>11s}" for v in VARIANTS) + "   gain(C/K, R/K)")
walls17 = []
for cu in C17:
    p, b, n, lam, G = cu
    k2b = math.isqrt(n) + 1
    ws = {}
    for variant in VARIANTS:
        wall = None
        for K1 in K1_17:
            w, t = rate(cu, M17, dsec, K1, variant, SEEDS[:3])
            if w == t:
                wall = K1
        ws[variant] = wall
    walls17.append((cu, ws))
    gk = ws['kannan']
    g1 = f"{ws['kannanC']/gk:.2f}x" if (gk and ws['kannanC']) else "n/a"
    g2 = f"{ws['reduced']/gk:.2f}x" if (gk and ws['reduced']) else "n/a"
    cells = []
    for variant in VARIANTS:
        wl = ws[variant]
        nu = nu_closed_form(variant, n, M17, wl, k2b) if wl else float('nan')
        cells.append(f"{str(wl):>4s}(nu{nu:.2f})")
    print(f"  {p:7d} {n:7d} {lam_star(lam, n):6.4f} | " + " ".join(cells)
          + f"   {g1}, {g2}")

ratios_c = [ws['kannanC'] / ws['kannan'] for _, ws in walls17
            if ws['kannan'] and ws['kannanC']]
ratios_r = [ws['reduced'] / ws['kannan'] for _, ws in walls17
            if ws['kannan'] and ws['reduced']]
nus = [nu_closed_form(v, cu[2], M17, ws[v], math.isqrt(cu[2]) + 1)
       for cu, ws in walls17 for v in VARIANTS if ws[v]]
print()
print(f"  kannanC/kannan wall ratio: {ratios_c}  (geo-mean "
      f"{math.exp(sum(map(math.log, ratios_c))/len(ratios_c)):.2f}x)" if ratios_c else "")
print(f"  reduced/kannan wall ratio: {ratios_r}  (geo-mean "
      f"{math.exp(sum(map(math.log, ratios_r))/len(ratios_r)):.2f}x)" if ratios_r else "")
print(f"  nu at wall over all cells: min={min(nus):.3f} max={max(nus):.3f} "
      f"mean={sum(nus)/len(nus):.3f}  (n={len(nus)} cells)")

# ===========================================================================
# E7 - the ceiling.  How far is the reformulated attack from the
#      information-theoretic (counting) limit?
# ===========================================================================
print()
print("=" * 78)
print("E7  ceiling: measured wall vs Gaussian-heuristic vs counting bound")
print("=" * 78)
print("""    Counting bound: for a WRONG d, all m signatures land in the
    K1 x K2 box with probability eff^m, and there are ~n wrong d, so the
    solution is unique only while  n*eff^m < 1, i.e.  eff < n^(-1/m).
    Gaussian heuristic (nu<1) gives  eff < (1/1.193^2)*n^(-1/m)
                                          = 0.703*n^(-1/m)  for L_R,
    so GH is a constant factor 0.70 inside the counting bound; both scale as
    n^(-1/m).  c_emp = eff_at_wall / n^(-1/m) measures where we actually are.""")
print()
print(f"  {'curve':>16s} {'m':>3s} {'variant':9s} {'K1w':>5s} {'eff_w':>7s} "
      f"{'n^-1/m':>7s} {'c_emp':>6s} {'K1_count':>9s}")
ROWS7 = []
for label, cur, _K1, _m in CURVES[1:]:
    ROWS7.append((label, cur, 12, {v: E2[(label, v)][1] for v in VARIANTS}))
for cu, ws in walls17:
    ROWS7.append((f"17b/{cu[0]}", cu, M17, ws))
c_emp_by_variant = {v: [] for v in VARIANTS}
for label, cu, mm, ws in ROWS7:
    n = cu[2]
    k2b = math.isqrt(n) + 1
    nb = n ** (-1.0 / mm)
    for v in VARIANTS:
        wl = ws.get(v)
        if not wl:
            continue
        eff_w = wl * k2b / n
        c_emp = eff_w / nb
        c_emp_by_variant[v].append(c_emp)
        print(f"  {label:>16s} {mm:3d} {v:9s} {wl:5d} {eff_w:7.3f} {nb:7.3f} "
              f"{c_emp:6.3f} {nb * n / k2b:9.1f}")
print()
for v in VARIANTS:
    xs = c_emp_by_variant[v]
    if not xs:
        continue
    print(f"  c_emp[{v:8s}]  n={len(xs):2d}  min={min(xs):.3f} "
          f"max={max(xs):.3f} mean={sum(xs)/len(xs):.3f}   "
          f"(GH predicts 0.703 for L_R / 0.176 for L_K, counting bound = 1.0)")

print()
print("Done.")
