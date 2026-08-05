"""
GLV-HNP Phase 2, Thread 23: drop the Kannan embedding, solve the BDD directly.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, result T5):
  The Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` (build_glv_lattice) is a
  (2m+2)-dimensional Kannan embedding with columns

      [ residue_0 .. residue_{m-1} | d | k2_0 .. k2_{m-1} | kannan ]

  and it contains the vector  n*S_D*e_m  (obtained as n*row_m - sum_i B_i*row_i,
  every residue column being a multiple of n*S_K1).  Its norm is n*S_D, while the
  planted vector has norm ~ n*sqrt(2m/3 + 4/3).  So the planted vector is NEVER
  lambda_1 -- measured on every curve tested, |sv[m]|/n = 1.0000 exactly and
  sv/pv in [0.337, 0.368].  Recovery was therefore never an SVP condition;
  it was a BDD/coset condition that LLL was being asked to solve by accident.

  The 2026-07-29 run proposed Thread 23: reformulate so that the target is
  lambda_1, i.e. quotient out the trivial direction and solve BDD/CVP directly.

THIS SCRIPT.  The reformulation used here removes the trivial direction by
construction rather than by projection: the d-column and the Kannan row are
simply deleted, and d is recovered arithmetically from (k1_0, k2_0) afterwards.

  Lattice L (dimension 2m, columns [residue_0..residue_{m-1} | k2_0..k2_{m-1}]):
      n*S_K1*e_i                          i = 0..m-1     (residue wraparound)
      b   = (B_i*S_K1)_i ; 0                             (coefficient d)
      r_i = -lam*S_K1 on residue col i, S_K2 on k2 col i (coefficient k2_i)

  Target  t_i = -A_i*S_K1  on the residue columns, 0 on the k2 columns.

  For the true secret,  d*B_i + A_i = k1_i + lam*k2_i  (mod n),  so the lattice
  point  v = d*b + sum_i k2_i*r_i  (reduced mod the n*S_K1*e_i) satisfies

      v - t  =  ( k1_i*S_K1 )_i ; ( k2_i*S_K2 )_i

  which is short.  Recovery: k1_i = (v-t)_i / S_K1, k2_i = (v-t)_{m+i} / S_K2,
  then  d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n).

  There is no d-column, hence no n*S_D*e_m vector, hence no trivial obstruction.

  det(L) = S_K1^m * n^(m-1) * S_K2^m      (the (B_i) generator has order n in
                                           Z^m / nZ^m, contributing index n)

Gaussian-heuristic prediction (derived here, tested in E3).  With S_K1 = n/K1,
S_K2 = n/K2, dim = 2m:

    GH(L) = sqrt(m/(pi*e)) * n^1.5 / sqrt(K1*K2)        (expected dist to a
                                                         random target)
    dist_uncentered = n*sqrt(2m/3)      (k1 ~ U[0,K1), so E[k1^2] = K1^2/3)
    dist_centered   = n*sqrt(m/6)       (target shifted to the box centre)

  dist < GH  <=>   eff := K1*K2/n  <  3/(2*pi*e) = 0.1757    (uncentered)
                   eff             <  6/(pi*e)   = 0.7027    (centred)

  So the reformulation is predicted to have an eff wall, NOT a K1 wall and NOT a
  lam wall -- and centring the target is predicted to move that wall by 4x.

  *** THIS GH PREDICTION IS FALSIFIED BY E3 -- DO NOT RE-DERIVE IT. ***
  Measured walls: uncentred ~0.05-0.10 (predicted 0.176), centred ~0.18-0.30
  (predicted 0.703).  GH overestimates by 2-4x because L is nowhere near a
  random lattice: lambda_1(L) comes from a single 2-D lambda block and sits well
  below GH(L) (E1: lambda_1/n = 1.08-1.93 against GH/n ~ 2.1-2.4).  The *ratio*
  centred:uncentred ~ 3x does survive, so the centring half of the derivation is
  sound; only the absolute constant is not.

EXPERIMENTS
  E1  Structural check.  For the three historical curves: is the trivial vector
      gone, and does the planted vector become lambda_1 (sv/pv > 1)?
  E2  Replay of the 2026-07-29 T4 grid (the stated falsifier).  Curves
      12-bit/2557 (lam*=0.340) and 12-bit/2677 (lam*=0.070), K2=52, m=12,
      K1 in {2,3,4,6,8,12,16,24}, 5 seeds, five methods:
        KANNAN-LLL     the 2026-07-29 baseline, verbatim
        KANNAN-1GUESS  ditto but scored on ONE candidate (see below)
        CVP-BABAI      LLL + Babai nearest-plane, uncentred target
        CVP-BABAI-C    LLL + Babai nearest-plane, centred target
        CVP-EXACT-C    exact CVP by enumeration, centred target
      CVP-EXACT-C answers the information-theoretic question (is the planted
      point the closest lattice point at all?); the Babai rows answer the
      algorithmic one.  KANNAN-1GUESS controls for a scoring asymmetry: the
      published recover_d is handed d_secret and allowed to match ANY of the
      2m+2 rows, whereas the CVP decoder gets exactly one candidate.  The
      control shows the asymmetry is immaterial (KANNAN-1GUESS ~ KANNAN-LLL).
  E3  eff sweep on fresh 17-bit curves, to locate the wall of the reformulated
      attack and compare it against the 0.176 / 0.703 predictions.
  E4  predictor study: exact lambda_1(L)/dist vs eff vs lam*, marginal AUC and
      within-eff-stratum AUC over 180 instances.
  E5  head-to-head KANNAN vs CVP-BABAI-C over 12 x 12-bit and 12 x 17-bit fresh
      curves at eff in {0.08,0.12,0.18,0.24,0.30}: does the reformulation
      dominate, or are the two formulations complementary?

HEADLINE RESULTS (2026-08-05 run; regenerate with the command below)
  * Thread 23's falsifier condition is MET.  sv/pv rises from 0.337-0.368 (old)
    to 0.690-1.925, and on the lam*=0.070 curve the K1 wall moves from ~4-6 out
    to ~16-24.  The 2026-07-29 K1 wall was a formulation artifact, not an
    information-theoretic ceiling.
  * The gain is general: CVP-BABAI-C beats KANNAN-LLL in every one of the ten
    (bit-size, eff) cells of E5, by up to 9x (17-bit, eff=0.18: 4/60 -> 37/60).
  * Centring the target is the essential ingredient, not the CVP framing:
    UNCENTRED CVP is worse than the Kannan baseline everywhere.
  * CVP-EXACT-C ~ CVP-BABAI-C everywhere: in this formulation nearest-plane is
    already optimal, so stronger reduction (the 2026-07-26 "BKZ(40) rescue")
    cannot help.  Residual failures are information-theoretic for this lattice.
  * lam* remains a non-predictor (E4 AUC 0.457), reconfirming 2026-07-29 T1/T3.

Run: python3 glv_hnp_phase2_cvp.py
Deps: fpylll, sympy  (pip install cysignals fpylll sympy)
Note: fpylll 0.6.4 wheels ship no pruning-strategies file, so SVP.shortest_vector
and BKZ.Param(..., strategies=...) raise FileNotFoundError; full-block BKZ works.
"""

import math
import random
import sys
import time

import sympy
from fpylll import BKZ, CVP, GSO, LLL, IntegerMatrix

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- copied verbatim from
# glv_hnp_phase2_20bit.py so this script stands alone.
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
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory for j=0 curves (copied from glv_hnp_phase2_lambda_threshold.py)
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
    """Both roots of x^2 + x + 1 = 0 mod n, or (None, None)."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    return min(r1, r2), max(r1, r2)

def lam_star(lam, n):
    """Root-convention-free lambda statistic (see 2026-07-29 T1)."""
    return min(lam, n - lam) / n

def build_curve(p, n, seed=12345):
    """Find twist b with #E = n and a generator. Returns (p,b,n,lam,G) or None."""
    lam, _ = glv_roots(n)
    if lam is None:
        return None
    for b_try in range(1, 400):
        rng = random.Random(seed)
        P = None
        for _ in range(500):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                break
        if P is None:
            continue
        if ec_mul(P, n, p) is None:
            G = find_generator(p, b_try, n, seed)
            if G is not None:
                return (p, b_try, n, lam, G)
    return None

def search_curves(lo_bits, hi_bits, want, seed_prime=None):
    """Prime-order j=0 GLV curves with n of the given bit length."""
    out = []
    p = sympy.nextprime(seed_prime if seed_prime else 2 ** lo_bits)
    while len(out) < want and p < 2 ** hi_bits:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 4 or n_cand.bit_length() < lo_bits:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        cur = build_curve(p, n_cand)
                        if cur is not None:
                            out.append(cur)
                            break
        p = sympy.nextprime(p)
    return out

# ---------------------------------------------------------------------------
# Signature generation -- identical convention to glv_hnp_phase2_20bit.py:236
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

def scales(n, k1_bound, k2_bound):
    """Column scales, identical to build_glv_lattice (glv_hnp_phase2_20bit.py:262)."""
    return n // k1_bound, 1, max(1, n // k2_bound), n

# ---------------------------------------------------------------------------
# BASELINE: the (2m+2)-dim Kannan lattice, verbatim from 2026-07-29
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

def kannan_planted_vector(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    v = [0] * dim
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
    v[m] = d_secret * S_D
    for i in range(m):
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[dim - 1] = S_KANNAN
    return v

def kannan_recover(M_reduced, m, n, S_KANNAN, d_secret):
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

def kannan_recover_strict(M_reduced, m, n, S_KANNAN, d_secret):
    """Single-candidate version of kannan_recover.

    The published recover_d (glv_hnp_phase2_20bit.py:288) scans ALL 2m+2 rows and
    declares success if ANY of them yields d_secret -- i.e. it is handed the answer
    and allowed ~26 guesses.  The CVP decoder below gets exactly one candidate, so
    comparing the two directly would flatter the Kannan baseline.  This variant
    takes only the min-norm row with |last| == S_KANNAN, making the comparison
    one-guess-vs-one-guess.
    """
    dim = 2 * m + 2
    cands = [row for row in M_reduced if abs(row[dim - 1]) == S_KANNAN]
    if not cands:
        return None
    row = min(cands, key=norm)
    sign = 1 if row[dim - 1] > 0 else -1
    d_cand = (sign * row[m]) % n
    return d_cand if d_cand == d_secret else None

def run_kannan(curve, m, d_secret, k1_bound, seed=42, strict=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KANNAN = scales(n, k1_bound, k2_bound)
    fn = kannan_recover_strict if strict else kannan_recover
    return fn(red, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# THREAD 23: the reformulated BDD lattice (dimension 2m, no d, no Kannan)
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Rows of L; 2m+1 generators of a rank-2m lattice (b has order n)."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                                   # residue wraparound
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                                        # coefficient d
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                   # coefficient k2_i
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    return rows

def cvp_target(sigs, n, k1_bound, k2_bound, centred):
    """t on the residue columns; centring shifts t into the middle of the box."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    t = [0] * dim
    off1 = (k1_bound // 2) if centred else 0
    off2 = (k2_bound // 2) if centred else 0
    for i in range(m):
        t[i] = (-sigs[i]['A'] + off1) * S_K1
        t[m + i] = off2 * S_K2
    return t

def cvp_planted_offset(sigs, n, k1_bound, k2_bound, centred):
    """The vector v - t for the true secret (what CVP must find)."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    off1 = (k1_bound // 2) if centred else 0
    off2 = (k2_bound // 2) if centred else 0
    return ([(s['k1'] - off1) * S_K1 for s in sigs] +
            [(s['k2'] - off2) * S_K2 for s in sigs])

def reduced_basis(rows, dim):
    """LLL-reduce the generating set and drop the zero rows (rank 2m of 2m+1)."""
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    keep = []
    for i in range(A.nrows):
        row = [A[i][j] for j in range(dim)]
        if any(row):
            keep.append(row)
    return IntegerMatrix.from_matrix(keep)

def babai(B, t):
    """Nearest-plane on an LLL-reduced B. Returns the lattice vector."""
    M = GSO.Mat(B, float_type="ld")
    M.update_gso()
    w = M.babai(list(t))
    out = [0] * B.ncols
    for i, c in enumerate(w):
        if c:
            for j in range(B.ncols):
                out[j] += c * B[i][j]
    return out

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def decode_cvp(v, t, sigs, n, k1_bound, k2_bound, centred, d_secret):
    """Turn a candidate lattice point into d; return d or None."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    off1 = (k1_bound // 2) if centred else 0
    off2 = (k2_bound // 2) if centred else 0
    diff = [v[i] - t[i] for i in range(2 * m)]
    if any(x % S_K1 for x in diff[:m]) or any(x % S_K2 for x in diff[m:]):
        return None
    k1s = [diff[i] // S_K1 + off1 for i in range(m)]
    k2s = [diff[m + i] // S_K2 + off2 for i in range(m)]
    for i in range(m):
        if not (0 <= k1s[i] < k1_bound and 0 <= k2s[i] < k2_bound):
            return None
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    try:
        d = (k1s[0] + lam_cache[0] * k2s[0] - A0) * modinv(B0, n) % n
    except ValueError:
        return None
    return d if d == d_secret else None

lam_cache = [0]   # set per-experiment; keeps decode_cvp's signature small

def run_cvp(curve, m, d_secret, k1_bound, seed=42, centred=True, exact=False):
    """Returns (recovered_bool, diag dict) or (None, None) if sig gen failed."""
    p, b, n, lam, G = curve
    lam_cache[0] = lam
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, None
    dim = 2 * m
    rows = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    B = reduced_basis(rows, dim)
    t = cvp_target(sigs, n, k1_bound, k2_bound, centred)
    pv = cvp_planted_offset(sigs, n, k1_bound, k2_bound, centred)

    if exact:
        v = list(CVP.closest_vector(B, list(t)))
    else:
        v = babai(B, t)

    d = decode_cvp(v, t, sigs, n, k1_bound, k2_bound, centred, d_secret)
    lam1 = norm([B[0][j] for j in range(dim)])
    diag = {
        'rank': B.nrows,
        'dim': dim,
        'lam1': lam1,
        'dist_planted': norm(pv),
        'dist_found': norm([v[i] - t[i] for i in range(dim)]),
        'svpv': lam1 / norm(pv) if norm(pv) else float('inf'),
    }
    return (d is not None), diag

# ---------------------------------------------------------------------------
# Heuristics
# ---------------------------------------------------------------------------

def gh_distance(n, m, k1_bound, k2_bound):
    """Gaussian-heuristic expected distance from a random target to L."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    log_det = m * math.log(S_K1) + (m - 1) * math.log(n) + m * math.log(S_K2)
    dim = 2 * m
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

EFF_PRED_UNCENTRED = 3.0 / (2 * math.pi * math.e)     # 0.1757
EFF_PRED_CENTRED = 6.0 / (math.pi * math.e)           # 0.7027

# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]
T0 = time.time()

def hdr(s):
    print("\n" + "=" * 74)
    print(s)
    print("=" * 74)
    sys.stdout.flush()

print("=" * 74)
print("GLV-HNP Phase 2 / Thread 23 — BDD reformulation (no Kannan, no d-column)")
print("=" * 74)
print(f"predicted eff wall: uncentred {EFF_PRED_UNCENTRED:.4f}, "
      f"centred {EFF_PRED_CENTRED:.4f}")

# The three historical curves (see 2026-07-29 T1/T5 table).
CURVE_8 = (211, 2, 199, glv_roots(199)[0], find_generator(211, 2, 199))
CURVE_2557 = (2557, 2, 2659, glv_roots(2659)[0], find_generator(2557, 2, 2659))
CURVE_2677 = (2677, 2, 2647, glv_roots(2647)[0], find_generator(2677, 2, 2647))
HIST = [("8-bit/199", CURVE_8, 2),
        ("12-bit/2557", CURVE_2557, 8),
        ("12-bit/2677", CURVE_2677, 8)]

# ---- E1: structural check --------------------------------------------------

hdr("E1  Is the trivial vector gone, and is the planted vector now lambda_1?")
print(f"{'curve':<14}{'lam*':>7}{'K1':>4}{'dim':>5}{'rank':>6}"
      f"{'lam_1/n':>10}{'dist/n':>9}{'sv/pv':>8}{'GH/dist':>9}{'eff':>8}")
for label, cur, k1b in HIST:
    p, b, n, lam, G = cur
    m = 12
    k2b = math.isqrt(n) + 1
    d_secret = random.Random(4242).randint(1, n - 1)
    ok, dg = run_cvp(cur, m, d_secret, k1b, seed=42, centred=True, exact=False)
    gh = gh_distance(n, m, k1b, k2b)
    eff = k1b * k2b / n
    print(f"{label:<14}{lam_star(lam, n):>7.3f}{k1b:>4}{dg['dim']:>5}{dg['rank']:>6}"
          f"{dg['lam1']/n:>10.3f}{dg['dist_planted']/n:>9.3f}"
          f"{dg['svpv']:>8.3f}{gh/dg['dist_planted']:>9.3f}{eff:>8.4f}")
print("\n  sv/pv = lambda_1(L) / ||v_planted - t||.  In the OLD Kannan lattice this")
print("  ratio was 0.337-0.368 on every curve (2026-07-29 T5) because of the")
print("  trivial vector n*S_D*e_m.  Values > 1 here mean the planted offset is")
print("  shorter than any nonzero lattice vector, i.e. the BDD solution is unique.")

# ---- E2: the stated falsifier — replay of the T4 grid ----------------------

hdr("E2  T4 grid replay: K1 wall for each method (m=12, K2=52, 5 seeds)")
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
METHODS = [
    ("KANNAN-LLL  ", lambda c, m, d, k, s: run_kannan(c, m, d, k, s)),
    ("KANNAN-1GUESS", lambda c, m, d, k, s: run_kannan(c, m, d, k, s, strict=True)),
    ("CVP-BABAI   ", lambda c, m, d, k, s: run_cvp(c, m, d, k, s, centred=False, exact=False)[0]),
    ("CVP-BABAI-C ", lambda c, m, d, k, s: run_cvp(c, m, d, k, s, centred=True, exact=False)[0]),
    ("CVP-EXACT-C ", lambda c, m, d, k, s: run_cvp(c, m, d, k, s, centred=True, exact=True)[0]),
]
e2 = {}
for label, cur, _ in HIST[1:]:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    print(f"\n  {label}  n={n}  lam={lam}  lam*={lam_star(lam, n):.3f}  K2={k2b}")
    print(f"  {'method':<13}" + "".join(f"{f'K1={k}':>9}" for k in K1_GRID))
    print(f"  {'eff=K1K2/n':<13}" + "".join(f"{k*k2b/n:>9.3f}" for k in K1_GRID))
    for mname, fn in METHODS:
        cells = []
        for k1b in K1_GRID:
            wins = 0
            for s in SEEDS:
                d_trial = random.Random(s + 7777).randint(1, n - 1)
                try:
                    wins += bool(fn(cur, 12, d_trial, k1b, s))
                except Exception:
                    pass
            cells.append(wins)
            e2[(label, mname.strip(), k1b)] = wins
        print(f"  {mname:<13}" + "".join(f"{f'{c}/5':>9}" for c in cells))
        sys.stdout.flush()

# ---- E3: where is the wall, and does it track eff? -------------------------

hdr("E3  eff sweep on fresh 17-bit curves (m=12, 5 seeds, CVP-EXACT-C)")
curves17 = search_curves(17, 18, 6)
print(f"  found {len(curves17)} fresh 17-bit prime-order j=0 GLV curves")
for c in curves17:
    print(f"    p={c[0]:<7} n={c[2]:<7} lam={c[3]:<7} lam*={lam_star(c[3], c[2]):.4f}")

EFFS = [0.05, 0.10, 0.18, 0.30, 0.45, 0.60, 0.70, 0.85, 1.00]
print(f"\n  {'eff':>6}{'K1(typ)':>9}" + "".join(f"{n:>14}" for n in
      ("BABAI-uncent", "BABAI-cent", "EXACT-cent")))
e3 = {}
for eff in EFFS:
    tallies = {'u': [0, 0], 'c': [0, 0], 'x': [0, 0]}
    k1_typ = None
    for cur in curves17:
        p, b, n, lam, G = cur
        k2b = math.isqrt(n) + 1
        k1b = max(2, round(eff * n / k2b))
        if k1_typ is None:
            k1_typ = k1b
        for s in SEEDS:
            d_trial = random.Random(s + 7777).randint(1, n - 1)
            for key, kw in (('u', dict(centred=False, exact=False)),
                            ('c', dict(centred=True, exact=False)),
                            ('x', dict(centred=True, exact=True))):
                try:
                    ok, _ = run_cvp(cur, 12, d_trial, k1b, s, **kw)
                except Exception:
                    ok = False
                tallies[key][1] += 1
                tallies[key][0] += bool(ok)
    e3[eff] = tallies
    print(f"  {eff:>6.2f}{k1_typ:>9}" +
          "".join(f"{f'{tallies[k][0]}/{tallies[k][1]}':>14}" for k in ('u', 'c', 'x')))
    sys.stdout.flush()

print(f"\n  prediction: uncentred wall at eff={EFF_PRED_UNCENTRED:.3f}, "
      f"centred wall at eff={EFF_PRED_CENTRED:.3f}")

# ---- E4: what actually predicts success in the reformulated lattice? -------
#
# E1+E2 show that at essentially IDENTICAL eff (0.156 vs 0.157) the two 12-bit
# curves land on opposite sides of the wall (1/5 vs 5/5), so eff alone cannot be
# the predictor and the Gaussian heuristic (which only sees eff) cannot be
# either.  The quantity that does separate them in E1 is
# lambda_1(L) / ||v_planted - t||.
#
# Note this is the SAME variable as the 2026-07-29 T2 hypothesis (rho = mu/||v||,
# mu = shortest vector of the 2-D lambda block), which was FALSIFIED there --
# for the Kannan lattice.  E4 tests it for the reformulated lattice, where mu is
# genuinely lambda_1 because the trivial vector n*S_D*e_m no longer exists.
# lambda_1 is computed by EXACT SVP enumeration, not by LLL's first row.

hdr("E4  predictor study: does lambda_1/dist separate success from failure?")

def exact_lam1(B):
    """Exact lambda_1 via full-block BKZ (= HKZ), whose b_0 is a shortest vector.

    NOT SVP.shortest_vector: the fpylll 0.6.4 wheel hardcodes its build-time path
    for the pruning strategies (/project/local/share/fplll/strategies/default.json)
    and raises FileNotFoundError.  Full-block BKZ needs no strategies file.
    Mutates B.
    """
    BKZ.reduction(B, BKZ.Param(B.nrows))
    return norm([B[0][j] for j in range(B.ncols)])

def instance(curve, m, d_secret, k1_bound, seed):
    """One CVP-EXACT-C instance; returns (ok, svpv_exact, eff, lam_star)."""
    p, b, n, lam, G = curve
    lam_cache[0] = lam
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    dim = 2 * m
    B = reduced_basis(build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound), dim)
    t = cvp_target(sigs, n, k1_bound, k2_bound, True)
    pv = cvp_planted_offset(sigs, n, k1_bound, k2_bound, True)
    v = list(CVP.closest_vector(B, list(t)))
    ok = decode_cvp(v, t, sigs, n, k1_bound, k2_bound, True, d_secret) is not None
    # copy the basis: SVP.shortest_vector mutates its argument
    Bc = IntegerMatrix.from_matrix([[B[i][j] for j in range(dim)]
                                    for i in range(B.nrows)])
    return ok, exact_lam1(Bc) / norm(pv), k1_bound * k2_bound / n, lam_star(lam, n)

def auc(rows, key):
    """Mann-Whitney AUC of `key` as a score for success. 0.5 = no signal."""
    pos = [r[key] for r in rows if r['ok']]
    neg = [r[key] for r in rows if not r['ok']]
    if not pos or not neg:
        return float('nan'), len(pos), len(neg)
    wins = sum((a > b) + 0.5 * (a == b) for a in pos for b in neg)
    return wins / (len(pos) * len(neg)), len(pos), len(neg)

curves17b = search_curves(17, 18, 12)
print(f"  {len(curves17b)} curves x effs {{0.10,0.18,0.30}} x {len(SEEDS)} seeds, m=12")
rows = []
for cur in curves17b:
    n = cur[2]
    k2b = math.isqrt(n) + 1
    for eff in (0.10, 0.18, 0.30):
        k1b = max(2, round(eff * n / k2b))
        for s in SEEDS:
            d_trial = random.Random(s + 7777).randint(1, n - 1)
            r = instance(cur, 12, d_trial, k1b, s)
            if r is None:
                continue
            ok, svpv, e, ls = r
            rows.append({'ok': ok, 'svpv': svpv, 'eff': e, 'lamstar': ls,
                         'n': n, 'K1': k1b})

nsucc = sum(r['ok'] for r in rows)
print(f"  n={len(rows)} instances, {nsucc} recovered ({nsucc/len(rows):.1%})")
print(f"\n  {'predictor':<12}{'AUC':>8}   interpretation")
for key, sign in (('svpv', +1), ('eff', -1), ('lamstar', +1)):
    a, npos, nneg = auc(rows, key)
    a_signed = a if sign > 0 else 1 - a
    print(f"  {key:<12}{a_signed:>8.3f}   ({npos} succ / {nneg} fail)"
          + ("   [higher predicts success]" if sign > 0
             else "   [lower predicts success]"))

cands = sorted({r['svpv'] for r in rows})
best = max(((sum((r['svpv'] >= c) == r['ok'] for r in rows) / len(rows), c)
            for c in cands), key=lambda z: z[0])
maj = max(nsucc, len(rows) - nsucc) / len(rows)
print(f"\n  best svpv threshold: {best[1]:.3f}  ->  accuracy {best[0]:.1%} "
      f"(majority baseline {maj:.1%})")

print(f"\n  {'svpv bin':<14}{'recovered':>12}")
for lo, hi in [(0, .8), (.8, 1.0), (1.0, 1.2), (1.2, 1.5), (1.5, 99)]:
    sub = [r for r in rows if lo <= r['svpv'] < hi]
    if sub:
        print(f"  [{lo:.1f},{hi:>4.1f})   {sum(r['ok'] for r in sub):>5}/{len(sub):<6}")

# The marginal AUCs of svpv and eff are confounded: svpv is a decreasing function
# of eff, and eff takes only 3 values here.  The discriminating test is whether
# svpv still separates WITHIN an eff stratum -- something eff cannot do by
# construction, and the situation E1/E2 exhibit (2557 vs 2677 at eff 0.156/0.157,
# 1/5 vs 5/5).
print(f"\n  within-eff-stratum AUC of svpv (eff is constant inside each stratum):")
print(f"  {'eff':>6}{'AUC':>8}{'succ':>7}{'fail':>7}   svpv range")
for eff in (0.10, 0.18, 0.30):
    sub = [r for r in rows if abs(r['eff'] - eff) < 0.02]
    a, npos, nneg = auc(sub, 'svpv')
    rng = (min(r['svpv'] for r in sub), max(r['svpv'] for r in sub)) if sub else (0, 0)
    print(f"  {eff:>6.2f}{a:>8.3f}{npos:>7}{nneg:>7}   [{rng[0]:.2f}, {rng[1]:.2f}]")

# ---- E5: are the two formulations complementary, and is that lambda-driven? --
#
# E2 shows the two formulations trade places: on 12-bit/2677 (lam*=0.070) the BDD
# reformulation moves the K1 wall from ~4-6 out to ~16-24, while on 12-bit/2557
# (lam*=0.340) it moves the wall backwards from ~12-16 in to ~3-4.  Taking the
# better of the two per curve puts BOTH walls at exactly K1=12 -- i.e. the
# apparent lam-dependence of the wall would be a formulation artifact, not a
# property of the curve.  On two curves that could be coincidence.  E5 tests it
# on the 12 fresh 17-bit curves at a fixed eff inside the transition zone.
#
# Hypothesis H23: sign(KANNAN - CVP) tracks lam*: Kannan wins at large lam*,
# the BDD reformulation wins at small lam*.

hdr("E5  complementarity of the two formulations: is it lam-driven or n-driven?")

curves12 = search_curves(11, 13, 12)
print(f"  curve sets: {len(curves12)} x 12-bit, {len(curves17b)} x 17-bit; m=12, "
      f"{len(SEEDS)} seeds")

def head_to_head(curves, eff):
    """Returns (kannan_wins, cvp_wins, union, total, per-curve rows)."""
    wk = wc = wb = tot = 0
    per = []
    for cur in curves:
        p, b, n, lam, G = cur
        k2b = math.isqrt(n) + 1
        k1b = max(2, round(eff * n / k2b))
        ck = cc = 0
        for s in SEEDS:
            d_trial = random.Random(s + 7777).randint(1, n - 1)
            try:
                a = bool(run_kannan(cur, 12, d_trial, k1b, s))
            except Exception:
                a = False
            try:
                c = bool(run_cvp(cur, 12, d_trial, k1b, s, centred=True,
                                 exact=False)[0])
            except Exception:
                c = False
            ck += a; cc += c; wb += (a or c); tot += 1
        wk += ck; wc += cc
        per.append((lam_star(lam, n), ck, cc))
    return wk, wc, wb, tot, per

EFFS5 = [0.08, 0.12, 0.18, 0.24, 0.30]
for bits, curves in (("12-bit", curves12), ("17-bit", curves17b)):
    print(f"\n  {bits} curves")
    print(f"  {'eff':>6}{'KANNAN':>10}{'CVP-B-C':>10}{'best-of':>10}"
          f"{'union gain':>12}   split by lam* (K/C)")
    for eff in EFFS5:
        wk, wc, wb, tot, per = head_to_head(curves, eff)
        lo = [r for r in per if r[0] < 0.20]
        hi = [r for r in per if r[0] >= 0.20]
        s_lo = (f"lo {sum(r[1] for r in lo)}/{sum(r[2] for r in lo)}"
                if lo else "lo -")
        s_hi = (f"hi {sum(r[1] for r in hi)}/{sum(r[2] for r in hi)}"
                if hi else "hi -")
        print(f"  {eff:>6.2f}{f'{wk}/{tot}':>10}{f'{wc}/{tot}':>10}"
              f"{f'{wb}/{tot}':>10}{wb - max(wk, wc):>12}   {s_lo}  {s_hi}")
        sys.stdout.flush()

print("\n  'union gain' = instances solved by the best-of-two but by neither")
print("  method alone being uniformly better; >0 means the formulations are")
print("  genuinely complementary.  'lo/hi' split tests whether lam* drives which")
print("  formulation wins (lo = lam* < 0.20, hi = lam* >= 0.20).")

print(f"\nelapsed {time.time() - T0:.1f}s")
