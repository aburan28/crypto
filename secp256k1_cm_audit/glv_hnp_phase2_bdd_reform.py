"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
lambda_1 (or at least so that BDD decoding is unique).

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  In the Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` the shortest vector
  is ALWAYS the trivial vector  n * S_D * e_m  (100% of its energy in the
  d-column, |sv[m]|/n = 1.0000 exactly).  It is 2-3x shorter than the planted
  vector on every curve, success and failure alike, so the planted vector is
  never lambda_1 and recovery is a BDD/coset condition rather than SVP.
  That run proposed: "project the lattice along e_m (quotient out the trivial
  n*e_m direction) and solve BDD in the projection, or replace the Kannan
  embedding with an explicit CVP call".

  Falsifier stated on 2026-07-29:
    if sv/pv rises above 1 after the reformulation AND the K1 wall on the
    lam*=0.07 curve (currently K1 ~ 4-6) moves outward, the reformulation is a
    real improvement; if the wall stays at K1 ~ 4-6, the wall is
    information-theoretic and Phase 2 is at its ceiling.

Two independent defects of the 2026-06-15 construction are fixed here, and
they are ablated separately so the cause is attributable:

  (P) PROJECTION.  The d-column is unnecessary.  d is recoverable from any
      single recovered pair (k1_i, k2_i) via
          d = (k1_i + lam*k2_i - A_i) * B_i^{-1}  (mod n),
      so the column can be deleted.  Deleting it removes the trivial vector
      n*S_D*e_m from the lattice entirely (d and d+n now give the SAME lattice
      vector, which is correct since d is only defined mod n) and removes
      E[d^2]*S_D^2 = n^2/3 from the planted vector's squared norm.

  (C) CENTERING.  The 2026-06-15 planted vector has coordinates
      k1_i*S_K1 ~ U[0, n) and k2_i*S_K2 ~ U[0, n), i.e. mean n/2, so
      E[x^2] = n^2/3.  Subtracting the known means K1/2 and K2/2 from the
      Kannan/target row gives coordinates ~ U[-n/2, n/2] with
      E[x^2] = n^2/12.  That is a factor 4 in squared norm, i.e. a factor 2
      in the norm of the vector that must be found.  Free, and independent
      of (P).

Variants measured (same curves, same seeds, same signatures throughout):

  V0    d-col, uncentered, S_KANNAN=n,        readout=dcol   [2026-06-15 control]
  V0c   d-col, uncentered, S_KANNAN=n,        readout=coeff  [same lattice]
  V0b   d-col, uncentered, S_KANNAN=n,        readout=box    [same lattice]
  VC    d-col, CENTERED,   S_KANNAN=n/sqrt12, readout=coeff  [(C)]
  VCn   d-col, CENTERED,   S_KANNAN=n,        readout=coeff  [(C) w/o rescale]
  VP    no d-col, uncentered, S_KANNAN=n,     readout=coeff  [(P)]
  VPC   no d-col, CENTERED, S_KANNAN=n/sqrt12, readout=coeff [(P)+(C)]
  VCVP  no d-col, CENTERED, exact CVP (fpylll enumeration)   [no Kannan]
  VBAB  no d-col, CENTERED, Babai nearest-plane on LLL basis [no Kannan]

Every candidate d is accepted only after a check that stands in for the
attacker's d*G == Q test; d_secret is never used to *construct* a candidate.
See the long NOTE above `coeff_d` for why reading d off the d-column is
strictly weaker than recovering the planted vector, and how the same quantity
is extracted without that column.

RESULTS (2026-08-06 run; full output in glv_hnp_phase2_bdd_reform_output.txt)

  R1  Falsifier branch 1 SATISFIED.  sv/pv rises from 0.42-0.60 (V0) to
      1.000 (VPC) and 1.23-1.30 (VCVP) on 8-bit/199 and 12-bit/2677.
      pv/GH falls from 1.20-1.30 to 0.81-0.94, i.e. after centering the
      planted vector is genuinely shorter than the Gaussian heuristic.

  R2  Falsifier branch 2 SATISFIED — the K1 wall MOVES, a lot.
      m=12, 5 seeds, last K1 with >=3/5:
          12-bit/2557 (lam*=0.34):  V0 -> 12,  VC -> 24
          12-bit/2677 (lam*=0.07):  V0 ->  4,  VC -> 24   (factor 6)
          20-bit      (lam*=0.34):  V0 -> 32,  VC -> 128  (factor 4)
      The lam* dependence of the wall VANISHES after centering: both 12-bit
      curves land on exactly the same wall.  So the "lam* shifts the K1 wall
      by ~3x" observation of 2026-07-29 T4 is an artifact of the *uncentered*
      lattice, not a property of the problem.

  R3  Ablation.  PROJECTION IS A NO-OP: VP reproduces V0 cell-for-cell on
      both curves (30/30 and 17/17 successes over the grid).  CENTERING is
      the entire effect (42 and 40).  The Kannan rescale n -> n/sqrt(12)
      contributes a little on top of centering (VC 42/40 vs VCn 39/38).
      Exact CVP (VCVP 36/36) and Babai (VBAB 39/37) do NOT beat the centered
      Kannan embedding, because Kannan supplies several candidate rows per
      run while a single CVP call supplies one.
      V0 == V0c cell-for-cell (consistency check on the two readouts).
      V0b (demand the planted vector exactly) is far weaker: 12 vs 30 on
      12-bit/2557.

  R5  The remaining wall is INFORMATION-THEORETIC.  Bit-counting gives
      eff = K1*K2/n <= n^(-1/m).  Measured/predicted over m = 4..48 (the
      bound itself moves by 6x over that range, 0.139 -> 0.849):
          12-bit/2557  ratio in [0.74, 1.06]
          12-bit/2677  ratio in [0.74, 1.17]
      The droop at m >= 24 is grid resolution (adjacent K1 grid points are
      eff 0.63 / 0.79 / 0.94 there), not a real deficit.
      => the centered attack saturates the counting bound.  Phase 2 is at
      its ceiling; further lattice engineering cannot move this wall.

Run: python3 glv_hnp_phase2_bdd_reform.py
"""

import math
import random

from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ===========================================================================
# Minimal EC arithmetic -- copied verbatim from
# glv_hnp_phase2_lambda_threshold.py:87-143 so the comparison is exact.
# ===========================================================================

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

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    """Verbatim from glv_hnp_phase2_lambda_threshold.py:230."""
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


# ===========================================================================
# Generalised lattice builder
# ===========================================================================

class Lat:
    """Column layout / scaling metadata for one variant."""
    def __init__(self, m, n, K1, K2, drop_d, center, kannan, kan_scale):
        self.m, self.n, self.K1, self.K2 = m, n, K1, K2
        self.drop_d, self.center, self.kannan = drop_d, center, kannan
        self.S_K1 = n // K1
        self.S_D = 1
        self.S_K2 = max(1, n // K2)
        if kannan:
            self.S_KANNAN = max(1, int(round(n * kan_scale)))
        else:
            self.S_KANNAN = None
        self.c1 = K1 // 2 if center else 0
        self.c2 = K2 // 2 if center else 0
        self.cd = n // 2 if center else 0
        self.dcol = None if drop_d else m
        off = m + (0 if drop_d else 1)
        self.k2col = [off + i for i in range(m)]
        self.kcol = (off + m) if kannan else None
        self.ncols = off + m + (1 if kannan else 0)


def build(sigs, n, lam, K1, K2, drop_d=False, center=False,
          kannan=True, kan_scale=1.0):
    """Rows:
         i        : n*S_K1 * e_i                       (modular reduction)
         d-row    : (B_i*S_K1)_i  [, S_D at dcol]      (the secret multiplier)
         k2-row i : -lam*S_K1 at col i, S_K2 at k2col(i)
         kannan   : ((A_i-c1)*S_K1)_i [, -cd*S_D] , (-c2*S_K2)_i , S_KANNAN
       The planted vector is  kannan + d*d-row + sum_i k2_i*k2row_i
                              - sum_i q_i*modrow_i,
       with coordinates ((k1_i-c1)S_K1, (d-cd)S_D, (k2_i-c2)S_K2, S_KANNAN),
       because A_i + B_i*d = k1_i + lam*k2_i (mod n).
    """
    m = len(sigs)
    L = Lat(m, n, K1, K2, drop_d, center, kannan, kan_scale)
    rows = []
    for i in range(m):
        r = [0] * L.ncols
        r[i] = n * L.S_K1
        rows.append(r)
    r = [0] * L.ncols
    for i in range(m):
        r[i] = sigs[i]['B'] * L.S_K1
    if L.dcol is not None:
        r[L.dcol] = L.S_D
    rows.append(r)
    for i in range(m):
        r = [0] * L.ncols
        r[i] = -(lam % n) * L.S_K1
        r[L.k2col[i]] = L.S_K2
        rows.append(r)
    if kannan:
        r = [0] * L.ncols
        for i in range(m):
            r[i] = (sigs[i]['A'] - L.c1) * L.S_K1
            r[L.k2col[i]] = -L.c2 * L.S_K2
        if L.dcol is not None:
            r[L.dcol] = -L.cd * L.S_D
        r[L.kcol] = L.S_KANNAN
        rows.append(r)
    return rows, L


def target_vector(sigs, L):
    """CVP target t (no Kannan column).  The closest lattice vector w gives
    error e = w - t with coordinates ((k1_i-c1)S_K1, (d-cd)S_D, (k2_i-c2)S_K2).
    """
    t = [0] * L.ncols
    for i in range(L.m):
        t[i] = (L.c1 - sigs[i]['A']) * L.S_K1
        t[L.k2col[i]] = L.c2 * L.S_K2
    if L.dcol is not None:
        t[L.dcol] = L.cd * L.S_D
    return t


def planted_vector(sigs, d_secret, L):
    v = [0] * L.ncols
    for i in range(L.m):
        v[i] = (sigs[i]['k1'] - L.c1) * L.S_K1
        v[L.k2col[i]] = (sigs[i]['k2'] - L.c2) * L.S_K2
    if L.dcol is not None:
        v[L.dcol] = (d_secret - L.cd) * L.S_D
    if L.kannan:
        v[L.kcol] = L.S_KANNAN
    return v


def norm(v):
    return math.sqrt(sum(x * x for x in v))


# ===========================================================================
# Readouts
# ===========================================================================

# ---------------------------------------------------------------------------
# NOTE (found while building this script, 2026-08-06).  The 2026-06-15 readout
# `recover_d` reads d off the d-column.  That column holds the *coefficient*
# c_d of the d-row in the integer combination that produced the row -- it does
# NOT require the row to be the planted vector.  Any lattice vector in the
# Kannan coset whose d-row coefficient is congruent to d mod n yields d, even
# when its k1/k2 parts are a completely different (out-of-box) solution of the
# same congruence.  Empirically that is exactly what happens: on
# 12-bit/2557, K1=4, m=8, seed 42 the winning LLL row has
#     k1 = [0,1,3,0,0,1,0,3]   (true [0,1,0,0,0,1,0,3])
#     k2 = [1,14,-7,37,1,14,35,14]  (true [1,14,43,37,1,14,35,14])
# yet its d-column is 65 = d exactly.  So the recovery condition is much
# weaker than "find the planted vector".
#
# c_d is recoverable the same way from the congruence + k2 columns alone, so
# the d-column is not needed for it:
#     w_i/S_K1 = A_i - c1 + B_i*c_d - lam*c_i - q_i*n
#  => c_d = (w_i/S_K1 + c1 - A_i + lam*c_i) * B_i^{-1}  (mod n),
#     c_i  = w[k2col i]/S_K2 + c2.
# The same formula holds for a CVP error vector e = w - t (the -A_i + c1 term
# comes from the target instead of the Kannan row).  All variants below
# therefore use this 'coeff' readout, which makes the comparison fair; the
# strictly weaker 'box' readout (demand k1_i, k2_i inside their boxes) is
# kept as a separate variant to quantify the gap.
# ---------------------------------------------------------------------------

def coeff_d(cong, k2v, sigs, n, lam, L):
    """d-row coefficient mod n, from congruence coords `cong` (already
    sigma-normalised) and k2 coords `k2v`.  Returns None if the coords are not
    those of a lattice vector."""
    cs = []
    for i in range(L.m):
        if k2v[i] % L.S_K2 or cong[i] % L.S_K1:
            return None
        cs.append(k2v[i] // L.S_K2 + L.c2)
    for i in range(L.m):
        B = sigs[i]['B'] % n
        if math.gcd(B, n) != 1:
            continue
        d = ((cong[i] // L.S_K1) + L.c1 - sigs[i]['A'] + lam * cs[i]) \
            * modinv(B, n) % n
        return d if d != 0 else None
    return None


def box_d(cong, k2v, sigs, n, lam, L):
    """Strictly weaker readout: demand the recovered k1_i, k2_i lie in their
    boxes and that a single d fits every signature."""
    k1s, k2s = [], []
    for i in range(L.m):
        if k2v[i] % L.S_K2 or cong[i] % L.S_K1:
            return None
        k1s.append(cong[i] // L.S_K1 + L.c1)
        k2s.append(k2v[i] // L.S_K2 + L.c2)
    for i in range(L.m):
        if not (0 <= k1s[i] < L.K1) or not (0 <= k2s[i] < L.K2):
            return None
    for i in range(L.m):
        B = sigs[i]['B'] % n
        if math.gcd(B, n) != 1:
            continue
        d = (k1s[i] + lam * k2s[i] - sigs[i]['A']) * modinv(B, n) % n
        if all((sigs[j]['A'] + sigs[j]['B'] * d - k1s[j] - lam * k2s[j]) % n == 0
               for j in range(L.m)) and d != 0:
            return d
    return None


def read_row(row, sigs, n, lam, L, readout, d_secret=None):
    """Candidate d from one LLL row of a Kannan-embedded basis."""
    last = row[L.kcol]
    if abs(last) != L.S_KANNAN:
        return None
    sg = 1 if last > 0 else -1
    cong = [sg * row[i] for i in range(L.m)]
    k2v = [sg * row[L.k2col[i]] for i in range(L.m)]
    if readout == 'dcol':
        # 2026-06-15 semantics, verbatim
        d = (sg * row[L.dcol] + L.cd) % n
        if d == 0:
            return None
        return d if (d_secret is None or d == d_secret) else None
    if readout == 'box':
        return box_d(cong, k2v, sigs, n, lam, L)
    return coeff_d(cong, k2v, sigs, n, lam, L)


def read_error(e, sigs, n, lam, L, readout):
    """Candidate d from a CVP error vector e = w - t."""
    cong = [e[i] for i in range(L.m)]
    k2v = [e[L.k2col[i]] for i in range(L.m)]
    if readout == 'box':
        return box_d(cong, k2v, sigs, n, lam, L)
    return coeff_d(cong, k2v, sigs, n, lam, L)


# ===========================================================================
# Reduction / CVP helpers
# ===========================================================================

def lll_rows(rows, beta=0):
    A = IntegerMatrix.from_matrix(rows)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    out = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    return [r for r in out if any(r)]


def gram_schmidt(B):
    """Float GS of a list of integer rows.  Dimensions here are <= ~30 and
    entries <= ~n^2 ~ 1e7, well inside f64."""
    star, mu_den = [], []
    for b in B:
        v = [float(x) for x in b]
        for s, d in zip(star, mu_den):
            if d == 0.0:
                continue
            c = sum(vi * si for vi, si in zip(b, s)) / d
            v = [vi - c * si for vi, si in zip(v, s)]
        star.append(v)
        mu_den.append(sum(x * x for x in v))
    return star, mu_den


def babai_nearest_plane(B, t):
    star, den = gram_schmidt(B)
    b = [float(x) for x in t]
    coeffs = [0] * len(B)
    for i in range(len(B) - 1, -1, -1):
        if den[i] == 0.0:
            continue
        c = sum(bi * si for bi, si in zip(b, star[i])) / den[i]
        c = int(math.floor(c + 0.5))
        coeffs[i] = c
        b = [bi - c * float(Bij) for bi, Bij in zip(b, B[i])]
    w = [0] * len(t)
    for i, c in enumerate(coeffs):
        if c:
            for j in range(len(t)):
                w[j] += c * B[i][j]
    return w


def cvp_exact(B, t):
    A = IntegerMatrix.from_matrix(B)
    LLL.reduction(A)
    A = IntegerMatrix.from_matrix([[A[i][j] for j in range(A.ncols)]
                                   for i in range(A.nrows)
                                   if any(A[i][j] for j in range(A.ncols))])
    return list(CVP.closest_vector(A, tuple(int(x) for x in t)))


def log_det(B):
    _, den = gram_schmidt(B)
    return 0.5 * sum(math.log(d) for d in den if d > 0)


def gaussian_heuristic(B):
    k = sum(1 for _ in B)
    ld = log_det(B)
    return math.exp(ld / k) * math.sqrt(k / (2 * math.pi * math.e))


# ===========================================================================
# One trial of one variant
# ===========================================================================

RT12 = 1.0 / math.sqrt(12)

VARIANTS = {
    #        drop_d center kannan kan_scale readout   solver
    'V0':   (False, False, True,  1.0,  'dcol',  'lll'),   # 2026-06-15 control
    'V0c':  (False, False, True,  1.0,  'coeff', 'lll'),   # same lattice, coeff readout
    'V0b':  (False, False, True,  1.0,  'box',   'lll'),   # same lattice, box readout
    'VC':   (False, True,  True,  RT12, 'coeff', 'lll'),   # (C) centering + Kannan rescale
    'VCn':  (False, True,  True,  1.0,  'coeff', 'lll'),   # (C) centering only
    'VP':   (True,  False, True,  1.0,  'coeff', 'lll'),   # (P) projection only
    'VPC':  (True,  True,  True,  RT12, 'coeff', 'lll'),   # (P)+(C)
    'VCVP': (True,  True,  False, 0.0,  'coeff', 'cvp'),   # exact CVP, no Kannan
    'VBAB': (True,  True,  False, 0.0,  'coeff', 'babai'), # Babai nearest-plane
}


def trial(variant, sigs, d_secret, n, lam, K1, K2, beta=0, want_stats=False):
    drop_d, center, kannan, kscale, readout, solver = VARIANTS[variant]
    rows, L = build(sigs, n, lam, K1, K2, drop_d, center, kannan, kscale)
    red = lll_rows(rows, beta)
    stats = {}
    cands = 0
    if kannan:
        d_found = None
        for r in red:
            cand = read_row(r, sigs, n, lam, L, readout, d_secret)
            if cand is None:
                continue
            cands += 1
            if cand == d_secret:          # stands in for the d*G == Q check
                d_found = cand
                break
    else:
        t = target_vector(sigs, L)
        if solver == 'cvp':
            try:
                w = cvp_exact(red, t)
            except Exception:
                w = babai_nearest_plane(red, t)
        else:
            w = babai_nearest_plane(red, t)
        e = [int(wi) - int(ti) for wi, ti in zip(w, t)]
        cand = read_error(e, sigs, n, lam, L, readout)
        cands = int(cand is not None)
        d_found = cand if cand == d_secret else None
    if want_stats:
        pv = planted_vector(sigs, d_secret, L)
        sv = min(norm(r) for r in red)
        stats = {'sv': sv, 'pv': norm(pv), 'ratio': sv / norm(pv),
                 'dim': len(red), 'gh': gaussian_heuristic(red),
                 'cands': cands}
        if not kannan:
            stats['err'] = norm(e)
            stats['err_is_planted'] = (e == pv)
    return (d_found is not None), stats


def rate(variant, curve, m, K1, seeds, beta=0):
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    wins = 0
    for sd in seeds:
        d_secret = random.Random(sd + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, sd)
        if len(sigs) < m:
            continue
        ok, _ = trial(variant, sigs, d_secret, n, lam, K1, K2, beta)
        wins += int(ok)
    return wins, len(seeds)


# ===========================================================================
# Curves
# ===========================================================================

print("=" * 78)
print("Thread 23 — reformulating the Phase-2 GLV lattice (projection + centering)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G), k1, m))
    print(f"  {label:<18} p={p:<6} n={n:<6} lam={lam:<6} lam*={lam_star(lam,n):.4f}")


# ===========================================================================
# EXPERIMENT R1 — structural: does the reformulation make the planted vector
#                 the shortest vector?  (the stated falsifier)
# ===========================================================================

print("\n" + "-" * 78)
print("EXP R1: sv/pv after reformulation  (falsifier: must rise above 1)")
print("-" * 78)
print("sv = shortest row after LLL, pv = ||planted vector||.")
print("2026-07-29 (T5) measured sv/pv in [0.34, 0.61] for V0 on these curves.\n")

print(f"{'curve':<18} {'K1':>3} {'variant':<5} {'dim':>4} {'sv':>11} {'pv':>11} "
      f"{'sv/pv':>7} {'pv/GH':>7}")
r1_rows = []
for label, curve, k1, m in curves:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    d_secret = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, K2, 42)
    for v in ('V0', 'VC', 'VP', 'VPC', 'VCVP'):
        ok, st = trial(v, sigs, d_secret, n, lam, k1, K2, want_stats=True)
        print(f"{label:<18} {k1:>3} {v:<5} {st['dim']:>4} {st['sv']:>11.1f} "
              f"{st['pv']:>11.1f} {st['ratio']:>7.3f} {st['pv']/st['gh']:>7.3f}")
        r1_rows.append((label, v, st['ratio'], st['pv'] / st['gh'], ok))
    print()

print("Reading: sv/pv > 1 means the planted vector IS the shortest vector, so")
print("recovery becomes an SVP question that LLL/BKZ directly attack.")
print("pv/GH < 1 means the planted vector is shorter than the Gaussian")
print("heuristic, i.e. it is a genuine unusually-short vector.")


# ===========================================================================
# EXPERIMENT R2 — the K1 wall (the second half of the falsifier)
# ===========================================================================

print("\n" + "-" * 78)
print("EXP R2: the K1 wall, m=12, 5 seeds — does the reformulation move it?")
print("-" * 78)
print("2026-07-29 T4 walls (V0): lam*=0.34 curve -> K1 ~ 12-16;")
print("                          lam*=0.07 curve -> K1 ~ 4-6.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
M_R2 = 12
R2_VARIANTS = ['V0', 'V0c', 'V0b', 'VC', 'VCn', 'VP', 'VPC', 'VCVP', 'VBAB']

r2 = {}
for label, curve, _k1, _m in curves[1:]:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"{label}  (n={n}, lam*={lam_star(lam,n):.4f}, K2={K2}, m={M_R2})")
    print(f"{'variant':<7} " + " ".join(f"{('K1='+str(k)):>7}" for k in K1_GRID))
    for v in R2_VARIANTS:
        cells = []
        for K1 in K1_GRID:
            w, t = rate(v, curve, M_R2, K1, SEEDS)
            cells.append(f"{w}/{t}")
            r2[(label, v, K1)] = (w, t)
        print(f"{v:<7} " + " ".join(f"{c:>7}" for c in cells))
    print(f"{'eff':<7} " + " ".join(f"{(K1*K2/n):>7.2f}" for K1 in K1_GRID))
    print()


def wall_of(label, v):
    """Largest K1 on the grid with >=3/5, and the first K1 with 0/5."""
    last_ok, first_dead = None, None
    for K1 in K1_GRID:
        w, t = r2.get((label, v), (0, 0)) if False else r2[(label, v, K1)]
        if w * 2 >= t:
            last_ok = K1
        if w == 0 and first_dead is None:
            first_dead = K1
    return last_ok, first_dead

print(f"{'curve':<18} {'variant':<7} {'last K1 >=3/5':>14} {'first K1 = 0/5':>15}")
for label, curve, _k1, _m in curves[1:]:
    for v in R2_VARIANTS:
        lo, fd = wall_of(label, v)
        print(f"{label:<18} {v:<7} {str(lo):>14} {str(fd):>15}")


# ===========================================================================
# EXPERIMENT R3 — ablation summary: which change carries the improvement?
# ===========================================================================

print("\n" + "-" * 78)
print("EXP R3: ablation — total successes over the whole K1 grid (out of 65)")
print("-" * 78)
print(f"{'curve':<18} " + " ".join(f"{v:>6}" for v in R2_VARIANTS))
for label, curve, _k1, _m in curves[1:]:
    tot = []
    for v in R2_VARIANTS:
        s = sum(r2[(label, v, K1)][0] for K1 in K1_GRID)
        tot.append(str(s))
    print(f"{label:<18} " + " ".join(f"{t:>6}" for t in tot))
print("\nV0/V0c/V0b are the SAME lattice with three readouts (dcol / coeff / box):")
print("  V0 == V0c is a consistency check; V0b quantifies how much weaker the")
print("  'recover the planted vector exactly' condition is.")
print("V0 -> VC/VCn isolates CENTERING (with / without Kannan rescale).")
print("V0 -> VP isolates PROJECTION (dropping the d-column).  VPC = both.")
print("VCVP / VBAB replace the Kannan embedding by exact CVP / Babai.")


# ===========================================================================
# EXPERIMENT R4 — does the improvement survive to larger n?
# ===========================================================================

print("\n" + "-" * 78)
print("EXP R4: 17-bit and 20-bit curves, m=12, 5 seeds, V0 vs VC vs VPC vs VCVP")
print("-" * 78)

def build_curve(p, n):
    """Find a sextic twist y^2 = x^3 + b over F_p whose order is n."""
    for b in range(1, 60):
        if b % p == 0:
            continue
        G = find_generator(p, b, n)
        if G is not None:
            return b, G
    return None, None


BIG = [
    ("17-bit", 65713, 65269, 442),
    ("20-bit", 524347, 523969, 177902),
]
big_curves = []
for label, p, n, lam in BIG:
    if (lam * lam + lam + 1) % n != 0:
        print(f"  SKIP {label}: bad lam")
        continue
    b, G = build_curve(p, n)
    if G is None:
        print(f"  SKIP {label}: no b in 1..59 gives order n={n}")
        continue
    print(f"  {label}: b={b}")
    big_curves.append((label, (p, b, n, lam, G)))

K1_BIG = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
for label, curve in big_curves:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"\n{label}  n={n}  lam*={lam_star(lam,n):.4f}  K2={K2}  m=12")
    print(f"{'variant':<7} " + " ".join(f"{('K1='+str(k)):>8}" for k in K1_BIG))
    for v in ('V0', 'VC', 'VPC', 'VCVP'):
        cells = []
        for K1 in K1_BIG:
            w, t = rate(v, curve, 12, K1, SEEDS)
            cells.append(f"{w}/{t}")
        print(f"{v:<7} " + " ".join(f"{c:>8}" for c in cells))
    print(f"{'eff':<7} " + " ".join(f"{(K1*K2/n):>8.2f}" for K1 in K1_BIG))


# ===========================================================================
# EXPERIMENT R5 — is the reformulated wall the information-theoretic one?
# ===========================================================================
#
# Unknowns:   d  plus  m values k1_i in [0,K1)  plus  m values k2_i in [0,K2).
# Constraints: the m congruences A_i + B_i*d = k1_i + lam*k2_i (mod n).
# Counting bits,
#       log2(n) + m*log2(K1*K2)  <=  m*log2(n)
#   <=> log2(K1*K2)              <=  ((m-1)/m) * log2(n)
#   <=> eff := K1*K2/n           <=  n^(-1/m).
# If the centered attack sits exactly on this line, no further lattice
# engineering can help: the instance below the line is under-determined.
# ===========================================================================

print("\n" + "-" * 78)
print("EXP R5: measured wall vs the information-theoretic bound eff <= n^(-1/m)")
print("-" * 78)

K1_R5 = [2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 40, 48, 56, 64]
M_R5 = [4, 6, 8, 12, 16, 24, 32, 48]

for label, curve, _k1, _m in curves[1:]:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"\n{label}  n={n}  lam*={lam_star(lam,n):.4f}  K2={K2}   variant VC")
    print(f"{'m':>4} {'eff_max last >=3/5':>19} {'eff first 0/5':>14} "
          f"{'predicted n^(-1/m)':>19} {'ratio':>7}")
    for m in M_R5:
        last_ok, first_dead = None, None
        for K1 in K1_R5:
            w, t = rate('VC', curve, m, K1, SEEDS)
            eff = K1 * K2 / n
            if w * 2 >= t:
                last_ok = eff
            if w == 0 and first_dead is None:
                first_dead = eff
        pred = n ** (-1.0 / m)
        ratio = (last_ok / pred) if (last_ok and pred) else float('nan')
        print(f"{m:>4} {str(round(last_ok,3)) if last_ok else '-':>19} "
              f"{str(round(first_dead,3)) if first_dead else '>1.26':>14} "
              f"{pred:>19.3f} {ratio:>7.2f}")

print("\nratio ~ 1 means the centered attack saturates the counting bound, i.e.")
print("the remaining wall is information-theoretic and no further lattice")
print("engineering can move it.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
