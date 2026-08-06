"""
GLV-HNP Phase 2 / Thread 23 — reformulate the lattice so the target is lambda_1,
and decide whether the K1 wall is algorithmic or information-theoretic.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 Kannan lattice always contains the trivial vector  n*S_D*e_m  (obtained
  as  n*row_m - sum_i B_i*row_i ), of norm n*S_D, while the planted vector has norm
  ~ n*sqrt(2m/3 + 4/3) > n*S_D.  So the planted vector is NEVER lambda_1, and
  recovery is a BDD/coset condition, not an SVP condition.  The 2026-07-29 run
  proposed Thread 23: quotient out the trivial direction and re-measure the K1 wall.

This script runs four experiments:

  U1  Confirm the trivial vector algebraically + measure sv/pv before and after
      projecting along e_m (the trivial direction).
  U2  Three decoders on a K1 grid, on the lam*=0.34 curve and the lam*=0.07 curve:
        (a) KANNAN   -- baseline: full lattice, LLL, scan rows for |last|=S_KANNAN
        (b) PROJ     -- drop column m, LLL with transformation tracking, read d
                        off the transformation row (d is recovered mod n, which is
                        exactly the ambiguity the trivial vector represents)
        (c) BABAI    -- no Kannan row: CVP by Babai nearest-plane on the LLL-reduced
                        (2m+1)-dim lattice, target t = -(A_i*S_K1, 0, ..., 0)
  U3  Exact-CVP oracle: is the planted vector actually THE closest vector to t?
      CAUTION (measured, see log 2026-08-06): d_ex < d_pl does NOT imply failure.
      When eff = K1*K2/n grows, the GLV decomposition k = k1 + lam*k2 stops being
      unique, so a strictly shorter lattice point can carry the SAME d.  The only
      sound success criterion is the d-coordinate mod n, reported separately as
      'CVP=d?'.  Compare that column, not the distance ratio.
  U4  Gaussian-heuristic prediction of the wall:  ||error|| / GH(L0)  as a function
      of eff = K1*K2/n, compared against the measured walls.

Helpers (EC arithmetic, curve search, lattice construction) are copied verbatim from
glv_hnp_phase2_lambda_threshold.py so the comparison to 2026-07-29 is exact.

Run: python3 glv_hnp_phase2_thread23_bdd.py
"""

import math
import random
import time

from fpylll import IntegerMatrix, LLL, GSO, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
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
        t = t * c % p
        r = r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and ec_mul(P, 1, p) is not None:
            return P
    return None

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    """Symmetric size of lam: min(lam, n-lam)/n in [0, 0.5]."""
    return min(lam % n, n - (lam % n)) / n

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

# ---------------------------------------------------------------------------
# Decoder (a): KANNAN baseline -- 2026-07-29 code path, unchanged
# ---------------------------------------------------------------------------

def decode_kannan(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 2
    S_KANNAN = scales(n, k1_bound, k2_bound)[3]
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[i][m]) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Decoder (b): PROJ -- quotient out the trivial direction e_m
#
# The trivial vector is n*S_D*e_m = n*row_m - sum_i B_i*row_i.  Deleting column m
# is exactly the orthogonal projection along e_m, which sends the trivial vector to
# 0.  The projected generating set is therefore rank-deficient by 1 (LLL puts the
# zero row first).  d is read off the transformation matrix: in the planted
# combination the coefficient of row_m is d, determined mod n -- which is precisely
# the ambiguity the trivial vector encodes, and d is only defined mod n anyway.
# ---------------------------------------------------------------------------

def decode_proj(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 2
    S_KANNAN = scales(n, k1_bound, k2_bound)[3]
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    Mp = [row[:m] + row[m + 1:] for row in M]      # delete column m
    A = IntegerMatrix.from_matrix(Mp)
    U = IntegerMatrix(A.nrows, A.nrows)
    LLL.reduction(A, U)
    last_col = dim - 2                             # last column after deletion
    for i in range(A.nrows):
        v = A[i][last_col]
        if abs(v) != S_KANNAN: continue
        sign = 1 if v > 0 else -1
        d_cand = (sign * U[i][m]) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Decoder (c): BABAI -- drop the Kannan row/column, solve CVP directly.
#
# L0 = rows 0..2m of M restricted to columns 0..2m (modular rows, d-row, lambda-rows).
# planted' = (k1_i*S_K1, d*S_D, k2_i*S_K2) = a + w  with a = (A_i*S_K1, 0, ..., 0)
# and w in L0.  So w is the lattice vector closest to t = -a, and its coordinate m
# is d*S_D.  The trivial vector n*e_m still lies in L0, but for a *coset* problem it
# is harmless: it shifts d by n, and d is only defined mod n.
# ---------------------------------------------------------------------------

def build_L0_and_target(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    d0 = 2 * m + 1
    L0 = [[0] * d0 for _ in range(d0)]
    for i in range(m):
        L0[i][i] = n * S_K1
    for i in range(m):
        L0[m][i] = sigs[i]['B'] * S_K1
    L0[m][m] = S_D
    for i in range(m):
        L0[m + 1 + i][i] = -lam * S_K1
        L0[m + 1 + i][m + 1 + i] = S_K2
    t = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * (m + 1)
    return L0, t

def planted_error(sigs, d_secret, n, k1_bound, k2_bound):
    """planted' = w - t, the BDD error vector (dimension 2m+1)."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    e = [0] * (2 * m + 1)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1
        e[m + 1 + i] = sigs[i]['k2'] * S_K2
    e[m] = d_secret * S_D
    return e

def decode_babai(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    L0, t = build_L0_and_target(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(L0)
    LLL.reduction(A)
    G = GSO.Mat(A)
    G.update_gso()
    try:
        coeffs = G.babai(list(t))
    except Exception:
        return False
    w = [0] * (2 * m + 1)
    for i, c in enumerate(coeffs):
        if c == 0: continue
        for j in range(2 * m + 1):
            w[j] += c * A[i][j]
    d_cand = w[m] % n                              # S_D = 1
    return d_cand == d_secret

# ---------------------------------------------------------------------------
# Decoder (d): CVPEX -- exact CVP (fplll enumeration) on the same L0/t as BABAI.
# This is the OPTIMAL decoder for this lattice formulation, so it separates an
# algorithmic wall (LLL/Babai too weak) from an information-theoretic one.
# Success criterion is d mod n, so a closest vector differing from the planted one
# by a multiple of the trivial vector n*e_m still counts as a recovery.
# ---------------------------------------------------------------------------

def decode_cvp_exact(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    L0, t = build_L0_and_target(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(L0)
    LLL.reduction(A)
    try:
        cv = CVP.closest_vector(A, tuple(int(x) for x in t))
    except Exception:
        return False
    return (cv[m] % n) == d_secret

# ---------------------------------------------------------------------------
# U3 helper: exact CVP oracle
# ---------------------------------------------------------------------------

def exact_cvp_verdict(sigs, n, lam, k1_bound, k2_bound, d_secret, budget_s=60.0):
    """Return (status, dist_exact, dist_planted, d_ok).

    status: 'ok'      -- exact CVP completed
            'timeout' -- exceeded budget / enumeration failed
    d_ok:   whether the exact closest vector yields the true d.
    If dist_exact < dist_planted the planted vector is not the closest vector and NO
    CVP-based decoder can recover d: the wall is information-theoretic.
    """
    m = len(sigs)
    L0, t = build_L0_and_target(sigs, n, lam, k1_bound, k2_bound)
    err = planted_error(sigs, d_secret, n, k1_bound, k2_bound)
    dist_planted = norm(err)
    A = IntegerMatrix.from_matrix(L0)
    LLL.reduction(A)
    t0 = time.time()
    try:
        cv = CVP.closest_vector(A, tuple(int(x) for x in t))
    except Exception:
        return ('timeout', None, dist_planted, False)
    if time.time() - t0 > budget_s:
        pass  # completed anyway; keep the result
    diff = [cv[j] - t[j] for j in range(2 * m + 1)]
    dist_exact = norm(diff)
    d_ok = (cv[m] % n) == d_secret
    return ('ok', dist_exact, dist_planted, d_ok)

# ---------------------------------------------------------------------------
# U4 helper: Gaussian-heuristic prediction
# ---------------------------------------------------------------------------

def gh_ratio(m, n, k1_bound, k2_bound):
    """||E[error]|| / GH(L0) for the (2m+1)-dim BDD lattice L0.

    det(L0) = (n*S_K1)^m * S_D * S_K2^m,  dim = 2m+1
    E[||error||^2] = m*(K1*S_K1)^2/3 + n^2/3 + m*(K2*S_K2)^2/3
    """
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    log_det = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2)
    gh = math.exp(log_det / dim) * math.sqrt(dim / (2 * math.pi * math.e))
    e2 = (m * (k1_bound * S_K1) ** 2 / 3.0
          + (n * S_D) ** 2 / 3.0
          + m * (k2_bound * S_K2) ** 2 / 3.0)
    return math.sqrt(e2) / gh

# ---------------------------------------------------------------------------
# Curves
# ---------------------------------------------------------------------------

HIST = [
    # label,             p,    b, n,    lam,  K1_hist, m_hist
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",      2677, 2, 2647, 185,  8,  10),
]

CURVES = []
for label, p, b, n, lam, k1h, mh in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    CURVES.append({'label': label, 'p': p, 'b': b, 'n': n, 'lam': lam, 'G': G,
                   'K1': k1h, 'm': mh, 'lstar': lam_star(lam, n)})

SEEDS = [42, 1234, 9999, 555, 31337]

print("=" * 78)
print("Thread 23 — is the Phase-2 K1 wall algorithmic or information-theoretic?")
print("=" * 78)
for c in CURVES:
    k2 = math.isqrt(c['n']) + 1
    print(f"  {c['label']:<14} p={c['p']:<5} n={c['n']:<5} lam={c['lam']:<5} "
          f"lam*={c['lstar']:.4f}  K2=sqrt(n)+1={k2}")

# ===========================================================================
# U1 — confirm the trivial vector; sv/pv before and after projection
# ===========================================================================
print("\n" + "-" * 78)
print("U1: trivial vector n*S_D*e_m, and sv/pv before vs after projecting along e_m")
print("-" * 78)
print("sv  = ||shortest row after LLL||,  pv = ||planted vector||")
print("A reformulation 'works' in the Thread-23 sense iff sv/pv >= 1 (planted is lambda_1).")
print()
print(f"{'curve':<14} {'K1':>3} {'m':>3} | {'sv/pv KANNAN':>13} {'sv/pv PROJ':>11} "
      f"| {'trivial in L':>12} {'|sv[m]|/n':>10}")

for c in CURVES:
    n, lam, G, p = c['n'], c['lam'], c['G'], c['p']
    k1, m = c['K1'], c['m']
    k2 = math.isqrt(n) + 1
    d_secret = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2, seed=42)
    if len(sigs) < m:
        print(f"{c['label']:<14} — signature generation failed"); continue
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    M = build_glv_lattice(sigs, n, lam, k1, k2)
    pv = planted_vector(sigs, d_secret, n, k1, k2)
    pv_norm = norm(pv)

    # explicit check that n*row_m - sum_i B_i*row_i == n*S_D*e_m
    dim = 2 * m + 2
    comb = [n * M[m][j] - sum(sigs[i]['B'] * M[i][j] for i in range(m)) for j in range(dim)]
    expect = [0] * dim; expect[m] = n * S_D
    trivial_ok = (comb == expect)

    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    nz = [r for r in rows if any(r)]
    sv_row = min(nz, key=norm)
    sv_kan = norm(sv_row)
    frac_d = abs(sv_row[m]) / n

    Mp = [row[:m] + row[m + 1:] for row in M]
    Ap = IntegerMatrix.from_matrix(Mp)
    LLL.reduction(Ap)
    rowsp = [[Ap[i][j] for j in range(dim - 1)] for i in range(Ap.nrows)]
    nzp = [r for r in rowsp if any(r)]
    sv_proj = norm(min(nzp, key=norm))
    pv_proj = math.sqrt(max(0.0, pv_norm ** 2 - float(pv[m]) ** 2))

    print(f"{c['label']:<14} {k1:>3} {m:>3} | {sv_kan / pv_norm:>13.3f} "
          f"{sv_proj / pv_proj:>11.3f} | {str(trivial_ok):>12} {frac_d:>10.4f}")

# ===========================================================================
# U2 — three decoders on the K1 grid
# ===========================================================================
print("\n" + "-" * 78)
print("U2: K1 wall, three decoders.  5 seeds each; entry = #recovered/5")
print("-" * 78)
print("KANNAN = 2026-07-29 baseline | PROJ = quotient e_m | BABAI = nearest-plane CVP")
print("CVPEX  = exact CVP (optimal decoder for this formulation)")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_FIXED = 12
DECODERS = [('KANNAN', decode_kannan), ('PROJ', decode_proj),
            ('BABAI', decode_babai), ('CVPEX', decode_cvp_exact)]

u2_table = {}
for c in CURVES[1:]:                                  # the two 12-bit curves
    n, lam, G, p = c['n'], c['lam'], c['G'], c['p']
    k2 = math.isqrt(n) + 1
    print(f"\ncurve {c['label']}  (lam*={c['lstar']:.4f}, K2={k2}, m={M_FIXED})")
    hdr = f"  {'decoder':<8}" + "".join(f"{('K1=' + str(k)):>8}" for k in K1_GRID)
    print(hdr)
    for dname, dfn in DECODERS:
        cells = []
        for k1 in K1_GRID:
            wins = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_trial, M_FIXED, n, lam, p, k1, k2, seed)
                if len(sigs) < M_FIXED: continue
                if dfn(sigs, n, lam, k1, k2, d_trial):
                    wins += 1
            cells.append(wins)
            u2_table[(c['label'], dname, k1)] = wins
        print(f"  {dname:<8}" + "".join(f"{str(v) + '/5':>8}" for v in cells))

# ===========================================================================
# U3 — exact-CVP oracle: is the planted vector THE closest vector?
# ===========================================================================
print("\n" + "-" * 78)
print("U3: exact CVP oracle.  d_ex/d_pl < 1 => planted is NOT the closest vector.")
print("    But that does NOT mean failure: once eff=K1*K2/n is large the GLV")
print("    decomposition k = k1 + lam*k2 is no longer unique, so a strictly shorter")
print("    lattice point can carry the SAME d.  Read the 'CVP=d?' column, not the ratio.")
print("-" * 78)

M_CVP = 6                                             # dim(L0) = 13, exact CVP feasible
print(f"m={M_CVP} (dim L0 = {2 * M_CVP + 1}), seed=42\n")
print(f"{'curve':<14} {'K1':>3} {'eff':>7} {'GH ratio':>9} | {'d_pl':>11} {'d_ex':>11} "
      f"{'d_ex/d_pl':>10} {'CVP=d?':>7} {'Babai=d?':>9}")

for c in CURVES[1:]:
    n, lam, G, p = c['n'], c['lam'], c['G'], c['p']
    k2 = math.isqrt(n) + 1
    for k1 in K1_GRID:
        d_trial = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_trial, M_CVP, n, lam, p, k1, k2, 42)
        if len(sigs) < M_CVP: continue
        status, d_ex, d_pl, d_ok = exact_cvp_verdict(sigs, n, lam, k1, k2, d_trial)
        bab = decode_babai(sigs, n, lam, k1, k2, d_trial)
        eff = k1 * k2 / n
        ghr = gh_ratio(M_CVP, n, k1, k2)
        if status != 'ok':
            print(f"{c['label']:<14} {k1:>3} {eff:>7.4f} {ghr:>9.3f} | "
                  f"{d_pl:>11.1f} {'TIMEOUT':>11} {'—':>10} {'—':>7} {str(bab):>9}")
            continue
        print(f"{c['label']:<14} {k1:>3} {eff:>7.4f} {ghr:>9.3f} | "
              f"{d_pl:>11.1f} {d_ex:>11.1f} {d_ex / d_pl:>10.4f} "
              f"{str(d_ok):>7} {str(bab):>9}")

# ===========================================================================
# U4 — Gaussian-heuristic wall prediction vs measurement
# ===========================================================================
print("\n" + "-" * 78)
print("U4: Gaussian-heuristic prediction of the wall (unique decoding needs ratio<1)")
print("-" * 78)
print("ratio(K1) = E||error|| / GH(L0);  asymptotically ratio -> sqrt(2*pi*e/3)*sqrt(eff)")
print(f"  => predicted critical eff = 3/(2*pi*e) = {3 / (2 * math.pi * math.e):.4f}")
print()
for c in CURVES[1:]:
    n = c['n']; k2 = math.isqrt(n) + 1
    print(f"curve {c['label']} (lam*={c['lstar']:.4f}), m={M_FIXED}:")
    print(f"  {'K1':>4} {'eff':>8} {'ratio':>8}   KANNAN/PROJ/BABAI/CVPEX")
    for k1 in K1_GRID:
        r = gh_ratio(M_FIXED, n, k1, k2)
        got = "/".join(str(u2_table.get((c['label'], d, k1), '-')) for d, _ in DECODERS)
        print(f"  {k1:>4} {k1 * k2 / n:>8.4f} {r:>8.3f}   {got}")
    k1c = 3.0 / (2 * math.pi * math.e) * n / k2
    print(f"  predicted critical K1 = eff_crit*n/K2 = {k1c:.2f}")

# ===========================================================================
# U5 — does more data rescue the optimal decoder past the wall?
# ===========================================================================
print("\n" + "-" * 78)
print("U5: m-sweep at and past the wall, LLL/Kannan vs the optimal decoder (CVPEX)")
print("-" * 78)
print("If CVPEX keeps winning where KANNAN dies, the wall is ALGORITHMIC.")
print("If both die at the same K1 regardless of m, it is INFORMATION-THEORETIC.\n")

M_SWEEP = [6, 8, 12, 16, 24, 32]
for c in CURVES[1:]:
    n, lam, G, p = c['n'], c['lam'], c['G'], c['p']
    k2 = math.isqrt(n) + 1
    for k1 in [6, 8, 12, 16, 24]:
        print(f"curve {c['label']} (lam*={c['lstar']:.4f}), K1={k1}, "
              f"eff={k1 * k2 / n:.4f}")
        print(f"  {'decoder':<8}" + "".join(f"{('m=' + str(mm)):>8}" for mm in M_SWEEP))
        for dname, dfn in [('KANNAN', decode_kannan), ('CVPEX', decode_cvp_exact)]:
            cells = []
            for mm in M_SWEEP:
                wins = 0
                for seed in SEEDS:
                    d_trial = random.Random(seed + 7777).randint(1, n - 1)
                    sigs = gen_signatures(G, d_trial, mm, n, lam, p, k1, k2, seed)
                    if len(sigs) < mm: continue
                    t0 = time.time()
                    ok = dfn(sigs, n, lam, k1, k2, d_trial)
                    if time.time() - t0 > 120:
                        cells.append(-1); break
                    wins += ok
                else:
                    cells.append(wins)
                    continue
                break
            padded = cells + ['skip'] * (len(M_SWEEP) - len(cells))
            print(f"  {dname:<8}" + "".join(
                f"{(str(v) + '/5' if isinstance(v, int) and v >= 0 else str(v)):>8}"
                for v in padded))
        print()

print("\n" + "=" * 78)
print("done")
print("=" * 78)
