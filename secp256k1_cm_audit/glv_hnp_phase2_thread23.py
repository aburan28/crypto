"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector
can actually be lambda_1.

Context (RESEARCH_AUTOLAB_LOG.md 2026-07-29, EXP T5):
  The Phase-2 lattice of glv_hnp_phase2_20bit.py:262 always has
      lambda_1 = n * S_D * e_m      (the "trivial" d-direction vector)
  because d is only defined mod n.  |sv[m]|/n = 1.0000 on every curve tested,
  and sv/pv sits in [0.337, 0.368] for successes and failures alike.
  That entry asserted "no choice of S_D removes it -- both vectors scale
  linearly in S_D".  EXP U1 below tests that assertion directly.

Experiments:
  U0  baseline reproduction of the 2026-07-29 T4 K1-grid (exact seeds)
  U1  S_D sweep on the ORIGINAL lattice          (tests the S_D claim)
  U2  d-ELIMINATED lattice (dim 2m+1, no d coordinate) on the same K1 grid
  U3  anatomy of lambda_1 in the d-eliminated lattice
  U4  quantitative model: the 2D lambda-block shortest vector mu vs ||pv||

Run: python3 glv_hnp_phase2_thread23.py
Deps: fpylll, sympy   (pip install cysignals fpylll sympy)
"""

import math
import random

from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + GLV helpers
# (copied verbatim from glv_hnp_phase2_lambda_threshold.py so that the
#  comparison against the 2026-07-29 numbers is exact)
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

def tonelli_shanks(nn, p):
    nn %= p
    if nn == 0: return 0
    if pow(nn, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(nn, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(nn, q, p), pow(nn, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0: break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u): break
        u, v = v, u
    return u

def lambda_block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1])

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- the 2026-06-15 column scaling."""
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
# FORMULATION A -- the original Phase-2 lattice, with S_D exposed
# dim 2m+2 : [k1-congruence cols 0..m-1 | d col m | k2 cols m+1..2m | kannan]
# ---------------------------------------------------------------------------

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound, S_D=None):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D_def, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if S_D is None:
        S_D = S_D_def
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

def planted_A(sigs, d_secret, n, k1_bound, k2_bound, S_D=None):
    m = len(sigs)
    S_K1, S_D_def, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if S_D is None:
        S_D = S_D_def
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_A(rows, m, n, S_KAN, S_D, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if S_D != 1 and val % S_D != 0: continue
        d_cand = (val // S_D) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# FORMULATION B -- d ELIMINATED via signature 0
#
#   d == B0^{-1} (k1_0 + lam*k2_0 - A0)               (mod n)
#   => k1_j == C_j + D_j*k1_0 + lam*D_j*k2_0 - lam*k2_j   (mod n),  j=1..m-1
#      with D_j = B_j/B_0,  C_j = A_j - D_j*A_0.
#
# generators : k1_0, k2_0, ..., k2_{m-1}          (m+1 of them)
# residuals  : k1_1..k1_{m-1} in m-1 congruence columns
# dim = (m-1) + 1 + m + 1 = 2m+1   -- exactly one less than formulation A,
# and the n*e_d direction is gone.
#
# columns: 0..m-2   congruence cols for j=1..m-1     (scale S_K1)
#          m-1      k1_0                             (scale S_K1)
#          m..2m-1  k2_0..k2_{m-1}                   (scale S_K2)
#          2m       kannan                           (scale S_KAN)
# ---------------------------------------------------------------------------

def build_lattice_B(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    B0i = modinv(B0, n)
    D = [0] * m
    C = [0] * m
    for j in range(1, m):
        D[j] = sigs[j]['B'] * B0i % n
        C[j] = (sigs[j]['A'] - D[j] * A0) % n

    M = [[0] * dim for _ in range(dim)]
    r = 0
    # congruence reduction rows
    for j in range(1, m):
        M[r][j - 1] = n * S_K1
        r += 1
    # generator row: k1_0
    for j in range(1, m):
        M[r][j - 1] = D[j] * S_K1
    M[r][m - 1] = S_K1
    r += 1
    # generator row: k2_0
    for j in range(1, m):
        M[r][j - 1] = lam * D[j] % n * S_K1
    M[r][m] = S_K2
    r += 1
    # generator rows: k2_l, l = 1..m-1
    for l in range(1, m):
        M[r][l - 1] = -lam * S_K1
        M[r][m + l] = S_K2
        r += 1
    # kannan row (constants)
    for j in range(1, m):
        M[r][j - 1] = C[j] * S_K1
    M[r][2 * m] = S_KAN
    r += 1
    assert r == dim, (r, dim)
    return M, (D, C, B0i, A0)

def planted_B(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * dim
    for j in range(1, m):
        v[j - 1] = sigs[j]['k1'] * S_K1
    v[m - 1] = sigs[0]['k1'] * S_K1
    for l in range(m):
        v[m + l] = sigs[l]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_B(rows, m, n, S_K1, S_K2, S_KAN, aux, d_secret):
    """Read k1_0, k2_0 off any row with |kannan| == S_KAN, then invert."""
    dim = 2 * m + 1
    D, C, B0i, A0 = aux
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        a = sign * row[m - 1]
        b = sign * row[m]
        if a % S_K1 or b % S_K2: continue
        k1_0 = a // S_K1
        k2_0 = b // S_K2
        d_cand = B0i * (k1_0 + LAM_CUR * k2_0 - A0) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

LAM_CUR = None   # set per experiment (recover_B needs lam)

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_A(curve, m, d_secret, k1_bound, seed, S_D=None, use_bkz=False, beta=20):
    global LAM_CUR
    p, b, n, lam, G = curve
    LAM_CUR = lam
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D_def, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if S_D is None:
        S_D = S_D_def
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_lattice_A(sigs, n, lam, k1_bound, k2_bound, S_D)
    dim = 2 * m + 2
    X = IntegerMatrix.from_matrix(M)
    if use_bkz: BKZ.reduction(X, BKZ.Param(beta))
    else:       LLL.reduction(X)
    rows = [[X[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_A(rows, m, n, S_KAN, S_D, d_secret) is not None
    pv = planted_A(sigs, d_secret, n, k1_bound, k2_bound, S_D)
    sv = min(rows, key=norm)
    return ok, norm(pv), norm(sv), sv, rows, sigs

def run_B(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    global LAM_CUR
    p, b, n, lam, G = curve
    LAM_CUR = lam
    k2_bound = math.isqrt(n) + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M, aux = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    X = IntegerMatrix.from_matrix(M)
    if use_bkz: BKZ.reduction(X, BKZ.Param(beta))
    else:       LLL.reduction(X)
    rows = [[X[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_B(rows, m, n, S_K1, S_K2, S_KAN, aux, d_secret) is not None
    pv = planted_B(sigs, n, k1_bound, k2_bound)
    sv = min(rows, key=norm)
    return ok, norm(pv), norm(sv), sv, rows, sigs

def rate(fn, curve, m, k1_bound, seeds, **kw):
    wins, ratios = 0, []
    for seed in seeds:
        n = curve[2]
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = fn(curve, m, d_trial, k1_bound, seed, **kw)
        if res is None: continue
        ok, pn, sn = res[0], res[1], res[2]
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================
print("=" * 78)
print("Thread 23 -- can the planted vector be made lambda_1?  (GLV-HNP Phase 2)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 12                      # same as EXP T4 (2026-07-29)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

HIST = [
    # label,              p,    b, n,    lam
    ("12-bit/2557",       2557, 2, 2659, 1755),
    ("12-bit/2677 FAIL",  2677, 2, 2647, 185),
]
curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G)))

for label, c in curves:
    p, b, n, lam, G = c
    k2 = math.isqrt(n) + 1
    print(f"  {label:<18} p={p} n={n} lam={lam} lam*={lam_star(lam,n):.4f} K2={k2}")
print(f"  m = {M_SIGS} signatures, seeds = {SEEDS}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U0: baseline -- reproduce the 2026-07-29 T4 grid (formulation A)")
print("-" * 78)
print(f"{'curve':<18} {'lam*':>7} " + " ".join(f"K1={k:<4}" for k in K1_GRID))
base = {}
for label, c in curves:
    cells = []
    for k1 in K1_GRID:
        w, t, _ = rate(run_A, c, M_SIGS, k1, SEEDS)
        base[(label, k1)] = w
        cells.append(f"{w}/{t} ")
    print(f"{label:<18} {lam_star(c[3],c[2]):>7.4f} " + " ".join(cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U1: S_D sweep -- is lambda_1 = n*S_D*e_m really unavoidable?")
print("-" * 78)
print("2026-07-29 claim: 'no choice of S_D removes it -- both vectors scale")
print("linearly in S_D'.  But ||n*S_D*e_m|| = n*S_D scales fully, while only")
print("ONE component of the planted vector does:")
print("   ||pv||^2 = (d*S_D)^2 + C^2,   C^2 = sum of the k1/k2/kannan blocks.")
print("Since d < n strictly, ||pv|| < n*S_D as soon as S_D > C/sqrt(n^2-d^2).")
print("Prediction: a large enough S_D makes the trivial vector NOT lambda_1.\n")

SD_GRID = [1, 2, 4, 8, 16, 32, 64, 256]
for label, c in curves:
    p, b, n, lam, G = c
    for k1 in (4, 8):
        print(f"\n  {label}, K1={k1}  (eff={k1*(math.isqrt(n)+1)/n:.3f})")
        print(f"    {'S_D':>5} {'succ':>6} {'sv/pv':>7} {'sv=triv?':>9} "
              f"{'|sv[m]|/(n*S_D)':>16} {'pv/triv':>8}")
        for sd in SD_GRID:
            w, t, mr = rate(run_A, c, M_SIGS, k1, SEEDS, S_D=sd)
            # anatomy on seed 42
            d0 = random.Random(42 + 7777).randint(1, n - 1)
            res = run_A(c, M_SIGS, d0, k1, 42, S_D=sd)
            ok, pn, sn, sv, rows, sigs = res
            m = M_SIGS
            triv_norm = n * sd
            is_triv = (abs(sv[m]) == triv_norm and
                       all(sv[i] == 0 for i in range(2 * m + 2) if i != m))
            print(f"    {sd:>5} {w}/{t:<4} {mr:>7.3f} {str(is_triv):>9} "
                  f"{abs(sv[m])/triv_norm:>16.4f} {pn/triv_norm:>8.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U2: formulation B (d eliminated, dim 2m+1) on the same K1 grid")
print("-" * 78)
print(f"{'curve':<18} {'form':>5} " + " ".join(f"K1={k:<4}" for k in K1_GRID))
for label, c in curves:
    cells_b = []
    for k1 in K1_GRID:
        w, t, _ = rate(run_B, c, M_SIGS, k1, SEEDS)
        cells_b.append(f"{w}/{t} ")
    print(f"{label:<18} {'A':>5} " +
          " ".join(f"{base[(label,k1)]}/{len(SEEDS)} " for k1 in K1_GRID))
    print(f"{'':<18} {'B':>5} " + " ".join(cells_b))

print("\nB + BKZ(beta=20) on the two grid cells nearest each wall:")
for label, c in curves:
    for k1 in (6, 8, 12, 16):
        w, t, _ = rate(run_B, c, M_SIGS, k1, SEEDS, use_bkz=True, beta=20)
        wa, _, _ = rate(run_A, c, M_SIGS, k1, SEEDS, use_bkz=True, beta=20)
        print(f"  {label:<18} K1={k1:<3} A-BKZ {wa}/{t}   B-BKZ {w}/{t}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U3: what IS lambda_1 in formulation B?")
print("-" * 78)
print("cols: 0..m-2 congruence(k1_1..k1_{m-1}), m-1 k1_0, m..2m-1 k2, 2m kannan\n")
print(f"{'curve':<18} {'K1':>3} {'sv/pv':>7} {'k1-blk':>7} {'k2-blk':>7} "
      f"{'kan':>7} {'kan=S?':>7} {'sv=pv?':>7}")
for label, c in curves:
    p, b, n, lam, G = c
    for k1 in (4, 8, 16):
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        res = run_B(c, M_SIGS, d0, k1, 42)
        ok, pn, sn, sv, rows, sigs = res
        m = M_SIGS
        dim = 2 * m + 1
        S_K1, _, S_K2, S_KAN = scales(n, k1, math.isqrt(n) + 1)
        tot = sum(x * x for x in sv) or 1
        e_k1 = sum(sv[i] ** 2 for i in range(m)) / tot
        e_k2 = sum(sv[m + i] ** 2 for i in range(m)) / tot
        e_kan = sv[dim - 1] ** 2 / tot
        pv = planted_B(sigs, n, k1, math.isqrt(n) + 1)
        same = all(abs(a) == abs(bb) for a, bb in zip(sv, pv))
        print(f"{label:<18} {k1:>3} {sn/pn:>7.3f} {e_k1:>7.3f} {e_k2:>7.3f} "
              f"{e_kan:>7.3f} {str(abs(sv[dim-1])==S_KAN):>7} {str(same):>7}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U4: model -- the 2D lambda-block vector mu is the real competitor")
print("-" * 78)
print("The block < (n*S_K1, 0), (-lam*S_K1, S_K2) > survives BOTH formulations")
print("(it has kannan coord 0, so it can never encode d).  Hermite bound:")
print("   mu <= sqrt(2/sqrt(3) * n*S_K1*S_K2) ~ n^1.5 * sqrt(1.155/(K1*K2))")
print("   ||pv_B|| ~ n * sqrt(2m/3 + 1)")
print("   => mu < ||pv|| iff eff = K1*K2/n > 1.155/(2m/3+1)\n")
eff_star = 1.155 / (2 * M_SIGS / 3.0 + 1.0)
print(f"  predicted wall at eff* = {eff_star:.4f}  (m={M_SIGS})\n")
print(f"{'curve':<18} {'K1':>3} {'eff':>6} {'mu':>12} {'||pv_B||':>12} "
      f"{'mu/pv':>7} {'A':>5} {'B':>5}")
for label, c in curves:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        S_K1, _, S_K2, S_KAN = scales(n, k1, k2b)
        mu = lambda_block_mu(n, lam, S_K1, S_K2)
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        res = run_B(c, M_SIGS, d0, k1, 42)
        pvb = res[1]
        wa = base[(label, k1)]
        wb, t, _ = rate(run_B, c, M_SIGS, k1, SEEDS)
        print(f"{label:<18} {k1:>3} {k1*k2b/n:>6.3f} {mu:>12.1f} {pvb:>12.1f} "
              f"{mu/pvb:>7.3f} {wa}/{len(SEEDS)} {wb}/{t}")


# ---------------------------------------------------------------------------
# EXP U5 -- the uSVP / Gaussian-heuristic predictor
#
# Both formulations are block-triangular, so their determinants are exact:
#   det_A = (n*S_K1)^m     * S_D * S_K2^m * S_KAN        (dim 2m+2)
#   det_B = (n*S_K1)^(m-1) * S_K1 * S_K2^m * S_KAN       (dim 2m+1)
# GH(L) = det^(1/dim) * sqrt(dim / (2*pi*e)) is the expected lambda_1 of a
# "random" lattice of that volume.  The planted vector is findable only if it
# is anomalously short, i.e. tau = ||pv|| / GH(L) < 1.
# ---------------------------------------------------------------------------

def det_A(n, m, S_K1, S_D, S_K2, S_KAN):
    return (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN

def det_B(n, m, S_K1, S_K2, S_KAN):
    return (n * S_K1) ** (m - 1) * S_K1 * S_K2 ** m * S_KAN

def gh(det, dim):
    """Gaussian heuristic for lambda_1 -- computed in logs (det is huge)."""
    log_det = math.log(det)
    return math.exp(log_det / dim) * math.sqrt(dim / (2 * math.pi * math.e))

print("\n" + "-" * 78)
print("EXP U5: uSVP predictor  tau = ||pv|| / GH(L)")
print("-" * 78)
print("If tau < 1 the planted vector is anomalously short and LLL/BKZ can")
print("isolate it; if tau > 1 it is a typical lattice vector and no amount of")
print("reduction distinguishes it.  Note det_B = det_A / (n*S_D) and")
print("dim_B = dim_A - 1, which predicts tau_A ~ tau_B -- i.e. EXP U2's")
print("cell-for-cell tie is not a coincidence.\n")
print(f"{'curve':<18} {'K1':>3} {'eff':>6} {'GH_A':>10} {'tau_A':>7} "
      f"{'GH_B':>10} {'tau_B':>7} {'A':>5} {'B':>5}")
u5_rows = []
for label, c in curves:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    m = M_SIGS
    for k1 in K1_GRID:
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
        gA = gh(det_A(n, m, S_K1, S_D, S_K2, S_KAN), 2 * m + 2)
        gB = gh(det_B(n, m, S_K1, S_K2, S_KAN), 2 * m + 1)
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        rA = run_A(c, m, d0, k1, 42)
        rB = run_B(c, m, d0, k1, 42)
        wa = base[(label, k1)]
        wb, t, _ = rate(run_B, c, m, k1, SEEDS)
        print(f"{label:<18} {k1:>3} {k1*k2b/n:>6.3f} {gA:>10.1f} "
              f"{rA[1]/gA:>7.3f} {gB:>10.1f} {rB[1]/gB:>7.3f} "
              f"{wa}/{len(SEEDS)} {wb}/{t}")
        u5_rows.append({'label': label, 'K1': k1, 'lam_star': lam_star(lam, n),
                        'eff': k1 * k2b / n, 'tau': rB[1] / gB,
                        'wins': wb, 'tot': t})

print("\nSeparation of the 5/5 cells from the 0/5 cells, both curves pooled:")
succ = [r['tau'] for r in u5_rows if r['wins'] == r['tot']]
fail = [r['tau'] for r in u5_rows if r['wins'] == 0]
if succ and fail:
    print(f"  tau on 5/5 cells : [{min(succ):.3f}, {max(succ):.3f}]  (n={len(succ)})")
    print(f"  tau on 0/5 cells : [{min(fail):.3f}, {max(fail):.3f}]  (n={len(fail)})")
    print(f"  separated = {max(succ) < min(fail)}")
    for key in ('eff', 'lam_star'):
        s = [r[key] for r in u5_rows if r['wins'] == r['tot']]
        f = [r[key] for r in u5_rows if r['wins'] == 0]
        print(f"  {key:<9} 5/5 [{min(s):.3f}, {max(s):.3f}]   "
              f"0/5 [{min(f):.3f}, {max(f):.3f}]   "
              f"separated = {max(s) < min(f)}")

# ---------------------------------------------------------------------------
# EXP U6 -- validate tau out of sample on freshly-searched curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                bb = num // 2
                if bb >= 0 and a * a - a * bb + bb * bb == p:
                    return (a, bb)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

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

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, want, nbins=8):
    """j=0 GLV curves with n prime, n = 1 mod 3, spread over lam* in [0,0.5]."""
    import sympy
    bins = {i: [] for i in range(nbins)}
    per = max(1, want // nbins)
    p = int(sympy.nextprime(lo))
    out = []
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1: continue
                    if not sympy.isprime(nc): continue
                    roots = glv_roots(nc)
                    if roots is None: continue
                    lam = roots[0]
                    ls = lam_star(lam, nc)
                    idx = min(nbins - 1, int(ls / (0.5 / nbins)))
                    if len(bins[idx]) >= per: continue
                    cur = build_curve(p, nc)
                    if cur is None: continue
                    bins[idx].append((p, cur[1], nc, lam, cur[4]))
        if sum(len(v) for v in bins.values()) >= want: break
        p = int(sympy.nextprime(p))
    for i in range(nbins):
        out.extend(bins[i])
    return out

print("\n" + "-" * 78)
print("EXP U6: out-of-sample test of tau on fresh 17-bit curves")
print("-" * 78)
FRESH = search_curves(2 ** 16, 2 ** 17, want=16, nbins=8)
print(f"found {len(FRESH)} fresh j=0 GLV curves, "
      f"lam* in [{min(lam_star(c[3],c[2]) for c in FRESH):.4f}, "
      f"{max(lam_star(c[3],c[2]) for c in FRESH):.4f}]")
print(f"\n{'p':>8} {'n':>8} {'lam*':>7} {'K1':>4} {'eff':>6} {'tau':>7} {'B':>5}")
u6 = []
for c in FRESH:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    m = M_SIGS
    for k1 in (4, 12, 40):
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
        gB = gh(det_B(n, m, S_K1, S_K2, S_KAN), 2 * m + 1)
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        rB = run_B(c, m, d0, k1, 42)
        if rB is None: continue
        tau = rB[1] / gB
        w, t, _ = rate(run_B, c, m, k1, SEEDS)
        print(f"{p:>8} {n:>8} {lam_star(lam,n):>7.4f} {k1:>4} "
              f"{k1*k2b/n:>6.3f} {tau:>7.3f} {w}/{t}")
        u6.append({'tau': tau, 'eff': k1 * k2b / n,
                   'lam_star': lam_star(lam, n), 'wins': w, 'tot': t})

if u6:
    print("\nPooled 17-bit cells (majority-vote outcome = recovered if wins>=3):")
    for key in ('tau', 'eff', 'lam_star'):
        s = [r[key] for r in u6 if r['wins'] >= 3]
        f = [r[key] for r in u6 if r['wins'] < 3]
        if not s or not f:
            print(f"  {key}: degenerate")
            continue
        # AUC by rank
        pairs = [(1 if a < bb else 0.5 if a == bb else 0) for a in s for bb in f]
        auc = sum(pairs) / len(pairs)      # small key -> success
        auc = max(auc, 1 - auc)
        print(f"  {key:<9} succ [{min(s):.3f}, {max(s):.3f}] (n={len(s)})  "
              f"fail [{min(f):.3f}, {max(f):.3f}] (n={len(f)})  "
              f"overlap={not (max(s) < min(f) or max(f) < min(s))}  AUC={auc:.3f}")

print("\nExact tau at the knife-edge cell K1=8 (both curves), full precision:")
for label, c in curves:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, 8, k2b)
    gB = gh(det_B(n, M_SIGS, S_K1, S_K2, S_KAN), 2 * M_SIGS + 1)
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    rB = run_B(c, M_SIGS, d0, 8, 42)
    print(f"  {label:<18} tau_B = {rB[1]/gB:.6f}   outcome "
          f"{base[(label,8)]}/{len(SEEDS)}")
print("  -> the two curves sit at the SAME tau with OPPOSITE outcomes.")
print("     tau is a reparametrisation of eff (det_B = n^m*S_K1^m*S_K2^m")
print("     depends on K1*K2 only, and ||pv|| ~ n*sqrt(2m/3+1) is K-free), so")
print("     it carries NO curve-specific information.  It predicts WHERE the")
print("     wall is from first principles; it does not explain lam*.")

# ---------------------------------------------------------------------------
# EXP U7 -- does MORE DATA cross the wall?  tau as a function of m.
#
# Asymptotically  det_B^(1/dim) -> sqrt(n^3/(K1*K2)) and GH ~ sqrt(m/(pi*e)),
# while ||pv|| ~ n*sqrt(2m/3):  both grow as sqrt(m), so
#        tau(m) -> sqrt(2*pi*e/3 * eff)     -- INDEPENDENT of m.
# Prediction: no number of signatures crosses the wall.  This is the
# quantitative explanation of T4b (2026-07-29: m=8..32 at K1=8 all failed).
# ---------------------------------------------------------------------------

print("\n" + "-" * 78)
print("EXP U7: tau(m) -- can more signatures cross the wall?")
print("-" * 78)
print("asymptote  tau_inf = sqrt(2*pi*e/3 * eff)\n")
M_GRID = [6, 8, 12, 16, 24, 32, 48]
for label, c in curves:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    for k1 in (4, 8):
        eff = k1 * k2b / n
        tinf = math.sqrt(2 * math.pi * math.e / 3.0 * eff)
        print(f"\n  {label}, K1={k1}, eff={eff:.3f}, tau_inf={tinf:.3f}")
        print(f"    {'m':>4} {'dim':>5} {'tau(m)':>8} {'succ':>6}")
        for m in M_GRID:
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
            gB = gh(det_B(n, m, S_K1, S_K2, S_KAN), 2 * m + 1)
            d0 = random.Random(42 + 7777).randint(1, n - 1)
            rB = run_B(c, m, d0, k1, 42)
            if rB is None:
                print(f"    {m:>4} {'-':>5} {'-':>8} {'(no sigs)':>6}")
                continue
            w, t, _ = rate(run_B, c, m, k1, SEEDS)
            print(f"    {m:>4} {2*m+1:>5} {rB[1]/gB:>8.3f} {w}/{t:<4}")


# ---------------------------------------------------------------------------
# EXP U8 -- is the lam*=0.07 failure information-theoretic or algorithmic?
#
# n is tiny (2647), so enumerate ALL candidate secrets d' and count how many
# are consistent with every signature, i.e. how many satisfy
#     (A_i + B_i*d') mod n  in  S = {k1 + lam*k2 mod n : k1<K1, k2<K2}
# for all i.  If the count is 1, the instance is uniquely determined and the
# failure is purely a lattice-reduction failure.  If it is >1, no algorithm
# can succeed and the "wall" is information-theoretic.
# ---------------------------------------------------------------------------

print("\n" + "-" * 78)
print("EXP U8: brute-force uniqueness of d (n is small enough to enumerate)")
print("-" * 78)
print("|S| = K1*K2 distinct nonces reachable; eff = |S|/n is the per-signature")
print("survival probability of a wrong d', so E[#wrong survivors] ~ n*eff^m.\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} {'|S|':>6} {'eff':>6} "
      f"{'#consistent':>12} {'E[wrong]':>10} {'LLL':>5}")
for label, c in curves:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    for k1 in (4, 8):
        S = set()
        for a in range(k1):
            for bb in range(k2b):
                S.add((a + lam * bb) % n)
        for m in (6, 12, 24):
            d0 = random.Random(42 + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
            if len(sigs) < m: continue
            cons = []
            for dp in range(1, n):
                good = True
                for s in sigs:
                    if (s['A'] + s['B'] * dp) % n not in S:
                        good = False
                        break
                if good:
                    cons.append(dp)
            assert d0 in cons, "true secret must be consistent"
            eff = len(S) / n
            w, t, _ = rate(run_B, c, m, k1, SEEDS)
            print(f"{label:<18} {k1:>3} {m:>3} {len(S):>6} {eff:>6.3f} "
                  f"{len(cons):>15} {n*eff**m:>10.2e} {w}/{t}")

print("\nIf #consistent == 1 while LLL reports 0/5, the instance IS uniquely")
print("solvable and the Phase-2 lattice is simply failing to find it.")

# ---------------------------------------------------------------------------
# EXP U9 -- if the instance is unique, is it a REDUCTION-QUALITY failure?
# Push beta up to the full dimension (beta = dim is HKZ, i.e. an exact SVP
# oracle on every projected block).  If even HKZ misses d, the planted vector
# is not the minimum of its coset and the ENCODING, not the reduction, is at
# fault.
# ---------------------------------------------------------------------------

print("\n" + "-" * 78)
print("EXP U9: BKZ beta sweep up to HKZ (beta = dim) on the failing cell")
print("-" * 78)
print(f"{'curve':<18} {'m':>3} {'dim_B':>6} " +
      " ".join(f"b={b:<4}" for b in (2, 10, 20, 25)) + "  HKZ")
for label, c in curves:
    for m in (8, 12):
        dim = 2 * m + 1
        cells = []
        for beta in (2, 10, 20, 25):
            bb = min(beta, dim)
            w, t, _ = rate(run_B, c, m, 8, SEEDS, use_bkz=True, beta=bb)
            cells.append(f"{w}/{t}  ")
        w, t, _ = rate(run_B, c, m, 8, SEEDS, use_bkz=True, beta=dim)
        print(f"{label:<18} {m:>3} {dim:>6} " + " ".join(cells) + f"  {w}/{t}")
print("\n(K1=8 throughout -- the cell where the two curves disagree.)")

print("\nDone.")
