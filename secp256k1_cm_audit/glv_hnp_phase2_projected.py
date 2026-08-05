"""
GLV-HNP Phase 2, Thread 23: make the planted vector reachable as lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, exp T5):
  The Phase-2 lattice L (dim 2m+2, built by build_glv_lattice) always
  contains the trivial vector  n*S_D*e_m  of norm n*S_D, because the d
  column carries a full n*Z sublattice (d is only defined mod n).  With
  ||v_planted|| ~ n*sqrt(2m/3 + 4/3) > n, the planted vector is NEVER
  lambda_1 of L: measured sv/pv in [0.34, 0.61] on every curve, success
  and failure alike.  Recovery is therefore a BDD/coset condition, not an
  SVP condition, which is why six consecutive curve-level separators
  (2026-06-21 .. 2026-06-29) all failed.

  Proposed fix (Thread 23): quotient out the trivial direction e_m and
  solve the resulting problem directly.

Two reformulations are implemented and compared against the original:

  FULL   the 2026-06-15 lattice, verbatim (control).

  PROJ   the projected lattice  L' = L / (L cap R*e_m),  realised by
         deleting column m from the generator matrix.  L' has rank 2m+1
         and 2m+2 generators (n*drow lies in the span of the modular
         rows), so fpylll reduces a rank-deficient input and emits one
         zero row.  d is no longer a lattice coordinate; it is recovered
         algebraically from signature 0 once (k1_0, k2_0) are known:
             d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  mod n.
         det(L') = det(L)/(n*S_D) = (n*S_K1)^m * S_K2^m * S_KANNAN / n.

  BABAI  drop the Kannan embedding as well.  L0 = <n*S_K1*e_i, drow,
         k2row_i> in dimension 2m (k1 cols + k2 cols), target
             t = (-A_i * S_K1)_{i<m} , 0 ... 0
         whose closest lattice vector differs from t by exactly
             (k1_i*S_K1 , k2_i*S_K2),   norm ~ n*sqrt(2m/3).
         Solved by Babai nearest-plane on an LLL-reduced basis.  This
         removes the S_KANNAN = n component that the embedding adds to
         the target norm.

Hypotheses under test:

  H23a  In L', the planted vector becomes lambda_1 (sv/pv -> 1.0 instead
        of 0.34-0.61).
  H23b  If H23a holds, the K1 wall measured in T4 (2026-07-29) moves
        outward: 12-bit/2677 currently fails for K1 >= 6 at m=10.
        FALSIFIER: if the wall stays at K1 ~ 4-6 the wall is
        information-theoretic and Phase 2 is at its ceiling.
  H23c  ||v_planted|| / GH(L') with GH = sqrt(N/2*pi*e) * det^(1/N)
        predicts success with a threshold near 1.

Experiments:
  P1  norms and lambda_1 status of FULL vs PROJ on the 3 historical curves
  P2  K1 wall: T4 grid re-run for FULL / PROJ / BABAI, 5 seeds
  P3  is ||pv||/||sv|| (or ||pv||/GH) a predictor across curves and K1?
  P4  does more data rescue the wall in PROJ?  (T4b said no for FULL)

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

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
    mm, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures and scales (verbatim)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- 2026-06-15 column-diagonal scaling."""
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
# FULL lattice (control) -- verbatim from glv_hnp_phase2_20bit.py:262
# ---------------------------------------------------------------------------

def build_full(sigs, n, lam, k1_bound, k2_bound):
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

def planted_full(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_full(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        if abs(row[dim - 1]) != S_KAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# PROJ lattice:  L' = L / (L cap R*e_m)   (delete column m)
# ---------------------------------------------------------------------------
# columns of L':  0..m-1  = k1 cols,  m..2m-1 = k2 cols,  2m = Kannan col.

def build_proj(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # modular rows
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim                            # d row (d coordinate deleted)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                       # k2 rows
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * dim                            # Kannan / A row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KAN
    rows.append(r)
    return rows                              # 2m+2 generators, rank 2m+1

def planted_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def d_from_k(k1_0, k2_0, sigs, n, lam):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    return (k1_0 + lam * k2_0 - A0) * modinv(B0, n) % n

def recover_proj(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        if abs(row[2 * m]) != S_KAN: continue
        sign = 1 if row[2 * m] > 0 else -1
        v = [sign * x for x in row]
        if v[0] % S_K1 or v[m] % S_K2: continue
        k1_0, k2_0 = v[0] // S_K1, v[m] // S_K2
        d_cand = d_from_k(k1_0, k2_0, sigs, n, lam)
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# BABAI:  CVP form, no Kannan embedding, no d coordinate
# ---------------------------------------------------------------------------
# L0 in dim 2m: cols 0..m-1 = k1 cols, m..2m-1 = k2 cols.
# target t = (-A_i*S_K1, 0);  closest vector w satisfies w - t = planted.

def build_babai(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    t = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    return rows, t                           # 2m+1 generators, rank 2m

def gram_schmidt(B):
    """Float GS.  Entries here are <= 2^40 for the curves used, so f64 is
    exact enough; a scaled variant would be needed for cryptographic sizes
    (cf. src/cryptanalysis/lattice.rs big_to_f64_scaled)."""
    n_rows = len(B)
    bstar, mu = [], [[0.0] * n_rows for _ in range(n_rows)]
    for i in range(n_rows):
        v = [float(x) for x in B[i]]
        for j in range(len(bstar)):
            nj = sum(y * y for y in bstar[j])
            if nj == 0:
                mu[i][j] = 0.0
                continue
            mu[i][j] = sum(float(B[i][k]) * bstar[j][k] for k in range(len(v))) / nj
            for k in range(len(v)):
                v[k] -= mu[i][j] * bstar[j][k]
        bstar.append(v)
    return bstar, mu

def babai_nearest_plane(B, t):
    """Babai nearest plane; B must be LLL-reduced and full rank."""
    bstar, _ = gram_schmidt(B)
    b = [float(x) for x in t]
    coeffs = [0] * len(B)
    for i in range(len(B) - 1, -1, -1):
        nj = sum(y * y for y in bstar[i])
        if nj == 0: continue
        c = sum(b[k] * bstar[i][k] for k in range(len(b))) / nj
        ci = int(math.floor(c + 0.5))
        coeffs[i] = ci
        for k in range(len(b)):
            b[k] -= ci * float(B[i][k])
    w = [0] * len(t)
    for i, ci in enumerate(coeffs):
        if ci:
            for k in range(len(t)):
                w[k] += ci * B[i][k]
    return w

def reduce_rows(rows, dim):
    """LLL-reduce a (possibly rank-deficient) generator set; returns the
    reduced rows with all-zero rows dropped."""
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    out = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    return [r for r in out if any(r)]

# ---------------------------------------------------------------------------
# One trial, all three methods
# ---------------------------------------------------------------------------

def gaussian_heuristic(log_det, N):
    """GH = sqrt(N/(2*pi*e)) * det^(1/N), computed in logs."""
    return math.sqrt(N / (2 * math.pi * math.e)) * math.exp(log_det / N)

def trial(curve, m, k1_bound, seed, d_secret=None):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    if d_secret is None:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    res = {'n': n, 'lam': lam, 'm': m, 'K1': k1_bound, 'K2': k2_bound,
           'eff': k1_bound * k2_bound / n}

    # --- FULL ---
    Mf = build_full(sigs, n, lam, k1_bound, k2_bound)
    rf = reduce_rows(Mf, 2 * m + 2)
    pvf = norm(planted_full(sigs, d_secret, n, k1_bound, k2_bound))
    svf = min(norm(r) for r in rf)
    res['full_ok'] = recover_full(rf, m, n, S_KAN, d_secret) is not None
    res['full_pv'] = pvf
    res['full_sv'] = svf
    res['full_sv_over_pv'] = svf / pvf
    # energy of the shortest vector in the d column
    svrow = min(rf, key=norm)
    res['full_sv_dcol'] = abs(svrow[m]) / n

    # --- PROJ ---
    Mp = build_proj(sigs, n, lam, k1_bound, k2_bound)
    rp = reduce_rows(Mp, 2 * m + 1)
    pvp = norm(planted_proj(sigs, n, k1_bound, k2_bound))
    svp = min(norm(r) for r in rp)
    res['proj_ok'] = recover_proj(rp, sigs, n, lam, k1_bound, k2_bound,
                                  d_secret) is not None
    res['proj_pv'] = pvp
    res['proj_sv'] = svp
    res['proj_sv_over_pv'] = svp / pvp
    N = 2 * m + 1
    log_det = m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_KAN) \
        - math.log(n)
    res['proj_gh'] = gaussian_heuristic(log_det, N)
    res['proj_pv_over_gh'] = pvp / res['proj_gh']

    # --- BABAI ---
    rows_b, t = build_babai(sigs, n, lam, k1_bound, k2_bound)
    rb = reduce_rows(rows_b, 2 * m)
    w = babai_nearest_plane(rb, t)
    diff = [w[k] - t[k] for k in range(2 * m)]
    ok_b = False
    if diff[0] % S_K1 == 0 and diff[m] % S_K2 == 0:
        k1_0, k2_0 = diff[0] // S_K1, diff[m] // S_K2
        ok_b = d_from_k(k1_0, k2_0, sigs, n, lam) == d_secret
    res['babai_ok'] = ok_b
    res['babai_dist'] = norm(diff)
    res['babai_true_dist'] = norm(
        [sigs[i]['k1'] * S_K1 for i in range(m)] +
        [sigs[i]['k2'] * S_K2 for i in range(m)])
    return res


def rates(curve, m, k1_bound, seeds):
    agg = {'full': 0, 'proj': 0, 'babai': 0, 'tot': 0}
    acc = {'full_sv_over_pv': [], 'proj_sv_over_pv': [], 'proj_pv_over_gh': [],
           'full_sv_dcol': []}
    for s in seeds:
        r = trial(curve, m, k1_bound, s)
        if r is None: continue
        agg['tot'] += 1
        agg['full'] += bool(r['full_ok'])
        agg['proj'] += bool(r['proj_ok'])
        agg['babai'] += bool(r['babai_ok'])
        for k in acc:
            acc[k].append(r[k])
    mean = {k: (sum(v) / len(v) if v else float('nan')) for k, v in acc.items()}
    return agg, mean


# ===========================================================================
print("=" * 78)
print("Thread 23 — reformulate the Phase-2 lattice so the target is lambda_1")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (identical to Thread 20's HIST table)
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
print("EXP P1: is the planted vector lambda_1 after projection?")
print("-" * 78)
print("sv/pv = ||shortest reduced vector|| / ||planted vector||.")
print("T5 (2026-07-29) measured sv/pv in [0.34,0.61] for FULL, with 100% of")
print("the shortest vector's energy in the d column (dcol = |sv[m]|/n = 1).")
print("If the projection works, PROJ sv/pv should rise to ~1.0.\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} {'lam*':>6} "
      f"{'FULL sv/pv':>11} {'dcol':>6} {'PROJ sv/pv':>11} "
      f"{'pv/GH':>7} {'FULL':>5} {'PROJ':>5} {'BAB':>5}")

p1_rows = []
for label, curve, k1, m in hist:
    agg, mean = rates(curve, m, k1, SEEDS)
    ls = lam_star(curve[3], curve[2])
    print(f"{label:<18} {k1:>3} {m:>3} {ls:>6.3f} "
          f"{mean['full_sv_over_pv']:>11.3f} {mean['full_sv_dcol']:>6.3f} "
          f"{mean['proj_sv_over_pv']:>11.3f} {mean['proj_pv_over_gh']:>7.3f} "
          f"{agg['full']}/{agg['tot']:<3} {agg['proj']}/{agg['tot']:<3} "
          f"{agg['babai']}/{agg['tot']:<3}")
    p1_rows.append((label, agg, mean))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P2: the K1 wall.  T4 (2026-07-29) grid, re-run for all 3 methods")
print("-" * 78)
print("T4 result for FULL (5 seeds):")
print("  12-bit/2557 (lam*=0.340): K1=2,3,4,6,8 -> 5/5;  12 -> 4/5; 16 -> 1/5; 24 -> 0/5")
print("  12-bit/2677 (lam*=0.070): K1=2,3,4 -> 5/5;  6 -> 2/5;  8,12,16,24 -> 0/5")
print("FALSIFIER: if PROJ/BABAI do not move the wall outward, the wall is")
print("information-theoretic and Phase 2 is at its ceiling.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
p2 = {}
for label, curve, _k1, m in hist[1:]:          # the two 12-bit curves
    print(f"\n  {label}   (n={curve[2]}, lam*={lam_star(curve[3], curve[2]):.3f}, m={m})")
    print(f"  {'K1':>4} {'eff':>7} {'FULL':>6} {'PROJ':>6} {'BABAI':>6} "
          f"{'FULLsv/pv':>10} {'PROJsv/pv':>10} {'pv/GH':>7}")
    p2[label] = []
    for k1 in K1_GRID:
        agg, mean = rates(curve, m, k1, SEEDS)
        eff = k1 * (math.isqrt(curve[2]) + 1) / curve[2]
        print(f"  {k1:>4} {eff:>7.3f} "
              f"{agg['full']}/{agg['tot']:<4} {agg['proj']}/{agg['tot']:<4} "
              f"{agg['babai']}/{agg['tot']:<4} "
              f"{mean['full_sv_over_pv']:>10.3f} {mean['proj_sv_over_pv']:>10.3f} "
              f"{mean['proj_pv_over_gh']:>7.3f}")
        p2[label].append((k1, eff, agg, mean))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P3: is pv/GH (a lattice-level quantity) a predictor?")
print("-" * 78)
print("Every (curve, K1) cell of P2 is one data point: predictor value vs")
print("observed PROJ success rate.  A usable predictor separates the 5/5")
print("cells from the 0/5 cells with a single threshold across BOTH curves.\n")
print(f"{'curve':<18} {'K1':>4} {'pv/GH':>7} {'PROJsv/pv':>10} {'PROJ':>6}")
pts = []
for label in p2:
    for k1, eff, agg, mean in p2[label]:
        rate = agg['proj'] / agg['tot'] if agg['tot'] else float('nan')
        print(f"{label:<18} {k1:>4} {mean['proj_pv_over_gh']:>7.3f} "
              f"{mean['proj_sv_over_pv']:>10.3f} {agg['proj']}/{agg['tot']:<4}")
        pts.append((label, k1, mean['proj_pv_over_gh'],
                    mean['proj_sv_over_pv'], rate))

def best_threshold(points, idx):
    """Best single threshold on points[i][idx]; predicts success iff val < thr.
    Returns (threshold, accuracy, baseline)."""
    lab = [(pt[idx], pt[4] >= 0.5) for pt in points]
    cands = sorted(set(v for v, _ in lab))
    best, bacc = None, -1.0
    for c in cands:
        for thr in (c - 1e-9, c + 1e-9):
            acc = sum((v < thr) == ok for v, ok in lab) / len(lab)
            if acc > bacc:
                best, bacc = thr, acc
    base = max(sum(1 for _, ok in lab if ok), sum(1 for _, ok in lab if not ok)) / len(lab)
    return best, bacc, base

for name, idx in (("pv/GH", 2), ("PROJ sv/pv", 3)):
    thr, acc, base = best_threshold(pts, idx)
    print(f"\n  best single threshold on {name}: {thr:.4f} "
          f"-> accuracy {acc:.3f}  (majority baseline {base:.3f})")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P4: does more data rescue the wall under PROJ?")
print("-" * 78)
print("T4b (FULL, 12-bit/2677, K1=8): m = 8/12/16/24/32 -> 0,0,1,0,1 of 5.\n")
label, curve, _k1, _m = hist[2]
print(f"  {label}, K1=8")
print(f"  {'m':>4} {'FULL':>6} {'PROJ':>6} {'BABAI':>6} {'PROJsv/pv':>10} {'pv/GH':>7}")
for m in (8, 12, 16, 24, 32):
    agg, mean = rates(curve, m, 8, SEEDS)
    print(f"  {m:>4} {agg['full']}/{agg['tot']:<4} {agg['proj']}/{agg['tot']:<4} "
          f"{agg['babai']}/{agg['tot']:<4} "
          f"{mean['proj_sv_over_pv']:>10.3f} {mean['proj_pv_over_gh']:>7.3f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P5: exact lambda_1 of the projected lattice as a separator")
print("-" * 78)
print("P3 shows pv/GH separates at ~1.01 with 87.5% accuracy, the two misses")
print("being the lam*=0.340 curve at K1=6,8 (pv/GH=1.33,1.46 yet 5/5).  So the")
print("residual lam* effect should be a lambda_1(L') vs GH(L') GAP.  Test it")
print("directly: lambda_1 via BKZ at block_size = dim (HKZ in these small")
print("dimensions; fpylll SVP.shortest_vector is unavailable, the fplll")
print("strategies json is not shipped with the wheel).\n")
print(f"{'curve':<18} {'K1':>4} {'pv':>10} {'lam1':>10} {'GH':>10} "
      f"{'lam1/GH':>8} {'pv/lam1':>8} {'PROJ':>6}")

def lambda1_proj(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    rows = build_proj(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    red = reduce_rows(rows, dim)
    A = IntegerMatrix.from_matrix(red)
    BKZ.reduction(A, BKZ.Param(block_size=min(dim, 40)))
    l1 = min(norm([A[i][j] for j in range(dim)])
             for i in range(A.nrows)
             if any(A[i][j] for j in range(dim)))
    pv = norm(planted_proj(sigs, n, k1_bound, k2_bound))
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    log_det = m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_KAN) \
        - math.log(n)
    return pv, l1, gaussian_heuristic(log_det, dim)

p5 = []
for label, curve, _k1, m in hist[1:]:
    for k1, eff, agg, mean in p2[label]:
        vals = [lambda1_proj(curve, m, k1, s) for s in SEEDS]
        pv = sum(v[0] for v in vals) / len(vals)
        l1 = sum(v[1] for v in vals) / len(vals)
        gh = sum(v[2] for v in vals) / len(vals)
        rate = agg['proj'] / agg['tot'] if agg['tot'] else float('nan')
        print(f"{label:<18} {k1:>4} {pv:>10.0f} {l1:>10.0f} {gh:>10.0f} "
              f"{l1/gh:>8.3f} {pv/l1:>8.3f} {agg['proj']}/{agg['tot']:<4}")
        p5.append((label, k1, pv / gh, pv / l1, rate))

for name, idx in (("pv/GH", 2), ("pv/lambda_1", 3)):
    thr, acc, base = best_threshold(p5, idx)
    print(f"\n  best single threshold on {name}: {thr:.4f} "
          f"-> accuracy {acc:.3f}  (majority baseline {base:.3f})")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P6: m-sensitivity of the K1 wall (reconciles P2 with T4)")
print("-" * 78)
print("T4 (2026-07-29) did not record m.  P2 above uses the HIST m (8 for")
print("2557, 10 for 2677) and differs from T4 in two cells: 2557/K1=12")
print("(1/5 here vs 4/5 in T4) and 2677/K1=6 (0/5 here vs 2/5 in T4).")
print("Sweep m to find which m reproduces T4.\n")
print(f"{'curve':<18} {'K1':>4} {'m=8':>6} {'m=10':>6} {'m=12':>6} {'m=16':>6}")
for label, curve, _k1, _m in hist[1:]:
    for k1 in ((12, 16) if '2557' in label else (6, 8)):
        cells = []
        for m in (8, 10, 12, 16):
            agg, _ = rates(curve, m, k1, SEEDS)
            cells.append(f"{agg['proj']}/{agg['tot']}")
        print(f"{label:<18} {k1:>4} " + " ".join(f"{c:>6}" for c in cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P7: what IS lambda_1(L') when it is not the planted vector?")
print("-" * 78)
print("Candidate: the 2D 'lambda block' B_i = <(n*S_K1,0), (-lam*S_K1,S_K2)>")
print("that Thread 20 (T2) already isolated.  B_i lives in coordinate pair")
print("(i, m+i) of L' and is independent of i, so mu = |shortest vec of B_i|")
print("is a lattice vector of L' for every i.  Prediction: lambda_1(L') =")
print("min(mu, ||pv||), and the planted vector is lambda_1 exactly when")
print("||pv|| < mu.\n")

def gauss_reduce_2d(u, v):
    """Exact shortest vector of a 2D integer lattice (Lagrange/Gauss)."""
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

def block_mu(n, lam, S_K1, S_K2):
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1])

print(f"{'curve':<18} {'K1':>4} {'mu':>10} {'lam1(BKZ)':>10} {'pv':>10} "
      f"{'min(mu,pv)':>11} {'== lam1':>8} {'pv<mu':>6} {'PROJ':>6}")
agree = 0
total = 0
for (label, curve, _k1, m), _ in zip(hist[1:], range(2)):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    for k1, eff, agg, mean in p2[label]:
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2_bound)
        mu = block_mu(n, lam, S_K1, S_K2)
        pv, l1, gh = [sum(x) / len(SEEDS) for x in
                      zip(*[lambda1_proj(curve, m, k1, s) for s in SEEDS])]
        pred = min(mu, pv)
        ok = abs(pred - l1) / l1 < 0.02
        agree += ok; total += 1
        print(f"{label:<18} {k1:>4} {mu:>10.0f} {l1:>10.0f} {pv:>10.0f} "
              f"{pred:>11.0f} {str(ok):>8} {str(pv < mu):>6} "
              f"{agg['proj']}/{agg['tot']:<4}")
print(f"\n  lambda_1(L') == min(mu, ||pv||) in {agree}/{total} cells (2% tol)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
