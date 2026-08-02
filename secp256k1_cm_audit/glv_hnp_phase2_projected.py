"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
short, and re-measure the K1 wall.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  In the Phase-2 lattice of glv_hnp_phase2_20bit.py:262 (build_glv_lattice)
  the shortest vector is ALWAYS the trivial vector n*S_D*e_m (|sv[m]|/n =
  1.0000 exactly, Kannan coordinate 0).  It is 2-3x shorter than the planted
  vector on every curve tested, success and failure alike.  Recovery is
  therefore a BDD/coset condition, not an SVP condition.  That run proposed:

    "Thread 23 — reformulate the Phase-2 lattice so the target is lambda_1.
     Project the lattice along e_m (quotient out the trivial n*e_m
     direction) ... Falsifier: if sv/pv rises above 1 after the
     reformulation and the K1 wall in T4 moves outward on the lam*=0.07
     curve (currently K1 ~ 4-6), the reformulation is a real improvement;
     if the wall stays at K1 ~ 4-6, then the wall is information-theoretic
     and Phase 2 is at its ceiling."

Two independent defects of the 2026-06-15 lattice are fixed here:

  (P) PROJECTION.  The d-column is a don't-care coordinate: d is only defined
      mod n, so ||v||^2 should not contain d^2.  Quotienting L by the kernel
      Z*(n*S_D*e_m) gives pi(L) of rank 2m+1 with
          covol(pi L) = covol(L) / (n*S_D)
      and the trivial vector maps to 0 -- it is gone, not merely demoted.
      d is still recoverable: the projected vector still carries every k1_i
      and k2_i, so d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n.  No
      transformation matrix is needed.

  (C) CENTERING.  k1_i is uniform on [0,K1) and k2_i on [0,K2), so the
      uncentered planted coordinates have E[x^2] = n^2/3 instead of n^2/12.
      Replacing A_i by A_i - lam*(K2//2) - (K1//2) recentres both blocks and
      keeps the target a lattice vector (no CVP needed).

Predicted effect (computed before running; see EXP P0):
      ||v_planted||^2 = n^2 * (2m/3 + 4/3)      V0  baseline
                      = n^2 * (m/6  + 1)        V3  projected + centered
  For m=12 that is 9.33 n^2 vs 3.00 n^2, a norm ratio of 1.76.  The Gaussian
  heuristic for this lattice scales as GH ~ K1^{-m/(2m+2)}, so a 1.76x
  shorter target should move the K1 wall outward by
      1.76^((2m+2)/m) = 1.76^2.167 = 3.5x,
  i.e. from K1 ~ 4-6 to K1 ~ 14-21 on the lam*=0.07 curve.

Experiments:
  P0  closed-form norm / Gaussian-heuristic table for V0..V3 (no LLL)
  P1  sanity: all four variants on the historical curves at their logged K1
  P2  the T4 K1 grid, all four variants, both 12-bit curves, 5 seeds
  P3  is the planted vector lambda_1 now?  sv/pv per variant
  P4  17-bit cross-curve sweep at eff = 0.15 and 0.25 (where V0 scored
      3/20 and 0/20 on 2026-07-29), V0 vs V3

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- same as
# glv_hnp_phase2_20bit.py / glv_hnp_phase2_lambda_threshold.py
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
        t, r = t * c % p, r * b % p
    return r

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(4000):
        x = rng.randrange(p)
        rhs = (x * x * x + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, sorted."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Lattice variants
#
#   V0  baseline      (2026-06-15; glv_hnp_phase2_20bit.py:262)   dim 2m+2
#   V1  centered                                                  dim 2m+2
#   V2  projected                                                 dim 2m+1
#   V3  projected + centered                                      dim 2m+1
#
# Column layout, baseline:   [ k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan ]
# Column layout, projected:  [ k1_0..k1_{m-1} |   | k2_0..k2_{m-1} | kannan ]
#
# The projected generating set has 2m+2 rows in 2m+1 columns; it is
# rank-deficient by exactly one (the trivial vector is in the kernel).
# fpylll's LLL handles that and emits a single zero row.
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def build_lattice(sigs, n, lam, k1_bound, k2_bound, project, center):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    k1h = k1_bound // 2 if center else 0
    k2h = k2_bound // 2 if center else 0

    ncols = (2 * m + 1) if project else (2 * m + 2)
    dcol = None if project else m
    # column index of k2 block / kannan
    k2c = m if project else m + 1
    kanc = ncols - 1

    M = [[0] * ncols for _ in range(2 * m + 2)]
    for i in range(m):
        M[i][i] = n * S_K1                        # modular rows
    for i in range(m):                            # the d row
        M[m][i] = sigs[i]['B'] * S_K1
    if dcol is not None:
        M[m][dcol] = S_D
    for i in range(m):                            # the k2 rows
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][k2c + i] = S_K2
    for i in range(m):                            # the Kannan row
        At = (sigs[i]['A'] - lam * k2h - k1h) % n
        M[2 * m + 1][i] = At * S_K1
    M[2 * m + 1][kanc] = S_KANNAN
    return M

def planted_vector(sigs, d_secret, n, lam, k1_bound, k2_bound, project, center):
    """The vector the attack is supposed to find, in the same coordinates."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    k1h = k1_bound // 2 if center else 0
    k2h = k2_bound // 2 if center else 0
    ncols = (2 * m + 1) if project else (2 * m + 2)
    k2c = m if project else m + 1
    v = [0] * ncols
    for i in range(m):
        v[i] = (sigs[i]['k1'] - k1h) * S_K1
        v[k2c + i] = (sigs[i]['k2'] - k2h) * S_K2
    if not project:
        v[m] = d_secret * S_D
    v[ncols - 1] = S_KANNAN
    return v

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Recovery: read k1_0, k2_0 straight out of the short vector and solve for d.
# Works for every variant, projected or not, so all four are scored the same
# way.  (Baseline glv_hnp_phase2_20bit.py:recover_d read d out of the d
# column; that column does not exist in V2/V3.)
# ---------------------------------------------------------------------------

def recover_d(reduced, sigs, n, lam, k1_bound, k2_bound, project, center,
              d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    k1h = k1_bound // 2 if center else 0
    k2h = k2_bound // 2 if center else 0
    ncols = (2 * m + 1) if project else (2 * m + 2)
    k2c = m if project else m + 1
    Binv = modinv(sigs[0]['B'], n)
    for row in reduced:
        last = row[ncols - 1]
        if abs(last) != S_KANNAN: continue
        r = [-x for x in row] if last < 0 else list(row)
        if r[0] % S_K1 != 0 or r[k2c] % S_K2 != 0: continue
        k1_0 = r[0] // S_K1 + k1h
        k2_0 = r[k2c] // S_K2 + k2h
        k_full = (k1_0 + lam * k2_0) % n
        d_cand = (k_full - sigs[0]['A']) * Binv % n
        if d_cand == d_secret:
            return d_cand
    return None

def shortest_nonzero(reduced):
    best = None
    for row in reduced:
        nz = norm(row)
        if nz == 0: continue
        if best is None or nz < best[0]:
            best = (nz, list(row))
    return best

# ---------------------------------------------------------------------------
# One trial
# ---------------------------------------------------------------------------

VARIANTS = [
    ("V0 base",      False, False),
    ("V1 center",    False, True),
    ("V2 project",   True,  False),
    ("V3 proj+cent", True,  True),
]

def run_trial(curve, m, k1_bound, seed, project, center,
              use_bkz=False, bkz_beta=20, want_ratio=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return (False, None) if want_ratio else False
    M = build_lattice(sigs, n, lam, k1_bound, k2_bound, project, center)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    ncols = A.ncols
    reduced = [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]
    ok = recover_d(reduced, sigs, n, lam, k1_bound, k2_bound, project, center,
                   d_secret) is not None
    if not want_ratio:
        return ok
    pv = planted_vector(sigs, d_secret, n, lam, k1_bound, k2_bound,
                        project, center)
    sv = shortest_nonzero(reduced)
    ratio = (sv[0] / norm(pv)) if sv else None
    return ok, ratio

def success_rate(curve, m, k1_bound, seeds, project, center,
                 use_bkz=False, bkz_beta=20):
    w = sum(1 for s in seeds
            if run_trial(curve, m, k1_bound, s, project, center,
                         use_bkz, bkz_beta))
    return w, len(seeds)

# ---------------------------------------------------------------------------
# Curve construction / search (same as glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def build_curve(p, n, seed=12345):
    lam, _ = glv_roots(n)
    if lam is None: return None
    for b_try in range(1, 400):
        G = find_generator(p, b_try, n, seed)
        if G is not None:
            return (p, b_try, n, lam, G)
    return None

def search_curves(lo, hi, count, seed=12345):
    """j=0 GLV curves with prime n = p+1-t, n = 1 mod 3, in [lo, hi)."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < count:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 8 or not sympy.isprime(nc) or nc % 3 != 1:
                        continue
                    cur = build_curve(p, nc, seed)
                    if cur is not None:
                        out.append(cur)
                        break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Closed-form predictions
# ---------------------------------------------------------------------------

def predict(n, m, k1_bound, k2_bound, project, center):
    """(E||v_planted||, Gaussian heuristic) under the uniform model."""
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    # per-coordinate second moment of S*x for x uniform on an interval of
    # width W: W^2/3 uncentered, W^2/12 centered.  S_K1*K1 ~ S_K2*K2 ~ n.
    f = 12.0 if center else 3.0
    e2 = 2 * m * (S_K1 * k1_bound) ** 2 / f          # k1 and k2 blocks
    # (uses S_K1*K1 ~ S_K2*K2 ~ n; exact enough for a prediction)
    e2 += 0.0 if project else n * n / 3.0            # d column
    e2 += S_KANNAN ** 2                              # kannan column
    pv = math.sqrt(e2)

    log_det = m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_KANNAN)
    if not project:
        log_det += math.log(S_D)
    else:
        log_det -= math.log(n * S_D)                 # covol(L)/covol(kernel)
    dim = (2 * m + 1) if project else (2 * m + 2)
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)
    return pv, gh

# ===========================================================================

print("=" * 78)
print("Thread 23 — project out the d-column and recentre the GLV-HNP lattice")
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
print("EXP P0: closed-form planted norm vs Gaussian heuristic (no LLL)")
print("-" * 78)
print("pv = E||v_planted||, gh = Gaussian heuristic for lambda_1.")
print("pv/gh < 1 is the (heuristic) condition for the target to be findable.\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} {'variant':<14} {'dim':>4} "
      f"{'pv':>10} {'gh':>10} {'pv/gh':>7}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for vname, proj, cent in VARIANTS:
        pv, gh = predict(n, m, k1, k2b, proj, cent)
        dim = (2 * m + 1) if proj else (2 * m + 2)
        print(f"{label:<18} {k1:>3} {m:>3} {vname:<14} {dim:>4} "
              f"{pv:>10.0f} {gh:>10.0f} {pv/gh:>7.3f}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P1: all four variants on the historical curves, logged K1 and m")
print("-" * 78)
print(f"{'curve':<18} {'lam*':>6} {'K1':>3} {'m':>3} "
      + " ".join(f"{v[0]:>14}" for v in VARIANTS))
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    cells = []
    for vname, proj, cent in VARIANTS:
        w, t = success_rate(curve, m, k1, SEEDS, proj, cent)
        cells.append(f"{w}/{t}")
    print(f"{label:<18} {lam_star(lam, n):>6.3f} {k1:>3} {m:>3} "
          + " ".join(f"{c:>14}" for c in cells))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P2: the T4 K1 grid (2026-07-29), all four variants, m=12")
print("-" * 78)
print("T4 baseline for reference (m as logged, V0):")
print("  12-bit/2557 (lam*=0.340): 5/5 up to K1=8, 4/5 at 12, 1/5 at 16, 0/5 at 24")
print("  12-bit/2677 (lam*=0.070): 5/5 up to K1=4, 2/5 at 6, 0/5 from K1=8\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
M_GRID = 12
walls = {}
for label, curve, _k1, _m in hist[1:]:
    p, b, n, lam, G = curve
    print(f"{label}  (n={n}, lam*={lam_star(lam, n):.3f}, m={M_GRID}, "
          f"K2={math.isqrt(n)+1})")
    print(f"  {'variant':<14} " + " ".join(f"{k:>5}" for k in K1_GRID))
    for vname, proj, cent in VARIANTS:
        row, wall = [], None
        for k1 in K1_GRID:
            w, t = success_rate(curve, M_GRID, k1, SEEDS, proj, cent)
            row.append(f"{w}/{t}")
            if w == t:
                wall = k1
        walls[(label, vname)] = wall
        print(f"  {vname:<14} " + " ".join(f"{c:>5}" for c in row))
    print(f"  last K1 with 5/5:  "
          + ",  ".join(f"{v[0]}={walls[(label, v[0])]}" for v in VARIANTS))
    print()


# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP P3: is the planted vector lambda_1 now?  (sv/pv, one seed)")
print("-" * 78)
print("sv = shortest nonzero vector after LLL, pv = planted vector.")
print("sv/pv >= 1 means nothing in the reduced basis is shorter than the target.\n")
print(f"{'curve':<18} {'K1':>3} " + " ".join(f"{v[0]:>14}" for v in VARIANTS))
for label, curve, k1, m in hist:
    for K1 in sorted({k1, 8, 16}):
        cells = []
        for vname, proj, cent in VARIANTS:
            ok, ratio = run_trial(curve, m, K1, 42, proj, cent, want_ratio=True)
            cells.append(f"{ratio:.3f}{'*' if ok else ' '}" if ratio else "n/a")
        print(f"{label:<18} {K1:>3} " + " ".join(f"{c:>14}" for c in cells))
print("(* = d recovered on that trial)")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P4: 17-bit cross-curve sweep at eff = K1*K2/n in {0.15, 0.25}")
print("-" * 78)
print("2026-07-29 T3 measured V0 at eff=0.15 -> 3/20 curves 5/5, "
      "eff=0.25 -> 0/20.\n")

curves17 = search_curves(2 ** 16, 2 ** 17, 20)
print(f"found {len(curves17)} 17-bit j=0 GLV curves\n")
for eff in (0.15, 0.25):
    tot = {v[0]: 0 for v in VARIANTS}
    per = {v[0]: 0 for v in VARIANTS}
    for cur in curves17:
        p, b, n, lam, G = cur
        k2b = math.isqrt(n) + 1
        k1 = max(2, int(round(eff * n / k2b)))
        for vname, proj, cent in VARIANTS:
            if vname not in ("V0 base", "V3 proj+cent"): continue
            w, t = success_rate(cur, 12, k1, SEEDS, proj, cent)
            tot[vname] += w
            per[vname] += 1 if w == t else 0
    print(f"eff={eff}:  "
          + ",  ".join(f"{v}: {per[v]}/{len(curves17)} curves 5/5 "
                       f"({tot[v]}/{len(curves17)*len(SEEDS)} trials)"
                       for v in ("V0 base", "V3 proj+cent")))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P5: predicted critical K1 (largest K1 with pv/gh < 1) vs observed")
print("-" * 78)
print("Prediction uses only the closed form of EXP P0 -- no LLL.  The observed")
print("column is the last K1 with 5/5 in EXP P2 (m=12).\n")

def critical_k1(n, m, k2_bound, project, center, hi=4096):
    best = None
    k1 = 2
    while k1 <= hi:
        pv, gh = predict(n, m, k1, k2_bound, project, center)
        if pv / gh < 1.0:
            best = k1
        else:
            break
        k1 += 1
    return best

print(f"{'curve':<18} {'variant':<14} {'K1_crit pred':>12} {'K1 obs (5/5)':>13}")
for label, curve, _k1, _m in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for vname, proj, cent in VARIANTS:
        kc = critical_k1(n, M_GRID, k2b, proj, cent)
        print(f"{label:<18} {vname:<14} {str(kc):>12} "
              f"{str(walls.get((label, vname))):>13}")

print("\nAsymptotic form.  With S_K1*K1 ~ S_K2*K2 ~ n and S_D = 1,")
print("  ||v||^2 = n^2 (2m/3 + 4/3)   V0      vs   n^2 (m/6 + 1)   V3")
print("  covol   = n^(3m+1)/(K1 K2)^m V0      vs   n^(3m)/(K1 K2)^m V3")
print("so pv/gh = 1 gives  K1*K2 = n^((m-1)/m) * C^((2m+2)/m)  and the V3/V0")
print("ratio of critical K1 is  (pv_V0/pv_V3)^((2m+2)/m):")
for label, curve, _k1, _m in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    pv0, _ = predict(n, M_GRID, 8, k2b, False, False)
    pv3, _ = predict(n, M_GRID, 8, k2b, True, True)
    r = (pv0 / pv3) ** ((2 * M_GRID + 2) / M_GRID)
    w0, w3 = walls.get((label, "V0 base")), walls.get((label, "V3 proj+cent"))
    obs = (w3 / w0) if (w0 and w3) else float('nan')
    print(f"  {label:<18} pv_V0/pv_V3 = {pv0/pv3:.3f}  ->  predicted "
          f"{r:.2f}x,  observed {obs:.2f}x")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
