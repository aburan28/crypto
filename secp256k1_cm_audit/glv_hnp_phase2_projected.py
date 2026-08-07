"""
GLV-HNP Phase 2, Thread 23: remove the trivial vector n*S_D*e_m from the
Phase-2 lattice and test whether that moves the K1 wall.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, "T5"):
  The shortest vector of the Phase-2 GLV lattice is always the trivial
  vector n*S_D*e_m (all its energy in the d column).  It is 2-3x shorter
  than the planted vector on every curve tested, so the planted vector is
  never lambda_1 and recovery is a BDD/coset condition.  That entry
  asserted:

      "no choice of S_D removes it -- both vectors scale linearly in S_D"

  and proposed Thread 23: project the lattice along e_m (quotient out the
  trivial direction) and re-measure the K1 wall.

Falsifier stated on 2026-07-29:
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4
   moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
   reformulation is a real improvement; if the wall stays at K1 ~ 4-6, then
   the wall is information-theoretic and Phase 2 is at its ceiling."

Experiments:
  T23.1  S_D sweep in the ORIGINAL lattice.  Tests the "both vectors scale
         linearly in S_D" claim directly.  (They do not: only the d
         coordinate of the planted vector scales, so sv/pv is a non-trivial
         function of S_D.)
  T23.2  PROJECTED lattice: delete column m from the basis.  This is exactly
         the coordinate projection pi: L -> Z^(dim-1) whose kernel on L is
         Z*(n*S_D*e_m).  The trivial vector maps to 0 and disappears; d is
         recovered algebraically from (k1_i, k2_i) instead of being read off
         a lattice coordinate.  Measure sv/pv.
  T23.3  K1 wall grid (the T4 grid of 2026-07-29) for original vs projected,
         on both historical 12-bit curves.  This is the falsifier.
  T23.4  m-sweep at the wall (K1=8, lam*=0.07) -- does the projection let
         extra signatures rescue the attack where T4b showed they do not?
  T23.5  eff-scaling on fresh 17-bit curves: where is the eff = K1*K2/n wall
         for each formulation, vs the information-theoretic bound
         eff_max = n^(-1/m)?

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- same as glv_hnp_phase2_20bit.py
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

def find_generator(p, b, n):
    rng = random.Random(12345)
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
# CM curve search (j=0), same routines as glv_hnp_phase2_20bit.py
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0, f"GLV eigenvalue check failed for n={n}"
    lam = min(r1, r2)
    return lam, n - 1 - lam

def find_curves(bits, want, seed_skip=0):
    """Find `want` prime-order j=0 GLV curves with n of ~`bits` bits."""
    out = []
    p = sympy.nextprime(2 ** (bits - 1))
    skipped = 0
    while p < 2 ** bits and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 8 or not sympy.isprime(n_cand) or n_cand % 3 != 1:
                        continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    if skipped < seed_skip:
                        skipped += 1
                        continue
                    found_b = None
                    for b_try in range(1, 60):
                        G = find_generator(p, b_try, n_cand)
                        if G is not None:
                            found_b = b_try
                            break
                    if found_b is None:
                        continue
                    out.append((p, found_b, n_cand, lam, G))
                    break
        p = sympy.nextprime(p)
    return out

# ---------------------------------------------------------------------------
# Signature generation (identical to glv_hnp_phase2_20bit.py:gen_signatures)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
# Lattices
# ---------------------------------------------------------------------------
#
# ORIGINAL (glv_hnp_phase2_20bit.py:build_glv_lattice), dim = 2m+2:
#   cols          0..m-1        m       m+1..2m       2m+1
#   meaning       congruence    d       k2            Kannan
#   row i<m       n*S_K1 e_i
#   row m         B_i*S_K1      S_D
#   row m+1+i     -lam*S_K1 e_i           S_K2 e_i
#   row 2m+1      A_i*S_K1                              S_KANNAN
#
# planted = (k1_i*S_K1, d*S_D, k2_i*S_K2, S_KANNAN)
# trivial = n*S_D*e_m   (norm n*S_D; = n*row_m - sum_i B_i*row_i)
#
# PROJECTED: delete column m.  dim = 2m+1.  The trivial vector maps to 0.
#   cols          0..m-1        m..2m-1     2m
#   planted = (k1_i*S_K1, k2_i*S_K2, S_KANNAN)

def scales(n, k1_bound, k2_bound, S_D=1):
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    return S_K1, S_D, S_K2, S_KANNAN

def build_original(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
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

def build_projected(sigs, n, lam, k1_bound, k2_bound):
    """Original lattice with column m (the d column) deleted, and the now
    redundant generator row m kept -- the resulting (2m+2) x (2m+1) matrix
    generates pi(L).  fpylll accepts a rectangular generating set and LLL
    returns a basis (with a zero row for the rank defect)."""
    M = build_original(sigs, n, lam, k1_bound, k2_bound, S_D=1)
    return [row[:len(sigs)] + row[len(sigs) + 1:] for row in M]

def planted_original(sigs, n, lam, k1_bound, k2_bound, d_secret, S_D=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def planted_projected(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_original(rows, m, n, S_KAN, d_secret, S_D=1):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % S_D:           # the d column carries the factor S_D
            continue
        d_cand = (val // S_D) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

def recover_projected(rows, sigs, m, n, lam, S_K1, S_K2, S_KAN, d_secret):
    """From a projected-lattice vector (k1_i*S_K1, k2_i*S_K2, +-S_KANNAN),
    read off (k1_0, k2_0) and solve A_0 + d*B_0 = k1_0 + lam*k2_0 (mod n)."""
    dim = 2 * m + 1
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        for i in range(m):
            a = sign * row[i]
            b = sign * row[m + i]
            if a % S_K1 or b % S_K2:
                continue
            k1 = a // S_K1
            k2 = b // S_K2
            k_full = (k1 + lam * k2) % n
            B = sigs[i]['B']
            if math.gcd(B, n) != 1:
                continue
            d_cand = (k_full - sigs[i]['A']) * modinv(B, n) % n
            if d_cand and d_cand == d_secret:
                return d_cand
    return None

# ---------------------------------------------------------------------------
# One trial
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, mode="orig", S_D=1,
          use_bkz=False, beta=20, want_ratio=False):
    """mode in {"orig", "proj"}.  Returns (recovered, sv/pv) ."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, None
    S_K1, S_Dv, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)

    if mode == "orig":
        M = build_original(sigs, n, lam, k1_bound, k2_bound, S_D)
        pv = planted_original(sigs, n, lam, k1_bound, k2_bound, d_secret, S_D)
        dim = 2 * m + 2
    else:
        M = build_projected(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_projected(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1

    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]

    if mode == "orig":
        ok = recover_original(rows, m, n, S_KAN, d_secret, S_Dv) is not None
    else:
        ok = recover_projected(rows, sigs, m, n, lam, S_K1, S_K2, S_KAN,
                               d_secret) is not None

    if want_ratio == "coset":
        # norm of the best Kannan-coset representative LLL actually produced,
        # relative to the planted vector.  Recovery needs this to be >= 1.
        cos = [r for r in rows if abs(r[dim - 1]) == S_KAN]
        return ok, (min(norm(r) for r in cos) / norm(pv)) if cos else None

    ratio = None
    if want_ratio:
        nz = [r for r in rows if any(r)]
        if nz:
            sv = min(nz, key=norm)
            ratio = norm(sv) / norm(pv)
            if want_ratio == "energy":
                # fraction of ||sv||^2 in each block
                if mode == "orig":
                    blocks = {"k1": sv[:m], "d": sv[m:m + 1],
                              "k2": sv[m + 1:2 * m + 1], "kan": sv[2 * m + 1:]}
                else:
                    blocks = {"k1": sv[:m], "d": [],
                              "k2": sv[m:2 * m], "kan": sv[2 * m:]}
                tot = norm(sv) ** 2 or 1.0
                energy = {k: norm(v) ** 2 / tot for k, v in blocks.items()}
                return ok, (ratio, energy)
    return ok, ratio

def rate(curve, m, k1_bound, seeds, mode="orig", S_D=1, use_bkz=False, beta=20):
    n = curve[2]
    wins = 0
    for seed in seeds:
        d = random.Random(seed + 7777).randint(1, n - 1)
        ok, _ = trial(curve, m, d, k1_bound, seed, mode, S_D, use_bkz, beta)
        wins += bool(ok)
    return wins

# ---------------------------------------------------------------------------
# Historical curves (2026-06-15 / 2026-07-29)
# ---------------------------------------------------------------------------

C_2557 = (2557, 2, 2659, 1755, find_generator(2557, 2, 2659))   # lam* = 0.340
C_2677 = (2677, 2, 2647, 185,  find_generator(2677, 2, 2647))   # lam* = 0.070
C_211  = (211,  2, 199,  106,  find_generator(211, 2, 199))     # lam* = 0.467

HIST = [("8-bit/199",   C_211,  0.467, 2),
        ("12-bit/2557", C_2557, 0.340, 8),
        ("12-bit/2677", C_2677, 0.070, 8)]

SEEDS5 = [42, 1234, 9999, 31337, 271828]

print("=" * 74)
print("Thread 23 -- projecting out the trivial vector n*S_D*e_m")
print("=" * 74)
for lbl, c, ls, _ in HIST:
    assert c[4] is not None, f"no generator for {lbl}"
    print(f"  {lbl}: p={c[0]} n={c[2]} lam={c[3]} lam*={ls}")

# ---------------------------------------------------------------------------
# T23.1  S_D sweep in the ORIGINAL lattice
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("T23.1  S_D sweep, ORIGINAL lattice, m=12")
print("  claim under test (2026-07-29): 'no choice of S_D removes the trivial")
print("  vector -- both vectors scale linearly in S_D'")
print("=" * 74)

SD_GRID = [1, 2, 4, 8, 16, 64, 256, 1024]
print(f"{'curve':<14}{'S_D':>6}  {'sv/pv':>7}  {'recov':>7}")
t231 = {}
for lbl, curve, ls, k1 in HIST[1:]:
    for S_D in SD_GRID:
        n = curve[2]
        d = random.Random(42 + 7777).randint(1, n - 1)
        _, ratio = trial(curve, 12, d, k1, 42, "orig", S_D, want_ratio=True)
        w = rate(curve, 12, k1, SEEDS5, "orig", S_D)
        t231[(lbl, S_D)] = (ratio, w)
        print(f"{lbl:<14}{S_D:>6}  {ratio:>7.3f}  {w:>4}/{len(SEEDS5)}")
    print()

# ---------------------------------------------------------------------------
# T23.2  PROJECTED lattice: sv/pv
# ---------------------------------------------------------------------------
print("=" * 74)
print("T23.2  sv/pv, ORIGINAL (S_D=1) vs PROJECTED, m=12")
print("=" * 74)
print(f"{'curve':<14}{'K1':>4}  {'sv/pv orig':>11}  {'sv/pv proj':>11}  "
      f"{'rec orig':>9}  {'rec proj':>9}")
for lbl, curve, ls, k1 in HIST:
    n = curve[2]
    d = random.Random(42 + 7777).randint(1, n - 1)
    _, r_o = trial(curve, 12, d, k1, 42, "orig", 1, want_ratio=True)
    _, r_p = trial(curve, 12, d, k1, 42, "proj", 1, want_ratio=True)
    w_o = rate(curve, 12, k1, SEEDS5, "orig")
    w_p = rate(curve, 12, k1, SEEDS5, "proj")
    print(f"{lbl:<14}{k1:>4}  {r_o:>11.3f}  {r_p:>11.3f}  "
          f"{w_o:>6}/{len(SEEDS5)}  {w_p:>6}/{len(SEEDS5)}")

# ---------------------------------------------------------------------------
# T23.2b  where does the shortest vector's energy sit?
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("T23.2b  energy of the shortest vector by block (fraction of ||sv||^2)")
print("  2026-07-29 T5 measured this for ORIGINAL only: d-block = 1.000")
print("=" * 74)
print(f"{'curve':<14}{'mode':>6}{'S_D':>5}  {'sv/pv':>7}  {'k1':>7}{'d':>8}"
      f"{'k2':>8}{'kan':>8}")
for lbl, curve, ls, k1 in HIST[1:]:
    n = curve[2]
    d = random.Random(42 + 7777).randint(1, n - 1)
    for mode, S_D in (("orig", 1), ("orig", 8), ("proj", 1)):
        _, (r, e) = trial(curve, 12, d, k1, 42, mode, S_D, want_ratio="energy")
        print(f"{lbl:<14}{mode:>6}{S_D:>5}  {r:>7.3f}  {e['k1']:>7.3f}"
              f"{e['d']:>8.3f}{e['k2']:>8.3f}{e['kan']:>8.3f}")

# ---------------------------------------------------------------------------
# T23.3  the falsifier: the T4 K1 grid, original vs projected
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("T23.3  K1 wall (the T4 grid of 2026-07-29), m=12, 5 seeds")
print("  FALSIFIER: does the wall on the lam*=0.07 curve move outward?")
print("=" * 74)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
t233 = {}
for lbl, curve, ls, _ in HIST[1:]:
    print(f"\n{lbl}  (lam*={ls})")
    hdr = "".join(f"{k:>5}" for k in K1_GRID)
    print(f"  {'K1':<10}{hdr}")
    for mode in ("orig", "proj"):
        cells = []
        for k1 in K1_GRID:
            w = rate(curve, 12, k1, SEEDS5, mode)
            t233[(lbl, mode, k1)] = w
            cells.append(f"{w:>5}")
        print(f"  {mode:<10}" + "".join(cells))
    # BKZ(20) on the projected lattice, at and beyond the LLL wall
    cells = []
    for k1 in K1_GRID:
        if t233[(lbl, "proj", k1)] == len(SEEDS5):
            cells.append(f"{'-':>5}")
        else:
            w = rate(curve, 12, k1, SEEDS5, "proj", use_bkz=True, beta=20)
            cells.append(f"{w:>5}")
    print(f"  {'proj+BKZ20':<10}" + "".join(cells))

# per-trial agreement: is "proj" the same predicate as "orig", trial by trial?
agree = disagree = 0
for lbl, curve, ls, _ in HIST:
    n = curve[2]
    for k1 in K1_GRID:
        for seed in SEEDS5:
            d = random.Random(seed + 7777).randint(1, n - 1)
            o, _ = trial(curve, 12, d, k1, seed, "orig")
            p_, _ = trial(curve, 12, d, k1, seed, "proj")
            if bool(o) == bool(p_):
                agree += 1
            else:
                disagree += 1
print(f"\n  per-trial agreement orig vs proj: {agree}/{agree + disagree} "
      f"({disagree} disagreements)")

# ---------------------------------------------------------------------------
# T23.4  m-sweep at the wall (T4b of 2026-07-29 was orig-only)
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("T23.4  m-sweep at K1=8 on the lam*=0.07 curve (T4b said orig: 0,0,1,0,1)")
print("=" * 74)
M_GRID = [8, 12, 16, 24, 32]
print(f"  {'m':<10}" + "".join(f"{m:>5}" for m in M_GRID))
for mode in ("orig", "proj"):
    cells = []
    for m in M_GRID:
        cells.append(f"{rate(C_2677, m, 8, SEEDS5, mode):>5}")
    print(f"  {mode:<10}" + "".join(cells))

# ---------------------------------------------------------------------------
# T23.5  eff wall on fresh 17-bit curves
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("T23.5  eff = K1*K2/n wall on fresh 17-bit curves, m=12, 3 seeds")
print("  information-theoretic ceiling: eff_max = n^(-1/m)")
print("=" * 74)
CURVES17 = find_curves(17, 6)
SEEDS3 = [42, 1234, 9999]
EFF_GRID = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
print(f"{'p':>7}{'n':>8}{'lam*':>8}   " + "".join(f"{e:>11.2f}" for e in EFF_GRID))
walls = {"orig": [], "proj": []}
for curve in CURVES17:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    lstar = min(lam, n - lam) / n
    eff_max = n ** (-1.0 / 12)
    for mode in ("orig", "proj"):
        cells, wall = [], None
        for e in EFF_GRID:
            k1 = max(2, round(e * n / k2))
            w = rate(curve, 12, k1, SEEDS3, mode)
            cells.append(f"{w:>4}/3 (K1={k1})".rjust(11))
            if w < len(SEEDS3) and wall is None:
                wall = e
        walls[mode].append(wall if wall is not None else 999)
        tag = f"{p:>7}{n:>8}{lstar:>8.3f} {mode[:4]}" if mode == "orig" \
              else f"{'':>7}{'':>8}{'':>8} {mode[:4]}"
        print(tag + "".join(cells))
    print(f"{'':>23} eff_max(IT) = {eff_max:.3f}")

# ---------------------------------------------------------------------------
# T23.6  the coset ratio kappa as a per-instance predictor
# ---------------------------------------------------------------------------
print("\n" + "=" * 74)
print("T23.6  kappa = ||best Kannan-coset row after LLL|| / ||v_planted||")
print("  T5 (2026-07-29) said recovery is a coset/BDD condition, not SVP.")
print("  If so, kappa >= 1 should be the exact success criterion.")
print("=" * 74)

obs = []   # (kappa, success)
for lbl, curve, ls, _ in HIST:
    n = curve[2]
    for k1 in K1_GRID:
        for seed in SEEDS5:
            d = random.Random(seed + 7777).randint(1, n - 1)
            ok, kap = trial(curve, 12, d, k1, seed, "proj", want_ratio="coset")
            if kap is not None:
                obs.append((kap, bool(ok)))
for curve in CURVES17:
    n = curve[2]
    k2 = math.isqrt(n) + 1
    for e in EFF_GRID:
        k1 = max(2, round(e * n / k2))
        for seed in SEEDS3:
            d = random.Random(seed + 7777).randint(1, n - 1)
            ok, kap = trial(curve, 12, d, k1, seed, "proj", want_ratio="coset")
            if kap is not None:
                obs.append((kap, bool(ok)))

pos = [k for k, s in obs if s]
neg = [k for k, s in obs if not s]
print(f"  n={len(obs)} trials, {len(pos)} success / {len(neg)} failure")
if pos and neg:
    print(f"  kappa | success: min={min(pos):.4f} max={max(pos):.4f}")
    print(f"  kappa | failure: min={min(neg):.4f} max={max(neg):.4f}")
    auc = sum((a > b) + 0.5 * (a == b) for a in pos for b in neg) / (len(pos) * len(neg))
    print(f"  AUC(kappa) = {auc:.4f}")
    # exact rule kappa >= 1
    tp = sum(1 for k in pos if k >= 1.0)
    fp = sum(1 for k in neg if k >= 1.0)
    print(f"  rule 'kappa >= 1': {tp}/{len(pos)} successes captured, "
          f"{fp}/{len(neg)} failures misfired  "
          f"(accuracy {(tp + len(neg) - fp) / len(obs):.4f})")
    best = max(((sum(1 for k in pos if k >= t) + sum(1 for k in neg if k < t)) / len(obs), t)
               for t in sorted(set(k for k, _ in obs)))
    print(f"  best threshold: kappa >= {best[1]:.4f} -> accuracy {best[0]:.4f}")
    print(f"  majority baseline: {max(len(pos), len(neg)) / len(obs):.4f}")

print("\n" + "=" * 74)
print("SUMMARY")
print("=" * 74)
for lbl, _, ls, _ in HIST[1:]:
    def wall_of(mode):
        for k1 in K1_GRID:
            if t233[(lbl, mode, k1)] < len(SEEDS5):
                return k1
        return f">{K1_GRID[-1]}"
    print(f"  {lbl} (lam*={ls}): K1 wall  orig={wall_of('orig')}  "
          f"proj={wall_of('proj')}")
ow = [w for w in walls['orig'] if w != 999]
pw = [w for w in walls['proj'] if w != 999]
if ow and pw:
    print(f"  17-bit mean eff wall: orig={sum(ow)/len(ow):.3f}  "
          f"proj={sum(pw)/len(pw):.3f}")
print("\nDone.")
