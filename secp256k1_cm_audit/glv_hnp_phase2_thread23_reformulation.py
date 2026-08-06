"""
GLV-HNP Phase 2, Thread 23: make the planted vector lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, exp T5):
  The Phase-2 lattice `build_glv_lattice` (glv_hnp_phase2_20bit.py:262) always
  has shortest vector  n * S_D * e_m  -- the "trivial" vector living entirely
  in the d column.  With S_D = 1 its norm is exactly n, while

      E||v_planted||^2 = n^2 * (2m/3 + S_D^2/3 + 1)                      (*)

  so the planted vector is ~2-3x longer and is NEVER lambda_1.  Recovery is
  therefore a BDD/coset condition, not SVP.  The 07-29 entry concluded:

      "no choice of S_D removes it -- both vectors scale linearly in S_D."

  That conclusion is WRONG, and (*) shows why: only the d COLUMN of the
  planted vector scales with S_D, not the whole vector.  Setting
  ||v_planted||^2 < ||n*S_D*e_m||^2 = n^2 S_D^2 gives

      2m/3 + S_D^2/3 + 1 < S_D^2   <=>   S_D^2 > m + 3/2                 (**)

  i.e. S_D > sqrt(m + 1.5) kills the trivial vector.  For m = 8 that is
  S_D >= 4; for m = 12, S_D >= 4.  This is a one-line change to `scales()`.

Two reformulations are tested here against the 07-29 baseline:

  V1(S_D)  the existing lattice with S_D swept (dim 2m+2).  Prediction (**).
  V2       d eliminated algebraically by pivoting on signature 0 (dim 2m+1).
           From k_full_i = A_i + B_i*d (mod n),
               t_i = B_i/B_0,   c_i = A_i - t_i*A_0   (mod n)
               k_full_i = t_i * k_full_0 + c_i        (mod n)
           with k_full_j = k1_j + lam*k2_j.  d never appears, so the trivial
           vector cannot exist.  d is read back off the recovered (k1_0,k2_0).

Falsifier stated on 07-29:
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  R1  S_D sweep: sv/pv and recovery rate vs prediction (**).
  R2  V2 correctness + sv/pv.
  R3  the T4 K1-wall grid re-run for V1(S_D=1), V1(S_D=opt), V2.
  R4  Gaussian-heuristic model for the eff = K1*K2/n wall; check against the
      07-29 T3/T4 numbers.

Run: python3 glv_hnp_phase2_thread23_reformulation.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic -- verbatim from glv_hnp_phase2_lambda_threshold.py
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
    m_, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m_ - i - 1), p)
        m_, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (norm form of Z[omega])."""
    for a in range(1, int(math.isqrt(4 * p // 3)) + 2):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        r = math.isqrt(disc)
        if r * r == disc and (a + r) % 2 == 0:
            b = (a + r) // 2
            if a * a - a * b + b * b == p:
                return (a, b)
    return None

def j0_traces(a, b):
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, n - 1 - r1), max(r1, n - 1 - r1))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures -- verbatim from glv_hnp_phase2_lambda_threshold.py:230
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# V1: the 07-29 lattice, with S_D exposed as a parameter
# ---------------------------------------------------------------------------

def scales_v1(n, k1_bound, k2_bound, s_d=1):
    return (n // k1_bound, s_d, max(1, n // k2_bound), n)

def build_v1(sigs, n, lam, k1_bound, k2_bound, s_d=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales_v1(n, k1_bound, k2_bound, s_d)
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

def planted_v1(sigs, d_secret, n, k1_bound, k2_bound, s_d=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales_v1(n, k1_bound, k2_bound, s_d)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_v1(reduced, m, n, S_D, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        dv = sign * row[m]
        if dv % S_D != 0: continue
        d_cand = (dv // S_D) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# V2: d eliminated by pivoting on signature 0.  dim = 2m+1.
#
# columns:  0 .. m-2      k1_i for i=1..m-1     (scale S_K1)
#           m-1           k1_0                  (scale S_K1)
#           m             k2_0                  (scale S_K2)
#           m+1 .. 2m-1   k2_i for i=1..m-1     (scale S_K2)
#           2m            Kannan                (scale S_KAN)
#
# rows:     0 .. m-2      n * S_K1 * e_{i-1}          (mod-n rows)
#           m-1           coefficient of k1_0
#           m             coefficient of k2_0
#           m+1 .. 2m-1   coefficient of k2_j, j=1..m-1
#           2m            Kannan row (constants c_i)
# ---------------------------------------------------------------------------

def v2_reduction(sigs, n):
    """t_i = B_i/B_0, c_i = A_i - t_i*A_0 (mod n), for i = 1..m-1."""
    B0inv = modinv(sigs[0]['B'] % n, n)
    t, c = [], []
    for i in range(1, len(sigs)):
        ti = sigs[i]['B'] * B0inv % n
        ci = (sigs[i]['A'] - ti * sigs[0]['A']) % n
        t.append(ti)
        c.append(ci)
    return t, c

def build_v2(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KAN = n
    t, c = v2_reduction(sigs, n)
    lam = lam % n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m - 1):                       # mod-n rows
        M[i][i] = n * S_K1
    for i in range(m - 1):                       # k1_0 row
        M[m - 1][i] = t[i] * S_K1
    M[m - 1][m - 1] = S_K1
    for i in range(m - 1):                       # k2_0 row
        M[m][i] = (lam * t[i] % n) * S_K1
    M[m][m] = S_K2
    for j in range(1, m):                        # k2_j rows
        M[m + j][j - 1] = -lam * S_K1
        M[m + j][m + j] = S_K2
    for i in range(m - 1):                       # Kannan row
        M[2 * m][i] = c[i] * S_K1
    M[2 * m][2 * m] = S_KAN
    return M, (S_K1, S_K2, S_KAN)

def planted_v2(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KAN = n
    v = [0] * (2 * m + 1)
    for i in range(1, m):
        v[i - 1] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[m - 1] = sigs[0]['k1'] * S_K1
    v[m] = sigs[0]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_v2(reduced, sigs, n, lam, S_K1, S_K2, S_KAN, d_secret):
    """Read (k1_0, k2_0) off the Kannan row, rebuild k_full_0, then d."""
    m = len(sigs)
    dim = 2 * m + 1
    B0inv = modinv(sigs[0]['B'] % n, n)
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        a, b = sign * row[m - 1], sign * row[m]
        if a % S_K1 or b % S_K2: continue
        k10, k20 = a // S_K1, b // S_K2
        k_full0 = (k10 + lam * k20) % n
        d_cand = (k_full0 - sigs[0]['A']) * B0inv % n
        if d_cand and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def reduce_rows(M, dim, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def trial(curve, m, d_secret, k1_bound, seed, variant, s_d=1,
          use_bkz=False, beta=20):
    """-> (recovered, ||planted||, ||shortest reduced row||, trivial_is_sv)"""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if variant == 'v1':
        M = build_v1(sigs, n, lam, k1_bound, k2_bound, s_d)
        dim = 2 * m + 2
        pv = planted_v1(sigs, d_secret, n, k1_bound, k2_bound, s_d)
        red = reduce_rows(M, dim, use_bkz, beta)
        S_K1, S_D, S_K2, S_KAN = scales_v1(n, k1_bound, k2_bound, s_d)
        ok = recover_v1(red, m, n, S_D, S_KAN, d_secret)
    else:
        M, (S_K1, S_K2, S_KAN) = build_v2(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1
        pv = planted_v2(sigs, n, k1_bound, k2_bound)
        red = reduce_rows(M, dim, use_bkz, beta)
        ok = recover_v2(red, sigs, n, lam, S_K1, S_K2, S_KAN, d_secret)
    nz = [r for r in red if any(r)]
    sv = min(nz, key=norm)
    # is the shortest vector the trivial n*S_D*e_m one?
    if variant == 'v1':
        triv = (abs(sv[m]) == n * s_d and all(x == 0 for j, x in enumerate(sv)
                                              if j != m))
    else:
        triv = False
    return ok, norm(pv), norm(sv), triv

def rate(curve, m, k1_bound, seeds, variant, s_d=1, use_bkz=False, beta=20):
    wins, ratios, trivs = 0, [], 0
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        r = trial(curve, m, d_trial, k1_bound, seed, variant, s_d, use_bkz, beta)
        if r is None:
            continue
        ok, pn, sn, tv = r
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
        trivs += bool(tv)
    mean = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mean, trivs

def expected_planted_norm(m, n, k1_bound, k2_bound, variant, s_d=1):
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    if variant == 'v1':
        return math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                         + (n * s_d) ** 2 / 3.0
                         + m * (k2_bound * S_K2) ** 2 / 3.0
                         + n ** 2)
    return math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                     + m * (k2_bound * S_K2) ** 2 / 3.0
                     + n ** 2)

def gaussian_heuristic(m, n, k1_bound, k2_bound, variant, s_d=1):
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    if variant == 'v1':
        dim = 2 * m + 2
        logdet = (m * math.log(n * S_K1) + math.log(s_d)
                  + m * math.log(S_K2) + math.log(n))
    else:
        dim = 2 * m + 1
        logdet = ((m - 1) * math.log(n * S_K1) + math.log(S_K1)
                  + m * math.log(S_K2) + math.log(n))
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)


# ===========================================================================
print("=" * 78)
print("Thread 23 -- reformulating the Phase-2 lattice so the target is lambda_1")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [                                 # label, p, b, n, lam, K1, m
    ("8-bit/199",   211,  2, 199,  106,  2,  6),
    ("12-bit/2557", 2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677", 2677, 2, 2647, 185,  8,  10),
]
curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    curves.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP R1: S_D sweep.  Prediction (**): trivial vector dies for S_D^2 > m+1.5")
print("-" * 78)
S_D_GRID = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]
for label, cur, k1, m in curves:
    n = cur[2]
    pred = math.sqrt(m + 1.5)
    print(f"\n{label}: n={n}, lam*={lam_star(cur[3], n):.4f}, K1={k1}, m={m}, "
          f"predicted S_D crit = sqrt(m+1.5) = {pred:.2f}")
    print(f"  {'S_D':>4} | {'sv/pv':>7} | {'trivial=sv':>10} | {'recovered':>9} "
          f"| {'E|pv|/n':>8} | {'GH/n':>7}")
    print("  " + "-" * 62)
    for s_d in S_D_GRID:
        w, tot, mean, trivs = rate(cur, m, k1, SEEDS, 'v1', s_d=s_d)
        k2b = math.isqrt(n) + 1
        epv = expected_planted_norm(m, n, k1, k2b, 'v1', s_d) / n
        gh = gaussian_heuristic(m, n, k1, k2b, 'v1', s_d) / n
        flag = "  <-- pred" if abs(s_d - pred) < 1.0 else ""
        print(f"  {s_d:>4} | {mean:>7.3f} | {trivs:>4}/{tot:<5} | {w:>4}/{tot:<4} "
              f"| {epv:>8.2f} | {gh:>7.2f}{flag}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP R2: V2 (d eliminated, dim 2m+1) -- correctness and sv/pv")
print("-" * 78)
print(f"  {'curve':<14} | {'variant':<10} | {'dim':>4} | {'sv/pv':>7} | "
      f"{'recovered':>9} | {'E|pv|/n':>8} | {'GH/n':>7}")
print("  " + "-" * 74)
for label, cur, k1, m in curves:
    n = cur[2]; k2b = math.isqrt(n) + 1
    for variant, s_d, dim in [('v1', 1, 2 * m + 2),
                              ('v1', int(math.ceil(math.sqrt(m + 1.5))), 2 * m + 2),
                              ('v2', 1, 2 * m + 1)]:
        w, tot, mean, trivs = rate(cur, m, k1, SEEDS, variant, s_d=s_d)
        epv = expected_planted_norm(m, n, k1, k2b, variant, s_d) / n
        gh = gaussian_heuristic(m, n, k1, k2b, variant, s_d) / n
        name = f"{variant}(S_D={s_d})" if variant == 'v1' else 'v2'
        print(f"  {label:<14} | {name:<10} | {dim:>4} | {mean:>7.3f} | "
              f"{w:>4}/{tot:<4} | {epv:>8.2f} | {gh:>7.2f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP R3: the T4 K1-wall grid, re-run for each variant")
print("       (07-29 falsifier: does the wall move outward?)")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
for label, cur, _k1, m in curves[1:]:              # the two 12-bit curves
    n = cur[2]
    print(f"\n{label}: lam*={lam_star(cur[3], n):.4f}, m={m}, n={n}")
    hdr = "  ".join(f"K1={k:<3}" for k in K1_GRID)
    print(f"  {'variant':<12} | {hdr}")
    print("  " + "-" * (14 + len(hdr)))
    sopt = int(math.ceil(math.sqrt(m + 1.5)))
    for variant, s_d in [('v1', 1), ('v1', sopt), ('v2', 1)]:
        cells = []
        for k1 in K1_GRID:
            w, tot, mean, _ = rate(cur, m, k1, SEEDS, variant, s_d=s_d)
            cells.append(f"{w}/{tot} ".ljust(7))
        name = f"v1(S_D={s_d})" if variant == 'v1' else 'v2'
        print(f"  {name:<12} | {'  '.join(cells)}")
    # eff for reference
    k2b = math.isqrt(n) + 1
    effs = "  ".join(f"{k1*k2b/n:<5.2f} " for k1 in K1_GRID)
    print(f"  {'eff=K1*K2/n':<12} | {effs}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP R4: Gaussian-heuristic model of the eff wall")
print("-" * 78)
print("""
  With S_K1 = n/K1, S_K2 = n/K2, S_KAN = n, S_D = s:
     det(V1)  = n^(3m+1) * s / (K1*K2)^m ,   dim = 2m+2
     E||pv||^2= n^2 * (2m/3 + s^2/3 + 1)
  Setting E||pv|| < GH = sqrt(dim/(2*pi*e)) * det^(1/dim) and letting m -> inf
  (so det^(1/dim) -> n * sqrt(n/(K1*K2)) = n/sqrt(eff), dim/(2*pi*e) -> 2m/(2*pi*e)):

     n*sqrt(2m/3) < sqrt(2m/(2*pi*e)) * n / sqrt(eff)
     <=>  eff  <  3/(2*pi*e)  =  %.4f
""" % (3.0 / (2 * math.pi * math.e)))
print(f"  {'curve':<14} | {'m':>3} | {'K1':>3} | {'eff':>6} | {'E|pv|/GH':>9} | "
      f"{'v1(S_D=1)':>10} | {'predict':>8}")
print("  " + "-" * 74)
EFF_CRIT = 3.0 / (2 * math.pi * math.e)
for label, cur, _k1, m in curves[1:]:
    n = cur[2]; k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        eff = k1 * k2b / n
        epv = expected_planted_norm(m, n, k1, k2b, 'v1', 1)
        gh = gaussian_heuristic(m, n, k1, k2b, 'v1', 1)
        w, tot, _, _ = rate(cur, m, k1, SEEDS, 'v1', s_d=1)
        pred = "OK" if eff < EFF_CRIT else "FAIL"
        print(f"  {label:<14} | {m:>3} | {k1:>3} | {eff:>6.3f} | {epv/gh:>9.3f} | "
              f"{w:>4}/{tot:<4}  | {pred:>8}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)


# ===========================================================================
# EXP R5 -- if the trivial vector is not lambda_1 any more, WHAT is?
#
# R1/R2 show sv/pv stays ~0.34-0.79 in every variant, so some OTHER vector is
# still shorter than the planted one.  Two candidates:
#   (a) the 2-D "lambda block"  < (n*S_K1, 0), (-lam*S_K1, S_K2) >  living in
#       the (col j, col m+1+j) plane -- present in V1 AND V2 alike;
#   (b) something spread across the C1 columns.
# Identify by dumping the per-column-group energy of the shortest reduced row.
# ===========================================================================

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0: break
        q = (2*num + den)//(2*den) if num >= 0 else -((-2*num + den)//(2*den))
        v = (v[0] - q*u[0], v[1] - q*u[1])
        if nrm2(v) >= nrm2(u): break
        u, v = v, u
    return u

print("\n" + "-" * 78)
print("EXP R5: identity of the surviving lambda_1 (energy profile of sv)")
print("-" * 78)
print(f"  {'curve':<12} | {'variant':<10} | {'|sv|/n':>7} | {'mu_blk/n':>8} | "
      f"{'k1cols':>7} | {'dcol':>6} | {'k2cols':>7} | {'kan':>6} | {'supp':>5}")
print("  " + "-" * 86)
for label, cur, k1, m in curves:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    S_K1, S_K2 = n // k1, max(1, n // k2b)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    mu = math.sqrt(w[0]**2 + w[1]**2)
    sopt = int(math.ceil(math.sqrt(m + 1.5)))
    for variant, s_d in [('v1', 1), ('v1', sopt), ('v2', 1)]:
        d_trial = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_trial, m, n, lam, p, k1, k2b, 42)
        if variant == 'v1':
            M = build_v1(sigs, n, lam, k1, k2b, s_d); dim = 2*m + 2
            pv = planted_v1(sigs, d_trial, n, k1, k2b, s_d)
            red = reduce_rows(M, dim)
            grp = {'k1cols': range(0, m), 'dcol': range(m, m+1),
                   'k2cols': range(m+1, 2*m+1), 'kan': range(2*m+1, 2*m+2)}
        else:
            M, _ = build_v2(sigs, n, lam, k1, k2b); dim = 2*m + 1
            pv = planted_v2(sigs, n, k1, k2b)
            red = reduce_rows(M, dim)
            grp = {'k1cols': list(range(0, m-1)) + [m-1], 'dcol': [],
                   'k2cols': [m] + list(range(m+1, 2*m)), 'kan': [2*m]}
        sv = min((r for r in red if any(r)), key=norm)
        tot = sum(x*x for x in sv) or 1
        e = {k: math.sqrt(sum(sv[j]**2 for j in v) / tot) for k, v in grp.items()}
        supp = sum(1 for x in sv if x)
        name = f"v1(S_D={s_d})" if variant == 'v1' else 'v2'
        print(f"  {label:<12} | {name:<10} | {norm(sv)/n:>7.3f} | {mu/n:>8.3f} | "
              f"{e['k1cols']:>7.3f} | {e['dcol']:>6.3f} | {e['k2cols']:>7.3f} | "
              f"{e['kan']:>6.3f} | {supp:>5}")

print("\n  (energies are fractions of ||sv||; they square-sum to 1."
      "  supp = #nonzero coords)")

# --- BKZ at the wall: is the wall a reduction-strength artefact? ------------
print("\n" + "-" * 78)
print("EXP R5b: BKZ(30) at the K1 wall -- reduction strength or information?")
print("-" * 78)
print(f"  {'curve':<12} | {'K1':>3} | {'variant':<10} | {'LLL':>6} | {'BKZ30':>7}")
print("  " + "-" * 50)
for label, cur, _k1, m in curves[1:]:
    n = cur[2]
    sopt = int(math.ceil(math.sqrt(m + 1.5)))
    walls = [8, 12, 16] if '2557' in label else [6, 8, 12]
    for k1 in walls:
        for variant, s_d in [('v1', 1), ('v2', 1)]:
            wl, tot, _, _ = rate(cur, m, k1, SEEDS, variant, s_d=s_d)
            wb, _, _, _ = rate(cur, m, k1, SEEDS, variant, s_d=s_d,
                               use_bkz=True, beta=30)
            name = f"v1(S_D={s_d})" if variant == 'v1' else 'v2'
            print(f"  {label:<12} | {k1:>3} | {name:<10} | {wl:>3}/{tot:<2} | "
                  f"{wb:>3}/{tot:<3}")
