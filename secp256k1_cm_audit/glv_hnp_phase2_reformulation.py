"""
GLV-HNP Phase 2 / Thread 23 — reformulate the lattice so the planted vector
can actually BE the shortest vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5):
  In the Phase-2 lattice of `glv_hnp_phase2_20bit.py:build_glv_lattice`, the
  shortest vector after LLL is ALWAYS the trivial vector n*S_D*e_m (100% of its
  energy in the d-column, |sv[m]|/n = 1.0000 exactly).  It carries no
  information (d is only defined mod n) and no choice of S_D removes it,
  because both it and the planted vector scale linearly in S_D.
  Measured sv/pv sits in [0.337, 0.368] for successes and failures alike.
  The 2026-07-29 entry proposed Thread 23: remove that direction and see
  whether the K1 wall on the lam*=0.070 curve (currently K1 ~ 4-6) moves.

This script runs a 2x2 factorial over two independent modifications:

  ELIM (d-elimination)
    Eliminate the secret d algebraically using signature 0 as a reference
    before building the lattice.  With k_i = A_i + B_i*d and t_i = B_i/B_0:
        k_i - t_i*k_0 = A_i - t_i*A_0 =: c_i   (mod n)
    so with the GLV split k_i = k1_i + lam*k2_i,
        k1_i = c_i + t_i*(k1_0 + lam*k2_0) - lam*k2_i   (mod n),  i = 1..m-1.
    The d-column is gone, hence so is the trivial vector n*S_D*e_m.
    Dimension drops 2m+2 -> 2m+1 and det drops by a factor n.
    d is recovered afterwards as d = (k1_0 + lam*k2_0 - A_0)/B_0 mod n.

  CENTER (recentring)
    k1 ~ U[0,K1) has E[k1^2] = K1^2/3; the shifted k1' = k1 - K1//2 has
    E[k1'^2] ~ K1^2/12.  Shifting k1, k2 (and d, when it is still present)
    to be centred on 0 only changes the Kannan row, and shrinks the planted
    vector by ~2x.  Pure BDD win, independent of ELIM.

Four formulations:
  R1  = baseline (verbatim 2026-06-15 construction)      dim 2m+2
  R1c = baseline + CENTER                                 dim 2m+2
  R2  = ELIM                                              dim 2m+1
  R2c = ELIM + CENTER                                     dim 2m+1

Falsifier (stated in the 2026-07-29 log entry):
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_reformulation.py
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py so the comparison is exact.
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
# CM helpers (j=0 curve search)
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
    """The two roots of x^2+x+1 = 0 mod n, sorted."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return tuple(sorted((r1, r2)))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def build_curve(p, n, seed=12345):
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
# Signatures
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- verbatim from the 2026-06-15 build."""
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
# R1 / R1c -- baseline lattice, optionally recentred
# ---------------------------------------------------------------------------
#
# Column layout (dim 2m+2):
#   col i        (i < m)      : residual column, value = k1_i * S_K1
#   col m                     : d column,        value = d     * S_D
#   col m+1+i    (i < m)      : k2 column,       value = k2_i  * S_K2
#   col 2m+1                  : Kannan column,   value = S_KANNAN
#
# With CENTER, the represented unknowns are k1' = k1 - h1, k2' = k2 - h2,
# d' = d - hd; only the Kannan row's entries change.

def build_R1(sigs, n, lam, k1_bound, k2_bound, center=False):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    h1 = k1_bound // 2 if center else 0
    h2 = k2_bound // 2 if center else 0
    hd = n // 2 if center else 0

    M = [[0] * dim for _ in range(dim)]
    for i in range(m):                       # n*z_i rows
        M[i][i] = n * S_K1
    for i in range(m):                       # d row
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):                       # k2_i rows
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):                       # Kannan row
        # k1'_i = (A_i - h1 + B_i*hd - lam*h2) + B_i*d' - lam*k2'_i - n*z_i
        const = (sigs[i]['A'] - h1 + sigs[i]['B'] * hd - lam * h2) % n
        M[2 * m + 1][i] = const * S_K1
    M[2 * m + 1][dim - 1] = S_KAN
    return M, (h1, h2, hd)

def planted_R1(sigs, d_secret, n, k1_bound, k2_bound, center=False):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    h1 = k1_bound // 2 if center else 0
    h2 = k2_bound // 2 if center else 0
    hd = n // 2 if center else 0
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - h1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - h2) * S_K2
    v[m] = (d_secret - hd) * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_R1(reduced, m, n, S_KAN, shifts, sigs, lam, d_secret):
    _, _, hd = shifts
    for row in reduced:
        last = row[2 * m + 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m] + hd) % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# R2 / R2c -- d eliminated via signature 0
# ---------------------------------------------------------------------------
#
# Congruence, for i = 1..m-1:
#   k1_i = c_i + t_i*k1_0 + t_i*lam*k2_0 - lam*k2_i - n*z_i
# with t_i = B_i/B_0 mod n and c_i = A_i - t_i*A_0 mod n.
#
# Column layout (dim 2m+1):
#   col i-1      (i = 1..m-1) : residual column, value = k1_i * S_K1
#   col m-1                   : k1_0 column,     value = k1_0 * S_K1
#   col m                     : k2_0 column,     value = k2_0 * S_K2
#   col m+i      (i = 1..m-1) : k2_i column,     value = k2_i * S_K2
#   col 2m                    : Kannan column,   value = S_KANNAN

def elim_coeffs(sigs, n):
    """(t_i, c_i) for i = 1..m-1, using signature 0 as the reference."""
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    B0inv = modinv(B0, n)
    out = []
    for i in range(1, len(sigs)):
        t = sigs[i]['B'] * B0inv % n
        c = (sigs[i]['A'] - t * A0) % n
        out.append((t, c))
    return out

def build_R2(sigs, n, lam, k1_bound, k2_bound, center=False):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    h1 = k1_bound // 2 if center else 0
    h2 = k2_bound // 2 if center else 0
    tc = elim_coeffs(sigs, n)

    M = [[0] * dim for _ in range(dim)]
    for j in range(m - 1):                   # n*z_i rows
        M[j][j] = n * S_K1
    for j, (t, c) in enumerate(tc):          # k1_0 row
        M[m - 1][j] = t % n * S_K1
    M[m - 1][m - 1] = S_K1
    for j, (t, c) in enumerate(tc):          # k2_0 row
        M[m][j] = t * lam % n * S_K1
    M[m][m] = S_K2
    for j, (t, c) in enumerate(tc):          # k2_i rows
        M[m + 1 + j][j] = (-lam) % n * S_K1
        M[m + 1 + j][m + 1 + j] = S_K2
    for j, (t, c) in enumerate(tc):          # Kannan row
        # k1'_i = [c_i + (t_i - 1)*h1 + (t_i - 1)*lam*h2] + t_i*k1'_0
        #                            + t_i*lam*k2'_0 - lam*k2'_i - n*z_i
        const = (c + (t - 1) * h1 + (t - 1) * lam * h2) % n
        M[2 * m][j] = const * S_K1
    M[2 * m][dim - 1] = S_KAN
    return M, (h1, h2, 0)

def planted_R2(sigs, n, k1_bound, k2_bound, center=False):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    h1 = k1_bound // 2 if center else 0
    h2 = k2_bound // 2 if center else 0
    v = [0] * (2 * m + 1)
    for j in range(m - 1):
        v[j] = (sigs[j + 1]['k1'] - h1) * S_K1
        v[m + 1 + j] = (sigs[j + 1]['k2'] - h2) * S_K2
    v[m - 1] = (sigs[0]['k1'] - h1) * S_K1
    v[m] = (sigs[0]['k2'] - h2) * S_K2
    v[2 * m] = S_KAN
    return v

def recover_R2(reduced, m, n, lam, S_K1, S_K2, S_KAN, shifts, sigs, d_secret):
    """d = (k1_0 + lam*k2_0 - A_0) / B_0 mod n."""
    h1, h2, _ = shifts
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    B0inv = modinv(B0, n)
    for row in reduced:
        last = row[2 * m]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        if (sign * row[m - 1]) % S_K1 or (sign * row[m]) % S_K2:
            continue
        k1_0 = sign * row[m - 1] // S_K1 + h1
        k2_0 = sign * row[m] // S_K2 + h2
        k0 = (k1_0 + lam * k2_0) % n
        d_cand = (k0 - A0) * B0inv % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Lattice-geometry diagnostics
# ---------------------------------------------------------------------------

def log_det(M):
    """log(|det|) via fraction-free Gaussian elimination on a copy (exact)."""
    from sympy import Matrix
    return float(sympy.log(abs(Matrix(M).det())))

def gaussian_heuristic(logdet, dim):
    """GH(L) = det^(1/dim) * sqrt(dim / (2*pi*e))."""
    return math.exp(logdet / dim) * math.sqrt(dim / (2 * math.pi * math.e))

# ---------------------------------------------------------------------------
# One experiment
# ---------------------------------------------------------------------------

FORMS = ('R1', 'R1c', 'R2', 'R2c')

def run_one(curve, m, d_secret, k1_bound, form, seed=42,
            use_bkz=False, bkz_beta=20, want_geom=False):
    """Returns dict with ok / pv / sv / dim / (optionally) gh."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    center = form.endswith('c')

    if form.startswith('R1'):
        M, shifts = build_R1(sigs, n, lam, k1_bound, k2_bound, center)
        pv = planted_R1(sigs, d_secret, n, k1_bound, k2_bound, center)
    else:
        M, shifts = build_R2(sigs, n, lam, k1_bound, k2_bound, center)
        pv = planted_R2(sigs, n, k1_bound, k2_bound, center)
    dim = len(M)

    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]

    if form.startswith('R1'):
        d_out = recover_R1(reduced, m, n, S_KAN, shifts, sigs, lam, d_secret)
    else:
        d_out = recover_R2(reduced, m, n, lam, S_K1, S_K2, S_KAN,
                           shifts, sigs, d_secret)

    res = {'ok': d_out is not None, 'pv': norm(pv),
           'sv': min(norm(r) for r in reduced), 'dim': dim}
    if want_geom:
        ld = log_det(M)
        res['gh'] = gaussian_heuristic(ld, dim)
        res['logdet'] = ld
    return res

def rate(curve, m, k1_bound, form, seeds, use_bkz=False, bkz_beta=20):
    """(wins, total, mean sv/pv)."""
    p, b, n, lam, G = curve
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = run_one(curve, m, d_trial, k1_bound, form, seed,
                    use_bkz=use_bkz, bkz_beta=bkz_beta)
        if r is None:
            continue
        wins += bool(r['ok'])
        ratios.append(r['sv'] / r['pv'] if r['pv'] else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================

print("=" * 78)
print("Thread 23 — Phase-2 lattice reformulation (d-elimination x recentring)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,             p,    b, n,    lam,  K1_hist, m
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
# EXP U0 — correctness: the planted vector is really in each lattice, and the
#          recovery routine really recovers d from it.
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U0: sanity — planted vector in lattice + recovery round-trip")
print("-" * 78)
print("Solves M^T x = v_planted over Q and checks x is integral (=> v in L),")
print("then feeds the planted vector alone through the recovery routine.\n")
print(f"{'curve':<18} {'form':>5} {'dim':>4} {'v in L':>7} {'d round-trip':>13}")

from sympy import Matrix

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sg = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    for form in FORMS:
        center = form.endswith('c')
        if form.startswith('R1'):
            M, shifts = build_R1(sg, n, lam, k1, k2b, center)
            v = planted_R1(sg, d0, n, k1, k2b, center)
            got = recover_R1([v], m, n, S_KAN, shifts, sg, lam, d0)
        else:
            M, shifts = build_R2(sg, n, lam, k1, k2b, center)
            v = planted_R2(sg, n, k1, k2b, center)
            got = recover_R2([v], m, n, lam, S_K1, S_K2, S_KAN, shifts, sg, d0)
        x = Matrix(M).T.solve(Matrix(v))
        integral = all(e.q == 1 for e in x)
        print(f"{label:<18} {form:>5} {len(M):>4} {str(integral):>7} "
              f"{str(got == d0):>13}")
    print()


# ---------------------------------------------------------------------------
# EXP U1 — lattice geometry: pv / GH for each formulation
# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP U1: planted norm vs Gaussian heuristic")
print("-" * 78)
print("pv/GH < 1 means the planted vector is plausibly the shortest vector.")
print("sv/pv < 1 means something ELSE is shorter (2026-07-29 EXP T5).\n")
print(f"{'curve':<18} {'form':>5} {'dim':>4} {'pv':>13} {'GH':>13} "
      f"{'pv/GH':>7} {'sv/pv':>7}")

u1 = {}
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    for form in FORMS:
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        r = run_one(curve, m, d0, k1, form, 42, want_geom=True)
        u1[(label, form)] = r
        print(f"{label:<18} {form:>5} {r['dim']:>4} {r['pv']:>13.1f} "
              f"{r['gh']:>13.1f} {r['pv']/r['gh']:>7.3f} "
              f"{r['sv']/r['pv']:>7.3f}")
    print()


# ---------------------------------------------------------------------------
# EXP U2 — THE FALSIFIER: K1 wall on the lam*=0.070 curve, all 4 formulations
# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP U2: K1 wall (the 2026-07-29 falsifier)")
print("-" * 78)
print("2026-07-29 EXP T4 measured, for R1: 2557 wall at K1~12-16,")
print("2677 wall at K1~4-6.  Does any reformulation move the wall outward?\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]

for label, curve, k1_hist, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"{label}  (lam*={lam_star(lam,n):.4f}, n={n}, m={m}, K2={k2b})")
    hdr = "  ".join(f"{k:>5}" for k in K1_GRID)
    print(f"{'form':>5}  {hdr}")
    for form in FORMS:
        cells = []
        for k1 in K1_GRID:
            w, t, _ = rate(curve, m, k1, form, SEEDS)
            cells.append(f"{w}/{t}")
        print(f"{form:>5}  " + "  ".join(f"{c:>5}" for c in cells))
    # sv/pv along the same grid, to see whether pv ever becomes shortest
    print(f"{'sv/pv':>5}")
    for form in FORMS:
        cells = []
        for k1 in K1_GRID:
            _, _, sp = rate(curve, m, k1, form, SEEDS)
            cells.append(f"{sp:.2f}")
        print(f"{form:>5}  " + "  ".join(f"{c:>5}" for c in cells))
    print()


# ---------------------------------------------------------------------------
# EXP U3 — does BKZ change the R2c picture at the wall?
# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP U3: BKZ(beta=40) at and beyond the wall, R1 vs R2c")
print("-" * 78)
print("If the wall is information-theoretic, stronger reduction cannot help.\n")
print(f"{'curve':<18} {'K1':>4} {'R1 LLL':>7} {'R1 BKZ':>7} "
      f"{'R2c LLL':>8} {'R2c BKZ':>8}")
for label, curve, k1_hist, m in hist:
    for k1 in (8, 12, 16, 24):
        row = []
        for form, bkz in (('R1', False), ('R1', True),
                          ('R2c', False), ('R2c', True)):
            w, t, _ = rate(curve, m, k1, form, SEEDS,
                           use_bkz=bkz, bkz_beta=40)
            row.append(f"{w}/{t}")
        print(f"{label:<18} {k1:>4} {row[0]:>7} {row[1]:>7} "
              f"{row[2]:>8} {row[3]:>8}")
    print()


# ---------------------------------------------------------------------------
# EXP U4 — information-theoretic reference line
# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP U4: where SHOULD the wall be?")
print("-" * 78)
print("Unknown bits = m*(log2 K1 + log2 K2); available bits = (m-1)*log2 n")
print("(one signature's worth is consumed by d).  K1_info = the largest K1")
print("with unknown bits < available bits.\n")
print(f"{'curve':<18} {'n':>7} {'m':>3} {'K2':>5} {'K1_info':>8} "
      f"{'K1_obs(R1)':>11} {'K1_obs(R2c)':>12}")

def k1_info_bound(n, m, k2):
    """Largest K1 with m*(lg K1 + lg K2) < (m-1)*lg n."""
    budget = (m - 1) * math.log2(n) / m - math.log2(k2)
    return 2 ** budget

def observed_wall(curve, m, form, seeds, grid):
    """Largest K1 in `grid` with a full sweep of successes."""
    best = None
    for k1 in grid:
        w, t, _ = rate(curve, m, k1, form, seeds)
        if w == t:
            best = k1
        else:
            break
    return best

for label, curve, k1_hist, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    ki = k1_info_bound(n, m, k2b)
    w1 = observed_wall(curve, m, 'R1', SEEDS, K1_GRID)
    w2 = observed_wall(curve, m, 'R2c', SEEDS, K1_GRID)
    print(f"{label:<18} {n:>7} {m:>3} {k2b:>5} {ki:>8.1f} "
          f"{str(w1):>11} {str(w2):>12}")


# ---------------------------------------------------------------------------
# EXP U5 — independent confirmation on fresh 17-bit curves
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP U5: fresh 17-bit curves, R1 vs R2c at the wall")
print("-" * 78)

M17 = 12
fresh = search_curves(2 ** 16, 2 ** 17, want=6)
print(f"Found {len(fresh)} fresh 17-bit j=0 GLV curves, m={M17}.")
print("K1 grid is chosen per-curve as fractions of the information bound")
print("K1_info, so the comparison sits ON the wall rather than far below it.\n")
print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K2':>5} {'K1inf':>6} {'K1':>5} "
      f"{'eff':>6} {'R1':>6} {'R1c':>6} {'R2':>6} {'R2c':>6}")

u5 = {f: [0, 0] for f in FORMS}
u5_wall = {f: [] for f in FORMS}
for (p, b, n, lam, G) in fresh:
    curve = (p, b, n, lam, G)
    k2b = math.isqrt(n) + 1
    ki = k1_info_bound(n, M17, k2b)
    grid = sorted({max(2, int(ki * f)) for f in (0.25, 0.5, 0.75, 1.0, 1.5)})
    for k1 in grid:
        cells = []
        for form in FORMS:
            w, t, _ = rate(curve, M17, k1, form, SEEDS[:3])
            u5[form][0] += w
            u5[form][1] += t
            cells.append(f"{w}/{t}")
        print(f"{p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k2b:>5} {ki:>6.1f} "
              f"{k1:>5} {k1*k2b/n:>6.3f} " + " ".join(f"{c:>6}" for c in cells))
    for form in FORMS:
        w = observed_wall(curve, M17, form, SEEDS[:3], grid)
        u5_wall[form].append((w / ki) if w else 0.0)
    print()

print("Totals across the 17-bit set:")
for form in FORMS:
    w, t = u5[form]
    print(f"  {form:>4}: {w}/{t}  ({100.0*w/t:.1f}%)" if t else f"  {form}: n/a")
print("\nMean observed wall as a fraction of the information bound K1_info:")
for form in FORMS:
    xs = u5_wall[form]
    print(f"  {form:>4}: {sum(xs)/len(xs):.3f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
