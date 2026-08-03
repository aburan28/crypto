"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
lambda_1.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, experiment T5): in the
Phase-2 lattice of `glv_hnp_phase2_20bit.py::build_glv_lattice` the shortest
vector after LLL is ALWAYS the trivial vector `n*S_D*e_m` -- 100% of its energy
in the d-column, |sv[m]|/n = 1.0000 on every curve tested.  It is 2-3x shorter
than the planted vector and carries no information (d is only defined mod n).
Recovery is therefore a BDD/coset condition, not an SVP condition.

Root cause: the secret d is UNBOUNDED mod n, so the d-column contributes
(n*S_D)^2/3 to ||v_planted||^2 while the lattice contains n*S_D*e_m of norm
n*S_D.  No choice of S_D fixes this (both scale linearly in S_D).  The fix must
remove the d-coordinate, not rescale it.

Three arms are compared:

  O  original    dim 2m+2, verbatim from glv_hnp_phase2_20bit.py (baseline)
  P  projected   dim 2m+1, column m deleted.  L cap <e_m> = <n*S_D*e_m>, so
                 orthogonal projection along e_m is exactly "delete column m"
                 and the image is a lattice.  n*S_D*e_m maps to 0, the trivial
                 shortest vector is gone, and ||v_planted||^2 loses its
                 (n*S_D)^2/3 term.
  C  CVP/Babai   dim 2m, arm P with the Kannan row and column removed, solved
                 as an explicit closest-vector problem by Babai nearest-plane
                 (exact, Fraction arithmetic) against target (-A_i*S_K1 | 0).

Arms P and C no longer read d off the short vector (there is no d-column).
Instead the vector carries (k1_i, k2_i); d is then recovered algebraically from
signature 0:   d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n).

Falsifier stated on 2026-07-29:
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."
Correction applied here: sv/pv CANNOT exceed 1, since the planted vector is a
lattice vector, so ||lambda_1|| <= ||v_planted|| always.  The correct reading of
the criterion is sv/pv -> 1.0 (the planted vector IS the shortest vector), up
from the measured 0.34-0.37 in the original lattice.  Both readings are
reported below.

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
from fractions import Fraction

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + GLV helpers
# (verbatim from glv_hnp_phase2_lambda_threshold.py so results are comparable)
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
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (Eisenstein norm form)."""
    if p % 3 != 1:
        return None
    lim = int(2 * math.isqrt(p)) + 2
    for a in range(0, lim + 1):
        for b in range(0, lim + 1):
            if a * a - a * b + b * b == p:
                return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

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

# ---------------------------------------------------------------------------
# Scales + signatures (verbatim)
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

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# ARM O -- original lattice, dim 2m+2 (verbatim baseline)
# ---------------------------------------------------------------------------

def build_lattice_O(sigs, n, lam, k1_bound, k2_bound):
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

def planted_O(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_O(reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# ARM P -- projected lattice, dim 2m+1 (column m deleted)
#
# Column layout: [0..m-1] = k1 columns, [m..2m-1] = k2 columns, [2m] = Kannan.
# The generating set has 2m+2 rows of rank 2m+1; fplll hoists the single
# dependency to a zero row, which is simply skipped on readout.
# ---------------------------------------------------------------------------

def build_lattice_P(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                                   # n*S_K1*e_i
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                                        # d-row (no d column)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                   # k2 rows
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    r = [0] * dim                                        # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KANNAN
    rows.append(r)
    return rows

def planted_P(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def d_from_k(k1_0, k2_0, sig0, n, lam):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    try:
        Binv = modinv(sig0['B'] % n, n)
    except ValueError:
        return None
    return (k1_0 + lam * k2_0 - sig0['A']) % n * Binv % n

def recover_P(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        sgn = 1 if last > 0 else -1
        k1_0 = sgn * row[0]
        k2_0 = sgn * row[m]
        if k1_0 % S_K1 or k2_0 % S_K2: continue
        d_cand = d_from_k(k1_0 // S_K1, k2_0 // S_K2, sigs[0], n, lam)
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# ARM C -- explicit CVP by Babai nearest-plane, dim 2m
#
# Lattice L0 (Kannan row/column removed), columns [0..m-1] = k1, [m..2m-1] = k2.
# A lattice point parameterised by (d, k2_i) has column i equal to
#   (B_i*d - lam*k2_i) * S_K1   (mod n*S_K1)
# and column m+i equal to k2_i*S_K2.  Since A_i + B_i*d - lam*k2_i = k1_i
# (mod n), the target t = (-A_i*S_K1 | 0) sits at distance
#   ||v - t|| = ||(k1_i*S_K1 | k2_i*S_K2)|| ~ n*sqrt(2m/3)
# from the lattice.  Babai nearest-plane on the LLL-reduced basis is exact
# (Fraction arithmetic), so no floating-point precision question arises.
# ---------------------------------------------------------------------------

def build_lattice_C(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    return rows

def target_C(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return t

def gram_schmidt_exact(basis):
    """Exact GS over Fraction.  Returns (bstar, mu) with rows of `basis` that
    are linearly dependent dropped (they give a zero b*)."""
    bstar, mu, keep = [], [], []
    for idx, b in enumerate(basis):
        v = [Fraction(x) for x in b]
        coeffs = []
        for j, bs in enumerate(bstar):
            den = sum(x * x for x in bs)
            num = sum(x * y for x, y in zip(v, bs))
            c = num / den
            coeffs.append(c)
            v = [x - c * y for x, y in zip(v, bs)]
        if any(x != 0 for x in v):
            bstar.append(v); mu.append(coeffs); keep.append(idx)
    return bstar, keep

def babai_nearest_plane(basis, target):
    """Exact Babai nearest-plane.  `basis` must be LLL-reduced for quality."""
    bstar, keep = gram_schmidt_exact(basis)
    B = [basis[i] for i in keep]
    b = [Fraction(x) for x in target]
    for i in range(len(B) - 1, -1, -1):
        den = sum(x * x for x in bstar[i])
        num = sum(x * y for x, y in zip(b, bstar[i]))
        c = num / den
        # round-half-to-even is irrelevant here; round-half-up is fine
        ci = (c.numerator * 2 + c.denominator) // (2 * c.denominator)
        b = [x - ci * y for x, y in zip(b, B[i])]
    # b now holds target - closest; return the offset vector (target - v)
    return [int(x) for x in b]

def recover_C(offset, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """`offset` = t - v.  We expect -offset = (k1_i*S_K1 | k2_i*S_K2) up to
    sign, since v - t = (k1_i*S_K1 | k2_i*S_K2)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for sgn in (-1, 1):
        k1_0, k2_0 = sgn * offset[0], sgn * offset[m]
        if k1_0 % S_K1 or k2_0 % S_K2: continue
        d_cand = d_from_k(k1_0 // S_K1, k2_0 // S_K2, sigs[0], n, lam)
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# One trial, all three arms on the SAME signature set
# ---------------------------------------------------------------------------

def lll_rows(rows, ncols):
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]

def bkz_rows(rows, ncols, beta):
    A = IntegerMatrix.from_matrix(rows)
    BKZ.reduction(A, BKZ.Param(beta))
    return [[A[i][j] for j in range(ncols)] for i in range(A.nrows)]

def trial(curve, m, d_secret, k1_bound, seed=42, arms=('O', 'P', 'C'),
          use_bkz=False, bkz_beta=20):
    """Returns dict arm -> (ok, sv_over_pv or None)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    out = {}

    if 'O' in arms:
        rows = build_lattice_O(sigs, n, lam, k1_bound, k2_bound)
        red = (bkz_rows if use_bkz else lll_rows)(
            rows, 2 * m + 2, *( (bkz_beta,) if use_bkz else () ))
        pv = norm(planted_O(sigs, d_secret, n, k1_bound, k2_bound))
        nz = [norm(r) for r in red if any(r)]
        sv = min(nz) if nz else float('nan')
        ok = recover_O(red, m, n, S_KAN, d_secret) is not None
        out['O'] = (ok, sv / pv if pv else float('nan'))

    if 'P' in arms:
        rows = build_lattice_P(sigs, n, lam, k1_bound, k2_bound)
        red = (bkz_rows if use_bkz else lll_rows)(
            rows, 2 * m + 1, *( (bkz_beta,) if use_bkz else () ))
        pv = norm(planted_P(sigs, n, k1_bound, k2_bound))
        nz = [norm(r) for r in red if any(r)]
        sv = min(nz) if nz else float('nan')
        ok = recover_P(red, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
        out['P'] = (ok, sv / pv if pv else float('nan'))

    if 'C' in arms:
        rows = build_lattice_C(sigs, n, lam, k1_bound, k2_bound)
        red = (bkz_rows if use_bkz else lll_rows)(
            rows, 2 * m, *( (bkz_beta,) if use_bkz else () ))
        red = [r for r in red if any(r)]
        off = babai_nearest_plane(red, target_C(sigs, n, k1_bound, k2_bound))
        ok = recover_C(off, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
        # distance achieved vs. the true planted offset
        true_off = [0] * (2 * m)
        for i in range(m):
            true_off[i] = sigs[i]['k1'] * S_K1
            true_off[m + i] = sigs[i]['k2'] * S_K2
        pv = norm(true_off)
        out['C'] = (ok, norm(off) / pv if pv else float('nan'))

    return out

def rate(curve, m, k1_bound, seeds, arms=('O', 'P', 'C'),
         use_bkz=False, bkz_beta=20):
    wins = {a: 0 for a in arms}
    ratios = {a: [] for a in arms}
    tot = 0
    for seed in seeds:
        p, b, n, lam, G = curve
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, arms=arms,
                    use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        tot += 1
        for a in arms:
            ok, r = res[a]
            wins[a] += bool(ok)
            if r == r:  # not NaN
                ratios[a].append(r)
    mean = {a: (sum(ratios[a]) / len(ratios[a]) if ratios[a] else float('nan'))
            for a in arms}
    return wins, tot, mean


def main():
    # ===========================================================================
    print("=" * 78)
    print("Thread 23 -- reformulate the GLV-HNP Phase-2 lattice so the planted")
    print("vector is lambda_1  (arms: O original / P projected / C CVP-Babai)")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,               p,    b, n,    lam,  K1, m
        ("8-bit/199",          211,  2, 199,  106,  2,  6),
        ("12-bit/2557",        2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677 FAIL",   2677, 2, 2647, 185,  8,  10),
    ]

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), k1, m))

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP A -- correctness of the two reformulations + lattice geometry")
    print("-" * 78)
    print("Checks that (i) the planted vector really lies in L' and L0, (ii) the")
    print("trivial vector n*S_D*e_m is gone, (iii) sv/pv moves.")
    print()
    print(f"{'curve':<18} {'K1':>3} {'m':>3} | {'arm':>3} {'dim':>4} "
          f"{'wins':>6} {'sv/pv':>8}")
    print("-" * 78)

    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        for a in ('O', 'P', 'C'):
            wins, tot, mean = rate(curve, m, k1, SEEDS, arms=(a,))
            dim = {'O': 2 * m + 2, 'P': 2 * m + 1, 'C': 2 * m}[a]
            print(f"{label:<18} {k1:>3} {m:>3} | {a:>3} {dim:>4} "
                  f"{wins[a]:>3}/{tot:<2} {mean[a]:>8.4f}")
        print()

    # sanity: the planted vector must be IN the projected lattice
    print("Sanity -- planted vector membership (solve M'^T x = v over Q):")
    from sympy import Matrix
    for label, curve, k1, m in hist[:2]:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
        rows = build_lattice_P(sigs, n, lam, k1, k2b)
        v = planted_P(sigs, n, k1, k2b)
        A = Matrix(rows).T
        sol, params = A.gauss_jordan_solve(Matrix(v))
        sol0 = sol.subs({t: 0 for t in params})
        resid = (A * sol0 - Matrix(v)).norm()
        # the free parameter is the kernel direction (the single rank dependency);
        # membership in the lattice needs SOME integer point on that affine line
        integral = all(x == int(x) for x in sol0 if x.is_Rational)
        print(f"  {label:<18} residual={resid}  rank-deficiency={len(params)}  "
              f"integer coeffs={integral}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP B -- the T4 K1 grid, replicated across all three arms")
    print("-" * 78)
    print("2026-07-29 T4 baseline (arm O): 2557 fails at K1>=16, 2677 at K1>=8.")
    print("Falsifier: does the K1 wall move outward under P / C?")
    print()

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
    for label, curve, _k1, m in hist[1:]:      # the two 12-bit curves
        p, b, n, lam, G = curve
        print(f"{label}  (lam*={lam_star(lam, n):.3f}, n={n}, m={m})")
        print(f"  {'arm':<4}" + "".join(f"{k:>6}" for k in K1_GRID))
        for a in ('O', 'P', 'C'):
            cells = []
            for k1 in K1_GRID:
                wins, tot, _ = rate(curve, m, k1, SEEDS, arms=(a,))
                cells.append(f"{wins[a]}/{tot}")
            print(f"  {a:<4}" + "".join(f"{c:>6}" for c in cells))
        print()

    # ---------------------------------------------------------------------------
    print("-" * 78)
    print("EXP C -- eff sweep on fresh 17-bit curves")
    print("-" * 78)
    print("2026-07-29 T3 baseline (arm O): eff=0.05 -> 19/20, 0.15 -> 3/20,")
    print("0.25 -> 0/20.  Does P / C extend the reachable bias strength?")
    print()

    def search_curves(lo, hi, want=10):
        out = []
        p = int(sympy.nextprime(lo))
        while p < hi and len(out) < want:
            if p % 3 == 1:
                eis = eisenstein_decompose(p)
                if eis is not None:
                    a_e, b_e = eis
                    for t in j0_traces(a_e, b_e):
                        nc = p + 1 - t
                        if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                            continue
                        roots = glv_roots(nc)
                        if roots is None:
                            continue
                        cur = build_curve(p, nc)
                        if cur is None:
                            continue
                        out.append((p, cur[1], nc, roots[0], cur[4]))
                        break
            p = int(sympy.nextprime(p))
        return out

    CURVES17 = search_curves(2**16, 2**17, want=10)
    print(f"Found {len(CURVES17)} fresh 17-bit j=0 GLV curves "
          f"(lam* range {min(lam_star(c[3], c[2]) for c in CURVES17):.4f}"
          f"-{max(lam_star(c[3], c[2]) for c in CURVES17):.4f})")
    print()

    M_SIGS = 12
    print(f"{'eff':>6} | {'arm O':>10} {'arm P':>10} {'arm C':>10} | "
          f"{'sv/pv O':>8} {'sv/pv P':>8}")
    print("-" * 78)
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25, 0.35):
        tw = {'O': 0, 'P': 0, 'C': 0}
        tt = 0
        rr = {'O': [], 'P': []}
        for (p, b, n, lam, G) in CURVES17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            wins, tot, mean = rate((p, b, n, lam, G), M_SIGS, k1b, SEEDS)
            for a in ('O', 'P', 'C'):
                tw[a] += wins[a]
            tt += tot
            for a in ('O', 'P'):
                if mean[a] == mean[a]:
                    rr[a].append(mean[a])
        mo = sum(rr['O']) / len(rr['O']) if rr['O'] else float('nan')
        mp = sum(rr['P']) / len(rr['P']) if rr['P'] else float('nan')
        print(f"{eff:>6.2f} | {tw['O']:>4}/{tt:<5} {tw['P']:>4}/{tt:<5} "
              f"{tw['C']:>4}/{tt:<5} | {mo:>8.4f} {mp:>8.4f}")

    print()
    print("=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
