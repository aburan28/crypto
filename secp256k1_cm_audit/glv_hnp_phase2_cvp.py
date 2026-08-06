"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so recovery is a CVP,
not an SVP -- and centre the target.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, experiment T5):
  The Phase-2 Kannan-embedded lattice ALWAYS contains the vector

        n * S_D * e_m          ( = n*row_m - sum_i B_i*row_i )

  of norm n*S_D, while the planted vector has expected norm
  n*S_D*sqrt(2m/3 + 4/3) > n*S_D for every m >= 1.  So the planted vector
  is NEVER lambda_1; measured sv/pv was in [0.34, 0.61] on every curve,
  success and failure alike.  Kannan embedding reduces CVP to SVP only
  when the embedded target IS the shortest vector -- here it provably is
  not.  Thread 23 (proposed 2026-07-29): drop the Kannan row/column and
  solve the underlying CVP directly.

Construction.  Let L' = the Phase-2 lattice with the Kannan row and column
removed (dim 2m+1, coordinates [k1-block | d | k2-block]):

    rows 0..m-1     :  n*S_K1 * e_i                        (modular)
    row  m          :  (B_i*S_K1)_i  at cols 0..m-1, S_D at col m
    rows m+1..2m    :  -lam*S_K1 at col i, S_K2 at col m+1+i

and let the target be  t = (A_i * S_K1)_i  at cols 0..m-1, 0 elsewhere.
Writing v in L' as sum_i c_i*row_i + x*row_m + sum_i y_i*row_{m+1+i}:

    (t - v)[i]     = (A_i - c_i*n - x*B_i + y_i*lam) * S_K1
    (t - v)[m]     = -x * S_D
    (t - v)[m+1+i] = -y_i * S_K2

Substituting  A_i + B_i*d = k1_i + lam*k2_i (mod n)  shows x = -d,
y_i = -k2_i, and

    t - v  =  (k1_i*S_K1, d*S_D, k2_i*S_K2)                  [*]

so d is read off coordinate m of the residual.  [*] is verified exactly in
EXP C0 below.

CENTRING.  [*] is the observation that makes the second half of this
thread necessary.  The true residual lives in the *positive orthant* box
[0, K1*S_K1) x [0, n*S_D) x [0, K2*S_K2), i.e. its coordinates are
uniform on [0, W) with W ~ n, not on [-W/2, W/2).  Its expected squared
norm is therefore (2m+1)*W^2/3 rather than (2m+1)*W^2/12: the uncentred
target is 2x longer than it needs to be, uniformly in m.  Both the
2026-07-26 Kannan lattice and the naive CVP above pay this factor.  The
fix is to translate the target by the box centre

    c = ( (K1-1)*S_K1/2, (n-1)*S_D/2, (K2-1)*S_K2/2 )

so that the residual becomes uniform on [-W/2, W/2) per coordinate, and
add c[m] back when reading off d.  This is free -- it changes no lattice,
only the coset representative.

The 2x2 design below therefore separates two independent effects:

              | uncentred      centred
    ----------+-------------------------------
    Kannan/SVP| KAN-LLL        KAN-LLL-C
    CVP on L' | BAB-LLL, CVP   BAB-LLL-C, CVP-C

where BAB-* is Babai nearest-plane on the LLL-reduced basis (polynomial
time) and CVP* is fplll's exact enumeration CVP (the ceiling of the
formulation, not a practical attack at scale).

Reference bound.  Expected number of spurious solutions in the box is
K1^m * K2^m * n / n^m, so a unique solution needs

    K1 * K2  <  n^((m-1)/m)                                  [IT]

Run: python3 glv_hnp_phase2_cvp.py
"""

import math
import random
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic  -- verbatim from glv_hnp_phase2_lambda_threshold.py:87
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
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m_ - i - 1), p)
        m_, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(4000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Phase-2 scales and signatures -- verbatim from
# glv_hnp_phase2_lambda_threshold.py:227 (so the comparison to the
# 2026-07-26 / 2026-07-29 numbers is exact)
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

# ---------------------------------------------------------------------------
# The de-embedded lattice L', its target, its centre
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """(2m+1)-dim lattice: the Kannan row AND column removed."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    return M

def cvp_target(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, _, _ = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m + 1)
    for i in range(m):
        t[i] = sigs[i]['A'] * S_K1
    return t

def box_centre(m, n, k1_bound, k2_bound):
    """Integer centre of the box the true residual [*] lives in."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    c = [0] * (2 * m + 1)
    c1 = ((k1_bound - 1) * S_K1) // 2
    c2 = ((k2_bound - 1) * S_K2) // 2
    for i in range(m):
        c[i] = c1
        c[m + 1 + i] = c2
    c[m] = ((n - 1) * S_D) // 2
    return c

def true_residual(sigs, d_secret, n, k1_bound, k2_bound):
    """The residual t - v induced by the true solution:  [*] in the docstring."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    e = [0] * (2 * m + 1)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1
        e[m + 1 + i] = sigs[i]['k2'] * S_K2
    e[m] = d_secret * S_D
    return e

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Kannan-embedded lattice (2026-07-26 baseline), uncentred and centred
# ---------------------------------------------------------------------------

def build_kannan_lattice(sigs, n, lam, k1_bound, k2_bound, centred=False):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    c = box_centre(m, n, k1_bound, k2_bound) if centred else [0] * (2 * m + 1)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    # last row = (t - c , S_KANNAN); reduces to (e - c, S_KANNAN) at the target
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1 - c[i]
        M[2 * m + 1][m + 1 + i] = -c[m + 1 + i]
    M[2 * m + 1][m] = -c[m]
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def recover_d_kannan(M_reduced, m, n, k1_bound, k2_bound, d_secret, centred):
    dim = 2 * m + 2
    _, S_D, _, S_KANNAN = scales(n, k1_bound, k2_bound)
    cd = box_centre(m, n, k1_bound, k2_bound)[m] if centred else 0
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m] + cd) % n          # S_D == 1
        if d_cand == d_secret % n:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Exact Babai nearest-plane (rational GS -- no float-precision question at
# these sizes, and reproducible across fpylll builds)
# ---------------------------------------------------------------------------

def _dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def babai_nearest_plane(B, t):
    """B: integer rows (LLL/BKZ reduced).  Returns v in L nearest-plane-close
    to t."""
    k = len(B)
    gs, gs_n2 = [], []
    for i in range(k):
        bi = [Fraction(x) for x in B[i]]
        v = list(bi)
        for j in range(i):
            if gs_n2[j] == 0:
                continue
            c = _dot(bi, gs[j]) / gs_n2[j]
            v = [v[q] - c * gs[j][q] for q in range(len(v))]
        gs.append(v)
        gs_n2.append(_dot(v, v))
    b = [Fraction(x) for x in t]
    for i in reversed(range(k)):
        if gs_n2[i] == 0:
            continue
        c = _dot(b, gs[i]) / gs_n2[i]
        num, den = c.numerator, c.denominator
        ci = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        b = [b[q] - ci * Fraction(B[i][q]) for q in range(len(b))]
    return [t[q] - int(b[q]) for q in range(len(t))]

# ---------------------------------------------------------------------------
# Methods -- each returns (recovered_ok, residual_norm)
# ---------------------------------------------------------------------------

def _reduce(M, beta=None):
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if beta is None:
        LLL.reduction(A)
    else:
        BKZ.reduction(A, BKZ.Param(beta))
    return A, [[A[i][j] for j in range(dim)] for i in range(dim)]

def run_kannan(sigs, curve, d_secret, k1, k2b, centred=False, beta=None):
    p, b, n, lam, G = curve
    m = len(sigs)
    M = build_kannan_lattice(sigs, n, lam, k1, k2b, centred=centred)
    _, red = _reduce(M, beta)
    ok = recover_d_kannan(red, m, n, k1, k2b, d_secret, centred) is not None
    return ok, min(norm(r) for r in red)

def _cvp_common(sigs, curve, k1, k2b, centred):
    p, b, n, lam, G = curve
    m = len(sigs)
    M = build_cvp_lattice(sigs, n, lam, k1, k2b)
    t = cvp_target(sigs, n, k1, k2b)
    cd = 0
    if centred:
        c = box_centre(m, n, k1, k2b)
        t = [t[i] - c[i] for i in range(len(t))]
        cd = c[m]
    return M, t, cd, m, n

def run_babai(sigs, curve, d_secret, k1, k2b, centred=False, beta=None):
    M, t, cd, m, n = _cvp_common(sigs, curve, k1, k2b, centred)
    _, red = _reduce(M, beta)
    v = babai_nearest_plane(red, t)
    e = [t[i] - v[i] for i in range(len(t))]
    d_cand = (e[m] + cd) % n                      # S_D == 1
    return (d_cand == d_secret % n), norm(e)

def run_cvp_exact(sigs, curve, d_secret, k1, k2b, centred=False):
    M, t, cd, m, n = _cvp_common(sigs, curve, k1, k2b, centred)
    A, _ = _reduce(M, None)
    try:
        v = CVP.closest_vector(A, tuple(t))
    except Exception:                             # pragma: no cover
        return None, float('nan')
    e = [t[i] - v[i] for i in range(len(t))]
    d_cand = (e[m] + cd) % n
    return (d_cand == d_secret % n), norm(e)

# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

METHODS = ("KAN-LLL", "KAN-LLL-C", "BAB-LLL", "BAB-LLL-C", "CVP", "CVP-C")

_DISPATCH = {
    "KAN-LLL":   lambda s, c, d, a, b: run_kannan(s, c, d, a, b, centred=False),
    "KAN-LLL-C": lambda s, c, d, a, b: run_kannan(s, c, d, a, b, centred=True),
    "BAB-LLL":   lambda s, c, d, a, b: run_babai(s, c, d, a, b, centred=False),
    "BAB-LLL-C": lambda s, c, d, a, b: run_babai(s, c, d, a, b, centred=True),
    "CVP":       lambda s, c, d, a, b: run_cvp_exact(s, c, d, a, b, centred=False),
    "CVP-C":     lambda s, c, d, a, b: run_cvp_exact(s, c, d, a, b, centred=True),
    "BAB-BKZ-C": lambda s, c, d, a, b: run_babai(s, c, d, a, b, centred=True, beta=20),
}

def trial(curve, m, k1, seed, methods=METHODS):
    """One instance; every method sees the SAME signature set and secret."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2b, seed)
    if len(sigs) < m:
        return None
    out = {k: _DISPATCH[k](sigs, curve, d_secret, k1, k2b) for k in methods}
    e = true_residual(sigs, d_secret, n, k1, k2b)
    c = box_centre(m, n, k1, k2b)
    out["_true"] = norm(e)
    out["_true_c"] = norm([e[i] - c[i] for i in range(len(e))])
    return out

def rates(curve, m, k1, seeds, methods=METHODS):
    acc = {k: 0 for k in methods}
    tot = 0
    gap = []          # ||t-v_CVP-C|| / ||t-v_true_centred||
    for s in seeds:
        r = trial(curve, m, k1, s, methods)
        if r is None:
            continue
        tot += 1
        for k in methods:
            if r[k][0]:
                acc[k] += 1
        if "CVP-C" in methods and not math.isnan(r["CVP-C"][1]) and r["_true_c"]:
            gap.append(r["CVP-C"][1] / r["_true_c"])
    return acc, tot, (sum(gap) / len(gap) if gap else float('nan'))

def it_bound_k1(n, m, k2_bound):
    """Largest K1 with K1*K2 < n^((m-1)/m)  -- the count bound [IT]."""
    return n ** ((m - 1) / m) / k2_bound

def fmt(acc, tot, methods):
    return " ".join(f"{str(acc[k]) + '/' + str(tot):>9}" for k in methods)

# ===========================================================================
# Experiments (guarded so the helpers above can be imported by
# glv_hnp_phase2_cvp_sweep.py without re-running them)
# ===========================================================================

def main():
    SEEDS = [42, 1234, 9999, 555, 31337]

    # Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26)
    HIST = [
        # label,             p,    b, n,    lam,  K1_hist, m_hist
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

    print("=" * 92)
    print("Thread 23 -- Phase-2 as a CVP instead of a Kannan/SVP instance, and centring the target")
    print("=" * 92)

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 92)
    print("EXP C0: sanity -- does the residual identity [*] hold exactly?")
    print("-" * 92)
    ok_all = True
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sg = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
        t = cvp_target(sg, n, k1, k2b)
        e = true_residual(sg, d0, n, k1, k2b)
        v = [t[i] - e[i] for i in range(len(t))]
        S_K1, S_D, S_K2, _ = scales(n, k1, k2b)
        y = [v[m + 1 + i] // S_K2 for i in range(m)]
        x = v[m] // S_D
        good = all(v[m + 1 + i] % S_K2 == 0 for i in range(m)) and v[m] % S_D == 0
        for i in range(m):
            resid = v[i] - (x * sg[i]['B'] * S_K1 - y[i] * lam * S_K1)
            good = good and (resid % (n * S_K1) == 0)
        good = good and (x % n == (-d0) % n) and all(
            y[i] % n == (-sg[i]['k2']) % n for i in range(m))
        ok_all = ok_all and good
        print(f"  {label:<18} m={m:<3} t-v is a lattice coset elt: {str(good):<6} "
              f"x == -d: {x % n == (-d0) % n}")
    print(f"  ==> [*] verified on all historical curves: {ok_all}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 92)
    print("EXP C0b: the centring factor -- ||e|| vs ||e - c||")
    print("-" * 92)
    print("Predicted ratio ||e-c||/||e|| -> sqrt((1/12)/(1/3)) = 0.500 for large m.\n")
    print(f"  {'curve':<18} {'m':>3} {'K1':>4} {'||e||':>11} {'||e-c||':>11} {'ratio':>7}")
    for label, curve, _k1, _m in hist:
        p, b, n, lam, G = curve
        for m_try in (6, 12, 24):
            r = trial(curve, m_try, 8 if n > 500 else 2, 42, methods=("CVP-C",))
            if r is None:
                continue
            print(f"  {label:<18} {m_try:>3} {8 if n > 500 else 2:>4} "
                  f"{r['_true']:>11.1f} {r['_true_c']:>11.1f} "
                  f"{r['_true_c'] / r['_true']:>7.4f}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 92)
    print("EXP C1: the 2x2 design on the 2026-07-29 T4 grid (m=12, 5 seeds)")
    print("-" * 92)
    print("T4 baseline (KAN-LLL) walls were: lam*=0.340 -> K1 ~ 12-16,")
    print("                                   lam*=0.070 -> K1 ~ 4-6.")
    print("Falsifier (2026-07-29): does the reformulation push the wall outward?\n")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    M_FIX = 12

    for label, curve, _k1, _m in hist:
        p, b, n, lam, G = curve
        if n < 2000:
            continue
        k2b = math.isqrt(n) + 1
        print(f"\ncurve {label}   n={n}  lam*={lam_star(lam, n):.4f}  K2={k2b}  m={M_FIX}")
        print(f"  [IT] count bound: K1*K2 < n^(11/12) = {n ** (11 / 12):.0f}"
              f"  ->  K1_max = {it_bound_k1(n, M_FIX, k2b):.1f}")
        print(f"  {'K1':>4} {'eff':>6} " + " ".join(f"{k:>9}" for k in METHODS)
              + f" {'cvpC/true':>10}")
        for k1 in K1_GRID:
            acc, tot, g = rates(curve, M_FIX, k1, SEEDS)
            print(f"  {k1:>4} {k1 * k2b / n:>6.3f} {fmt(acc, tot, METHODS)} {g:>10.4f}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 92)
    print("EXP C2: does MORE data now help?  (T4b re-run at K1=8, lam*=0.07 curve)")
    print("-" * 92)
    print("T4b (KAN-LLL): m = 8/12/16/24/32 -> 0,0,1,0,1 of 5, from which 2026-07-29")
    print("concluded the K1 wall is 'genuine'.  If a centred method rises monotonically")
    print("in m, that conclusion was an artefact of the embedding, not an info limit.\n")

    fail_curve = [c for lbl, c, _, _ in hist if c[2] == 2647][0]
    k2b_f = math.isqrt(fail_curve[2]) + 1
    print(f"  {'m':>4} {'IT K1max':>9} " + " ".join(f"{k:>9}" for k in METHODS))
    for m_try in (4, 6, 8, 12, 16, 24, 32):
        acc, tot, _ = rates(fail_curve, m_try, 8, SEEDS)
        print(f"  {m_try:>4} {it_bound_k1(fail_curve[2], m_try, k2b_f):>9.1f} "
              f"{fmt(acc, tot, METHODS)}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 92)
    print("EXP C3: where is the ceiling?  K1 sweep at large m on the lam*=0.07 curve")
    print("-" * 92)
    print("Practical (poly-time) method BAB-LLL-C vs the formulation ceiling CVP-C")
    print("vs the counting bound [IT].  A gap CVP-C < [IT] means the *lattice*")
    print("(not the information) is the binding constraint.\n")

    SUB = ("KAN-LLL", "BAB-LLL-C", "BAB-BKZ-C", "CVP-C")
    for m_try in (16, 24):
        print(f"\n  m = {m_try}   [IT] K1_max = {it_bound_k1(2647, m_try, k2b_f):.1f}")
        print(f"  {'K1':>4} {'eff':>6} " + " ".join(f"{k:>9}" for k in SUB))
        for k1 in (8, 12, 16, 24, 32, 48, 64):
            acc, tot, _ = rates(fail_curve, m_try, k1, SEEDS, methods=SUB)
            print(f"  {k1:>4} {k1 * k2b_f / 2647:>6.3f} {fmt(acc, tot, SUB)}")

    print("\n" + "=" * 92)
    print("done")
    print("=" * 92)


if __name__ == "__main__":
    main()
