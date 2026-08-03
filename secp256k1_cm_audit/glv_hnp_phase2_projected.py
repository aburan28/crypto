"""
GLV-HNP -- Thread 23: make the planted vector the target, not a runner-up.

Background (2026-07-29, log T5)
-------------------------------
The Phase-2 lattice `build_glv_lattice()` in
`glv_hnp_phase2_lambda_threshold.py:254` is (2m+2)-dimensional with column
scales (S_K1, S_D, S_K2, S_KANNAN) = (n//K1, 1, n//K2, n).  Its d-column is
scaled by S_D = 1 while d ~ U[0,n) is *not* small, and the lattice therefore
contains the trivial vector

    n * S_D * e_m          (= n*row_m reduced by the B_i*row_i)

of norm n*S_D, whereas ||v_planted|| ~ n*sqrt(2m/3 + 4/3).  The trivial vector
is shorter for every m >= 1, it is the LLL shortest vector on 100% of curves
tested, and it carries no information (d is only defined mod n).  Rescaling
S_D does not help: both vectors scale linearly in S_D.

So Phase-2 recovery was never an SVP condition.  This script tests the
reformulation proposed on 2026-07-29: remove the d-direction entirely.

Variant P (this run's proposal) -- eliminate d algebraically
------------------------------------------------------------
With A_i + B_i*d = k_i (mod n) and k_i = k1_i + lam*k2_i, pick signature 0 as
reference and set

    C_i = B_i * B_0^{-1}  mod n,      T_i = A_i - C_i*A_0  mod n     (i=1..m-1)

Then d cancels:

    k_i - C_i*k_0 = T_i   (mod n)

Unknowns are the 2m small values (k1_0, k2_0, k1_i, k2_i), d is gone.  Fusing
the congruence column with the k1_i value column exactly as the original
construction does gives a (2m+1)-dimensional lattice with NO under-scaled
column, hence no trivial short vector.  d is read back afterwards from
d = (k_0 - A_0) * B_0^{-1} mod n.

Variant B -- Babai / CVP on the un-embedded baseline lattice
-------------------------------------------------------------
Drop the Kannan row from the baseline lattice, LLL-reduce, and solve CVP for
the target (-A_i*S_K1, 0, ..., 0) by exact-rational nearest-plane.  This tests
whether the Kannan embedding (rather than the d-column per se) is what costs
the attack.

Falsifier (stated 2026-07-29)
-----------------------------
  * if sv/pv rises above 1 AND the K1 wall on the lam*=0.070 curve
    (12-bit/2677, currently K1 ~ 4-6) moves outward -> the reformulation is a
    real improvement;
  * if the wall stays at K1 ~ 4-6 -> the wall is information-theoretic and
    Phase 2 is at its ceiling.

Run: python3 glv_hnp_phase2_projected.py
"""

import importlib.util
import math
import random
import sys
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, BKZ

_spec = importlib.util.spec_from_file_location(
    "_t20a", __file__.rsplit("/", 1)[0] + "/glv_hnp_phase2_lambda_threshold.py")
_t20a = importlib.util.module_from_spec(_spec)
# that module runs its own T1..T5 experiments at import time; mute them
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    _spec.loader.exec_module(_t20a)

modinv = _t20a.modinv
find_generator = _t20a.find_generator
gen_signatures = _t20a.gen_signatures
build_glv_lattice = _t20a.build_glv_lattice
planted_vector = _t20a.planted_vector
scales = _t20a.scales
norm = _t20a.norm
recover_d = _t20a.recover_d
lam_star = _t20a.lam_star
eisenstein_decompose = _t20a.eisenstein_decompose
j0_traces = _t20a.j0_traces
glv_roots = _t20a.glv_roots
build_curve = _t20a.build_curve

# ---------------------------------------------------------------------------
# Variant P: d-eliminated lattice, dimension 2m+1
# ---------------------------------------------------------------------------
#
# column layout (dim = 2m+1):
#   0 .. m-2      congruence column for signature i = j+1, fused with k1_i value
#   m-1           k1_0 value
#   m             k2_0 value
#   m+1 .. 2m-1   k2_i value for i = j+1
#   2m            Kannan column
#
# row layout (dim = 2m+1):
#   0 .. m-2      n*S_K1 * e_j                       (mod-n rows)
#   m-1 .. 2m-3   k2_i rows: -lam*S_K1 at col j, S_K2 at col m+1+j
#   2m-2          k1_0 row:  C_i*S_K1 at col j,     S_K1 at col m-1
#   2m-1          k2_0 row:  (C_i*lam mod n)*S_K1 at col j, S_K2 at col m
#   2m            Kannan row: T_i*S_K1 at col j,    S_KANNAN at col 2m

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    lam_r = lam % n

    B0inv = modinv(sigs[0]['B'], n)
    C = [0] * m
    T = [0] * m
    for i in range(1, m):
        C[i] = sigs[i]['B'] * B0inv % n
        T[i] = (sigs[i]['A'] - C[i] * sigs[0]['A']) % n

    M = [[0] * dim for _ in range(dim)]
    # mod-n rows
    for j in range(m - 1):
        M[j][j] = n * S_K1
    # k2_i rows
    for j in range(m - 1):
        r = m - 1 + j
        M[r][j] = -lam_r * S_K1
        M[r][m + 1 + j] = S_K2
    # k1_0 row
    r = 2 * m - 2
    for j in range(m - 1):
        M[r][j] = C[j + 1] * S_K1
    M[r][m - 1] = S_K1
    # k2_0 row
    r = 2 * m - 1
    for j in range(m - 1):
        M[r][j] = C[j + 1] * lam_r % n * S_K1
    M[r][m] = S_K2
    # Kannan row
    r = 2 * m
    for j in range(m - 1):
        M[r][j] = T[j + 1] * S_K1
    M[r][2 * m] = S_KAN
    return M, C, T


def projected_planted_vector(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * dim
    for j in range(m - 1):
        v[j] = sigs[j + 1]['k1'] * S_K1
        v[m + 1 + j] = sigs[j + 1]['k2'] * S_K2
    v[m - 1] = sigs[0]['k1'] * S_K1
    v[m] = sigs[0]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def recover_d_projected(M_reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Scan reduced rows for |last| == S_KANNAN, read (k1_0, k2_0), rebuild d."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        a, b = sign * row[m - 1], sign * row[m]
        if a % S_K1 or b % S_K2:
            continue
        k1_0, k2_0 = a // S_K1, b // S_K2
        k_0 = (k1_0 + lam * k2_0) % n
        d_cand = (k_0 - sigs[0]['A']) * B0inv % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None



# ---------------------------------------------------------------------------
# Variant PC: projected lattice, RECENTERED
# ---------------------------------------------------------------------------
#
# Both the baseline and Variant P are *uncentered*: k1_i in [0,K1) maps to the
# coordinate k1_i*S_K1 in [0,n), so the planted vector sits in the positive
# orthant with E[||pv||^2] = n^2*(2m/3 + 1).  Subtracting the box centre from
# the Kannan row (which enters the planted combination with coefficient
# exactly 1) moves every coordinate into [-n/2, n/2) and gives
# E[||pv||^2] = n^2*(m/6 + 1) -- shorter by ~sqrt(2m/3)/sqrt(m/6) -> 2 for
# large m.  This is the standard Boneh-Venkatesan / Nguyen-Shparlinski
# recentring, and it is absent from the Phase-2 construction.

def build_projected_lattice_centered(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M, C, T = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    c1, c2 = k1_bound // 2, k2_bound // 2
    r = 2 * m
    for j in range(m - 1):
        M[r][j] -= c1 * S_K1
        M[r][m + 1 + j] -= c2 * S_K2
    M[r][m - 1] -= c1 * S_K1
    M[r][m] -= c2 * S_K2
    return M, C, T


def projected_planted_vector_centered(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = projected_planted_vector(sigs, n, k1_bound, k2_bound)
    c1, c2 = k1_bound // 2, k2_bound // 2
    for j in range(m - 1):
        v[j] -= c1 * S_K1
        v[m + 1 + j] -= c2 * S_K2
    v[m - 1] -= c1 * S_K1
    v[m] -= c2 * S_K2
    return v


def recover_d_projected_centered(M_reduced, sigs, n, lam, k1_bound, k2_bound,
                                 d_secret):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1, c2 = k1_bound // 2, k2_bound // 2
    B0inv = modinv(sigs[0]['B'], n)
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        a, b = sign * row[m - 1], sign * row[m]
        if a % S_K1 or b % S_K2:
            continue
        k1_0, k2_0 = a // S_K1 + c1, b // S_K2 + c2
        k_0 = (k1_0 + lam * k2_0) % n
        d_cand = (k_0 - sigs[0]['A']) * B0inv % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None


# ---------------------------------------------------------------------------
# Variant B: exact-rational Babai nearest-plane CVP on the baseline lattice
# ---------------------------------------------------------------------------

def _gso_exact(basis):
    """Gram-Schmidt over Q.  Returns (B_star, mu) with exact Fractions."""
    bstar, mu = [], []
    for i, b in enumerate(basis):
        v = [Fraction(x) for x in b]
        row_mu = []
        for j in range(i):
            bj = bstar[j]
            den = sum(x * x for x in bj)
            c = Fraction(0) if den == 0 else sum(Fraction(b[k]) * bj[k]
                                                 for k in range(len(b))) / den
            row_mu.append(c)
            if c:
                v = [v[k] - c * bj[k] for k in range(len(v))]
        bstar.append(v)
        mu.append(row_mu)
    return bstar, mu


def babai_nearest_plane(basis, target):
    """Exact nearest-plane.  basis: list of int rows (LLL-reduced)."""
    bstar, _ = _gso_exact(basis)
    w = [Fraction(x) for x in target]
    coeffs = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        den = sum(x * x for x in bstar[i])
        if den == 0:
            continue
        c = sum(w[k] * bstar[i][k] for k in range(len(w))) / den
        q = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = q
        if q:
            w = [w[k] - q * basis[i][k] for k in range(len(w))]
    approx = [0] * len(target)
    for i, q in enumerate(coeffs):
        if q:
            for k in range(len(target)):
                approx[k] += q * basis[i][k]
    return approx, coeffs


def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Baseline lattice minus the Kannan row: dim 2m+1, columns 0..2m."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    lam_r = lam % n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam_r * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    return M


def run_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Target t = (-A_i*S_K1, 0, ..., 0); closest vector should be
    (k1_i*S_K1 - A_i*S_K1, d*S_D, k2_i*S_K2)."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    target = [0] * dim
    for i in range(m):
        target[i] = -sigs[i]['A'] * S_K1
    approx, _ = babai_nearest_plane(red, target)
    d_cand = (approx[m] // S_D) % n if S_D else None
    err = [approx[k] - target[k] for k in range(dim)]
    return (d_cand == d_secret), norm(err)


# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False,
          bkz_beta=20):
    """variant in {'base','proj','cvp'} -> (ok, planted_norm, shortest_norm)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None

    if variant == 'cvp':
        ok, errn = run_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret)
        pn = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
        return (ok, pn, errn)

    if variant == 'base':
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    elif variant == 'projc':
        M, _, _ = build_projected_lattice_centered(sigs, n, lam, k1_bound,
                                                   k2_bound)
        pv = projected_planted_vector_centered(sigs, n, k1_bound, k2_bound)
    else:
        M, _, _ = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = projected_planted_vector(sigs, n, k1_bound, k2_bound)
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]

    if variant == 'base':
        _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
        ok = recover_d(red, m, n, S_KAN, d_secret) is not None
    elif variant == 'projc':
        ok = recover_d_projected_centered(red, sigs, n, lam, k1_bound,
                                          k2_bound, d_secret) is not None
    else:
        ok = recover_d_projected(red, sigs, n, lam, k1_bound, k2_bound,
                                 d_secret) is not None
    sn = min(norm(r) for r in red)
    return (ok, norm(pv), sn)


def rate(curve, m, k1_bound, seeds, variant, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, variant,
                    use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        ok, pn, sn = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr


# ---------------------------------------------------------------------------

SEEDS = [1, 2, 3, 4, 5]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

# (label, p, b, n, lam, m) -- the historical Phase-2 curves.
# m = 12 on the two 12-bit curves so the 'base' rows reproduce EXP T4
# (2026-07-29) exactly; the 8-bit curve has K2 = 15 and is not grid-comparable.
HIST = [
    ("8-bit/199",   211,  2, 199,  106,  6),
    ("12-bit/2557", 2557, 2, 2659, 1755, 12),
    ("12-bit/2677", 2677, 2, 2647, 185,  12),
]


def load_hist():
    out = []
    for label, p, b, n, lam, m in HIST:
        G = find_generator(p, b, n)
        if G is None:
            print(f"  !! no generator for {label}", file=sys.stderr)
            continue
        out.append((label, (p, b, n, lam, G), m, lam_star(lam, n)))
    return out


def main():
    print("=" * 74)
    print("Thread 23 -- reformulated GLV-HNP Phase-2 lattice")
    print("=" * 74)

    curves = load_hist()

    # -- E1: sanity, does the projected lattice reproduce the planted vector? --
    print("\nE1  construction check (planted vector is in the lattice, "
          "d reads back)")
    print(f"  {'curve':14s} {'variant':6s} {'dim':>4s} {'||pv||/n':>9s} "
          f"{'||sv||/||pv||':>13s} {'rec':>4s}")
    for label, curve, m, mu in curves:
        p, b, n, lam, G = curve
        for variant in ('base', 'proj'):
            k1 = 4
            res = trial(curve, m, 12345 % n or 7, k1, 1, variant)
            if res is None:
                print(f"  {label:14s} {variant:6s}  signature generation failed")
                continue
            ok, pn, sn = res
            dim = 2 * m + 2 if variant == 'base' else 2 * m + 1
            print(f"  {label:14s} {variant:6s} {dim:4d} {pn / n:9.3f} "
                  f"{sn / pn:13.4f} {str(ok):>4s}")

    # -- E2: the K1 wall, baseline vs projected --------------------------------
    print("\nE2  K1 wall  (wins out of 5 seeds; K1*K2/n effective bias grows "
          "left to right)")
    for label, curve, m, mu in curves:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        print(f"\n  {label}  (n={n}, lam*={mu:.4f}, m={m}, K2={k2})")
        hdr = "  " + " ".join(f"{k:>5d}" for k in K1_GRID)
        print(f"    {'K1':10s}{hdr}")
        eff = "  " + " ".join(f"{k * k2 / n:5.2f}" for k in K1_GRID)
        print(f"    {'eff':10s}{eff}")
        for variant in ('base', 'proj', 'projc', 'cvp'):
            cells = []
            for k1 in K1_GRID:
                w, t, _ = rate(curve, m, k1, SEEDS, variant)
                cells.append(f"{w}/{t}")
            print(f"    {variant:10s}  " + " ".join(f"{c:>5s}" for c in cells))

    # -- E3: sv/pv across the grid (is the planted vector now lambda_1?) -------
    print("\nE3  mean ||sv||/||pv|| across the K1 grid  "
          "(>1 means planted is the shortest)")
    for label, curve, m, mu in curves:
        print(f"\n  {label}")
        for variant in ('base', 'proj', 'projc'):
            cells = []
            for k1 in K1_GRID:
                _, _, mr = rate(curve, m, k1, SEEDS, variant)
                cells.append(f"{mr:5.2f}")
            print(f"    {variant:10s}  " + " ".join(f"{c:>5s}" for c in cells))

    # -- E4: BKZ on the projected lattice at the wall -------------------------
    print("\nE4  BKZ(beta=30) vs LLL on the projected lattice at the wall")
    print(f"  {'curve':14s} {'K1':>4s} {'proj-LLL':>9s} {'proj-BKZ30':>11s}")
    for label, curve, m, mu in curves:
        for k1 in (8, 12, 16, 24):
            wl, t, _ = rate(curve, m, k1, SEEDS, 'proj')
            wb, _, _ = rate(curve, m, k1, SEEDS, 'proj',
                            use_bkz=True, bkz_beta=30)
            print(f"  {label:14s} {k1:4d} {f'{wl}/{t}':>9s} {f'{wb}/{t}':>11s}")

    # -- E5: the L2 parasite -- mu is m-independent, ||pv|| grows like sqrt(m) -
    print("\nE5  the lambda-block parasite L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)>")
    print("    mu = lambda_1(L2) does NOT depend on m; ||pv|| ~ n*sqrt(2m/3+1).")
    print("    So the planted vector cannot be lambda_1 for m above m_cross.")
    print(f"  {'curve':14s} {'K1':>4s} {'mu':>10s} {'copies':>7s} "
          f"{'||pv||(m=12)':>13s} {'m_cross':>8s}")
    for label, curve, m, mu_star in curves:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in (4, 8):
            S_K1, _, S_K2, S_KAN = scales(n, k1, k2)
            mu_l2, _ = _t20a.lambda_block_mu(n, lam, S_K1, S_K2)
            pv12 = _t20a.planted_norm_expected(12, n, k1, k2, S_K1, 1, S_K2,
                                               S_KAN)
            # ||pv(m)||^2 = m*((K1*S_K1)^2+(K2*S_K2)^2)/3 + n^2/3 + S_KAN^2
            per = ((k1 * S_K1) ** 2 + (k2 * S_K2) ** 2) / 3.0
            const = n ** 2 / 3.0 + S_KAN ** 2
            m_cross = (mu_l2 ** 2 - const) / per if per else float('nan')
            print(f"  {label:14s} {k1:4d} {mu_l2:10.1f} {m - 1:7d} "
                  f"{pv12:13.1f} {m_cross:8.2f}")

    # -- E6: does shrinking m below the parasite crossover rescue recovery? ---
    print("\nE6  m-sweep at the wall (K1=8, eff~0.157; info-theoretic floor "
          "m~4.3)")
    print("    T4b tested m = 8..32 only.  If mu is the obstruction, small m "
          "should win.")
    M_GRID = [3, 4, 5, 6, 7, 8, 10, 12]
    for label, curve, _m, mu_star in curves:
        p, b, n, lam, G = curve
        if n < 2000:
            continue
        print(f"\n  {label} (lam*={mu_star:.4f})")
        print("    m         " + " ".join(f"{v:>5d}" for v in M_GRID))
        for variant in ('base', 'proj'):
            cells = []
            for mm in M_GRID:
                w, t, _ = rate(curve, mm, 8, SEEDS, variant)
                cells.append(f"{w}/{t}")
            print(f"    {variant:10s}" + " ".join(f"{c:>5s}" for c in cells))

    run_e7()
    run_e8()
    run_e9()
    print("\nDone.")




# ---------------------------------------------------------------------------
# E7: is the wall information-theoretic?  Exhaustive uniqueness check.
# ---------------------------------------------------------------------------
#
# For a toy-size n we can settle the question directly instead of inferring it
# from reduction failures.  A candidate secret d is *consistent* with the
# signature set iff for every i the value k = A_i + B_i*d (mod n) admits a
# decomposition k = k1 + lam*k2 (mod n) with 0 <= k1 < K1 and 0 <= k2 < K2.
# If more than one d is consistent, no algorithm can recover the true one --
# the wall is information-theoretic, not a failure of lattice reduction.

def decomposable(k, n, lam, k1_bound, k2_bound):
    for k2 in range(k2_bound):
        if (k - lam * k2) % n < k1_bound:
            return True
    return False


def consistent_secrets(sigs, n, lam, k1_bound, k2_bound, cap=None):
    out = []
    for d in range(1, n):
        for s in sigs:
            k = (s['A'] + s['B'] * d) % n
            if not decomposable(k, n, lam, k1_bound, k2_bound):
                break
        else:
            out.append(d)
            if cap and len(out) >= cap:
                break
    return out


# ---------------------------------------------------------------------------
# Variant BC: baseline lattice (d-column KEPT), recentred.  Attribution control.
# ---------------------------------------------------------------------------
# Comparing BC against PC isolates the two independent changes:
#   proj  = d-column removed, uncentred
#   basec = d-column kept,    recentred
# If basec ~ projc then recentring is the whole effect and the trivial
# n*S_D*e_m vector diagnosed on 2026-07-29 is a red herring.

def build_base_lattice_centered(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    c1, c2 = k1_bound // 2, k2_bound // 2
    r = 2 * m + 1                      # Kannan row: coefficient 1 in planted
    for i in range(m):
        M[r][i] -= c1 * S_K1
        M[r][m + 1 + i] -= c2 * S_K2
    return M


def trial_basec(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    dim = 2 * m + 2
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_base_lattice_centered(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for i in range(dim):
        if abs(A[i][dim - 1]) != S_KAN:
            continue
        sg = 1 if A[i][dim - 1] > 0 else -1
        if (sg * A[i][m]) % n == d_secret:
            return True
    return False


def rate_basec(curve, m, k1_bound, seeds):
    n = curve[2]
    w = 0
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        if trial_basec(curve, m, d, k1_bound, s):
            w += 1
    return w


def run_e8():
    """Attribution: which of the two changes moves the wall?"""
    print("\nE8  attribution -- recentring vs d-elimination (m=12, 5 seeds)")
    for label, curve, m, mu_star in load_hist():
        if curve[2] < 2000:
            continue
        print(f"\n  {label} lam*={mu_star:.4f}")
        print("    K1        " + " ".join(f"{k:>5d}" for k in K1_GRID))
        for v in ('base', 'proj', 'projc'):
            cells = [f"{rate(curve, m, k, SEEDS, v)[0]}/5" for k in K1_GRID]
            print(f"    {v:10s}" + " ".join(f"{c:>5s}" for c in cells))
        cells = [f"{rate_basec(curve, m, k, SEEDS)}/5" for k in K1_GRID]
        print(f"    {'basec':10s}" + " ".join(f"{c:>5s}" for c in cells))


def run_e9():
    """Replication of the recentring gain on fresh 17-bit curves."""
    print("\nE9  replication on fresh 17-bit j=0 GLV curves (m=12, 5 seeds)")
    with contextlib.redirect_stdout(io.StringIO()):
        found = _t20a.search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=6)
    grid = (8, 16, 24, 32)
    print(f"  {'p':>7s} {'n':>7s} {'lam*':>6s} | "
          + " ".join(f"K1={k:<4d}" for k in grid))
    for item in list(found)[:6]:
        p, b, n, lam, G = item[1] if isinstance(item[0], str) else item
        row = []
        for k1 in grid:
            wb = rate((p, b, n, lam, G), 12, k1, SEEDS, 'base')[0]
            wc = rate((p, b, n, lam, G), 12, k1, SEEDS, 'projc')[0]
            row.append(f"{wb}->{wc}")
        print(f"  {p:7d} {n:7d} {lam_star(lam, n):6.3f} | "
              + " ".join(f"{c:<7s}" for c in row))
    print("\n  (cells are base -> projc, wins out of 5)")


def run_e7():
    print("\nE7  uniqueness of the secret (exhaustive over all d in [1,n))")
    print("    #consistent d = 1  -> wall is algorithmic (reduction too weak)")
    print("    #consistent d > 1  -> wall is information-theoretic")
    curves = load_hist()
    print(f"  {'curve':14s} {'K1':>4s} {'m':>3s} {'eff':>6s} {'seed':>5s} "
          f"{'#consistent d':>14s} {'true d in set':>14s} {'LLL':>5s}")
    for label, curve, m, mu_star in curves:
        p, b, n, lam, G = curve
        if n < 2000:
            continue
        k2b = math.isqrt(n) + 1
        for k1b in (4, 6, 8, 12):
            for seed in (1, 2, 3):
                d_true = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_true, m, n, lam, p, k1b, k2b, seed)
                if len(sigs) < m:
                    continue
                cons = consistent_secrets(sigs, n, lam, k1b, k2b)
                res = trial(curve, m, d_true, k1b, seed, 'proj')
                ok = res[0] if res else None
                print(f"  {label:14s} {k1b:4d} {m:3d} {k1b * k2b / n:6.2f} "
                      f"{seed:5d} {len(cons):14d} "
                      f"{str(d_true in cons):>14s} {str(ok):>5s}")

if __name__ == "__main__":
    main()
