"""
GLV-HNP Phase 2, Thread 23: reformulate the Phase-2 lattice so that the
planted vector is not dominated by trivial short vectors.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 / commit e845207, exp T5):
  The shortest vector of the Phase-2 lattice was found to be, on every curve
  tested, the trivial vector  n*S_D*e_m  (all its energy in the d column).
  That entry concluded:
     "no choice of S_D removes it -- both vectors scale linearly in S_D"
  and proposed Thread 23: quotient out the trivial direction, or replace the
  Kannan embedding by an explicit CVP call, then check whether the K1 wall
  measured in T4 moves outward.

  The "both scale linearly" statement is checked directly here (E1) and is
  found to be FALSE: only the d-coordinate of v_planted carries S_D, so
      ||v_planted||^2 = (d*S_D)^2 + R^2,   ||trivial||^2 = (n*S_D)^2
  with R^2 = sum(k1_i*S_K1)^2 + sum(k2_i*S_K2)^2 + S_KANNAN^2  independent of
  S_D.  Since d < n, the ratio ||v_planted||/||trivial|| -> d/n < 1 as
  S_D -> infinity.  So a large S_D *does* demote the trivial vector.  The
  experiments below measure whether that buys any recovery.

Experiments
  E1  S_D sweep: norm ratio ||pv||/||trivial||, rank of the planted vector by
      norm inside the LLL basis, and recovery rate, for S_D = 2^0..2^24.
  E2  d-column-dropped Kannan lattice (dim 2m+1 instead of 2m+2).  The d
      column is deleted entirely and d is read off the LLL transformation
      matrix U instead (the coefficient of the d-row).  This removes the
      trivial vector from the lattice by construction.
  E3  Explicit CVP: drop the d column AND the Kannan row (dim 2m), Babai
      nearest-plane against target t = (-A_i*S_K1, 0, ..., 0), d recovered
      from the Babai coefficients composed with U.
  E4  K1 wall for baseline / E2 / E3 side by side on the historical
      success/failure pair, against the information-theoretic bound
      m_min = log(n) / log(n/(K1*K2)).

Run: python3 glv_hnp_phase2_bdd.py
"""

import math
import random
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fpylll import IntegerMatrix, LLL, BKZ

import glv_hnp_lib as L
from glv_hnp_lib import (
    find_generator, gen_signatures, scales, build_glv_lattice,
    planted_vector, norm, recover_d, glv_roots, lam_star,
)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (identical table to the 2026-07-29 run).
#   label,             p,    b, n,    lam,  K1, m
HIST = [
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]


def load_hist():
    out = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, (p, b, n, lam, G), k1, m))
    return out


# ===========================================================================
# Variant A -- baseline lattice with a settable S_D
# ===========================================================================

def build_lattice_SD(sigs, n, lam, k1_bound, k2_bound, S_D):
    """Baseline Phase-2 lattice with S_D overridden.  S_D=1 reproduces
    build_glv_lattice() exactly (asserted in E1)."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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


def planted_SD(sigs, d_secret, n, k1_bound, k2_bound, S_D):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v


def run_A(curve, m, d_secret, k1_bound, seed, S_D):
    """Baseline with S_D override.  Returns dict of diagnostics or None."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_lattice_SD(sigs, n, lam, k1_bound, k2_bound, S_D)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

    # d is read off the d column; with S_D != 1 divide it back out.
    ok = False
    for row in red:
        if abs(row[dim - 1]) != S_KAN:
            continue
        sgn = 1 if row[dim - 1] > 0 else -1
        val = sgn * row[m]
        if val % S_D != 0:
            continue
        if (val // S_D) % n == d_secret:
            ok = True
            break

    pv = planted_SD(sigs, d_secret, n, k1_bound, k2_bound, S_D)
    pn = norm(pv)
    triv = n * S_D                      # ||n*S_D*e_m||
    norms = sorted(norm(r) for r in red)
    sn = norms[0]
    rank = sum(1 for x in norms if x < pn) + 1   # rank of pv among basis norms
    return {'ok': ok, 'pn': pn, 'triv': triv, 'sn': sn, 'rank': rank, 'dim': dim}


# ===========================================================================
# Variant B (E2) -- d column DROPPED, d read from the LLL transform U
# ===========================================================================

def build_lattice_nod(sigs, n, lam, k1_bound, k2_bound):
    """Same generators as the baseline but with column m (the d column)
    deleted.  Columns: 0..m-1 = k1 block, m..2m-1 = k2 block, 2m = Kannan.
    2m+2 generators spanning a rank-(2m+1) lattice: the relation
        n*row_m - sum_i B_i*row_i = 0
    is exactly the trivial vector, which is annihilated by the deletion."""
    m = len(sigs)
    ncols = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = [[0] * ncols for _ in range(2 * m + 2)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][2 * m] = S_KAN
    return M


def run_B(curve, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_lattice_nod(sigs, n, lam, k1_bound, k2_bound)
    nrows, ncols = 2 * m + 2, 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    U = IntegerMatrix.identity(nrows)
    if use_bkz:
        LLL.reduction(A, U)
        BKZ.reduction(A, BKZ.Param(bkz_beta))   # BKZ without transform
        # BKZ invalidates U, so only the LLL transform is usable; re-run LLL
        # semantics: fall back to the LLL result for coefficient recovery.
        A = IntegerMatrix.from_matrix(M)
        U = IntegerMatrix.identity(nrows)
        LLL.reduction(A, U)
    else:
        LLL.reduction(A, U)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

    ok = False
    for j in range(nrows):
        last = A[j][ncols - 1]
        if abs(last) != S_KAN:
            continue
        sgn = 1 if last > 0 else -1
        d_cand = (sgn * U[j][m]) % n
        if d_cand == d_secret:
            ok = True
            break

    # planted vector in this lattice = baseline planted with the d entry removed
    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    pv_nod = pv[:m] + pv[m + 1:]
    pn = norm(pv_nod)
    norms = sorted(norm([A[i][j] for j in range(ncols)]) for i in range(nrows))
    nz = [x for x in norms if x > 0]
    sn = nz[0] if nz else 0.0
    rank = sum(1 for x in nz if x < pn) + 1
    return {'ok': ok, 'pn': pn, 'sn': sn, 'rank': rank, 'dim': ncols}


# ===========================================================================
# Variant C (E3) -- explicit CVP (Babai nearest plane), no d column, no Kannan
# ===========================================================================

def build_lattice_cvp(sigs, n, lam, k1_bound, k2_bound):
    """Generators rows 0..2m (no Kannan row) over columns 0..2m-1
    (k1 block then k2 block, no d column).  Rank 2m."""
    m = len(sigs)
    ncols = 2 * m
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * ncols for _ in range(2 * m + 1)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    return M


def babai_nearest_plane(basis, target):
    """Babai nearest plane on an (already LLL-reduced) list-of-lists basis.
    Returns the integer coefficient vector c with sum c_i*basis_i ~ target."""
    k = len(basis)
    if k == 0:
        return []
    # float Gram-Schmidt
    gs, mu = [], [[0.0] * k for _ in range(k)]
    for i in range(k):
        v = [float(x) for x in basis[i]]
        for j in range(i):
            d = sum(basis[i][t] * gs[j][t] for t in range(len(v)))
            nj = sum(x * x for x in gs[j])
            mu[i][j] = d / nj if nj else 0.0
            for t in range(len(v)):
                v[t] -= mu[i][j] * gs[j][t]
        gs.append(v)
    b = [float(x) for x in target]
    c = [0] * k
    for i in range(k - 1, -1, -1):
        nj = sum(x * x for x in gs[i])
        if nj == 0:
            continue
        ci = round(sum(b[t] * gs[i][t] for t in range(len(b))) / nj)
        c[i] = int(ci)
        for t in range(len(b)):
            b[t] -= c[i] * basis[i][t]
    return c


def run_C(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_lattice_cvp(sigs, n, lam, k1_bound, k2_bound)
    nrows, ncols = 2 * m + 1, 2 * m
    A = IntegerMatrix.from_matrix(M)
    U = IntegerMatrix.identity(nrows)
    LLL.reduction(A, U)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)

    rows = [[A[i][j] for j in range(ncols)] for i in range(nrows)]
    keep = [i for i in range(nrows) if any(rows[i])]
    basis = [rows[i] for i in keep]

    # Target: the lattice point closest to t should be at distance ||planted||.
    # u - t = (k1_i*S_K1, k2_i*S_K2) when u = d*row_m + sum k2_i*row_{m+1+i} + red.
    t = [0] * ncols
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1

    c = babai_nearest_plane(basis, t)
    # coefficient of the original d-row (generator index m)
    coef_m = sum(c[a] * U[keep[a]][m] for a in range(len(keep)))
    d_cand = coef_m % n

    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    pv_cvp = pv[:m] + pv[m + 1:2 * m + 1]
    pn = norm(pv_cvp)
    # actual Babai residual
    u = [sum(c[a] * basis[a][j] for a in range(len(basis))) for j in range(ncols)]
    resid = norm([u[j] - t[j] for j in range(ncols)])
    return {'ok': d_cand == d_secret, 'pn': pn, 'resid': resid, 'dim': ncols}


# ===========================================================================
# Success-rate drivers
# ===========================================================================

def rate(fn, curve, m, k1_bound, seeds, **kw):
    p, b, n, lam, G = curve
    wins, extra = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = fn(curve, m, d_trial, k1_bound, seed, **kw)
        if r is None:
            continue
        wins += bool(r['ok'])
        extra.append(r)
    return wins, len(seeds), extra


def fmt(w, t):
    return f"{w}/{t}"


# ===========================================================================
print("=" * 78)
print("Thread 23 -- Phase-2 lattice reformulation (BDD / d-column / CVP)")
print("=" * 78)

hist = load_hist()

print("\nCurves (identical to the 2026-07-29 table):")
print(f"{'curve':<18} {'p':>6} {'n':>6} {'bits':>5} {'lam':>6} {'lam*':>7} {'K2':>5}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    print(f"{label:<18} {p:>6} {n:>6} {n.bit_length():>5} {lam:>6} "
          f"{lam_star(lam, n):>7.4f} {math.isqrt(n)+1:>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1: does S_D demote the trivial vector n*S_D*e_m ?")
print("-" * 78)
print("2026-07-29 claimed 'both vectors scale linearly in S_D'.  Only the")
print("d-coordinate of v_planted does; the k1/k2/Kannan coordinates do not.")
print("Prediction: ||pv||/||triv|| -> d/n < 1 as S_D grows, so the trivial")
print("vector stops being shortest.  Question: does recovery follow?\n")

# sanity: S_D=1 must reproduce the archived construction exactly
_lbl, _cur, _k1, _m = hist[1]
_p, _b, _n, _lam, _G = _cur
_k2b = math.isqrt(_n) + 1
_sg = gen_signatures(_G, 12345, _m, _n, _lam, _p, _k1, _k2b, 42)
assert build_lattice_SD(_sg, _n, _lam, _k1, _k2b, 1) == \
       build_glv_lattice(_sg, _n, _lam, _k1, _k2b), "S_D=1 must match archive"
print("[check] build_lattice_SD(...,S_D=1) == build_glv_lattice(...)  OK\n")

SD_GRID = [1, 4, 16, 64, 256, 1 << 10, 1 << 14, 1 << 18, 1 << 24]
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    print(f"{label}   (K1={k1}, m={m})")
    print(f"  {'S_D':>10} {'||pv||/||triv||':>16} {'||sv||/||pv||':>14} "
          f"{'rank(pv)':>9} {'LLL':>6}")
    for sd in SD_GRID:
        w, t, ex = rate(run_A, curve, m, k1, SEEDS, S_D=sd)
        if not ex:
            continue
        r_pt = sum(e['pn'] / e['triv'] for e in ex) / len(ex)
        r_sp = sum(e['sn'] / e['pn'] for e in ex) / len(ex)
        rk = sum(e['rank'] for e in ex) / len(ex)
        print(f"  {sd:>10} {r_pt:>16.4f} {r_sp:>14.4f} {rk:>9.1f} {fmt(w,t):>6}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E2/E3: d-column-dropped Kannan lattice, and explicit Babai CVP")
print("-" * 78)
print("E2 deletes the d column (dim 2m+2 -> 2m+1); d is read from the LLL")
print("transformation matrix U instead.  The trivial vector is not in the")
print("lattice at all.  E3 additionally deletes the Kannan row and solves")
print("CVP directly by Babai nearest plane (dim 2m).\n")

print(f"{'curve':<18} {'K1':>3} {'m':>3} {'base':>6} {'E2':>6} {'E3':>6} "
      f"{'dim b/E2/E3':>13}")
for label, curve, k1, m in hist:
    wb, tb, _ = rate(lambda c, mm, d, k, s: (lambda r: r and
                     {'ok': r[0], 'pn': r[1], 'sn': r[2]})(
                        L.run_experiment(c, mm, d, k, s)),
                     curve, m, k1, SEEDS)
    w2, t2, _ = rate(run_B, curve, m, k1, SEEDS)
    w3, t3, _ = rate(run_C, curve, m, k1, SEEDS)
    print(f"{label:<18} {k1:>3} {m:>3} {fmt(wb,tb):>6} {fmt(w2,t2):>6} "
          f"{fmt(w3,t3):>6} {str(2*m+2)+'/'+str(2*m+1)+'/'+str(2*m):>13}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E4: does the K1 wall move?  (m=12, the 2026-07-29 T4 grid)")
print("-" * 78)
print("T4 (2026-07-29) measured, at m=12: 2557 walls at K1~12-16, 2677 at")
print("K1~4-6.  If a reformulation is a real improvement the wall moves out.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
M_FIX = 12
for label, curve, k1_0, m0 in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    print(f"{label}  (m={M_FIX}, K2={k2b}, n={n})")
    hdr = "  " + f"{'method':<8}" + "".join(f"{k:>7}" for k in K1_GRID)
    print(hdr)
    rows_out = {}
    for name, fn in (("base", None), ("E2", run_B), ("E3", run_C)):
        cells = []
        for k1 in K1_GRID:
            if fn is None:
                w, t, _ = rate(lambda c, mm, d, k, s: (lambda r: r and
                               {'ok': r[0]})(L.run_experiment(c, mm, d, k, s)),
                               curve, M_FIX, k1, SEEDS)
            else:
                w, t, _ = rate(fn, curve, M_FIX, k1, SEEDS)
            cells.append(f"{w}/{t}")
        rows_out[name] = cells
        print("  " + f"{name:<8}" + "".join(f"{c:>7}" for c in cells))
    # information-theoretic reference
    info = []
    for k1 in K1_GRID:
        eff = k1 * k2b / n
        if eff >= 1.0:
            info.append("inf")
        else:
            info.append(f"{math.log(n)/math.log(1.0/eff):.1f}")
    print("  " + f"{'m_min':<8}" + "".join(f"{c:>7}" for c in info))
    print()

print("-" * 78)
print("m_min = log(n)/log(n/(K1*K2)) is the information-theoretic minimum")
print("number of signatures; the sweep uses m=12, so any column with")
print("m_min < 12 is information-theoretically satisfiable and a failure")
print("there is a reduction/geometry failure, not a data failure.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E5: Gaussian-heuristic crossing -- is the K1 wall a geometry wall?")
print("-" * 78)
print("For the E2 lattice (dim D = 2m+1) compute")
print("    GH(L) = sqrt(D/(2*pi*e)) * det(L)^(1/D)")
print("and the ratio  gamma = ||v_planted|| / GH(L).  gamma < 1 means the")
print("planted vector is (heuristically) unique-shortest-ish and a good")
print("reduction should find it; gamma > 1 means exponentially many lattice")
print("vectors are shorter than the target and NO amount of reduction helps.\n")


def log_det_from_basis(basis):
    """sum log||b*_i|| via float Gram-Schmidt on a list-of-lists basis."""
    gs = []
    tot = 0.0
    for i in range(len(basis)):
        v = [float(x) for x in basis[i]]
        for g in gs:
            ng = sum(x * x for x in g)
            if ng == 0:
                continue
            c = sum(basis[i][t] * g[t] for t in range(len(v))) / ng
            for t in range(len(v)):
                v[t] -= c * g[t]
        gs.append(v)
        nv = math.sqrt(sum(x * x for x in v))
        if nv > 0:
            tot += math.log(nv)
    return tot


def gamma_E2(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_lattice_nod(sigs, n, lam, k1_bound, k2_bound)
    nrows, ncols = 2 * m + 2, 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(ncols)] for i in range(nrows)]
    basis = [r for r in rows if any(r)]
    D = len(basis)
    ld = log_det_from_basis(basis)
    gh = math.sqrt(D / (2 * math.pi * math.e)) * math.exp(ld / D)
    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    pn = norm(pv[:m] + pv[m + 1:])
    return {'ok': True, 'gamma': pn / gh, 'gh': gh, 'pn': pn, 'D': D}


for label, curve, k1_0, m0 in hist[1:]:
    p, b, n, lam, G = curve
    print(f"{label}  (m={M_FIX})")
    print("  " + f"{'K1':<8}" + "".join(f"{k:>8}" for k in K1_GRID))
    gam = []
    for k1 in K1_GRID:
        vals = []
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = gamma_E2(curve, M_FIX, d_trial, k1, seed)
            if r:
                vals.append(r['gamma'])
        gam.append(sum(vals) / len(vals) if vals else float('nan'))
    print("  " + f"{'gamma':<8}" + "".join(f"{g:>8.3f}" for g in gam))
    cells = []
    for k1 in K1_GRID:
        w, t, _ = rate(run_B, curve, M_FIX, k1, SEEDS)
        cells.append(f"{w}/{t}")
    print("  " + f"{'E2':<8}" + "".join(f"{c:>8}" for c in cells))
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E6: BKZ-beta sweep at the wall -- can better reduction cross it?")
print("-" * 78)
print("If the wall is an approximation-factor wall, raising beta moves it.")
print("If it is the gamma>1 geometry ceiling from E5, beta does nothing.\n")

BETAS = [2, 20, 30, 40, 50]
for label, curve, k1_0, m0 in hist[1:]:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    # pick the two K1 values straddling that curve's wall
    probe = [8, 12, 16] if '2557' in label else [4, 6, 8]
    print(f"{label}  (m={M_FIX})")
    print("  " + f"{'K1':<6}" + "".join(f"{'b='+str(bb):>8}" for bb in BETAS)
          + f"{'gamma':>9}")
    for k1 in probe:
        cells = []
        for bb in BETAS:
            wins = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_trial, M_FIX, n, lam, p, k1, k2b, seed)
                if len(sigs) < M_FIX:
                    continue
                Mx = build_glv_lattice(sigs, n, lam, k1, k2b)
                dim = 2 * M_FIX + 2
                Ax = IntegerMatrix.from_matrix(Mx)
                if bb <= 2:
                    LLL.reduction(Ax)
                else:
                    BKZ.reduction(Ax, BKZ.Param(bb))
                red = [[Ax[i][j] for j in range(dim)] for i in range(dim)]
                _, _, _, S_KAN = scales(n, k1, k2b)
                wins += recover_d(red, M_FIX, n, S_KAN, d_trial) is not None
            cells.append(f"{wins}/{len(SEEDS)}")
        gv = []
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = gamma_E2(curve, M_FIX, d_trial, k1, seed)
            if r:
                gv.append(r['gamma'])
        g = sum(gv) / len(gv) if gv else float('nan')
        print("  " + f"{k1:<6}" + "".join(f"{c:>8}" for c in cells)
              + f"{g:>9.3f}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E7: gamma as a function of m -- does more data ever cross back?")
print("-" * 78)
print("2026-07-29 T4b found that at K1=8 the lam*=0.07 curve is NOT rescued by")
print("m = 8/12/16/24/32 (0,0,1,0,1 of 5).  If gamma(m) decreases in m there")
print("must be an m where it works, and T4b just did not go far enough; if")
print("gamma(m) increases, the K1 wall is permanent at every m and Phase 2 is")
print("at its ceiling.\n")

M_GRID = [6, 8, 12, 16, 24, 32, 48, 64]
for label, curve, k1_0, m0 in hist[1:]:
    p, b, n, lam, G = curve
    for k1 in (8,):
        print(f"{label}  (K1={k1}, K2={math.isqrt(n)+1})")
        print("  " + f"{'m':<8}" + "".join(f"{mm:>8}" for mm in M_GRID))
        gam, suc = [], []
        for mm in M_GRID:
            vals = []
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = gamma_E2(curve, mm, d_trial, k1, seed)
                if r:
                    vals.append(r['gamma'])
            gam.append(sum(vals) / len(vals) if vals else float('nan'))
            w, t, _ = rate(run_B, curve, mm, k1, SEEDS)
            suc.append(f"{w}/{t}")
        print("  " + f"{'gamma':<8}" + "".join(f"{g:>8.3f}" for g in gam))
        print("  " + f"{'E2':<8}" + "".join(f"{c:>8}" for c in suc))
        print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP E8: what IS lambda_1 once the d column is gone, and where does the")
print("        planted vector sit in the reduced basis?")
print("-" * 78)
print("gamma is nearly curve-independent (E5: 1.288 vs 1.291) yet the two")
print("curves behave completely differently, and at m>=48 gamma<1 for both")
print("while 2677 still fails.  So the obstruction is the SHAPE of the GS")
print("profile, not det/dim.  Measure it directly.\n")

print(f"{'curve':<18} {'m':>4} {'||pv||':>10} {'lambda_1':>10} {'lam_1 = mu?':>12} "
      f"{'rank(pv)':>10} {'#|last|=SK':>11} {'d found':>8}")
for label, curve, k1_0, m0 in hist[1:]:
    p, b, n, lam, G = curve
    K1, K2 = 8, math.isqrt(n) + 1
    S_K1, _, S_K2, S_KAN = scales(n, K1, K2)
    mu, wvec = L.lambda_block_mu(n, lam, S_K1, S_K2)
    for mm in (32, 64):
        d = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d, mm, n, lam, p, K1, K2, 42)
        M = build_lattice_nod(sigs, n, lam, K1, K2)
        nr, nc = 2 * mm + 2, 2 * mm + 1
        A = IntegerMatrix.from_matrix(M)
        U = IntegerMatrix.identity(nr)
        LLL.reduction(A, U)
        rows = [[A[i][j] for j in range(nc)] for i in range(nr)]
        pv = planted_vector(sigs, d, n, K1, K2)
        pn = norm(pv[:mm] + pv[mm + 1:])
        nz = sorted(norm(r) for r in rows if any(r))
        kan = [i for i in range(nr) if abs(A[i][nc - 1]) == S_KAN]
        cands = set()
        for i in kan:
            sg = 1 if A[i][nc - 1] > 0 else -1
            cands.add((sg * U[i][mm]) % n)
        rank = sum(1 for x in nz if x < pn) + 1
        print(f"{label:<18} {mm:>4} {pn:>10.3e} {nz[0]:>10.3e} "
              f"{str(abs(nz[0]-mu) < 1.0):>12} {str(rank)+'/'+str(len(nz)):>10} "
              f"{len(kan):>11} {str(d in cands):>8}")

print("\nlambda-block invariants at K1=8 (exact Gauss reduction of the 2D block):")
print(f"{'curve':<18} {'lam*':>8} {'mu':>9} {'blockdet/mu':>12} {'skew':>7} "
      f"{'wrap=lam*K2/n':>14}")
for label, curve, k1_0, m0 in hist[1:]:
    p, b, n, lam, G = curve
    K1, K2 = 8, math.isqrt(n) + 1
    S_K1, _, S_K2, _ = scales(n, K1, K2)
    mu, _w = L.lambda_block_mu(n, lam, S_K1, S_K2)
    bdet = n * S_K1 * S_K2
    print(f"{label:<18} {lam_star(lam, n):>8.4f} {mu:>9.0f} {bdet/mu:>12.0f} "
          f"{(bdet/mu)/mu:>7.2f} {lam*K2/n:>14.2f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E9: reconciling 2026-07-29 T3 (lam* irrelevant) with T4 (lam* shifts")
print("        the wall by 3x) -- the two-regime picture")
print("-" * 78)
print("T3 swept 20 curves at eff = K1*K2/n in {0.05, 0.15, 0.25}.  Map each")
print("eff to the gamma it induces at m=12 on a representative 12-bit curve.\n")

_lbl, _cur, _k10, _m0 = hist[1]
_p, _b, _n, _lam, _G = _cur
_K2 = math.isqrt(_n) + 1
print(f"{'eff':>6} {'K1':>5} {'gamma(m=12)':>12}   T3 outcome (2026-07-29)")
T3 = {0.05: "19/20 curves 5/5, lam* range [0.007,0.491] = the WHOLE range",
      0.15: "3/20 curves,     lam* range [0.318,0.482] = high-lam* only",
      0.25: "0/20 curves"}
for eff in (0.05, 0.15, 0.25):
    k1 = max(2, round(eff * _n / _K2))
    vals = []
    for seed in SEEDS:
        dt = random.Random(seed + 7777).randint(1, _n - 1)
        r = gamma_E2(_cur, 12, dt, k1, seed)
        if r:
            vals.append(r['gamma'])
    g = sum(vals) / len(vals) if vals else float('nan')
    print(f"{eff:>6.2f} {k1:>5} {g:>12.3f}   {T3[eff]}")
print("\nReading: lam* is irrelevant when gamma is small (everything works) and")
print("irrelevant when gamma is large (nothing works); it decides ONLY inside")
print("the transition band gamma ~ 1.1-1.8.  T3 and T4 are both correct and")
print("were measured on opposite sides of that band.")

print("=" * 78)
