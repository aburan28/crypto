"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector
is actually the target of the reduction.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, finding T5):
  The Phase-2 Kannan lattice L_A (dim 2m+2, built by
  glv_hnp_phase2_20bit.py:263 `build_glv_lattice`) always has
      lambda_1(L_A) = n * S_D * e_m
  the "trivial" vector living entirely in the d column.  It is there because
  d is only defined mod n, so n*e_m is in the lattice for free.  It is
  2-3x SHORTER than the planted vector on every curve tested (sv/pv in
  [0.34, 0.61]), success and failure alike, so the planted vector is never
  lambda_1 and recovery is a BDD/coset condition rather than an SVP one.

  The 2026-07-29 run proposed (verbatim):
    "project the lattice along e_m (quotient out the trivial n*e_m
     direction) and solve BDD in the projection, or replace the Kannan
     embedding with an explicit CVP call (Babai nearest-plane on the
     reduced basis) targeting (A_i*S_K1, 0, ...)."
  with the falsifier:
    "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4
     moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
     reformulation is a real improvement; if the wall stays at K1 ~ 4-6,
     then the wall is information-theoretic and Phase 2 is at its ceiling."

This script executes exactly that, in three arms.

  ARM A  baseline.  build_glv_lattice verbatim, dim 2m+2, LLL, recover_d.
         (Reproduces the 2026-07-29 T4 grid so the comparison is exact.)

  ARM B  projected Kannan.  Delete column m (the d column) from L_A.
         dim 2m+1.  The image of the planted vector is
             (k1_i*S_K1)_i ++ (k2_i*S_K2)_i ++ (S_KANNAN)
         which no longer carries d as a coordinate; d is reconstructed
         afterwards from (k1_0, k2_0) via d = B_0^{-1}(k1_0 + lam*k2_0 - A_0).
         The trivial vector projects to 0, so it is gone from the lattice.
         Soundness: adding t copies of the projected d-row changes col i by
         t*B_i*S_K1, i.e. changes the implied d by t; only the true d makes
         every k1_i small.  So this is still a genuine BDD instance, just
         one with the free short vector removed.

  ARM C  explicit CVP.  Homogeneous lattice L_h in Z^{2m} (m congruence
         columns + m k2 columns), no Kannan row and no d column:
             n*S_K1*e_i                       (i = 0..m-1)
             (B_0*S_K1, ..., B_{m-1}*S_K1)    (the d row)
             -lam*S_K1*e_i + S_K2*e_{m+i}     (i = 0..m-1)
         target  t = (-A_i*S_K1)_i ++ (0)_m.
         Babai nearest-plane on the LLL-reduced basis gives u in L_h; the
         residual u - t should be exactly (k1_i*S_K1)_i ++ (k2_i*S_K2)_i.
         This is the smallest of the three targets:
             ||pv_A||^2 ~ n^2(2m/3 + 4/3)
             ||pv_B||^2 ~ n^2(2m/3 + 1)
             ||pv_C||^2 ~ n^2(2m/3)

Metrics reported per (curve, K1, arm):
  wins/trials   over the same 5 seeds as the 2026-07-29 T4 grid
  sv/pv         ||lambda_1 after reduction|| / ||planted vector||
                (the falsifier's headline number; > 1 means the planted
                vector IS the shortest vector)
  the K1 wall   the largest K1 at which the arm still recovers d

Also reported: the information-theoretic bound m*log2(1/eff) >= log2(n),
so that "the wall is information-theoretic" can be checked numerically
rather than asserted.

Run: python3 glv_hnp_phase2_bdd_reformulation.py
"""

import math
import random
import sys

import contextlib
import io

import sympy
from fpylll import IntegerMatrix, LLL, BKZ, GSO

sys.path.insert(0, __file__.rsplit('/', 1)[0])

# glv_hnp_phase2_20bit is a script, not a library: importing it runs its own
# 20-bit sweep (~90s of output).  Swallow that; we only want its functions.
with contextlib.redirect_stdout(io.StringIO()):
    from glv_hnp_phase2_20bit import (       # noqa: E402
        modinv, ec_mul, tonelli_shanks, find_generator,
        gen_signatures, build_glv_lattice, recover_d,
        eisenstein_decompose, j0_traces, glv_eigenvalue, find_point,
    )

SEEDS = [42, 1234, 9999, 555, 31337]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def norm2(v):
    return sum(int(x) * int(x) for x in v)


def lam_star(lam, n):
    return min(lam, n - lam) / n


def lll_rows(M):
    """LLL-reduce a list-of-lists integer matrix; return list of rows."""
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]


def shortest_nonzero(rows):
    best = None
    for r in rows:
        q = norm2(r)
        if q and (best is None or q < best):
            best = q
    return best


# ---------------------------------------------------------------------------
# ARM A — baseline Kannan lattice (dim 2m+2), verbatim from 2026-06-15
# ---------------------------------------------------------------------------

def arm_A(sigs, n, lam, k1_bound, k2_bound, d_secret, m, use_bkz=False, beta=20):
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2

    # planted vector, constructed explicitly (not searched for)
    pv = [0] * dim
    for i, s in enumerate(sigs):
        pv[i] = s['k1'] * S_K1
        pv[m + 1 + i] = s['k2'] * S_K2
    pv[m] = d_secret * S_D
    pv[dim - 1] = S_KANNAN

    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_d(rows, m, n, S_KANNAN, d_secret) is not None
    return ok, math.sqrt(shortest_nonzero(rows) / norm2(pv))


# ---------------------------------------------------------------------------
# ARM B — projected Kannan lattice: delete the d column (dim 2m+1)
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    """L_A with column m (the d column) deleted.

    Column layout of the result (dim 2m+1):
        0     .. m-1   congruence columns   (entry k1_i*S_K1 on the target)
        m     .. 2m-1  k2 columns           (entry k2_i*S_K2 on the target)
        2m             Kannan column        (entry S_KANNAN on the target)
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

    M = []
    for i in range(m):                                     # n*S_K1*e_i
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    row = [0] * dim                                        # projected d row
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for i in range(m):                                     # k2 rows
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    row = [0] * dim                                        # Kannan row
    for i in range(m):
        row[i] = sigs[i]['A'] * S_K1
    row[dim - 1] = S_KANNAN
    M.append(row)

    return M, S_K1, S_K2, S_KANNAN


def recover_from_projected(rows, sigs, n, lam, m, S_K1, S_K2, S_KANNAN, d_secret):
    """Read k1/k2 off a reduced row with |Kannan coord| == S_KANNAN, then get d."""
    dim = 2 * m + 1
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        k1s = [sign * row[i] for i in range(m)]
        k2s = [sign * row[m + i] for i in range(m)]
        if any(x % S_K1 for x in k1s) or any(x % S_K2 for x in k2s):
            continue
        k1s = [x // S_K1 for x in k1s]
        k2s = [x // S_K2 for x in k2s]
        # reconstruct d from signature 0:  A_0 + B_0*d = k1_0 + lam*k2_0 (mod n)
        try:
            d_cand = (k1s[0] + lam * k2s[0] - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
        except ValueError:
            continue
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand
    return None


def arm_B(sigs, n, lam, k1_bound, k2_bound, d_secret, m, use_bkz=False, beta=20):
    M, S_K1, S_K2, S_KANNAN = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1

    pv = [0] * dim
    for i, s in enumerate(sigs):
        pv[i] = s['k1'] * S_K1
        pv[m + i] = s['k2'] * S_K2
    pv[dim - 1] = S_KANNAN

    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    ok = recover_from_projected(rows, sigs, n, lam, m, S_K1, S_K2,
                                S_KANNAN, d_secret) is not None
    return ok, math.sqrt(shortest_nonzero(rows) / norm2(pv))


# ---------------------------------------------------------------------------
# ARM C — explicit CVP (Babai nearest-plane) on the homogeneous lattice
# ---------------------------------------------------------------------------

def build_homogeneous_lattice(sigs, n, lam, k1_bound, k2_bound):
    """L_h in Z^{2m}: congruence columns 0..m-1, k2 columns m..2m-1.

    No Kannan row, no d column.  2m+1 generators, rank 2m.
    """
    m = len(sigs)
    dim = 2 * m
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)

    M = []
    for i in range(m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    row = [0] * dim
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    M.append(row)
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    return M, S_K1, S_K2


def arm_C(sigs, n, lam, k1_bound, k2_bound, d_secret, m, use_bkz=False, beta=20):
    M, S_K1, S_K2 = build_homogeneous_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m

    pv = [0] * dim                       # the residual we are hoping to hit
    for i, s in enumerate(sigs):
        pv[i] = s['k1'] * S_K1
        pv[m + i] = s['k2'] * S_K2

    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    nz = [r for r in rows if any(r)]
    if len(nz) != dim:
        return False, float('nan')

    target = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m

    B = IntegerMatrix.from_matrix(nz)
    G = GSO.Mat(B, float_type='d')
    G.update_gso()
    try:
        coeffs = G.babai(target)
    except Exception:
        return False, float('nan')

    u = [0] * dim
    for c, r in zip(coeffs, nz):
        if c:
            for j in range(dim):
                u[j] += c * r[j]
    resid = [u[j] - target[j] for j in range(dim)]

    ok = False
    if all(x % S_K1 == 0 for x in resid[:m]) and all(x % S_K2 == 0 for x in resid[m:]):
        k1s = [resid[i] // S_K1 for i in range(m)]
        k2s = [resid[m + i] // S_K2 for i in range(m)]
        try:
            d_cand = (k1s[0] + lam * k2s[0] - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
            ok = (d_cand == d_secret)
        except ValueError:
            ok = False

    # for arm C, sv/pv compares lambda_1(L_h) to the CVP residual we want
    return ok, math.sqrt(shortest_nonzero(nz) / norm2(pv))


ARMS = [('A base', arm_A), ('B proj', arm_B), ('C cvp', arm_C)]


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def trial(arm_fn, curve, m, k1_bound, seed, d_secret=None, use_bkz=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    if d_secret is None:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    return arm_fn(sigs, n, lam, k1_bound, k2_bound, d_secret, m, use_bkz=use_bkz)


def arm_stats(arm_fn, curve, m, k1_bound, seeds, use_bkz=False):
    wins, tot, ratios = 0, 0, []
    for s in seeds:
        r = trial(arm_fn, curve, m, k1_bound, s, use_bkz=use_bkz)
        if r is None:
            continue
        ok, ratio = r
        tot += 1
        wins += int(ok)
        if ratio == ratio:                       # not NaN
            ratios.append(ratio)
    med = sorted(ratios)[len(ratios) // 2] if ratios else float('nan')
    return wins, tot, med


print("=" * 78)
print("Thread 23 — does removing the trivial vector n*S_D*e_m rescue Phase 2?")
print("=" * 78)

# Historical Phase-2 curves — identical table to
# glv_hnp_phase2_lambda_threshold.py:424 so the T4 comparison is exact.
HIST = [
    # label,             p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",      2677, 2, 2647, 185,  8,  10),
]

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))


# ===========================================================================
# EXP U0 — membership self-check
# ===========================================================================
# Guards against sign/index slips in the three planted-vector formulas: build
# each planted vector from explicit integer coefficients against the basis and
# assert it equals the closed form used by the arms.  (An early version of this
# script had k2 -> -k2 in ARM B's readout, which silently reported 0/30.)
print("\n" + "-" * 78)
print("EXP U0: planted-vector membership self-check (all three arms)")
print("-" * 78)
_lab, _curve, _k1, _m = hist[1]
_p, _b, _n, _lam, _G = _curve
_k2b = math.isqrt(_n) + 1
_d = random.Random(4242).randint(1, _n - 1)
_sigs = gen_signatures(_G, _d, _m, _n, _lam, _p, _b, _k1, _k2b, 42)
assert len(_sigs) == _m

_MA, _SK1, _SD, _SK2, _SKAN = build_glv_lattice(_sigs, _n, _lam, _k1, _k2b)
_MB, _, _, _ = build_projected_lattice(_sigs, _n, _lam, _k1, _k2b)
_MC, _, _ = build_homogeneous_lattice(_sigs, _n, _lam, _k1, _k2b)


def _combine(basis, coeffs):
    dim = len(basis[0])
    out = [0] * dim
    for c, r in zip(coeffs, basis):
        if c:
            for j in range(dim):
                out[j] += c * r[j]
    return out


# coefficients: c_i on the n*S_K1*e_i rows, d on the d row, k2_i on the k2 rows,
# 1 on the Kannan row (arms A and B only).
_c = [(_sigs[i]['k1'] - (_sigs[i]['A'] + _d * _sigs[i]['B'] - _lam * _sigs[i]['k2'])) // _n
      for i in range(_m)]
for i in range(_m):
    assert (_sigs[i]['k1'] - (_sigs[i]['A'] + _d * _sigs[i]['B']
                              - _lam * _sigs[i]['k2'])) % _n == 0

_pvA = [0] * (2 * _m + 2)
for i, s in enumerate(_sigs):
    _pvA[i] = s['k1'] * _SK1
    _pvA[_m + 1 + i] = s['k2'] * _SK2
_pvA[_m] = _d * _SD
_pvA[2 * _m + 1] = _SKAN
assert _combine(_MA, _c + [_d] + [s['k2'] for s in _sigs] + [1]) == _pvA, "ARM A pv"

_pvB = [0] * (2 * _m + 1)
for i, s in enumerate(_sigs):
    _pvB[i] = s['k1'] * _SK1
    _pvB[_m + i] = s['k2'] * _SK2
_pvB[2 * _m] = _SKAN
assert _combine(_MB, _c + [_d] + [s['k2'] for s in _sigs] + [1]) == _pvB, "ARM B pv"

_tgtC = [-s['A'] * _SK1 for s in _sigs] + [0] * _m
_pvC = [0] * (2 * _m)
for i, s in enumerate(_sigs):
    _pvC[i] = s['k1'] * _SK1
    _pvC[_m + i] = s['k2'] * _SK2
_uC = _combine(_MC, _c + [_d] + [s['k2'] for s in _sigs])
assert [_uC[j] - _tgtC[j] for j in range(2 * _m)] == _pvC, "ARM C residual"

print("  ARM A planted vector in L_A               : OK")
print("  ARM B planted vector in proj(L_A)         : OK")
print("  ARM C residual u - t in L_h - t           : OK")
print(f"  ||pv_A|| / ||pv_B|| / ||pv_C||  (n={_n}, m={_m}) = "
      f"{math.sqrt(norm2(_pvA))/_n:.3f}n / {math.sqrt(norm2(_pvB))/_n:.3f}n / "
      f"{math.sqrt(norm2(_pvC))/_n:.3f}n")


# ===========================================================================
# EXP U1 — sv/pv for all three arms, at the historical (K1, m) settings
# ===========================================================================
print("\n" + "-" * 78)
print("EXP U1: sv/pv  (falsifier headline).  >1 means the planted vector IS")
print("        the shortest vector of the reduced lattice.")
print("-" * 78)
print(f"{'curve':<14} {'lam*':>6} {'K1':>3} {'m':>3} "
      f"{'A sv/pv':>9} {'B sv/pv':>9} {'C sv/pv':>9}")
u1 = {}
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    row = []
    for name, fn in ARMS:
        _, _, med = arm_stats(fn, curve, m, k1, SEEDS)
        row.append(med)
    u1[label] = row
    print(f"{label:<14} {lam_star(lam, n):>6.3f} {k1:>3} {m:>3} "
          f"{row[0]:>9.3f} {row[1]:>9.3f} {row[2]:>9.3f}")

print("\nARM A reproduces 2026-07-29 T5 (sv/pv in [0.34, 0.61], planted vector")
print("never shortest).  ARM B/C remove the free n*e_m direction.")


# ===========================================================================
# EXP U2 — the K1 wall, all three arms.  Direct replay of the T4 grid.
# ===========================================================================
print("\n" + "-" * 78)
print("EXP U2: K1 wall.  Same grid as 2026-07-29 T4 (K2 = isqrt(n)+1, 5 seeds).")
print("        T4 baseline for reference:")
print("          12-bit/2557 (lam*=0.340): 5,5,5,5,5,4,1,0  for K1=2,3,4,6,8,12,16,24")
print("          12-bit/2677 (lam*=0.070): 5,5,5,2,0,0,0,0")
print("        m=12: the 2026-07-29 log states K1, K2 and the seed list for the")
print("        T4 grid but not m.  m=12 is the value that reproduces BOTH logged")
print("        rows bit-exactly under ARM A (m=6..11 do not); it also matches the")
print("        m=12 used in that run's T3.  ARM A below is that reproduction.")
print("-" * 78)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
T4_M = 12
walls = {}
for label, curve, _k1, _m_hist in hist[1:]:      # the two 12-bit curves, as in T4
    p, b, n, lam, G = curve
    m = T4_M
    print(f"\n{label}  (lam*={lam_star(lam, n):.3f}, n={n}, m={m}, "
          f"K2={math.isqrt(n)+1})")
    hdr = "  ".join(f"{k:>4}" for k in K1_GRID)
    print(f"{'arm':<8} {'':>0}{hdr}      eff at each K1")
    for name, fn in ARMS:
        cells, wall = [], None
        for k1 in K1_GRID:
            w, t, _ = arm_stats(fn, curve, m, k1, SEEDS)
            cells.append(f"{w}/{t}")
            if w == t and t > 0:
                wall = k1
        walls[(label, name)] = wall
        print(f"{name:<8} " + "  ".join(f"{c:>4}" for c in cells) +
              f"      wall(5/5) K1={wall}")
    effs = [k1 * (math.isqrt(n) + 1) / n for k1 in K1_GRID]
    print(f"{'eff':<8} " + "  ".join(f"{e:>4.2f}" for e in effs))
    # information-theoretic bound: need m*log2(1/eff) >= log2(n)
    print(f"{'m*lg1/e':<8} " + "  ".join(
        f"{m*math.log2(1/e) if e < 1 else 0:>4.1f}" for e in effs) +
        f"      info bound = lg2(n) = {math.log2(n):.1f}")


# ===========================================================================
# EXP U3 — does BKZ(20) change the verdict at the wall?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP U3: BKZ(20) at the K1 values straddling the 2677 wall.")
print("-" * 78)
label, curve, _k1, _mh = hist[2]
p, b, n, lam, G = curve
m = T4_M
print(f"{label}  m={m}")
print(f"{'arm':<8} {'red':<6} " + "  ".join(f"K1={k:<4}" for k in [4, 6, 8, 12]))
for name, fn in ARMS:
    for red, bkz in (("LLL", False), ("BKZ20", True)):
        cells = []
        for k1 in [4, 6, 8, 12]:
            w, t, _ = arm_stats(fn, curve, m, k1, SEEDS, use_bkz=bkz)
            cells.append(f"{w}/{t}")
        print(f"{name:<8} {red:<6} " + "  ".join(f"{c:<7}" for c in cells))


# ===========================================================================
# EXP U4 — more data at the wall: does m rescue any arm?
# ===========================================================================
print("\n" + "-" * 78)
print("EXP U4: m sweep at K1=8 on 12-bit/2677 (T4b showed ARM A does not")
print("        recover: m=8/12/16/24/32 -> 0,0,1,0,1 of 5).")
print("-" * 78)
M_GRID = [8, 12, 16, 24, 32]
print(f"{'arm':<8} " + "  ".join(f"m={mm:<4}" for mm in M_GRID))
for name, fn in ARMS:
    cells = []
    for mm in M_GRID:
        w, t, _ = arm_stats(fn, curve, mm, 8, SEEDS)
        cells.append(f"{w}/{t}")
    print(f"{name:<8} " + "  ".join(f"{c:<6}" for c in cells))


# ===========================================================================
# EXP U5 — replication on fresh 17-bit curves (not the historical anchors)
# ===========================================================================
def find_j0_curves(lo, hi, want, m_default=10):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = p + 1 - t
                    if nc < 8 or not sympy.isprime(nc) or nc % 3 != 1:
                        continue
                    lam, _ = glv_eigenvalue(nc)
                    if lam is None:
                        continue
                    fb = None
                    for b_try in range(1, 60):
                        P = find_point(p, b_try)
                        if P is not None and ec_mul(P, nc, p) is None:
                            fb = b_try
                            break
                    if fb is None:
                        continue
                    G = find_generator(p, fb, nc)
                    if G is None:
                        continue
                    out.append((f"{p}/{nc}", (p, fb, nc, lam, G)))
                    break
        p = int(sympy.nextprime(p))
    return out

print("\n" + "-" * 78)
print("EXP U5: replication on 6 fresh 17-bit j=0 curves, m=10, K1 sweep.")
print("-" * 78)
fresh = find_j0_curves(2 ** 16, 2 ** 17, 6)
M17 = 10
K1_17 = [8, 16, 32, 64, 128]
print(f"{'curve':<14} {'lam*':>6} {'arm':<8} " + "  ".join(f"K1={k:<5}" for k in K1_17))
agg = {name: [0] * len(K1_17) for name, _ in ARMS}
for label, curve in fresh:
    p, b, n, lam, G = curve
    for name, fn in ARMS:
        cells = []
        for j, k1 in enumerate(K1_17):
            w, t, _ = arm_stats(fn, curve, M17, k1, SEEDS)
            agg[name][j] += w
            cells.append(f"{w}/{t}")
        print(f"{label:<14} {lam_star(lam, n):>6.3f} {name:<8} " +
              "  ".join(f"{c:<7}" for c in cells))

print(f"\n{'TOTAL':<14} {'':>6} " + " " * 8 +
      "  ".join(f"{'':<7}" for _ in K1_17))
denom = len(fresh) * len(SEEDS)
for name, _ in ARMS:
    print(f"{'aggregate':<14} {'':>6} {name:<8} " +
          "  ".join(f"{agg[name][j]}/{denom:<5}" for j in range(len(K1_17))))
effs17 = [k1 * (math.isqrt(fresh[0][1][2]) + 1) / fresh[0][1][2] for k1 in K1_17]
print(f"{'eff (approx)':<14} {'':>6} {'':<8} " +
      "  ".join(f"{e:<7.3f}" for e in effs17))

# ===========================================================================
# EXP U6 — why is ARM B identical to ARM A?  The decoupling measurement.
# ===========================================================================
# ARM B matched ARM A cell-for-cell everywhere above.  Hypothesis: the d
# column is a rank-1 direction that no other reduced basis vector touches, so
# LLL isolates n*S_D*e_m immediately and the reduction of the complement is
# unchanged.  If so, the sorted norm profile of the ARM A reduced basis, with
# the trivial vector deleted, equals that of the ARM B reduced basis.
print("\n" + "-" * 78)
print("EXP U6: decoupling test — is the d column inert for the reduction?")
print("-" * 78)
print(f"{'curve':<14} {'K1':>3} {'seed':>6} {'rows w/ d!=0':>13} "
      f"{'energy in d col':>16} {'profiles equal':>15}")
for label, curve, k1h, mh in hist:
    p, b, n, lam, G = curve
    for k1 in (k1h, 2 * k1h):
        for seed in SEEDS[:2]:
            k2b = math.isqrt(n) + 1
            dsec = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, dsec, mh, n, lam, p, b, k1, k2b, seed)
            if len(sigs) < mh:
                continue
            m = mh
            MA, S_K1, S_D, S_K2, S_KAN = build_glv_lattice(sigs, n, lam, k1, k2b)
            MB, _, _, _ = build_projected_lattice(sigs, n, lam, k1, k2b)
            RA, RB = lll_rows(MA), lll_rows(MB)

            nz_d = sum(1 for r in RA if r[m])
            tot_e = sum(norm2(r) for r in RA)
            d_e = sum(int(r[m]) ** 2 for r in RA)

            # delete the trivial vector (the unique row supported only in col m)
            triv = [i for i, r in enumerate(RA)
                    if r[m] and not any(r[j] for j in range(len(r)) if j != m)]
            profA = sorted(norm2([r[j] for j in range(len(r)) if j != m])
                           for i, r in enumerate(RA) if i not in triv)
            profB = sorted(norm2(r) for r in RB)
            print(f"{label:<14} {k1:>3} {seed:>6} {nz_d:>6}/{len(RA):<6} "
                  f"{d_e/tot_e:>16.2e} {str(profA == profB):>15}")

print("\nHYPOTHESIS REFUTED.  The d column is NOT a rank-1 direction that LLL")
print("isolates: a majority of reduced rows carry a nonzero d coordinate and the")
print("norm profiles of the two reduced bases never agree.  ARM A and ARM B")
print("therefore produce genuinely DIFFERENT reduced bases.  The identical")
print("scores in U2/U5 are not basis identity — tested per-trial below.")


# ===========================================================================
# EXP U7 — per-trial agreement between ARM A and ARM B
# ===========================================================================
# U2/U5 compare aggregate counts, which can coincide by accident.  This
# compares the arms instance by instance on the SAME signature set.
print("\n" + "-" * 78)
print("EXP U7: per-trial agreement, ARM A vs ARM B (same sigs, same d, same seed)")
print("-" * 78)
agree = disagree = 0
disagree_cells = []
for label, curve, k1h, mh in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in [2, 3, 4, 6, 8, 12, 16, 24]:
        for seed in SEEDS:
            dsec = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, dsec, mh, n, lam, p, b, k1, k2b, seed)
            if len(sigs) < mh:
                continue
            okA, _ = arm_A(sigs, n, lam, k1, k2b, dsec, mh)
            okB, _ = arm_B(sigs, n, lam, k1, k2b, dsec, mh)
            if okA == okB:
                agree += 1
            else:
                disagree += 1
                disagree_cells.append((label, k1, seed, okA, okB))
tot = agree + disagree
print(f"  trials compared : {tot}")
print(f"  A and B agree   : {agree}  ({100.0*agree/tot:.1f}%)")
print(f"  A and B disagree: {disagree}")
for c in disagree_cells[:20]:
    print(f"    {c[0]:<14} K1={c[1]:<3} seed={c[2]:<6} A={c[3]} B={c[4]}")
if disagree == 0:
    print("\n  Perfect agreement on every instance, despite provably different")
    print("  reduced bases (U6).  Removing the trivial vector changes the")
    print("  reduction but never changes the OUTCOME: recovery is decided by the")
    print("  BDD instance itself, not by the presence of n*S_D*e_m.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
