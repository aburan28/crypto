"""
GLV-HNP Phase 2, Thread 23 (part 2): the eff wall.

Part 1 (`glv_hnp_phase2_reformulation.py`) showed:
  * F1/F2 fix the lattice geometry exactly as Thread 23 predicted -- the
    parasite n*S_D*e_m is gone and the planted vector becomes lambda_1
    (pv/b1 = 1.000) resp. a BDD error shorter than lambda_1 (0.55-0.77) --
    but the K1 wall does NOT move.  Pre-registered falsifier: second branch.
  * Exp D: the Gaussian-heuristic ratio R does not separate success from
    failure.  R in [0.73, 1.10] over winners, [0.89, 2.19] over losers.
    Instead the wall looked like a pure eff = K1*K2/n threshold sitting
    between eff = 0.079 (wins at every m) and eff = 0.118 (loses at every m),
    with NO m-dependence out to m = 48.

Two things Part 1 left open, answered here:

  F  Harness cross-check.  Exp A did not reproduce the 2026-07-29 Exp T4 table
     at its two margin cells (2557 K1=12: got 1/5 vs 4/5 logged; 2677 K1=6: got
     0/5 vs 2/5 logged).  Call the ORIGINAL module's own `success_rate` on the
     same grid to find out whether the discrepancy is in this session's
     reimplementation or in the logged table.

     *** SUPERSEDED by glv_hnp_phase2_t4_m_identification.py (same session). ***
     Exp F holds m fixed at the module's HIST values (8 and 10).  Sweeping m
     shows the logged T4 table is reproduced exactly, and only, at m = 12 --
     a parameter the 2026-07-29 entry omitted.  The logged numbers are right;
     Exp A and Exp F simply used the wrong m.  Exp F is retained because its
     BKZ(20) column is still a valid measurement.

  G  Universality of the eff wall.  Part 1's Exp E only probed eff <= 0.063 on
     the 17- and 20-bit curves, i.e. entirely BELOW the wall, so it could not
     test whether the wall is at the same eff for every n.  Sweep eff through
     the wall on 12-, 17- and 20-bit curves at two values of m.

  H  Does stronger reduction move the wall?  BKZ(20) and BKZ(40) on F0 and F1
     at the wall.  The 2026-07-26 entry claimed BKZ(40) also fails, but that
     was measured at one (curve, K1, m) point and on F0 only -- i.e. with the
     parasite still present.

Run: python3 glv_hnp_phase2_eff_wall.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ, GSO, FPLLL

# Reuse Part 1's helpers verbatim (everything above its experiment banner).
_SRC = open(__file__.replace('eff_wall', 'reformulation')).read()
_HEAD = _SRC.split('# ' + '=' * 75)[0]
exec(compile(_HEAD, 'glv_hnp_phase2_reformulation.py', 'exec'), globals())

SEEDS = [42, 1234, 9999, 555, 31337]

print("=" * 78)
print("Thread 23 part 2 -- the eff wall (F / G / H)")
print("=" * 78)

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP F -- harness cross-check against glv_hnp_phase2_lambda_threshold.py")
print("-" * 78)
print("Runs the ORIGINAL 2026-07-29 module's own success_rate() on the T4 grid.")
print("If it agrees with Part 1's Exp A, the logged T4 table is what differs.\n")

_ORIG = open('glv_hnp_phase2_lambda_threshold.py').read()
_ORIG_HEAD = _ORIG.split('# ' + '=' * 75)[0]
orig = {}
exec(compile(_ORIG_HEAD, 'glv_hnp_phase2_lambda_threshold.py', 'exec'), orig)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
T4_LOGGED = {
    "12-bit/2557": [5, 5, 5, 5, 5, 4, 1, 0],
    "12-bit/2677": [5, 5, 5, 2, 0, 0, 0, 0],
}
HIST = [("12-bit/2557", 2557, 2, 2659, 1755, 8),
        ("12-bit/2677", 2677, 2, 2647, 185, 10)]

print(f"{'curve':<14} {'source':<26} " + " ".join(f"{'K1='+str(k):>6}" for k in K1_GRID))
for label, p, b, n, lam, m in HIST:
    G = orig['find_generator'](p, b, n)
    cur = (p, b, n, lam, G)
    orig_row, bkz_row = [], []
    for k1 in K1_GRID:
        w, t, _ = orig['success_rate'](cur, m, k1, SEEDS)
        orig_row.append(w)
        wb, tb, _ = orig['success_rate'](cur, m, k1, SEEDS,
                                         use_bkz=True, bkz_beta=20)
        bkz_row.append(wb)
    print(f"{label:<14} {'ORIGINAL module, LLL':<26} "
          + " ".join(f"{v:>6}" for v in orig_row))
    print(f"{label:<14} {'ORIGINAL module, BKZ(20)':<26} "
          + " ".join(f"{v:>6}" for v in bkz_row))
    print(f"{label:<14} {'logged T4 table (07-29)':<26} "
          + " ".join(f"{v:>6}" for v in T4_LOGGED[label]))
    print(f"{label:<14} {'agrees with logged?':<26} "
          f"LLL={orig_row == T4_LOGGED[label]}  "
          f"BKZ20={bkz_row == T4_LOGGED[label]}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP G -- is the eff wall at the same eff for every n?")
print("-" * 78)
print("Part 1 bracketed the 12-bit wall to eff in (0.079, 0.118), m-independent")
print("out to m=48.  Sweep eff THROUGH that band on three curve sizes.\n")

EFF_TARGETS = [0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.24]

curves = []
G0 = find_generator(2677, 2, 2647)
curves.append(("12-bit", (2677, 2, 2647, 185, G0)))
for lo, hi, tag in [(2 ** 16, 2 ** 17, "17-bit"), (2 ** 19, 2 ** 20, "20-bit")]:
    c = find_curve_in_range(lo, hi)
    if c is not None:
        curves.append((tag, c))

for tag, c in curves:
    p, b, n, lam, G = c
    k2b = math.isqrt(n) + 1
    ft = 'd' if n.bit_length() <= 18 else 'mpfr'
    if ft == 'mpfr':
        FPLLL.set_precision(212)
    print(f"{tag}  p={p} n={n} ({n.bit_length()}b) lam*={lam_star(lam,n):.3f} "
          f"K2={k2b}  (Babai {ft})")
    print(f"  {'m':>3} | " + " | ".join(f"{'eff~'+format(e,'.2f'):>16}"
                                        for e in EFF_TARGETS))
    for m in (12, 24):
        cells = []
        for e in EFF_TARGETS:
            k1 = max(2, round(e * n / k2b))
            w, t, _ = rate(c, m, k1, SEEDS, kinds=('F0', 'F1', 'F2'),
                           float_type=ft)
            eff = k1 * k2b / n
            cells.append(f"K1={k1} {w['F0']}{w['F1']}{w['F2']} e={eff:.3f}")
        print(f"  {m:>3} | " + " | ".join(f"{s:>16}" for s in cells))
    print("  (cell = K1, then F0/F1/F2 wins out of 5, then realised eff)\n")

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP H -- does stronger reduction move the wall, once the parasite is gone?")
print("-" * 78)
print("12-bit/2677, m=12, K1 across the wall.  LLL vs BKZ(20) vs BKZ(40),")
print("on F0 (parasite present) and F1 (parasite removed).\n")


def rate_bkz(curve, m, k1_bound, seeds, kind, beta):
    """Success rate for F0/F1 under BKZ(beta) instead of LLL."""
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2b)
    wins = 0
    for seed in seeds:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
        if len(sigs) < m:
            continue
        rows = build_F0(sigs, n, lam, k1_bound, k2b) if kind == 'F0' \
            else build_F1(sigs, n, lam, k1_bound, k2b)
        A = IntegerMatrix.from_matrix(rows)
        if beta:
            BKZ.reduction(A, BKZ.Param(beta))
        else:
            LLL.reduction(A)
        red = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
        red = [r for r in red if any(r)]
        if kind == 'F0':
            ok = recover_F0(red, m, n, S_KAN, d_secret) is not None
        else:
            ok = recover_F1(red, sigs, m, n, lam, k1_bound, k2b,
                            d_secret) is not None
        wins += bool(ok)
    return wins


cur = (2677, 2, 2647, 185, G0)
K1_H = [4, 5, 6, 7, 8, 10, 12]
print(f"  {'K1':>4} {'eff':>7} | " + " ".join(
    f"{k + '/' + str(bb):>9}" for k in ('F0', 'F1') for bb in (0, 20, 40)))
for k1 in K1_H:
    vals = [rate_bkz(cur, 12, k1, SEEDS, k, bb)
            for k in ('F0', 'F1') for bb in (0, 20, 40)]
    print(f"  {k1:>4} {k1*52/2647:>7.4f} | " + " ".join(f"{v:>9}" for v in vals))
print("  (columns: formulation / BKZ block size; 0 = plain LLL; wins out of 5)")

print("\nDone.")
