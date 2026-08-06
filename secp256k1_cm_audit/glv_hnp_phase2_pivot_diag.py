"""
GLV-HNP Phase 2, Thread 23 follow-up: two diagnostics on the pivot lattice.

Companion to glv_hnp_phase2_pivot_lattice.py (same run, 2026-08-06).  That
script showed the pivot lattice removes the trivial vector n*S_D*e_d but does
NOT move the K1 wall.  Two loose ends it left:

  E4  The 2026-07-29 T4 grid reported 5/5 for 12-bit/2557 at K1=8 and 4/5 at
      K1=12; the E2 grid here reports 3/5 and 1/5 for the same curve, m and
      seeds.  The only uncontrolled variable is d_secret (T4 did not record
      it).  E4 sweeps d over 8 values to test whether d alone explains the
      gap, i.e. whether the wall position is an instance property rather than
      a curve property.

  E5  In the pivot lattice the shortest reduced vector is no longer the
      trivial vector.  What is it?  E5 decomposes its energy over the
      coordinate groups
          res   = k1_j columns, j = 1..m-1      (m-1 coords)
          k1_0  = the pivot k1 column           (1 coord)
          k2_0  = the pivot k2 column           (1 coord)
          k2    = k2_j columns, j = 1..m-1      (m-1 coords)
          kan   = the Kannan column             (1 coord)
      and reports the Kannan coordinate.  A vector with kan = 0 lies in the
      homogeneous sublattice and carries no information about d.

Run: python3 glv_hnp_phase2_pivot_diag.py
"""

import math

from glv_hnp_phase2_pivot_lattice import (
    find_generator, gen_signatures, scales, lam_star, norm,
    build_glv_lattice, planted_vector, recover_d,
    build_pivot_lattice, pivot_planted_vector, pivot_recover_d,
    reduce_matrix, run_orig, run_pivot, rate,
    SEEDS, HIST,
)

print("=" * 78)
print("Thread 23 follow-up: d-sensitivity of the wall, and what lambda_1 is")
print("=" * 78)

curves = {}
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    curves[label] = ((p, b, n, lam, G), k1, m)

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E4: is the K1-wall position a property of the curve or of d?")
print("-" * 78)
print("Same curve, same m, same 5 seeds; only d_secret varies.")
print("If the wall moves with d, the 2026-07-29 T4 grid and the E2 grid in")
print("glv_hnp_phase2_pivot_lattice.py are both correct -- they used different d.")
print()
K1_GRID = [4, 6, 8, 12, 16]
for label in ("12-bit/2557", "12-bit/2677"):
    cv, k1_hist, m = curves[label]
    p, b, n, lam, G = cv
    print(f"\n  {label}  (m={m}, orig lattice, LLL; entries = #recovered/5)")
    print("    " + f"{'d_secret':>10}" + "".join(f"{'K1='+str(k):>8}" for k in K1_GRID))
    tot = {k: 0 for k in K1_GRID}
    for frac in (0.0501, 0.1237, 0.2500, 0.3319, 0.4213, 0.5771, 0.7182, 0.9033):
        d_sec = max(1, int(round(frac * n)))
        cells = []
        for K in K1_GRID:
            h, _ = rate(lambda s, K=K, d=d_sec: run_orig(cv, m, d, K, seed=s), SEEDS)
            tot[K] += h
            cells.append(f"{h}/5")
        print("    " + f"{d_sec:>10}" + "".join(f"{c:>8}" for c in cells))
    print("    " + f"{'TOTAL':>10}" + "".join(f"{str(tot[k])+'/40':>8}" for k in K1_GRID))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E5: energy decomposition of the shortest pivot-lattice vector")
print("-" * 78)
print("frac_X = (sum of squares over group X) / ||sv||^2.  kan = Kannan coord")
print("in units of S_KANNAN; kan = 0 means the vector says nothing about d.")
print()
print(f"{'curve':<13} {'K1':>3} {'m':>3} {'sv/pv':>7} | {'res':>6} {'k1_0':>6} {'k2_0':>6} "
      f"{'k2':>6} {'kan':>6} | {'kan/S':>6} {'ok':>5}")
for label in ("8-bit/199", "12-bit/2557", "12-bit/2677"):
    cv, k1, m = curves[label]
    p, b, n, lam, G = cv
    for K in (k1, 4, 2):
        K2 = math.isqrt(n) + 1
        d_sec = max(1, int(round(0.4213 * n)))
        sigs = gen_signatures(G, d_sec, m, n, lam, p, K, K2, seed=42)
        S_K1, _, S_K2, S_KAN = scales(n, K, K2)
        M = build_pivot_lattice(sigs, n, lam, K, K2)
        pv = pivot_planted_vector(sigs, n, K, K2)
        red = reduce_matrix(M, use_bkz=False)
        dim = 2 * m + 1
        sv, svn = None, None
        for r in red:
            nr = norm(r)
            if nr > 0 and (svn is None or nr < svn):
                sv, svn = r, nr
        e_res = sum(x * x for x in sv[0:m - 1])
        e_k10 = sv[m - 1] ** 2
        e_k20 = sv[m] ** 2
        e_k2 = sum(x * x for x in sv[m + 1:2 * m])
        e_kan = sv[dim - 1] ** 2
        tot = float(svn * svn)
        ok = pivot_recover_d(red, sigs, n, lam, K, K2, d_sec) is not None
        print(f"{label:<13} {K:>3} {m:>3} {svn/norm(pv):>7.4f} | "
              f"{e_res/tot:>6.3f} {e_k10/tot:>6.3f} {e_k20/tot:>6.3f} "
              f"{e_k2/tot:>6.3f} {e_kan/tot:>6.3f} | "
              f"{sv[dim-1]/S_KAN:>6.2f} {str(ok):>5}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
