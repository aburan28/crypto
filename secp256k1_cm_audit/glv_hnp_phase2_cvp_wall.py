"""
Thread 23c — does exact CVP move the Phase-2 K1 wall?

EXP P3 (glv_hnp_phase2_projected.py) showed one grid cell where exact CVP beat
LLL-Kannan on the lam*=0.070 curve (K1=6: CVP 4/5 vs KAN 0/5) at 5 seeds.  This
script re-measures K1 = 4..8 on that curve at 15 seeds, to decide whether the
cell was noise or a real shift of the wall.

Four methods (all defined in glv_hnp_phase2_projected.py):
  KAN    LLL on the historical Kannan embedding, dim 2m+2
  PROJ   LLL on the same lattice projected along e_m (trivial vector removed)
  BABAI  nearest-plane on the LLL-reduced BDD lattice, dim 2m
  CVP    exact closest-vector (fplll enumeration) on the same BDD lattice

Run: python3 glv_hnp_phase2_cvp_wall.py
"""

import math

from glv_hnp_phase2_projected import (
    find_generator, sweep, gh_margin, lam_star, nu_hat, scales, METHODS,
)

if __name__ == "__main__":
    P, B, N, LAM, M_SIGS = 2677, 2, 2647, 185, 10
    G = find_generator(P, B, N)
    assert G is not None
    assert (LAM * LAM + LAM + 1) % N == 0
    curve = (P, B, N, LAM, G)
    K2 = math.isqrt(N) + 1
    s = scales(N, 2, K2)

    SEEDS = [42, 1234, 9999, 555, 31337, 101, 202, 303, 404, 505,
             606, 707, 808, 909, 111]

    print("=" * 70)
    print("Thread 23c — exact CVP vs LLL-Kannan at the K1 wall")
    print("=" * 70)
    print(f"curve 12-bit/2677: p={P} n={N} K2={K2} m={M_SIGS} "
          f"lam*={lam_star(LAM, N):.3f} nu_hat@K1=2={nu_hat(N, LAM, s[0], s[2]):.3f}")
    print(f"{len(SEEDS)} seeds per cell\n")
    print(f"{'K1':>4}{'eff':>8}{'gamma':>7}   "
          f"{'KAN':>7}{'PROJ':>7}{'BABAI':>7}{'CVP':>7}")
    for k1 in (4, 5, 6, 7, 8):
        agg, used, _ = sweep(curve, M_SIGS, k1, SEEDS)
        cells = "".join(f"{agg[k]:>4}/{used:<3}" for k in METHODS)
        print(f"{k1:>4}{k1*K2/N:>8.4f}{gh_margin(N, M_SIGS, k1, K2):>7.2f}   {cells}")
    print("\nDone.")
