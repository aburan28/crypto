"""
Thread 23 (part 2) — the K1 wall is a BDD condition in the Kannan-free sublattice.

Part 1 (glv_hnp_phase2_projected.py) showed:
  * projecting out the trivial vector n*S_D*e_m does NOT move the K1 wall
    (identical grids for 12-bit/2557; marginally worse for 12-bit/2677),
  * so "the planted vector is not lambda_1" (2026-07-29 T5) was a true
    observation but a false diagnosis -- recovery never required lambda_1,
    only that the planted vector appear among the LLL basis rows,
  * and the Gaussian-heuristic criterion tau < 1 under-predicts the wall by a
    factor 2-4 in K1, with no common tau value at the two observed walls
    (tau = 2.05 at the 2557 wall vs 1.35 at the 2677 wall).

The correct frame is BDD, not SVP.  Write

    L'   = projected lattice, rank 2m+1
    L'_0 = { v in L' : last coordinate 0 }, rank 2m     (the "noise" sublattice;
           basis = all generators EXCEPT the Kannan row)

Recovery asks for the lattice point of L'_0 nearest the target, with error
exactly v_planted.  Unique decoding holds when

    gamma = ||v_planted|| / lambda_1(L'_0)  <  1/2.

This script measures gamma directly (lambda_1(L'_0) via BKZ) and tests whether a
single gamma threshold predicts the wall on BOTH curves, controlling for m --
which the 2026-07-29 T4 table did not (2557 was run at m=8, 2677 at m=10).

  E5  gamma vs observed success, both curves, K1 grid
  E6  m-controlled K1 wall: m in {8,10,12,16}, both curves
  E7  does BKZ-40 move the wall in the projected lattice?

Run: python3 glv_hnp_phase2_bdd.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

from glv_hnp_phase2_projected import (
    find_generator, lam_star, scales, gen_signatures,
    build_projected_lattice, planted_vector_projected, recover_d_projected,
    run_projected, norm, rate, gh_tau,
)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755, 8),
    ("12-bit/2677", 2677, 2, 2647, 185,  10),
]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

curves = []
for label, p, b, n, lam, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0
    curves.append((label, (p, b, n, lam, G), m))


def lambda1_noise(curve, m, d_secret, k1_bound, seed, beta=30):
    """
    lambda_1(L'_0) where L'_0 = projected lattice minus the Kannan row,
    together with ||v_planted|| for the same signature set.
    """
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m:
        return None
    rows = build_projected_lattice(sigs, n, lam, k1_bound, k2b)
    noise = rows[:-1]                      # drop the Kannan / target row
    cols = 2 * m + 1
    A = IntegerMatrix.from_matrix(noise)
    BKZ.reduction(A, BKZ.Param(beta))
    red = [[A[i][j] for j in range(cols)] for i in range(A.nrows)]
    mu = min((norm(r) for r in red if any(r)), default=float('nan'))
    pn = norm(planted_vector_projected(sigs, n, k1_bound, k2b))
    return mu, pn


def mean_gamma(curve, m, k1_bound, seeds):
    gs = []
    for s in seeds:
        d = random.Random(s + 7777).randint(1, curve[2] - 1)
        r = lambda1_noise(curve, m, d, k1_bound, s)
        if r is None:
            continue
        mu, pn = r
        gs.append(pn / mu if mu else float('nan'))
    return sum(gs) / len(gs) if gs else float('nan')


print("=" * 78)
print("Thread 23 part 2 — BDD criterion for the GLV-HNP Phase-2 K1 wall")
print("=" * 78)
for label, cur, m in curves:
    print(f"  {label}: n={cur[2]} lam={cur[3]} lam*={lam_star(cur[3], cur[2]):.4f} m={m}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("E5 — gamma = ||v_planted|| / lambda_1(L'_0)   vs observed success")
print("     (unique-decoding BDD predicts success when gamma < 0.5)")
print("=" * 78)
e5 = {}
for label, cur, m in curves:
    gammas, wins = [], []
    for k1 in K1_GRID:
        gammas.append(mean_gamma(cur, m, k1, SEEDS))
        w, t, _ = rate(run_projected, cur, m, k1, SEEDS)
        wins.append(w)
    e5[label] = (gammas, wins)
    print(f"\n{label}  m={m}")
    print("  K1      " + " ".join(f"{k:>7}" for k in K1_GRID))
    print("  gamma   " + " ".join(f"{g:>7.3f}" for g in gammas))
    print("  wins    " + " ".join(f"{str(w)+'/5':>7}" for w in wins))
    last_ok = max((g for g, w in zip(gammas, wins) if w == 5), default=float('nan'))
    first_bad = min((g for g, w in zip(gammas, wins) if w < 5), default=float('nan'))
    print(f"  gamma at last 5/5: {last_ok:.3f}   gamma at first failure: {first_bad:.3f}")

print("\n  --- cross-curve consistency of the gamma threshold ---")
for label in e5:
    g, w = e5[label]
    lo = max((x for x, y in zip(g, w) if y == 5), default=float('nan'))
    hi = min((x for x, y in zip(g, w) if y < 5), default=float('nan'))
    print(f"  {label}: threshold bracketed in ({lo:.3f}, {hi:.3f})")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("E6 — m-controlled K1 wall (projected lattice, LLL, 5 seeds)")
print("     the 2026-07-29 T4 table compared 2557@m=8 against 2677@m=10;")
print("     this removes that confound.")
print("=" * 78)
M_GRID = [8, 10, 12, 16]
for label, cur, _m in curves:
    print(f"\n{label}   lam*={lam_star(cur[3], cur[2]):.4f}")
    print("  m\\K1  " + " ".join(f"{k:>6}" for k in K1_GRID))
    for m in M_GRID:
        row = []
        for k1 in K1_GRID:
            w, t, _ = rate(run_projected, cur, m, k1, SEEDS)
            row.append(w)
        wall = next((k for k, w in zip(K1_GRID, row) if w < 5), None)
        print(f"  m={m:<4} " + " ".join(f"{str(w)+'/5':>6}" for w in row) +
              f"   wall K1={wall}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("E7 — does BKZ-40 move the wall in the projected lattice?")
print("=" * 78)
print(f"{'curve':<14} {'m':<4} {'K1':<5} {'LLL':<7} {'BKZ-40':<8}")
for label, cur, m in curves:
    for k1 in [4, 6, 8, 12, 16]:
        wl, t, _ = rate(run_projected, cur, m, k1, SEEDS)
        wb, t, _ = rate(run_projected, cur, m, k1, SEEDS, use_bkz=True, beta=40)
        print(f"{label:<14} {m:<4} {k1:<5} {str(wl)+'/5':<7} {str(wb)+'/5':<8}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
