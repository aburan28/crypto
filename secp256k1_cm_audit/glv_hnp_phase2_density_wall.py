"""
GLV-HNP Phase 2, Thread 23b: the K1 wall is a lattice-DENSITY wall at
eff = 3/(2*pi*e) ~ 0.1757, independent of m and of n.

Thread 23a (glv_hnp_phase2_projected.py) removed the trivial n*S_D*e_m vector
by quotienting the lattice along e_m and found the success rate unchanged in
every single trial (eff=0.05: 96/100 vs 96/100; 0.15: 18/100 vs 18/100;
0.25: 9/100 vs 9/100).  Masking by the trivial vector was therefore not the
obstruction.  This script identifies what is.

Closed form
-----------
Write eff = K1*K2/n.  With S_K1 = n/K1, S_K2 = n/K2, S_D = 1, S_KANNAN = n:

    det(L)   = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN = n^(2m+1) / eff^m
    dim(L)   = 2m+2
    E||pv||^2 = m*(K1*S_K1)^2/3 + n^2/3 + m*(K2*S_K2)^2/3 + n^2
              = n^2 * (2m/3 + 4/3)

    det(L_proj) = det(L) / (n*S_D) = n^(2m) / eff^m
    dim(L_proj) = 2m+1
    E||pv_proj||^2 = n^2 * (2m/3 + 1)

Gaussian heuristic  GH = sqrt(D/(2*pi*e)) * det^(1/D).  For either lattice
det^(1/D) -> n / sqrt(eff) as m -> infinity, so

    ||pv|| / GH  ->  sqrt(2m/3) * sqrt(2*pi*e / (2m)) * sqrt(eff)
                 =  sqrt(2*pi*e*eff/3)
                 =  2.3861 * sqrt(eff)

-- independent of m and of n.  The planted vector is the unique shortest
vector (so unique-SVP / BDD is well posed) exactly when

    eff  <  3 / (2*pi*e)  =  0.175656...

Above that, the ball of radius ||pv|| already contains many lattice vectors
by the Gaussian heuristic, so no amount of reduction strength and no
reformulation of the same lattice can single the planted vector out.  Note
this is NOT an information-theoretic bound: the solution is still uniquely
determined whenever m*log2(1/eff) > log2(n), which holds far above 0.176.
It is a bound on the lattice formulation.

Predictions tested here
-----------------------
  E1  numerical: ||pv||/GH converges to 2.3861*sqrt(eff) for L and L_proj.
  E2  the empirical success curve crosses 50% near eff = 0.176.
  E3  m-independence: at fixed eff, adding signatures does not move the wall
      (retro-explains T4b of 2026-07-29: m = 8/12/16/24/32 at K1=8 gave
      0,0,1,0,1 of 5).
  E4  n-independence: same crossing at 14, 17 and 20 bits.

Run: python3 glv_hnp_phase2_density_wall.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL

from glv_hnp_phase2_projected import (  # noqa: E402
    ec_mul, tonelli_shanks, eisenstein_decompose, j0_traces, glv_roots,
    lam_star, scales, gen_signatures, build_glv_lattice,
    build_projected_lattice, planted_vector, planted_vector_proj, norm,
    recover_d_kannan, recover_d_proj, build_curve, search_curves, rate,
)

TWO_PI_E = 2 * math.pi * math.e
EFF_STAR = 3.0 / TWO_PI_E
CONST = math.sqrt(TWO_PI_E / 3.0)

print("=" * 78)
print("Thread 23b -- the Phase-2 K1 wall is a lattice-density wall")
print("=" * 78)
print(f"predicted threshold  eff* = 3/(2*pi*e) = {EFF_STAR:.6f}")
print(f"predicted ratio      ||pv||/GH -> {CONST:.4f} * sqrt(eff)")

# ---------------------------------------------------------------------------
# E1 -- numerical check of the closed form
# ---------------------------------------------------------------------------

def gh(dim, log_det):
    """Gaussian heuristic length from log(det)."""
    return math.sqrt(dim / TWO_PI_E) * math.exp(log_det / dim)

def geometry(n, k1_bound, k2_bound, m, projected):
    """Exact log-det and expected planted norm for the actual integer scales."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if projected:
        dim = 2 * m + 1
        log_det = (m * math.log(n * S_K1) + m * math.log(S_K2)
                   + math.log(S_KAN) - math.log(n * S_D))
        pv2 = (m * (k1_bound * S_K1) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    else:
        dim = 2 * m + 2
        log_det = (m * math.log(n * S_K1) + math.log(S_D)
                   + m * math.log(S_K2) + math.log(S_KAN))
        pv2 = (m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0 + S_KAN ** 2)
    return dim, log_det, math.sqrt(pv2), gh(dim, log_det)

print("\n" + "-" * 78)
print("E1: ||pv||/GH vs the closed form  (n = 17-bit, K2 = isqrt(n)+1)")
print("-" * 78)
n_demo = 131101 if sympy.isprime(131101) else int(sympy.nextprime(2 ** 17))
k2_demo = math.isqrt(n_demo) + 1
print(f"n = {n_demo} ({n_demo.bit_length()} bits), K2 = {k2_demo}")
print(f"\n{'eff':>6}{'pred':>9} | " + "".join(f"{'m=' + str(m):>16}" for m in (6, 12, 24, 48)))
print(f"{'':>6}{'':>9} | " + "".join(f"{'L':>8}{'L_proj':>8}" for _ in range(4)))
for eff in (0.05, 0.10, 0.15, EFF_STAR, 0.20, 0.25, 0.40):
    k1 = max(2, round(eff * n_demo / k2_demo))
    eff_actual = k1 * k2_demo / n_demo
    cells = []
    for m in (6, 12, 24, 48):
        for proj in (False, True):
            _d, ld, pv, g = geometry(n_demo, k1, k2_demo, m, proj)
            cells.append(f"{pv / g:>8.3f}")
    print(f"{eff_actual:>6.3f}{CONST * math.sqrt(eff_actual):>9.3f} | " + "".join(cells))
print("\n-> both lattices converge to the same limit; the projection cannot help.")
sys.stdout.flush()

# ---------------------------------------------------------------------------
# E2 -- empirical crossing
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E2: empirical success curve vs eff  (17-bit, 20 curves, m=12, 5 seeds)")
print("-" * 78)
sys.stdout.flush()

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 12
curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
print(f"{len(curves17)} curves, lam* in "
      f"[{min(lam_star(c[3], c[2]) for c in curves17):.4f}, "
      f"{max(lam_star(c[3], c[2]) for c in curves17):.4f}]\n")
sys.stdout.flush()

EFF_GRID = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.24, 0.30]
print(f"{'eff':>6}{'pv/GH':>8}{'trials':>9}{'rate':>8}{'curves 5/5':>12}")
e2 = []
for eff in EFF_GRID:
    wins = tot = full = 0
    ratios = []
    for cur in curves17:
        p, b, n, lam, G = cur
        k2 = math.isqrt(n) + 1
        k1 = max(2, round(eff * n / k2))
        _d, ld, pv, g = geometry(n, k1, k2, M_SIGS, False)
        ratios.append(pv / g)
        w, t, _ = rate(cur, M_SIGS, k1, SEEDS, 'kannan')
        wins += w; tot += t; full += (w == t)
    r = wins / tot
    e2.append((eff, sum(ratios) / len(ratios), r))
    print(f"{eff:>6.2f}{sum(ratios)/len(ratios):>8.3f}{tot:>9}{r:>8.2f}"
          f"{full:>8}/{len(curves17)}")
    sys.stdout.flush()

# locate the 50% crossing by linear interpolation on the measured curve
cross = None
for i in range(1, len(e2)):
    if e2[i - 1][2] >= 0.5 > e2[i][2]:
        (x0, _, y0), (x1, _, y1) = e2[i - 1], e2[i]
        cross = x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
        break
print(f"\nempirical 50% crossing: eff = {cross:.4f}" if cross
      else "\nno 50% crossing in grid")
print(f"predicted:              eff* = {EFF_STAR:.4f}")

# ---------------------------------------------------------------------------
# E3 -- m-independence
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E3: m-independence of the wall  (17-bit, 10 curves, 3 seeds)")
print("-" * 78)
print("if the wall were a data-shortage, success would rise with m.")
sys.stdout.flush()

M_GRID = [6, 8, 12, 16, 24, 32]
print(f"\n{'eff':>6} | " + "".join(f"{'m=' + str(m):>8}" for m in M_GRID))
for eff in (0.08, 0.14, 0.22):
    cells = []
    for m in M_GRID:
        wins = tot = 0
        for cur in curves17[:10]:
            p, b, n, lam, G = cur
            k2 = math.isqrt(n) + 1
            k1 = max(2, round(eff * n / k2))
            w, t, _ = rate(cur, m, k1, SEEDS[:3], 'kannan')
            wins += w; tot += t
        cells.append(f"{wins/tot:>8.2f}")
    print(f"{eff:>6.2f} | " + "".join(cells))
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# E4 -- n-independence
# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E4: n-independence of the wall  (m=12, 3 seeds, 8 curves per size)")
print("-" * 78)
sys.stdout.flush()

print(f"\n{'bits':>5} | " + "".join(f"{e:>8.2f}" for e in (0.08, 0.14, 0.22)))
for bits in (14, 17, 20):
    cs = search_curves(2 ** (bits - 1), 2 ** bits, per_bin=1, nbins=8)[:8]
    cells = []
    for eff in (0.08, 0.14, 0.22):
        wins = tot = 0
        for cur in cs:
            p, b, n, lam, G = cur
            k2 = math.isqrt(n) + 1
            k1 = max(2, round(eff * n / k2))
            w, t, _ = rate(cur, M_SIGS, k1, SEEDS[:3], 'kannan')
            wins += w; tot += t
        cells.append(f"{wins/tot:>8.2f}")
    print(f"{bits:>5} | " + "".join(cells) + f"   ({len(cs)} curves)")
    sys.stdout.flush()

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("Information-theoretic margin (for contrast -- NOT the binding constraint)")
print("-" * 78)
print("unknowns m*log2(K1*K2) + log2(n) bits vs constraints m*log2(n) bits;")
print("solution unique when m*log2(1/eff) > log2(n).\n")
print(f"{'eff':>6}{'m*log2(1/eff)':>16}{'log2(n), 17b':>14}{'unique?':>9}")
for eff in (0.08, 0.14, 0.22, 0.40, 0.70):
    lhs = M_SIGS * math.log2(1 / eff)
    print(f"{eff:>6.2f}{lhs:>16.1f}{17.0:>14.1f}{'yes' if lhs > 17 else 'no':>9}")
print("\n-> at m=12 the solution stays information-theoretically unique up to")
print("   eff ~ 0.37; the lattice wall at 0.176 is strictly tighter.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
