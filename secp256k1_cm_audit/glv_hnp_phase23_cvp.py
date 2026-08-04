"""
GLV-HNP Thread 23: reformulate the Phase-2 lattice as an explicit CVP.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 Kannan lattice (dim 2m+2, `build_glv_lattice` in
  glv_hnp_phase2_20bit.py:263) always has lambda_1 = the trivial vector
  n*S_D*e_m -- 2-3x SHORTER than the planted vector, on every curve tested
  (|sv[m]|/n = 1.0000 exactly).  So Phase-2 recovery was never an SVP
  condition; it is a BDD/coset condition.  That run proposed:

    "Thread 23 -- reformulate the lattice so the target is lambda_1.
     Project along e_m, or replace the Kannan embedding with an explicit
     CVP call (Babai nearest-plane) targeting (A_i*S_K1, 0, ...).
     Falsifier: if the K1 wall on the lam*=0.07 curve (currently K1 ~ 4-6)
     moves outward, the reformulation is a real improvement; if the wall
     stays, the wall is information-theoretic and Phase 2 is at its ceiling."

This script executes exactly that, and sharpens the falsifier: instead of
only Babai (an approximation), it also runs *exact* CVP (fplll enumeration).
Exact CVP is the strongest possible attack in this formulation, so:

    exact CVP recovers d   -> any residual failure is reduction quality
    exact CVP fails        -> the planted vector is not the closest vector,
                              i.e. the wall is INFORMATION-THEORETIC and no
                              amount of lattice reduction can move it.

The CVP lattice L23 (dim 2m+1 -- the Kannan coordinate is gone):

  cols:  0..m-1  = k1 columns      (scale S_K1 = n // K1)
         m       = d column        (scale S_D  = 1)
         m+1..2m = k2 columns      (scale S_K2 = n // K2)

  rows:  i          :  n*S_K1  at col i                         (i < m)
         m          :  B_i*S_K1 at col i (all i),  S_D at col m
         m+1+i      : -lam*S_K1 at col i,  S_K2 at col m+1+i

  target t[i] = -A_i*S_K1,  t[m] = 0,  t[m+1+i] = 0
  planted lattice vector v* = t + (k1_i*S_K1, d, k2_i*S_K2),
  using  A_i + B_i*d = k1_i + lam*k2_i  (mod n).

  The trivial vector n*S_D*e_m is STILL in L23 (n*row_m - sum_i B_i*row_i),
  but it is now harmless: adding it maps d -> d+n, and d is only defined
  mod n.  That is the whole point of the reformulation.

  "Centered" variant: shift t by half the range of each unknown
  (t[i] += (K1//2)*S_K1, t[m] += n//2, t[m+1+i] += (K2//2)*S_K2), so the
  error vector is centered on 0 and ||e|| drops by ~2x.

Run: python3 glv_hnp_phase23_cvp.py
Requires: fpylll (+ cysignals), sympy.
"""

import math
import random
from fpylll import IntegerMatrix, LLL, GSO, CVP

# Lattice construction, EC arithmetic, signature generation and both attacks
# live in the shared library so that glv_hnp_phase23_cvp_scaling.py uses the
# identical code path.
from glv_hnp_phase23_cvp_lib import (
    find_generator, gen_signatures,
    build_glv_lattice, attack_kannan,
    build_cvp_lattice, planted_vector, babai_nearest_plane, attack_cvp,
    norm2, dist2, lam_star,
)

# ---------------------------------------------------------------------------
# Sweep driver
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

def trial_set(curve, k1_bound, m, seeds, method):
    """method in {'kannan','kannan_bkz','babai','babai_c','cvp','cvp_c'}"""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    wins = 0
    ratios = []      # sqrt(dist(v_cvp,t)^2 / dist(v*,t)^2); <1 means a
                     # non-planted vector is strictly closer than the truth
    planted_closest = 0
    total = 0
    for seed in seeds:
        d_secret = random.Random(seed ^ 0xABCDEF).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, b,
                              k1_bound, k2_bound, seed)
        if len(sigs) < m:
            continue
        total += 1
        if method == 'kannan':
            wins += bool(attack_kannan(sigs, n, lam, k1_bound, k2_bound,
                                       d_secret))
        elif method == 'kannan_bkz':
            wins += bool(attack_kannan(sigs, n, lam, k1_bound, k2_bound,
                                       d_secret, use_bkz=True, bkz_beta=30))
        else:
            centered = method.endswith('_c')
            exact = method.startswith('cvp')
            ok, dv, dp, is_pl = attack_cvp(sigs, n, lam, k1_bound, k2_bound,
                                           d_secret, centered=centered,
                                           exact=exact)
            wins += bool(ok)
            planted_closest += bool(is_pl)
            ratios.append(math.sqrt(dv / dp) if dp else float('nan'))
    mean_ratio = (sum(ratios) / len(ratios)) if ratios else float('nan')
    return wins, total, mean_ratio, planted_closest

# ---------------------------------------------------------------------------
# Historical curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
# ---------------------------------------------------------------------------

HIST = [
    # label,          p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",     211,  2, 199,  106,  2,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",   2677, 2, 2647, 185,  8,  10),
]

curves = []
for label, p, b, n, lam, k1h, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G), k1h, m))

print("=" * 78)
print("Thread 23 - Phase-2 lattice reformulated as explicit CVP")
print("=" * 78)
for label, (p, b, n, lam, G), k1h, m in curves:
    k2 = math.isqrt(n) + 1
    print(f"  {label:14s} p={p:5d} n={n:5d} lam={lam:5d} "
          f"lam*={lam_star(lam,n):.4f} K2={k2} m={m}")

# ===========================================================================
# EXP A - sanity: does the CVP formulation work at all?
# ===========================================================================

print("\n" + "-" * 78)
print("EXP A: CVP formulation sanity check (all three historical curves)")
print("-" * 78)
print("Methods: kannan = 2026-07-29 baseline (dim 2m+2, SVP/Kannan, LLL)")
print("         babai  = nearest-plane on L23 (dim 2m+1), _c = centered target")
print("         cvp    = exact CVP (fplll enumeration) on L23")
print()
hdr = (f"{'curve':14s} {'K1':>3s} {'m':>3s} | {'kannan':>7s} {'babai':>7s} "
       f"{'babai_c':>7s} {'cvp':>7s} {'cvp_c':>7s} | {'d_cvp/d_pl':>10s} "
       f"{'pl=closest':>10s}")
print(hdr)
print("-" * len(hdr))
for label, curve, k1h, m in curves:
    row = {}
    for meth in ('kannan', 'babai', 'babai_c', 'cvp', 'cvp_c'):
        row[meth] = trial_set(curve, k1h, m, SEEDS, meth)
    _, _, ratio_c, pl_c = row['cvp_c']
    tot = row['cvp_c'][1]
    print(f"{label:14s} {k1h:3d} {m:3d} | "
          + " ".join(f"{row[me][0]}/{row[me][1]:<5d}"[:7].rjust(7)
                     for me in ('kannan', 'babai', 'babai_c', 'cvp', 'cvp_c'))
          + f" | {ratio_c:10.4f} {pl_c}/{tot}".rjust(11))

# ===========================================================================
# EXP B - the falsifier: does the K1 wall move?
# ===========================================================================

print("\n" + "-" * 78)
print("EXP B: K1 sweep (the 2026-07-29 T4 grid), Kannan vs CVP")
print("-" * 78)
print("2026-07-29 T4 recorded, with Kannan+LLL:")
print("  12-bit/2557 (lam*=0.340): 5/5 up to K1=6, 4/5 at K1=12, 0/5 at K1=24")
print("  12-bit/2677 (lam*=0.070): 5/5 up to K1=4, 0/5 at K1>=8")
print("Falsifier: if the CVP wall sits further right, the reformulation is a")
print("real improvement; if it sits in the same place, the wall is a data wall.")
print()

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
wall_summary = {}
for label, curve, k1h, m in curves[1:]:      # the two 12-bit curves
    p, b, n, lam, G = curve
    print(f"\n  {label}  (lam*={lam_star(lam,n):.4f}, m={m}, n={n})")
    hdr2 = (f"    {'K1':>3s} {'eff':>7s} | {'kannan':>7s} {'bkz30':>7s} "
            f"{'babai_c':>7s} {'cvp_c':>7s} | {'d_cvp/d_pl':>10s} "
            f"{'pl=closest':>10s}")
    print(hdr2)
    print("    " + "-" * (len(hdr2) - 4))
    walls = {}
    for k1 in K1_GRID:
        k2 = math.isqrt(n) + 1
        eff = k1 * k2 / n
        res = {}
        for meth in ('kannan', 'kannan_bkz', 'babai_c', 'cvp_c'):
            res[meth] = trial_set(curve, k1, m, SEEDS, meth)
        ratio = res['cvp_c'][2]
        pl = res['cvp_c'][3]
        tot = res['cvp_c'][1]
        print(f"    {k1:3d} {eff:7.4f} | "
              + " ".join(f"{res[me][0]}/{res[me][1]}".rjust(7)
                         for me in ('kannan', 'kannan_bkz', 'babai_c', 'cvp_c'))
              + f" | {ratio:10.4f} {pl}/{tot}".rjust(11))
        for meth in ('kannan', 'babai_c', 'cvp_c'):
            w, t, _, _ = res[meth]
            if w == t and t > 0:
                walls[meth] = k1          # last fully-successful K1
    wall_summary[label] = walls
    print(f"    wall (largest K1 with {len(SEEDS)}/{len(SEEDS)}): "
          + ", ".join(f"{k}={v}" for k, v in walls.items()))

# ===========================================================================
# EXP C - is the wall information-theoretic?  m-sweep at the wall
# ===========================================================================

print("\n" + "-" * 78)
print("EXP C: does more data rescue CVP past the wall?  (T4b analog)")
print("-" * 78)
print("2026-07-29 T4b: at K1=8 on 12-bit/2677, Kannan+LLL scored 0,0,1,0,1 of 5")
print("for m = 8,12,16,24,32.  Re-run with exact CVP.")
print()
label, curve, _, _ = curves[2]      # 12-bit/2677
p, b, n, lam, G = curve
hdr3 = (f"    {'m':>3s} {'dim':>4s} | {'kannan':>7s} {'cvp_c':>7s} | "
        f"{'d_cvp/d_pl':>10s} {'pl=closest':>10s}")
print(f"  {label}  K1=8  (eff={8*(math.isqrt(n)+1)/n:.4f})")
print(hdr3)
print("    " + "-" * (len(hdr3) - 4))
for m in (8, 12, 16, 24, 32):
    res = {}
    for meth in ('kannan', 'cvp_c'):
        res[meth] = trial_set(curve, 8, m, SEEDS, meth)
    ratio = res['cvp_c'][2]
    pl = res['cvp_c'][3]
    tot = res['cvp_c'][1]
    print(f"    {m:3d} {2*m+1:4d} | "
          + " ".join(f"{res[me][0]}/{res[me][1]}".rjust(7)
                     for me in ('kannan', 'cvp_c'))
          + f" | {ratio:10.4f} {pl}/{tot}".rjust(11))

# ===========================================================================
# EXP D - the structural claim: lambda_1 vs the planted vector, before/after
# ===========================================================================

print("\n" + "-" * 78)
print("EXP D: is the planted vector now the target?  (T5 re-measured on L23)")
print("-" * 78)
print("On the Kannan lattice, 2026-07-29 T5 found sv = n*S_D*e_m always, with")
print("sv/pv in [0.34, 0.61].  On L23 the same trivial vector still exists, but")
print("the question is now whether v* is the CLOSEST vector to t, not the")
print("shortest vector.  Reported: ||v_cvp - t|| / ||v* - t||.")
print()
hdr4 = (f"{'curve':14s} {'K1':>3s} | {'||v*-t||/n':>11s} {'lam1(L23)/n':>12s} "
        f"{'sv/pv':>7s} | {'d_cvp/d_pl':>10s} {'pl=closest':>11s} {'recov':>7s}")
print(hdr4)
print("-" * len(hdr4))
for label, curve, k1h, m in curves:
    p, b, n, lam, G = curve
    for k1 in sorted({k1h, 2, 8}):
        k2 = math.isqrt(n) + 1
        seed = SEEDS[0]
        d_secret = random.Random(seed ^ 0xABCDEF).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1, k2, seed)
        if len(sigs) < m:
            continue
        B, t, scales, dim = build_cvp_lattice(sigs, n, lam, k1, k2,
                                              centered=True)
        A = IntegerMatrix.from_matrix(B)
        LLL.reduction(A)
        sv = min((norm2([A[i][j] for j in range(dim)]) for i in range(dim)))
        vp = planted_vector(sigs, n, lam, k1, k2, d_secret)
        dp = dist2(vp, t)
        ok, dv, dp2, is_pl = attack_cvp(sigs, n, lam, k1, k2, d_secret,
                                        centered=True, exact=True)
        print(f"{label:14s} {k1:3d} | {math.sqrt(dp)/n:11.4f} "
              f"{math.sqrt(sv)/n:12.4f} {math.sqrt(sv/dp):7.4f} | "
              f"{math.sqrt(dv/dp):10.4f} {str(is_pl):>11s} {str(bool(ok)):>7s}")

print("\nDone.")
