"""
GLV-HNP Phase 2, Thread 23 follow-up (E4): what blocks the projected lattice?

E1-E3 (glv_hnp_phase2_projected.py) established:
  * the d-free lattice L' (dim 2m+1) is correct: det == closed form, planted
    vector is a member, and the trivial vector n*S_D*e_m is GONE;
  * sv/pv rose (0.417 -> 0.813 on the lam*=0.07 curve) but stayed BELOW 1, so
    the planted vector is still not lambda_1;
  * the K1 success grid is bit-identical between OLD and NEW (both curves,
    all 9 K1 values, 5 seeds).  The pre-registered falsifier therefore says
    the wall is NOT caused by the trivial vector.
  * BKZ(20) on L' rescued the lam*=0.07 curve at K1=4 (4/5 -> 5/5) and
    K1=6 (3/5 -> 5/5).  That gain needs checking against BKZ on L.

E4 asks the three questions those results raise.

  E4a  IDENTITY OF THE NEW BLOCKER.  Hypothesis H23a: the shortest vector of
       L' is the 2-D "lambda-block" vector, i.e. the Gauss-reduced shortest
       vector mu of
             B = < (n*S_K1, 0),  (-lam*S_K1, S_K2) >
       embedded in the (col j, col m+j) plane for some j >= 1.  (j = 0 is
       excluded: column 0 of the HNF basis carries the unit row (1,c_1,...),
       not n*S_K1*e_0.)  This is the same mu whose RATIO to ||v_planted|| was
       falsified as a *predictor* by Thread 20/T2 -- being the identity of the
       shortest vector and being a predictor of success are different claims,
       and T2 only tested the second.
       Falsifier: if min_j |mu| != observed sv, or the observed sv has support
       outside a single (j, m+j) plane, H23a is wrong.

  E4b  IS THE BKZ GAIN REAL?  Head-to-head LLL / BKZ-20 / BKZ-40 on OLD and
       NEW, 10 seeds (not 5), on the K1 values where E3 saw a difference.
       Falsifier: if NEW-BKZ <= OLD-BKZ at every K1, the reformulation buys
       nothing at all and Thread 23 is a pure negative result.

  E4c  DOES MORE DATA RESCUE L'?  T4b showed that in the OLD lattice, m =
       8/12/16/24/32 at K1=8 gave 0,0,1,0,1 of 5 -- more signatures do not
       help.  Re-test in L'.  If L' also flatlines, the K1 wall is
       information-theoretic in both formulations.

Run: python3 glv_hnp_phase2_projected_e4.py
Deps: fpylll, sympy
"""

import math

from glv_hnp_phase2_projected import (
    scales, gen_signatures, find_generator, lam_star, norm,
    build_old, build_new, planted_old, planted_new,
    recover_old, recover_new, rate, run_once,
    HIST, SEEDS, D_SECRET_FRAC,
)
from fpylll import IntegerMatrix, LLL, BKZ


def gauss_reduce_2d(u, v):
    """Lagrange-Gauss reduction of a rank-2 integer lattice; returns the
    exact shortest nonzero vector."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        # v <- v - round(<u,v>/<u,u>) * u
        q = (u[0] * v[0] + u[1] * v[1]) / n2(u)
        q = int(math.floor(q + 0.5))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if n2(v) >= n2(u):
            return u
        u, v = v, u


def lambda_block_mu(n, lam, S_K1, S_K2):
    return gauss_reduce_2d((n * S_K1, 0), (-lam * S_K1, S_K2))


SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 20260804, 65537, 123456789, 314159]

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    hist.append((label, (p, b, n, lam, G), k1, m))

CURVES12 = [(lb, cv, m) for (lb, cv, k1, m) in hist if '12-bit' in lb]

print("=" * 78)
print("Thread 23 / E4 - what blocks the projected lattice L'?")
print("=" * 78)

# ===========================================================================
# E4a - identity of the shortest vector in L'
# ===========================================================================
print()
print("=" * 78)
print("E4a  H23a: is sv(L') the 2-D lambda-block vector mu?")
print("=" * 78)
SVLP = "sv(L')"
print(f"{'curve':<16} {'K1':>3} {'sd':>3} | {'|mu|':>12} {SVLP:>12} "
      f"{'sv/|mu|':>8} | {'supp k1':>7} {'supp k2':>7} {'kan':>4} {'planes':>6}")
print("-" * 78)

a_rows = []
for label, curve, m in CURVES12:
    p, b, n, lam, G = curve
    d_secret = max(1, int(D_SECRET_FRAC * n))
    k2_bound = math.isqrt(n) + 1
    for k1_bound in (4, 8, 16):
        S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
        mu = lambda_block_mu(n, lam, S_K1, S_K2)
        mu_norm = math.sqrt(mu[0] ** 2 + mu[1] ** 2)
        for seed in SEEDS[:3]:
            r = run_once(curve, m, d_secret, k1_bound, seed, 'new')
            rows = [row for row in r['rows'] if any(row)]
            sv_row = min(rows, key=norm)
            sv = norm(sv_row)
            k1_supp = sum(1 for i in range(m) if sv_row[i])
            k2_supp = sum(1 for i in range(m) if sv_row[m + i])
            kan = sv_row[2 * m]
            # how many distinct (j, m+j) planes does the support touch?
            planes = sorted({i for i in range(m)
                             if sv_row[i] or sv_row[m + i]})
            single = (len(planes) == 1 and kan == 0)
            print(f"{label:<16} {k1_bound:>3} {seed % 1000:>3} | {mu_norm:>12.1f} "
                  f"{sv:>12.1f} {sv/mu_norm:>8.4f} | {k1_supp:>7} {k2_supp:>7} "
                  f"{kan:>4} {len(planes):>6}{'  <- single plane' if single else ''}")
            a_rows.append((label, k1_bound, seed, mu_norm, sv, len(planes), kan))
    print("-" * 78)

exact = sum(1 for r in a_rows if abs(r[4] - r[3]) < 1e-6)
print(f"H23a exact matches (sv == |mu|): {exact}/{len(a_rows)}")
print(f"sv/|mu| range: [{min(r[4]/r[3] for r in a_rows):.4f}, "
      f"{max(r[4]/r[3] for r in a_rows):.4f}]")
single_plane = sum(1 for r in a_rows if r[5] == 1 and r[6] == 0)
print(f"single-plane support: {single_plane}/{len(a_rows)}")

# ===========================================================================
# E4b - is the BKZ gain real?
# ===========================================================================
print()
print("=" * 78)
print("E4b  LLL / BKZ-20 / BKZ-40, OLD vs NEW, 10 seeds")
print("=" * 78)
print("(cell = #seeds recovering d, out of 10)")
print()
print(f"{'curve':<16} {'K1':>3} | {'old LLL':>7} {'old B20':>7} {'old B40':>7} | "
      f"{'new LLL':>7} {'new B20':>7} {'new B40':>7}")
print("-" * 78)

b_rows = []
for label, curve, m in CURVES12:
    p, b, n, lam, G = curve
    d_secret = max(1, int(D_SECRET_FRAC * n))
    for k1_bound in (4, 6, 8, 12):
        cells = []
        for variant in ('old', 'new'):
            for bkz, beta in ((False, 0), (True, 20), (True, 40)):
                g, _ = rate(curve, m, d_secret, k1_bound, SEEDS10, variant,
                            use_bkz=bkz, beta=beta)
                cells.append(g)
        print(f"{label:<16} {k1_bound:>3} | {cells[0]:>7} {cells[1]:>7} "
              f"{cells[2]:>7} | {cells[3]:>7} {cells[4]:>7} {cells[5]:>7}")
        b_rows.append((label, k1_bound, cells))
    print("-" * 78)

gain = [(lb, k, c) for (lb, k, c) in b_rows if c[4] > c[1] or c[5] > c[2]]
print(f"cells where NEW beats OLD at equal reduction strength: {len(gain)}/{len(b_rows)}")
for lb, k, c in gain:
    print(f"   {lb} K1={k}: B20 {c[1]}->{c[4]},  B40 {c[2]}->{c[5]}")
loss = [(lb, k, c) for (lb, k, c) in b_rows if c[4] < c[1] or c[5] < c[2]]
print(f"cells where NEW loses to OLD: {len(loss)}/{len(b_rows)}")
for lb, k, c in loss:
    print(f"   {lb} K1={k}: B20 {c[1]}->{c[4]},  B40 {c[2]}->{c[5]}")

# ===========================================================================
# E4c - does more data rescue L' at the wall?
# ===========================================================================
print()
print("=" * 78)
print("E4c  more signatures at the wall (K1=8), OLD vs NEW, 5 seeds")
print("=" * 78)
print("T4b (old lattice, 2026-07-29): m = 8/12/16/24/32 -> 0,0,1,0,1 of 5")
print()
label, curve, _m0 = CURVES12[1]      # the lam*=0.07 failure curve
p, b, n, lam, G = curve
d_secret = max(1, int(D_SECRET_FRAC * n))
print(f"  curve {label}  n={n}  lam*={lam_star(lam, n):.4f}  K1=8")
print(f"  {'m':>4} | {'old LLL':>7} {'new LLL':>7} {'new B20':>7} {'new B40':>7}")
print("  " + "-" * 42)
for m in (8, 12, 16, 24, 32):
    go, _ = rate(curve, m, d_secret, 8, SEEDS, 'old')
    gn, _ = rate(curve, m, d_secret, 8, SEEDS, 'new')
    g20, _ = rate(curve, m, d_secret, 8, SEEDS, 'new', use_bkz=True, beta=20)
    g40, _ = rate(curve, m, d_secret, 8, SEEDS, 'new', use_bkz=True, beta=40)
    print(f"  {m:>4} | {go:>7} {gn:>7} {g20:>7} {g40:>7}")

print()
print("=" * 78)
print("done")
print("=" * 78)
