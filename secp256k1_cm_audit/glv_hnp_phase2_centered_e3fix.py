"""
Thread 23, corrected E3 + the high-m test.

The addendum's check B showed the 2026-07-29 T4 table labels its curves by p,
not n: "12-bit/2557" is (p=2557, n=2659, lam=903, lam*=0.3396) and
"12-bit/2677" is (p=2677, n=2647, lam=185, lam*=0.0699).  E3 in the main run
resolved those numbers as group orders and therefore swept two DIFFERENT
curves.  This re-runs the K1 sweep on the curves the log actually used —
in particular the lam*=0.0699 curve named in the Thread-23 falsifier.

E6 then tests the model's sharpest prediction: gap(m) climbs to
sqrt(a)/sqrt(2*pi*e*eff) from below, so the centred lattice should cross
gap = 1 at eff = 0.55 somewhere near m = 64, an operating point where the
uncentred lattice is flatly dead at every m.

Run: python3 glv_hnp_phase2_centered_e3fix.py
"""

import math

from glv_hnp_phase2_centered import (
    VARIANTS, SEEDS, build_curve, glv_roots, lam_star, predict, rate,
    collect_curves,
)

M_SIGS = 12
K1S = (2, 3, 4, 6, 8, 12, 16, 24, 32, 48)

# (p, n) exactly as used by the 2026-07-29 T4 table.
T4_CURVES = [(2557, 2659), (2677, 2647)]

print("=" * 78)
print("E3-fixed — K1 sweep on the ACTUAL 2026-07-29 T4 curves")
print("=" * 78)

for p, n in T4_CURVES:
    cur = build_curve(p, n)
    lam = glv_roots(n)[0]
    curve = (p, cur[1], n, lam, cur[4])
    k2 = math.isqrt(n) + 1
    print(f"\ncurve p={p} n={n} lam={lam} lam*={lam_star(lam, n):.4f} K2={k2}")
    print(f"{'K1':>4} {'eff':>6} " + " ".join(f"{v[0]:>11}" for v in VARIANTS))
    for k1 in K1S:
        cells = []
        for name, cen, proj in VARIANTS:
            w, t = rate(curve, M_SIGS, k1, SEEDS, cen, proj)
            g = predict(M_SIGS, n, k1, k2, cen, proj)[0]
            cells.append(f"{w}/{t} {g:>5.2f}")
        print(f"{k1:>4} {k1*k2/n:>6.3f} " + " ".join(f"{c:>11}" for c in cells))

print("\n" + "=" * 78)
print("E6 — high-m test: does the centred lattice cross gap=1 at eff=0.55?")
print("=" * 78)
CURVES = collect_curves(2 ** 16, 2 ** 17, 4)
print(f"{'eff':>6} {'m':>4} {'base':>12} {'cent':>12}   (wins over "
      f"{len(CURVES)} curves x {len(SEEDS)} seeds, gap in parens)")
for eff_t in (0.40, 0.55):
    for m in (32, 48, 64, 96):
        cells = []
        for cen in (False, True):
            w = t = 0
            gs = []
            for curve in CURVES:
                n = curve[2]
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff_t * n / k2))
                gs.append(predict(m, n, k1, k2, cen, False)[0])
                ww, tt = rate(curve, m, k1, SEEDS, cen, False)
                w += ww; t += tt
            cells.append(f"{w:>2}/{t:<3}({sum(gs)/len(gs):.2f})")
        print(f"{eff_t:>6.2f} {m:>4} {cells[0]:>12} {cells[1]:>12}")

print("\ndone.")
