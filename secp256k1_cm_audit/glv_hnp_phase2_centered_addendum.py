"""
Thread 23 addendum — two checks on the main run.

A. WHY is the e_m projection a bit-exact no-op?  Hypothesis: LLL immediately
   isolates the trivial vector n*S_D*e_m (axis-aligned, norm n, shorter than
   everything else) and reduces the rest inside its orthogonal complement,
   which IS the projected lattice.  Test: reduce the unprojected basis, delete
   column m from the result, and compare the row-norm multiset against the
   independently-reduced projected basis.

B. Curve-label check.  The 2026-07-29 T4 table calls its two curves
   "12-bit/2557" and "12-bit/2677", and the 2026-06-15 session says
   "p=2677, lam=185".  E3 of the main run resolved 2557/2677 as GROUP ORDERS n.
   Print both readings so future runs stop conflating p and n.

Run: python3 glv_hnp_phase2_centered_addendum.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

from glv_hnp_phase2_centered import (
    build_lattice, gen_signatures, scales, collect_curves, eisenstein_decompose,
    j0_traces, glv_roots, lam_star, build_curve, predict,
)

print("=" * 78)
print("A — is the e_m projection the same lattice LLL already works in?")
print("=" * 78)

CURVES = collect_curves(2 ** 16, 2 ** 17, 3)
M_SIGS = 12

print(f"{'n':>7} {'eff':>6} {'cent':>5} {'trivial row present':>20} "
      f"{'norm multisets equal':>22}")
for curve in CURVES:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    for eff_t in (0.10, 0.30):
        k1 = max(2, int(eff_t * n / k2))
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
        d = random.Random(7).randint(1, n - 1)
        sigs = gen_signatures(G, d, M_SIGS, n, lam, p, k1, k2, 42)
        for cen in (False, True):
            Mu, _ = build_lattice(sigs, n, lam, k1, k2, cen, False)
            Mp, _ = build_lattice(sigs, n, lam, k1, k2, cen, True)
            Au = IntegerMatrix.from_matrix(Mu); LLL.reduction(Au)
            Ap = IntegerMatrix.from_matrix(Mp); LLL.reduction(Ap)
            ru = [[Au[i][j] for j in range(Au.ncols)] for i in range(Au.nrows)]
            rp = [[Ap[i][j] for j in range(Ap.ncols)] for i in range(Ap.nrows)]

            m = M_SIGS
            triv = [r for r in ru
                    if abs(r[m]) == n * S_D and not any(r[j] for j in range(len(r))
                                                        if j != m)]
            # delete column m from the unprojected reduced basis
            ru_del = [r[:m] + r[m + 1:] for r in ru if r not in triv]
            nu = sorted(sum(x * x for x in r) for r in ru_del)
            np_ = sorted(sum(x * x for x in r) for r in rp if any(r))
            print(f"{n:>7} {eff_t:>6.2f} {str(cen):>5} "
                  f"{('yes' if triv else 'NO'):>20} "
                  f"{str(nu == np_):>22}")

print("\n" + "=" * 78)
print("B — p vs n for the two '12-bit' curves of the 2026-07-29 T4 table")
print("=" * 78)
for val in (2557, 2677):
    print(f"\nreading '{val}' as p:")
    if sympy.isprime(val) and val % 3 == 1:
        eis = eisenstein_decompose(val)
        if eis:
            for t in sorted(set(j0_traces(*eis))):
                nc = val + 1 - t
                if nc > 1 and sympy.isprime(nc) and nc % 3 == 1:
                    r = glv_roots(nc)
                    if r:
                        print(f"   p={val} -> n={nc:>6} lam={r[0]:>6} "
                              f"lam*={lam_star(r[0], nc):.4f}")
    else:
        print(f"   p={val}: prime={sympy.isprime(val)}, "
              f"p mod 3={val % 3} (need 1)")
    print(f"reading '{val}' as n:")
    pp = int(sympy.nextprime(2000))
    hits = []
    while pp < 4000:
        if pp % 3 == 1:
            eis = eisenstein_decompose(pp)
            if eis and val in [pp + 1 - t for t in j0_traces(*eis)]:
                hits.append(pp)
        pp = int(sympy.nextprime(pp))
    for pp in hits:
        r = glv_roots(val)
        print(f"   p={pp:>6} -> n={val} lam={r[0]:>6} "
              f"lam*={lam_star(r[0], val):.4f}")

print("\n" + "=" * 78)
print("C — asymptote vs finite-m: gap(m) for base and cent at fixed eff")
print("=" * 78)
n = CURVES[0][2]
k2 = math.isqrt(n) + 1
print(f"n={n}  (gap -> sqrt(a)/sqrt(2*pi*e*eff) as m -> infinity)")
print(f"{'eff':>6} " + " ".join(f"m={m:<5}" for m in (8, 12, 20, 32, 64, 128,
                                                      512)) + "   asympt")
for eff_t in (0.25, 0.55):
    k1 = max(2, int(eff_t * n / k2))
    eff = k1 * k2 / n
    for cen, a in ((False, 3.0), (True, 12.0)):
        gs = [predict(m, n, k1, k2, cen, False)[0]
              for m in (8, 12, 20, 32, 64, 128, 512)]
        asym = math.sqrt(a) / math.sqrt(2 * math.pi * math.e * eff)
        tag = "cent" if cen else "base"
        print(f"{eff:>6.3f} " + " ".join(f"{g:<7.3f}" for g in gs) +
              f"  {asym:>6.3f}  {tag}")

print("\ndone.")
