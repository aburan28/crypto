"""
GLV-HNP Phase 2, Thread 23 follow-up (E5): is the centered eff-ceiling lambda*-independent?

glv_hnp_phase2_projected.py showed that CENTERING (not the d-column projection) is what
moves the K1 wall, and that it moves it by ~4x -- exactly the factor predicted by the
second-moment argument (E[x^2]=K^2/3 uncentered vs Var=K^2/12 centered, so ||pv||^2 drops
4x and the eff ceiling rises 4x):

    uncentered ceiling  eff < 3/(4*pi*e) = 0.0878
    centered   ceiling  eff < 3/(pi*e)   = 0.351

This script tests that against the recorded 2026-07-29 T3 sweep, which used 20 fresh
17-bit j=0 GLV curves at eff in {0.05, 0.15, 0.25} and found (UNCENTERED):

    eff=0.05 -> 19/20 curves recover 5/5
    eff=0.15 ->  3/20
    eff=0.25 ->  0/20

Prediction under the ceiling formula: centering should lift eff=0.15 and eff=0.25 to
near-total recovery (both are < 0.351), while eff=0.45 should start to fail. And the
lambda* spread should have no residual effect at all in the centered variant.

Run: python3 glv_hnp_phase2_eff_ceiling.py
"""

import math
import random
import sympy

# Reuse the exact constructions from the Thread 23 script (helpers only, no experiments).
_src = open('glv_hnp_phase2_projected.py').read()
exec(compile(_src[:_src.index('# Reference curves')], 'glv_hnp_phase2_projected.py', 'exec'))

print("=" * 78)
print("E5: is the centered eff-ceiling lambda*-independent?  (17-bit curves)")
print("=" * 78)

# ---------------------------------------------------------------------------
# Find 17-bit j=0 GLV curves with lambda* spread over (0, 0.5)  [T3 protocol]
# ---------------------------------------------------------------------------

def find_17bit_curves(target=12):
    """17-bit prime-order j=0 curves, chosen to spread lambda* across (0, 0.5)."""
    found = []
    p = sympy.nextprime(2 ** 16)
    while p < 2 ** 17 and len(found) < 400:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_c = p + 1 - t
                    if n_c < 2 or not sympy.isprime(n_c) or n_c % 3 != 1:
                        continue
                    lam, _ = glv_eigenvalue(n_c)
                    if lam is None:
                        continue
                    b_ok = None
                    for b_try in range(1, 60):
                        P = find_point(p, b_try)
                        if P is not None and ec_mul(P, n_c, p) is None:
                            b_ok = b_try
                            break
                    if b_ok is None:
                        continue
                    G = find_generator(p, b_ok, n_c)
                    if G is None:
                        continue
                    found.append((p, b_ok, n_c, lam, G))
                    break
        p = sympy.nextprime(p)
    # spread over lambda* deciles
    found.sort(key=lambda c: min(c[3], c[2] - c[3]) / c[2])
    if len(found) <= target:
        return found
    idx = [round(i * (len(found) - 1) / (target - 1)) for i in range(target)]
    return [found[i] for i in sorted(set(idx))]

def find_point(p, b):
    rng = random.Random(42)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

CURVES = find_17bit_curves(12)
print(f"\n  {len(CURVES)} curves, lambda* in "
      f"[{min(min(c[3], c[2]-c[3])/c[2] for c in CURVES):.4f}, "
      f"{max(min(c[3], c[2]-c[3])/c[2] for c in CURVES):.4f}]")

SEEDS5 = [42, 1234, 9999, 31337, 271828]
EFFS = [0.05, 0.15, 0.25, 0.35, 0.45]
M = 12

def k1_for_eff(n, eff):
    return max(2, int(round(eff * n / (math.isqrt(n) + 1))))

print("\n  Per-curve recovery (out of 5 seeds), m=12.  U = V0 uncentered, C = V1 centered.")
header = "  " + f"{'n':>7s} {'lam*':>7s} " + " ".join(f"eff={e:<4.2f}    " for e in EFFS)
print(header)
print("  " + " " * 15 + " ".join("  U    C     " for _ in EFFS))

tallies = {e: {'U': 0, 'C': 0} for e in EFFS}
rows = []
for c in CURVES:
    p, b, n, lam, G = c
    lstar = min(lam, n - lam) / n
    cells = []
    for e in EFFS:
        k1b = k1_for_eff(n, e)
        wu, t, _ = sweep('V0 base', c, M, k1b, SEEDS5)
        wc, _, _ = sweep('V1 base+ctr', c, M, k1b, SEEDS5)
        tallies[e]['U'] += (wu == t)
        tallies[e]['C'] += (wc == t)
        cells.append(f"{wu}/5  {wc}/5   ")
    rows.append((n, lstar, cells))
    print(f"  {n:7d} {lstar:7.4f} " + " ".join(cells))

print("\n  Curves recovering 5/5, out of %d:" % len(CURVES))
print("    eff:        " + " ".join(f"{e:>6.2f}" for e in EFFS))
print("    uncentered: " + " ".join(f"{tallies[e]['U']:>6d}" for e in EFFS))
print("    centered:   " + " ".join(f"{tallies[e]['C']:>6d}" for e in EFFS))
print("    [2026-07-29 uncentered, 20 curves: eff=0.05 -> 19, 0.15 -> 3, 0.25 -> 0]")

# ---------------------------------------------------------------------------
# lambda* correlation within each eff level
# ---------------------------------------------------------------------------
print("\n  Does lambda* still predict success at fixed eff?")
print(f"    {'eff':>6s} {'variant':>10s} {'lam* of successes':>26s} {'lam* of failures':>26s}")
for e in EFFS:
    ei = EFFS.index(e)
    for var, off in (('uncentered', 0), ('centered', 1)):
        succ = [f"{r[1]:.3f}" for r in rows
                if r[2][ei].split()[off].startswith('5/')]
        fail = [f"{r[1]:.3f}" for r in rows
                if not r[2][ei].split()[off].startswith('5/')]
        def rng(xs):
            return f"[{min(xs)}, {max(xs)}] n={len(xs)}" if xs else "-- none --"
        print(f"    {e:6.2f} {var:>10s} {rng(succ):>26s} {rng(fail):>26s}")

# ---------------------------------------------------------------------------
# Finite-m ceiling formula vs. observation
# ---------------------------------------------------------------------------
print("\n  Finite-m eff ceiling (the asymptote 3/(pi*e) is NOT attained at m=12):")
print("""
    GH_0 > ||pv||  <=>  eff  <  [2m/(2*pi*e)] / [C(m) * n^(1/m)]
        C(m) = m/6 + 1     (centered)
        C(m) = 2m/3 + 1    (uncentered)
    ratio centered/uncentered = (2m/3+1)/(m/6+1) -> 4 as m->inf, = 3.0 at m=12.
""")

def eff_ceiling(n, m, centered):
    C = (m / 6 + 1) if centered else (2 * m / 3 + 1)
    return (2 * m / (2 * math.pi * math.e)) / (C * n ** (1.0 / m))

print(f"    {'n':>8s} {'m':>3s} {'pred U':>8s} {'pred C':>8s} {'ratio':>6s}")
for n_ref in [2659, 2647, 75000]:
    print(f"    {n_ref:8d} {M:3d} {eff_ceiling(n_ref, M, False):8.4f} "
          f"{eff_ceiling(n_ref, M, True):8.4f} "
          f"{eff_ceiling(n_ref, M, True)/eff_ceiling(n_ref, M, False):6.2f}")
print("""
    Observed crossovers (last eff with majority 5/5 -> first eff with none):
      12-bit n=2659/2647, m=12   uncentered  0.079 -> 0.157   (pred 0.081)
      12-bit n=2659/2647, m=12   centered    0.313 -> 0.626   (pred 0.243)
      17-bit n~75000,     m=12   uncentered  0.05  -> 0.15    (pred 0.061)
      17-bit n~75000,     m=12   centered    0.15  -> 0.35    (pred 0.184)
""")

print("\n" + "=" * 78)
print("DONE")
print("=" * 78)
