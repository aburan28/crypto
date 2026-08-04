"""
GLV-HNP Phase 2, Thread 23 / E5: is the K1 wall information-theoretic or
algorithmic?  Decided by exhaustive search over d.

The 2026-07-29 entry framed Thread 23 as a dichotomy:
    "if the wall stays at K1 ~ 4-6, then the wall is information-theoretic
     and Phase 2 is at its ceiling."
E1-E4 showed the wall does stay put (OLD and NEW grids bit-identical, 10
seeds, LLL/BKZ-20/BKZ-40 alike).  But "the reformulation didn't help" does
NOT imply "the information isn't there" -- that inference needs its own test.

E5 runs it directly.  The historical failure curve has n = 2647, so the whole
key space is enumerable.  For each candidate d in [1, n) compute
    k_full_i = A_i + B_i*d  mod n
and ask whether every signature admits a GLV decomposition with the promised
bounds, i.e. whether there is k2 in [0,K2) with (k_full_i - lam*k2) mod n < K1.
d_secret always passes.  The question is how many OTHER d also pass.

  #survivors == 1  ->  the instance is information-theoretically determined;
                       LLL/BKZ failing on it is an ALGORITHMIC limit, and the
                       "information-theoretic ceiling" reading is FALSIFIED.
  #survivors >> 1  ->  the information really is absent and Phase 2 is at its
                       ceiling, as the 2026-07-29 entry proposed.

E5b compares the counting bound already coded into `sweep_curve`
(m_thresh = ceil(log n / log(1/eff)), glv_hnp_phase2_20bit.py) against both
the measured survivor count and the measured lattice wall.

Run: python3 glv_hnp_phase2_projected_e5.py
Deps: fpylll, sympy
"""

import math

from glv_hnp_phase2_projected import (
    scales, gen_signatures, find_generator, lam_star, rate,
    HIST, SEEDS, D_SECRET_FRAC,
)


def decomposable(k_full, n, lam, K1, K2):
    """Is k_full = k1 + lam*k2 mod n with 0 <= k1 < K1, 0 <= k2 < K2?"""
    for k2 in range(K2):
        if (k_full - lam * k2) % n < K1:
            return True
    return False


def survivors(sigs, n, lam, K1, K2, d_secret, cap=None):
    """All d in [1,n) consistent with every signature's GLV bound."""
    out = []
    for d in range(1, n):
        ok = True
        for s in sigs:
            k_full = (s['A'] + s['B'] * d) % n
            if not decomposable(k_full, n, lam, K1, K2):
                ok = False
                break
        if ok:
            out.append(d)
            if cap and len(out) >= cap:
                break
    assert d_secret in out, "d_secret must always survive"
    return out


hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    hist.append((label, (p, b, n, lam, G), k1, m))
CURVES12 = [(lb, cv, m) for (lb, cv, k1, m) in hist if '12-bit' in lb]

print("=" * 78)
print("Thread 23 / E5 - information-theoretic or algorithmic?  (exhaustive d)")
print("=" * 78)

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]

for label, curve, m in CURVES12:
    p, b, n, lam, G = curve
    d_secret = max(1, int(D_SECRET_FRAC * n))
    K2 = math.isqrt(n) + 1
    print()
    print(f"  {label}  n={n}  lam*={lam_star(lam, n):.4f}  m={m}  K2={K2}  "
          f"d_secret={d_secret}")
    print(f"  {'K1':>3} {'eff':>6} {'m_thresh':>9} | {'#survivors':>11} "
          f"{'unique?':>8} | {'LLL':>4} {'B40':>4} | {'verdict':>26}")
    print("  " + "-" * 82)
    for K1 in K1_GRID:
        eff = K1 * K2 / n
        m_thresh = (math.ceil(math.log(n) / math.log(1.0 / eff))
                    if eff < 1.0 else float('inf'))
        sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, 42)
        surv = survivors(sigs, n, lam, K1, K2, d_secret)
        uniq = (len(surv) == 1)
        g_lll, _ = rate(curve, m, d_secret, K1, SEEDS, 'new')
        g_b40, _ = rate(curve, m, d_secret, K1, SEEDS, 'new',
                        use_bkz=True, beta=40)
        if uniq and g_b40 == 0:
            verdict = "ALGORITHMIC failure"
        elif uniq:
            verdict = "solved, info unique"
        elif g_b40 == 0:
            verdict = "info genuinely absent"
        else:
            verdict = "solved despite ambiguity"
        print(f"  {K1:>3} {eff:>6.3f} {str(m_thresh):>9} | {len(surv):>11} "
              f"{str(uniq):>8} | {g_lll:>4} {g_b40:>4} | {verdict:>26}")

print()
print("=" * 78)
print("done")
print("=" * 78)
