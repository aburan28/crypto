"""
GLV-HNP Phase 2, Thread 23 part 2: the K1 wall is an m-wall, and there is a
closed formula for the m you need.

Derivation (from T23-A of glv_hnp_phase2_thread23_bdd.py).  With
S_K1 = n/K1, S_K2 = n/K2, S_D = 1, S_KANNAN = n, dim = 2m+2:

    ||v_planted||  =  n * sqrt(2m/3 + 4/3)
    det(L)         =  (n*S_K1)^m * S_K2^m * n  =  n^(3m+1) / (K1*K2)^m
    GH(L)          =  sqrt(dim/(2*pi*e)) * det^(1/dim)

    gamma(m) = ||v_planted|| / GH(L)

    log gamma(m) = 1/2 * log( 2*pi*e*(m+2) / (3*(m+1)) )
                   + ( log n + m*log(eff) ) / (2m+2),      eff = K1*K2/n

Two consequences:

  (1) ASYMPTOTE.  gamma(m) -> sqrt(2*pi*e/3) * sqrt(eff) = 2.3862*sqrt(eff).
      So gamma can only ever drop below 1 if

          eff  <  3/(2*pi*e)  =  0.17564...

      This is a LATTICE-GEOMETRY wall, strictly stronger than the
      information-theoretic wall eff < 1 that the existing code uses via
      m_thresh = log(n)/log(1/eff) (glv_hnp_phase2_20bit.py:330).  Above
      eff_crit, NO amount of signatures helps.

  (2) FINITE m.  Setting gamma(m) < 1 and dropping the O(1/m) term:

          m_lat  =  ( log n + log(2*pi*e/3) ) / ( log(1/eff) - log(2*pi*e/3) )

      with log(2*pi*e/3) = 1.73937.  The convergence to the asymptote is
      governed by (n/eff)^(1/(2m+2)), i.e. it decays only like
      exp(log n / 2m) -- to get within 5% of the asymptote needs m ~ 10*log n.
      THIS is why the 2026-07-29 T4b sweep (m = 8/12/16/24/32 at K1=8, all
      failures) saw no improvement: it needed m ~ 87.

FALSIFIER (pre-registered, this is the whole point of the script):

  curve 12-bit/2677 (n=2647, K2=52, lam*=0.070), the historical hard curve.

    K1=6  -> eff=0.1179, m_lat = 24.1   PREDICT: fails at m=12, succeeds m>=24
    K1=8  -> eff=0.1572, m_lat = 86.7   PREDICT: fails at m<=32 (matches T4b),
                                                 succeeds at m>=87
    K1=12 -> eff=0.2357 > eff_crit      PREDICT: fails at EVERY m, incl. 128
                                                 (the sharp control)

  If K1=8 flips to success near m~87 while K1=12 stays dead at m=128, the
  formula is confirmed and the "K1 wall" of 2026-07-26/07-29 is dissolved into
  an m-wall with a closed form.  If K1=8 stays dead at m=128, the formula is
  wrong and the wall is something else.

Run: python3 glv_hnp_phase2_thread23_mwall.py
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL

from glv_hnp_phase2_thread23_bdd import (            # noqa: E402
    find_generator, gen_signatures, scales, build_glv_lattice,
    recover_kannan, gamma_ratio, planted_norm, gaussian_heuristic,
    lam_star, GAMMA_ASYMPTOTIC,
)

LOG_C = math.log(GAMMA_ASYMPTOTIC ** 2)        # log(2*pi*e/3) = 1.73937
EFF_CRIT = 3.0 / (2 * math.pi * math.e)        # 0.175639


def m_lat(n, eff):
    """Smallest m with gamma(m) < 1, from the closed form (O(1/m) dropped)."""
    denom = math.log(1.0 / eff) - LOG_C
    if denom <= 0:
        return float('inf')
    return (math.log(n) + LOG_C) / denom


def m_lat_exact(n, k1, k2, m_max=4096):
    """Smallest integer m with the exact gamma(m) < 1 (no asymptotics)."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    for m in range(2, m_max):
        if gamma_ratio(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN) < 1.0:
            return m
    return float('inf')


def m_info(n, eff):
    """Information-theoretic threshold used by the existing code."""
    return math.log(n) / math.log(1.0 / eff) if eff < 1 else float('inf')


def trial(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_glv_lattice(sigs, n, lam, S_K1, S_D, S_K2, S_KAN)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_kannan(reduced, m, n, S_D, S_KAN, d_secret)


# ===========================================================================
print("=" * 92)
print("Thread 23 part 2 — the K1 wall is an m-wall (closed-form m_lat)")
print("=" * 92)
print(f"eff_crit = 3/(2*pi*e) = {EFF_CRIT:.6f}   log(2*pi*e/3) = {LOG_C:.5f}")

HIST = [
    ("12-bit/2557",      2557, 2, 2659, 1755),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185),
]
curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0
    curves.append((label, (p, b, n, lam, G)))
    print(f"  {label}: n={n} lam*={lam_star(lam,n):.3f} G={G}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 92)
print("Predicted thresholds: m_info (information) vs m_lat (lattice geometry)")
print("-" * 92)
print(f"{'curve':>16} {'K1':>4} {'K2':>4} {'eff':>8} {'m_info':>8} "
      f"{'m_lat(asy)':>11} {'m_lat(exact)':>13}")
for label, (p, b, n, lam, G) in curves:
    k2 = math.isqrt(n) + 1
    for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
        eff = k1 * k2 / n
        ml = m_lat(n, eff)
        me = m_lat_exact(n, k1, k2)
        print(f"{label:>16} {k1:>4} {k2:>4} {eff:>8.4f} {m_info(n,eff):>8.2f} "
              f"{ml:>11.1f} {str(me):>13}")
    print()

# ---------------------------------------------------------------------------
print("-" * 92)
print("THE FALSIFIER — curve 12-bit/2677 (the historical hard curve), m sweep")
print("  2026-07-29 T4b baseline at K1=8: m = 8/12/16/24/32 -> 0,0,1,0,1 of 5")
print("-" * 92)

SEEDS = [42, 1234, 9999, 555, 31337]
hard = curves[1][1]
n_hard = hard[2]
k2_hard = math.isqrt(n_hard) + 1

PLAN = [
    (6,  [12, 16, 20, 24, 32, 48]),
    (8,  [12, 24, 32, 48, 64, 80, 96, 128]),
    (12, [12, 32, 64, 128]),           # eff > eff_crit: control, must all fail
]

results = {}
for k1, ms in PLAN:
    eff = k1 * k2_hard / n_hard
    ml = m_lat(n_hard, eff)
    me = m_lat_exact(n_hard, k1, k2_hard)
    print(f"\n=== K1={k1}  eff={eff:.4f}  m_info={m_info(n_hard,eff):.1f}  "
          f"m_lat(asy)={ml:.1f}  m_lat(exact)={me} ===")
    if eff > EFF_CRIT:
        print("    eff > eff_crit -> gamma(m) >= 1 for ALL m; predict total failure")
    print(f"{'m':>5} {'dim':>5} {'gamma':>8} {'pred':>6} {'LLL':>7} {'secs':>7}")
    S_K1, S_D, S_K2, S_KAN = scales(n_hard, k1, k2_hard)
    for m in ms:
        g = gamma_ratio(m, n_hard, k1, k2_hard, S_K1, S_D, S_K2, S_KAN)
        t0 = time.time()
        wins = tot = 0
        for s in SEEDS:
            r = trial(hard, m, k1, s)
            if r is None: continue
            tot += 1; wins += bool(r)
        dt = time.time() - t0
        print(f"{m:>5} {2*m+2:>5} {g:>8.4f} {str(g<1.0):>6} "
              f"{str(wins)+'/'+str(tot):>7} {dt:>7.1f}")
        results[(k1, m)] = (g, wins, tot)

# ---------------------------------------------------------------------------
print("\n" + "-" * 92)
print("VERDICT")
print("-" * 92)
for k1, ms in PLAN:
    eff = k1 * k2_hard / n_hard
    flip = None
    for m in ms:
        g, w, t = results[(k1, m)]
        if t and w == t:
            flip = m; break
    me = m_lat_exact(n_hard, k1, k2_hard)
    print(f"K1={k1:>3} eff={eff:.4f}  m_lat(exact)={str(me):>6}  "
          f"first m with 5/5 = {str(flip):>6}")
print("\nFormula confirmed iff (a) K1=6 and K1=8 flip at m near m_lat, and")
print("(b) K1=12 (eff > eff_crit) never flips, including at m=128.")
print("=" * 92)
