"""
GLV-HNP Phase 2, Thread 23 part 3: why MORE DATA HURTS, and the optimal m.

Part 2 (glv_hnp_phase2_thread23_mwall.py) FALSIFIED the pure Gaussian-heuristic
criterion gamma(m) < 1: on the hard curve 12-bit/2677 at K1=8, gamma drops below
1 at m=96 and keeps falling, yet recovery is 0/5 at m=96 AND at m=128.

The missing factor is that LLL's approximation quality degrades exponentially in
the lattice DIMENSION, and dim = 2m+2 grows with the data.  Under the standard
Gama-Nguyen unique-SVP heuristic, LLL/BKZ recovers the planted vector when

    lambda_2 / lambda_1  >~  tau * delta^dim        (delta ~ 1.021 for LLL)

Taking lambda_2 ~ GH(L) and lambda_1 = ||v_planted||, this is 1/gamma >~ tau*delta^dim,
i.e. the quantity to MINIMISE is

    Psi(m) = gamma(m) * delta^(2m+2)

gamma(m) decreases in m but saturates at sqrt(2*pi*e/3)*sqrt(eff); delta^(2m+2)
grows without bound.  So Psi has an INTERIOR MINIMUM.  Differentiating

    log Psi(m) = 1/2*log(2*pi*e*(m+2)/(3*(m+1))) + (L + m*E)/(2m+2) + (2m+2)*log(delta)
                 with L = log n, E = log(eff)

and dropping the O(1/m^2) term gives  d/dm = 2(E-L)/(2m+2)^2 + 2*log(delta) = 0, so

    dim*  =  2m* + 2  =  sqrt( log(n/eff) / log(delta) )                       (*)

    log Psi_min  =  1/2*log(2*pi*e/3) + 1/2*log(eff) + 2*sqrt( log(n/eff)*log(delta) )

PREDICTIONS (pre-registered):

  P1  On 12-bit/2677 (n=2647) at K1=6: dim* = sqrt(log(2647/0.1179)/log 1.021)
      = 22.0, i.e. m* = 10.0.  The success rate as a function of m must RISE to
      a peak near m ~ 10 and FALL thereafter -- not increase monotonically.
      The part-2 data (m = 12,16,20,24,32,48 -> 2,1,2,1,1,0 of 5) is consistent
      with the falling half; this script measures the rising half (m = 4..12),
      which has never been measured.

  P2  Psi_min, not eff, is the cross-K1 wall.  Predicted Psi_min on 2677:
      K1=2 -> 1.24, K1=4 -> 1.70, K1=6 -> 2.04, K1=8 -> 2.33.
      Observed 2026-07-29 T4 at m=12: 5/5, 5/5, 2/5, 0/5.  So the threshold
      should sit near Psi ~ 1.9.

  P3  If Psi is the right variable, then a STRONGER reduction (BKZ, smaller
      delta) must move dim* OUTWARD, i.e. the optimal m must grow.  With
      delta_BKZ20 ~ 1.0128, dim* grows by sqrt(0.02078/0.01272) = 1.28.

Run: python3 glv_hnp_phase2_thread23_psi.py
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL, BKZ

from glv_hnp_phase2_thread23_bdd import (            # noqa: E402
    find_generator, gen_signatures, scales, build_glv_lattice,
    recover_kannan, gamma_ratio, lam_star,
)

DELTA_LLL = 1.021        # Gama-Nguyen experimental Hermite factor for LLL
DELTA_BKZ20 = 1.0128     # ditto for BKZ-20


def psi(m, n, k1, k2, delta=DELTA_LLL):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    g = gamma_ratio(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN)
    return g * delta ** (2 * m + 2)


def dim_star(n, eff, delta=DELTA_LLL):
    return math.sqrt(math.log(n / eff) / math.log(delta))


def psi_min_closed(n, eff, delta=DELTA_LLL):
    return math.exp(0.5 * math.log(2 * math.pi * math.e / 3.0)
                    + 0.5 * math.log(eff)
                    + 2 * math.sqrt(math.log(n / eff) * math.log(delta)))


def trial(curve, m, k1_bound, seed, bkz_beta=0):
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
    if bkz_beta:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_kannan(reduced, m, n, S_D, S_KAN, d_secret)


def rate(curve, m, k1, seeds, bkz_beta=0):
    w = t = 0
    for s in seeds:
        r = trial(curve, m, k1, s, bkz_beta)
        if r is None: continue
        t += 1; w += bool(r)
    return w, t


# ===========================================================================
print("=" * 92)
print("Thread 23 part 3 — Psi(m) = gamma(m)*delta^dim  has an interior minimum")
print("=" * 92)

HIST = [
    ("12-bit/2557",      2557, 2, 2659, 1755),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185),
]
curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0
    curves.append((label, (p, b, n, lam, G)))
    print(f"  {label}: n={n} lam*={lam_star(lam,n):.3f}")

SEEDS = [42, 1234, 9999, 555, 31337, 7, 271828, 141421, 100003, 65537]
print(f"  {len(SEEDS)} seeds per cell (doubled vs part 2 to resolve the peak)")

# ---------------------------------------------------------------------------
print("\n" + "-" * 92)
print("P2: closed-form Psi_min vs the 2026-07-29 T4 outcomes at m=12")
print("-" * 92)
print(f"{'curve':>16} {'K1':>4} {'eff':>7} {'dim*':>7} {'m*':>6} "
      f"{'Psi_min':>8} {'Psi(12)':>8} {'T4@m=12':>9}")
T4 = {('12-bit/2557', 2): '5/5', ('12-bit/2557', 3): '5/5',
      ('12-bit/2557', 4): '5/5', ('12-bit/2557', 6): '5/5',
      ('12-bit/2557', 8): '5/5', ('12-bit/2557', 12): '4/5',
      ('12-bit/2557', 16): '1/5', ('12-bit/2557', 24): '0/5',
      ('12-bit/2677 FAIL', 2): '5/5', ('12-bit/2677 FAIL', 3): '5/5',
      ('12-bit/2677 FAIL', 4): '5/5', ('12-bit/2677 FAIL', 6): '2/5',
      ('12-bit/2677 FAIL', 8): '0/5', ('12-bit/2677 FAIL', 12): '0/5',
      ('12-bit/2677 FAIL', 16): '0/5', ('12-bit/2677 FAIL', 24): '0/5'}
for label, (p, b, n, lam, G) in curves:
    k2 = math.isqrt(n) + 1
    for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
        eff = k1 * k2 / n
        ds = dim_star(n, eff)
        print(f"{label:>16} {k1:>4} {eff:>7.4f} {ds:>7.1f} {(ds-2)/2:>6.1f} "
              f"{psi_min_closed(n, eff):>8.3f} {psi(12, n, k1, k2):>8.3f} "
              f"{T4.get((label,k1),'-'):>9}")
    print()

# ---------------------------------------------------------------------------
print("-" * 92)
print("P1: THE PEAK.  m sweep across dim* on the hard curve (12-bit/2677).")
print("    Part 2 measured only m >= 12 (the falling half).  m = 4..12 is new.")
print("-" * 92)

hard = curves[1][1]
n_h = hard[2]
k2_h = math.isqrt(n_h) + 1
M_SWEEP = [4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24, 32, 48]

for k1 in (6, 8):
    eff = k1 * k2_h / n_h
    ds = dim_star(n_h, eff)
    print(f"\n=== K1={k1}  eff={eff:.4f}  dim*={ds:.1f}  m*={(ds-2)/2:.1f}  "
          f"Psi_min={psi_min_closed(n_h, eff):.3f} ===")
    print(f"{'m':>4} {'dim':>5} {'gamma':>8} {'Psi':>8} {'LLL':>7} {'bar':>22}")
    best = (-1, None)
    for m in M_SWEEP:
        S_K1, S_D, S_K2, S_KAN = scales(n_h, k1, k2_h)
        g = gamma_ratio(m, n_h, k1, k2_h, S_K1, S_D, S_K2, S_KAN)
        w, t = rate(hard, m, k1, SEEDS)
        mark = " <-- dim*" if abs(2*m+2 - ds) <= 2 else ""
        print(f"{m:>4} {2*m+2:>5} {g:>8.4f} {psi(m,n_h,k1,k2_h):>8.3f} "
              f"{str(w)+'/'+str(t):>7} {'#'*w:>22}{mark}")
        if w > best[0]:
            best = (w, m)
    print(f"  peak: m={best[1]} with {best[0]}/{len(SEEDS)}   "
          f"(predicted m* = {(ds-2)/2:.1f})")

# ---------------------------------------------------------------------------
print("\n" + "-" * 92)
print("P3: stronger reduction must move the peak OUTWARD.")
print(f"    dim*(BKZ20)/dim*(LLL) = sqrt(log{DELTA_LLL}/log{DELTA_BKZ20}) = "
      f"{math.sqrt(math.log(DELTA_LLL)/math.log(DELTA_BKZ20)):.3f}")
print("-" * 92)

k1 = 8
eff = k1 * k2_h / n_h
ds_l = dim_star(n_h, eff, DELTA_LLL)
ds_b = dim_star(n_h, eff, DELTA_BKZ20)
print(f"K1={k1}: predicted m*(LLL) = {(ds_l-2)/2:.1f}, "
      f"m*(BKZ20) = {(ds_b-2)/2:.1f}")
print(f"{'m':>4} {'dim':>5} {'LLL':>7} {'BKZ20':>7} {'secs':>7}")
for m in (6, 8, 10, 12, 16, 20, 24, 32):
    t0 = time.time()
    wl, tl = rate(hard, m, k1, SEEDS)
    wb, tb = rate(hard, m, k1, SEEDS, bkz_beta=20)
    print(f"{m:>4} {2*m+2:>5} {str(wl)+'/'+str(tl):>7} "
          f"{str(wb)+'/'+str(tb):>7} {time.time()-t0:>7.1f}")

print("\n" + "=" * 92)
print("done")
print("=" * 92)
