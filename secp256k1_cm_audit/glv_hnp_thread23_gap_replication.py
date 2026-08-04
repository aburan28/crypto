"""
Thread 23, part 5 — does the reduction gap appear on FRESH curves at larger m?

Part 4 (glv_hnp_thread23_mlever.py) found on 12-bit/2677 (lam*=0.070) that the
exact-CVP wall moves outward with m while the Kannan+LLL wall does not:

    K1=8:  m=6 -> LLL 1/5 CVP 0/5 ;  m=18 -> LLL 0/5 CVP 4/5

Part 3 (U9) found NO gap on 6 fresh 17-bit curves — but it only tested m=10.
So "the gap is idiosyncratic to 2677" and "the gap only opens at larger m" are
both still live.  This script separates them: two fresh 17-bit curves (one
small lam*, one mid lam*), m in {10, 18}, at the eff step just past each
curve's m=10 LLL wall.

Deliberately small: 3 seeds, 2 curves, 2 m values, 2 eff steps, CVP capped at
CAP_SEC per call.  Timeouts are counted, not hidden.

Run: python3 glv_hnp_thread23_gap_replication.py
"""

import math
import random
import signal

from fpylll import IntegerMatrix, LLL, BKZ, CVP

from glv_hnp_thread23_bdd import (
    gen_signatures, scales, build_kannan, recover_kannan,
    build_bdd, error_vector, recover_bdd, norm, lam_star,
)
from glv_hnp_thread23_replication import search_curves, nu_hat

SEEDS = [42, 1234, 9999]
CAP_SEC = 150


class Timeout(Exception):
    pass


signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(Timeout()))

print("=" * 78)
print("Thread 23 part 5 — gap replication at larger m")
print("=" * 78)
print("Searching fresh 17-bit j=0 GLV curves ...")
FRESH = search_curves(2 ** 16, 2 ** 17, nbins=6, per_bin=1)
# keep the smallest-lam* and a mid-lam* curve
FRESH.sort(key=lambda c: lam_star(c[3], c[2]))
PICK = [FRESH[0], FRESH[len(FRESH) // 2]]
for (p, b, n, lam, G) in PICK:
    k2b = math.isqrt(n) + 1
    print(f"  p={p} n={n} lam*={lam_star(lam,n):.4f} "
          f"nu_hat={nu_hat(n, lam, n//8, max(1, n//k2b)):.4f}")

EFFS = [0.08, 0.12]
M_GRID = [10, 18]
print(f"\n{len(SEEDS)} seeds, CVP cap {CAP_SEC}s\n")
print("  p         lam*    m   eff   dim   K+LLL  CVP    decoy-closer  timeouts")
print("  " + "-" * 72)
for (p, b, n, lam, G) in PICK:
    k2b = math.isqrt(n) + 1
    ls = lam_star(lam, n)
    for m in M_GRID:
        for eff in EFFS:
            k1 = max(2, int(round(eff * n / k2b)))
            lll_w, cvp_w, closer, touts = 0, 0, 0, 0
            for s in SEEDS:
                d = random.Random(s + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
                if len(sigs) < m:
                    continue
                MK = build_kannan(sigs, n, lam, k1, k2b, 1)
                dk = 2 * m + 2
                A = IntegerMatrix.from_matrix(MK)
                LLL.reduction(A)
                rows = [[A[i][j] for j in range(dk)] for i in range(dk)]
                S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b, 1)
                lll_w += recover_kannan(rows, m, n, S_KAN, S_D, d) is not None

                MB, t = build_bdd(sigs, n, lam, k1, k2b, 1)
                db = 2 * m + 1
                AB = IntegerMatrix.from_matrix(MB)
                BKZ.reduction(AB, BKZ.Param(20))
                e1 = error_vector(sigs, d, n, k1, k2b, 1, reduced=True)
                e0 = error_vector(sigs, d, n, k1, k2b, 1, reduced=False)
                en = min(norm(e1), norm(e0))
                try:
                    signal.alarm(CAP_SEC)
                    w = list(CVP.closest_vector(AB, tuple(t)))
                    signal.alarm(0)
                except Exception:
                    signal.alarm(0)
                    touts += 1
                    continue
                cvp_w += recover_bdd(w, t, m, n, S_D, d) is not None
                if norm([w[j] - t[j] for j in range(db)]) < en - 1e-9:
                    closer += 1
            print(f"  {p:<9d} {ls:.4f}  {m:<3d} {eff:.2f}  {db:<4d}  "
                  f"{lll_w}/{len(SEEDS)}    {cvp_w}/{len(SEEDS)}    "
                  f"{closer}/{len(SEEDS)}           {touts}")
    print()

print("=" * 78)
print("done")
print("=" * 78)
