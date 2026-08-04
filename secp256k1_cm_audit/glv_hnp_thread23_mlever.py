"""
Thread 23, part 4 — is the wall movable by more signatures?

Parts 1-3 established:
  * the trivial vector n*S_D*e_m is NOT the obstruction (S_D sweep and the
    embedding-free Babai/CVP arms all fail to beat the baseline);
  * beyond the wall, exact CVP returns a decoy that is STRICTLY closer to the
    target than the planted error (U5), so the failure is decoding-uniqueness,
    not reduction quality (U9: LLL wall == CVP wall on 6/6 fresh curves).

The remaining lever in the Phase-2 lattice is m.  Adding a signature adds one
constraint but also two dimensions (one k1 column, one k2 column), so it is not
obvious that decoy density falls.  EXP T4b (2026-07-29) tried m = 8..32 with
Kannan+LLL at K1=8 on 12-bit/2677 and got 0,0,1,0,1 of 5 — but LLL there could
not see past its own wall.  Re-run the m lever with the exact decoder.

U10  m in {6,10,14,18} x K1 in {6,8,12} on 12-bit/2677, arms Kannan+LLL and
     exact CVP; also report how often a strictly-closer decoy exists.

Enumeration cost grows fast with dim = 2m+1, so each CVP call is capped by
SIGALRM (CAP_SEC) and timeouts are reported rather than silently dropped.

Run: python3 glv_hnp_thread23_mlever.py
"""

import math
import random
import signal

from fpylll import IntegerMatrix, LLL, BKZ, CVP

from glv_hnp_thread23_bdd import (
    find_generator, lam_star, gen_signatures, scales,
    build_kannan, recover_kannan, build_bdd, error_vector,
    recover_bdd, norm,
)

SEEDS = [42, 1234, 9999, 555, 31337]
CAP_SEC = 180


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


signal.signal(signal.SIGALRM, _alarm)

print("=" * 78)
print("Thread 23 part 4 — does m move the wall?")
print("=" * 78)

p, b, n, lam = 2677, 2, 2647, 185
G = find_generator(p, b, n)
assert G is not None and (lam * lam + lam + 1) % n == 0
curve = (p, b, n, lam, G)
print(f"  12-bit/2677  n={n}  lam={lam}  lam*={lam_star(lam,n):.4f}")
print(f"  CVP cap: {CAP_SEC}s per call\n")

M_GRID = [6, 10, 14, 18]
K1_GRID = [6, 8, 12]

print("  m   dim  K1   K+LLL   CVP     decoy-closer  timeouts")
print("  " + "-" * 62)
for m in M_GRID:
    for k1 in K1_GRID:
        k2b = math.isqrt(n) + 1
        lll_w, cvp_w, closer, touts = 0, 0, 0, 0
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, s)
            if len(sigs) < m:
                continue
            # arm 1: Kannan + LLL
            MK = build_kannan(sigs, n, lam, k1, k2b, 1)
            dk = 2 * m + 2
            A = IntegerMatrix.from_matrix(MK)
            LLL.reduction(A)
            rows = [[A[i][j] for j in range(dk)] for i in range(dk)]
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b, 1)
            lll_w += recover_kannan(rows, m, n, S_KAN, S_D, d) is not None
            # arm 2: exact CVP
            MB, t = build_bdd(sigs, n, lam, k1, k2b, 1)
            db = 2 * m + 1
            AB = IntegerMatrix.from_matrix(MB)
            BKZ.reduction(AB, BKZ.Param(20))     # cheap pre-reduction
            e1 = error_vector(sigs, d, n, k1, k2b, 1, reduced=True)
            e0 = error_vector(sigs, d, n, k1, k2b, 1, reduced=False)
            en = min(norm(e1), norm(e0))
            try:
                signal.alarm(CAP_SEC)
                w = list(CVP.closest_vector(AB, tuple(t)))
                signal.alarm(0)
            except Timeout:
                signal.alarm(0)
                touts += 1
                continue
            except Exception:
                signal.alarm(0)
                touts += 1
                continue
            cvp_w += recover_bdd(w, t, m, n, S_D, d) is not None
            if norm([w[j] - t[j] for j in range(db)]) < en - 1e-9:
                closer += 1
        print(f"  {m:<3d} {2*m+1:<4d} {k1:<4d} {lll_w}/5     {cvp_w}/5     "
              f"{closer}/5           {touts}")
    print()

print("=" * 78)
print("done")
print("=" * 78)
