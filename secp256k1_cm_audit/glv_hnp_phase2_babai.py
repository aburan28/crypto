"""
Thread 23 (part 3) — the K1 wall is a Gram-Schmidt-profile condition.

Part 2 (glv_hnp_phase2_bdd.py) falsified the BDD-radius framing:

  * gamma = ||v_planted|| / lambda_1(L'_0) has NO common threshold across the
    two curves -- 12-bit/2557 flips between gamma 2.13 and 2.31, while
    12-bit/2677 flips between 0.71 and 0.86, a factor ~3 apart;
  * worse, 2557 recovers 5/5 at gamma = 2.13, i.e. FAR beyond the unique-decoding
    radius gamma < 1/2, and even beyond gamma = 1.  There exist noise vectors
    shorter than the planted error and recovery still succeeds.

So recovery is not "find the closest lattice point" either.  The mechanism must
be the coset marker: the Kannan column forces exactly one LLL basis vector into
the target coset, and LLL size-reduces it against the noise sublattice L'_0.
That succeeds exactly when Babai nearest-plane on the reduced basis of L'_0
returns the planted representative, i.e. when

    beta = max_i |<v_planted, b*_i>| / ||b*_i||^2  <  1/2

over the Gram-Schmidt vectors b*_i of the reduced L'_0 basis.  This is an exact
condition on the GS PROFILE, not on lambda_1 -- which is why lambda_1-based and
determinant-based (Gaussian-heuristic) predictors both failed.

  E8  does more signal (larger m) move the wall?  m in {8,...,40}, at each
      curve's wall K1.  Part 2 hinted 2557 improves with m while 2677 does not.
  E9  beta vs observed success on both curves -- is beta < 1/2 the universal
      predictor that gamma and tau were not?

Run: python3 glv_hnp_phase2_babai.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

from glv_hnp_phase2_projected import (
    find_generator, lam_star, scales, gen_signatures,
    build_projected_lattice, planted_vector_projected,
    run_projected, norm, rate,
)

SEEDS = [42, 1234, 9999, 555, 31337]
SEEDS10 = [42, 1234, 9999, 555, 31337, 8675309, 271828, 141421, 1618033, 6022140]

HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755, 0.3400),
    ("12-bit/2677", 2677, 2, 2647, 185,  0.0699),
]
curves = []
for label, p, b, n, lam, ls in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0
    curves.append((label, (p, b, n, lam, G)))


def gs_profile(rows):
    """Gram-Schmidt b*_i of the given (independent) integer rows, as floats."""
    bstar = []
    for v in rows:
        w = [float(x) for x in v]
        for bs in bstar:
            nb = sum(x * x for x in bs)
            if nb <= 0:
                continue
            c = sum(float(a) * b for a, b in zip(v, bs)) / nb
            w = [wi - c * bi for wi, bi in zip(w, bs)]
        if sum(x * x for x in w) > 1e-9:
            bstar.append(w)
    return bstar


def babai_beta(curve, m, d_secret, k1_bound, seed, beta_bkz=30):
    """
    beta = max_i |<v_planted, b*_i>| / ||b*_i||^2  over the reduced basis of the
    Kannan-free sublattice L'_0.  Babai nearest-plane recovers the planted coset
    representative exactly when beta < 1/2.
    """
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m:
        return None
    rows = build_projected_lattice(sigs, n, lam, k1_bound, k2b)
    cols = 2 * m + 1
    A = IntegerMatrix.from_matrix(rows[:-1])          # drop Kannan row
    BKZ.reduction(A, BKZ.Param(beta_bkz))
    red = [[A[i][j] for j in range(cols)] for i in range(A.nrows)]
    red = [r for r in red if any(r)]
    bstar = gs_profile(red)
    t = planted_vector_projected(sigs, n, k1_bound, k2b)
    worst = 0.0
    for bs in bstar:
        nb = sum(x * x for x in bs)
        if nb <= 0:
            continue
        c = abs(sum(float(a) * b for a, b in zip(t, bs))) / nb
        worst = max(worst, c)
    return worst


def mean_beta(curve, m, k1_bound, seeds):
    vals = []
    for s in seeds:
        d = random.Random(s + 7777).randint(1, curve[2] - 1)
        v = babai_beta(curve, m, d, k1_bound, s)
        if v is not None:
            vals.append(v)
    return sum(vals) / len(vals) if vals else float('nan')


def _main():
    print("=" * 78)
    print("Thread 23 part 3 — GS-profile (Babai) criterion for the Phase-2 K1 wall")
    print("=" * 78)

    # ---------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("E8 — does more signal move the wall?  (projected lattice, LLL, 10 seeds)")
    print("=" * 78)
    M_GRID = [8, 12, 16, 24, 32, 40]
    E8_K1 = {"12-bit/2557": [8, 12, 16], "12-bit/2677": [4, 6, 8]}
    for label, cur in curves:
        print(f"\n{label}   lam*={lam_star(cur[3], cur[2]):.4f}")
        print("  K1\\m   " + " ".join(f"{m:>7}" for m in M_GRID))
        for k1 in E8_K1[label]:
            row = []
            for m in M_GRID:
                w, t, _ = rate(run_projected, cur, m, k1, SEEDS10)
                row.append(w)
            print(f"  K1={k1:<4} " + " ".join(f"{str(w)+'/10':>7}" for w in row))

    # ---------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("E9 — beta = max_i |<v_planted,b*_i>|/||b*_i||^2   vs observed success")
    print("     Babai nearest-plane recovers exactly when beta < 0.5")
    print("=" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    brackets = {}
    for label, cur in curves:
        m = 8 if label == "12-bit/2557" else 10
        betas, wins = [], []
        for k1 in K1_GRID:
            betas.append(mean_beta(cur, m, k1, SEEDS))
            w, t, _ = rate(run_projected, cur, m, k1, SEEDS)
            wins.append(w)
        print(f"\n{label}  m={m}")
        print("  K1      " + " ".join(f"{k:>7}" for k in K1_GRID))
        print("  beta    " + " ".join(f"{x:>7.3f}" for x in betas))
        print("  wins    " + " ".join(f"{str(w)+'/5':>7}" for w in wins))
        lo = max((x for x, w in zip(betas, wins) if w == 5), default=float('nan'))
        hi = min((x for x, w in zip(betas, wins) if w < 5), default=float('nan'))
        brackets[label] = (lo, hi)
        print(f"  beta at last 5/5: {lo:.3f}   beta at first failure: {hi:.3f}")

    print("\n  --- is there a single beta threshold across both curves? ---")
    for label, (lo, hi) in brackets.items():
        print(f"  {label}: flip bracketed in ({lo:.3f}, {hi:.3f})")
    los = [v[0] for v in brackets.values()]
    his = [v[1] for v in brackets.values()]
    if max(los) < min(his):
        print(f"  CONSISTENT: a common threshold exists in "
              f"({max(los):.3f}, {min(his):.3f});  Babai predicts 0.5 -> "
              f"{'CONTAINS 0.5' if max(los) < 0.5 < min(his) else 'does NOT contain 0.5'}")
    else:
        print(f"  INCONSISTENT: brackets ({max(los):.3f},{min(his):.3f}) do not overlap; "
              f"beta is not a universal predictor either")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    _main()
