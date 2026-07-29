"""
GLV-HNP Phase 2 / Thread 20 - robustness confirmation for the lambda/n falsification.

The main study (glv_hnp_phase2_lambda_threshold.py) uses 3 seeds per (curve, m).
This script re-tests the decisive curves with 12 seeds so the falsification of the
2026-07-26 "lam/n threshold" claim does not rest on seed luck.

Decisive pairs, all at eff = K1*K2/n ~ 0.157:

  A. n=2551, lam/n=0.0196  - FAR below the claimed (0.07,0.34) threshold.
     Claim predicts FAILURE.  Main study measured 3/3 at m=7.
  B. n=2647, lam/n=0.0699  - the single failing curve the claim was built on.
     Claim predicts FAILURE.  Main study reproduced the failure.
  C. n=2251, lam/n=0.3145  - ABOVE the claimed threshold.
     Claim predicts SUCCESS.  Main study measured failure at every m<=12.
  D. n=2293, lam/n=0.4313  - above threshold; main study measured 3/3 at m=8.

If A succeeds and C fails under 12 seeds, lam/n is falsified in both directions.

Run: python3 glv_hnp_phase2_lambda_confirm.py
"""

import math

import sympy

from glv_hnp_phase2_lambda_threshold import (
    find_generator, glv_eigenvalues, first_full_recovery_m, predictors,
    k1_for_eff, eisenstein_decompose, j0_traces,
)


def find_p_for_n(n, p_lo, p_hi):
    """Recover the base prime p of the j=0 CM curve with group order n.

    n = p + 1 - t for one of the six sextic-twist traces of p, so scan the
    primes in [p_lo, p_hi) and check whether n appears among them.
    """
    p = int(sympy.nextprime(p_lo - 1))
    while p < p_hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None and any(p + 1 - t == n for t in j0_traces(*eis)):
                return p
        p = int(sympy.nextprime(p))
    return None

SEEDS12 = [42, 1234, 9999, 7, 31337, 555, 90210, 4242, 13, 271828, 161803, 8675309]

# (label, n, expected verdict of the 2026-07-26 lam/n claim)
CASES = [
    ("A n=2551  lam/n=0.020", 2551, "FAILURE"),
    ("B n=2647  lam/n=0.070", 2647, "FAILURE"),
    ("C n=2251  lam/n=0.315", 2251, "SUCCESS"),
    ("D n=2293  lam/n=0.431", 2293, "SUCCESS"),
]


def main():
    print("=" * 78)
    print("Thread 20 confirmation - 12 seeds, 12/12 required for SUCCESS")
    print("=" * 78)
    print(f"  {'case':<24} {'K1':>4} {'eff':>6} {'lam/n':>7} {'m*':>5} "
          f"{'claim':>8} {'measured':>9}  verdict")

    rows = []
    for label, n, claim in CASES:
        lam, _ = glv_eigenvalues(n)
        if lam is None:
            print(f"  {label:<24} no GLV eigenvalue - skipped")
            continue
        p = find_p_for_n(n, 2**11, 2**12)
        if p is None:
            print(f"  {label:<24} no base prime found - skipped")
            continue
        G, b_found = None, None
        for b_try in range(1, 60):
            G = find_generator(p, b_try, n)
            if G is not None:
                b_found = b_try
                break
        if G is None:
            print(f"  {label:<24} no generator found - skipped")
            continue
        curve = (p, b_found, n, lam, G)
        k1 = k1_for_eff(n, 0.157)
        hit, tally, pn = first_full_recovery_m(curve, k1, range(4, 15),
                                               seeds=SEEDS12)
        pr = predictors(n, lam, k1, 14, pn)
        measured = "SUCCESS" if hit is not None else "FAILURE"
        agree = "claim ok" if measured == claim else "CLAIM WRONG"
        print(f"  {label:<24} {k1:>4} {pr['eff']:>6.3f} {pr['lam_n']:>7.4f} "
              f"{str(hit):>5} {claim:>8} {measured:>9}  {agree}")
        print(f"      tally (m: wins/12): "
              + ' '.join(f"{m}:{w}" for m, w in sorted(tally.items())))
        rows.append((label, claim, measured))

    wrong = [r for r in rows if r[1] != r[2]]
    print(f"\n  {len(wrong)}/{len(rows)} decisive cases contradict the "
          f"2026-07-26 lam/n claim.")
    for label, claim, measured in wrong:
        print(f"    {label}: claim said {claim}, measured {measured}")


if __name__ == '__main__':
    main()
