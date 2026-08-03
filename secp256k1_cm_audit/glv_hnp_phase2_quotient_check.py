"""
Thread 23, structural check: L0 / L_lambda = Z/n, with d as the coset label.

This is the statement that explains every negative result of the session:

  * L_lam := direct sum over i of the 2D lambda-block L2^(i)
             = < (n*S_K1, 0), (-lam*S_K1, S_K2) >  in columns (i, m+i)
    is a FINITE-INDEX sublattice of L0 (the Kannan-free dim-2m sublattice),
    of index exactly n.
  * Every vector of L_lam leaves the recovered d unchanged.
  * Hence L0 / L_lam = Z/nZ and the coset label IS the candidate d.

Consequences:
  - the m lambda-block rivals that keep the planted vector off lambda_1 are
    HARMLESS: they never change the answer, so removing them (Thread 23's
    proposed reformulation) cannot change recovery.  Confirmed empirically:
    300/300 instance-level agreement between the original and projected arms.
  - block-wise disagreement is impossible: either all m blocks report the true
    d or none do.  Confirmed empirically: clean-block counts are bimodal at
    {0, m}, and majority voting gains exactly nothing.
  - recovery is precisely "does the closest point of L0 to the target lie in
    the correct one of n cosets", i.e. a decoding problem over Z/n in a direct
    sum of m two-dimensional lattices.

Run: python3 glv_hnp_phase2_quotient_check.py
"""

import math
import random

from sympy import Matrix

from glv_hnp_phase2_projected import scales, gen_signatures
from glv_hnp_phase2_exact_cvp import l0_basis, search_curves


def main():
    print("=" * 78)
    print("Thread 23 -- L0 / L_lambda = Z/n, with d as the coset label")
    print("=" * 78)
    print()
    print("Algebra.  A lambda-block vector in block i is")
    print("  c1*(n*S_K1, 0) + c2*(-lam*S_K1, S_K2)")
    print("     = ((c1*n - c2*lam)*S_K1,  c2*S_K2),")
    print("so  dk1_i = c1*n - c2*lam  and  dk2_i = c2, hence")
    print("  dk1_i + lam*dk2_i = c1*n - c2*lam + lam*c2 = c1*n = 0  (mod n).")
    print()
    print("The readout d_i = (k1_i + lam*k2_i - A_i) * B_i^-1  (mod n) depends")
    print("on (k1_i, k2_i) only through k1_i + lam*k2_i mod n, so d_i is")
    print("INVARIANT under all of L_lam.  QED.")
    print()
    print("Numerical check of the index [L0 : L_lam]  (expect exactly n):")
    print(f"{'n':>7} {'eff':>5} {'m':>3} | {'index':>10} {'index == n':>11}")
    print("-" * 46)

    ok = True
    for (p, b, n, lam, G) in search_curves(2 ** 12, 2 ** 13, want=3):
        for eff, m in ((0.10, 4), (0.20, 5)):
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
            d0 = random.Random(42).randint(1, n - 1)
            sigs = gen_signatures(G, d0, m, n, lam, p, k1b, k2b, 42)
            if len(sigs) < m:
                continue
            det0 = int(abs(Matrix(l0_basis(sigs, n, lam, k1b, k2b)).det()))
            detl = (n * S_K1 * S_K2) ** m
            idx = detl // det0
            good = (detl % det0 == 0 and idx == n)
            ok &= good
            print(f"{n:>7} {eff:>5.2f} {m:>3} | {idx:>10} {str(good):>11}")

    print()
    print(f"all indices equal n: {ok}")
    print("=" * 78)


if __name__ == "__main__":
    main()
