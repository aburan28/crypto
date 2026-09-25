#!/usr/bin/env python3
"""Exact Frobenius-invariant linear-subspace admission check for ECC2K-130.

For prime extension degree n, ord_n(2) determines irreducible degrees of
(X^n-1)/(X-1) over F_2. This checks n=131 without a CAS.
"""
from math import gcd, isqrt

N = 131
PHI = N - 1
PRIME_DIVISORS = (2, 5, 13)

assert all(N % divisor for divisor in range(2, isqrt(N) + 1))
assert PHI == 2 * 5 * 13
assert all(pow(2, PHI // divisor, N) != 1 for divisor in PRIME_DIVISORS)
assert pow(2, PHI, N) == 1
assert gcd(N, 2) == 1
order = PHI

# X^131-1 has the linear factor X-1 and one degree-130 irreducible factor.
# Since the characteristic does not divide 131, each invariant subspace
# is a direct sum of some of these two irreducible Frobenius modules.
dimensions = {0}
for degree in (1, order):
    dimensions |= {existing + degree for existing in tuple(dimensions)}
assert dimensions == {0, 1, 130, 131}
print(f"ord_{N}(2)={order}; invariant F_2 subspace dimensions={sorted(dimensions)}")
