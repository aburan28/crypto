#!/usr/bin/env python3
"""Exact Frobenius-invariant linear-subspace admission check for ECC2K-130.

For prime extension degree n, ord_n(2) determines irreducible degrees of
(X^n-1)/(X-1) over F_2. This checks the n=131 case without a CAS.
"""
from math import gcd

N = 131
PHI = N - 1
PRIME_DIVISORS = (2, 5, 13)

assert all(PHI % p == 0 for p in PRIME_DIVISORS)
assert all(pow(2, PHI // p, N) != 1 for p in PRIME_DIVISORS)
assert pow(2, PHI, N) == 1
assert gcd(N, 2) == 1
ORDER = PHI
INVARIANT_DIMENSIONS = (0, 1, 130, 131)
print(f"ord_{N}(2)={ORDER}; invariant F_2 subspace dimensions={INVARIANT_DIMENSIONS}")
