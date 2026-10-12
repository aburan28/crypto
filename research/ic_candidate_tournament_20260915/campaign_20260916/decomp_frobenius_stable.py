#!/usr/bin/env python3
"""Which binary Koblitz degrees admit a Frobenius-stable factor base, and of what size.

An F_2-subspace V of F_{2^n} is stable under x -> x^2 exactly when V = ker g(sigma)
for a divisor g of t^n - 1 over F_2 (research/notes/index-calculus/RESEARCH_QUASI_SUBFIELD.md
section 3).  For n prime, t^n - 1 = (t - 1) * Phi_n(t), and Phi_n splits over F_2 into
(n - 1)/o irreducible factors of degree o = ord_n(2).  So the achievable stable dimensions
are exactly { e + k*o : e in {0, 1}, 0 <= k <= (n-1)/o }.

For each degree this prints o, the achievable dimensions, and -- for m summands -- the
stable dimension closest to the two model optima of DECOMPOSITION-SURVEY.md section 2:
  * l_bal  = n/(m+1)   (free-oracle balance of relation collection against linear algebra)
  * l_sat  = n/m       (the largest l before every target decomposes, m*l ~ n)
A '-' means the only stable dimensions are the trivial ones {0, 1, n-1, n}.

Exact integer arithmetic; no randomness; runs in well under a second.
"""


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def order_of_two(n):
    k, x = 1, 2 % n
    while x != 1:
        x = x * 2 % n
        k += 1
    return k


def stable_dims(n):
    o = order_of_two(n)
    blocks = (n - 1) // o
    return o, sorted({e + k * o for e in (0, 1) for k in range(blocks + 1)})


def nearest_nontrivial(dims, n, target):
    useful = [d for d in dims if 2 <= d <= n - 2]
    if not useful:
        return None
    return min(useful, key=lambda d: (abs(d - target), d))


TOURNAMENT = [13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
CRYPTO = [73, 127, 131, 151, 163, 233, 239, 257, 283, 409, 571]

assert all(is_prime(n) for n in TOURNAMENT + CRYPTO)

print("n    ord_n(2)  #Phi_n factors  stable dims (nontrivial, 2..n-2)            "
      "m=3 bal/sat   m=4 bal/sat   m=5 bal/sat")
for label, ns in (("tournament degrees", TOURNAMENT), ("cryptographic / reference degrees", CRYPTO)):
    print(f"-- {label}")
    for n in ns:
        o, dims = stable_dims(n)
        nontriv = [d for d in dims if 2 <= d <= n - 2]
        shown = ",".join(map(str, nontriv)) if nontriv else "-"
        if len(shown) > 44:
            shown = shown[:41] + "..."
        cells = []
        for m in (3, 4, 5):
            b = nearest_nontrivial(dims, n, n / (m + 1))
            s = nearest_nontrivial(dims, n, n / m)
            cells.append(f"{'-' if b is None else b:>3}/{'-' if s is None else s:<3}")
        print(f"{n:<4} {o:>8}  {(n - 1) // o:>14}  {shown:<44}  " + "     ".join(cells))

print()
print("Reading: a stable base folds the factor base, the relation count and the linear algebra by n")
print("(the orbit length), a poly(n) factor that does not move the exponent.  Where the nearest stable")
print("dimension is far from n/(m+1) or n/m, the fold costs more than it buys; where it is '-', it")
print("does not exist and only an orbit-union base (materialised, not a subspace) is available.")
