#!/usr/bin/env python3
"""Embedding degree of the large prime factor of #E(F_{2^m}) for

    E_0: y^2 + xy = x^3 + 1

over GF(2). The trace recurrence is the one in
`koblitz_point_count` (`src/cryptanalysis/koblitz_index_calculus.rs`):
s_0 = 2, s_1 = t = -1, s_k = t*s_{k-1} - 2*s_{k-2}, and
#E = 2^m + 1 - s_m. The embedding degree is the multiplicative
order of 2^m modulo that prime. MOV/Frey-Rueck only helps when this
order is tiny.
"""

from __future__ import annotations


def factor(n: int) -> list[tuple[int, int]]:
    n = abs(n)
    factors: list[tuple[int, int]] = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            c = 0
            while n % p == 0:
                n //= p
                c += 1
            factors.append((p, c))
        p += 1 if p == 2 else 2
    if n > 1:
        factors.append((n, 1))
    return factors


def trace(m: int) -> int:
    s_prev, s_cur = 2, -1
    if m == 0:
        return 2
    for _ in range(1, m):
        s_prev, s_cur = s_cur, -s_cur - 2 * s_prev
    return s_cur


def group_order(m: int) -> int:
    return (1 << m) + 1 - trace(m)


def multiplicative_order(a: int, r: int) -> int:
    k = r - 1
    for p, _ in factor(r - 1):
        while k % p == 0 and pow(a, k // p, r) == 1:
            k //= p
    return k


def rows(degrees: range | list[int]):
    out = []
    for m in degrees:
        order = group_order(m)
        factors = factor(order)
        prime = max(p for p, _ in factors)
        degree = multiplicative_order(pow(2, m, prime), prime)
        assert pow(pow(2, m, prime), degree, prime) == 1
        out.append((m, order, factors, prime, degree))
    return out


def main() -> None:
    assert group_order(1) == 4
    assert group_order(2) == 8
    data = rows([5, 7, 11, 13, 15, 17, 19, 23, 29, 31, 41])
    print(f"{'m':>3} {'bits(r)':>7} {'k':>12} {'log2(k)':>8}")
    for m, _order, _factors, prime, degree in data:
        bits = prime.bit_length()
        logk = degree.bit_length() - 1
        print(f"{m:3d} {bits:7d} {degree:12d} {logk:8d}")
    by_m = {m: (prime, degree) for m, _o, _f, prime, degree in data}
    # Past the toy rungs, k is far above the supersingular range k <= 6.
    for m in (19, 23, 29, 31, 41):
        _prime, degree = by_m[m]
        assert degree > 6, (m, degree)
    assert by_m[41][1] > 2**30
    assert by_m[31][1] > 10**4


if __name__ == "__main__":
    main()
