"""What Cheon's algorithm could use on standard curves: the divisors of
n − 1 and n + 1 (n the prime subgroup order) that a q-SDH-style leak of
`[α^i]G`, i ≤ q, would make available.

For each curve: trial division of n ∓ 1 to 2^22, then Pollard rho with a
bounded budget on the cofactor.  Whatever is left unsplit is reported as
`unsplit` and excluded, so every "largest usable d ≤ q" below is a lower
bound on the true value (a larger prime factor could only add divisors
above the budgets shown, and only if it is itself ≤ q).

    python3 curve_divisors.py > results/curve_divisors.md
"""

from __future__ import annotations

import math
import random
from typing import Dict, List, Tuple

from ec import divisors, is_probable_prime

CURVES: Dict[str, int] = {
    # prime subgroup orders
    "secp256k1 (n)": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
    "P-256 (n)": 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    "P-384 (n)": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973,
    "Curve25519 (ℓ)": 2**252 + 27742317777372353535851937790883648493,
    "BN254 (r)": 21888242871839275222246405745257275088548364400416034343698204186575808495617,
    "BLS12-381 (r)": 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001,
    "ECC2K-130 (r)": 680564733841876926932320129493409985129,
}

BUDGET = 3_000_000  # Pollard-rho iterations per cofactor

# Factors beyond the rho budget that are known from the literature; each
# list is verified below (every entry prime, product divides the cofactor)
# before it is used, so a wrong entry is reported, not silently trusted.
KNOWN_FACTORS: Dict[Tuple[str, str], List[int]] = {
    ("secp256k1 (n)", "n−1"): [107361793816595537, 174723607534414371449,
                               341948486974166000522343609283189],
}


def trial(n: int, bound: int) -> Tuple[Dict[int, int], int]:
    fac: Dict[int, int] = {}
    for q in range(2, bound):
        if q * q > n:
            break
        while n % q == 0:
            fac[q] = fac.get(q, 0) + 1
            n //= q
    return fac, n


def rho_split(n: int, budget: int, rng: random.Random) -> int | None:
    if n % 2 == 0:
        return 2
    for _ in range(6):
        c = rng.randrange(1, n)
        x = y = rng.randrange(0, n)
        g = 1
        it = 0
        while g == 1 and it < budget // 6:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            g = math.gcd(abs(x - y), n)
            it += 1
        if 1 < g < n:
            return g
    return None


def factor_bounded(n: int, rng: random.Random, known: List[int] = ()) -> Tuple[Dict[int, int], List[int]]:
    fac, rest = trial(n, 1 << 22)
    for q in known:
        if is_probable_prime(q) and rest % q == 0:
            while rest % q == 0:
                fac[q] = fac.get(q, 0) + 1
                rest //= q
        else:
            print(f"<!-- known factor {q} rejected: not prime or does not divide -->")
    unsplit: List[int] = []
    stack = [rest] if rest > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_probable_prime(m):
            fac[m] = fac.get(m, 0) + 1
            continue
        g = rho_split(m, BUDGET, rng)
        if g is None:
            unsplit.append(m)
        else:
            stack += [g, m // g]
    return fac, unsplit


def largest_divisor_below(fac: Dict[int, int], bound: int) -> int:
    best = 1
    for d in divisors(fac):
        if d <= bound and d > best:
            best = d
    return best


def main():
    rng = random.Random(1)
    print("# Divisors of n ± 1 on standard curves, and the Cheon gain a q-SDH leak would buy\n")
    print("Trial division to 2^22, then Pollard rho with a bounded budget; `unsplit` cofactors are")
    print("excluded, so each `largest d ≤ q` is a lower bound.  Gain is √d for the p−1 case (needs")
    print("`[α^d]G`, d | n−1) and √d/… for the p+1 case (needs 2d inputs, d | n+1, cost √(n/d)+d).\n")
    print("| curve | side | known factorisation | unsplit | largest d ≤ 2^20 | ≤ 2^32 | ≤ 2^48 | ≤ 2^64 | √d at 2^32 | √d at 2^64 |")
    print("|:--|:--|:--|:--|---:|---:|---:|---:|---:|---:|")
    for name, n in CURVES.items():
        for side, m in (("n−1", n - 1), ("n+1", n + 1)):
            fac, unsplit = factor_bounded(m, rng, KNOWN_FACTORS.get((name, side), []))
            fs = " · ".join(f"{q}" + (f"^{e}" if e > 1 else "") for q, e in sorted(fac.items()))
            us = ", ".join(f"c{len(str(u))}" for u in unsplit) or "—"
            ds = [largest_divisor_below(fac, 1 << k) for k in (20, 32, 48, 64)]
            print(f"| {name} | {side} | {fs} | {us} | {ds[0]:,} | {ds[1]:,} | {ds[2]:,} | {ds[3]:,} | "
                  f"2^{math.log2(ds[1]) / 2:.1f} | 2^{math.log2(ds[3]) / 2:.1f} |")


if __name__ == "__main__":
    main()
