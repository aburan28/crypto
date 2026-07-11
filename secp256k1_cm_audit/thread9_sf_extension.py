"""
thread9_sf_extension.py
Thread 9: Extended squarefree discriminant sweep for naive-cover Jacobians.

Goals:
  1. Extend to all primes p≡1 mod 6 up to p=211.
  2. Track sf(D_new) = sf(a1²-4(a2-2p)) for each pair of sextic twists.
  3. Test conjecture: for pair (0,1) [1-indexed (1,2)], is sf=73 structural?
  4. Detect if a2-2p takes constant value across multiple primes for any fixed pair.

Run: python3 thread9_sf_extension.py
"""

import math
from itertools import combinations


def legendre(a, p):
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else r


def is_square_mod(a, p):
    return legendre(a, p) == 1


def curve_order_j0(p, b):
    count = 1
    for x in range(p):
        rhs = (pow(x, 3, p) + b) % p
        if rhs == 0:
            count += 1
        elif is_square_mod(rhs, p):
            count += 2
    return count


def genus2_point_count(p, b_i, b_j):
    """Count affine points on y²=(x³+b_i)(x³+b_j) over F_p, plus 2 points at infinity."""
    count = 2
    for x in range(p):
        xi3 = pow(x, 3, p)
        fx = ((xi3 + b_i) % p * (xi3 + b_j)) % p
        if fx == 0:
            count += 1
        elif is_square_mod(fx, p):
            count += 2
    return count


def genus2_point_count_f_p2(p, b_i, b_j):
    """Count F_{p^2} points on y²=(x³+b_i)(x³+b_j)."""
    ns = next(r for r in range(2, p) if not is_square_mod(r, p))

    def fq2_mul(x1, x2):
        a, b = x1
        c, d = x2
        return ((a * c + b * d * ns) % p, (a * d + b * c) % p)

    def fq2_pow(x, n):
        r = (1, 0)
        base = x
        while n > 0:
            if n & 1:
                r = fq2_mul(r, base)
            base = fq2_mul(base, base)
            n >>= 1
        return r

    def is_square_fq2(a):
        # a is in F_{p^2}; check if it's a square there
        # #F_{p^2}* = p^2-1; a is a square iff a^((p^2-1)/2) = 1
        exp = (p * p - 1) // 2
        r = fq2_pow(a, exp)
        return r == (1, 0)

    count = 2  # two points at infinity (same as over F_p for genus 2)
    for a in range(p):
        for b in range(p):
            x = (a, b)
            # compute x^3 over F_{p^2}
            x2 = fq2_mul(x, x)
            x3 = fq2_mul(x2, x)
            # add b_i and b_j (elements of F_p ⊂ F_{p^2})
            lft = ((x3[0] + b_i) % p, x3[1])
            rgt = ((x3[0] + b_j) % p, x3[1])
            fx = fq2_mul(lft, rgt)
            if fx == (0, 0):
                count += 1
            elif is_square_fq2(fx):
                count += 2
    return count


def charpoly_genus2(p, b_i, b_j):
    """Compute Weil polynomial coefficients (a1, a2) for y²=(x³+b_i)(x³+b_j)/F_p."""
    N1 = genus2_point_count(p, b_i, b_j)
    N2 = genus2_point_count_f_p2(p, b_i, b_j)
    a1 = p + 1 - N1
    sum_sq = p * p + 1 - N2   # = a1² - 2*a2
    assert (a1 * a1 - sum_sq) % 2 == 0, f"Parity error: a1={a1}, sum_sq={sum_sq}"
    a2 = (a1 * a1 - sum_sq) // 2
    return a1, a2


def squarefree_part(n):
    if n == 0:
        return 0, 0
    n = abs(n)
    sf = 1
    sq = 1
    d = 2
    tmp = n
    while d * d <= tmp:
        cnt = 0
        while tmp % d == 0:
            tmp //= d
            cnt += 1
        if cnt % 2 == 1:
            sf *= d
        sq *= d ** (cnt // 2 * 2)
        d += 1 if d == 2 else 2
    if tmp > 1:
        sf *= tmp
    return sf, n // sf


def is_perfect_square(n):
    if n < 0:
        return False
    s = int(math.isqrt(n))
    return s * s == n


def primitive_root(p):
    phi = p - 1
    factors = set()
    n = phi
    for d in range(2, int(math.isqrt(n)) + 1):
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
    if n > 1:
        factors.add(n)
    for g in range(2, p):
        if all(pow(g, phi // f, p) != 1 for f in factors):
            return g
    return None


def primes_1mod6_up_to(limit):
    """Return primes p≡1 mod 6 up to limit."""
    def is_prime(n):
        if n < 2:
            return False
        if n == 2:
            return True
        if n % 2 == 0:
            return False
        for d in range(3, int(math.isqrt(n)) + 1, 2):
            if n % d == 0:
                return False
        return True
    return [p for p in range(7, limit + 1) if p % 6 == 1 and is_prime(p)]


def analyze_prime(p):
    """Return list of (pair_idx, a1, a2, D, sf, category) for all 15 pairs."""
    g = primitive_root(p)
    twists = [pow(g, k, p) for k in range(6)]
    results = []
    for i, j in combinations(range(6), 2):
        bi, bj = twists[i], twists[j]
        a1, a2 = charpoly_genus2(p, bi, bj)
        D = a1 * a1 - 4 * (a2 - 2 * p)
        if D == 0:
            cat = 'SPLIT-REP'
            sf = 0
        elif is_perfect_square(D):
            cat = 'SPLIT'
            sf = 0
        elif D < 0:
            cat = 'IRRED'
            sf, _ = squarefree_part(abs(D))
            sf = -sf  # negative to indicate irreducible
        else:
            cat = 'NONSPLIT-Q'
            sf, _ = squarefree_part(D)
        results.append((i, j, a1, a2, D, sf, cat))
    return results


def main():
    primes = primes_1mod6_up_to(211)
    print(f"Primes p≡1 mod 6, up to 211: {primes}")
    print()

    # Track sf values by pair ratio class
    # Pairs with ratio g^1: (0,1), (3,4) — same D by symmetry
    # Pairs with ratio g^3: (0,3), (1,4) — same D by symmetry
    # Pairs with ratio g^4: (0,4), (1,5)? — TBD
    # Pairs with ratio g^2: (0,2), (1,3) — TBD
    # Pairs with ratio g^5: (0,5), (2,3)? — TBD

    pair_01_data = []  # (1,2) in 1-indexed, i.e., (0,1) in 0-indexed
    pair_03_data = []  # (1,4) in 1-indexed

    print("=" * 90)
    print(f"{'p':>5} | {'pair':>6} | {'a1':>5} | {'a2':>7} | {'D':>8} | {'sf':>6} | category")
    print("=" * 90)

    all_nonsplit = []
    all_irred = []

    for p in primes:
        results = analyze_prime(p)
        for (i, j, a1, a2, D, sf, cat) in results:
            if cat in ('NONSPLIT-Q', 'IRRED'):
                pair_label = f"({i+1},{j+1})"
                print(f"{p:>5} | {pair_label:>6} | {a1:>5} | {a2:>7} | {D:>8} | {abs(sf):>6} | {cat}")
                if cat == 'NONSPLIT-Q':
                    all_nonsplit.append((p, i, j, a1, a2, D, sf))
                    if (i, j) == (0, 1):
                        pair_01_data.append((p, a1, a2, D, sf))
                    if (i, j) == (0, 3):
                        pair_03_data.append((p, a1, a2, D, sf))
                else:
                    all_irred.append((p, i, j, a1, a2, D, sf))

    print()
    print("=" * 90)
    print(f"Total NONSPLIT-Q cases: {len(all_nonsplit)}")
    print(f"Total IRRED cases:      {len(all_irred)}")
    print()

    print("Pair (0,1) = 1-indexed (1,2): sf(D_new) progression")
    print(f"  {'p':>5} | {'a1':>5} | {'a2':>7} | {'D':>8} | {'2p-a2':>7} | sf")
    for (p, a1, a2, D, sf) in pair_01_data:
        print(f"  {p:>5} | {a1:>5} | {a2:>7} | {D:>8} | {2*p-a2:>7} | {sf}")
    print()

    print("Pair (0,3) = 1-indexed (1,4): sf(D_new) progression")
    print(f"  {'p':>5} | {'a1':>5} | {'a2':>7} | {'D':>8} | sf")
    for (p, a1, a2, D, sf) in pair_03_data:
        print(f"  {p:>5} | {a1:>5} | {a2:>7} | {D:>8} | {sf}")
    print()

    # Check if any sf value recurs across multiple primes
    from collections import defaultdict
    sf_to_primes = defaultdict(list)
    for (p, i, j, a1, a2, D, sf) in all_nonsplit:
        sf_to_primes[sf].append((p, i, j))
    print("Recurring squarefree values (sf appearing at ≥2 primes):")
    for sf, pts in sorted(sf_to_primes.items()):
        if len(pts) >= 2:
            print(f"  sf={sf}: {pts}")
    print()

    # Look for constant a2-2p pattern
    val_to_cases = defaultdict(list)
    for (p, i, j, a1, a2, D, sf) in all_nonsplit:
        val_to_cases[a2 - 2*p].append((p, i+1, j+1, a1, a2, sf))
    print("Constant a2-2p values (appearing at ≥2 primes):")
    for val, cases in sorted(val_to_cases.items()):
        if len(cases) >= 2:
            print(f"  a2-2p={val}: {[(c[0], f'({c[1]},{c[2]})', f'a1={c[3]}', f'sf={c[5]}') for c in cases]}")

    print()
    print("Structural check: do pairs (0,1)=(1,2) and (3,4)=(4,5) always give same D?")
    for p in primes:
        results = analyze_prime(p)
        by_pair = {(r[0],r[1]): r for r in results}
        for (a,b),(c,d) in [((0,1),(3,4)), ((0,3),(1,4))]:
            if (a,b) in by_pair and (c,d) in by_pair:
                D_ab = by_pair[(a,b)][4]
                D_cd = by_pair[(c,d)][4]
                if D_ab != D_cd:
                    print(f"  MISMATCH at p={p}: ({a},{b}) D={D_ab}, ({c},{d}) D={D_cd}")
    print("  Symmetry check complete (no mismatch = symmetry holds).")


if __name__ == '__main__':
    main()
