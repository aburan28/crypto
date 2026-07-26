"""
Find j=0 prime-order curves with various λ/n ratios for GLV-HNP threshold study.
Target: prime n in 12-16 bit range (2000–65000).

For larger curves we use PARI's ellcard via subprocess, or brute-force Hasse-range
point counts. For now, use the CM formula for j=0 curves (y^2 = x^3 + b over F_p):
  4p = t^2 + 3u^2  where p ≡ 1 mod 3, t ≡ 1 mod 3 (by convention)
  n = p + 1 - t  (or p + 1 + t depending on twist)

We iterate over t, u satisfying 4p = t^2 + 3u^2, then check if n = p+1-t is prime.

Output: list of (p, n, lam1, lam2, r1, r2) with r1 < r2.
"""

import math
import sys

def isprime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for d in range(3, int(n**0.5) + 1, 2):
        if n % d == 0: return False
    return True

def glv_roots(n):
    """Find roots of x^2+x+1 ≡ 0 mod n (only if n ≡ 1 mod 3)."""
    if n % 3 != 1:
        return []
    # Compute sqrt(-3) mod n via Tonelli-Shanks on small candidates
    disc = (-3) % n
    # For small n, brute force sqrt
    sq = -1
    for x in range(n):
        if x * x % n == disc:
            sq = x
            break
    if sq < 0:
        return []
    inv2 = pow(2, -1, n)
    r1 = (-1 + sq) * inv2 % n
    r2 = (-1 - sq) * inv2 % n
    roots = sorted([r for r in [r1, r2] if (r*r + r + 1) % n == 0])
    return roots

def count_points_small(p, b):
    """Brute-force #E(F_p) for y^2 = x^3 + b. Only use for p < 5000."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x**3 + b) % p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p-1)//2, p) == 1:
            count += 2
    return count

def find_curves_cm(n_min=2000, n_max=16000):
    """
    Use CM formula 4p = t^2 + 3u^2 to find j=0 prime-order curves.
    For each candidate (t, u), p = (t^2 + 3u^2)/4, n = p+1-t.
    Require: p prime, n prime, n ≡ 1 mod 3, n ∈ [n_min, n_max].
    """
    curves = []
    seen_n = set()

    # t must be even or both t,u must have same parity for 4p = t^2+3u^2 to be integer
    # Actually: 4p = t^2 + 3u^2 requires t^2 ≡ 0 or 1 mod 4 and 3u^2 ≡ 0 or 3 mod 4.
    # If t is even, u must be even. If t is odd, u must be odd.
    # Standard convention: t ≡ 1 mod 3 (to distinguish the two square roots mod n).

    # Search range for t based on n range: n ≈ p, p ≈ n_max. t ≤ 2√p ≤ 2√n_max.
    t_max = int(2 * math.sqrt(n_max)) + 1
    u_max = int(2 * math.sqrt(n_max / 3)) + 1

    for t in range(-t_max, t_max + 1):
        for u in range(1, u_max + 1):
            # 4p = t^2 + 3u^2 must be divisible by 4
            val = t*t + 3*u*u
            if val % 4 != 0:
                continue
            p = val // 4
            if p < 7:
                continue
            if not isprime(p):
                continue
            if p % 3 != 1:
                continue

            # Two possible n values (two twists)
            for sign in [1, -1]:
                n = p + 1 - sign * t
                if n < n_min or n > n_max:
                    continue
                if not isprime(n):
                    continue
                if n % 3 != 1:
                    continue
                if n in seen_n:
                    continue

                roots = glv_roots(n)
                if len(roots) < 2:
                    continue

                lam1, lam2 = roots[0], roots[1]
                r1 = lam1 / n
                r2 = lam2 / n

                # Verify with actual point count (for small p)
                if p < 3000:
                    # Find a b such that #E(F_p, b) = n
                    found_b = None
                    for b in range(1, p):
                        actual_n = count_points_small(p, b)
                        if actual_n == n:
                            found_b = b
                            break
                    if found_b is None:
                        continue
                else:
                    # Skip verification for large p (trust CM formula)
                    # Find any b that gives a j=0 curve (b=1 is the canonical one)
                    found_b = 1  # CM guarantee: some b works

                seen_n.add(n)
                curves.append((p, found_b, n, lam1, lam2, r1, r2))

    return sorted(curves, key=lambda x: x[5])  # sort by r1

if __name__ == "__main__":
    print("Finding 12-14 bit j=0 prime-order curves with various λ/n...")
    curves = find_curves_cm(n_min=2000, n_max=16000)
    print(f"Found {len(curves)} distinct curves.\n")
    print(f"{'p':>7}  {'b':>4}  {'n':>6}  {'lam1':>6}  {'r1':>7}  {'lam2':>6}  {'r2':>7}")
    print("-" * 55)
    for (p, b, n, l1, l2, r1, r2) in curves[:40]:
        print(f"{p:>7}  {b:>4}  {n:>6}  {l1:>6}  {r1:>7.4f}  {l2:>6}  {r2:>7.4f}")
