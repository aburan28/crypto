"""
scaling.py -- measure how the cost of finding ONE vanishing minor (= one ECDLP
relation) scales with the group order n, for the Abdullah-Mahalanobis-Mallick
minors method.

The method finds a relation when a k-subset of the points R_i = a_iP + b_iQ
sums to O.  The density of such subsets is

        rho_k = Pr[ k random curve points sum to O ] ~ 1/n,

independent of k for k >= 2 (the sum of k>=2 independent ~uniform points is
~uniform on the group, hence O with probability 1/n).  So the number of minors
that must be *examined* before a vanishing one appears is ~ n.

This script measures E[subsets-examined-to-first-relation] across a wide range
of n and checks whether it tracks c*n (exponential in log n, i.e. NO
subexponential gain) or grows sub-exponentially.
"""

import random
import time
from ecmin import Curve, curve_order, legendre, _tonelli
from verify_correspondence import is_prime, random_prime, find_point


def make_curve(target_bits, tries=400):
    lo, hi = 1 << target_bits, 1 << (target_bits + 1)
    for _ in range(tries):
        p = random_prime(lo, hi)
        for _ in range(30):
            a, b = random.randrange(p), random.randrange(p)
            if (4 * a * a * a + 27 * b * b) % p == 0:
                continue
            n = curve_order(a, b, p)
            if is_prime(n):
                C = Curve(a, b, p)
                P = find_point(C, n)
                if P is not None:
                    return p, a, b, n, C, P
    raise RuntimeError("no prime-order curve found")


def build_pool(C, P, Q, n, size):
    pool, coeffs = [], []
    while len(pool) < size:
        aa, bb = random.randrange(n), random.randrange(n)
        R = C.add(C.mul(aa, P), C.mul(bb, Q))
        if R is not None:
            pool.append(R)
            coeffs.append((aa, bb))
    return pool, coeffs


def cost_to_first_relation(C, pool, coeffs, n, k, cap):
    """Sample random k-subsets of the pool until one sums to O; return the count
    of subsets examined (capped)."""
    M = len(pool)
    for cnt in range(1, cap + 1):
        idx = random.sample(range(M), k)
        S = None
        for i in idx:
            S = C.add(S, pool[i])
        if S is None:
            # genuine relation; sanity-recover the discrete log
            A = sum(coeffs[i][0] for i in idx) % n
            B = sum(coeffs[i][1] for i in idx) % n
            return cnt
    return None  # hit cap without a relation


def main():
    random.seed(7)
    k = 4
    bit_sizes = [8, 10, 12, 14, 16, 18, 20]
    print(f"# minors method: cost to first vanishing minor (relation), k={k}")
    print(f"{'bits':>4} {'n':>10} {'reps':>5} {'E[subsets]':>12} "
          f"{'E[subsets]/n':>13} {'wall(s)':>8}")
    rows = []
    for bits in bit_sizes:
        t0 = time.time()
        p, a, b, n, C, P = make_curve(bits)
        m = random.randrange(2, n)
        Q = C.mul(m, P)
        pool, coeffs = build_pool(C, P, Q, n, size=min(3000, max(200, 4 * n // 1)))
        # choose reps so total work stays bounded; cap search at ~12n subsets
        cap = max(2000, 16 * n)
        reps = max(2, min(300, 2_000_000 // max(1, n)))
        counts = []
        for _ in range(reps):
            c = cost_to_first_relation(C, pool, coeffs, n, k, cap)
            if c is not None:
                counts.append(c)
        if counts:
            avg = sum(counts) / len(counts)
            rows.append((bits, n, len(counts), avg, avg / n))
            print(f"{bits:>4} {n:>10} {len(counts):>5} {avg:>12.1f} "
                  f"{avg / n:>13.3f} {time.time() - t0:>8.1f}")
        else:
            print(f"{bits:>4} {n:>10}   no relations within cap")
    # summary: is E[subsets]/n roughly constant?
    if len(rows) >= 2:
        ratios = [r[4] for r in rows]
        print(f"\n# E[subsets]/n over the range: "
              f"min={min(ratios):.3f} max={max(ratios):.3f} "
              f"mean={sum(ratios)/len(ratios):.3f}")
        print("# A constant ratio across a >1000x range of n means cost = "
              "Theta(n):")
        print("# the natural minor search is NOT subexponential.")


if __name__ == "__main__":
    main()
