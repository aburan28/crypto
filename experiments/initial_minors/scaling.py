"""
scaling.py -- measure how the cost of finding ONE vanishing minor (= one ECDLP
relation) scales with the group order n, for the Abdullah-Mahalanobis-Mallick
minors method.

The method builds M[i][j] = phi_j(R_i) for random points R_i = a_iP + b_iQ and
finds a relation when a leading principal ("initial") k x k minor of M vanishes,
i.e. when the first k rows sum to O.  The density of such prefixes is

        rho_k = Pr[ k random curve points sum to O ] ~ 1/n,

independent of k for k >= 2 (the sum of k>=2 independent ~uniform points is
~uniform on the group, hence O with probability 1/n).  So the number of initial
minors that must be *examined* before a vanishing one appears is ~ n.

This script draws random N-row matrices, reads their initial minors off one
pivoting-free Gaussian elimination (ecmin.leading_minor_zero_pivots), and
measures E[initial-minors-examined-to-first-relation] across a wide range of n,
checking whether it tracks c*n (exponential in log n, i.e. NO subexponential
gain) or grows sub-exponentially.  Every zero pivot is cross-checked against the
group law and the recovered discrete log is checked against the planted secret.
"""

import random
import time
from ecmin import Curve, curve_order, build_matrix, leading_minor_zero_pivots
from verify_correspondence import is_prime, random_prime, find_point, sum_points


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
    """Distinct affine points a P + b Q.  Distinctness matters: a repeated row
    makes a minor vanish without any relation, so it must not be counted."""
    pool, coeffs, seen = [], [], set()
    while len(pool) < size:
        aa, bb = random.randrange(n), random.randrange(n)
        R = C.add(C.mul(aa, P), C.mul(bb, Q))
        if R is not None and R not in seen:
            seen.add(R)
            pool.append(R)
            coeffs.append((aa, bb))
    return pool, coeffs


def cost_to_first_relation(C, pool, coeffs, rows, m_secret, n, N, cap):
    """Draw random N-row matrices M[i][j] = phi_j(R_i) from the pool and scan
    their leading principal minors (LU pivots, no row swaps) until one vanishes.
    Returns (initial minors examined, discrete log recovered correctly), or None
    if the cap was hit.  Minors of size 2..N are counted: the 1x1 minor is
    phi_1 = 1 identically and can never vanish."""
    M = len(pool)
    examined = 0
    while examined < cap:
        idx = random.sample(range(M), N)
        k = leading_minor_zero_pivots([rows[i] for i in idx], C.p)
        if k is None:
            examined += N - 1
            continue
        examined += k - 1
        S = idx[:k]
        if sum_points(C, [pool[i] for i in S]) is not None:
            raise AssertionError("zero pivot but prefix does not sum to O")
        A = sum(coeffs[i][0] for i in S) % n
        B = sum(coeffs[i][1] for i in S) % n
        recovered = (-A * pow(B, n - 2, n)) % n if B else None
        return examined, recovered == m_secret
    return None  # hit cap without a relation


def main():
    random.seed(7)
    N = 16
    bit_sizes = [8, 10, 12, 14, 16, 18, 20]
    print(f"# minors method: cost to first vanishing initial minor (relation), "
          f"N={N} rows per matrix")
    print(f"{'bits':>4} {'n':>10} {'reps':>5} {'ok':>5} {'E[minors]':>12} "
          f"{'E[minors]/n':>13} {'wall(s)':>8}")
    rows = []
    for bits in bit_sizes:
        t0 = time.time()
        p, a, b, n, C, P = make_curve(bits)
        m = random.randrange(2, n)
        Q = C.mul(m, P)
        pool, coeffs = build_pool(C, P, Q, n, size=min(3000, n // 2))
        mrows = build_matrix(pool, p, ncols=N)
        # choose reps so total work stays bounded; cap search at 16n minors
        cap = max(2000, 16 * n)
        reps = max(2, min(300, 2_000_000 // max(1, n)))
        counts, ok = [], 0
        for _ in range(reps):
            c = cost_to_first_relation(C, pool, coeffs, mrows, m, n, N, cap)
            if c is not None:
                counts.append(c[0])
                ok += c[1]
        if counts:
            avg = sum(counts) / len(counts)
            rows.append((bits, n, len(counts), avg, avg / n))
            print(f"{bits:>4} {n:>10} {len(counts):>5} {ok:>5} {avg:>12.1f} "
                  f"{avg / n:>13.3f} {time.time() - t0:>8.1f}")
        else:
            print(f"{bits:>4} {n:>10}   no relations within cap")
    # summary: is E[minors]/n roughly constant?
    if len(rows) >= 2:
        ratios = [r[4] for r in rows]
        print(f"\n# E[minors]/n over the range: "
              f"min={min(ratios):.3f} max={max(ratios):.3f} "
              f"mean={sum(ratios)/len(ratios):.3f}")
        print("# A constant ratio across a >1000x range of n means cost = "
              "Theta(n):")
        print("# the initial-minor search is NOT subexponential.")


if __name__ == "__main__":
    main()
