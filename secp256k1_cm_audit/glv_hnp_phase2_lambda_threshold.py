"""
GLV-HNP Phase 2: lambda/n threshold bisection study.

Goal: empirically determine the critical lambda/n ratio where the Phase 2
GLV-aware attack transitions from success to failure.

From 2026-07-26 log:
  lambda/n = 0.07 (p=2677, n=2647, lam=185): FAIL (LLL and BKZ-40 fail)
  lambda/n = 0.34 (p=524347, n=523969, lam=177902): SUCCESS at m=9

This script:
  1. Sweeps primes p in [2^12, 2^15], finds j=0 prime-order curves.
  2. Groups found curves by lambda/n bucket (width 0.05).
  3. For each bucket, runs Phase 2 attack with 3 seeds and reports success rate.
  4. Identifies the approximate threshold.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Arithmetic helpers
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def isqrt_exact(n):
    if n < 0: return None
    x = int(math.isqrt(n))
    return x if x * x == n else None

def tonelli_shanks(n, p):
    """Compute sqrt(n) mod p. Returns None if n is not a QR mod p."""
    n %= p
    if n == 0: return 0
    if p == 2: return n % 2
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m_ts, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 1, t * t % p
        while tmp != 1:
            tmp = tmp * tmp % p
            i += 1
        b = pow(c, 1 << (m_ts - i - 1), p)
        m_ts, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def miller_rabin(n):
    """Deterministic Miller-Rabin for n < 3.2e18."""
    if n < 2: return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n == p: return True
        if n % p == 0: return False
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

# ---------------------------------------------------------------------------
# Elliptic curve arithmetic (j=0, y^2 = x^3 + b, a=0)
# ---------------------------------------------------------------------------

def ec_add(P, Q, p, b=None):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p, b=None):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def find_point(p, b, rng):
    """Find a point on y^2 = x^3 + b over F_p (not the point at infinity)."""
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n, rng):
    """Find a generator of order n on y^2 = x^3 + b over F_p (n prime)."""
    for _ in range(5000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

# ---------------------------------------------------------------------------
# CM theory for j=0 curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """
    Find (a, b) with a^2 - a*b + b^2 = p, a, b > 0, using the identity
    4p = (2a - b)^2 + 3b^2.
    """
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        sq = isqrt_exact(disc)
        if sq is None:
            continue
        for num in [a + sq, a - sq]:
            if num % 2 == 0:
                bval = num // 2
                if bval > 0 and a * a - a * bval + bval * bval == p:
                    return (a, bval)
    return None

def j0_traces(a, bval):
    """
    Return the 6 Frobenius traces for the 6 sextic twists of y^2=x^3+1 over F_p,
    given Eisenstein norm a^2 - a*b + b^2 = p.
    """
    return [2*a - bval, -(2*a - bval), a + bval, -(a + bval), 2*bval - a, -(2*bval - a)]

def glv_eigenvalue(n):
    """
    Compute the GLV eigenvalue lambda: the smaller root of x^2+x+1 ≡ 0 mod n.
    Returns (lam_small, lam_large) or (None, None) if n is not 1 mod 3.
    """
    if n % 3 != 1:
        return None, None
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    assert (r1 * r1 + r1 + 1) % n == 0, f"GLV check failed r1={r1} n={n}"
    assert (r2 * r2 + r2 + 1) % n == 0, f"GLV check failed r2={r2} n={n}"
    lam_small = min(r1, r2)
    lam_large = max(r1, r2)
    return lam_small, lam_large

def count_points_direct(p, b):
    """Direct point count on y^2 = x^3 + b over F_p. O(p) — only for small p."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x * x * x + b) % p
        if rhs == 0:
            count += 1
        else:
            leg = pow(rhs, (p - 1) // 2, p)
            if leg == 1:
                count += 2
    return count

# ---------------------------------------------------------------------------
# Curve discovery: sweep small primes for j=0 curves with prime group order
# ---------------------------------------------------------------------------

def find_curves_by_lambda_ratio(p_min=4096, p_max=32768, verbose=False):
    """
    Sweep primes p in [p_min, p_max] with p ≡ 1 mod 3.
    For each p, compute the 6 Frobenius traces and check which give prime n.
    For prime n with n ≡ 1 mod 3, compute lambda/n.
    Returns list of (p, b_coeff, n, lam, lam_ratio) sorted by lam_ratio.
    """
    curves = []
    rng = random.Random(42)

    p = p_min | 1  # start odd
    while p <= p_max:
        # Advance to next prime ≡ 1 mod 3
        while not (miller_rabin(p) and p % 3 == 1):
            p += 2
        if p > p_max:
            break

        # Eisenstein decomposition
        dec = eisenstein_decompose(p)
        if dec is None:
            p += 2
            continue
        a_eis, b_eis = dec

        traces = j0_traces(a_eis, b_eis)
        for t in set(traces):
            n = p + 1 - t
            if n < 100 or not miller_rabin(n) or n % 3 != 1:
                continue

            lam_small, lam_large = glv_eigenvalue(n)
            if lam_small is None:
                continue

            # Find curve parameter b such that #E(F_p) = n.
            # Try b = 1, 2, ..., 6 (representatives of the 6 sextic twist classes).
            b_found = None
            for b_try in range(1, p):
                # Quick heuristic: check one random point
                pt = find_point(p, b_try, rng)
                if pt is None:
                    continue
                # Check if n * pt = O (necessary for #E = n or a multiple of n)
                if ec_mul(pt, n, p) is None:
                    # Verify more carefully: count points directly for small p
                    n_actual = count_points_direct(p, b_try)
                    if n_actual == n:
                        b_found = b_try
                        break
                # Only try first 6 candidates
                if b_try >= 6:
                    break

            if b_found is None:
                continue

            ratio = lam_small / n
            curves.append((p, b_found, n, lam_small, ratio))
            if verbose:
                print(f"  Found: p={p} b={b_found} n={n} lam={lam_small} ratio={ratio:.4f}")

        p += 2

    curves.sort(key=lambda x: x[4])
    return curves

# ---------------------------------------------------------------------------
# Phase 2 attack (generic)
# ---------------------------------------------------------------------------

def gen_sigs(d_secret, m, n, lam, k1_bound, k2_bound, G, p, seed=42):
    """Generate m GLV-biased signatures."""
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 100000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0:
            continue
        R = ec_mul(G, k_full, p)
        if R is None:
            continue
        r = R[0] % n
        if r == 0:
            continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0:
            continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Build the (2m+2) x (2m+2) GLV-HNP Phase 2 lattice."""
    S_K1 = n // k1_bound
    S_K2 = n // k2_bound
    S_D = 1
    S_KAN = n
    m = len(sigs)
    dim = 2 * m + 2
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m+1+i][i] = -lam * S_K1
        M[m+1+i][m+1+i] = S_K2
    for i in range(m):
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KAN
    return M, S_D, S_KAN

def try_recover(reduced, m, n, S_D, S_KAN):
    """Scan LLL-reduced rows for d. Returns d or None."""
    dim = 2 * m + 2
    for row in reduced:
        if abs(row[dim-1]) != S_KAN:
            continue
        sign = 1 if row[dim-1] > 0 else -1
        d_cand = sign * row[m] * modinv(S_D, n) % n
        if 1 <= d_cand < n:
            return d_cand
    return None

def run_phase2_attack(p, b_coeff, n, lam, m, d_secret, seed=42, use_bkz=False, bkz_beta=20):
    """
    Run Phase 2 attack on curve y^2=x^3+b_coeff over F_p with group order n, GLV eigenvalue lam.
    Returns True if d_secret is recovered.
    """
    K1_BOUND = 2
    K2_BOUND = math.isqrt(n) + 1

    rng_gen = random.Random(seed + 777)
    G = find_generator(p, b_coeff, n, rng_gen)
    if G is None:
        return False

    sigs = gen_sigs(d_secret, m, n, lam, K1_BOUND, K2_BOUND, G, p, seed)
    if len(sigs) < m:
        return False

    M, S_D, S_KAN = build_lattice(sigs, n, lam, K1_BOUND, K2_BOUND)
    dim = 2 * m + 2
    A = IntegerMatrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            A[i, j] = M[i][j]

    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)

    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d_rec = try_recover(reduced, m, n, S_D, S_KAN)
    return d_rec == d_secret

# ---------------------------------------------------------------------------
# Main experiment: threshold sweep
# ---------------------------------------------------------------------------

def lambda_ratio_bucket(ratio, width=0.05):
    return math.floor(ratio / width) * width

def main():
    print("=" * 70)
    print("GLV-HNP Phase 2 — lambda/n Threshold Bisection Study")
    print("=" * 70)
    print()

    # Step 1: Collect known data points from previous runs.
    #   p=2677, n=2647, lam=185, ratio=0.070: FAIL (from 07-26 log)
    #   p=199 (toy), n=199, lam=106, ratio≈0.533: SUCCESS (from 07-26 log, but lam=106)
    # We also include p=524347, n=523969, lam=177902, ratio≈0.34: SUCCESS.
    # NOTE: For the sweep we use j=0 curves with b=1 or found b.

    # Step 2: Discover curves in [2^12, 2^15].
    print("Discovering j=0 prime-order curves in [4096, 32768]...")
    curves = find_curves_by_lambda_ratio(p_min=4096, p_max=32768, verbose=False)
    print(f"  Found {len(curves)} curves with prime group order and n≡1 mod 3")
    print()

    if not curves:
        print("No curves found — widening search range to [256, 65536]")
        curves = find_curves_by_lambda_ratio(p_min=256, p_max=65536, verbose=False)

    # Step 3: Pick one representative curve per bucket.
    buckets = {}
    for (p, b, n, lam, ratio) in curves:
        bk = round(lambda_ratio_bucket(ratio), 2)
        if bk not in buckets:
            buckets[bk] = (p, b, n, lam, ratio)

    # Also add known data from 07-26 log:
    #   8-bit toy: p=211, n=199, lam=106, ratio=0.53
    #   small-lam fail: p=2677, n=2647, lam=185, ratio=0.07
    known_extra = [
        (211, 2, 199, 106, 106/199),      # 8-bit toy (SUCCESS at m=4)
        (2677, 1, None, None, None),       # small-lam (FAIL) — needs real n/lam
    ]

    # For p=2677, count directly:
    n_2677 = count_points_direct(2677, 1)
    if miller_rabin(n_2677) and n_2677 % 3 == 1:
        lam_s, _ = glv_eigenvalue(n_2677)
        if lam_s is not None:
            ratio_2677 = lam_s / n_2677
            bk = round(lambda_ratio_bucket(ratio_2677), 2)
            if bk not in buckets:
                buckets[bk] = (2677, 1, n_2677, lam_s, ratio_2677)

    # Add 8-bit toy
    n_211 = count_points_direct(211, 2)
    if miller_rabin(n_211) and n_211 % 3 == 1:
        lam_s, _ = glv_eigenvalue(n_211)
        if lam_s is not None:
            ratio_211 = lam_s / n_211
            bk = round(lambda_ratio_bucket(ratio_211), 2)
            if bk not in buckets:
                buckets[bk] = (211, 2, n_211, lam_s, ratio_211)
            else:
                # If there's already an entry in this bucket but from a known curve, prefer known
                pass

    print("Bucket representatives:")
    print(f"  {'bucket':>7}  {'p':>7}  {'n':>7}  {'lam':>7}  {'ratio':>7}")
    for bk in sorted(buckets):
        p, b, n, lam, ratio = buckets[bk]
        print(f"  [{bk:.2f}, {bk+0.05:.2f})  p={p:6d}  n={n:6d}  lam={lam:6d}  ratio={ratio:.4f}")
    print()

    # Step 4: For each bucket, sweep m from info-theoretic threshold to threshold+4.
    #         Test 3 seeds each.
    SEEDS = [42, 1234, 9999]
    NUM_SIGS_MIN = 3
    NUM_SIGS_MAX = 12

    print("=" * 70)
    print("Attack results (K1_BOUND=2, K2_BOUND=ceil(sqrt(n))+1, LLL)")
    print("Format: m=X: wins/3")
    print("=" * 70)
    print()

    results_summary = {}  # ratio -> (m_thresh_found, success_at_thresh)
    rng_main = random.Random(314159)

    for bk in sorted(buckets):
        p, b, n, lam, ratio = buckets[bk]
        K1 = 2
        K2 = math.isqrt(n) + 1
        eff = K1 * K2 / n
        if eff >= 1:
            m_thresh = 3
        else:
            m_thresh = max(3, math.ceil(math.log(n) / math.log(1.0 / eff)))

        print(f"--- lambda/n ≈ {ratio:.3f} (bucket [{bk:.2f}, {bk+0.05:.2f})) ---")
        print(f"    p={p}  n={n}  lam={lam}  K1={K1}  K2={K2}")
        print(f"    eff=K1*K2/n={eff:.4f}  info-theoretic m_thresh≈{m_thresh}")

        bucket_results = {}
        for m in range(NUM_SIGS_MIN, min(NUM_SIGS_MAX, m_thresh + 5) + 1):
            wins = 0
            for seed in SEEDS:
                d_test = rng_main.randint(1, n - 1)
                ok = run_phase2_attack(p, b, n, lam, m, d_test, seed=seed)
                wins += ok
            bucket_results[m] = wins
            tag = "SUCCESS" if wins >= 2 else ("FAIL" if wins == 0 else "PARTIAL")
            print(f"    m={m}: {wins}/{len(SEEDS)}  [{tag}]")

            # Early stop: if 3/3 success at this m, no need to go higher for threshold
            if wins == len(SEEDS) and m >= m_thresh:
                print(f"    -> 3/3 SUCCESS at m={m}. Threshold found.")
                results_summary[ratio] = (m, True)
                break
        else:
            # Check if we ever got 3/3
            max_m = max(bucket_results, key=lambda x: bucket_results[x])
            best = bucket_results[max_m]
            if best < 2:
                print(f"    -> FAIL: never 3/3 (best {best}/3 at m={max_m})")
                results_summary[ratio] = (max_m, False)
            else:
                results_summary[ratio] = (max_m, best >= 2)
        print()

    # Step 5: Summary table
    print("=" * 70)
    print("SUMMARY: lambda/n threshold results")
    print(f"{'ratio':>8}  {'m_thresh':>10}  {'outcome':>10}")
    print("-" * 32)
    for ratio in sorted(results_summary):
        m_used, success = results_summary[ratio]
        print(f"{ratio:8.4f}  {'m='+str(m_used):>10}  {'SUCCESS' if success else 'FAIL':>10}")
    print()

    # Identify threshold
    success_ratios = sorted(r for r, (_, ok) in results_summary.items() if ok)
    fail_ratios = sorted(r for r, (_, ok) in results_summary.items() if not ok)
    if success_ratios and fail_ratios:
        threshold_lo = max(fail_ratios)
        threshold_hi = min(success_ratios)
        if threshold_lo < threshold_hi:
            print(f"THRESHOLD ESTIMATE: lambda/n in ({threshold_lo:.4f}, {threshold_hi:.4f})")
        else:
            print(f"THRESHOLD: overlapping bands — needs finer grained search")
    elif success_ratios:
        print(f"All tested ratios SUCCEED. Min tested: {min(success_ratios):.4f}")
    else:
        print(f"All tested ratios FAIL. Max tested: {max(fail_ratios):.4f}")

    print()
    print("Done.")

if __name__ == "__main__":
    main()
