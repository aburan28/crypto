"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Prior findings (2026-07-26):
  - λ/n ≥ 0.34 → LLL succeeds (8-bit n=199, 12-bit n=2659, 20-bit n=523969)
  - λ/n = 0.07 → LLL fails at all m ≤ 12 (12-bit n=2647)
  - BKZ(β=40) does NOT rescue the small-λ failure

This script searches for j=0 curves at specific λ/n ratio targets,
then runs the Phase 2 LLL attack on each to map the success/failure boundary.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
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

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition for j=0 curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

# ---------------------------------------------------------------------------
# Curve database: find one curve per λ/n ratio bucket
# ---------------------------------------------------------------------------

RATIO_TARGETS = [0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48]
BUCKET_WIDTH = 0.025   # each bucket is ±0.025 around the target

def find_curves_by_ratio(p_max=20000, verbose=True):
    """
    Search j=0 curves over F_p for p up to p_max, collect one curve per ratio bucket.
    Returns dict: ratio_target -> (p, b, n, lam, G).
    """
    buckets = {}  # target -> (p, b, n, lam, G)
    remaining = set(RATIO_TARGETS)

    if verbose:
        print(f"Searching for curves at λ/n targets: {RATIO_TARGETS}")
        print(f"  p range: [7, {p_max}], bucket width ±{BUCKET_WIDTH}")

    p = 7
    while remaining and p < p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n = p + 1 - t
                    if n < 7:
                        continue
                    if not sympy.isprime(n):
                        continue
                    if n % 3 != 1:
                        continue
                    lam, _ = glv_eigenvalue(n)
                    if lam is None:
                        continue
                    ratio = lam / n
                    # Check which buckets this ratio falls into
                    for target in list(remaining):
                        if abs(ratio - target) <= BUCKET_WIDTH:
                            # Find the curve twist
                            G = None
                            for b_try in range(1, min(p, 200)):
                                G = find_generator(p, b_try, n)
                                if G is not None:
                                    b_found = b_try
                                    break
                            if G is not None:
                                buckets[target] = (p, b_found, n, lam, G)
                                remaining.discard(target)
                                if verbose:
                                    print(f"  Bucket {target:.2f}: p={p}, n={n} "
                                          f"({n.bit_length()}b), lam={lam}, "
                                          f"lam/n={ratio:.4f}")
                            break
        p = sympy.nextprime(p)

    if remaining and verbose:
        print(f"  WARNING: no curve found for targets: {sorted(remaining)}")

    return buckets

# ---------------------------------------------------------------------------
# GLV-HNP Phase 2 lattice (column-scaled, from 2026-07-26)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = max(1, n // k1_bound)
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1         # S_D = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_KANNAN

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_lll(curve_params, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def run_bkz(curve_params, m, d_secret, k1_bound, seed, beta):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    BKZ.reduction(A, BKZ.Param(beta))
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def sweep_attack(label, curve_params, k1_bound, m_range, seeds,
                 use_bkz=False, bkz_beta=20, verbose=True):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n

    algo = f"BKZ(β={bkz_beta})" if use_bkz else "LLL"
    if verbose:
        print(f"\n  [{algo}] {label}")
        print(f"    p={p}, n={n} ({n.bit_length()}b), λ={lam}, λ/n={lam/n:.4f}")
        print(f"    K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}")

    first_success_m = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 99991).randint(1, n - 1)
            ok = (run_bkz(curve_params, m, d, k1_bound, seed, bkz_beta)
                  if use_bkz else
                  run_lll(curve_params, m, d, k1_bound, seed))
            wins += ok
        frac = wins / len(seeds)
        marker = " ← FIRST SUCCESS" if (wins == len(seeds) and first_success_m is None) else ""
        if verbose:
            print(f"    m={m}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first_success_m is None:
            first_success_m = m
    return first_success_m  # None = never succeeded

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 68)
print("GLV-aware HNP Phase 2 — λ/n threshold bisection")
print("Known: λ/n≥0.34 → LLL OK; λ/n=0.07 → LLL fail (BKZ(40) also fail)")
print("Goal: find the approximate λ/n threshold where LLL starts failing")
print("=" * 68)

SEEDS = [42, 1234, 9999, 5678, 1111]

# Step 1: Collect one curve per ratio bucket
print("\n--- Step 1: Finding curves at ratio targets ---")
buckets = find_curves_by_ratio(p_max=25000)

print(f"\nFound {len(buckets)}/{len(RATIO_TARGETS)} target buckets")

# Step 2: For each curve, determine the K1_BOUND and sweep
print("\n--- Step 2: LLL attack sweep per ratio ---")
print("(K1_BOUND = max(2, int(0.05*sqrt(n)))")
print()

results = {}
sorted_targets = sorted(buckets.keys())

for target in sorted_targets:
    curve = buckets[target]
    p, b, n, lam, G = curve
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    k2_bound = math.isqrt(n) + 1
    ratio = lam / n

    label = f"λ/n≈{ratio:.3f} (n={n}, {n.bit_length()}b)"
    m_range = range(3, 18)

    first_m = sweep_attack(label, curve, k1_bound, m_range, SEEDS, verbose=True)
    results[target] = (ratio, first_m, n)

# Step 3: For ambiguous/failed cases, try BKZ(40)
print("\n--- Step 3: BKZ(β=40) on failed LLL cases ---")
bkz_results = {}
for target in sorted_targets:
    _, first_m, n = results[target]
    if first_m is None:
        curve = buckets[target]
        p, b, n2, lam, G = curve
        k1_bound = max(2, int(0.05 * math.sqrt(n2)))
        ratio = lam / n2
        label = f"λ/n≈{ratio:.3f} (n={n2}, BKZ(40))"
        first_m_bkz = sweep_attack(label, curve, k1_bound, range(3, 18), SEEDS,
                                   use_bkz=True, bkz_beta=40, verbose=True)
        bkz_results[target] = first_m_bkz
    else:
        bkz_results[target] = None  # LLL already succeeded, no need for BKZ

# Step 4: Summary table
print("\n" + "=" * 68)
print("SUMMARY TABLE")
print(f"{'λ/n target':>12} {'actual λ/n':>12} {'n bits':>7} "
      f"{'LLL m*':>8} {'BKZ(40)':>10} {'verdict':>10}")
print("-" * 68)

for target in sorted_targets:
    ratio, first_m, n = results[target]
    n_bits = n.bit_length()
    lll_str = str(first_m) if first_m is not None else "fail"
    bkz_m = bkz_results.get(target)
    bkz_str = str(bkz_m) if bkz_m is not None else ("n/a" if first_m is not None else "fail")
    verdict = "LLL OK" if first_m is not None else ("BKZ OK" if bkz_m is not None else "FAIL")
    print(f"{target:>12.2f} {ratio:>12.4f} {n_bits:>7} {lll_str:>8} {bkz_str:>10} {verdict:>10}")

print("=" * 68)

# Step 5: Estimate threshold
successes = [(r, fm) for (_, (r, fm, _)) in results.items() if fm is not None]
failures = [(r, fm) for (_, (r, fm, _)) in results.items() if fm is None]

# Also include BKZ rescues
bkz_success_ratios = [results[t][0] for t in sorted_targets
                      if results[t][1] is None and bkz_results.get(t) is not None]

if successes and failures:
    max_fail_ratio = max(r for r, _ in failures)
    min_success_ratio = min(r for r, _ in successes)
    print(f"\nThreshold estimate:")
    print(f"  Max failing λ/n  = {max_fail_ratio:.4f}")
    print(f"  Min succeeding λ/n = {min_success_ratio:.4f}")
    print(f"  Threshold in [{max_fail_ratio:.4f}, {min_success_ratio:.4f}]")
elif not failures:
    min_r = min(r for r, _ in successes)
    print(f"\nAll tested ratios succeed (min λ/n tested = {min_r:.4f})")
    print(f"Threshold is below {min_r:.4f} (or no threshold exists in this range)")
elif not successes:
    max_r = max(r for r, _ in failures)
    print(f"\nAll tested ratios fail (max λ/n tested = {max_r:.4f})")
    print(f"Threshold is above {max_r:.4f}")

if bkz_success_ratios:
    print(f"\n  BKZ(40) rescues LLL failures at λ/n: {sorted(bkz_success_ratios)}")

# Also report known reference points
print("\nKnown reference points (from 2026-07-26 log):")
print("  λ/n=0.53 (8-bit, n=199):      LLL 3/3 at m=4")
print("  λ/n=0.66 (12-bit, n=2659):    LLL 3/3 at m=7")
print("  λ/n=0.34 (20-bit, n=523969):  LLL 3/3 at m=9")
print("  λ/n=0.07 (12-bit, n=2647):    LLL fail, BKZ(40) fail")

print("\nDone.")
