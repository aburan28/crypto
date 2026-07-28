"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Goal: find the boundary ratio λ/n below which LLL fails to recover d.
Known data points (from 2026-07-26):
  lam/n ≈ 0.07  →  FAIL (p=2677, even BKZ-40 fails)
  lam/n ≈ 0.34  →  SUCCESS (p≈524347)

Strategy: scan 10–13-bit j=0 curves (y²=x³+b), compute lam/n for each,
group by ratio, and find min m for 3/3 LLL recovery per bucket.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys
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
    if k < 0: return ec_mul((P[0], (-P[1]) % p), -k, p)
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
    mc, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mc - i - 1), p)
        mc, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(1234)
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
# CM theory for j=0 curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -(2*a - b), a + b, -(a + b), 2*b - a, -(2*b - a)]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

# ---------------------------------------------------------------------------
# Curve collection
# ---------------------------------------------------------------------------

def collect_curves(bit_min=9, bit_max=13, max_per_bucket=2):
    """
    Scan j=0 prime-order curves in [2^bit_min, 2^bit_max].
    Returns dict: ratio_bucket -> list of (p, b, n, lam, G, ratio).
    Buckets at 0.025 width: 0.000, 0.025, ..., 0.500.
    """
    buckets = {}

    p = sympy.nextprime(2**bit_min - 1)
    found = 0
    while p < 2**bit_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                traces = j0_traces(a_e, b_e)
                for t in traces:
                    n_cand = p + 1 - t
                    if n_cand < 4: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    # Only keep ratio in (0.04, 0.50] — the interesting range
                    if ratio < 0.04 or ratio > 0.50: continue
                    bucket = round(int(ratio / 0.025) * 0.025, 3)
                    if bucket not in buckets:
                        buckets[bucket] = []
                    if len(buckets[bucket]) >= max_per_bucket:
                        continue
                    # Find a curve twist with this group order
                    found_b = None
                    for b_try in range(1, p):
                        G = find_generator(p, b_try, n_cand)
                        if G is not None:
                            found_b = b_try
                            break
                    if found_b is None: continue
                    G = find_generator(p, found_b, n_cand)
                    if G is None: continue
                    buckets[bucket].append((p, found_b, n_cand, lam, G, ratio))
                    found += 1

        p = sympy.nextprime(p)

    return buckets

# ---------------------------------------------------------------------------
# Signature generation
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

# ---------------------------------------------------------------------------
# Lattice
# ---------------------------------------------------------------------------

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
    M[m][m] = 1  # S_D = 1
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
        if d_cand == d_secret: return True
    return False

def run_lll(p, b, n, lam, G, m, d_secret, k1_bound, seed=42):
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

SEEDS = [42, 1234, 9999]

def min_m_for_success(p, b, n, lam, G, k1_bound, m_max=18, required=3):
    """Return minimum m such that LLL recovers d for `required` out of len(SEEDS) seeds."""
    for m in range(3, m_max + 1):
        wins = 0
        for seed in SEEDS:
            d = random.Random(seed * 31337).randint(1, n - 1)
            if run_lll(p, b, n, lam, G, m, d, k1_bound, seed):
                wins += 1
        if wins >= required:
            return m, wins
    return None, 0  # never 3/3 in range

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 72)
print("GLV-HNP Phase 2: λ/n Threshold Bisection Study")
print("=" * 72)

# --- Step 1: Collect curves grouped by lam/n bucket ----------------------
print("\nStep 1: Scanning 9–13-bit j=0 prime-order curves...")
sys.stdout.flush()

buckets = collect_curves(bit_min=9, bit_max=13, max_per_bucket=2)

all_ratios = sorted(buckets.keys())
total_curves = sum(len(v) for v in buckets.values())
print(f"  Found {total_curves} curves across {len(all_ratios)} ratio buckets.")
print(f"  Ratio range: [{min(all_ratios):.3f}, {max(all_ratios):.3f}]")
print(f"  Buckets with ≥1 curve: {[f'{r:.3f}' for r in all_ratios]}")
sys.stdout.flush()

# --- Step 2: For each bucket, test LLL ------------------------------------
print("\nStep 2: Testing LLL recovery per bucket (m_max=18, 3/3 threshold)...")
print(f"  {'Bucket':>8}  {'p':>6}  {'n':>6}  {'lam/n':>7}  {'bits':>4}  "
      f"{'K1':>5}  {'m*':>4}  {'status':}")
print("  " + "-" * 65)
sys.stdout.flush()

threshold_data = []  # (ratio, min_m_or_None)

for bucket_ratio in all_ratios:
    curves = buckets[bucket_ratio]
    # Use first curve in bucket
    p, b, n, lam, G, ratio = curves[0]
    n_bits = n.bit_length()
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n

    min_m, wins = min_m_for_success(p, b, n, lam, G, k1_bound, m_max=18, required=3)
    status = f"m*={min_m} (3/{len(SEEDS)})" if min_m else f"FAIL (best {wins}/{len(SEEDS)} up to m=18)"
    print(f"  {bucket_ratio:8.3f}  {p:6d}  {n:6d}  {ratio:7.4f}  {n_bits:4d}  "
          f"{k1_bound:5d}  {min_m if min_m else '--':>4}  {status}")
    sys.stdout.flush()
    threshold_data.append((ratio, min_m))

# --- Step 3: Identify threshold -------------------------------------------
print("\n" + "=" * 72)
print("Step 3: Threshold analysis")
print("=" * 72)

success_ratios = [(r, m) for r, m in threshold_data if m is not None]
fail_ratios    = [(r, m) for r, m in threshold_data if m is None]

if success_ratios and fail_ratios:
    max_fail_ratio  = max(r for r, _ in fail_ratios)
    min_succ_ratio  = min(r for r, _ in success_ratios)
    print(f"\n  Highest FAIL ratio:  {max_fail_ratio:.4f}")
    print(f"  Lowest SUCCESS ratio: {min_succ_ratio:.4f}")
    print(f"  Threshold bracket:   ({max_fail_ratio:.3f}, {min_succ_ratio:.3f})")
    mid = (max_fail_ratio + min_succ_ratio) / 2
    print(f"  Midpoint estimate:   {mid:.3f}")
elif not fail_ratios:
    print(f"\n  All buckets SUCCEEDED — threshold is below {min(r for r, _ in success_ratios):.3f}")
elif not success_ratios:
    print(f"\n  All buckets FAILED — threshold is above {max(r for r, _ in fail_ratios):.3f}")

# --- Step 4: Print summary table ------------------------------------------
print("\n" + "=" * 72)
print("Summary table (lam/n → min m for 3/3 LLL success)")
print("=" * 72)
print(f"  {'lam/n':>8}  {'m*':>6}  {'outcome':}")
print("  " + "-" * 35)

# Include known prior data points
known = [
    (0.070, None,  "FAIL (prior, p=2677, BKZ-40 also fails)"),
    (0.340, 9,     "SUCCESS (prior, p≈524347, from 2026-07-26)"),
    (0.440, 4,     "SUCCESS (secp256k1 proxy lam/n, from 2026-06-15)"),
    (0.530, 4,     "SUCCESS (prior, p=211)"),
    (0.660, 3,     "SUCCESS (prior, p=2557)"),
]

all_data = sorted(threshold_data + [(r, m) for r, m, _ in known if r not in dict(threshold_data)],
                  key=lambda x: x[0])

# Merge prior data for display only
for ratio, min_m in sorted(threshold_data, key=lambda x: x[0]):
    status = f"m*={min_m}" if min_m else "FAIL"
    print(f"  {ratio:8.4f}  {str(min_m) if min_m else '--':>6}  {status}")

print("\nPrior data points (not re-run this session):")
for r, m, note in known:
    print(f"  {r:8.4f}  {str(m) if m else '--':>6}  {note}")

print("\nDone.")
