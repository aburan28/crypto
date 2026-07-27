"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Goal: find the threshold in λ/n ∈ [0.07, 0.34] below which LLL fails
to recover d. Known boundary points:
  λ/n = 0.07  (p=2677): FAIL at m≤12, BKZ(40) also fails
  λ/n = 0.34  (p=523969, 20-bit): LLL 3/3 at m=9

Method:
  1. Enumerate 12-bit j=0 primes p ≡ 1 mod 3, collect curves with
     varying λ/n buckets [0.07..0.34].
  2. For each bucket, pick one representative curve.
  3. Run LLL sweep (m=6..14, 3 seeds) and record first m with 3/3.
  4. Also test BKZ(beta=20) at m=10 for each curve.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (a=0, short Weierstrass)
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
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory helpers (Eisenstein for j=0)
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
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
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
# Curve search: find j=0 prime-order curves bucketed by λ/n
# ---------------------------------------------------------------------------

BUCKETS = [
    (0.07, 0.09),   # known fail zone
    (0.09, 0.12),
    (0.12, 0.16),
    (0.16, 0.21),
    (0.21, 0.27),
    (0.27, 0.34),   # known success zone boundary
    (0.34, 0.42),   # known success
    (0.42, 0.50),   # near 0.44 (secp256k1's ratio)
]

def find_curves_by_lambda_ratio(bit_lo=11, bit_hi=13, target_count=2):
    """
    Enumerate j=0 curves in [2^bit_lo, 2^bit_hi], group by λ/n bucket.
    Returns dict: bucket_label -> list of (p, b, n, lam, G, ratio) (up to target_count each).
    """
    bucket_curves = {f"{lo:.2f}-{hi:.2f}": [] for (lo, hi) in BUCKETS}
    full_buckets = set()

    p = sympy.nextprime(2**bit_lo)
    p_max = 2**bit_hi
    examined = 0

    while p < p_max and len(full_buckets) < len(BUCKETS):
        examined += 1
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

                    for (lo, hi) in BUCKETS:
                        key = f"{lo:.2f}-{hi:.2f}"
                        if lo <= ratio < hi and len(bucket_curves[key]) < target_count:
                            # Find a b value for this curve
                            found_b = None
                            for b_try in range(1, p):
                                rng_tmp = random.Random(42)
                                for _ in range(100):
                                    x = rng_tmp.randint(0, p - 1)
                                    rhs = (pow(x, 3, p) + b_try) % p
                                    y = tonelli_shanks(rhs, p)
                                    if y is not None and y != 0:
                                        P_test = (x, y)
                                        if ec_mul(P_test, n_cand, p) is None:
                                            found_b = b_try
                                            break
                                if found_b is not None:
                                    break
                            if found_b is None:
                                continue
                            G = find_generator(p, found_b, n_cand)
                            if G is None:
                                continue
                            bucket_curves[key].append((p, found_b, n_cand, lam, G, ratio))
                            if len(bucket_curves[key]) >= target_count:
                                full_buckets.add(key)
                            break
        p = sympy.nextprime(p)

    print(f"  [scan] examined {examined} primes up to {p_max}")
    return bucket_curves

# ---------------------------------------------------------------------------
# Lattice experiment (same structure as glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_K1, S_D, S_K2, S_KANNAN

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
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
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return d_cand
    return None

def run_lll_experiment(curve_params, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G, ratio = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep_one_curve(label, curve_params, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G, ratio = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"
    print(f"  {label} [{algo}]  lam/n={ratio:.4f}  K1={k1_bound}  eff={eff:.4f}")
    results = {}
    first_33 = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 9999).randint(1, n - 1)
            ok = run_lll_experiment(curve_params, m, d_trial, k1_bound, seed,
                                    use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
        status = f"{wins}/{len(seeds)}"
        print(f"    m={m}: {status}")
        if wins == len(seeds) and first_33 is None:
            first_33 = m
    return results, first_33

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection study")
print("Prior data: fail at λ/n=0.07, success at λ/n=0.34")
print("=" * 70)

SEEDS = [42, 1234, 9999]

# Step 1: Enumerate curves across λ/n buckets (12-bit range)
print("\nStep 1: Scanning 12-bit j=0 curves for λ/n buckets...")
bucket_curves = find_curves_by_lambda_ratio(bit_lo=11, bit_hi=14, target_count=1)

print("\nCurves found per bucket:")
all_curves = []
for (lo, hi) in BUCKETS:
    key = f"{lo:.2f}-{hi:.2f}"
    curves = bucket_curves[key]
    if curves:
        p, b, n, lam, G, ratio = curves[0]
        print(f"  λ/n∈[{lo:.2f},{hi:.2f}): p={p}, n={n}, lam={lam}, ratio={ratio:.4f}")
        all_curves.append((key, lo, hi, curves[0]))
    else:
        print(f"  λ/n∈[{lo:.2f},{hi:.2f}): (no curve found in scan)")

# Step 2: LLL sweep for each found curve
print("\nStep 2: LLL sweep per λ/n bucket (m=6..16, 3 seeds, K1=8 fixed)")
print("NOTE: K1=8 fixed to match known-failure experiment (p=2677, lam/n=0.07).")
print("K1=2 trivially succeeds always; K1=8 tests whether lam/n is causal.")
print("-" * 70)

# K1=8 FIXED: matches the known-failure experiment. With K1=2, every curve
# trivially succeeds since k1 in {0,1} gives an essentially binary constraint.
K1_FIXED = 8
threshold_table = []  # (ratio, first_33_m or None)

for key, lo, hi, curve_params in all_curves:
    p, b, n, lam, G, ratio = curve_params
    k1_bound = K1_FIXED
    label = f"λ/n∈[{lo:.2f},{hi:.2f}) p={p}"
    res, first_33 = sweep_one_curve(label, curve_params, k1_bound,
                                    m_range=range(6, 17), seeds=SEEDS)
    threshold_table.append((ratio, first_33, key, p, n, lam))

# Step 3: BKZ(20) on the curves that failed LLL
print("\nStep 3: BKZ(20) rescue on LLL-failure curves")
print("-" * 70)

bkz_results = {}
for ratio, first_33, key, p, n, lam, in threshold_table:
    if first_33 is None:
        # Find the curve again
        for (k2, lo, hi, cp) in all_curves:
            if k2 == key:
                k1_bound = K1_FIXED
                label = f"λ/n∈[{lo:.2f},{hi:.2f}) BKZ(20) p={p}"
                res_bkz, first_33_bkz = sweep_one_curve(
                    label, cp, k1_bound,
                    m_range=range(8, 17), seeds=SEEDS,
                    use_bkz=True, bkz_beta=20
                )
                bkz_results[key] = first_33_bkz

# Step 4: Summary table
print("\n" + "=" * 70)
print("THRESHOLD SUMMARY TABLE")
print(f"{'λ/n range':<22} {'p':>8} {'n':>6} {'λ':>6} {'LLL 3/3 at m':<18} {'BKZ(20) 3/3'}")
print("-" * 70)

for ratio, first_33, key, p, n, lam in sorted(threshold_table, key=lambda x: x[0]):
    lll_str = f"m={first_33}" if first_33 else "NEVER"
    bkz_str = ""
    if first_33 is None and key in bkz_results:
        bkz_str = f"m={bkz_results[key]}" if bkz_results[key] else "NEVER"
    print(f"  {key:<20}  {p:>8}  {n:>6}  {lam:>6}  {lll_str:<18} {bkz_str}")

# Step 5: Interpret threshold
print("\nINTERPRETATION:")
threshold_candidates = [(ratio, first_33) for ratio, first_33, *_ in threshold_table
                        if first_33 is not None]
failure_ratios = [ratio for ratio, first_33, *_ in threshold_table if first_33 is None]

if threshold_candidates and failure_ratios:
    min_success_ratio = min(r for r, _ in threshold_candidates)
    max_failure_ratio = max(failure_ratios)
    print(f"  LLL failure zone:  λ/n ≤ {max_failure_ratio:.4f}")
    print(f"  LLL success zone:  λ/n ≥ {min_success_ratio:.4f}")
    print(f"  Threshold bracket: ({max_failure_ratio:.4f}, {min_success_ratio:.4f})")
    midpoint = (max_failure_ratio + min_success_ratio) / 2
    print(f"  Estimated threshold: ~{midpoint:.4f}")
    print(f"  secp256k1 GLV: λ/n ≈ 0.44 → {'ABOVE' if 0.44 > min_success_ratio else 'BELOW'} threshold")
elif not failure_ratios:
    print(f"  All curves succeeded — threshold is BELOW {min(r for r, *_ in threshold_table):.4f}")
elif not threshold_candidates:
    print(f"  All curves failed — threshold is ABOVE {max(r for r, *_ in threshold_table):.4f}")

print("\nDone.")
