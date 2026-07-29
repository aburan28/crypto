"""
GLV-HNP Phase 2 — λ/n threshold bisection study.

Goal: Find the critical λ/n ratio below which LLL fails to recover d,
bisecting the known bracket [0.07 (fails), 0.34 (succeeds)].

Strategy:
  1. Search for j=0 GLV curves (p ≡ 1 mod 3, n prime, n ≡ 1 mod 3) in
     12–16 bit range.
  2. Bucket by λ/n in bands of width 0.05.
  3. For each bucket, pick one representative curve and run LLL at 3 seeds,
     sweeping m from 3 to m_thresh+5.
  4. Record minimum m giving 3/3 success (or FAIL if never).

Known data points going in:
  λ/n ≈ 0.07 (p=2677, n=2647): LLL FAIL at all m≤12, BKZ(40) FAIL
  λ/n ≈ 0.44 (secp256k1 proxy, lam/n=0.44): SUCCESS (from Phase 2 analysis)
  λ/n ≈ 0.53 (p=211, n=199): SUCCESS at m=6
  λ/n ≈ 0.66 (p=2557, n=2659): SUCCESS at m=7

Run: python3 glv_hnp_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass a=0)
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
# CM / Eisenstein
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
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

# ---------------------------------------------------------------------------
# Signature generation and lattice
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
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

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

def run_lll_experiment(curve_params, m, d_secret, k1_bound, seed=42):
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

# ---------------------------------------------------------------------------
# Curve finder: collect curves across λ/n spectrum
# ---------------------------------------------------------------------------

SEEDS_TEST = [42, 1234, 9999]

def find_curves_by_lambda_ratio(bit_lo=10, bit_hi=16, max_per_bucket=1, bucket_width=0.05):
    """
    Search primes p in [2^bit_lo, 2^bit_hi] with p ≡ 1 mod 3.
    For each, compute Eisenstein traces → candidate group orders n.
    Filter n prime and n ≡ 1 mod 3.
    Bucket by λ/n in bands of width bucket_width.
    Return dict: bucket_id -> list of (p, b, n, lam, G, lam/n).
    """
    buckets = {}
    n_buckets = int(1.0 / bucket_width)
    for i in range(n_buckets):
        lo = i * bucket_width
        hi = (i + 1) * bucket_width
        buckets[i] = {'range': (lo, hi), 'curves': []}

    p = sympy.nextprime(2**bit_lo - 1)
    p_limit = 2**bit_hi
    checked = 0
    while p < p_limit:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis:
                a_e, b_e = eis
                traces = j0_traces(a_e, b_e)
                for t in traces:
                    n_cand = p + 1 - t
                    if n_cand < 8: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    bkt = min(int(ratio / bucket_width), n_buckets - 1)
                    if len(buckets[bkt]['curves']) >= max_per_bucket:
                        continue
                    # Find a twist b such that #E(F_p) = n_cand
                    found_b = None
                    for b_try in range(1, min(p, 500)):
                        Pt = None
                        rng2 = random.Random(b_try)
                        for _ in range(200):
                            x = rng2.randint(0, p - 1)
                            rhs = (pow(x, 3, p) + b_try) % p
                            y = tonelli_shanks(rhs, p)
                            if y and y != 0:
                                Pt = (x, y)
                                break
                        if Pt is None: continue
                        if ec_mul(Pt, n_cand, p) is None:
                            found_b = b_try
                            break
                    if found_b is None: continue
                    G_cand = find_generator(p, found_b, n_cand)
                    if G_cand is None: continue
                    buckets[bkt]['curves'].append((p, found_b, n_cand, lam, G_cand, ratio))
        p = sympy.nextprime(p)
        checked += 1
    return buckets

# ---------------------------------------------------------------------------
# Per-curve threshold sweep
# ---------------------------------------------------------------------------

def sweep_threshold(curve_entry, k1_bound, m_max=20):
    """
    Run LLL at m=3..m_max with 3 seeds.
    Return (min_m_for_3of3, results_dict) where results_dict maps m -> wins.
    """
    p, b, n, lam, G, ratio = curve_entry
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1.0 else 9999

    results = {}
    min_m_3of3 = None
    for m in range(3, m_max + 1):
        wins = 0
        for seed in SEEDS_TEST:
            d_trial = random.Random(seed + 9999).randint(1, n - 1)
            ok = run_lll_experiment((p, b, n, lam, G), m, d_trial, k1_bound, seed)
            wins += ok
        results[m] = wins
        if wins == len(SEEDS_TEST) and min_m_3of3 is None:
            min_m_3of3 = m
        # Early stop: if 3/3 at two consecutive m, we're done
        if min_m_3of3 is not None and m >= min_m_3of3 + 1:
            break
    return min_m_3of3, results, m_thresh

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection")
print("Searching for j=0 GLV curves across the λ/n spectrum [0, 0.50)")
print("Known: λ/n=0.07 → FAIL; λ/n=0.34+ → PASS")
print("=" * 70)

# K1_BOUND for the sweep: set so the information-theoretic threshold is ~5-8
# We use K1 = max(2, int(0.05 * sqrt(n))) as in glv_hnp_phase2_20bit.py
# but we'll compute per-curve.

BUCKET_WIDTH = 0.05
BIT_LO = 10
BIT_HI = 16
MAX_PER_BUCKET = 2

print(f"\nPhase 1: Collecting curves in [{BIT_LO},{BIT_HI})-bit range...")
buckets = find_curves_by_lambda_ratio(
    bit_lo=BIT_LO, bit_hi=BIT_HI,
    max_per_bucket=MAX_PER_BUCKET, bucket_width=BUCKET_WIDTH
)

# Print bucket census
print("\nBucket census:")
print(f"  {'Band':12s}  {'Curves':6s}  {'Representatives'}")
for bkt_id in sorted(buckets.keys()):
    lo, hi = buckets[bkt_id]['range']
    curves = buckets[bkt_id]['curves']
    if curves:
        reps = ", ".join(f"n={c[2]}(λ/n={c[5]:.3f})" for c in curves)
    else:
        reps = "(none found)"
    print(f"  [{lo:.2f},{hi:.2f})   {len(curves):6d}  {reps}")

# Also add known failure curve (p=2677, lam/n=0.07)
known_fail_p = 2677
known_fail_b = 2
known_fail_n = 2647
known_fail_lam = 185
G_fail = find_generator(known_fail_p, known_fail_b, known_fail_n)
if G_fail:
    ratio_fail = known_fail_lam / known_fail_n
    known_fail_entry = (known_fail_p, known_fail_b, known_fail_n,
                        known_fail_lam, G_fail, ratio_fail)
    print(f"\n  + Known failure injected: p={known_fail_p}, n={known_fail_n}, "
          f"λ/n={ratio_fail:.4f}")
else:
    known_fail_entry = None
    print("  (Warning: could not find generator for known failure curve)")

print("\nPhase 2: LLL threshold sweep per bucket...")
print("=" * 70)

# Collect results: ratio -> (min_m or None, detailed results)
all_results = []  # list of (ratio, min_m, n, p)

for bkt_id in sorted(buckets.keys()):
    lo, hi = buckets[bkt_id]['range']
    curves = buckets[bkt_id]['curves']
    if not curves:
        continue
    # Use first representative
    entry = curves[0]
    p, b, n, lam, G, ratio = entry
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    k2_bound = math.isqrt(n) + 1

    print(f"\nBucket [{lo:.2f},{hi:.2f}): p={p}, n={n}, λ={lam}, "
          f"λ/n={ratio:.4f}, K1={k1_bound}, K2={k2_bound}")
    min_m, res, m_thresh = sweep_threshold(entry, k1_bound, m_max=18)
    all_results.append((ratio, min_m, n, p))

    for m, wins in sorted(res.items()):
        thresh_marker = " ← thresh" if m == m_thresh else ""
        ok_marker = " ← FIRST 3/3" if m == min_m else ""
        print(f"  m={m:2d}: {wins}/{len(SEEDS_TEST)}{thresh_marker}{ok_marker}")

    if min_m is not None:
        print(f"  → PASS: min m for 3/3 = {min_m}  (m_thresh={m_thresh:.0f})")
    else:
        print(f"  → FAIL: never 3/3 in m=3..18  (m_thresh={m_thresh:.0f})")

# Add known failure curve
if known_fail_entry is not None:
    p, b, n, lam, G, ratio = known_fail_entry
    k1_bound = 8  # same as in 20-bit script
    k2_bound = math.isqrt(n) + 1
    print(f"\nKnown failure curve: p={p}, n={n}, λ={lam}, "
          f"λ/n={ratio:.4f}, K1={k1_bound}, K2={k2_bound}")
    min_m, res, m_thresh = sweep_threshold(known_fail_entry, k1_bound, m_max=14)
    all_results.append((ratio, min_m, n, p))
    for m, wins in sorted(res.items()):
        print(f"  m={m:2d}: {wins}/{len(SEEDS_TEST)}")
    if min_m is not None:
        print(f"  → PASS: min m for 3/3 = {min_m}")
    else:
        print(f"  → FAIL: never 3/3 in m=3..14")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("THRESHOLD SUMMARY TABLE")
print("=" * 70)
print(f"  {'λ/n band':12s}  {'n':8s}  {'p':8s}  {'min m (3/3)':>12s}  {'status'}")
print(f"  {'-'*12}  {'-'*8}  {'-'*8}  {'-'*12}  {'-'*8}")

all_results.sort(key=lambda x: x[0])
pass_ratios = []
fail_ratios = []
for ratio, min_m, n, p in all_results:
    status = f"PASS (m={min_m})" if min_m is not None else "FAIL"
    print(f"  {ratio:.4f}        {n:8d}  {p:8d}  {str(min_m) if min_m else 'never':>12s}  {status}")
    if min_m is not None:
        pass_ratios.append(ratio)
    else:
        fail_ratios.append(ratio)

if fail_ratios and pass_ratios:
    max_fail = max(fail_ratios)
    min_pass = min(pass_ratios)
    print(f"\nBracket: FAIL at λ/n ≤ {max_fail:.4f}, PASS at λ/n ≥ {min_pass:.4f}")
    print(f"Threshold estimated in [{max_fail:.4f}, {min_pass:.4f}]")
elif pass_ratios:
    print(f"\nAll tested curves PASS.  Min ratio: {min(pass_ratios):.4f}")
elif fail_ratios:
    print(f"\nAll tested curves FAIL.  Max ratio: {max(fail_ratios):.4f}")

print("\nDone.")
