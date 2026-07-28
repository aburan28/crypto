"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Prior runs found:
  - λ/n = 0.07 (n=2647): LLL and BKZ(40) both FAIL at m≤12
  - λ/n = 0.34 (n=523969, 20-bit): LLL succeeds at m=9
  - λ/n = 0.53/0.66 (8-bit, 12-bit): LLL succeeds at m=4–7

Goal: find the transition threshold in the 12-bit range (n ≈ 2^11..2^14)
at which LLL transitions from always-failing to eventually-succeeding.

Method:
  1. Enumerate primes p ≡ 1 (mod 3) in [2^11, 2^14], Eisenstein-decompose,
     collect (p, n, λ) for all 6 twists where n is prime.
  2. Sort by λ/n; pick ~10 representative curves covering 0.05..0.50.
  3. For each, run LLL sweep m=3..14 with 3 seeds.
  4. Report minimum m for 3/3 success (or "never" if m>14).

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (j=0, a=0 Weierstrass)
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
# CM / Eisenstein helpers
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
    return [2*a - b, -2*a + b, -(a+b), a+b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return min(r1, r2)

def is_prime_simple(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

# ---------------------------------------------------------------------------
# Collect curves by λ/n bucket
# ---------------------------------------------------------------------------

def collect_curves(p_min, p_max, bucket_width=0.05):
    """
    Enumerate primes p ≡ 1 (mod 3) in [p_min, p_max].
    For each, compute all 6 twists. Keep ones where n=p+1-t is prime and
    n ≡ 1 (mod 3) (so GLV eigenvalue exists).
    Returns list of (lam_ratio, p, n, lam, trace) sorted by lam_ratio.
    Also returns dict: bucket -> best (p,n,lam,lam_ratio).
    """
    candidates = []
    p = p_min | 1  # start odd
    while p <= p_max:
        if p % 3 == 1 and is_prime_simple(p):
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n = p + 1 - t
                    if n < 2 or not is_prime_simple(n): continue
                    if n % 3 != 1: continue
                    lam = glv_eigenvalue(n)
                    if lam is None: continue
                    ratio = lam / n
                    # canonical: ratio in (0, 0.5]
                    if ratio > 0.5: ratio = 1.0 - ratio - 1.0/n
                    candidates.append((ratio, p, n, lam, t))
        p += 2

    candidates.sort()

    # Bucket by bucket_width, pick one representative per bucket
    buckets = {}
    for ratio, p, n, lam, t in candidates:
        key = int(ratio / bucket_width)
        if key not in buckets:
            buckets[key] = (ratio, p, n, lam, t)

    return candidates, buckets

# ---------------------------------------------------------------------------
# Attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, seed=42):
    k2_bound = math.isqrt(n) + 1
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_lattice(sigs, n, lam, k1_bound):
    m = len(sigs)
    k2_bound = math.isqrt(n) + 1
    dim = 2 * m + 2
    S_K1 = max(1, n // k1_bound)
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
        M[m+1+i][i] = -lam * S_K1
        M[m+1+i][m+1+i] = S_K2
    for i in range(m):
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_lll, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_lll:
        last = row[dim-1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret: return True
    return False

def run_once(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam = curve
    G = find_generator(p, b, n)
    if G is None: return False
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep(curve, k1_bound, m_range, seeds):
    """
    Returns dict: m -> (wins, total).
    Stops early if we find the first m with all seeds winning.
    """
    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 31337).randint(1, curve[2] - 1)
            ok = run_once(curve, m, d, k1_bound, seed)
            wins += ok
        results[m] = (wins, len(seeds))
        if wins == len(seeds):
            break  # found threshold
    return results

# ---------------------------------------------------------------------------
# Find b (twist) for a given (p, n)
# ---------------------------------------------------------------------------

def find_b_for_n(p, n):
    """Find the twist b ∈ {1..p-1} such that #E(y²=x³+b, F_p) = n."""
    rng = random.Random(99999)
    # Try sextic twist representatives via small b values
    for b in list(range(1, 50)) + [rng.randint(50, p-1) for _ in range(200)]:
        G = find_generator(p, b, n)
        if G is not None:
            return b
    return None

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection")
print("=" * 70)
print()
print("Collecting 12–14-bit j=0 GLV curves...")

P_MIN = 2**11
P_MAX = 2**14
all_curves, buckets = collect_curves(P_MIN, P_MAX, bucket_width=0.05)

print(f"Found {len(all_curves)} (p,n,λ) triples in [{P_MIN}, {P_MAX}]")
print(f"  λ/n range: [{all_curves[0][0]:.3f}, {all_curves[-1][0]:.3f}]")
print(f"  Buckets (0.05-wide): {sorted(buckets.keys())}")
print()

# Select representative curves from buckets spanning 0.05..0.50
TARGET_RATIOS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
selected = []
for target in TARGET_RATIOS:
    key = int(target / 0.05)
    # Try the target bucket and adjacent ones
    for k in [key, key-1, key+1, key-2, key+2]:
        if k in buckets:
            ratio, p, n, lam, trace = buckets[k]
            b = find_b_for_n(p, n)
            if b is not None:
                selected.append((ratio, p, n, lam, b))
                break
    else:
        print(f"  WARNING: no curve found near λ/n≈{target:.2f}")

# Deduplicate by (p, n)
seen_pn = set()
unique = []
for item in selected:
    ratio, p, n, lam, b = item
    if (p, n) not in seen_pn:
        seen_pn.add((p, n))
        unique.append(item)
selected = sorted(unique, key=lambda x: x[0])

print(f"Selected {len(selected)} representative curves:")
for ratio, p, n, lam, b in selected:
    print(f"  λ/n={ratio:.3f}  p={p:5d}  n={n:5d} ({n.bit_length()}b)  "
          f"λ={lam:5d}  b={b}")
print()

# Attack parameters (fixed across all curves for fair comparison)
K1_BOUND = 8
SEEDS = [42, 1234, 9999]
M_RANGE = range(3, 16)

print(f"Attack parameters: K1_BOUND={K1_BOUND}, seeds={SEEDS}")
print(f"  (K2_BOUND = isqrt(n)+1 ≈ sqrt(n))")
print()
print("-" * 70)

results_summary = []

for ratio, p, n, lam, b in selected:
    k2_bound = math.isqrt(n) + 1
    eff = K1_BOUND * k2_bound / n
    curve = (p, b, n, lam)
    print(f"λ/n={ratio:.3f}  n={n} ({n.bit_length()}b)  eff={eff:.4f}")
    sys.stdout.flush()

    res = sweep(curve, K1_BOUND, M_RANGE, SEEDS)
    m_thresh = None
    row_parts = []
    for m, (wins, total) in sorted(res.items()):
        row_parts.append(f"m={m}:{wins}/{total}")
        if wins == total and m_thresh is None:
            m_thresh = m

    print(f"  {' '.join(row_parts)}")
    status = f"3/3 at m={m_thresh}" if m_thresh is not None else f"never 3/3 (m≤{max(res.keys())})"
    print(f"  => {status}")
    results_summary.append((ratio, n, m_thresh, status))
    print()
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print("=" * 70)
print("THRESHOLD SUMMARY")
print("=" * 70)
print(f"{'λ/n':>8} {'n (bits)':>10} {'min_m (3/3)':>15}")
print("-" * 40)
threshold_found = False
prev_success = None
for ratio, n, m_thresh, status in results_summary:
    bits = n.bit_length()
    m_str = str(m_thresh) if m_thresh is not None else ">15"
    print(f"{ratio:8.3f} {f'{n} ({bits}b)':>10} {m_str:>15}   {status}")
    if m_thresh is not None and not threshold_found and prev_success is None:
        prev_success = ratio
    if m_thresh is None and prev_success is None:
        pass  # still in failure region
    if m_thresh is not None and not threshold_found:
        threshold_found = True

# Find the transition (first success after a failure)
failures = [(r, n, m) for r, n, m, s in results_summary if m is None]
successes = [(r, n, m) for r, n, m, s in results_summary if m is not None]

print()
if failures and successes:
    last_fail = max(failures, key=lambda x: x[0])
    first_success = min(successes, key=lambda x: x[0])
    print(f"TRANSITION: failure at λ/n={last_fail[0]:.3f}, success at λ/n={first_success[0]:.3f}")
    print(f"  => threshold λ/n ≈ {(last_fail[0]+first_success[0])/2:.3f}  "
          f"[bounds: ({last_fail[0]:.3f}, {first_success[0]:.3f})]")
elif not failures:
    print(f"No failures found — threshold may be below λ/n={min(r for r,_,_,_ in results_summary):.3f}")
else:
    print(f"No successes found — threshold may be above λ/n={max(r for r,_,_,_ in results_summary):.3f}")

print()
print("Done.")
