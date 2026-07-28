"""
GLV-HNP Phase 2: lambda/n threshold bisection study.

Goal: identify the lambda/n ratio at which the GLV-aware HNP lattice attack
(column-scaled (2m+2)×(2m+2) lattice, LLL + BKZ) transitions from failure
to success.

Prior data:
  lambda/n = 0.07 (p=2677, n=2647):  LLL and BKZ(40) BOTH FAIL at all m<=12
  lambda/n = 0.34 (p=524347, n=523969): LLL succeeds at m=9

Strategy:
1. Scan primes p in [2^11, 2^16] for j=0 CM curves (p ≡ 1 mod 3).
   Use Eisenstein decomposition to get all 6 sextic-twist group orders
   without point-counting. Keep one representative per 0.02-wide
   lambda/n bucket in [0.06, 0.38].
2. For each bucket representative, run:
     (a) LLL sweep  m=4..18, 3 seeds, criterion: first m with 3/3
     (b) BKZ(40) at m=failure+3 for curves near the threshold
3. Report the boundary: largest lambda/n that never 3/3 at m<=18.

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

# ---------------------------------------------------------------------------
# CM theory helpers
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
    """Six Frobenius traces for the 6 sextic twists of a j=0 curve over F_p."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Smaller root of x^2 + x + 1 = 0 mod n (requires n ≡ 1 mod 3)."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def find_generator(p, b, n):
    """Find a generator of E(F_p): y^2=x^3+b of prime order n."""
    rng = random.Random(p ^ b)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def find_b_for_order(p, n):
    """Find b in 1..p-1 such that #E(F_p) with y^2=x^3+b equals n (prime)."""
    for b in range(1, min(p, 1000)):
        rng = random.Random(b)
        for _ in range(30):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b
                break
    # fallback: brute scan
    for b in range(1, p):
        rng = random.Random(b + 999)
        for _ in range(5):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y_sq = rhs
            y = tonelli_shanks(y_sq, p)
            if y is not None and y != 0:
                if ec_mul((x, y), n, p) is None:
                    return b
                break
    return None

# ---------------------------------------------------------------------------
# Scan for bucket representatives
# ---------------------------------------------------------------------------

BUCKET_WIDTH = 0.02
BUCKET_MIN   = 0.06
BUCKET_MAX   = 0.38

def bucket_index(ratio):
    return int((ratio - BUCKET_MIN) / BUCKET_WIDTH)

def bucket_center(idx):
    return BUCKET_MIN + (idx + 0.5) * BUCKET_WIDTH

N_BUCKETS = int((BUCKET_MAX - BUCKET_MIN) / BUCKET_WIDTH)  # 16 buckets

def scan_for_curves(p_min=2**11, p_max=2**16, verbose=True):
    """
    Scan primes p in [p_min, p_max] for j=0 CM curves.
    Collect one (p, n, lam, bucket_idx) representative per lambda/n bucket.
    Returns dict: bucket_idx -> (p, n, lam, ratio).
    """
    buckets = {}
    p = sympy.nextprime(p_min - 1)
    total_checked = 0
    curves_found = 0

    while p <= p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 4 or n_cand > p + 2:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    bidx = bucket_index(ratio)
                    if 0 <= bidx < N_BUCKETS and bidx not in buckets:
                        buckets[bidx] = (p, n_cand, lam, ratio)
                        curves_found += 1
                total_checked += 1
        if len(buckets) == N_BUCKETS:
            break
        p = sympy.nextprime(p)

    if verbose:
        print(f"  Scanned {total_checked} primes p ≡ 1 (mod 3), found {curves_found} bucket reps")
    return buckets

# ---------------------------------------------------------------------------
# GLV attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 300000:
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
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
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_KANNAN

def recover_d(M_reduced, m, n, S_KANNAN):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) == S_KANNAN:
            sign = 1 if last > 0 else -1
            d_cand = (sign * row[m]) % n
            if d_cand != 0:
                return d_cand
    return None

def run_once(p, b, n, lam, G, d_secret, m, k1_bound, seed, use_bkz=False, bkz_beta=40):
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d_cand = recover_d(reduced, m, n, S_KANNAN)
    return d_cand == d_secret

SEEDS = [42, 1234, 9999]
M_RANGE = list(range(4, 19))  # m=4..18

def sweep_lll(p, b, n, lam, G, k1_bound, label, verbose=True, bkz_fallback=False):
    """
    Run LLL sweep over m_range. Return (first_m_with_3of3_or_None, full_table).
    If bkz_fallback=True and LLL never succeeds, also try BKZ(40) at m=12,15,18.
    """
    d_secret = random.Random(p * 31 + n).randint(1, n - 1)
    if verbose:
        print(f"  {label}  [lam/n={lam/n:.4f}]")

    table = {}
    first_m = None
    for m in M_RANGE:
        wins = sum(run_once(p, b, n, lam, G, d_secret, m, k1_bound, seed)
                   for seed in SEEDS)
        table[m] = (wins, len(SEEDS))
        if wins == len(SEEDS) and first_m is None:
            first_m = m
        marker = " ← 3/3!" if wins == len(SEEDS) else ""
        if verbose:
            print(f"    LLL m={m}: {wins}/{len(SEEDS)}{marker}")
        if first_m is not None:
            break  # early exit on success

    # If LLL failed throughout, try BKZ(40) at m=12,15,18
    bkz_result = None
    if first_m is None and bkz_fallback:
        if verbose:
            print(f"    LLL failed at all m. Trying BKZ(40)...")
        for m in [12, 15, 18]:
            wins = sum(run_once(p, b, n, lam, G, d_secret, m, k1_bound, seed,
                                use_bkz=True, bkz_beta=40)
                       for seed in SEEDS)
            if verbose:
                print(f"    BKZ(40) m={m}: {wins}/{len(SEEDS)}")
            if wins == len(SEEDS):
                bkz_result = m
                break

    return first_m, table, bkz_result

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — lambda/n threshold bisection study")
print("=" * 70)
print(f"Scanning for j=0 CM curves in {N_BUCKETS} buckets "
      f"covering lambda/n ∈ [{BUCKET_MIN:.2f}, {BUCKET_MAX:.2f}]")
print(f"  Known: fails at lambda/n=0.07, succeeds at lambda/n=0.34\n")

# Step 1: Find bucket representatives
print("Step 1: Curve discovery scan (p ∈ [2^11, 2^16])...")
buckets = scan_for_curves(p_min=2**11, p_max=2**16, verbose=True)

# Also try a wider scan if some buckets are empty (low-ratio buckets are rare)
if len(buckets) < N_BUCKETS:
    missing = [i for i in range(N_BUCKETS) if i not in buckets]
    print(f"  Missing {len(missing)} buckets: {[f'{bucket_center(i):.2f}' for i in missing]}")
    print("  Extending scan to p ≤ 2^20...")
    extra = scan_for_curves(p_min=2**16 + 1, p_max=2**20, verbose=True)
    for k, v in extra.items():
        if k not in buckets:
            buckets[k] = v
    missing2 = [i for i in range(N_BUCKETS) if i not in buckets]
    if missing2:
        print(f"  Still missing: {[f'{bucket_center(i):.2f}' for i in missing2]}")

print(f"\nFound {len(buckets)} bucket representatives.")

# Step 2: Build curves (find b and G for each bucket rep)
print("\nStep 2: Building generators for each bucket representative...")
curves = {}  # bucket_idx -> (p, b, n, lam, G, k1_bound)
for bidx in sorted(buckets):
    p, n, lam, ratio = buckets[bidx]
    b = find_b_for_order(p, n)
    if b is None:
        print(f"  bucket {bidx} ({ratio:.4f}): SKIP — could not find b for n")
        continue
    G = find_generator(p, b, n)
    if G is None:
        print(f"  bucket {bidx} ({ratio:.4f}): SKIP — could not find generator")
        continue
    # k1_bound targeting eff ≈ 0.05
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    eff = k1_bound * k2_bound / n
    curves[bidx] = (p, b, n, lam, G, k1_bound)
    print(f"  bucket {bidx}: lam/n={ratio:.4f}  p={p}, n={n} ({n.bit_length()}b), "
          f"b={b}, k1={k1_bound}, eff={eff:.4f}")

# Step 3: LLL sweep for each curve
print("\nStep 3: LLL sweep (m=4..18, 3 seeds, 3/3 criterion)...")
print("  Near-threshold buckets (lam/n ∈ [0.10, 0.32]) get BKZ(40) fallback.\n")

results = {}  # bucket_idx -> {'ratio': float, 'first_m': int or None, 'bkz_first_m': int or None}
for bidx in sorted(curves):
    p, b, n, lam, G, k1_bound = curves[bidx]
    ratio = lam / n
    near_threshold = 0.10 <= ratio <= 0.32
    label = f"bucket {bidx:2d} lam/n={ratio:.4f}  p={p}, n={n}"
    first_m, table, bkz_m = sweep_lll(
        p, b, n, lam, G, k1_bound, label,
        verbose=True, bkz_fallback=near_threshold
    )
    results[bidx] = {'ratio': ratio, 'first_m': first_m, 'bkz_first_m': bkz_m}
    status = f"LLL 3/3 at m={first_m}" if first_m else (
             f"LLL FAIL, BKZ(40) 3/3 at m={bkz_m}" if bkz_m else "FAIL (LLL + BKZ)")
    print(f"  → {status}\n")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

print("=" * 70)
print("SUMMARY: lambda/n threshold for GLV-HNP Phase 2 LLL attack")
print("=" * 70)
print(f"  {'bucket':>6}  {'lam/n':>7}  {'first_m (LLL)':>14}  {'BKZ(40)':>9}  {'result':>10}")
print(f"  {'-'*6}  {'-'*7}  {'-'*14}  {'-'*9}  {'-'*10}")

last_fail_ratio = None
first_success_ratio = None

for bidx in sorted(results):
    r = results[bidx]
    ratio = r['ratio']
    lll_m = r['first_m']
    bkz_m = r['bkz_first_m']
    if lll_m is not None:
        result_str = f"LLL@m={lll_m}"
        if first_success_ratio is None:
            first_success_ratio = ratio
    elif bkz_m is not None:
        result_str = f"BKZ@m={bkz_m}"
        if first_success_ratio is None:
            first_success_ratio = ratio
    else:
        result_str = "FAIL"
        last_fail_ratio = ratio
    bkz_str = f"m={bkz_m}" if bkz_m else ("—" if lll_m is not None else "FAIL")
    print(f"  {bidx:6d}  {ratio:7.4f}  {str(lll_m) if lll_m else 'NEVER':>14}  "
          f"{bkz_str:>9}  {result_str}")

# Also show reference points from prior runs
print()
print("  Reference points (from prior runs):")
print(f"     lam/n=0.07  (p=2677):    LLL FAIL, BKZ(40) FAIL   [2026-07-26]")
print(f"     lam/n=0.34  (p=524347):  LLL 3/3 at m=9            [2026-07-26]")
print(f"     lam/n=0.53  (p=211):     LLL 3/3 at m=4            [2026-07-26]")

print()
if last_fail_ratio is not None and first_success_ratio is not None:
    print(f"  THRESHOLD INTERVAL: ({last_fail_ratio:.4f}, {first_success_ratio:.4f})")
    mid = (last_fail_ratio + first_success_ratio) / 2
    print(f"  MIDPOINT ESTIMATE:  lam/n ≈ {mid:.4f}")
elif first_success_ratio is not None and last_fail_ratio is None:
    print(f"  LLL succeeds at ALL tested ratios ≥ {first_success_ratio:.4f}")
elif last_fail_ratio is not None and first_success_ratio is None:
    print(f"  LLL fails at ALL tested ratios ≤ {last_fail_ratio:.4f}")
else:
    print("  No clear threshold found (insufficient data)")

print("\nDone.")
