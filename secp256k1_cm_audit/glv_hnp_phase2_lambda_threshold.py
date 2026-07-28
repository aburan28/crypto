"""
Thread 20: lambda/n threshold sweep for GLV-HNP Phase 2 attack.

Known data points (from prior runs):
  lam/n = 0.07  (p=2677,  n=2647)   → LLL/BKZ(40) FAIL for all m≤12
  lam/n = 0.34  (p=524347, n=523969) → LLL 3/3 at m=9
  lam/n = 0.53  (p=211,   n=199)    → LLL 3/3 at m=4
  lam/n = 0.66  (p=2557,  n=2659)   → LLL 3/3 at m=7

Objective: bisect the critical lam/n below which LLL consistently fails.
Method: scan 12-14 bit j=0 prime-order curves; for each, compute lam/n,
        run LLL with fixed K1_BOUND=4, find minimum m for 3/3 success.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (a=0, y^2 = x^3 + b over F_p)
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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_point(p, b):
    rng = random.Random(b * p % 99991 + 1)
    for _ in range(p):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n):
    rng = random.Random(n * b % 99991 + 7)
    for _ in range(10000):
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
    """Find (a,b) with a^2 - a*b + b^2 = p."""
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
    """Six Frobenius traces for the sextic twists of j=0 over F_p."""
    return [2*a - b, -(2*a - b), a + b, -(a + b), 2*b - a, -(2*b - a)]

# ---------------------------------------------------------------------------
# GLV eigenvalue
# ---------------------------------------------------------------------------

def glv_eigenvalue(n):
    """
    Compute GLV eigenvalue: smaller root of x^2 + x + 1 = 0 mod n.
    Returns (lam, n-1-lam) with lam <= n//2, or (None, None) if no root.
    """
    if n % 3 != 1:
        return None, None
    neg3 = (n - 3) % n
    if pow(neg3, (n - 1) // 2, n) != 1:
        return None, None
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
# Find b such that #E(F_p) = target_n
# ---------------------------------------------------------------------------

def find_b_for_trace(p, target_n):
    """Brute-force: try b=1..min(200,p) until a point has order dividing target_n."""
    for b in range(1, min(200, p)):
        P = find_point(p, b)
        if P is None:
            continue
        if ec_mul(P, target_n, p) is None:
            return b
    return None

# ---------------------------------------------------------------------------
# Curve collection
# ---------------------------------------------------------------------------

def collect_curves(p_min=2000, p_max=16000, verbose=False):
    """
    Scan primes p in [p_min, p_max] with p ≡ 1 mod 6 (j=0 CM structure).
    For each prime-order j=0 curve, compute (p, b, n, lam, lam/n).
    Returns sorted list by lam/n.
    """
    curves = []
    seen_n = set()

    p = sympy.nextprime(p_min - 1)
    while p <= p_max:
        if p % 6 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n = p + 1 - t
                    if n < 4 or n in seen_n:
                        continue
                    if n % 3 != 1:
                        continue
                    if not sympy.isprime(n):
                        continue
                    lam, lam2 = glv_eigenvalue(n)
                    if lam is None:
                        continue
                    ratio = lam / n
                    b_coeff = find_b_for_trace(p, n)
                    if b_coeff is None:
                        continue
                    G = find_generator(p, b_coeff, n)
                    if G is None:
                        continue
                    seen_n.add(n)
                    curves.append((ratio, p, b_coeff, n, lam))
                    if verbose:
                        print(f"  p={p} n={n} lam={lam} lam/n={ratio:.4f}")
        p = sympy.nextprime(p)

    curves.sort()
    return curves

# ---------------------------------------------------------------------------
# GLV-HNP Phase 2 attack (from glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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
    return M, S_K1, S_D, S_K2, S_KANNAN

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return d_cand
    return None

def run_once(p, b_coeff, n, lam, G, k1_bound, m, seed):
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 13579).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b_coeff, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def find_min_m(p, b_coeff, n, lam, G, k1_bound, seeds, m_max):
    """Return minimum m such that all seeds succeed, or None if > m_max."""
    for m in range(3, m_max + 1):
        wins = sum(run_once(p, b_coeff, n, lam, G, k1_bound, m, s) for s in seeds)
        if wins == len(seeds):
            return m
    return None

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 68)
print("Thread 20: lambda/n threshold sweep (GLV-HNP Phase 2)")
print("=" * 68)

K1_BOUND = 8    # Match prior failing experiment (p=2677, K1=8 was the condition)
SEEDS = [42, 1234, 9999]
M_MAX = 18

# -- Phase A: collect curves ------------------------------------------------
print(f"\nPhase A: scanning primes p in [2000, 16000] for j=0 GLV curves...")
all_curves = collect_curves(p_min=2000, p_max=16000)
print(f"  Found {len(all_curves)} prime-order j=0 curves with n ≡ 1 mod 3 and GLV eigenvalue.")

# -- Phase B: select representatives by lam/n bin ---------------------------
BIN_WIDTH = 0.05
bins = {}
for ratio, p, b_coeff, n, lam in all_curves:
    bin_key = math.floor(ratio / BIN_WIDTH) * BIN_WIDTH
    if bin_key not in bins:
        bins[bin_key] = (ratio, p, b_coeff, n, lam)

print(f"\nPhase B: binning by lam/n (bin width {BIN_WIDTH})")
print(f"  Bins found: {sorted(bins.keys())}")

# Also add the known failing curve (p=2677, n=2647, lam=185) if not already present
KNOWN_CURVES = [
    (0.07, 2677, 2, 2647, 185, "KNOWN-FAIL"),
    (0.53, 211, 2, 199, 106, "KNOWN-PASS"),
    (0.66, 2557, 2, 2659, 1755, "KNOWN-PASS"),
]

# -- Phase C: run LLL attacks -----------------------------------------------
print(f"\nPhase C: LLL attack (K1_BOUND={K1_BOUND}, seeds={SEEDS}, m_max={M_MAX})")
print(f"{'lam/n':>8} {'p':>7} {'n':>7} {'lam':>7} {'min_m':>7} {'result':>10}")
print("-" * 55)

results = []

# Test known curves first
for ratio, p, b_coeff, n, lam, label in KNOWN_CURVES:
    G = find_generator(p, b_coeff, n)
    if G is None:
        print(f"  {ratio:.4f}: no generator found for known curve n={n}")
        continue
    min_m = find_min_m(p, b_coeff, n, lam, G, K1_BOUND, SEEDS, M_MAX)
    status = f"m={min_m}" if min_m else "FAIL"
    print(f"{ratio:>8.4f} {p:>7} {n:>7} {lam:>7} {min_m if min_m else '—':>7} {status:>10}  ({label})")
    results.append((ratio, p, n, lam, min_m, label))

print("  --- new curves ---")

# Test binned curves
for bin_key in sorted(bins.keys()):
    ratio, p, b_coeff, n, lam = bins[bin_key]
    G = find_generator(p, b_coeff, n)
    if G is None:
        print(f"  {ratio:.4f}: no generator for p={p}, n={n}")
        continue
    min_m = find_min_m(p, b_coeff, n, lam, G, K1_BOUND, SEEDS, M_MAX)
    status = f"m={min_m}" if min_m else "FAIL"
    print(f"{ratio:>8.4f} {p:>7} {n:>7} {lam:>7} {min_m if min_m else '—':>7} {status:>10}")
    results.append((ratio, p, n, lam, min_m, "scan"))

# -- Phase D: threshold summary --------------------------------------------
print("\n" + "=" * 68)
print("THRESHOLD SUMMARY")
print("=" * 68)

successes = [(r, m) for r, _, _, _, m, _ in results if m is not None]
failures  = [(r,)   for r, _, _, _, m, _ in results if m is None]

if successes and failures:
    min_success_ratio = min(r for r, _ in successes)
    max_fail_ratio    = max(r for (r,) in failures)
    print(f"\n  Lowest lam/n with LLL success:  {min_success_ratio:.4f}")
    print(f"  Highest lam/n with LLL failure: {max_fail_ratio:.4f}")
    if min_success_ratio > max_fail_ratio:
        print(f"\n  --> Threshold bracket: ({max_fail_ratio:.4f}, {min_success_ratio:.4f})")
        mid = (max_fail_ratio + min_success_ratio) / 2
        print(f"  --> Best estimate: lam/n_crit ≈ {mid:.4f}")
    else:
        print(f"\n  --> Threshold is not cleanly separated; success and failure overlap.")
        print(f"      Overlapping range: ({min_success_ratio:.4f}, {max_fail_ratio:.4f})")
elif successes:
    print(f"\n  All curves succeeded. Lowest lam/n tested: {min(r for r,_ in successes):.4f}")
    print(f"  Cannot establish lower bound from this data.")
elif failures:
    print(f"\n  All curves failed. Highest lam/n tested: {max(r for (r,) in failures):.4f}")
    print(f"  Cannot establish upper bound from this data.")

print("\n  lam/n sweep (all results):")
results.sort()
for ratio, p, n, lam, min_m, label in results:
    k2 = math.isqrt(n) + 1
    eff = K1_BOUND * k2 / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    status = f"3/3 at m={min_m}" if min_m else f"FAIL (m_thresh_theory≈{m_thresh:.0f})"
    print(f"    lam/n={ratio:.4f}  n={n:6d}  lam={lam:6d}  {status}")

print("\nDone.")
