"""
GLV-aware HNP Phase 2 — λ/n threshold bisection (Thread 20).

Prior data (from 2026-07-26 log):
  p=2677,   n=2647,   λ=185,    λ/n=0.070 → FAIL (LLL+BKZ-40)
  p=524347, n=523969, λ=177902, λ/n=0.340 → PASS (3/3 at m=9)

Goal: bisect threshold by finding y²=x³+2 curves over F_p with
λ_min/n in [0.07, 0.34] and testing LLL recovery rate.

λ_min = min(λ, n-1-λ) is the smaller root of x²+x+1≡0 mod n.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import random
import math
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Number theory helpers
# ---------------------------------------------------------------------------
def is_prime_miller_rabin(n, witnesses=(2, 3, 5, 7, 11, 13, 17, 19, 23)):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in witnesses:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def tonelli_shanks(n, p):
    """Square root of n mod p (p prime). Returns r s.t. r^2=n mod p, or None."""
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, rx = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return rx
        i, tmp = 1, t * t % p
        while tmp != 1:
            tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, rx = i, b * b % p, t * b * b % p, rx * b % p

def count_points_brute(p, b):
    """Count points on y^2 = x^3 + b over F_p using Legendre symbols."""
    half = (p - 1) // 2
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x * x * x + b) % p
        if rhs == 0:
            count += 1
        else:
            leg = pow(rhs, half, p)
            if leg == 1:
                count += 2
    return count

def find_lambda_roots(n):
    """
    Roots of x^2+x+1=0 mod n (assumes n prime).
    Returns (lam1, lam2) or None if no roots exist.
    """
    # discriminant = -3 mod n
    disc = (-3) % n
    sq = tonelli_shanks(disc, n)
    if sq is None:
        return None
    inv2 = pow(2, -1, n)
    lam1 = ((-1 + sq) * inv2) % n
    lam2 = ((-1 - sq) * inv2) % n
    return (lam1, lam2)

def find_point(p, b):
    """
    Find any non-identity point on y^2=x^3+b over F_p.
    Since n is prime (checked in search_curves), any non-identity point is a generator.
    """
    for x in range(p):
        rhs = (x * x * x + b) % p
        if rhs == 0:
            return (x, 0)
        sq = tonelli_shanks(rhs, p)
        if sq is not None:
            return (x, sq)
    return None

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass a=0)
# ---------------------------------------------------------------------------
def modinv(a, m): return pow(a, -1, m)

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
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

# ---------------------------------------------------------------------------
# Signature generation (GLV model)
# ---------------------------------------------------------------------------
def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 10000:
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

# ---------------------------------------------------------------------------
# GLV lattice and recovery
# ---------------------------------------------------------------------------
def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1
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
            return d_cand
    return None

def run_attack(p, n, lam, G, k1_bound, m, d_secret, seed=42, use_bkz=False):
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, 2, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False

    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(20))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Search: find y^2=x^3+2 curves with λ_min/n in target ranges
# ---------------------------------------------------------------------------
def search_curves(p_max=12000, b=2):
    """
    Enumerate primes p ≡ 1 (mod 3), p ≤ p_max.
    For each, compute n = #E(F_p) and find λ roots of x²+x+1=0 mod n.
    Return list of (p, n, lam_min, ratio) where ratio = lam_min/n.
    """
    results = []
    for p in range(5, p_max + 1, 2):
        if p % 3 != 1: continue
        if not is_prime_miller_rabin(p): continue
        n = count_points_brute(p, b)
        if not is_prime_miller_rabin(n): continue
        if n == p: continue  # trivial case
        roots = find_lambda_roots(n)
        if roots is None: continue
        lam1, lam2 = roots
        # Verify: both should satisfy x^2+x+1=0 mod n
        assert (lam1 * lam1 + lam1 + 1) % n == 0
        # λ_min = smaller of the two roots
        lam_min = min(lam1, lam2)
        ratio = lam_min / n
        results.append((p, n, lam_min, ratio))
    return results

# ---------------------------------------------------------------------------
# Main: search + attack for threshold bisection
# ---------------------------------------------------------------------------
print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection (Thread 20)")
print("Known: λ/n=0.070 FAILS, λ/n=0.340 PASSES")
print("=" * 70)

print("\n[Step 1] Searching for curves with λ_min/n in [0.07, 0.35] ...")
print("  (brute-force point counting, p ≤ 12000, b=2)")

curves = search_curves(p_max=12000, b=2)

# Add known data points
known = [
    (2677, 2647, 185, 185/2647),    # λ/n=0.070, FAIL (from log)
    (2557, 2659, 903, 903/2659),    # λ/n=0.340, PASS (from scaling test)
    (524347, 523969, 177902, 177902/523969),  # λ/n=0.340, PASS (from 20-bit)
]

# Filter curves in target range
target = [(p, n, lam, r) for (p, n, lam, r) in curves if 0.07 <= r <= 0.35]

print(f"\n  Found {len(target)} curves with λ_min/n in [0.07, 0.35]:")
print(f"  {'p':>8} {'n':>8} {'λ_min':>8} {'λ/n':>8}")
for p, n, lam, r in sorted(target, key=lambda x: x[3]):
    print(f"  {p:>8} {n:>8} {lam:>8} {r:>8.4f}")

# Bin curves into 5 ratio windows
bins = {
    "0.07-0.12": [],
    "0.12-0.17": [],
    "0.17-0.22": [],
    "0.22-0.27": [],
    "0.27-0.35": [],
}
for p, n, lam, r in target:
    if r < 0.12:    bins["0.07-0.12"].append((p, n, lam, r))
    elif r < 0.17:  bins["0.12-0.17"].append((p, n, lam, r))
    elif r < 0.22:  bins["0.17-0.22"].append((p, n, lam, r))
    elif r < 0.27:  bins["0.22-0.27"].append((p, n, lam, r))
    else:           bins["0.27-0.35"].append((p, n, lam, r))

# Pick one representative from each non-empty bin
representatives = []
for bin_name, entries in bins.items():
    if entries:
        # Pick entry closest to bin midpoint
        lo, hi = float(bin_name.split("-")[0]), float(bin_name.split("-")[1])
        mid = (lo + hi) / 2
        best = min(entries, key=lambda x: abs(x[3] - mid))
        representatives.append((bin_name, *best))

print(f"\n  Representatives (one per λ/n window):")
print(f"  {'bin':>12} {'p':>8} {'n':>8} {'λ_min':>8} {'λ/n':>8}")
for row in representatives:
    print(f"  {row[0]:>12} {row[1]:>8} {row[2]:>8} {row[3]:>8} {row[4]:>8.4f}")

# ---------------------------------------------------------------------------
print("\n[Step 2] Running GLV-HNP LLL attack on each representative ...")
print("  Parameters: K1_BOUND=4, m=6, 5 seeds")
print("  (K1_BOUND=4 means k1 is 2-bit biased, small relative to sqrt(n))")
print()

SEEDS = [42, 1234, 9999, 7, 31415]
M_VAL = 6
K1_BOUND = 4

print(f"  {'λ/n':>8} {'p':>7} {'n':>7} {'lam':>7} {'LLL':>6} {'BKZ20':>6}")

results_by_ratio = []

for row in representatives:
    bin_name, p, n, lam, ratio = row

    # Find a point on the curve (= generator since n is prime)
    G = find_point(p, 2)
    if G is None:
        print(f"  {ratio:>8.4f} {p:>7} {n:>7} {lam:>7}  [no point found]")
        continue
    if ec_mul(G, n, p) is not None:
        print(f"  {ratio:>8.4f} {p:>7} {n:>7} {lam:>7}  [bad order]")
        continue

    lll_wins = 0
    bkz_wins = 0
    for seed in SEEDS:
        d_secret = random.Random(seed + 1111).randint(1, n - 1)
        if run_attack(p, n, lam, G, K1_BOUND, M_VAL, d_secret, seed, use_bkz=False):
            lll_wins += 1
        if run_attack(p, n, lam, G, K1_BOUND, M_VAL, d_secret, seed, use_bkz=True):
            bkz_wins += 1

    lll_str = f"{lll_wins}/{len(SEEDS)}"
    bkz_str = f"{bkz_wins}/{len(SEEDS)}"
    print(f"  {ratio:>8.4f} {p:>7} {n:>7} {lam:>7}  {lll_str:>6}  {bkz_str:>6}")
    results_by_ratio.append((ratio, lll_wins, bkz_wins, len(SEEDS)))

# ---------------------------------------------------------------------------
print("\n[Step 3] Summary")
print("  λ/n    LLL     BKZ-20  verdict")
for ratio, lll_w, bkz_w, total in sorted(results_by_ratio, key=lambda x: x[0]):
    lll_ok = "PASS" if lll_w == total else ("PART" if lll_w > 0 else "FAIL")
    bkz_ok = "PASS" if bkz_w == total else ("PART" if bkz_w > 0 else "FAIL")
    print(f"  {ratio:.4f}  {lll_w}/{total}  {bkz_w}/{total}  LLL={lll_ok} BKZ={bkz_ok}")

print("\nKnown anchors:")
print("  0.070  (p=2677, n=2647): FAIL (LLL+BKZ-40, from 2026-07-26 log)")
print("  0.340  (p=2557, n=2659 / p=524347, n=523969): PASS")
