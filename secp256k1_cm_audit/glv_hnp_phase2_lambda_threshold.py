"""
GLV-HNP Phase 2: lambda/n threshold bisection study.

From 2026-07-26: LLL succeeds at lam/n >= 0.34, fails at lam/n = 0.07.
This script bisects the threshold by finding j=0 curves at lam/n ∈ {0.07..0.42}
and testing LLL recovery at each ratio.

Curve family: y^2 = x^3 + b over F_p, j=0, p ≡ 1 (mod 3).
GLV: k_full = k1 + lam*k2 (mod n), k1 < K1_BOUND, k2 < sqrt(n).

Run: python3 glv_hnp_phase2_lambda_threshold.py
Expected runtime: ~4-8 minutes.
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
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
    """Find a generator of order n on y^2 = x^3 + b."""
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
# Eisenstein CM for j=0 curves
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
# Find a curve with lam/n in target_range = (lo, hi)
# Search up to 20-bit primes; bias toward ~12-bit for speed.
# ---------------------------------------------------------------------------

def find_curve_in_range(target_lo, target_hi, bit_min=10, bit_max=16, seed_offset=0):
    """
    Find a j=0 prime-order curve with lam/n in (target_lo, target_hi).
    Returns (p, b, n, lam, G) or None.
    """
    p_lo = max(2**bit_min, 2**bit_min - 1)
    p = sympy.nextprime(p_lo + seed_offset)
    searched = 0
    while p < 2**bit_max and searched < 4000:
        searched += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 100: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    if target_lo <= ratio <= target_hi:
                        # Find the b for this twist
                        for b_try in range(1, 200):
                            P = None
                            for _ in range(100):
                                x = random.randint(0, p - 1)
                                rhs = (pow(x, 3, p) + b_try) % p
                                y = tonelli_shanks(rhs, p)
                                if y is not None and y != 0:
                                    P = (x, y)
                                    break
                            if P is None: continue
                            if ec_mul(P, n_cand, p) is None:
                                G = find_generator(p, b_try, n_cand)
                                if G is not None:
                                    return (p, b_try, n_cand, lam, G)
        p = sympy.nextprime(p)
    return None

# ---------------------------------------------------------------------------
# Signature generation with GLV nonces
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

# ---------------------------------------------------------------------------
# GLV lattice + recovery
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

def run_one(curve_params, m, d_secret, k1_bound, seed=42):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def test_curve(label, curve_params, seeds, m_max=14):
    """
    Run LLL at m=4..m_max with 3 seeds per m.
    Return the smallest m where all 3 seeds succeed (or None).
    """
    p, b, n, lam, G = curve_params
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    ratio = lam / n
    print(f"  Testing {label}: p={p}, n={n} ({n.bit_length()}b), "
          f"lam={lam}, lam/n={ratio:.3f}, K1={k1_bound}")
    for m in range(4, m_max + 1):
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 9999).randint(1, n - 1)
            if run_one(curve_params, m, d, k1_bound, seed):
                wins += 1
        print(f"    m={m}: {wins}/{len(seeds)}", end="", flush=True)
        if wins == len(seeds):
            print(f"  ← 3/3 PASS")
            return m
        print()
    print(f"    FAIL (never 3/3 up to m={m_max})")
    return None

# ---------------------------------------------------------------------------
# Main: bisect over lam/n bins
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-aware HNP Phase 2 — λ/n threshold bisection study")
print("Known: success at lam/n >= 0.34, failure at lam/n = 0.07")
print("Goal:  find the critical threshold in (0.07, 0.34)")
print("=" * 70)

SEEDS = [42, 1234, 9999]

# Target λ/n bins — 9 bins covering [0.07, 0.45]
# We use slightly wider windows (0.04 wide) to improve findability
BINS = [
    (0.07, 0.10, "bin_007_010"),
    (0.10, 0.13, "bin_010_013"),
    (0.13, 0.16, "bin_013_016"),
    (0.16, 0.19, "bin_016_019"),
    (0.19, 0.22, "bin_019_022"),
    (0.22, 0.26, "bin_022_026"),
    (0.26, 0.30, "bin_026_030"),
    (0.30, 0.34, "bin_030_034"),
    (0.34, 0.40, "bin_034_040"),   # known-good reference
]

results_by_bin = {}

for lo, hi, tag in BINS:
    mid = (lo + hi) / 2
    print(f"\n{'='*70}")
    print(f"Bin [{lo:.2f}, {hi:.2f}] (target lam/n ≈ {mid:.2f})")
    # Try progressively larger search ranges
    curve = None
    for seed_off in [0, 500, 2000, 8000]:
        curve = find_curve_in_range(lo, hi, bit_min=10, bit_max=17,
                                    seed_offset=seed_off)
        if curve is not None:
            break
    if curve is None:
        print(f"  SEARCH FAILED — no curve found in bin [{lo:.2f},{hi:.2f}]")
        results_by_bin[tag] = None
        continue
    m_success = test_curve(tag, curve, SEEDS, m_max=14)
    results_by_bin[tag] = (curve, m_success)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY TABLE — λ/n threshold bisection")
print(f"{'Bin [lo,hi]':<18} {'lam/n':<8} {'n bits':<8} {'result':<20}")
print("-" * 70)
for lo, hi, tag in BINS:
    entry = results_by_bin.get(tag)
    mid = (lo + hi) / 2
    if entry is None:
        print(f"[{lo:.2f},{hi:.2f}]         {mid:.3f}    —      CURVE NOT FOUND")
    else:
        curve, m_suc = entry
        p, b, n, lam, G = curve
        ratio = lam / n
        if m_suc is not None:
            status = f"LLL SUCCESS at m={m_suc}"
        else:
            status = "LLL FAIL (m≤14)"
        print(f"[{lo:.2f},{hi:.2f}]         {ratio:.3f}    "
              f"{n.bit_length():<6}   {status}")

# Identify threshold
print("\nThreshold analysis:")
prev_ok = True
threshold_lo, threshold_hi = None, None
for lo, hi, tag in BINS:
    entry = results_by_bin.get(tag)
    if entry is None:
        continue
    _, m_suc = entry
    ok = m_suc is not None
    if prev_ok and not ok:
        threshold_lo = lo
    if not prev_ok and ok:
        threshold_hi = hi
    prev_ok = ok

if threshold_lo is not None:
    print(f"  Threshold crossing: lam/n drops below {threshold_lo:.2f} → LLL fails")
elif threshold_hi is not None:
    print(f"  Threshold crossing: lam/n rises above {threshold_hi:.2f} → LLL succeeds")
else:
    # Check if we got clean monotone
    all_results = [(lo, hi, results_by_bin.get(tag))
                   for lo, hi, tag in BINS if results_by_bin.get(tag) is not None]
    all_pass = all(r[2][1] is not None for r in all_results)
    all_fail = all(r[2][1] is None for r in all_results)
    if all_pass:
        print("  All bins passed — threshold is below the lowest bin tested")
    elif all_fail:
        print("  All bins failed — threshold is above the highest bin tested")
    else:
        print("  Non-monotone result — likely noise; see per-bin data above")

print("\nDone.")
