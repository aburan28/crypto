"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Goal: find the critical λ/n ratio below which the column-scaled GLV lattice
attack (Phase 2, LLL) fails to recover d.

Prior data points:
  λ/n = 0.53 (n=199):    LLL 3/3 at m=6   ✓
  λ/n = 0.34 (n=523969): LLL 3/3 at m=9   ✓
  λ/n = 0.07 (n=2647):   LLL FAIL at all m≤12, BKZ(40) also fails ✗

This script:
1. Finds j=0 GLV curves with λ/n in bands {0.08..0.10, 0.10..0.13,
   0.13..0.17, 0.17..0.22, 0.22..0.28, 0.28..0.35, 0.35..0.50}.
2. For each band, runs LLL sweep over m=3..14, records min m for 3/3 success.
3. Produces a threshold table.

Run: python3 glv_hnp_phase2_lambda_threshold.py
Expected runtime: ~5-15 minutes (dominated by curve search and LLL).
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ── EC arithmetic ──────────────────────────────────────────────────────────

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
    m_v, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_v - i - 1), p)
        m_v, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(42)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ── CM / Eisenstein ────────────────────────────────────────────────────────

def eisenstein_decompose(p):
    """Find (a,b) with a²−ab+b²=p, a,b>=0 via O(sqrt(p)) search."""
    a = int(math.isqrt(p))
    while a >= 0:
        disc = 4 * p - 3 * a * a
        if disc < 0:
            a -= 1
            continue
        sq = int(math.isqrt(disc))
        if sq * sq == disc:
            for sign in (1, -1):
                b_num = a + sign * sq
                if b_num >= 0 and b_num % 2 == 0:
                    b = b_num // 2
                    if a * a - a * b + b * b == p and b >= 0:
                        return (a, b)
        a -= 1
    return None

def j0_traces(a, b):
    """The six possible traces for j=0 curves over F_p with Eisenstein norm a²-ab+b²=p."""
    return [2*a - b, -(2*a - b), 2*b - a, -(2*b - a), a + b, -(a + b)]

def glv_eigenvalue(n):
    """Smaller root of x²+x+1=0 mod n (requires n≡1 mod 3)."""
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

def find_curve_with_ratio(target_lo, target_hi, search_start=100, search_limit=50000,
                           require_lam_ratio_lo=None, require_lam_ratio_hi=None):
    """
    Search for a j=0 GLV prime-order curve with λ/n in [target_lo, target_hi].
    Returns (p, b_coeff, n, lam) or None.
    """
    lo = require_lam_ratio_lo if require_lam_ratio_lo is not None else target_lo
    hi = require_lam_ratio_hi if require_lam_ratio_hi is not None else target_hi
    p = sympy.nextprime(search_start)
    checked = 0
    while checked < search_limit:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 7:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam = glv_eigenvalue(n_cand)
                        if lam is None:
                            continue
                        ratio = lam / n_cand
                        if lo <= ratio <= hi:
                            # Find a b such that #E(F_p) = n_cand
                            for b_try in range(1, min(p, 300)):
                                P = find_generator(p, b_try, n_cand)
                                if P is not None:
                                    return (p, b_try, n_cand, lam, P)
        p = sympy.nextprime(p)
        checked += 1
    return None

# ── Lattice / LLL attack ───────────────────────────────────────────────────

def gen_glv_sigs(G, p, n, lam, d_secret, k1_bound, k2_bound, m, seed=42):
    """Generate m GLV-biased signatures with k1 ∈ [0,k1_bound), k2 ∈ [0,k2_bound)."""
    rng = random.Random(seed)
    sigs = []
    for _ in range(200 * m):
        if len(sigs) >= m:
            break
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
    return sigs if len(sigs) == m else None

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Build the (2m+2)×(2m+2) column-scaled GLV lattice."""
    m = len(sigs)
    dim = 2 * m + 2
    s_k1 = n // k1_bound
    s_d = 1
    s_k2 = n // max(k2_bound, 1)
    s_kan = n
    M = [[0] * dim for _ in range(dim)]
    # Mod-n rows
    for i in range(m):
        M[i][i] = n * s_k1
    # d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * s_k1
    M[m][m] = s_d
    # k2 rows
    for i in range(m):
        M[m + 1 + i][i] = -lam * s_k1
        M[m + 1 + i][m + 1 + i] = s_k2
    # Kannan row
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * s_k1
    M[2 * m + 1][dim - 1] = s_kan
    return M, s_d, s_kan

def try_recover_d(sigs, n, lam, k1_bound, k2_bound, d_secret, reduction='LLL', bkz_beta=20):
    """Build lattice, reduce, scan for d. Returns True if recovered."""
    M, s_d, s_kan = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    m = len(sigs)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if reduction == 'LLL':
        LLL.reduction(A)
    elif reduction == 'BKZ':
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    # Scan for Kannan row and extract d
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != s_kan:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[i][m]) % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return True
        # Flip sign (sometimes lattice finds negative Kannan row)
        d_cand2 = (-sign * A[i][m]) % n
        if d_cand2 == d_secret:
            return True
    return False

def sweep_m(G, p, n, lam, k1_bound, k2_bound, m_range, n_seeds=3,
            reduction='LLL', bkz_beta=20, label=""):
    """
    For each m in m_range, run n_seeds trials. Returns dict m -> wins.
    Stops early at first m achieving n_seeds/n_seeds.
    """
    rng = random.Random(777)
    results = {}
    for m in m_range:
        wins = 0
        for seed_idx in range(n_seeds):
            d_secret = rng.randint(1, n - 1)
            sigs = gen_glv_sigs(G, p, n, lam, d_secret, k1_bound, k2_bound, m,
                                 seed=seed_idx * 1000 + m)
            if sigs is None:
                continue
            if try_recover_d(sigs, n, lam, k1_bound, k2_bound, d_secret,
                              reduction=reduction, bkz_beta=bkz_beta):
                wins += 1
        results[m] = wins
        status = f"{wins}/{n_seeds}"
        print(f"  {label}m={m}: {status}")
        if wins == n_seeds:
            # Early exit once we hit perfect recovery
            for m2 in m_range:
                if m2 > m and m2 not in results:
                    results[m2] = -1  # not tested
            break
    return results

# ── Target ratio bands ─────────────────────────────────────────────────────
# We want curves with λ/n in each band. Note: for a j=0 curve,
# min(λ, n-1-λ)/n is the canonical ratio (always ≤ 0.5).
# Bands to test (lo, hi):
RATIO_BANDS = [
    (0.07, 0.09),   # near known failure
    (0.09, 0.12),
    (0.12, 0.16),
    (0.16, 0.21),
    (0.21, 0.27),
    (0.27, 0.34),   # just below known success at 0.34
    (0.34, 0.42),   # known success territory
    (0.42, 0.50),   # secp256k1 target (~0.44)
]

K1_BOUND = 2   # minimal bias: k1 ∈ {0, 1}
M_RANGE   = list(range(3, 15))   # m=3..14

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection")
print(f"K1_BOUND={K1_BOUND}, sweeping m={M_RANGE[0]}..{M_RANGE[-1]}, 3 seeds each")
print("=" * 70)
print()

# Results table: band_label -> (ratio, n, lam, min_m_for_3_3 or 'FAIL')
table = []

# ── Known prior data (no recomputation needed) ─────────────────────────────
KNOWN = [
    (0.07, 2647, 185, "FAIL (BKZ-40 also fails)"),
    (0.34, 523969, 177902, "9"),
    (0.53, 199, 106, "6"),
]

print("── Prior data points ──────────────────────────────────────────────────")
print(f"{'λ/n':>6}  {'n':>8}  {'λ':>8}  {'min_m':>10}")
for ratio, n, lam, min_m in KNOWN:
    print(f"{ratio:6.3f}  {n:8d}  {lam:8d}  {min_m:>10}")
print()

# ── Sweep each band ─────────────────────────────────────────────────────────
band_results = {}
for (lo, hi) in RATIO_BANDS:
    mid = (lo + hi) / 2
    print(f"── Band [{lo:.2f}, {hi:.2f}] (target≈{mid:.2f}) ──────────────────────────────")
    # Search with a smaller start for low-ratio bands (need more primes in range)
    search_start = 500
    result = find_curve_with_ratio(lo, hi, search_start=search_start, search_limit=80000)
    if result is None:
        print(f"  NO CURVE FOUND in [{lo:.2f},{hi:.2f}] — SKIP")
        band_results[(lo, hi)] = None
        print()
        continue
    p, b_coeff, n, lam, G = result
    ratio = lam / n
    k2_bound = math.isqrt(n) + 2
    eff = K1_BOUND * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else 999
    print(f"  Curve: p={p}, b={b_coeff}, n={n}, λ={lam}, λ/n={ratio:.4f}")
    print(f"  K2_BOUND={k2_bound}, eff={eff:.4f}, info-th. m≈{m_thresh:.1f}")

    res = sweep_m(G, p, n, lam, K1_BOUND, k2_bound, M_RANGE, n_seeds=3,
                  label=f"[{ratio:.3f}] ")
    # Find minimum m achieving 3/3
    min_m_full = None
    for m in M_RANGE:
        if res.get(m, 0) == 3:
            min_m_full = m
            break
    if min_m_full is None:
        outcome = f"FAIL (LLL ≤ m={M_RANGE[-1]})"
        print(f"  → LLL FAIL at all tested m")
    else:
        outcome = str(min_m_full)
        print(f"  → LLL succeeds at m≥{min_m_full}")

    band_results[(lo, hi)] = {
        'p': p, 'n': n, 'lam': lam, 'ratio': ratio,
        'k2_bound': k2_bound, 'min_m': min_m_full, 'outcome': outcome,
        'raw': res
    }
    print()

# ── Summary table ────────────────────────────────────────────────────────────
print("=" * 70)
print("SUMMARY TABLE")
print(f"{'band':>14}  {'λ/n':>6}  {'n':>8}  {'λ':>8}  {'min_m (3/3)':>14}  {'LLL status':>14}")
print("-" * 70)

for (lo, hi), d in band_results.items():
    if d is None:
        print(f"  [{lo:.2f},{hi:.2f}]  {'—':>6}  {'—':>8}  {'—':>8}  {'NO CURVE':>14}  {'SKIP':>14}")
    else:
        status = "✓ PASS" if d['min_m'] is not None else "✗ FAIL"
        min_m_str = str(d['min_m']) if d['min_m'] is not None else ">"+ str(M_RANGE[-1])
        print(f"  [{lo:.2f},{hi:.2f}]  {d['ratio']:6.3f}  {d['n']:8d}  {d['lam']:8d}  {min_m_str:>14}  {status:>14}")

print("-" * 70)
print("Prior known:")
for ratio, n, lam, min_m in KNOWN:
    status = "✗ FAIL" if "FAIL" in min_m else "✓ PASS"
    print(f"              {ratio:6.3f}  {n:8d}  {lam:8d}  {min_m:>14}  {status:>14}")
print()

# ── Threshold diagnosis ───────────────────────────────────────────────────────
print("=" * 70)
print("THRESHOLD DIAGNOSIS")
print()

passing_ratios = [d['ratio'] for d in band_results.values()
                  if d is not None and d['min_m'] is not None]
failing_ratios = [d['ratio'] for d in band_results.values()
                  if d is not None and d['min_m'] is None]
# Add prior known
failing_ratios.append(0.07)
passing_ratios.extend([0.34, 0.53])

if passing_ratios and failing_ratios:
    max_fail = max(failing_ratios)
    min_pass = min(passing_ratios)
    print(f"Largest confirmed FAIL: λ/n ≈ {max_fail:.3f}")
    print(f"Smallest confirmed PASS: λ/n ≈ {min_pass:.3f}")
    if max_fail < min_pass:
        print(f"Threshold bracket: ({max_fail:.3f}, {min_pass:.3f})")
        print(f"Threshold midpoint estimate: λ/n ≈ {(max_fail + min_pass)/2:.3f}")
    else:
        print("Non-monotone result — possible: threshold depends on n as well as λ/n.")
elif failing_ratios:
    print(f"All tested bands FAIL; threshold may be above max tested ({max(failing_ratios):.3f})")
elif passing_ratios:
    print(f"All tested bands PASS; threshold may be below min tested ({min(passing_ratios):.3f})")
print()
print("Structural note: when λ/n is small, the -λ·S_k1 entries in the k2-rows")
print("  have magnitude ~λ·n/K1 ≈ (λ/n)·n²/K1, comparable to the modular rows")
print("  (n·S_k1 = n²/K1). After LLL the k2-row cannot be 'pulled apart' from the")
print("  modular rows, destroying the bias signal. Above the threshold the k2-rows")
print("  are large enough (relative to modular rows) to be separated by LLL.")
