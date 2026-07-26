"""
GLV-HNP Phase 2: Lambda/n threshold bisection study.

Sweeps j=0 curves in the 12-16 bit range with varying lambda/n ratios
and tests LLL recovery success as a function of lambda/n.

Goal: find the critical lambda/n ratio below which LLL fails.

Prior data:
  lam/n = 0.07 (n=2647,  p=2677):   FAIL (LLL + BKZ(40))
  lam/n = 0.34 (n=523969, p=524347): SUCCESS
  lam/n = 0.53 (n=199,   p=211):    SUCCESS
  lam/n = 0.66 (n=2659,  p=2557):   SUCCESS

Expected finding: there is a threshold tau in (0.07, 0.34) below which
the lambda-row contribution is too small for LLL to isolate the k2 constraint.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys

try:
    import sympy
    from fpylll import IntegerMatrix, LLL
except ImportError as e:
    sys.exit(f"Missing dependency: {e}. Install with: pip install sympy fpylll")

# ---------------------------------------------------------------------------
# Tonelli-Shanks square root mod p
# ---------------------------------------------------------------------------

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m_ts, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_ts - i - 1), p)
        m_ts, c, t, r = i, b * b % p, t * b * b % p, r * b % p

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0, y^2 = x^3 + b)
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
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def find_generator(p, b, n, max_tries=5000):
    rng = random.Random(42)
    for _ in range(max_tries):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition and GLV eigenvalue
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - a*b + b^2 = p, a,b >= 0."""
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
    """Six Frobenius traces for j=0 curve over F_p given Eisenstein (a,b)."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Smaller root of x^2 + x + 1 = 0 mod n; returns None if n%3 != 1."""
    if n % 3 != 1: return None
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return min(r1, r2)

# ---------------------------------------------------------------------------
# Curve search: find j=0 curves with specific lambda/n ratios
# ---------------------------------------------------------------------------

def search_curves(bit_lo, bit_hi, target_ratios, per_ratio=3, max_p=None):
    """
    Search for j=0 prime-order curves in [2^bit_lo, 2^bit_hi] with
    lambda/n close to each target_ratio +/- bin_width/2.

    Returns dict: ratio_bin -> list of (p, b, n, lam, G, lam/n)
    """
    BIN = 0.025   # half-width of each target bin
    buckets = {r: [] for r in target_ratios}
    found_total = 0

    p = sympy.nextprime(2**bit_lo - 1)
    p_limit = min(2**bit_hi, max_p) if max_p else 2**bit_hi

    while p < p_limit:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 100: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    for r in target_ratios:
                        if abs(ratio - r) < BIN:
                            if len(buckets[r]) < per_ratio:
                                # find curve parameter b
                                found_b = None
                                for b_try in range(1, min(100, p)):
                                    P = None
                                    rng_tmp = random.Random(7)
                                    for _ in range(200):
                                        x = rng_tmp.randint(0, p - 1)
                                        rhs = (pow(x, 3, p) + b_try) % p
                                        y = tonelli_shanks(rhs, p)
                                        if y is not None and y != 0:
                                            P = (x, y)
                                            break
                                    if P is None: continue
                                    if ec_mul(P, n_cand, p) is None:
                                        found_b = b_try
                                        break
                                if found_b is None: continue
                                G = find_generator(p, found_b, n_cand)
                                if G is None: continue
                                buckets[r].append((p, found_b, n_cand, lam, G, ratio))
                                found_total += 1
                                break  # only add to one bucket per (p, trace)
        p = sympy.nextprime(p)

    return buckets, found_total

# ---------------------------------------------------------------------------
# Signature generation
# ---------------------------------------------------------------------------

def gen_sigs(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs = []
    for _ in range(200000):
        if len(sigs) >= m: break
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

# ---------------------------------------------------------------------------
# GLV lattice + LLL
# ---------------------------------------------------------------------------

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
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

def try_lll(curve, m, d, k1_bound, seed):
    p, b, n, lam, G, ratio = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_sigs(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) == S_KANNAN:
            sign = 1 if last > 0 else -1
            d_cand = (sign * A[i][m]) % n
            if d_cand == d:
                return True
    return False

# ---------------------------------------------------------------------------
# Main experiment: sweep lambda/n and record success at fixed m and seeds
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999]
FIXED_EFF = 0.05   # target effective bias per sig: K1*K2/n ≈ FIXED_EFF

def test_curve(curve, m_max=14):
    """Return first m where all 3 seeds succeed, or None if not in [3, m_max]."""
    p, b, n, lam, G, ratio = curve
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(FIXED_EFF * n / k2_bound))

    eff = k1_bound * k2_bound / n
    if eff >= 1.0:
        m_thresh = 2
    else:
        m_thresh = max(2, math.ceil(math.log(n) / math.log(1.0 / eff)))

    d = random.Random(7).randint(1, n - 1)

    for m in range(max(2, m_thresh - 1), m_max + 1):
        wins = 0
        for seed in SEEDS:
            if try_lll(curve, m, d, k1_bound, seed):
                wins += 1
        if wins == len(SEEDS):
            return m, k1_bound, eff, m_thresh
    return None, k1_bound, eff, m_thresh

# ---------------------------------------------------------------------------
# Hard-coded known-fail curve (p=2677, lam/n=0.07) for reference
# ---------------------------------------------------------------------------

def make_reference_fail_curve():
    p, b, n, lam = 2677, 2, 2647, 185
    G = find_generator(p, b, n)
    return (p, b, n, lam, G, lam / n)

def make_reference_success_curves():
    """Return the three known-success curves from prior runs."""
    curves = []
    # 8-bit: p=211, n=199, lam=106, lam/n=0.53
    G1 = find_generator(211, 2, 199)
    if G1: curves.append((211, 2, 199, 106, G1, 106/199))
    # 12-bit: p=2557, n=2659, lam=1755, lam/n=0.66
    G2 = find_generator(2557, 2, 2659)
    if G2: curves.append((2557, 2, 2659, 1755, G2, 1755/2659))
    return curves

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: lambda/n threshold bisection study")
print(f"Fixed effective bias: eff = K1*K2/n ≈ {FIXED_EFF}")
print(f"Seeds: {SEEDS}")
print("=" * 70)

# First, test the known-fail curve and known-success curves to calibrate
print("\n--- Calibration: known-fail (lam/n=0.07) ---")
fail_curve = make_reference_fail_curve()
if fail_curve[4] is not None:
    m_first, k1_b, eff, m_thr = test_curve(fail_curve, m_max=12)
    status = f"3/3 at m={m_first}" if m_first else "NEVER (m<=12)"
    print(f"  p=2677, n=2647, lam/n={fail_curve[5]:.3f}: LLL {status}"
          f"  (K1={k1_b}, eff={eff:.3f}, m_thresh={m_thr})")
else:
    print("  ERROR: could not build fail curve")

print("\n--- Calibration: known-success curves ---")
for c in make_reference_success_curves():
    m_first, k1_b, eff, m_thr = test_curve(c, m_max=12)
    status = f"3/3 at m={m_first}" if m_first else "NEVER (m<=12)"
    print(f"  p={c[0]}, n={c[2]}, lam/n={c[5]:.3f}: LLL {status}"
          f"  (K1={k1_b}, eff={eff:.3f}, m_thresh={m_thr})")

# Target lambda/n values for the sweep: 0.08 to 0.45 in steps of 0.05
TARGET_RATIOS = [0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]

print(f"\n--- Searching for curves in 12-16 bit range ---")
print(f"    Target lambda/n values: {TARGET_RATIOS}")
print(f"    Bin half-width: 0.025 (so range is target ± 0.025)")
print(f"    Curves per target: 3")
print()

# Search over 12-16 bit range; use up to 2^14 to keep search fast
buckets, found_total = search_curves(
    bit_lo=12, bit_hi=15, target_ratios=TARGET_RATIOS, per_ratio=3, max_p=2**15
)
print(f"Search complete: {found_total} curves found across {len(TARGET_RATIOS)} targets")
print()

# Test each found curve
results = {}   # ratio -> list of (lam_over_n, m_first_3of3 or None)
all_ratios_found = []

print("--- Sweep results ---")
print(f"{'lam/n':>10}  {'p':>8}  {'n':>8}  {'K1':>4}  {'eff':>6}  {'m_thresh':>9}  {'LLL result':>18}")
print("-" * 75)

for target_r in sorted(buckets.keys()):
    curves = buckets[target_r]
    for c in curves:
        p, b, n, lam, G, ratio = c
        m_first, k1_b, eff, m_thr = test_curve(c, m_max=14)
        status = f"3/3 at m={m_first}" if m_first else "FAIL (m<=14)"
        print(f"{ratio:>10.4f}  {p:>8}  {n:>8}  {k1_b:>4}  {eff:>6.3f}  {m_thr:>9}  {status:>18}")
        if target_r not in results:
            results[target_r] = []
        results[target_r].append((ratio, m_first))
        all_ratios_found.append((ratio, m_first, p, n))

print()
print("=" * 70)
print("THRESHOLD SUMMARY")
print("=" * 70)
print(f"\n{'lam/n range':>20}  {'outcome':>25}  {'notes'}")
print("-" * 70)

# Determine transition
all_tested = sorted(all_ratios_found, key=lambda x: x[0])
success_ratios = [(r, m) for r, m, p, n in all_tested if m is not None]
fail_ratios    = [(r, m) for r, m, p, n in all_tested if m is None]

if fail_ratios:
    max_fail = max(r for r, m in fail_ratios)
    print(f"  Max FAIL ratio: {max_fail:.4f}")
if success_ratios:
    min_success = min(r for r, m in success_ratios)
    print(f"  Min SUCCESS ratio: {min_success:.4f}")

if fail_ratios and success_ratios:
    max_f = max(r for r, m in fail_ratios)
    min_s = min(r for r, m in success_ratios)
    if min_s > max_f:
        print(f"\n  THRESHOLD BRACKET: ({max_f:.4f}, {min_s:.4f})")
        print(f"  Midpoint estimate: tau ≈ {(max_f + min_s) / 2:.4f}")
    else:
        print(f"\n  Transition is non-monotone or overlapping; bracket unclear.")

print("\n  Full data:")
for r, m, p, n in all_tested:
    status = f"3/3 at m={m}" if m else "FAIL"
    print(f"    lam/n={r:.4f} (p={p}, n={n}): {status}")

# Reference values from prior runs
print("\n  Reference points (from prior runs, different bit lengths):")
print(f"    lam/n=0.070 (n=2647,   12b): FAIL (LLL + BKZ(40))")
print(f"    lam/n=0.340 (n=523969, 19b): SUCCESS (3/3 at m=9)")
print(f"    lam/n=0.530 (n=199,     8b): SUCCESS (3/3 at m=4)")
print(f"    lam/n=0.660 (n=2659,   12b): SUCCESS (3/3 at m=7)")

print("\nDone.")
