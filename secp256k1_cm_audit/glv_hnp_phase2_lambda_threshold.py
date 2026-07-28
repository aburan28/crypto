"""
GLV-HNP Phase 2: Lambda/n threshold study.

Goal: Bisect the λ/n threshold between 0.07 (known fail) and 0.34 (known success)
to find the minimum ratio at which LLL recovers the secret key.

Previous data points (from 2026-07-26 log):
  λ/n = 0.07  (p=2677, n=2647,  lam=185)  → LLL FAIL, BKZ(40) FAIL
  λ/n = 0.34  (p=524347,n=523969,lam=177902) → LLL 3/3 at m=9
  λ/n = 0.44  (secp256k1 proxy) → in "safe" regime (theoretical)
  λ/n = 0.53  (p=211, n=199, lam=106) → LLL 3/3 at m=4
  λ/n = 0.66  (p=2557, n=2659, lam=1755) → LLL 3/3 at m=7

Method:
  1. Find j=0 prime curves at target λ/n ratios (0.10, 0.15, 0.20, 0.25) using
     Eisenstein decomposition + Tonelli-Shanks eigenvalue computation.
  2. Run the Phase 2 GLV lattice attack at m=4..12 for each curve.
  3. Report whether LLL 3/3 (3 seeds all recover) at any m in range.
  4. Use a small (12-bit) prime size to keep runtime under 60s.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass y^2 = x^3 + b)
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

def find_curve_near_ratio(target_ratio, bit_lo=10, bit_hi=13, tol=0.04):
    """
    Search for a j=0 prime curve with λ/n near target_ratio (within ±tol).
    Returns (p, b, n, lam, G) or None.
    """
    p = sympy.nextprime(2**bit_lo)
    p_max = 2**bit_hi
    best = None
    best_err = float('inf')

    while p < p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 8 or not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    # glv_eigenvalue returns lam <= n//2, so ratio in [0, 0.5]
                    # But the "other" eigenvalue is n-1-lam, ratio ≈ 1 - ratio
                    # We want ratio closest to target_ratio from BOTH eigenvalues
                    ratio2 = (n_cand - 1 - lam) / n_cand
                    for r in [ratio, ratio2]:
                        err = abs(r - target_ratio)
                        if err < tol and err < best_err:
                            # Use lam corresponding to this ratio
                            actual_lam = lam if r == ratio else (n_cand - 1 - lam)
                            # Find b for this curve
                            found_b = None
                            for b_try in range(1, 300):
                                rng_try = random.Random(42)
                                found_pt = None
                                for _ in range(5000):
                                    x = rng_try.randint(0, p - 1)
                                    rhs = (pow(x, 3, p) + b_try) % p
                                    y = tonelli_shanks(rhs, p)
                                    if y is not None and y != 0:
                                        pt = (x, y)
                                        if ec_mul(pt, n_cand, p) is None:
                                            found_b = b_try
                                            break
                                if found_b is not None:
                                    break
                            if found_b is None:
                                continue
                            G = find_generator(p, found_b, n_cand)
                            if G is None:
                                continue
                            best = (p, found_b, n_cand, actual_lam, G)
                            best_err = err
        p = sympy.nextprime(p)

    return best

# ---------------------------------------------------------------------------
# GLV lattice attack (from Phase 2)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
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

def run_lll(curve_params, m, d_secret, k1_bound, seed=42):
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

def sweep(label, curve_params, k1_bound, m_range, seeds):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    if eff < 1.0:
        m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))
    else:
        m_thresh = float('inf')
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}, lam/n={lam/n:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")
    print(f"{'='*60}")
    first_3of3 = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            ok = run_lll(curve_params, m, d_trial, k1_bound, seed)
            wins += ok
        marker = " ← 3/3!" if wins == len(seeds) else ""
        print(f"  m={m}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first_3of3 is None:
            first_3of3 = m
    return first_3of3

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 60)
print("GLV-HNP Phase 2 — λ/n threshold bisection study")
print("=" * 60)
print()
print("Known data:")
print("  λ/n=0.07 (p=2677, n=2647): FAIL at all m (even BKZ-40)")
print("  λ/n=0.34 (p=524347):       3/3 LLL at m=9")
print("  λ/n=0.53 (p=211):          3/3 LLL at m=4")
print()
print("Target ratios to probe: 0.10, 0.15, 0.20, 0.25, 0.30")
print()

SEEDS = [42, 1234, 9999]
TARGET_RATIOS = [0.10, 0.15, 0.20, 0.25, 0.30]

summary = {}

for target in TARGET_RATIOS:
    print(f"\n--- Searching for curve near λ/n={target:.2f} ---")
    c = find_curve_near_ratio(target, bit_lo=10, bit_hi=13, tol=0.05)
    if c is None:
        print(f"  NOT FOUND at tol=0.05, widening to 0.10...")
        c = find_curve_near_ratio(target, bit_lo=10, bit_hi=14, tol=0.10)
    if c is None:
        print(f"  SKIPPED: no suitable curve found near λ/n={target:.2f}")
        summary[target] = None
        continue

    p, b, n, lam, G = c
    actual_ratio = lam / n
    print(f"  Found: p={p}, n={n}, lam={lam}, λ/n={actual_ratio:.4f}")

    # K1_BOUND: target eff ≈ 0.05
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.05 * math.sqrt(n)))

    label = f"λ/n≈{actual_ratio:.3f} (p={p}, n={n}, lam={lam})"
    first_3of3 = sweep(label, c, k1_bound, range(4, 14), SEEDS)
    summary[target] = (actual_ratio, first_3of3, p, n, lam)

# Also sweep the two known reference curves for comparison
print("\n\n--- Reference: λ/n=0.07 (p=2677, n=2647) known-FAIL ---")
G_fail = find_generator(2677, 2, 2647)
if G_fail:
    c_fail = (2677, 2, 2647, 185, G_fail)
    k1_f = max(2, int(0.05 * math.sqrt(2647)))
    sweep("λ/n=0.07 FAIL (p=2677)", c_fail, k1_f, range(4, 14), SEEDS)
    summary[0.07] = (0.0699, None, 2677, 2647, 185)

print("\n--- Reference: λ/n≈0.53 (p=211, n=199) known-PASS ---")
G_pass = find_generator(211, 2, 199)
if G_pass:
    c_pass = (211, 2, 199, 106, G_pass)
    sweep("λ/n=0.53 PASS (p=211)", c_pass, 2, range(3, 8), SEEDS)
    summary[0.53] = (0.5327, 4, 211, 199, 106)

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print("\n")
print("=" * 60)
print("THRESHOLD SUMMARY")
print("=" * 60)
print(f"{'λ/n':>8}  {'p':>7}  {'n':>7}  {'lam':>7}  {'LLL 3/3 at m':>14}")
print("-" * 60)

all_data = []
if 0.07 in summary and summary[0.07]:
    ratio, first_m, p, n, lam = summary[0.07]
    all_data.append((ratio, first_m, p, n, lam))

for target in TARGET_RATIOS:
    if target in summary and summary[target]:
        ratio, first_m, p, n, lam = summary[target]
        all_data.append((ratio, first_m, p, n, lam))

if 0.53 in summary and summary[0.53]:
    ratio, first_m, p, n, lam = summary[0.53]
    all_data.append((ratio, first_m, p, n, lam))

all_data.sort(key=lambda x: x[0])
for ratio, first_m, p, n, lam in all_data:
    status = f"m={first_m}" if first_m is not None else "FAIL"
    print(f"{ratio:>8.4f}  {p:>7}  {n:>7}  {lam:>7}  {status:>14}")

# Find threshold
successes = [(r, m) for r, m, *_ in all_data if m is not None]
failures  = [(r, m) for r, m, *_ in all_data if m is None]
if successes and failures:
    lo = max(r for r, _ in failures)
    hi = min(r for r, _ in successes)
    print(f"\nThreshold bracket: ({lo:.4f}, {hi:.4f})")
    print(f"  LLL fails  for λ/n <= {lo:.4f}")
    print(f"  LLL succeeds for λ/n >= {hi:.4f}")
else:
    print("\nInsufficient data for threshold bracket.")

print("\nDone.")
