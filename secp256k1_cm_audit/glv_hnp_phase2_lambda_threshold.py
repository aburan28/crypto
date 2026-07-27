"""
GLV-HNP Phase 2: lambda/n threshold bisection study.

Goal: Find the threshold lam/n* such that the GLV-aware LLL attack succeeds
for lam/n > lam/n* and fails for lam/n < lam/n*.

Known bounds from 2026-07-26 log:
  - lam/n = 0.07 (p=2677, n=2647): LLL fails, BKZ(beta=40) fails
  - lam/n = 0.34 (p=524347, n=523969): LLL 3/3 at m=9

This script sweeps over j=0 CM curves with various lam/n values in [0.05, 0.50]
and records the minimum m for 3/3 LLL success (or "fail" up to m_max=20).

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ


# ---------------------------------------------------------------------------
# Pure-Python EC arithmetic (short Weierstrass a=0)
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


# ---------------------------------------------------------------------------
# Eisenstein CM machinery (from glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - a*b + b^2 = p, a,b >= 0. O(sqrt(p))."""
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
    """Six sextic-twist Frobenius traces for p = a^2-ab+b^2."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue_small(n):
    """
    Smaller root of x^2+x+1 = 0 mod n (needs n prime, n ≡ 1 mod 3).
    Returns lam_small in [1, n//2], or None if no root exists.
    """
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    # pick the one that satisfies the equation
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)


# ---------------------------------------------------------------------------
# Find twist coefficient b for a given (p, n_target)
# ---------------------------------------------------------------------------

def count_points_on_curve(p, b):
    """Count #E(F_p) for y^2 = x^3 + b. O(p)."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (pow(x, 3, p) + b) % p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2
    return count

def find_twist_b(p, n_target, b_limit=300):
    """Find b in [1, b_limit] such that #(y^2=x^3+b over F_p) = n_target."""
    for b in range(1, b_limit + 1):
        if count_points_on_curve(p, b) == n_target:
            return b
    return None

def find_generator(p, b, n):
    """Find generator of E(F_p): y^2=x^3+b of prime order n."""
    rng = random.Random(99991)
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
# Curve database: search for representatives at each lam/n bucket
# ---------------------------------------------------------------------------

def build_curve_database(p_max=8000, ratio_bins=None):
    """
    Search primes p ≡ 1 (mod 3) up to p_max.
    For each p, compute all 6 group orders n_k via Eisenstein decomposition.
    For prime n_k ≡ 1 (mod 3), compute lam/n and bin.
    Returns dict: ratio_bucket -> (p, n, lam) list.
    """
    if ratio_bins is None:
        # Bins of width 0.05 covering (0, 0.50)
        ratio_bins = [round(0.05 * i, 2) for i in range(1, 10)]
        # Also add finer bins in the critical zone [0.07, 0.34]
        ratio_bins = sorted(set(ratio_bins + [0.07, 0.10, 0.13, 0.16, 0.20, 0.24, 0.28, 0.32]))

    db = {rb: [] for rb in ratio_bins}

    for p in sympy.primerange(100, p_max):
        if p % 3 != 1:
            continue
        eis = eisenstein_decompose(p)
        if eis is None:
            continue
        a_e, b_e = eis
        traces = j0_traces(a_e, b_e)
        for t in traces:
            n_k = p + 1 - t
            if n_k < 50 or n_k == p:
                continue
            if not sympy.isprime(n_k):
                continue
            if n_k % 3 != 1:
                continue
            lam = glv_eigenvalue_small(n_k)
            if lam is None:
                continue
            ratio = lam / n_k
            # Find closest bin
            best_bin = min(ratio_bins, key=lambda rb: abs(ratio - rb))
            if abs(ratio - best_bin) < 0.025:  # within 0.025 of bin center
                db[best_bin].append((p, n_k, lam, ratio))

    # Sort each bin by n (ascending = faster to attack)
    for rb in ratio_bins:
        db[rb].sort(key=lambda x: x[1])

    return db, ratio_bins


# ---------------------------------------------------------------------------
# LLL/BKZ attack machinery (from glv_hnp_phase2_attack.py)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def run_lll_attack(curve, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def find_min_m(curve, k1_bound, seeds, m_range, use_bkz=False, bkz_beta=20):
    """
    Find the minimum m in m_range such that LLL (or BKZ) recovers d
    for ALL seeds. Returns m or None if never 3/3 in range.
    """
    p, b, n, lam, G = curve
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_secret = random.Random(seed + 3141).randint(1, n - 1)
            ok = run_lll_attack(curve, m, d_secret, k1_bound, seed,
                                use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        if wins == len(seeds):
            return m
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: lambda/n threshold bisection study")
print("=" * 70)

# Seeds and K1_BOUND
SEEDS = [42, 1234, 9999]
K1_BOUND = 2   # bias: k1 in {0, 1} — same as 8-bit reference
M_RANGE = range(3, 21)  # sweep m = 3..20

# ---- Step 1: Build curve database -----------------------------------------
print("\n[Step 1] Searching for j=0 CM curves with diverse lam/n ratios...")
print("  (primes p ≡ 1 mod 3 up to 8000, prime group orders n ≡ 1 mod 3)")

TARGET_BINS = sorted(set([0.05, 0.07, 0.10, 0.13, 0.16, 0.20, 0.24, 0.28, 0.32,
                           0.35, 0.40, 0.45, 0.48]))
db, bins = build_curve_database(p_max=8000, ratio_bins=TARGET_BINS)

print(f"\n  Candidates per bin:")
for rb in bins:
    n_found = len(db[rb])
    if n_found > 0:
        best = db[rb][0]
        print(f"    bin {rb:.2f}: {n_found} curves; smallest p={best[0]}, n={best[1]}, lam={best[2]}, ratio={best[3]:.4f}")
    else:
        print(f"    bin {rb:.2f}: (none found)")

# ---- Step 2: Select representative per bin and find twist b ---------------
print("\n[Step 2] Selecting representative per bin and finding twist coefficient b...")

selected = []  # list of (ratio_bin, p, b, n, lam, actual_ratio)
for rb in bins:
    if not db[rb]:
        continue
    # Pick smallest n for speed
    for (p, n, lam, ratio) in db[rb]:
        if n < 150:  # too small
            continue
        print(f"  bin {rb:.2f}: trying p={p}, n={n}, lam={lam} (ratio={ratio:.4f})... ", end='', flush=True)
        b_coeff = find_twist_b(p, n, b_limit=200)
        if b_coeff is None:
            print("b not found, trying next")
            continue
        G = find_generator(p, b_coeff, n)
        if G is None:
            print("generator not found, trying next")
            continue
        # Verify lambda
        assert (lam * lam + lam + 1) % n == 0, f"lambda check failed n={n}"
        print(f"b={b_coeff}, G={G}")
        selected.append((rb, p, b_coeff, n, lam, ratio))
        break

print(f"\n  Selected {len(selected)} curves for LLL sweep.")

# ---- Step 3: LLL sweep for each curve -------------------------------------
print("\n[Step 3] LLL sweep: find min m for 3/3 success at each lam/n ratio...")
print("  (K1_BOUND=2, seeds=42/1234/9999, m in 3..20)")
print()

results_lll = []
for (rb, p, b, n, lam, ratio) in selected:
    curve = (p, b, n, lam, find_generator(p, b, n))
    k2_bound = math.isqrt(n) + 1
    eff = K1_BOUND * k2_bound / n
    m_thresh_theory = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')

    # Quick per-m sweep with status
    min_m = None
    row_str = []
    for m in M_RANGE:
        wins = 0
        for seed in SEEDS:
            d_secret = random.Random(seed + 3141).randint(1, n - 1)
            ok = run_lll_attack(curve, m, d_secret, K1_BOUND, seed)
            wins += ok
        row_str.append(f"m{m}:{wins}/{len(SEEDS)}")
        if wins == len(SEEDS) and min_m is None:
            min_m = m
            break  # found the minimum

    status = f"3/3 at m={min_m}" if min_m is not None else f"FAIL (up to m={max(M_RANGE)-1})"
    print(f"  ratio={ratio:.4f} (bin {rb:.2f})  p={p}  n={n}  lam={lam}")
    print(f"    K2={k2_bound}  eff={eff:.4f}  m_thresh≈{m_thresh_theory:.1f}")
    print(f"    Sweep: {', '.join(row_str)}")
    print(f"    LLL result: {status}")
    print()

    results_lll.append({
        'ratio': ratio,
        'bin': rb,
        'p': p,
        'n': n,
        'lam': lam,
        'min_m_lll': min_m,
        'eff': eff,
        'm_thresh_theory': m_thresh_theory,
    })

# ---- Step 4: BKZ test on failing curves -----------------------------------
print("\n[Step 4] BKZ rescue test on curves where LLL fails...")

fail_curves = [(rb, p, b, n, lam, ratio) for (rb, p, b, n, lam, ratio) in selected
               if next((r['min_m_lll'] for r in results_lll
                        if r['p'] == p and r['n'] == n), None) is None]

if not fail_curves:
    print("  No LLL-failing curves — threshold may be below smallest tested ratio.")
else:
    for (rb, p, b, n, lam, ratio) in fail_curves[:3]:  # test at most 3
        curve = (p, b, n, find_generator(p, b, n))
        print(f"  ratio={ratio:.4f}  n={n}  lam={lam}")
        for bkz_beta in [20, 40, 60]:
            min_m_bkz = None
            for m in M_RANGE:
                wins = 0
                for seed in SEEDS:
                    d_secret = random.Random(seed + 3141).randint(1, n - 1)
                    c_full = (p, b, n, lam, find_generator(p, b, n))
                    ok = run_lll_attack(c_full, m, d_secret, K1_BOUND, seed,
                                       use_bkz=True, bkz_beta=bkz_beta)
                    wins += ok
                if wins == len(SEEDS):
                    min_m_bkz = m
                    break
            status = f"3/3 at m={min_m_bkz}" if min_m_bkz else "FAIL"
            print(f"    BKZ(beta={bkz_beta}): {status}")
        print()

# ---- Summary table ----------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY: lam/n threshold study")
print("=" * 70)
print(f"{'lam/n':>8}  {'n':>6}  {'lam':>6}  {'eff':>7}  {'m_thresh':>8}  {'LLL result':>15}")
print("-" * 70)

threshold_lo = None  # highest ratio where LLL fails
threshold_hi = None  # lowest ratio where LLL succeeds

for r in sorted(results_lll, key=lambda x: x['ratio']):
    m_result = r['min_m_lll']
    if m_result is None:
        lll_str = "FAIL"
        threshold_lo = max(threshold_lo or 0, r['ratio'])
    else:
        lll_str = f"3/3 at m={m_result}"
        if threshold_hi is None or r['ratio'] < threshold_hi:
            threshold_hi = r['ratio']
    print(f"{r['ratio']:>8.4f}  {r['n']:>6}  {r['lam']:>6}  {r['eff']:>7.4f}"
          f"  {r['m_thresh_theory']:>8.1f}  {lll_str:>15}")

print("-" * 70)
if threshold_lo is not None and threshold_hi is not None:
    print(f"\nThreshold bracket: lam/n ∈ ({threshold_lo:.4f}, {threshold_hi:.4f})")
    print(f"Best estimate: lam/n* ≈ {(threshold_lo + threshold_hi)/2:.4f}")
elif threshold_lo is None:
    print(f"\nLLL succeeds at ALL tested ratios (lam/n >= {min(r['ratio'] for r in results_lll):.4f})")
    print("Threshold may be below smallest tested ratio. Try curves with lam/n < 0.07.")
elif threshold_hi is None:
    print(f"\nLLL FAILS at ALL tested ratios (up to lam/n = {max(r['ratio'] for r in results_lll):.4f})")

print("\nDone.")
