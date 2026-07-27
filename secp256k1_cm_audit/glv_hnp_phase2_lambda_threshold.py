"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Prior data (from 2026-07-26 autolab run):
  λ/n = 0.07  (p=2677,   n=2647):   LLL/BKZ(40) FAIL — small-λ failure
  λ/n = 0.34  (p=524347, n=523969): LLL 3/3 at m=9
  λ/n = 0.44  (secp256k1 proxy):    expected success
  λ/n = 0.53  (p=211, n=199):       LLL 3/3 at m=4
  λ/n = 0.66  (p=2557, n=2659):     LLL 3/3 at m=7

Goal: find curves covering λ/n ∈ [0.07, 0.40] and determine the minimum viable λ/n
where LLL can recover d (with a reasonable m budget, say m ≤ 15).

Method:
  Enumerate small primes p ≡ 1 (mod 3), use Eisenstein CM decomposition to compute
  all 6 traces, derive λ for each n = p+1-t, and collect curves in λ/n bins.
  For each candidate curve run LLL sweep m=4..12 with 3 seeds and BKZ(40) as fallback.

Run: python3 glv_hnp_phase2_lambda_threshold.py
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
    m_v, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_v - i - 1), p)
        m_v, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(100000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory
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
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

def find_curve_b(p, n):
    """Find b such that #E(y^2=x^3+b over F_p) = n."""
    rng = random.Random(999)
    for _ in range(200):
        b = rng.randint(1, p - 1)
        P = None
        for x in range(min(p, 500)):
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                break
        if P is None: continue
        if ec_mul(P, n, p) is None:
            return b
    # systematic search
    for b in range(1, min(p, 10000)):
        for x in range(min(p, 1000)):
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b
                break
    return None

# ---------------------------------------------------------------------------
# Curve discovery: find representative curves per λ/n bin
# ---------------------------------------------------------------------------

def discover_curves(target_ratios, p_max=200000):
    """
    For each target ratio in target_ratios (e.g. 0.10, 0.15, ..., 0.35),
    find the first small prime-order j=0 curve whose λ/n is within 0.025 of target.

    Returns dict: {ratio_bin: (p, b, n, lam, G)}
    """
    found = {}  # ratio_bin -> curve
    remaining = set(target_ratios)

    p = 7
    while p < p_max and remaining:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                traces = j0_traces(a_e, b_e)
                for t in traces:
                    n_cand = p + 1 - t
                    if n_cand < 10: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand

                    # Check which bins this could fill
                    for target in list(remaining):
                        if abs(ratio - target) < 0.025:
                            b_param = find_curve_b(p, n_cand)
                            if b_param is None: continue
                            G = find_generator(p, b_param, n_cand)
                            if G is None: continue
                            found[target] = (p, b_param, n_cand, lam, G)
                            remaining.discard(target)
                            print(f"  [bin {target:.2f}] p={p}, n={n_cand}, λ={lam}, λ/n={ratio:.4f}")
                            break
        p = sympy.nextprime(p)

    if remaining:
        print(f"  WARNING: could not find curves for bins: {sorted(remaining)}")
    return found

# ---------------------------------------------------------------------------
# GLV lattice / signatures / recovery
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
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return True
    return False

def run_experiment(curve_params, m, d_secret, k1_bound, seed=42,
                   use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def sweep_curve(label, curve_params, k1_bound, m_range, seeds,
                use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')

    algo = f"BKZ(beta={bkz_beta})" if use_bkz else "LLL"
    print(f"\n  [{algo}] {label}")
    print(f"    n={n} ({n.bit_length()}b), lam={lam}, lam/n={lam/n:.4f}, "
          f"K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")

    results = {}
    first_full = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            ok = run_experiment(curve_params, m, d_trial, k1_bound, seed,
                                use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
        marker = "*" if m >= m_thresh else " "
        print(f"    {marker}m={m}: {wins}/{len(seeds)}")
        if wins == len(seeds) and first_full is None:
            first_full = m
    return results, first_full

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection study")
print("=" * 70)
print()

SEEDS = [42, 1234, 9999]

# Target λ/n bins: from 0.08 up to 0.38 in steps of 0.03
# (bracketing the unknown threshold between 0.07 [fail] and 0.34 [pass])
TARGET_RATIOS = [0.10, 0.13, 0.16, 0.20, 0.24, 0.28, 0.32]

print("Step 1: Discover one representative curve per λ/n bin")
print("-" * 70)
curves_by_bin = discover_curves(TARGET_RATIOS, p_max=300000)

# Also include the two known data points
print(f"\nKnown reference curves:")
print(f"  [fail] p=2677, n=2647, lam=185, lam/n=0.0699  (BKZ(40) fails)")
print(f"  [pass] p=524347, n=523969, lam=177902, lam/n=0.3396  (LLL 3/3 m=9)")

print()
print("Step 2: LLL sweep for each discovered curve")
print("=" * 70)

# k1_bound target: want eff ≈ 0.05 (balance between enough signatures and useful bias)
# eff = K1 * sqrt(n) / n = K1 / sqrt(n)
# K1 ≈ 0.05 * sqrt(n) but at least 2

summary = []  # [(lam_ratio, label, lll_first_full, bkz_first_full)]

for target in sorted(curves_by_bin.keys()):
    curve = curves_by_bin[target]
    p, b, n, lam, G = curve
    lam_ratio = lam / n

    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    k2_bound = math.isqrt(n) + 1

    label = f"bin {target:.2f}: p={p}, n={n}, lam/n={lam_ratio:.4f}"
    print(f"\n{'='*70}")
    print(f"Curve: {label}")

    # LLL sweep m=4..14
    m_range = range(4, 15)
    res_lll, first_lll = sweep_curve(label, curve, k1_bound, m_range, SEEDS,
                                     use_bkz=False)

    # If LLL failed, try BKZ(40)
    first_bkz = None
    if first_lll is None:
        print(f"\n  LLL failed all m — trying BKZ(40)")
        res_bkz, first_bkz = sweep_curve(label, curve, k1_bound, m_range, SEEDS,
                                          use_bkz=True, bkz_beta=40)

    summary.append((lam_ratio, target, label, first_lll, first_bkz))

print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'λ/n':>8}  {'LLL first 3/3':>14}  {'BKZ(40) first 3/3':>18}  {'Status'}")
print("-" * 70)

# Include known endpoints
print(f"{'0.0699':>8}  {'never':>14}  {'never':>18}  FAIL (known)")
for lam_ratio, target, label, first_lll, first_bkz in sorted(summary):
    lll_str = f"m={first_lll}" if first_lll is not None else "never"
    bkz_str = f"m={first_bkz}" if first_bkz is not None else (
        "n/a" if first_lll is not None else "never"
    )
    status = "PASS" if first_lll is not None else (
        "PASS(BKZ)" if first_bkz is not None else "FAIL"
    )
    print(f"{lam_ratio:>8.4f}  {lll_str:>14}  {bkz_str:>18}  {status}")
print(f"{'0.3396':>8}  {'m=9 (known)':>14}  {'n/a':>18}  PASS (known)")
print(f"{'0.5302':>8}  {'m=4 (known)':>14}  {'n/a':>18}  PASS (known)")

print()
print("Threshold interpretation:")
# Find first PASS in sorted summary
threshold_upper = None
threshold_lower = 0.07
for lam_ratio, target, label, first_lll, first_bkz in sorted(summary):
    if first_lll is not None or first_bkz is not None:
        threshold_upper = lam_ratio
        break
    else:
        threshold_lower = lam_ratio

if threshold_upper is not None:
    print(f"  LLL threshold bracketed: ({threshold_lower:.4f}, {threshold_upper:.4f})")
    if threshold_lower == 0.07:
        print(f"  Lower bracket from known failure at λ/n = 0.0699")
else:
    print(f"  LLL still failing at all tested ratios up to {max(t for t,_,_,_,_ in summary):.4f}")
    print(f"  Need to test higher λ/n or increase m budget.")

print("\nDone.")
