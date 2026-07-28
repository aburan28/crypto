"""
GLV-HNP Phase 2: lambda/n threshold bisection.

Prior results:
  lam/n = 0.53 (8-bit/199):   LLL 3/3 at m=4
  lam/n = 0.66 (12-bit/2659): LLL 3/3 at m=7
  lam/n = 0.34 (20-bit):      LLL 3/3 at m=9
  lam/n = 0.07 (12-bit/2647): LLL never 3/3 (all m<=12, all BKZ betas)

Goal: bisect the failure threshold between 0.07 and 0.34 by finding
j=0 curves with lam/n ∈ {0.10, 0.13, 0.17, 0.20, 0.25, 0.30} and
testing LLL recovery across m=4..14 with 3 seeds each.

Strategy: sweep 12-20 bit primes p≡1(mod 3), Eisenstein decompose,
compute all 6 twist traces, keep candidates with the desired lam/n bucket.
Then run the Phase 2 GLV attack.

Run: python3 glv_hnp_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass, a=0)
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
# CM / Eisenstein utilities
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

def find_curve_for_b_in_range(p, n, lam_min_frac, lam_max_frac):
    """Find b for curve y^2=x^3+b over F_p with order n (trial b=1..200)."""
    for b in range(1, 200):
        rhs = (pow(1, 3, p) + b) % p  # x=1
        # Quick check: try b until we find a point of correct order
        rng = random.Random(b * 7)
        for _ in range(20):
            x = rng.randint(0, p - 1)
            rhs_x = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs_x, p)
            if y is not None and y != 0:
                P = (x, y)
                Q = ec_mul(P, n, p)
                if Q is None:
                    return b
                break
    return None

# ---------------------------------------------------------------------------
# Build & attack lattice
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
            return d_cand
    return None

def run_lll(curve_params, m, d_secret, k1_bound, seed=42):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep_m(label, curve_params, k1_bound, m_range, seeds, verbose=True):
    p, b, n, lam, G = curve_params
    ratio = lam / n
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    if verbose:
        print(f"\n  [{label}]  lam/n={ratio:.4f}  n={n}({n.bit_length()}b)  K1={k1_bound}  eff={eff:.4f}")
    first_full = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 3333).randint(1, n - 1)
            if run_lll(curve_params, m, d_trial, k1_bound, seed):
                wins += 1
        result = f"{wins}/{len(seeds)}"
        marker = " *** THRESHOLD" if wins == len(seeds) and first_full is None else ""
        if verbose:
            print(f"    m={m:2d}: {result}{marker}")
        if wins == len(seeds) and first_full is None:
            first_full = m
    return first_full, ratio

# ---------------------------------------------------------------------------
# Curve search: find j=0 curves with lam/n in target bucket
# ---------------------------------------------------------------------------

TARGET_BUCKETS = [
    (0.08, 0.12, "bucket_0.10"),
    (0.12, 0.16, "bucket_0.13"),
    (0.16, 0.20, "bucket_0.18"),
    (0.20, 0.25, "bucket_0.22"),
    (0.25, 0.30, "bucket_0.27"),
    (0.30, 0.36, "bucket_0.33"),
]

def find_curves_per_bucket(bit_min=12, bit_max=20, max_per_bucket=1):
    """
    Sweep primes in [2^bit_min, 2^bit_max], Eisenstein-decompose, compute
    traces for all 6 twists, sort candidates into lam/n buckets.
    Returns dict: bucket_label -> list of (p, b, n, lam, G, ratio).
    """
    buckets = {label: [] for (_, _, label) in TARGET_BUCKETS}
    found_count = {label: 0 for (_, _, label) in TARGET_BUCKETS}

    print(f"Searching for j=0 curves in [{bit_min}, {bit_max}]-bit range...")
    p_start = sympy.nextprime(2**(bit_min - 1))
    p_end = 2**bit_max
    p = p_start
    total_checked = 0
    while p < p_end and any(c < max_per_bucket for c in found_count.values()):
        p = sympy.nextprime(p)
        if p % 3 != 1:
            continue
        total_checked += 1
        eis = eisenstein_decompose(p)
        if eis is None:
            continue
        a_e, b_e = eis
        traces = j0_traces(a_e, b_e)
        for t in traces:
            n_cand = p + 1 - t
            if n_cand < 2 or not sympy.isprime(n_cand):
                continue
            if n_cand % 3 != 1:
                continue
            lam, _ = glv_eigenvalue(n_cand)
            if lam is None:
                continue
            ratio = lam / n_cand
            for (lo, hi, label) in TARGET_BUCKETS:
                if lo <= ratio < hi and found_count[label] < max_per_bucket:
                    # Try to find b
                    found_b = None
                    for b_try in range(1, 300):
                        rng_b = random.Random(b_try * 13)
                        ok = False
                        for _ in range(30):
                            x = rng_b.randint(0, p - 1)
                            rhs_x = (pow(x, 3, p) + b_try) % p
                            y_x = tonelli_shanks(rhs_x, p)
                            if y_x is not None and y_x != 0:
                                Pt = (x, y_x)
                                if ec_mul(Pt, n_cand, p) is None:
                                    found_b = b_try
                                    ok = True
                                    break
                            if ok: break
                        if found_b is not None:
                            break
                    if found_b is None:
                        continue
                    G_cand = find_generator(p, found_b, n_cand)
                    if G_cand is None:
                        continue
                    buckets[label].append((p, found_b, n_cand, lam, G_cand, ratio))
                    found_count[label] += 1
                    print(f"  Found {label}: p={p}, n={n_cand}, lam={lam}, ratio={ratio:.4f}")
    print(f"Checked {total_checked} primes p≡1(mod 3).")
    return buckets

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — lambda/n threshold bisection")
print("Prior: lam/n=0.34 succeeds, lam/n=0.07 fails")
print("Goal:  find the critical threshold in (0.07, 0.34)")
print("=" * 70)

SEEDS = [42, 1234, 9999]
M_RANGE = range(4, 16)  # m=4..15

# ---- Step 1: find one curve per bucket ------------------------------------
buckets = find_curves_per_bucket(bit_min=11, bit_max=22, max_per_bucket=1)

# ---- Step 2: for each bucket, run sweep -----------------------------------
print("\n" + "=" * 70)
print("Threshold sweep: LLL m=4..15, 3 seeds")
print("=" * 70)

summary = []  # (ratio, first_full_m, label, curve_params)

# Also include known reference points
# Reference: lam/n=0.07 (known FAIL)
ref_fail = (2677, 2, 2647, 185, find_generator(2677, 2, 2647), 185/2647)
# Reference: lam/n=0.34 (known SUCCESS from 2026-07-26)
# We'll re-verify this with the 20-bit curve found in 20bit.py; use approximate stand-in
# p=524347 b=2 n=523969 lam=177902 from that run
_G_ref20 = find_generator(524347, 2, 523969)
ref_ok = (524347, 2, 523969, 177902, _G_ref20, 177902/523969)

print("\n[REFERENCE: lam/n=0.07, known FAIL]")
first_ref_fail, _ = sweep_m("ref_fail lam/n=0.07", ref_fail[:5], k1_bound=8,
                             m_range=M_RANGE, seeds=SEEDS)
summary.append((0.0701, first_ref_fail, "ref_FAIL_p2677", ref_fail[:5]))

if ref_ok[4] is not None:
    print("\n[REFERENCE: lam/n=0.34, known SUCCESS]")
    k1_ref20 = max(2, int(0.05 * math.sqrt(ref_ok[2])))
    first_ref_ok, _ = sweep_m("ref_ok lam/n=0.34", ref_ok[:5], k1_bound=k1_ref20,
                               m_range=M_RANGE, seeds=SEEDS)
    summary.append((ref_ok[5], first_ref_ok, "ref_OK_p524347", ref_ok[:5]))

for (lo, hi, label) in TARGET_BUCKETS:
    curves = buckets.get(label, [])
    if not curves:
        print(f"\n  [SKIP] {label}: no curve found in search.")
        summary.append(((lo + hi) / 2, None, label, None))
        continue
    p, b, n, lam, G, ratio = curves[0]
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    first_full, _ = sweep_m(label, (p, b, n, lam, G), k1_bound=k1_bound,
                             m_range=M_RANGE, seeds=SEEDS)
    summary.append((ratio, first_full, label, (p, b, n, lam, G)))

# ---- Step 3: Report -------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY TABLE — lambda/n threshold study")
print("=" * 70)
print(f"{'lam/n':>10}  {'m*':>5}  {'label'}")
print("-" * 40)
for ratio, m_star, label, _ in sorted(summary, key=lambda x: x[0]):
    m_str = str(m_star) if m_star is not None else "FAIL"
    print(f"{ratio:10.4f}  {m_str:>5}  {label}")

# Find threshold: last FAIL then first SUCCESS
sorted_summary = sorted(summary, key=lambda x: x[0])
threshold_lo = None
threshold_hi = None
for i, (ratio, m_star, label, _) in enumerate(sorted_summary):
    if m_star is None:
        threshold_lo = ratio
    elif threshold_lo is not None and threshold_hi is None:
        threshold_hi = ratio
        break

print()
if threshold_lo is not None and threshold_hi is not None:
    print(f"** Threshold located in ({threshold_lo:.4f}, {threshold_hi:.4f})")
elif threshold_lo is not None:
    print(f"** All curves above {threshold_lo:.4f} succeed; lower bound not improved")
else:
    print("** All curves succeeded (threshold < lowest tested ratio)")

print("\nDone.")
