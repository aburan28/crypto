"""
GLV-HNP Phase 2 — Thread 20: λ/n threshold bisection.

Known from prior runs (2026-07-26):
  lam/n = 0.07  (p=2677, n=2647):   LLL *and* BKZ(40) always FAIL
  lam/n = 0.34  (p=524347, n=523969): LLL succeeds at m=9

Goal: find the threshold lam/n* ∈ (0.07, 0.34) where LLL transitions from
fail to success.  Secondary: does BKZ(beta=40) push the threshold lower?

Method:
  - Search 12-bit j=0 curves across 8 λ/n buckets
  - For each representative curve run LLL at m ∈ {5,7,9,11,13}, 3 seeds
  - Also run BKZ(40) for curves with lam/n < 0.22 (potential rescue zone)
  - Report first-success m as a function of lam/n

Run: python3 glv_hnp_phase2_lambda_threshold.py
Date: 2026-07-28 (autolab Thread 20)
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic  (generic p, b)
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

def tonelli(n, p):
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

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM / GLV utilities  (same as glv_hnp_phase2_20bit.py)
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
    sq = tonelli(neg3, n)
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

# ---------------------------------------------------------------------------
# Curve finder: search 12-bit primes, bucket by lam/n
# ---------------------------------------------------------------------------

BUCKETS = [
    ('0.07-0.10', 0.07, 0.10),
    ('0.10-0.14', 0.10, 0.14),
    ('0.14-0.18', 0.14, 0.18),
    ('0.18-0.23', 0.18, 0.23),
    ('0.23-0.29', 0.23, 0.29),
    ('0.29-0.36', 0.29, 0.36),
    ('0.36-0.44', 0.36, 0.44),
    ('0.44-0.50', 0.44, 0.50),
]

def find_curves_by_lambda_ratio(n_bits=12, verbose=True):
    """
    Sweep n_bits-prime primes p (j=0 candidates), find one representative
    curve per BUCKET.  Returns dict: bucket_name -> (p, b, n, lam, G).
    """
    found = {}
    lo, hi = 2**(n_bits - 1), 2**n_bits
    p = sympy.nextprime(lo)
    searched = 0

    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 4 or n_cand >= hi:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam, _ = glv_eigenvalue(n_cand)
                        if lam is None:
                            continue
                        ratio = lam / n_cand
                        for (bname, lo_r, hi_r) in BUCKETS:
                            if lo_r <= ratio < hi_r and bname not in found:
                                # Find b such that #E(F_p)=n_cand
                                found_b = None
                                for b_try in range(1, 200):
                                    P = None
                                    rng_s = random.Random(999 + b_try)
                                    for _ in range(100):
                                        x = rng_s.randint(0, p - 1)
                                        rhs = (pow(x, 3, p) + b_try) % p
                                        y = tonelli(rhs, p)
                                        if y is not None and y != 0:
                                            P = (x, y)
                                            break
                                    if P is None:
                                        continue
                                    Q = ec_mul(P, n_cand, p)
                                    if Q is None:
                                        found_b = b_try
                                        break
                                if found_b is None:
                                    continue
                                G = find_generator(p, found_b, n_cand)
                                if G is None:
                                    continue
                                found[bname] = (p, found_b, n_cand, lam, G)
                                if verbose:
                                    print(f"  [{bname}] p={p}, n={n_cand}, "
                                          f"lam={lam}, lam/n={ratio:.4f}, b={found_b}")
                                break
        searched += 1
        if len(found) == len(BUCKETS):
            break
        p = sympy.nextprime(p)

    if verbose:
        print(f"  Searched {searched} primes, found {len(found)}/{len(BUCKETS)} buckets.")
    return found

# ---------------------------------------------------------------------------
# Signature generation  (planted k = k1 + lam*k2)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    for _ in range(200000):
        if len(sigs) == m:
            break
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
# GLV lattice builder  (same scaling as phase2_20bit.py)
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = max(1, n // k1_bound)
    S_D  = 1
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

def recover_d(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Single experiment
# ---------------------------------------------------------------------------

def run_one(curve, m, seed, use_bkz=False, bkz_beta=40):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.06 * math.sqrt(n)))  # eff ≈ 0.06

    d_secret = random.Random(seed + 7777).randint(1, n - 1)
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
    M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(M_red, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Sweep a single curve over m values and seeds
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999]
M_VALS = [5, 7, 9, 11, 13]

def sweep_one_curve(label, curve, use_bkz=False, bkz_beta=40):
    p, b, n, lam, G = curve
    algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.06 * math.sqrt(n)))
    eff = k1_bound * k2_bound / n
    print(f"  [{algo}] {label}  lam/n={lam/n:.4f}  n={n}  K1={k1_bound}  eff={eff:.4f}")
    first_success = None
    for m in M_VALS:
        wins = sum(run_one(curve, m, s, use_bkz=use_bkz, bkz_beta=bkz_beta)
                   for s in SEEDS)
        status = f"{wins}/{len(SEEDS)}"
        marker = " ← FIRST SUCCESS" if (wins == len(SEEDS) and first_success is None) else ""
        print(f"    m={m:2d}: {status}{marker}")
        if wins == len(SEEDS) and first_success is None:
            first_success = m
    return first_success  # None = never succeeded

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — Thread 20: λ/n threshold bisection")
print("=" * 70)
print()

# ---- Step 1: find one curve per bucket ------------------------------------
print("Step 1: Searching for representative curves per λ/n bucket (12-bit)...")
curves_by_bucket = find_curves_by_lambda_ratio(n_bits=12, verbose=True)
print()

# ---- Step 2: inject known anchor curves ----------------------------------
# Anchors from 2026-07-26 logs:
#   lam/n=0.07 (p=2677, n=2647, lam=185) — KNOWN FAIL
#   lam/n=0.53 (p=211,  n=199,  lam=106) — KNOWN SUCCESS (but lam > n/2 → use min)
print("Step 2: Injecting known anchor curves from prior sessions...")

G_2677 = find_generator(2677, 2, 2647)
if G_2677 is not None:
    curves_by_bucket['anchor_0.07'] = (2677, 2, 2647, 185, G_2677)
    print("  anchor_0.07: p=2677, n=2647, lam=185, lam/n=0.0699 (known fail)")

G_211 = find_generator(211, 2, 199)
if G_211 is not None:
    lam_199, _ = glv_eigenvalue(199)
    if lam_199 is None:
        lam_199 = 92   # min(106, 92)=92; 92/199=0.462
    curves_by_bucket['anchor_0.46'] = (211, 2, 199, lam_199, G_211)
    print(f"  anchor_0.46: p=211, n=199, lam={lam_199}, lam/n={lam_199/199:.4f} (known success)")
print()

# ---- Step 3: LLL sweep over all buckets ----------------------------------
print("Step 3: LLL sweep (m ∈ {5,7,9,11,13}, 3 seeds each)")
print("=" * 70)

results_lll = {}
sorted_buckets = sorted(
    curves_by_bucket.items(),
    key=lambda kv: kv[1][2] and kv[1][3] / kv[1][2]   # sort by lam/n
)
for bname, curve in sorted_buckets:
    p, b, n, lam, G = curve
    label = f"{bname}  n={n}  lam={lam}"
    first_m = sweep_one_curve(label, curve, use_bkz=False)
    results_lll[bname] = (lam / n, first_m)
    print()

# ---- Step 4: BKZ(40) on potential rescue zone (lam/n < 0.22) -----------
print("=" * 70)
print("Step 4: BKZ(40) rescue attempt on lam/n < 0.22 curves")
print("=" * 70)

results_bkz = {}
for bname, curve in sorted_buckets:
    p, b, n, lam, G = curve
    ratio = lam / n
    if ratio < 0.22:
        label = f"{bname}  n={n}  lam={lam}"
        first_m = sweep_one_curve(label, curve, use_bkz=True, bkz_beta=40)
        results_bkz[bname] = (ratio, first_m)
        print()

# ---- Step 5: Summary table -----------------------------------------------
print("=" * 70)
print("SUMMARY — λ/n threshold bisection")
print("=" * 70)
print(f"{'Bucket':<20} {'lam/n':>6} {'LLL first-m':>12} {'BKZ(40)':>10}")
print("-" * 54)
for bname, curve in sorted_buckets:
    p, b, n, lam, G = curve
    ratio = lam / n
    lll_m  = results_lll.get(bname, (ratio, None))[1]
    bkz_m  = results_bkz.get(bname, (ratio, None))[1]
    lll_s  = f"m={lll_m}" if lll_m else "FAIL"
    bkz_s  = f"m={bkz_m}" if bkz_m else ("FAIL" if bname in results_bkz else "not tested")
    print(f"  {bname:<18} {ratio:>6.4f} {lll_s:>12} {bkz_s:>10}")

print()
# Identify threshold
lll_success = [(ratio, m) for _, (ratio, m) in results_lll.items() if m is not None]
lll_fail    = [(ratio, m) for _, (ratio, m) in results_lll.items() if m is None]
if lll_success and lll_fail:
    thr_lo = max(r for r, _ in lll_fail)
    thr_hi = min(r for r, _ in lll_success)
    print(f"LLL threshold: lam/n ∈ ({thr_lo:.4f}, {thr_hi:.4f})")
    print(f"Estimate: lam/n* ≈ {(thr_lo + thr_hi) / 2:.4f}")
elif not lll_fail:
    print("LLL succeeded for ALL buckets tested.")
elif not lll_success:
    print("LLL FAILED for all buckets — threshold > 0.50 (unexpected).")

bkz_success = [(r, m) for _, (r, m) in results_bkz.items() if m is not None]
bkz_fail    = [(r, m) for _, (r, m) in results_bkz.items() if m is None]
if bkz_success:
    bkz_thr_hi = min(r for r, _ in bkz_success)
    print(f"BKZ(40) rescues down to lam/n ≥ {bkz_thr_hi:.4f}")
elif bkz_fail:
    print("BKZ(40) did NOT rescue any additional curve in the lam/n < 0.22 zone.")

print("\nDone.")
