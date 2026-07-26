"""
GLV-HNP Phase 2: lambda/n threshold study.

Bisects the LLL success threshold between lambda/n = 0.07 (known fail)
and lambda/n = 0.34 (known pass) by sweeping j=0 CM curves bucketed by
lambda/n ratio.  Each bucket uses a ~10-12 bit prime-order curve with
K1_BOUND=2 (1-bit k1 bias) for comparability.

Prior data:
  lam/n = 0.07 (p=2677, n=2647): LLL & BKZ(40) both FAIL
  lam/n = 0.34 (p=524347, n=523969): LLL passes at m=9
  lam/n = 0.53 (p=211, n=199): LLL passes at m=4
  lam/n = 0.66 (p=2557, n=2659): LLL passes at m=7

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---- Minimal EC arithmetic ------------------------------------------------

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

def tonelli_shanks(n_val, p):
    n_val %= p
    if n_val == 0: return 0
    if pow(n_val, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n_val, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    mm, c, t, r = s, pow(z, q, p), pow(n_val, q, p), pow(n_val, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        bb = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, bb * bb % p, t * bb * bb % p, r * bb % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---- CM arithmetic -------------------------------------------------------

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
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    return min(r1, r2)

def find_b_for_order(p, n_target, max_b=300):
    """Find curve parameter b so that #E(F_p): y^2=x^3+b equals n_target."""
    rng = random.Random(42)
    for b_try in range(1, max_b):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_try) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            Q = ec_mul((x, y), n_target, p)
            if Q is None:
                return b_try
    return None

# ---- Lattice attack -------------------------------------------------------

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
        r_val = R[0] % n
        if r_val == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r_val) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r_val * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

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

def recover_d(A_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for i in range(dim):
        last = A_reduced[i, dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A_reduced[i, m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_experiment(curve_params, m, d_secret, k1_bound, seed=42,
                   use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    return recover_d(A, m, n, S_KANNAN, d_secret)

# ---- Bucket search -------------------------------------------------------

BUCKETS = [
    (0.05, 0.10),
    (0.10, 0.15),
    (0.15, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.35),
    (0.35, 0.40),
    (0.40, 0.45),
    (0.45, 0.50),
]
N_LO, N_HI = 200, 8000   # target n range for consistent m_thresh

def collect_bucket_curves(p_max=2**16):
    """Find one j=0 curve per lambda/n bucket with n in [N_LO, N_HI]."""
    bucket_curves = [None] * len(BUCKETS)
    p = 5
    total_checked = 0
    while p <= p_max and not all(v is not None for v in bucket_curves):
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if not (N_LO <= n_cand <= N_HI): continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    for i, (lo, hi) in enumerate(BUCKETS):
                        if lo <= ratio < hi and bucket_curves[i] is None:
                            b_param = find_b_for_order(p, n_cand)
                            if b_param is None: continue
                            G = find_generator(p, b_param, n_cand)
                            if G is None: continue
                            bucket_curves[i] = (p, b_param, n_cand, lam, G)
                            total_checked += 1
        p = sympy.nextprime(p)
    return bucket_curves

# ---- Main ----------------------------------------------------------------

K1_BOUND = 2
MAX_M = 16
SEEDS = [42, 1234, 9999]

print("=" * 70)
print("GLV-HNP Phase 2 — lambda/n threshold study")
print("K1_BOUND=2 (1-bit k1 bias), n in [200,8000], seeds=3")
print("=" * 70)

print("\nSearching for representative curves per lambda/n bucket...")
bucket_curves = collect_bucket_curves(p_max=2**16)

# Print found curves
print(f"\n{'Bucket':18s}  {'p':8s}  {'n':8s}  {'lam':8s}  {'lam/n':6s}  {'n bits':6s}  {'found':5s}")
print("-" * 72)
for i, (lo, hi) in enumerate(BUCKETS):
    if bucket_curves[i] is None:
        print(f"[{lo:.2f},{hi:.2f})       NONE FOUND")
    else:
        p, b, n, lam, G = bucket_curves[i]
        print(f"[{lo:.2f},{hi:.2f})  {p:8d}  {n:8d}  {lam:8d}  {lam/n:.4f}  {n.bit_length():6d}  ✓")

# ---- Sweep: find minimum m for 3/3 LLL success ---------------------------

print(f"\n{'='*70}")
print("LLL sweep (3/3 success across seeds {42, 1234, 9999})")
print(f"{'='*70}")

lll_results = {}   # bucket_idx -> min_m or None
bkz_results = {}   # bucket_idx -> min_m or None (for FAIL buckets)

for i, (lo, hi) in enumerate(BUCKETS):
    cp = bucket_curves[i]
    if cp is None:
        lll_results[i] = None
        print(f"\nBucket [{lo:.2f},{hi:.2f}): NO CURVE — skipping")
        continue
    p, b, n, lam, G = cp
    k2_bound = math.isqrt(n) + 1
    eff = K1_BOUND * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    print(f"\nBucket [{lo:.2f},{hi:.2f}): p={p}, n={n}, lam={lam}, lam/n={lam/n:.4f}, "
          f"eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")

    min_m_lll = None
    for m in range(3, MAX_M + 1):
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok = run_experiment(cp, m, d_trial, K1_BOUND, seed, use_bkz=False)
            wins += ok
        print(f"  m={m:2d}: LLL {wins}/{len(SEEDS)}", end="", flush=True)
        if wins == len(SEEDS):
            min_m_lll = m
            print(f"  ← 3/3 ✓  STOP", flush=True)
            break
        print(flush=True)
    lll_results[i] = min_m_lll

# BKZ(40) rescue for the FAIL buckets
fail_indices = [i for i, v in lll_results.items() if v is None and bucket_curves[i] is not None]
if fail_indices:
    print(f"\n{'='*70}")
    print(f"BKZ(beta=40) rescue on FAIL buckets: {fail_indices}")
    print(f"{'='*70}")
    for i in fail_indices:
        cp = bucket_curves[i]
        lo, hi = BUCKETS[i]
        p, b, n, lam, G = cp
        print(f"\nBucket [{lo:.2f},{hi:.2f}): p={p}, n={n}, lam/n={lam/n:.4f}")
        min_m_bkz = None
        for m in range(3, MAX_M + 1):
            wins = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                ok = run_experiment(cp, m, d_trial, K1_BOUND, seed,
                                    use_bkz=True, bkz_beta=40)
                wins += ok
            print(f"  m={m:2d}: BKZ(40) {wins}/{len(SEEDS)}", end="", flush=True)
            if wins == len(SEEDS):
                min_m_bkz = m
                print(f"  ← 3/3 ✓  STOP", flush=True)
                break
            print(flush=True)
        bkz_results[i] = min_m_bkz

# ---- Summary table -------------------------------------------------------

print(f"\n{'='*70}")
print("SUMMARY — lambda/n threshold for LLL success")
print(f"{'='*70}")
print(f"{'Bucket':18s}  {'lam/n':6s}  {'n':6s}  {'LLL min-m':10s}  {'BKZ(40) min-m':14s}")
print("-" * 68)
for i, (lo, hi) in enumerate(BUCKETS):
    cp = bucket_curves[i]
    if cp is None:
        print(f"[{lo:.2f},{hi:.2f})      —       —  NO CURVE      —")
        continue
    p, b, n, lam, G = cp
    lll_m = lll_results.get(i)
    bkz_m = bkz_results.get(i, "—")
    lll_str = f"m={lll_m}" if lll_m is not None else "FAIL"
    bkz_str = f"m={bkz_m}" if isinstance(bkz_m, int) else ("FAIL" if bkz_m is None else "—")
    print(f"[{lo:.2f},{hi:.2f})  {lam/n:.4f}  {n:6d}  {lll_str:10s}  {bkz_str:14s}")

# Identify threshold
pass_buckets = [i for i, v in lll_results.items() if v is not None]
fail_buckets = [i for i, v in lll_results.items() if v is None and bucket_curves[i] is not None]
if pass_buckets and fail_buckets:
    min_pass_ratio = min(bucket_curves[i][3] / bucket_curves[i][2] for i in pass_buckets)
    max_fail_ratio = max(bucket_curves[i][3] / bucket_curves[i][2] for i in fail_buckets)
    print(f"\nConclusion:")
    print(f"  LLL PASS: lam/n >= {min_pass_ratio:.4f} (bucket {BUCKETS[min(pass_buckets)]})")
    print(f"  LLL FAIL: lam/n <= {max_fail_ratio:.4f} (bucket {BUCKETS[max(fail_buckets)]})")
    print(f"  Threshold: {max_fail_ratio:.4f} < lam/n* <= {min_pass_ratio:.4f}")
elif not fail_buckets:
    print("\nConclusion: LLL succeeds for all tested lambda/n values.")
else:
    print("\nConclusion: LLL fails for all tested lambda/n values.")

print("\nDone.")
