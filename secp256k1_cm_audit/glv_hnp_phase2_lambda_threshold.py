"""
GLV-HNP Phase 2 — λ/n threshold study.

Goal: empirically determine the λ/n ratio below which LLL fails to recover
the secret key d in the GLV-aware HNP attack.

From 2026-07-26 run:
  - λ/n = 0.53 (n=199):  LLL succeeds 5/5 at m=6
  - λ/n = 0.34 (n~524k): LLL succeeds 3/3 at m=9
  - λ/n = 0.07 (n=2647): LLL AND BKZ(40) both fail

This script:
  1. Enumerates all j=0 prime-order curves with p in [2^7, 2^13].
  2. Bins them by λ/n ratio into 9 buckets: [0.05,0.10), ..., [0.45,0.50].
  3. For each non-empty bucket, picks the best representative
     (largest n in range, λ/n closest to bucket center).
  4. Tests LLL recovery at m=6, 5 seeds per bucket.
  5. Also tests BKZ(β=20) on failed buckets.
  6. Reports success rate vs λ/n.

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
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_point_on_curve(p, b, seed=42):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n):
    for seed in range(100):
        P = find_point_on_curve(p, b, seed)
        if P is None: continue
        if ec_mul(P, n, p) is None:
            return P
    return None

# ---------------------------------------------------------------------------
# CM / Eisenstein theory for j=0 curves
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

def find_twist_b(p, n):
    """Find b in {1,...,6} such that #E(F_p) = n for y^2 = x^3 + b."""
    for b in range(1, 7):
        P = find_point_on_curve(p, b)
        if P is None: continue
        if ec_mul(P, n, p) is None:
            return b
    # Extended search if needed
    for b in range(7, 50):
        P = find_point_on_curve(p, b)
        if P is None: continue
        if ec_mul(P, n, p) is None:
            return b
    return None

# ---------------------------------------------------------------------------
# Enumerate j=0 curves with prime n across a range of p-values
# Returns list of (p, b, n, lam, lam_over_n)
# ---------------------------------------------------------------------------

def enumerate_j0_curves(p_min, p_max):
    curves = []
    p = sympy.nextprime(p_min - 1)
    while p <= p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 5:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam, _ = glv_eigenvalue(n_cand)
                        if lam is None:
                            continue
                        ratio = lam / n_cand
                        # Find twist b
                        b = find_twist_b(p, n_cand)
                        if b is None:
                            continue
                        G = find_generator(p, b, n_cand)
                        if G is None:
                            continue
                        curves.append((p, b, n_cand, lam, ratio, G))
        p = sympy.nextprime(p)
    return curves

# ---------------------------------------------------------------------------
# GLV-HNP Phase 2 lattice attack — parameterized by curve
# ---------------------------------------------------------------------------

def gen_signatures(G, p, n, lam, d_secret, k1_bound, k2_bound, m, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 100000:
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

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = n // k2_bound
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]
    # mod-n rows
    for i in range(m):
        M[i][i] = n * S_K1
    # d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    # k2 rows
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    # Kannan row
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_K1, S_D, S_K2, S_KANNAN

def recover_d(reduced_rows, m, n, S_KANNAN, S_D, sigs, k1_bound, k2_bound, lam):
    dim = 2 * m + 2
    for row in reduced_rows:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        # Verify: A + B*d decompose as k1 + lam*k2 with k1 < k1_bound, k2 < k2_bound
        ok = True
        for s in sigs:
            kf = (s['A'] + s['B'] * d_cand) % n
            found = False
            for k2 in range(k2_bound):
                k1 = (kf - lam * k2) % n
                if 0 <= k1 < k1_bound:
                    found = True
                    break
            if not found:
                ok = False
                break
        if ok:
            return d_cand
    return None

def run_attack(G, p, n, lam, k1_bound, k2_bound, m, d_secret, seed, use_bkz=False, bkz_beta=20):
    sigs = gen_signatures(G, p, n, lam, d_secret, k1_bound, k2_bound, m, seed)
    if len(sigs) < m:
        return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    dim = len(M)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d_rec = recover_d(reduced, m, n, S_KANNAN, S_D, sigs, k1_bound, k2_bound, lam)
    return (d_rec is not None and d_rec == d_secret)

# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold study")
print("Enumerating j=0 curves with p in [2^7, 2^13] and binning by λ/n")
print("=" * 70)
print()

# Enumerate curves
print("Step 1: Enumerating j=0 curves...")
all_curves = enumerate_j0_curves(128, 8192)
print(f"  Found {len(all_curves)} curves with prime n and λ/n ∈ (0, 0.5]")

# Compute λ/n range
ratios = [c[4] for c in all_curves]
if ratios:
    print(f"  λ/n range: [{min(ratios):.4f}, {max(ratios):.4f}]")
print()

# Bin into 9 buckets: [0.05,0.10), ..., [0.45,0.50]
BUCKETS = [(0.05 * i, 0.05 * (i + 1)) for i in range(1, 10)]  # (0.05,0.10) up to (0.45,0.50)
BUCKETS = [(0.00, 0.05)] + BUCKETS  # also [0.00, 0.05)

bins = {b: [] for b in BUCKETS}
for curve in all_curves:
    r = curve[4]
    for b in BUCKETS:
        if b[0] <= r < b[1]:
            bins[b].append(curve)
            break

print("Step 2: Bin contents:")
for b, curves in sorted(bins.items()):
    if curves:
        best = max(curves, key=lambda c: c[2])  # largest n
        print(f"  λ/n ∈ [{b[0]:.2f},{b[1]:.2f}): {len(curves)} curves, "
              f"representative: p={best[0]}, n={best[2]}, λ={best[3]}, λ/n={best[4]:.4f}")
    else:
        print(f"  λ/n ∈ [{b[0]:.2f},{b[1]:.2f}): EMPTY")
print()

# For each non-empty bin, pick one representative and test LLL
M_TEST = 6      # signatures per test (from prior results: 8-bit curves need m=6)
N_SEEDS = 5     # seeds to test

print("Step 3: LLL recovery test (m=6, 5 seeds) per λ/n bucket:")
print(f"{'λ/n range':>20}  {'n':>8}  {'λ/n':>7}  {'LLL':>8}  {'BKZ-20':>8}  note")
print("-" * 75)

threshold_results = []

for b in sorted(bins.keys()):
    curves = bins[b]
    if not curves:
        print(f"  [{b[0]:.2f},{b[1]:.2f}): EMPTY — skip")
        continue

    # Pick representative: largest n in bucket
    rep = max(curves, key=lambda c: c[2])
    p, b_coeff, n, lam, ratio, G = rep

    # Bounds for this curve.
    # K1_BOUND=2 is the Phase 2 threat model: k1 is extremely small (k1 in {0,1}).
    # K2_BOUND = sqrt(n) + 1 is the GLV domain constraint.
    k1_bound = 2
    k2_bound = math.isqrt(n) + 1

    # Effective bias eff = k1_bound * k2_bound / n
    eff = k1_bound * k2_bound / n

    # Generate test secrets
    rng = random.Random(0xDEADBEEF)
    secrets = [rng.randint(1, n - 1) for _ in range(N_SEEDS)]
    seeds_list = [42, 1234, 9999, 7, 31415]

    lll_wins = 0
    bkz_wins = 0
    lll_needed = False

    for d_sec, seed in zip(secrets, seeds_list):
        ok = run_attack(G, p, n, lam, k1_bound, k2_bound, M_TEST, d_sec, seed, use_bkz=False)
        if ok:
            lll_wins += 1

    # If LLL failed any, try BKZ-20
    if lll_wins < N_SEEDS:
        lll_needed = True
        for d_sec, seed in zip(secrets, seeds_list):
            ok = run_attack(G, p, n, lam, k1_bound, k2_bound, M_TEST, d_sec, seed, use_bkz=True, bkz_beta=20)
            if ok:
                bkz_wins += 1

    bkz_str = f"{bkz_wins}/{N_SEEDS}" if lll_needed else "n/a"
    note = ""
    if lll_wins == N_SEEDS:
        note = "LLL OK"
    elif bkz_wins == N_SEEDS:
        note = "BKZ-20 OK"
    elif bkz_wins > 0:
        note = "BKZ-20 partial"
    else:
        note = "BOTH FAIL"

    print(f"  [{b[0]:.2f},{b[1]:.2f}): {n:>8}  {ratio:>7.4f}  {lll_wins}/{N_SEEDS}  "
          f"{bkz_str:>8}  {note}  "
          f"(k1_bd={k1_bound}, k2_bd={k2_bound}, eff={eff:.4f})")

    threshold_results.append({
        'bucket': b, 'ratio': ratio, 'n': n, 'lam': lam,
        'lll_wins': lll_wins, 'bkz_wins': bkz_wins,
        'k1_bound': k1_bound, 'k2_bound': k2_bound, 'eff': eff
    })

print()
print("=" * 70)
print("Step 4: Summary — LLL threshold estimate")
print()

# Find the threshold: largest ratio where LLL fully fails
lll_full_pass = [r for r in threshold_results if r['lll_wins'] == N_SEEDS]
lll_full_fail = [r for r in threshold_results if r['lll_wins'] == 0]
lll_partial = [r for r in threshold_results if 0 < r['lll_wins'] < N_SEEDS]

if lll_full_pass:
    min_lam_pass = min(r['ratio'] for r in lll_full_pass)
    print(f"  LLL fully succeeds (5/5): λ/n >= {min_lam_pass:.4f} (min observed pass ratio)")

if lll_full_fail:
    max_lam_fail = max(r['ratio'] for r in lll_full_fail)
    print(f"  LLL fully fails  (0/5):  λ/n <= {max_lam_fail:.4f} (max observed fail ratio)")

if lll_partial:
    partial_wins = [r['lll_wins'] for r in lll_partial]
    partial_ratios = [f"{r['ratio']:.4f}" for r in lll_partial]
    print(f"  LLL partial ({partial_wins}): ratios {partial_ratios}")

print()
print("Step 5: BKZ-20 rescue capability")
for r in threshold_results:
    if r['lll_wins'] < N_SEEDS and r['bkz_wins'] > r['lll_wins']:
        print(f"  λ/n={r['ratio']:.4f}: BKZ-20 improves {r['lll_wins']}/{N_SEEDS} "
              f"→ {r['bkz_wins']}/{N_SEEDS}")
    elif r['lll_wins'] < N_SEEDS and r['bkz_wins'] == r['lll_wins']:
        print(f"  λ/n={r['ratio']:.4f}: BKZ-20 gives NO improvement "
              f"({r['bkz_wins']}/{N_SEEDS})")

print()
print("Step 6: Relevance to secp256k1")
print("  secp256k1: λ/n ≈ 0.44 (from prior log entries)")
print("  => This ratio is in bucket [0.40,0.45) or [0.45,0.50)")
for r in threshold_results:
    if 0.40 <= r['ratio'] < 0.50:
        status = "LLL OK" if r['lll_wins'] == N_SEEDS else "LLL FAILS"
        print(f"  Bucket λ/n=[{r['bucket'][0]:.2f},{r['bucket'][1]:.2f}): "
              f"{status} ({r['lll_wins']}/{N_SEEDS})")

print()
print("Done.")
