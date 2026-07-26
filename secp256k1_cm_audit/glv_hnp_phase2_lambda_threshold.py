"""
GLV-HNP Phase 2: lambda/n threshold bisection study.

Identifies the minimum lam/n at which the GLV-aware LLL attack succeeds.
Uses the SMALLER GLV root consistently: lam = min(root1, root2) of
  x^2 + x + 1 = 0 mod n,  so lam/n < 0.5 for all curves.

Prior data points (from 2026-07-26 run):
  lam/n = 0.07 (n=2647, larger root 185):   LLL FAILS (all m, BKZ-40 also)
  lam/n = 0.34 (n=523969, smaller root):    LLL 3/3 at m=9
  lam/n = 0.53 (n=199, larger root 106):    LLL 3/3 at m=4  [used larger root]

This study fixes the SMALLER root and sweeps lam/n in eight bins:
  [0.05,0.09], [0.09,0.13], [0.13,0.17], [0.17,0.22],
  [0.22,0.27], [0.27,0.32], [0.32,0.37], [0.37,0.43]

For each bin: find one representative curve (p < 2^16, n prime, n≡1 mod 3),
then run LLL at m=5..13 with 3 seeds.  Report first m achieving 3/3.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ─────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────

def modinv(a, m):
    return pow(a, -1, m)

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    mm, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

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

def glv_eigenvalue_min(n):
    """Smaller root of x^2+x+1=0 mod n.  Returns lam in [1, (n-2)//2]."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    lam = min(r1, r2)
    if (lam * lam + lam + 1) % n != 0:
        return None
    return lam

# ─────────────────────────────────────────────────────────────
# Curve catalogue: enumerate all small j=0 prime-order curves
# ─────────────────────────────────────────────────────────────

def build_curve_catalogue(p_max=2**16):
    """
    Return list of (p, n, lam, ratio) for all prime p < p_max
    with p≡1(mod 3), prime group order n≡1(mod 3), valid lam.
    """
    cat = []
    p = 521
    while p < p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_c = p + 1 - t
                    if n_c < 8 or not sympy.isprime(n_c) or n_c % 3 != 1:
                        continue
                    lam = glv_eigenvalue_min(n_c)
                    if lam is None: continue
                    cat.append((p, n_c, lam, lam / n_c))
        p = sympy.nextprime(p)
    return cat

def pick_representative(cat, lo, hi):
    """
    From catalogue, pick a curve with ratio in [lo, hi] and n as large
    as possible (for a more realistic test).  Return (p, n, lam) or None.
    """
    candidates = [(p, n, lam, r) for (p, n, lam, r) in cat if lo <= r <= hi]
    if not candidates: return None
    # prefer larger n for a more interesting test
    candidates.sort(key=lambda x: x[1], reverse=True)
    return candidates[0][:3]  # (p, n, lam)

# ─────────────────────────────────────────────────────────────
# GLV lattice construction (column-diagonal scaling)
# ─────────────────────────────────────────────────────────────

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
    return M, S_K1, S_D, S_K2, S_KANNAN

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

def run_lll(curve_triple, b_coeff, m, d_secret, k1_bound, seed=42):
    p, n, lam = curve_triple
    G = find_generator(p, b_coeff, n)
    if G is None: return None  # no generator found
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b_coeff, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M, _, _, _, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

# ─────────────────────────────────────────────────────────────
# Main: sweep
# ─────────────────────────────────────────────────────────────

BINS = [
    (0.05, 0.09),
    (0.09, 0.13),
    (0.13, 0.17),
    (0.17, 0.22),
    (0.22, 0.27),
    (0.27, 0.32),
    (0.32, 0.37),
    (0.37, 0.43),
]
SEEDS = [42, 1234, 9999]
M_RANGE = range(5, 14)
K1_FRAC = 0.06  # K1_BOUND = ceil(K1_FRAC * sqrt(n)); eff ≈ K1_FRAC^2 ≈ 0.004

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection (smaller root)")
print("Curve search: p < 2^16, n prime, n≡1 mod 3")
print("=" * 70)

print("\nBuilding curve catalogue (p < 2^16)...", flush=True)
cat = build_curve_catalogue(p_max=2**16)
print(f"  Found {len(cat)} prime-order j=0 GLV curves.")

# Show distribution of lam/n
buckets = {}
for (p, n, lam, r) in cat:
    key = int(r * 10) / 10
    buckets[key] = buckets.get(key, 0) + 1
print("  Distribution of lam/n (smaller root):")
for k in sorted(buckets):
    print(f"    [{k:.1f}, {k+0.1:.1f}): {buckets[k]} curves")

# Pick representatives and cache generators
print("\nPicking representatives and finding generators...")
reps = {}  # bin_lo -> (p, n, lam, b_coeff, ratio)
for (lo, hi) in BINS:
    r3 = pick_representative(cat, lo, hi)
    if r3 is None:
        print(f"  [{lo:.2f},{hi:.2f}]: NOT FOUND")
        reps[(lo, hi)] = None
        continue
    p, n, lam = r3
    ratio = lam / n
    # Find b_coeff
    b_found = None
    for b_try in range(1, min(p, 100)):
        G = find_generator(p, b_try, n)
        if G is not None:
            b_found = b_try
            break
    if b_found is None:
        print(f"  [{lo:.2f},{hi:.2f}]: no generator for p={p}, n={n}")
        reps[(lo, hi)] = None
        continue
    reps[(lo, hi)] = (p, n, lam, b_found, ratio)
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(4, math.ceil(K1_FRAC * math.sqrt(n)))
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    print(f"  [{lo:.2f},{hi:.2f}]: p={p}, n={n} ({n.bit_length()}b), "
          f"lam={lam}, lam/n={ratio:.4f}, K1={k1_bound}, eff={eff:.5f}, "
          f"m_thresh≈{m_thresh:.0f}")

# Sweep
print("\n" + "=" * 70)
print("LLL sweep (m=5..13, 3 seeds each)")
print("=" * 70)

summary = []  # (ratio, first_success_m)

for (lo, hi) in BINS:
    rep = reps.get((lo, hi))
    if rep is None:
        summary.append(((lo+hi)/2, None, None))
        continue
    p, n, lam, b_coeff, ratio = rep
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(4, math.ceil(K1_FRAC * math.sqrt(n)))
    print(f"\n--- lam/n={ratio:.4f}  n={n} ({n.bit_length()}b)  lam={lam} ---")

    first_m = None
    for m in M_RANGE:
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed ^ 0xDEAD).randint(1, n - 1)
            # cache generator per-curve (rebuild each call is safe for small n)
            ok = run_lll((p, n, lam), b_coeff, m, d_trial, k1_bound, seed)
            if ok is None:
                print(f"  m={m}: ERROR (gen/sig failure)")
                break
            wins += ok
        marker = " ★" if wins == len(SEEDS) else ""
        print(f"  m={m}: {wins}/{len(SEEDS)}{marker}", flush=True)
        if wins == len(SEEDS) and first_m is None:
            first_m = m

    summary.append((ratio, first_m, n))

# Results table
print("\n" + "=" * 70)
print("RESULTS — λ/n threshold sweep (smaller root, LLL only)")
print(f"{'lam/n':>8}  {'n bits':>7}  {'first 3/3 at m':>15}  status")
print("-" * 55)

threshold_candidates = []
for (ratio, first_m, n_val) in summary:
    n_bits = n_val.bit_length() if n_val else "N/A"
    if n_val is None:
        status = "NOT FOUND"
        row = f"{ratio:>8.4f}  {'N/A':>7}  {'N/A':>15}  {status}"
    elif first_m is None:
        status = "NEVER 3/3 (m≤13)"
        row = f"{ratio:>8.4f}  {n_bits:>7}  {'never':>15}  {status}"
    else:
        status = f"SUCCESS"
        row = f"{ratio:>8.4f}  {n_bits:>7}  {first_m:>15}  {status}"
        threshold_candidates.append(ratio)
    print(row)

print("-" * 55)

if threshold_candidates:
    thr = min(threshold_candidates)
    print(f"\n  Minimum lam/n achieving LLL 3/3: ~{thr:.4f}")

    # Estimate threshold from the last failing bin and first succeeding bin
    fail_ratios = [r for (r, m, _) in summary if m is None and r is not None]
    if fail_ratios:
        last_fail = max(fail_ratios)
        print(f"  Last failing lam/n: ~{last_fail:.4f}")
        print(f"  Threshold estimate: lam/n in ({last_fail:.3f}, {thr:.3f})")
    else:
        print(f"  No failures detected — all tested lam/n succeeded.")
else:
    print("\n  All tested lam/n failed.  Threshold > max tested ratio.")

print("\nDone.")
