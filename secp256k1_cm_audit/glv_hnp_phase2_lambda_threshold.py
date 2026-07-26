"""
GLV-HNP Phase 2: λ/n threshold study.

Goal: bisect the critical λ/n ratio below which the column-scaled GLV lattice
attack fails, given the known bracket:
  - lam/n = 0.07 (p=2677, n=2647): LLL and BKZ(40) FAIL  (2026-07-26)
  - lam/n = 0.34 (20-bit curve):   LLL SUCCEEDS at m=9    (2026-07-26)

Method:
  1. Enumerate j=0 prime-order curves for small primes p in [100, 6000].
  2. Compute lam/n for each, bucket into intervals of width 0.05 in [0.05, 0.50].
  3. For each bucket, select up to 3 representative curves.
  4. Run LLL at a fixed m (m_fixed = ceil(log(n)/log(1/eff)) + 3) on 3 seeds.
  5. Report success rate per bucket → locate the threshold.

Known reference points (appended as annotations):
  lam/n=0.53 (p=211, n=199):   LLL 3/3 at m=4
  lam/n=0.66 (p=2557, n=2659): LLL 3/3 at m=7
  lam/n=0.34 (20-bit):         LLL 3/3 at m=9
  lam/n=0.07 (p=2677, n=2647): LLL 0/3, BKZ(40) 0/3

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
# ---------------------------------------------------------------------------

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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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
    """Find a generator of E: y^2=x^3+b/F_p of prime order n."""
    rng = random.Random(42)
    for _ in range(5000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM / Eisenstein
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
    """Return the smaller root of x^2+x+1=0 mod n (requires n≡1 mod 3)."""
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
    return min(r1, r2)

# ---------------------------------------------------------------------------
# Curve enumeration
# ---------------------------------------------------------------------------

def enumerate_curves(p_min, p_max):
    """
    Enumerate j=0 prime-order curves for p in [p_min, p_max].
    Returns list of dicts: {p, b, n, lam, lam_ratio}.
    """
    curves = []
    p = sympy.nextprime(p_min - 1)
    while p <= p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n = p + 1 - t
                    if n < 7 or n > 2 * p:
                        continue
                    if not sympy.isprime(n):
                        continue
                    if n % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n)
                    if lam is None:
                        continue
                    ratio = lam / n
                    # Find the twist b
                    found_b = None
                    for b_try in range(1, min(p, 300)):
                        rng = random.Random(b_try * 17 + 3)
                        x = rng.randint(0, p - 1)
                        rhs = (pow(x, 3, p) + b_try) % p
                        y = tonelli_shanks(rhs, p)
                        if y is not None and y != 0:
                            P = (x, y)
                            if ec_mul(P, n, p) is None:
                                found_b = b_try
                                break
                    if found_b is not None:
                        curves.append({
                            'p': p, 'b': found_b, 'n': n,
                            'lam': lam, 'ratio': ratio
                        })
        p = sympy.nextprime(p)
    return curves

# ---------------------------------------------------------------------------
# GLV lattice and attack (from glv_hnp_phase2_attack.py)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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
    return M, S_KANNAN

def recover_d(M_reduced, m, n, S_KANNAN):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        return d_cand
    return None

def run_lll_attack(curve, m, d_secret, seed=42):
    p, b, n, lam = curve['p'], curve['b'], curve['n'], curve['lam']
    k2_bound = math.isqrt(n) + 1
    # k1_bound: target eff ~ 0.05 (5% of n covered per sig)
    k1_bound = max(2, int(0.05 * math.sqrt(n)))

    G = find_generator(p, b, n)
    if G is None:
        return None  # can't test

    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None

    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d_rec = recover_d(reduced, m, n, S_KANNAN)
    return d_rec is not None and d_rec == d_secret

def test_curve(curve, m_fixed, seeds=(42, 1234, 9999)):
    wins = 0
    for seed in seeds:
        d = random.Random(seed + 5555).randint(1, curve['n'] - 1)
        result = run_lll_attack(curve, m_fixed, d, seed)
        if result:
            wins += 1
    return wins

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold study")
print("Known bracket: lam/n=0.07 FAILS, lam/n=0.34 SUCCEEDS")
print("=" * 70)
print()

# --- Step 1: enumerate curves ---------------------------------------------
print("Enumerating j=0 prime-order curves, p in [100, 6000]...")
all_curves = enumerate_curves(100, 6000)
print(f"  Found {len(all_curves)} eligible curves (j=0, prime n, n≡1 mod 3, has GLV λ)")
print()

# --- Step 2: bin by lam/n in intervals of 0.05 ---------------------------
BIN_WIDTH = 0.05
N_BINS = 9   # [0.05, 0.10), [0.10, 0.15), ..., [0.45, 0.50)
bins = {i: [] for i in range(N_BINS)}
for c in all_curves:
    r = c['ratio']
    if r < 0.05 or r >= 0.50:
        continue
    idx = int((r - 0.05) / BIN_WIDTH)
    if 0 <= idx < N_BINS:
        bins[idx].append(c)

print("Curve counts per λ/n bin:")
for i in range(N_BINS):
    lo = 0.05 + i * BIN_WIDTH
    hi = lo + BIN_WIDTH
    print(f"  [{lo:.2f}, {hi:.2f}): {len(bins[i])} curves")
print()

# --- Step 3: per-bin LLL test ---------------------------------------------
M_FIXED = 12  # Fixed number of signatures per test
SEEDS = (42, 1234, 9999)
REPS_PER_BIN = 3  # Test up to 3 representative curves per bin

print(f"LLL attack test: m={M_FIXED}, seeds={SEEDS}, up to {REPS_PER_BIN} curves per bin")
print("=" * 70)

results = {}  # bin_idx -> list of (ratio, wins, n_curves_tested)

for i in range(N_BINS):
    lo = 0.05 + i * BIN_WIDTH
    hi = lo + BIN_WIDTH
    bucket = bins[i]
    if not bucket:
        print(f"  [{lo:.2f}, {hi:.2f}): NO CURVES — skip")
        results[i] = []
        continue

    # Select diverse representatives (by ratio spread within the bin)
    bucket_sorted = sorted(bucket, key=lambda c: c['ratio'])
    step = max(1, len(bucket_sorted) // REPS_PER_BIN)
    reps = bucket_sorted[::step][:REPS_PER_BIN]

    bin_results = []
    for c in reps:
        wins = test_curve(c, M_FIXED, SEEDS)
        bin_results.append((c['ratio'], c['n'], wins))
        status = "✓" if wins == len(SEEDS) else ("△" if wins > 0 else "✗")
        print(f"  [{lo:.2f}, {hi:.2f}) lam/n={c['ratio']:.4f} n={c['n']:6d} "
              f"wins={wins}/{len(SEEDS)} {status}")

    results[i] = bin_results

print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'λ/n range':<18} {'ratio':>8} {'n':>8} {'success':>10} {'verdict':>10}")
print("-" * 60)

threshold_lo = None
threshold_hi = None

for i in range(N_BINS):
    lo = 0.05 + i * BIN_WIDTH
    hi = lo + BIN_WIDTH
    bin_res = results.get(i, [])
    if not bin_res:
        print(f"  [{lo:.2f}, {hi:.2f})  {'(no curves)':>38}")
        continue
    for (ratio, n_val, wins) in bin_res:
        verdict = "PASS" if wins == len(SEEDS) else ("PARTIAL" if wins > 0 else "FAIL")
        print(f"  [{lo:.2f}, {hi:.2f})  {ratio:>8.4f} {n_val:>8d} {wins:>3}/{len(SEEDS):<6} {verdict:>10}")
        if verdict == "FAIL" and threshold_lo is None:
            threshold_lo = ratio
        if verdict in ("PASS", "PARTIAL") and threshold_lo is not None:
            threshold_hi = ratio

print()
# Include known reference points
print("Reference points (from 2026-07-26 run):")
ref = [
    (0.07, "p=2677,n=2647", "FAIL (LLL and BKZ(40))"),
    (0.34, "20-bit CM",     "PASS (LLL 3/3 at m=9)"),
    (0.44, "secp256k1 λ",   "expected PASS (λ/n structure)"),
    (0.53, "p=211,n=199",   "PASS (LLL 3/3 at m=4)"),
    (0.66, "p=2557,n=2659", "PASS (LLL 3/3 at m=7)"),
]
for (r, label, status) in ref:
    print(f"  lam/n={r:.2f}  {label:<25}  {status}")

print()
if threshold_lo is not None and threshold_hi is not None:
    print(f"Estimated threshold bracket: lam/n in [{threshold_lo:.4f}, {threshold_hi:.4f})")
elif threshold_lo is not None:
    print(f"Lowest known FAIL: lam/n={threshold_lo:.4f} (no PASS found above it in this run)")
else:
    print("No FAIL observed in this run — threshold < lowest tested ratio or all PASS")

print()
print("Done.")
