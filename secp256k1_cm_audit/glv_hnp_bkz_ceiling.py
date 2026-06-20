"""
GLV-HNP Phase 2: BKZ ceiling for Curve A (2026-06-20).

Established facts (2026-06-19):
  - Curve A (p=524347, lam/n=0.34): LLL fails at K1≥55 (eff≥0.076)
  - BKZ-20 rescues K1=55 at m=14 (first 3/3)
  - Root cause: LLL reduction quality, not lattice structure
  - The planted short vector EXISTS (BKZ finds it)

Today's question:
  Does BKZ-20's ceiling match Curve C's LLL ceiling (eff≈0.13+)?
  If yes: confirms reduction-quality interpretation for ALL K1 values.
  If no (BKZ also fails at some K1): stronger reduction (BKZ-40) might help,
         or we have a genuine lattice obstruction at high K1.

Experiments:
  1. BKZ-20 on Curve A at K1=72, 90, 109 with m=10..24, 3 seeds.
     Compare to Curve C's LLL success (all K1 up to 109 at m≤13).
  2. If BKZ-20 fails at K1=72+: test BKZ-40 as escalation.
  3. Final comparison table: Curve A BKZ vs Curve C LLL.

Reference curves:
  Curve A: p=524347, n=523969, lam=177902  (lam/n=0.3395, lam≈n/3)
  Curve C: p=624517, n=622957, lam=178530  (lam/n=0.2867)
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic
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
    mv, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mv - i - 1), p)
        mv, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(200000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            pt = (x, y)
            if ec_mul(pt, n, p) is None:
                return pt
    return None

# ---------------------------------------------------------------------------
# GLV-HNP lattice (same structure as prior scripts)
# ---------------------------------------------------------------------------

def gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
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
        sigs.append({
            'A': h * modinv(s, n) % n,
            'B': r * modinv(s, n) % n,
            'k1': k1, 'k2': k2,
        })
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
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
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        if abs(row[dim-1]) != S_KANNAN: continue
        sign = 1 if row[dim-1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_reduction(M, dim, method, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if method == 'lll':
        LLL.reduction(A)
    else:
        par = BKZ.Param(block_size=beta, strategies=[], auto_abort=True)
        BKZ.reduction(A, par)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def experiment(curve, m, d_secret, k1_bound, seed, method='lll', beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    reduced = run_reduction(M, dim, method, beta)
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

# ---------------------------------------------------------------------------
# Setup curves
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: BKZ ceiling for Curve A (2026-06-20)")
print("=" * 70)

# Curve A: lam ≈ n/3 (the problematic curve)
P_A, B_A, N_A, LAM_A = 524347, 2, 523969, 177902
G_A = find_generator(P_A, B_A, N_A)
assert G_A is not None and ec_mul(G_A, N_A, P_A) is None
CURVE_A = (P_A, B_A, N_A, LAM_A, G_A)

# Curve C: reference (LLL succeeds at all K1 up to 109)
P_C, B_C, N_C, LAM_C = 624517, 15, 622957, 178530
G_C = find_generator(P_C, B_C, N_C)
assert G_C is not None and ec_mul(G_C, N_C, P_C) is None
CURVE_C = (P_C, B_C, N_C, LAM_C, G_C)

print(f"Curve A: p={P_A}, n={N_A}, lam={LAM_A} (lam/n={LAM_A/N_A:.4f})")
print(f"  3*lam mod n = {3*LAM_A % N_A}  δ/n = {(3*LAM_A % N_A)/N_A:.4f}")
print(f"Curve C: p={P_C}, n={N_C}, lam={LAM_C} (lam/n={LAM_C/N_C:.4f})")
print()

SEEDS_3 = [42, 1234, 9999]
K1_VALUES = [55, 72, 90, 109]

# Reference: Curve C LLL results (from 2026-06-18 log)
# K1=55→m=12, K1=72→m=13, K1=90→m=13, K1=109→m=13
curve_c_lll_ref = {55: 12, 72: 13, 90: 13, 109: 13}

# ---------------------------------------------------------------------------
# Experiment 1: BKZ-20 on Curve A at K1=55,72,90,109 (m=10..24)
# Key question: does BKZ-20's ceiling reach Curve C's LLL ceiling?
# ---------------------------------------------------------------------------

print("--- Experiment 1: Curve A + BKZ-20, K1∈{55,72,90,109}, m=10..24, 3 seeds ---")
print()

bkz20_results = {}  # K1 → first m at 3/3

for k1 in K1_VALUES:
    k2_bound = math.isqrt(N_A) + 1
    eff = k1 * k2_bound / N_A
    m_thresh = math.ceil(math.log(N_A) / (-math.log(k1 * k2_bound / N_A)))
    print(f"  K1={k1}, eff={eff:.4f}, m_thresh={m_thresh}")
    first_33 = None
    t0 = time.time()
    for m in range(10, 25):
        wins = sum(
            experiment(CURVE_A, m,
                       random.Random(s + 9001).randint(1, N_A - 1),
                       k1, s, method='bkz', beta=20)
            for s in SEEDS_3
        )
        marker = ""
        if wins == 3 and first_33 is None:
            first_33 = m
            marker = " ← first 3/3"
        print(f"    m={m:2d}: {wins}/3{marker}")
    elapsed = time.time() - t0
    bkz20_results[k1] = first_33
    ref = curve_c_lll_ref.get(k1, '?')
    status = f"m={first_33}" if first_33 else "FAIL (m≤24)"
    print(f"  BKZ-20 result: {status}  |  Curve C LLL ref: m={ref}  [{elapsed:.1f}s]")
    print()

# ---------------------------------------------------------------------------
# Experiment 2: BKZ-40 escalation for K1 values where BKZ-20 fails
# Only run if BKZ-20 failed for any K1
# ---------------------------------------------------------------------------

bkz20_failed = [k1 for k1 in K1_VALUES if bkz20_results[k1] is None]

if bkz20_failed:
    print(f"--- Experiment 2: BKZ-40 escalation for K1∈{bkz20_failed}, m=10..22, 3 seeds ---")
    print()
    bkz40_results = {}
    for k1 in bkz20_failed:
        k2_bound = math.isqrt(N_A) + 1
        eff = k1 * k2_bound / N_A
        print(f"  K1={k1}, eff={eff:.4f}")
        first_33 = None
        t0 = time.time()
        for m in range(10, 23):
            wins = sum(
                experiment(CURVE_A, m,
                           random.Random(s + 9001).randint(1, N_A - 1),
                           k1, s, method='bkz', beta=40)
                for s in SEEDS_3
            )
            marker = ""
            if wins == 3 and first_33 is None:
                first_33 = m
                marker = " ← first 3/3"
            print(f"    m={m:2d}: {wins}/3 (BKZ-40){marker}")
        elapsed = time.time() - t0
        bkz40_results[k1] = first_33
        status = f"m={first_33}" if first_33 else "FAIL (m≤22)"
        ref = curve_c_lll_ref.get(k1, '?')
        print(f"  BKZ-40 result: {status}  |  Curve C LLL ref: m={ref}  [{elapsed:.1f}s]")
        print()
else:
    bkz40_results = {}
    print("--- Experiment 2: BKZ-40 escalation skipped (BKZ-20 succeeded at all K1) ---")
    print()

# ---------------------------------------------------------------------------
# Experiment 3: Curve C LLL verification (confirm reference numbers)
# Run K1=72, 90, 109 on Curve C to verify the 2026-06-18 claims.
# ---------------------------------------------------------------------------

print("--- Experiment 3: Curve C LLL reference verification (K1=72,90,109, m=10..16) ---")
print()

curve_c_actual = {}
for k1 in [72, 90, 109]:
    k2_bound_c = math.isqrt(N_C) + 1
    eff_c = k1 * k2_bound_c / N_C
    m_thresh_c = math.ceil(math.log(N_C) / (-math.log(k1 * k2_bound_c / N_C)))
    print(f"  K1={k1}, eff={eff_c:.4f}, m_thresh={m_thresh_c}")
    first_33 = None
    for m in range(10, 17):
        wins = sum(
            experiment(CURVE_C, m,
                       random.Random(s + 8001).randint(1, N_C - 1),
                       k1, s, method='lll')
            for s in SEEDS_3
        )
        marker = ""
        if wins == 3 and first_33 is None:
            first_33 = m
            marker = " ← first 3/3"
        print(f"    m={m:2d}: {wins}/3{marker}")
    curve_c_actual[k1] = first_33
    print(f"  Curve C LLL: {f'm={first_33}' if first_33 else 'FAIL'}")
    print()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print("=" * 70)
print("SUMMARY: Curve A BKZ ceiling vs Curve C LLL ceiling")
print("=" * 70)
print()
print(f"{'K1':>4} {'eff':>6} | {'A LLL(prev)':>12} | {'A BKZ-20':>10} | {'A BKZ-40':>10} | {'C LLL':>8}")
print("-" * 65)

curve_a_lll_prev = {36: 13, 55: None, 72: None, 90: None, 109: None}  # from 2026-06-19
for k1 in K1_VALUES:
    k2_bound = math.isqrt(N_A) + 1
    eff = k1 * k2_bound / N_A
    a_lll = curve_a_lll_prev.get(k1)
    a_lll_str = f"m={a_lll}" if a_lll else "FAIL"
    a_b20 = bkz20_results.get(k1)
    a_b20_str = f"m={a_b20}" if a_b20 else "FAIL"
    a_b40 = bkz40_results.get(k1, '—')
    a_b40_str = (f"m={a_b40}" if a_b40 and a_b40 != '—' else ("FAIL" if a_b40 == None else "—"))
    c_lll = curve_c_actual.get(k1, curve_c_lll_ref.get(k1))
    c_str = f"m={c_lll}" if c_lll else "FAIL"
    print(f"{k1:>4} {eff:>6.4f} | {a_lll_str:>12} | {a_b20_str:>10} | {a_b40_str:>10} | {c_str:>8}")

print()
print("Interpretation:")

all_bkz20_succeed = all(v is not None for k, v in bkz20_results.items() if k in K1_VALUES)
if all_bkz20_succeed:
    print("  BKZ-20 succeeds at ALL K1 values for Curve A.")
    print("  => Reduction quality (not lattice structure) is the root cause.")
    print("  => The near-rational lam/n does NOT create an unsurmountable obstruction.")
    print("  => For secp256k1-scale attacks, BKZ (not LLL) should be used when lam≈n/3.")
else:
    failed_k1 = [k for k in K1_VALUES if bkz20_results.get(k) is None]
    print(f"  BKZ-20 FAILS for K1∈{failed_k1}.")
    all_bkz40_succeed = all(bkz40_results.get(k) is not None for k in failed_k1)
    if all_bkz40_succeed:
        print("  BKZ-40 rescues the remaining cases.")
        print("  => Higher block-size BKZ can overcome the near-rational lam/n challenge.")
    else:
        still_failed = [k for k in failed_k1 if bkz40_results.get(k) is None]
        print(f"  BKZ-40 also fails for K1∈{still_failed}.")
        print("  => Possible genuine structural obstruction from lam≈n/3 at high K1.")
        print("  => Would need BKZ(β>40) or a different lattice construction.")

print()
print("Done.")
