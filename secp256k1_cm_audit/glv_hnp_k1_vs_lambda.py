"""
GLV-HNP Phase 2 — K1_BOUND vs λ/n: targeted decoupling experiment.

From lambda_threshold.py we found:
  n=1303, λ/n=0.073, K1=2 → PASS at m=4
  n=2647, λ/n=0.070, K1=8 → FAIL at all m≤14

The λ/n ratios are nearly equal. We test whether the failure is due to:
  (a) K1_BOUND=8 (not λ/n=0.07)
  (b) The specific n=2647 magnitude
  (c) Some combination

Experiments:
  E1: n=2647, K1=2  → if PASS, K1 is the culprit
  E2: n=1303, K1=8  → if FAIL, K1_BOUND is the culprit for small-n too
  E3: n=2647, K1=4  → binary search on K1
  E4: n=2647, K1=6  → further bisect
  E5: n=1303, K1=4  → same for small n

Secondary: does a large-K1 curve with large λ/n still fail?
  E6: n=1447, λ/n=0.487, K1=8 → if PASS, confirms λ/n threshold within K1=8 class
"""

import math
import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass a=0)
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
# Signatures and lattice
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
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret: return True
    return False

SEEDS = [42, 1234, 9999]

def sweep(label, curve, k1_bound, m_range):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = (math.ceil(math.log(n) / math.log(1.0/eff)) if eff < 1 else 9999)
    print(f"\n{'='*65}")
    print(f"{label}")
    print(f"  p={p}, n={n}, λ={lam}, λ/n={lam/n:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")
    first_3of3 = None
    for m in m_range:
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 9999).randint(1, n - 1)
            wins += run_lll(curve, m, d_trial, k1_bound, seed)
        thresh_mark = " ← thresh" if m == m_thresh else ""
        ok_mark = " ← FIRST 3/3" if wins == 3 and first_3of3 is None else ""
        print(f"  m={m:2d}: {wins}/3{thresh_mark}{ok_mark}")
        if wins == 3 and first_3of3 is None:
            first_3of3 = m
        # Stop after confirming 3/3 twice
        if first_3of3 and m >= first_3of3 + 1:
            break
    status = f"PASS at m={first_3of3}" if first_3of3 else f"FAIL (never 3/3 in m={m_range.start}..{m_range.stop-1})"
    print(f"  → {status}")
    return first_3of3

# ---------------------------------------------------------------------------
# Build the test curves
# ---------------------------------------------------------------------------

# Known failure: p=2677, n=2647, λ=185, λ/n=0.0699
G_2647 = find_generator(2677, 2, 2647)
curve_2647 = (2677, 2, 2647, 185, G_2647)

# Known near-match: n=1303, λ=95, λ/n=0.0729 (from lambda_threshold.py)
# p=1249 has Eisenstein trace giving n=1303, b=?
# Let's find b for p=1249, n=1303
G_1303 = None
b_1303 = None
lam_1303 = 95  # from lambda_threshold.py

# Verify λ: (95^2 + 95 + 1) % 1303
assert (95*95 + 95 + 1) % 1303 == 0, "λ check failed for n=1303"

# Find a generator of order 1303 on some j=0 curve over F_1249
# From the bucket output: p=1249, n=1303
for b_try in range(1, 300):
    rng_tmp = random.Random(b_try)
    found_pt = None
    for _ in range(500):
        x = rng_tmp.randint(0, 1248)
        rhs = (pow(x, 3, 1249) + b_try) % 1249
        y = tonelli_shanks(rhs, 1249)
        if y and y != 0:
            found_pt = (x, y)
            break
    if found_pt and ec_mul(found_pt, 1303, 1249) is None:
        b_1303 = b_try
        G_1303 = find_generator(1249, b_try, 1303)
        break

if G_1303 is None:
    print("ERROR: Could not build n=1303 curve")
else:
    print(f"Built n=1303 curve: p=1249, b={b_1303}, G={G_1303}")

curve_1303 = (1249, b_1303, 1303, lam_1303, G_1303)

# n=1447, λ=704, λ/n=0.487  (from bucket [0.45,0.50))
# p=1483 from lambda_threshold output
lam_1447 = 704
assert (704*704 + 704 + 1) % 1447 == 0, "λ check failed for n=1447"
b_1447 = None
G_1447 = None
for b_try in range(1, 300):
    rng_tmp = random.Random(b_try)
    found_pt = None
    for _ in range(500):
        x = rng_tmp.randint(0, 1482)
        rhs = (pow(x, 3, 1483) + b_try) % 1483
        y = tonelli_shanks(rhs, 1483)
        if y and y != 0:
            found_pt = (x, y)
            break
    if found_pt and ec_mul(found_pt, 1447, 1483) is None:
        b_1447 = b_try
        G_1447 = find_generator(1483, b_try, 1447)
        break

if G_1447 is None:
    print("ERROR: Could not build n=1447 curve")
else:
    print(f"Built n=1447 curve: p=1483, b={b_1447}, G={G_1447}")

curve_1447 = (1483, b_1447, 1447, lam_1447, G_1447)

# ---------------------------------------------------------------------------
# Run decoupling experiments
# ---------------------------------------------------------------------------

print("\n" + "=" * 65)
print("DECOUPLING EXPERIMENT: K1_BOUND vs λ/n")
print("Known: n=2647, λ/n=0.070, K1=8 → FAIL")
print("Known: n=1303, λ/n=0.073, K1=2 → PASS at m=4")
print("=" * 65)

results = {}

# E1: n=2647, K1=2 (same λ/n=0.070, reduce K1 to 2)
if G_2647:
    r = sweep("E1: n=2647, λ/n=0.070, K1=2 (lower K1 on failing curve)",
              curve_2647, k1_bound=2, m_range=range(3, 16))
    results['E1'] = r

# E2: n=1303, K1=8 (same λ/n=0.073, raise K1 to 8)
if G_1303:
    r = sweep("E2: n=1303, λ/n=0.073, K1=8 (raise K1 on passing curve)",
              curve_1303, k1_bound=8, m_range=range(3, 16))
    results['E2'] = r

# E3: n=2647, K1=4 (binary search between K1=2 pass and K1=8 fail)
if G_2647:
    r = sweep("E3: n=2647, λ/n=0.070, K1=4 (bisect)",
              curve_2647, k1_bound=4, m_range=range(3, 16))
    results['E3'] = r

# E4: n=2647, K1=6 (tighten bisect)
if G_2647:
    r = sweep("E4: n=2647, λ/n=0.070, K1=6",
              curve_2647, k1_bound=6, m_range=range(3, 16))
    results['E4'] = r

# E5: n=1303, K1=4
if G_1303:
    r = sweep("E5: n=1303, λ/n=0.073, K1=4",
              curve_1303, k1_bound=4, m_range=range(3, 16))
    results['E5'] = r

# E6: n=1447, λ/n=0.487, K1=8 (large λ/n, large K1 — expect PASS)
if G_1447:
    r = sweep("E6: n=1447, λ/n=0.487, K1=8 (large λ/n, large K1)",
              curve_1447, k1_bound=8, m_range=range(3, 16))
    results['E6'] = r

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print("\n" + "=" * 65)
print("DECOUPLING SUMMARY")
print("=" * 65)

rows = [
    ("Known-FAIL", "n=2647", "λ/n=0.070", 8, None, "FAIL baseline"),
    ("Known-PASS", "n=1303", "λ/n=0.073", 2, 4,    "PASS baseline"),
    ("E1",         "n=2647", "λ/n=0.070", 2, results.get('E1'), "K1↓ on FAIL curve"),
    ("E2",         "n=1303", "λ/n=0.073", 8, results.get('E2'), "K1↑ on PASS curve"),
    ("E3",         "n=2647", "λ/n=0.070", 4, results.get('E3'), "K1 bisect"),
    ("E4",         "n=2647", "λ/n=0.070", 6, results.get('E4'), "K1 bisect"),
    ("E5",         "n=1303", "λ/n=0.073", 4, results.get('E5'), "K1 bisect"),
    ("E6",         "n=1447", "λ/n=0.487", 8, results.get('E6'), "large λ/n, K1=8"),
]

print(f"  {'Exp':12s}  {'n':8s}  {'λ/n':10s}  {'K1':4s}  {'min m':7s}  {'Desc'}")
print(f"  {'-'*12}  {'-'*8}  {'-'*10}  {'-'*4}  {'-'*7}  {'-'*25}")
for exp, n_str, ratio_str, k1, min_m, desc in rows:
    status = f"m={min_m}" if min_m is not None else "FAIL"
    print(f"  {exp:12s}  {n_str:8s}  {ratio_str:10s}  {k1:4d}  {status:7s}  {desc}")

print("\nConclusion:")
if results.get('E1') and not results.get('E2'):
    print("  → K1_BOUND is the primary driver: K1=2 rescues the failing curve,")
    print("    K1=8 breaks the passing curve. λ/n is NOT the threshold variable.")
elif results.get('E1') and results.get('E2'):
    print("  → Both E1 and E2 pass: K1 alone doesn't explain failure.")
    print("    Something else (n magnitude, or K1 * λ/n interaction) drives it.")
elif not results.get('E1') and results.get('E2'):
    print("  → E1 fails and E2 fails: consistent with K1 being the culprit.")
    print("    With K1=2 the n=2647 curve ALSO fails → n magnitude matters too.")
else:
    print("  → E1 fails, E2 passes: λ/n is the driver, not K1.")

print("\nDone.")
