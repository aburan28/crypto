"""
GLV-HNP Phase 2: lambda/n threshold bisection — controlled for bit length.

Problem with first attempt (glv_hnp_lambda_threshold.py):
  - FAIL curve was 12-bit (n=2647) but all SUCCESS curves were 11-bit.
  - Need to compare same-size curves to isolate lam/n as the variable.

This script:
  1. Finds 12-bit j=0 prime-order curves with lam/n in [0.07, 0.40].
  2. Groups by lam/n bucket, sweeps LLL, reports m* per bucket.
  3. Also tests the specific p=2677 curve at more seeds to confirm it fails.

Run: python3 glv_hnp_lambda_bisect.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic
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

# ---------------------------------------------------------------------------
# Lattice attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    att = 0
    while len(sigs) < m and att < 500000:
        att += 1
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
        sigs.append({'A': A, 'B': B})
    return sigs

def run_lll(p, b, n, lam, G, k1_bound, m, d_secret, seed):
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    dim = 2 * m + 2
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1  # S_D = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for row_i in range(dim):
        last = A[row_i][dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[row_i][m]) % n
        if d_cand == d_secret:
            return True
    return False

def sweep_curve(label, p, b, n, lam, G, k1_bound, m_range, seeds):
    ratio = lam / n
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    print(f"\n[{label}]  lam/n={ratio:.4f}  n={n}({n.bit_length()}b)  K1={k1_bound}"
          f"  K2={k2_bound}  eff={eff:.4f}")
    first_full = None
    last_result = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 5555).randint(1, n - 1)
            if run_lll(p, b, n, lam, G, k1_bound, m, d_trial, seed):
                wins += 1
        last_result[m] = (wins, len(seeds))
        marker = " *** FIRST 3/3" if wins == len(seeds) and first_full is None else ""
        print(f"  m={m:2d}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first_full is None:
            first_full = m
    return first_full, ratio

# ---------------------------------------------------------------------------
# Search for 12-bit curves across lam/n buckets
# ---------------------------------------------------------------------------

BUCKETS_12BIT = [
    (0.07, 0.10, "12b_0.07-0.10"),
    (0.10, 0.14, "12b_0.10-0.14"),
    (0.14, 0.18, "12b_0.14-0.18"),
    (0.18, 0.23, "12b_0.18-0.23"),
    (0.23, 0.30, "12b_0.23-0.30"),
    (0.30, 0.40, "12b_0.30-0.40"),
]

def find_12bit_curves_per_bucket():
    # 12-bit: p in [2048, 4096]
    buckets = {label: None for (_, _, label) in BUCKETS_12BIT}
    found = {label: False for (_, _, label) in BUCKETS_12BIT}

    print("Searching for 12-bit j=0 GLV prime-order curves by lam/n bucket...")
    p = sympy.nextprime(2047)
    limit = 4096
    while p < limit and not all(found.values()):
        p = sympy.nextprime(p)
        if p % 3 != 1: continue
        eis = eisenstein_decompose(p)
        if eis is None: continue
        a_e, b_e = eis
        traces = j0_traces(a_e, b_e)
        for t in traces:
            n_cand = p + 1 - t
            if n_cand < 2 or not sympy.isprime(n_cand): continue
            if n_cand.bit_length() != 12: continue  # strict 12-bit
            if n_cand % 3 != 1: continue
            lam, _ = glv_eigenvalue(n_cand)
            if lam is None: continue
            ratio = lam / n_cand
            for (lo, hi, label) in BUCKETS_12BIT:
                if not found[label] and lo <= ratio < hi:
                    # find b
                    found_b = None
                    for b_try in range(1, 300):
                        rng_b = random.Random(b_try * 17)
                        for _ in range(40):
                            x = rng_b.randint(0, p - 1)
                            rhs_x = (pow(x, 3, p) + b_try) % p
                            y_x = tonelli_shanks(rhs_x, p)
                            if y_x is not None and y_x != 0:
                                if ec_mul((x, y_x), n_cand, p) is None:
                                    found_b = b_try
                                    break
                        if found_b is not None: break
                    if found_b is None: continue
                    G_cand = find_generator(p, found_b, n_cand)
                    if G_cand is None: continue
                    buckets[label] = (p, found_b, n_cand, lam, G_cand, ratio)
                    found[label] = True
                    print(f"  {label}: p={p}, b={found_b}, n={n_cand}, lam={lam},"
                          f" ratio={ratio:.4f}")
    return buckets

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — lambda/n threshold bisection (controlled: 12-bit)")
print("=" * 70)

SEEDS = [42, 1234, 9999]
# For 12-bit curves, K1=8 (as in failing reference) is a natural choice
# But the failing curve also has K1=8. Let's use same K1 for fairness.
K1 = 8
M_RANGE = range(4, 16)

# Step 1: find curves
buckets = find_12bit_curves_per_bucket()

# Step 2: also include the known failing reference (p=2677, n=2647, lam=185, lam/n=0.07)
G_fail = find_generator(2677, 2, 2647)
known_fail = (2677, 2, 2647, 185, G_fail, 185/2647)
print(f"  [known_fail] p=2677, b=2, n=2647, lam=185, ratio=0.0699 (12b)")

# Step 3: sweep
print("\n" + "=" * 70)
print(f"Sweep: K1={K1}, K2=sqrt(n)+1, m=4..15, 3 seeds")
print("=" * 70)

summary = []

# Known fail
m_star, ratio = sweep_curve("known_FAIL p=2677", 2677, 2, 2647, 185, G_fail,
                            K1, M_RANGE, SEEDS)
summary.append((ratio, m_star, "known_FAIL"))

for (lo, hi, label) in BUCKETS_12BIT:
    entry = buckets.get(label)
    if entry is None:
        print(f"\n[SKIP] {label}: no curve found")
        summary.append(((lo+hi)/2, None, label + "_notfound"))
        continue
    p, b, n, lam, G, ratio = entry
    m_star, ratio = sweep_curve(label, p, b, n, lam, G, K1, M_RANGE, SEEDS)
    summary.append((ratio, m_star, label))

# Step 4: report
print("\n" + "=" * 70)
print("SUMMARY — 12-bit curves, K1=8, 3 seeds")
print("=" * 70)
print(f"{'lam/n':>10}  {'m*':>5}  {'label'}")
print("-" * 45)
for ratio, m_star, label in sorted(summary, key=lambda x: x[0]):
    m_str = str(m_star) if m_star is not None else "FAIL(>15)"
    print(f"{ratio:10.4f}  {m_str:>9}  {label}")

sorted_s = sorted(summary, key=lambda x: x[0])
threshold_lo = None
threshold_hi = None
for ratio, m_star, label in sorted_s:
    if m_star is None:
        threshold_lo = ratio
    elif threshold_lo is not None and threshold_hi is None:
        threshold_hi = ratio
        break

print()
if threshold_lo is not None and threshold_hi is not None:
    print(f"** Threshold in ({threshold_lo:.4f}, {threshold_hi:.4f}) at fixed 12-bit / K1={K1}")
elif threshold_lo is not None:
    print(f"** All curves above {threshold_lo:.4f} succeed for 12-bit / K1={K1}")
else:
    print(f"** All tested curves succeed (threshold is below {sorted_s[0][0]:.4f})")

print("\nDone.")
