"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Prior data (from 2026-07-26 run):
  λ/n = 0.07  (p=2677, n=2647):  LLL + BKZ(40) FAIL at all m≤12
  λ/n = 0.34  (p=524347, n=523969): LLL succeeds at m=9 (3/3)
  λ/n = 0.53  (p=211, n=199):    LLL succeeds at m=4 (3/3)

Goal: bisect the failure threshold by scanning 12-bit curves with
λ/n ∈ [0.05, 0.45] in steps of ~0.05, running LLL at m=4..14.

Strategy:
  1. Scan primes p ≡ 1 (mod 3), p ≤ 2^13, using Eisenstein decomposition
     to compute all 6 group orders n without point counting.
  2. For each (p, n) with n prime and n ≡ 1 (mod 3), compute λ/n.
  3. Bucket curves by λ/n into bins [0.05,0.10), [0.10,0.15), ..., [0.45,0.50).
     Keep one representative per bin.
  4. For each representative, run LLL sweep m=4..14 with 3 seeds.
  5. Report threshold: smallest λ/n where LLL 3/3 succeeds.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (copied from glv_hnp_phase2_20bit.py)
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

def find_generator(p, b, n):
    """Find a generator of E(F_p): y^2=x^3+b of order n (n prime)."""
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

# ---------------------------------------------------------------------------
# CM theory for j=0 curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Smaller root of x^2 + x + 1 = 0 mod n. Requires n ≡ 1 (mod 3)."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def find_curve_b(p, n):
    """Find parameter b (1..p-1) such that #E(F_p): y^2=x^3+b = n."""
    rng = random.Random(999)
    for b in range(1, min(p, 200)):
        for _ in range(20):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                Q = ec_mul(P, n, p)
                if Q is None:
                    return b
                break
    return None

# ---------------------------------------------------------------------------
# Curve scan: collect one curve per λ/n bin
# ---------------------------------------------------------------------------

BIN_WIDTH = 0.05
BIN_EDGES = [i * BIN_WIDTH for i in range(1, 10)]  # 0.05, 0.10, ..., 0.45

def bin_index(ratio):
    """Return bin index 0..8 for ratio in [0.05, 0.50)."""
    idx = int(ratio / BIN_WIDTH)
    return max(0, min(idx, len(BIN_EDGES) - 1))

def scan_curves(p_max=2**14):
    """
    Scan primes p ≡ 1 (mod 3), p ≤ p_max.
    For each, compute all 6 group orders n via Eisenstein decomposition.
    Collect one (p, b, n, lam) per λ/n bin.
    """
    bins = {}  # bin_idx -> (p, b, n, lam, ratio)
    p = 7  # first prime ≡ 1 (mod 3)
    while p <= p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 7 or not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    if ratio < 0.05 or ratio >= 0.50:
                        continue
                    bi = bin_index(ratio)
                    if bi in bins:
                        continue  # already have this bin
                    # Find curve parameter b
                    b_param = find_curve_b(p, n_cand)
                    if b_param is None:
                        continue
                    G = find_generator(p, b_param, n_cand)
                    if G is None:
                        continue
                    bins[bi] = (p, b_param, n_cand, lam, ratio, G)
                    if len(bins) == len(BIN_EDGES):
                        return bins  # all bins filled
        p = sympy.nextprime(p)
    return bins

# ---------------------------------------------------------------------------
# Signature generation and GLV lattice (from glv_hnp_phase2_20bit.py)
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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret=None):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_secret is not None and d_cand == d_secret:
            return d_cand
    return None

def run_experiment(p, b, n, lam, G, m, d_secret, k1_bound, seed=42):
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep_one_curve(p, b, n, lam, ratio, G, m_range, seeds):
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok = run_experiment(p, b, n, lam, G, m, d_trial, k1_bound, seed)
            wins += ok
        results[m] = (wins, len(seeds))
    first_full = None
    for m, (w, t) in sorted(results.items()):
        if w == t:
            first_full = m
            break
    return results, first_full

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection")
print("Prior data:  λ/n=0.07 → FAIL;  λ/n=0.34 → succeed at m=9")
print("Goal:        find smallest λ/n where LLL 3/3 succeeds")
print("=" * 70)

SEEDS = [42, 1234, 9999]
M_RANGE = range(4, 15)

print("\nPhase 1: Scanning for one curve per λ/n bin (width=0.05)...")
bins = scan_curves(p_max=2**14)
print(f"  Found {len(bins)}/9 bins filled.")
for bi in sorted(bins):
    p, bparam, n, lam, ratio, G = bins[bi]
    print(f"  Bin {bi} (λ/n≈{BIN_EDGES[bi]:.2f}): p={p}, n={n} ({n.bit_length()}b), "
          f"lam={lam}, λ/n={ratio:.4f}")

print("\nPhase 2: LLL sweep for each bin (m=4..14, 3 seeds each)...")
print("-" * 70)

THRESHOLD_DATA = []
for bi in sorted(bins):
    p, bparam, n, lam, ratio, G = bins[bi]
    results, first_full = sweep_one_curve(p, bparam, n, lam, ratio, G,
                                          M_RANGE, SEEDS)
    THRESHOLD_DATA.append((ratio, n, first_full, results))
    status = f"LLL 3/3 at m={first_full}" if first_full else f"LLL FAIL (never 3/3 up to m={max(M_RANGE)-1})"
    row = "  ".join(f"m={m}:{w}/{t}" for m, (w, t) in sorted(results.items()))
    print(f"  λ/n={ratio:.4f} (p={p}, n={n}): {status}")
    print(f"    Detail: {row}")

# Also add known data points from prior runs
known = [
    (0.07, 2647, None, {}),   # p=2677, confirmed fail in 2026-07-26 run
    (0.34, 523969, 9, {}),    # p=524347, confirmed m=9 in 2026-07-26 run
    (0.53, 199, 4, {}),       # p=211, confirmed m=4 in 2026-07-26 run
]

print("\n" + "=" * 70)
print("SUMMARY TABLE (all curves, including prior-run reference points)")
print("=" * 70)
print(f"{'λ/n':>8}  {'n':>8}  {'n bits':>6}  {'LLL result':>30}")
print("-" * 70)

all_data = sorted(THRESHOLD_DATA + known, key=lambda x: x[0])
for ratio, n, first_full, _ in all_data:
    status = f"3/3 at m={first_full}" if first_full else "FAIL (never 3/3)"
    n_bits = n.bit_length()
    print(f"  {ratio:>6.4f}  {n:>8}  {n_bits:>6}b  {status:>30}")

# Find threshold
successes = [(r, mf) for r, n, mf, _ in all_data if mf is not None]
failures  = [(r, n) for r, n, mf, _ in all_data if mf is None]

if successes and failures:
    max_fail_ratio = max(r for r, _ in failures)
    min_succ_ratio = min(r for r, _ in successes)
    print(f"\nThreshold estimate: λ/n ∈ ({max_fail_ratio:.4f}, {min_succ_ratio:.4f})")
    print(f"  → Midpoint: {(max_fail_ratio + min_succ_ratio) / 2:.4f}")
elif successes:
    print(f"\nAll curves succeed. Min successful λ/n = {min(r for r,_ in successes):.4f}")
else:
    print("\nAll curves fail — threshold is above all tested λ/n values.")

print("\nDone.")
