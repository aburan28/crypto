"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Goal: Empirically locate the λ/n threshold below which the GLV-aware
HNP lattice attack (LLL, column-diagonal scaling) fails to recover d.

Prior data (from 2026-07-26 runs):
  λ/n=0.07 (p=2677, n=2647): LLL+BKZ(40) never 3/3 → FAIL
  λ/n=0.34 (p=524347, n=523969): LLL 3/3 at m=9 → PASS
  λ/n=0.53 (p=211,   n=199):    LLL 3/3 at m=4 → PASS (8-bit ref)
  λ/n=0.66 (p=2557,  n=2659):   LLL 3/3 at m=7 → PASS (12-bit ref)

Strategy: find curves with λ/n in bins [0.08,0.12), [0.12,0.17),
[0.17,0.22), [0.22,0.27), [0.27,0.32), all 12-14 bit for speed.
Run LLL sweep m=4..14. Report pass/fail per bin.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ─────────────────────────────────────────────────────────────────────────────
# EC arithmetic (short Weierstrass a=0)
# ─────────────────────────────────────────────────────────────────────────────

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

# ─────────────────────────────────────────────────────────────────────────────
# Eisenstein CM machinery for j=0 curves over F_p
# ─────────────────────────────────────────────────────────────────────────────

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                bv = num // 2
                if bv >= 0 and a * a - a * bv + bv * bv == p:
                    return (a, bv)
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

# ─────────────────────────────────────────────────────────────────────────────
# Curve finder: search 12-15 bit primes for target λ/n ∈ [lo, hi)
# ─────────────────────────────────────────────────────────────────────────────

def find_point_fast(p, b, seed=42):
    """Find a random non-identity point on y^2 = x^3 + b over F_p."""
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_curve_in_ratio_bin(lo, hi, bit_lo=11, bit_max=15, max_tries=200000):
    """Return (p, b, n, lam, G) with λ/n ∈ [lo, hi) and n in [2^bit_lo, 2^bit_max)."""
    p = sympy.nextprime(2**bit_lo)
    tried = 0
    while p < 2**bit_max and tried < max_tries:
        tried += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                ae, be = eis
                for t in j0_traces(ae, be):
                    nc = p + 1 - t
                    if nc < 2 or nc.bit_length() < bit_lo:
                        continue
                    if sympy.isprime(nc) and nc % 3 == 1:
                        lam, _ = glv_eigenvalue(nc)
                        if lam is None: continue
                        ratio = lam / nc
                        if lo <= ratio < hi:
                            # Find twist b using random point sampling (O(1) per b)
                            for b_try in range(1, 200):
                                pt = find_point_fast(p, b_try, seed=b_try * 1337)
                                if pt is None: continue
                                if ec_mul(pt, nc, p) is None:
                                    G = find_generator(p, b_try, nc)
                                    if G is not None:
                                        return (p, b_try, nc, lam, G)
                            break
        p = sympy.nextprime(p)
    return None

# ─────────────────────────────────────────────────────────────────────────────
# Lattice attack machinery (from glv_hnp_phase2_20bit.py)
# ─────────────────────────────────────────────────────────────────────────────

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
    M[m][m] = 1             # S_D = 1
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
            return True
    return False

def run_lll_trial(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def sweep(label, curve, k1_bound, m_range, seeds):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    if eff < 1.0:
        m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))
    else:
        m_thresh = 0
    print(f"\n  [{label}]  p={p}, n={n} ({n.bit_length()}b), lam={lam}, lam/n={lam/n:.4f}")
    print(f"           K1={k1_bound}, K2={k2_bound}, eff={eff:.5f}, m_thresh≈{m_thresh:.1f}")
    found_m = None
    for m in m_range:
        wins = sum(
            run_lll_trial(curve, m, random.Random(seed + 9999).randint(1, n-1), k1_bound, seed)
            for seed in seeds
        )
        total = len(seeds)
        print(f"    m={m}: {wins}/{total}", end="", flush=True)
        if wins == total and found_m is None:
            found_m = m
            print(" ← PASS", end="")
        print()
    return found_m   # None = never passed

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection study")
print("Prior data: 0.07 FAIL, 0.34 PASS → bisect [0.07, 0.34]")
print("=" * 70)

# Bins to probe (lo, hi, label)
BINS = [
    (0.07, 0.11, "bin_A [0.07,0.11)"),
    (0.11, 0.15, "bin_B [0.11,0.15)"),
    (0.15, 0.20, "bin_C [0.15,0.20)"),
    (0.20, 0.25, "bin_D [0.20,0.25)"),
    (0.25, 0.30, "bin_E [0.25,0.30)"),
    (0.30, 0.35, "bin_F [0.30,0.35)"),
]

SEEDS = [42, 1234, 9999]
M_RANGE = range(4, 18)

results = {}

for lo, hi, label in BINS:
    print(f"\n{'─'*70}")
    print(f"Searching for curve with λ/n ∈ [{lo:.2f}, {hi:.2f})...")
    curve = find_curve_in_ratio_bin(lo, hi)
    if curve is None:
        print(f"  WARNING: no curve found in bin [{lo:.2f},{hi:.2f}). Skipping.")
        results[label] = None
        continue
    p, b, n, lam, G = curve
    # K1_BOUND=8 to match the original 12-bit experiment regime:
    #   (p=2557/2677, n≈2600, K1=8, K2≈51, eff≈0.157)
    # The failure at lam/n=0.07 and pass at lam/n=0.66 both used K1=8.
    k1_bound = 8
    found_m = sweep(label, curve, k1_bound, M_RANGE, SEEDS)
    results[label] = found_m

# Known reference points
print(f"\n{'─'*70}")
print("Reference from prior runs (fixed, not re-run here):")
print("  lam/n=0.07 (p=2677, n=2647): LLL never 3/3 even at m=12 → FAIL")
print("  lam/n=0.34 (p=524347, n=523969): LLL 3/3 at m=9 → PASS")
print("  lam/n=0.53 (p=211, n=199): LLL 3/3 at m=4 → PASS")

print(f"\n{'='*70}")
print("SUMMARY — λ/n threshold bisection")
print(f"{'='*70}")
for lo, hi, label in BINS:
    fm = results.get(label)
    if fm is None:
        verdict = "NO CURVE / skipped"
    elif fm == "no_curve":
        verdict = "NO CURVE"
    else:
        verdict = f"PASS at m={fm}" if fm is not None else f"FAIL (never 3/3 at m≤{max(M_RANGE)-1})"
    print(f"  {label}: {verdict}")

print()
pass_bins = [(lo, hi, label) for lo, hi, label in BINS
             if results.get(label) is not None and not isinstance(results.get(label), str)]
fail_bins = [(lo, hi, label) for lo, hi, label in BINS
             if results.get(label) is None or results.get(label) == "no_curve"]

# Find the transition
all_sorted = [(lo, hi, label, results.get(label)) for lo, hi, label in BINS]
threshold_lo, threshold_hi = None, None
for i in range(len(all_sorted) - 1):
    r_cur = all_sorted[i][3]
    r_nxt = all_sorted[i+1][3]
    if r_cur is None and r_nxt is not None:
        threshold_lo = all_sorted[i][0]
        threshold_hi = all_sorted[i+1][1]

if threshold_lo is not None:
    print(f"Threshold estimate: λ/n ∈ [{threshold_lo:.2f}, {threshold_hi:.2f})")
else:
    print("Threshold: cannot determine from this sweep (all pass or all fail).")
    # Report what we have
    for lo, hi, label, fm in all_sorted:
        status = "PASS" if (fm is not None and fm != "no_curve") else "FAIL/skip"
        print(f"  [{lo:.2f},{hi:.2f}): {status}")

print("\nDone.")
