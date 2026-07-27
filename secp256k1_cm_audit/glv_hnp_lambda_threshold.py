"""
GLV-HNP λ/n threshold study (Thread 20).

Bisects the failure threshold between λ/n=0.07 (known FAIL, p=2677)
and λ/n=0.25 (known SUCCESS, 20-bit curves).

Finds j=0 curves in 7 λ/n bins and tests LLL attack for each,
reporting the first m at which 3/3 seeds recover d.

Usage: python3 glv_hnp_lambda_threshold.py
"""

import math
import random
import sys
import sympy
from fpylll import IntegerMatrix, LLL

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
    mv, cv, tv, rv = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if tv == 0: return 0
        if tv == 1: return rv
        i, tmp = 0, tv
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(cv, 1 << (mv - i - 1), p)
        mv, cv, tv, rv = i, b * b % p, tv * b * b % p, rv * b % p

def find_generator(p, b, n):
    rng = random.Random(42)
    for _ in range(60000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM helpers: Eisenstein decomposition and j0 traces
# ---------------------------------------------------------------------------

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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

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

def find_curve_for_n(p, n):
    """Find b such that #(y^2=x^3+b / F_p) = n. Try b=1..100."""
    rng = random.Random(77)
    for b in range(1, 101):
        for _ in range(300):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b
                break  # this b doesn't give order n; move to next b
    return None

# ---------------------------------------------------------------------------
# GLV lattice (same column-diagonal scaling as Phase 2 experiments)
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs); dim = 2 * m + 2
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

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 400000:
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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def run_lll(p, b, n, lam, G, k1_bound, m, seed=42):
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 4444).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Curve discovery: find one j=0 curve per λ/n bin
# ---------------------------------------------------------------------------

BINS = [
    (0.04, 0.09, 'λ/n∈[0.04,0.09)'),   # known FAIL anchor
    (0.09, 0.13, 'λ/n∈[0.09,0.13)'),
    (0.13, 0.17, 'λ/n∈[0.13,0.17)'),
    (0.17, 0.21, 'λ/n∈[0.17,0.21)'),
    (0.21, 0.25, 'λ/n∈[0.21,0.25)'),
    (0.25, 0.30, 'λ/n∈[0.25,0.30)'),
    (0.30, 0.50, 'λ/n∈[0.30,0.50)'),   # known SUCCESS anchor
]

print("=" * 70)
print("GLV-HNP λ/n threshold study")
print("Bisecting the LLL failure boundary in [0.07, 0.25]")
print("=" * 70)
print()

# Search 10..14 bit primes for one representative per bin
bin_curves = {}
p = sympy.nextprime(2**9)
p_max = 2**14
print(f"Searching for j=0 curves in λ/n bins (10-14 bit, p ≤ {p_max})...")
sys.stdout.flush()

found_count = 0
target_count = len(BINS)

while p <= p_max and found_count < target_count:
    if p % 3 == 1:
        eis = eisenstein_decompose(p)
        if eis is not None:
            ae, be_eis = eis
            traces = j0_traces(ae, be_eis)
            for t in traces:
                n_cand = p + 1 - t
                if n_cand < 64 or n_cand > 2 * p: continue
                if not sympy.isprime(n_cand): continue
                if n_cand % 3 != 1: continue
                lam, _ = glv_eigenvalue(n_cand)
                if lam is None: continue
                ratio = lam / n_cand
                for lo, hi, label in BINS:
                    if lo not in bin_curves and lo <= ratio < hi:
                        b_found = find_curve_for_n(p, n_cand)
                        if b_found is None: continue
                        G = find_generator(p, b_found, n_cand)
                        if G is None: continue
                        bin_curves[lo] = {
                            'p': p, 'b': b_found, 'n': n_cand,
                            'lam': lam, 'G': G, 'label': label, 'ratio': ratio
                        }
                        print(f"  [FOUND] {label}: p={p}, n={n_cand}, "
                              f"λ={lam}, λ/n={ratio:.4f}")
                        sys.stdout.flush()
                        found_count += 1
                        break
    p = sympy.nextprime(p)

print(f"\nFound {found_count}/{target_count} bins.\n")

if found_count < target_count:
    missing = [label for lo, hi, label in BINS if lo not in bin_curves]
    print(f"WARNING: missing bins: {missing}")
    print()

# ---------------------------------------------------------------------------
# LLL sweep for each bin
# ---------------------------------------------------------------------------

K1_BOUND = 2    # bias: k1 ∈ {0, 1}
SEEDS = [42, 1234, 9999]
M_RANGE = list(range(3, 12))

print(f"Running LLL sweep: K1_BOUND={K1_BOUND}, seeds={SEEDS}, m={M_RANGE[0]}..{M_RANGE[-1]}")
print()

col_hdr = "λ/n bin               | lam/n | λ     |  n    |  p    " + \
          "".join(f"| m={m} " for m in M_RANGE) + "| status"
print(col_hdr)
print("-" * len(col_hdr))

threshold_data = []

for lo, hi, label in BINS:
    if lo not in bin_curves:
        print(f"{label:22s} | NOT FOUND")
        continue

    cd = bin_curves[lo]
    p_c, b_c, n_c, lam_c, G_c, ratio = (
        cd['p'], cd['b'], cd['n'], cd['lam'], cd['G'], cd['ratio']
    )

    wins_per_m = {}
    for m in M_RANGE:
        wins = 0
        for seed in SEEDS:
            ok = run_lll(p_c, b_c, n_c, lam_c, G_c, K1_BOUND, m, seed)
            wins += ok
        wins_per_m[m] = wins
        sys.stdout.flush()

    first_full = next((m for m in M_RANGE if wins_per_m[m] == len(SEEDS)), None)
    any_success = any(w > 0 for w in wins_per_m.values())
    if first_full is not None:
        status = f"SUCCESS (3/3 first m={first_full})"
    elif any_success:
        status = "PARTIAL"
    else:
        status = "FAIL"

    threshold_data.append((ratio, label, status, first_full, wins_per_m))

    row = (f"{label:22s} | {ratio:.3f} | {lam_c:5d} | {n_c:5d} | {p_c:5d} " +
           "".join(f"| {wins_per_m[m]}/3 " for m in M_RANGE) +
           f"| {status}")
    print(row)
    sys.stdout.flush()

# ---------------------------------------------------------------------------
# Summary and threshold determination
# ---------------------------------------------------------------------------

print()
print("=" * 70)
print("THRESHOLD ANALYSIS SUMMARY")
print("=" * 70)
print(f"K1_BOUND={K1_BOUND}, K2_BOUND=ceil(sqrt(n)), 3 seeds per (m, curve)")
print()

fail_ratios = [r for r, lbl, s, _, _ in threshold_data if s == "FAIL"]
success_ratios = [r for r, lbl, s, _, _ in threshold_data if s == "SUCCESS (3/3 first m=3)" or "SUCCESS" in s]
partial_ratios = [r for r, lbl, s, _, _ in threshold_data if s == "PARTIAL"]

for ratio, label, status, first_full, _ in sorted(threshold_data):
    print(f"  λ/n={ratio:.4f}  {label:22s}  →  {status}")

print()
if fail_ratios and success_ratios:
    max_fail = max(fail_ratios)
    min_success = min(success_ratios)
    print(f"Threshold bracket: FAIL ≤ {max_fail:.3f}, SUCCESS ≥ {min_success:.3f}")
    if partial_ratios:
        print(f"Partial results at: {sorted(partial_ratios)}")
    print(f"Critical boundary: ~{(max_fail + min_success) / 2:.3f} (midpoint)")
elif fail_ratios:
    print(f"All failures at ratios: {sorted(fail_ratios)}")
elif success_ratios:
    print(f"All successes — threshold is below min tested ({min(success_ratios):.3f})")

print()
print("NOTE: Using λ = min-root of x²+x+1=0 (mod n). Failure may be rescued")
print("      by using the OTHER root (λ'=n-1-λ) which has λ'/n≥0.5 always.")
print()
print("Done.")
