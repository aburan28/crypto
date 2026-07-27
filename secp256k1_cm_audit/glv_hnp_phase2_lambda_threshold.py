"""
GLV-HNP Phase 2: λ/n threshold bisection.

From 2026-07-26 Phase 2 analysis (Thread 5 / Thread 20 proposal):
  - λ/n ≈ 0.44 (p=211, n=199, K1=2): LLL succeeds at m=6
  - λ/n ≈ 0.44 (20-bit curves, K1=8): LLL succeeds
  - λ/n = 0.07 (p=2677, n=2647, lam=185, K1=8): LLL+BKZ(40) both fail
    → structural failure (not a reduction-strength issue)

Hypothesis: there exists a critical λ*/n below which the Phase 2 lattice attack
always fails (for fixed K1). The mechanism (from log 2026-07-26):
  The -λ·S_K1 entries (k2-rows) can't be distinguished from the n·S_K1 entries
  (modular rows) when λ/n is small; the GS geometry degenerates.

Goal: sweep λ/n ∈ [0.05, 0.50] with K1=8, find the success/failure boundary.

Method:
  1. Scan j=0 CM primes (p ≡ 1 mod 3, Eisenstein form) in [2^11, 2^17]
  2. Bucket by λ/n into 0.05-wide bins
  3. Attack each representative with Phase 2 lattice (K1=8) at m=6..14, 3 seeds
  4. Report table: λ/n bin → first_m_with_3/3 or "fail"

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (a=0 short Weierstrass)
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
    mm, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b_param, n):
    rng = random.Random(12345)
    for _ in range(100000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_param) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def find_b_for_n(p, n):
    rng = random.Random(99)
    for b_try in range(1, min(p, 200)):
        for _ in range(20):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b_try
                break
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition + GLV eigenvalue
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a, b2) with a^2 - a*b2 + b2^2 = p, a > 0, b2 >= 0."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b2 = num // 2
                if b2 >= 0 and a * a - a * b2 + b2 * b2 == p:
                    return (a, b2)
    return None

def j0_traces(a, b2):
    """Six possible Frobenius traces for a j=0 CM curve with Eisenstein params."""
    return [2*a - b2, -(2*a - b2), b2 - 2*a,
            a + b2, -(a + b2), 2*b2 - a]

def glv_eigenvalue(n):
    """Smaller root of x^2 + x + 1 ≡ 0 (mod n). Requires n ≡ 1 (mod 3)."""
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

# ---------------------------------------------------------------------------
# Curve scan: bucket j=0 CM curves by λ/n
# ---------------------------------------------------------------------------

BIN_WIDTH = 0.05
BINS = [round(lo, 2) for lo in [i * BIN_WIDTH for i in range(1, 10)]]  # 0.05..0.45
MAX_PER_BIN = 2

def scan_curves(p_lo, p_hi, max_per_bin=MAX_PER_BIN, verbose=False):
    """
    Scan primes p in [p_lo, p_hi], find j=0 CM curves, bucket by λ/n.
    Returns dict: bin_lo -> list of (p, b, n, lam, lam_n).
    """
    buckets = {b: [] for b in BINS}
    p = sympy.nextprime(p_lo - 1)
    total_scanned = 0
    while p < p_hi:
        total_scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for tr in j0_traces(a_e, b_e):
                    n_cand = p + 1 - tr
                    if n_cand < 50:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    lam_n = lam / n_cand
                    bin_lo = math.floor(lam_n / BIN_WIDTH) * BIN_WIDTH
                    bin_lo = round(bin_lo, 2)
                    if bin_lo not in buckets:
                        continue
                    if len(buckets[bin_lo]) >= max_per_bin:
                        continue
                    b_param = find_b_for_n(p, n_cand)
                    if b_param is None:
                        continue
                    buckets[bin_lo].append((p, b_param, n_cand, lam, lam_n))
                    if verbose:
                        print(f"  [bin {bin_lo:.2f}] p={p}, n={n_cand}, "
                              f"lam={lam}, lam/n={lam_n:.4f}")
        p = sympy.nextprime(p)
        if all(len(v) >= max_per_bin for v in buckets.values()):
            break
    print(f"  Scanned {total_scanned} primes in [{p_lo}, {p_hi}]")
    filled = sum(1 for v in buckets.values() if v)
    print(f"  Filled {filled}/{len(BINS)} bins")
    return buckets

# ---------------------------------------------------------------------------
# Phase 2 attack with K1 bias
# ---------------------------------------------------------------------------

K1_BOUND = 8   # k1 in {0, ..., K1_BOUND-1}

def gen_sigs(p, b_param, n, lam, G, m, seed):
    K2_BOUND = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 300000:
        attempts += 1
        k1 = rng.randint(0, K1_BOUND - 1)
        k2 = rng.randint(0, K2_BOUND - 1)
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
        sigs.append({'A': A, 'B': B})
    return d_secret, sigs

def build_lattice(sigs, n, lam):
    m = len(sigs)
    K2_BOUND = math.isqrt(n) + 1
    dim = 2 * m + 2
    S_K1 = max(1, n // K1_BOUND)
    S_K2 = max(1, n // K2_BOUND)
    S_KANNAN = n

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
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_one(p, b_param, n, lam, G, m, seed):
    d_secret, sigs = gen_sigs(p, b_param, n, lam, G, m, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep_curve(p, b_param, n, lam, G, m_range, seeds):
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_one(p, b_param, n, lam, G, m, s))
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 72)
print("GLV-HNP Phase 2: λ/n threshold bisection (Thread 20)")
print(f"K1_BOUND={K1_BOUND}, m_range=6..16, 3 seeds per point")
print("=" * 72)
print()

# ---- Reference data points -------------------------------------------------
REFERENCE = [
    # (p, b, n, lam, lam_n, source)
    (211,  2,  199,   106, 106/199,  "toy 8-bit  [SUCCESS k/n≈0.53]"),
    (2677, 2, 2647,   185, 185/2647, "12-bit p=2677 [FAIL k/n≈0.07]"),
]

# The 8-bit toy uses the LARGER eigenvalue (λ=106 > 199/2=99.5 → wait, 106 > 99.5 so
# this is actually the larger root). For comparison we use the same convention.
# glv_eigenvalue returns min root; we note both for the toy.
# For p=211, n=199: roots are 92 and 106 (92+106=198=n-1 ✓). min=92, 92/199≈0.462.
# The toy attack uses λ=106. We use WHICHEVER eigenvalue corresponds to the attack.
# For this sweep we always use lam = glv_eigenvalue(n) = min root.

print("Reference points:")
for p, b, n, lam, lam_n, label in REFERENCE:
    print(f"  p={p}, n={n}, lam={lam}, lam/n={lam_n:.4f}  [{label}]")
print()

# ---- Curve scan ------------------------------------------------------------
print("Scanning j=0 CM curves in [2^11, 2^17] for λ/n bins...")
P_LO, P_HI = 2**11, 2**17
buckets = scan_curves(P_LO, P_HI, verbose=False)
print()

# ---- Print what we found ---------------------------------------------------
print("Curves found per bin:")
for bin_lo in sorted(buckets.keys()):
    curves = buckets[bin_lo]
    entries = [(f"n={n},lam={lam},lam/n={lam_n:.4f}") for (p,b,n,lam,lam_n) in curves]
    print(f"  [{bin_lo:.2f}, {bin_lo+BIN_WIDTH:.2f}): {len(curves)} curve(s)  "
          + ("  ".join(entries) if entries else "NONE"))
print()

# ---- LLL sweep per bin -----------------------------------------------------
M_RANGE = list(range(6, 17, 2))   # m = 6, 8, 10, 12, 14, 16
SEEDS = [42, 1234, 9999]
N_SEEDS = len(SEEDS)

print("=" * 72)
print(f"LLL sweep: K1={K1_BOUND}, m={M_RANGE[0]}..{M_RANGE[-1]} step 2, {N_SEEDS} seeds")
print("=" * 72)
print()

summary = []   # (bin_lo, lam_n, n, first_win_m or None)

for bin_lo in sorted(buckets.keys()):
    curves = buckets[bin_lo]
    if not curves:
        print(f"[λ/n ∈ [{bin_lo:.2f}, {bin_lo+BIN_WIDTH:.2f})]: NO CURVES FOUND — skip")
        summary.append((bin_lo, None, None, None))
        continue

    # Take only the first curve per bin for speed
    p, b_param, n, lam, lam_n = curves[0]
    G = find_generator(p, b_param, n)
    if G is None:
        print(f"[λ/n ∈ [{bin_lo:.2f}, {bin_lo+BIN_WIDTH:.2f})]: generator fail — skip")
        summary.append((bin_lo, lam_n, n, None))
        continue

    print(f"[λ/n ∈ [{bin_lo:.2f}, {bin_lo+BIN_WIDTH:.2f})]: "
          f"p={p}, n={n}, lam={lam}, lam/n={lam_n:.4f}")

    res = sweep_curve(p, b_param, n, lam, G, M_RANGE, SEEDS)
    first_win = next((m for m in M_RANGE if res[m][0] == N_SEEDS), None)

    row = "  "
    for m in M_RANGE:
        w, t = res[m]
        row += f"m={m}:{w}/{t} "
    print(row)
    if first_win:
        print(f"  → SUCCESS: first 3/3 at m={first_win}")
    else:
        max_w = max(w for w, t in res.values())
        print(f"  → FAIL: never 3/3 (max {max_w}/{N_SEEDS} at any m)")
    print()
    summary.append((bin_lo, lam_n, n, first_win))

# ---- Also run the two reference points -------------------------------------
print("=" * 72)
print("Reference point sweep (for calibration):")
print()

ref_summary = []
for (p, b, n, lam, lam_n, label) in REFERENCE:
    # For ref points use the stored lam (may be larger root, as in toy)
    G = find_generator(p, b, n)
    if G is None:
        print(f"  {label}: generator not found")
        continue
    print(f"  {label}  [p={p}, n={n}, lam={lam}, lam/n={lam_n:.4f}]")
    res = sweep_curve(p, b, n, lam, G, M_RANGE, SEEDS)
    first_win = next((m for m in M_RANGE if res[m][0] == N_SEEDS), None)
    row = "  "
    for m in M_RANGE:
        w, t = res[m]
        row += f"m={m}:{w}/{t} "
    print(row)
    if first_win:
        print(f"  → SUCCESS: first 3/3 at m={first_win}")
    else:
        max_w = max(w for w, t in res.values())
        print(f"  → FAIL: never 3/3 (max {max_w}/{N_SEEDS})")
    print()
    ref_summary.append((lam_n, label, first_win))

# ---- Summary table ---------------------------------------------------------
print("=" * 72)
print("THRESHOLD SUMMARY (K1=8, max m tested = 16, 3 seeds)")
print(f"{'λ/n bin':<20} {'lam/n':<8} {'n':>8}  {'first 3/3 m':<14} outcome")
print("-" * 65)

for (bin_lo, lam_n, n, first_win) in summary:
    if lam_n is None:
        print(f"  [{bin_lo:.2f},{bin_lo+BIN_WIDTH:.2f})  {'—':>8}  {'—':>8}  {'NO DATA':<14}")
        continue
    if first_win:
        outcome = f"SUCCESS m={first_win}"
    else:
        outcome = "FAIL (m≤16)"
    print(f"  [{bin_lo:.2f},{bin_lo+BIN_WIDTH:.2f})  {lam_n:.4f}   {n:>8}  {str(first_win) if first_win else 'never':<14} {outcome}")

print()
print("Reference calibration:")
for (lam_n, label, first_win) in ref_summary:
    outcome = f"SUCCESS m={first_win}" if first_win else "FAIL (m≤16)"
    print(f"  lam/n={lam_n:.4f}  {label}")
    print(f"    → {outcome}")

# ---- Threshold interpretation ----------------------------------------------
print()
print("=" * 72)
print("THRESHOLD ANALYSIS:")
successes = [(lam_n, bin_lo) for (bin_lo, lam_n, n, fw) in summary if fw is not None and lam_n is not None]
failures  = [(lam_n, bin_lo) for (bin_lo, lam_n, n, fw) in summary if fw is None and lam_n is not None]

if successes and failures:
    max_fail  = max(lam_n for lam_n, _ in failures)
    min_succ  = min(lam_n for lam_n, _ in successes)
    print(f"  Largest failing  λ/n: {max_fail:.4f}")
    print(f"  Smallest success λ/n: {min_succ:.4f}")
    if max_fail < min_succ:
        print(f"  → Clean boundary: λ/n* ∈ ({max_fail:.4f}, {min_succ:.4f})")
    else:
        print(f"  → Overlapping region: λ/n* is ambiguous or stochastic")
elif not failures:
    print("  All tested bins SUCCEED — threshold is below lowest λ/n tested")
elif not successes:
    print("  All tested bins FAIL — threshold is above highest λ/n tested (or K1/m insufficient)")
print()
print("Done.")
