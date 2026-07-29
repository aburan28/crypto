"""
GLV-HNP Phase 2 — λ/n threshold bisection (Thread 20, 2026-07-29)

Bisects the critical λ/n ratio below which the GLV-aware LLL attack fails.

Prior data (from 2026-07-26 Phase 2 validation):
  λ/n = 0.53  (p=211, n=199, 8-bit):    LLL 3/3 at m=6  ✓
  λ/n = 0.34  (p=524347, n=523969, 20b): LLL 3/3 at m=9  ✓
  λ/n = 0.07  (p=2677, n=2647, 12-bit):  LLL FAIL all m  ✗
              (BKZ-40 also fails — structural, not strength issue)

Goal: find threshold T* such that LLL succeeds iff λ/n > T*.

Approach:
  1. Enumerate j=0 curves y²=x³+b over primes p≤600
  2. For each curve compute n=#E(F_p), check n prime, compute GLV eigenvalue λ
  3. Bucket curves by λ/n into 0.05-wide bins [0.05,0.10), [0.10,0.15), ...
  4. For each bin pick one representative, run LLL sweep m=3..18 with 3 seeds
  5. Record minimum m for 3/3 recovery (or FAIL)

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

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
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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

def count_points_brute(p, b):
    """Count points on y^2 = x^3 + b over F_p (including point at infinity)."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None:
            continue
        elif y == 0:
            count += 1
        else:
            count += 2
    return count

def is_prime_miller(n):
    """Miller-Rabin primality test (deterministic for n < 3,215,031,751)."""
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0: r += 1; d //= 2
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def sieve_primes(limit):
    """Return list of primes up to limit."""
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, limit+1, i):
                sieve[j] = False
    return [i for i in range(2, limit+1) if sieve[i]]

def glv_eigenvalue(n):
    """
    Find the GLV eigenvalue for a j=0 curve with group order n (prime).
    Returns (lam, lam/n) where lam is the smallest root of x^2+x+1 ≡ 0 mod n.
    Returns None if no such root exists (n ≢ 1 mod 3, so D=-3 CM unavailable).
    """
    if n % 3 != 1:
        return None
    # x^2+x+1 ≡ 0 mod n iff -3 is QR mod n
    neg3 = (-3) % n
    if pow(neg3, (n - 1) // 2, n) != 1:
        return None
    sqrt_neg3 = tonelli_shanks(neg3, n)
    if sqrt_neg3 is None:
        return None
    # Two roots: (-1 + sqrt_neg3)/2 and (-1 - sqrt_neg3)/2 mod n
    inv2 = modinv(2, n)
    r1 = (-1 + sqrt_neg3) * inv2 % n
    r2 = (-1 - sqrt_neg3) * inv2 % n
    # Verify
    assert (r1 * r1 + r1 + 1) % n == 0, f"r1={r1} not a root mod {n}"
    assert (r2 * r2 + r2 + 1) % n == 0, f"r2={r2} not a root mod {n}"
    lam = min(r1, r2)
    return lam

def find_generator(p, b, n):
    """Find a generator of the prime-order subgroup of y^2=x^3+b over F_p."""
    rng = random.Random(42)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        G = (x, y)
        # Since n is prime and G is on the curve, G is either order n or 1.
        # #E(F_p) = n means every non-identity point has order n.
        if ec_mul(G, n, p) is None:
            return G
    return None

# ---------------------------------------------------------------------------
# Phase 2 lattice (same structure as glv_hnp_phase2_attack.py)
# ---------------------------------------------------------------------------

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    """
    Build the (2m+2) x (2m+2) GLV-aware HNP lattice.
    sigs: list of (r_i, h_i, s_inv_i) from ECDSA signing.

    Planted short vector: (k1_1*S_K1, ..., k1_m*S_K1, d*S_D,
                           k2_1*S_K2, ..., k2_m*S_K2, S_KANNAN)

    Rows:
      0..m-1:  n*S_K1  (mod-n constraint for each k1_i)
      m:       d-row:  [a_i*S_K1, ..., S_D, 0, ...]  (a_i = r*s_inv)
      m+1..2m: k2 rows: [-lam*S_K1 in col i, S_K2 in col m+1+i]
      2m+1:    Kannan: [b_i*S_K1, ..., 0, S_KANNAN]  (b_i = h*s_inv)

    Planted combination: -q_i*(mod row i) + d*(d-row) + k2_i*(k2 row i) + 1*(Kannan)
    k1-col_i = S_K1*(d*a_i + b_i - lam*k2_i - q_i*n) = S_K1*k1_i  (ECDSA: h*s_inv + d*r*s_inv = k_full)
    """
    m = len(sigs)
    dim = 2 * m + 2

    S_K1     = n // k1_bound if k1_bound > 0 else 1
    S_D      = 1
    S_K2     = n // k2_bound if k2_bound > 0 else 1
    S_KANNAN = n

    rows = [[0] * dim for _ in range(dim)]

    # Rows 0..m-1: mod-n reduction for k1_i
    for i in range(m):
        rows[i][i] = n * S_K1

    # Row m: d-row (a_i = r*s_inv is the d-coefficient)
    rows[m][m] = S_D
    for i in range(m):
        r_i, h_i, s_inv_i = sigs[i]
        a_i = r_i * s_inv_i % n
        rows[m][i] = a_i * S_K1   # no modular reduction: a_i < n, so a_i*S_K1 < n*S_K1 always

    # Rows m+1..2m: k2 rows (b_i = h*s_inv is the constant term)
    for i in range(m):
        rows[m + 1 + i][m + 1 + i] = S_K2      # k2 diagonal (NOT n*S_K2)
        rows[m + 1 + i][i] = -lam * S_K1        # negative, no modular reduction

    # Row 2m+1: Kannan embedding row (b_i constant term)
    rows[2 * m + 1][2 * m + 1] = S_KANNAN
    for i in range(m):
        r_i, h_i, s_inv_i = sigs[i]
        b_i = h_i * s_inv_i % n
        rows[2 * m + 1][i] = b_i * S_K1

    # Convert to fpylll IntegerMatrix (use M[i, j] not M[i][j])
    M = IntegerMatrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            M[i, j] = rows[i][j]

    return M, S_K1, S_D, S_K2, S_KANNAN

def gen_signatures(G, n, lam, p, b_coeff, d_secret, k1_bound, k2_bound, m, seed):
    """Generate m ECDSA signatures with GLV-biased nonces."""
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 100000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0:
            continue
        R = ec_mul(G, k_full, p)
        if R is None:
            continue
        r = R[0] % n
        if r == 0:
            continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0:
            continue
        s_inv = modinv(s, n)
        sigs.append((r, h, s_inv))
    return sigs if len(sigs) == m else None

def recover_d(reduced, n, m):
    """Try to recover d from LLL-reduced basis."""
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) == 0:
            continue
        # Try extracting d from col m (d slot)
        d_cand = row[m] % n
        if d_cand == 0:
            continue
        yield d_cand

def run_once(p, b_coeff, n, lam, G, k1_bound, k2_bound, m, d_secret, seed):
    """Run one attack instance. Returns True if d recovered."""
    sigs = gen_signatures(G, n, lam, p, b_coeff, d_secret, k1_bound, k2_bound, m, seed)
    if sigs is None:
        return False

    M, S_K1, S_D, S_K2, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    LLL.reduction(M)

    dim = 2 * m + 2
    for d_cand in recover_d([[M[i][j] for j in range(dim)] for i in range(dim)], n, m):
        if d_cand == d_secret:
            return True
        # Also check negation
        if (n - d_cand) % n == d_secret:
            return True
    return False

def sweep_m(label, p, b_coeff, n, lam, G, k1_bound, k2_bound,
            d_secret, seeds, m_range, verbose=True):
    """Sweep over m values. Return first m where all seeds succeed, or None."""
    if verbose:
        print(f"\n  [{label}]  n={n}, λ={lam}, λ/n={lam/n:.4f}, "
              f"K1={k1_bound}, K2={k2_bound}")
    results = {}
    first_success = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            rng_d = random.Random(seed + 777)
            d = rng_d.randint(1, n - 1)
            ok = run_once(p, b_coeff, n, lam, G, k1_bound, k2_bound, m, d, seed)
            wins += ok
        results[m] = wins
        total = len(seeds)
        if wins == total and first_success is None:
            first_success = m
        if verbose:
            status = "✓" if wins == total else ("✗" if wins == 0 else f"{wins}/{total}")
            print(f"    m={m:2d}: {status}")
        # Early-stop: once 3/3 confirmed, we're done (can still check a few more)
        if wins == total and m >= min(m_range) + 2:
            break
    return first_success, results

# ---------------------------------------------------------------------------
# Step 1: Enumerate j=0 curves over small primes, collect (n, λ, λ/n)
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n Threshold Bisection (Thread 20)")
print("=" * 70)
print()
print("Step 1: Enumerating j=0 curves y²=x³+b over primes p≤600 ...")

PRIMES = sieve_primes(600)
candidates = []   # (λ/n, p, b, n, λ)

for p in PRIMES:
    if p < 7:
        continue  # too small
    for b in range(1, min(p, 30)):
        n = count_points_brute(p, b)
        if not is_prime_miller(n):
            continue
        if n < 30:
            continue  # too small for meaningful lattice
        lam = glv_eigenvalue(n)
        if lam is None:
            continue
        ratio = lam / n
        candidates.append((ratio, p, b, n, lam))

candidates.sort()
print(f"  Found {len(candidates)} candidate (p, b, n, λ) tuples with prime n and GLV eigenvalue.")
print()

# Show distribution of λ/n
BINS = [(0.05*i, 0.05*(i+1)) for i in range(1, 10)]  # [0.05,0.10), ..., [0.45,0.50)
print("  λ/n distribution:")
for lo, hi in BINS:
    in_bin = [c for c in candidates if lo <= c[0] < hi]
    print(f"    [{lo:.2f}, {hi:.2f}): {len(in_bin)} curves")
print()

# ---------------------------------------------------------------------------
# Step 2: Pick one representative per bin
# ---------------------------------------------------------------------------

print("Step 2: Selecting one representative per bin ...")
SEEDS = [42, 1234, 9999]
EFF_TARGET = 0.15   # target efficiency eff = K1*K2/n (same as 8-bit baseline)
M_MAX = 18

# Add the two known reference curves
REFERENCE_CURVES = [
    # (λ/n, p, b_coeff, n, λ, label)
    (0.073, 2677, 2, 2647, 185, "known-fail  (p=2677, λ/n=0.07)"),
    (0.532, 211,  2,  199, 106, "known-pass  (p=211,  λ/n=0.53)"),
]

selected = {}  # bin_key -> (ratio, p, b, n, lam)
for lo, hi in BINS:
    in_bin = [c for c in candidates if lo <= c[0] < hi]
    if in_bin:
        # Prefer larger n for stronger signal
        best = max(in_bin, key=lambda c: c[3])
        selected[(lo, hi)] = best

print(f"  Selected {len(selected)} representative curves:")
for (lo, hi), (ratio, p, b, n, lam) in sorted(selected.items()):
    print(f"    [{lo:.2f}, {hi:.2f}): p={p}, b={b}, n={n}, λ={lam}, λ/n={ratio:.4f}")
print()

# ---------------------------------------------------------------------------
# Step 3: Run the attack on each representative
# ---------------------------------------------------------------------------

print("Step 3: Running LLL attack on each representative (3 seeds, m=3..18) ...")
print("  Target efficiency eff=K1*K2/n ≈ 0.15")
print()

all_results = []   # (λ/n, label, first_success, results_dict)

# First run reference curves
for (ref_ratio, ref_p, ref_b, ref_n, ref_lam, ref_label) in REFERENCE_CURVES:
    k2 = math.isqrt(ref_n) + 1
    k1 = max(2, int(EFF_TARGET * ref_n / k2))
    actual_eff = k1 * k2 / ref_n
    print(f"  Reference: {ref_label}")
    print(f"    K1={k1}, K2={k2}, eff={actual_eff:.3f}")
    G = find_generator(ref_p, ref_b, ref_n)
    if G is None:
        print(f"    ERROR: no generator found for p={ref_p}, b={ref_b}, n={ref_n}")
        all_results.append((ref_ratio, ref_label, None, {}))
        continue
    m_range = range(3, M_MAX + 1)
    first_ok, res = sweep_m(ref_label, ref_p, ref_b, ref_n, ref_lam, G,
                            k1, k2, 0, SEEDS, m_range, verbose=True)
    all_results.append((ref_ratio, ref_label, first_ok, res))

# Then run each selected bin
for (lo, hi), (ratio, p, b, n, lam) in sorted(selected.items()):
    label = f"bin[{lo:.2f},{hi:.2f}) p={p} n={n} λ/n={ratio:.3f}"
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(EFF_TARGET * n / k2))
    actual_eff = k1 * k2 / n
    print(f"\n  {label}")
    print(f"    K1={k1}, K2={k2}, eff={actual_eff:.3f}")
    G = find_generator(p, b, n)
    if G is None:
        print(f"    ERROR: no generator found")
        all_results.append((ratio, label, None, {}))
        continue
    m_range = range(3, M_MAX + 1)
    first_ok, res = sweep_m(label, p, b, n, lam, G, k1, k2, 0,
                            SEEDS, m_range, verbose=True)
    all_results.append((ratio, label, first_ok, res))

# ---------------------------------------------------------------------------
# Step 4: Summary table
# ---------------------------------------------------------------------------

print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'λ/n':>8}  {'n':>8}  {'K1':>4}  {'K2':>5}  {'LLL result':>20}")
print("-" * 70)

all_results.sort(key=lambda x: x[0])
threshold_lo = None
threshold_hi = None

for (ratio, label, first_ok, res) in all_results:
    _, p_s, b_s, n_s, lam_s = next(
        ((r, p, b, n, lam) for r, p, b, n, lam in candidates if abs(r-ratio) < 1e-6),
        (None, None, None, 0, 0)
    )
    # Get from reference if not in candidates
    for ref_ratio, ref_p, ref_b, ref_n, ref_lam, _ in REFERENCE_CURVES:
        if abs(ratio - ref_ratio) < 0.01:
            n_s = ref_n
            lam_s = ref_lam
            ref_k2 = math.isqrt(ref_n) + 1
            ref_k1 = max(2, int(EFF_TARGET * ref_n / ref_k2))
            p_s, b_s = ref_p, ref_b
    k2_s = math.isqrt(n_s) + 1 if n_s else 0
    k1_s = max(2, int(EFF_TARGET * n_s / k2_s)) if n_s else 0

    if first_ok is not None:
        lll_str = f"3/3 at m={first_ok:2d}"
        if threshold_hi is None or ratio > threshold_hi:
            threshold_hi = ratio
    else:
        lll_str = f"FAIL (all m≤{M_MAX})"
        if threshold_lo is None or ratio > threshold_lo:
            threshold_lo = ratio

    print(f"{ratio:>8.4f}  {n_s:>8d}  {k1_s:>4d}  {k2_s:>5d}  {lll_str:>20}  | {label}")

print()
print("=" * 70)
print("FINDINGS")
print("=" * 70)
if threshold_lo is not None and threshold_hi is not None:
    print(f"  Failure region:  λ/n ≤ {threshold_lo:.4f}")
    print(f"  Success region:  λ/n ≥ {threshold_hi:.4f}")
    print(f"  Threshold T* lies in ({threshold_lo:.4f}, {threshold_hi:.4f})")
elif threshold_lo is None:
    print(f"  All tested curves SUCCEED at λ/n ≥ {threshold_hi:.4f}")
    print(f"  Lower threshold not found in this sweep")
elif threshold_hi is None:
    print(f"  All tested curves FAIL at λ/n ≤ {threshold_lo:.4f}")
    print(f"  Upper threshold not found in this sweep")
else:
    print("  Inconclusive — need more data points")

print()
print("INTERPRETATION")
print("-" * 70)
print("  When λ/n is small, the λ-rows in the (2m+2)×(2m+2) lattice carry")
print("  small entries (-λ·S_K1 ≈ -λ·n/K1 << n²/K1). After LLL, the basis")
print("  mixes λ-rows with modular rows, losing the k2 separation signal.")
print("  The threshold T* quantifies how large λ must be relative to n for")
print("  LLL to cleanly isolate the short planted vector.")
print()
print("Done.")
