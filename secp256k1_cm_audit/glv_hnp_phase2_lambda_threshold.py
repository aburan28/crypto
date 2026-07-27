"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Bisects the LLL success threshold between λ/n = 0.07 (known fail)
and λ/n = 0.34 (known success).  Uses 12–14-bit j=0 CM curves
with fixed efficiency eff ≈ 0.10 to isolate the geometric effect
of λ/n on the lattice.

For each of 9 bins covering [0.05, 0.50] with width 0.05, we:
  1. Find 1–2 representative j=0 curves with λ/n in that bin.
  2. Sweep m ∈ {5, 6, 7, 8, 9, 10} at 3 independent seeds.
  3. Report the minimum m at which 3/3 seeds succeed (or "never").

Invariant: K1_BOUND chosen so that eff = K1*K2/n ≈ 0.10 for all
curves.  K2_BOUND = ⌊√n⌋ + 1.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0, reused from phase2 scripts)
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
# CM theory: Eisenstein decomposition + trace enumeration
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - a*b + b^2 = p, a,b >= 0."""
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
    """Six Frobenius traces for the six sextic twists of j=0."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def glv_eigenvalue_min(n):
    """
    Compute the smaller GLV eigenvalue: min root of x^2+x+1≡0 mod n.
    Returns (lam_min, lam_max) with lam_min < n//2.
    """
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    # Verify both roots
    if (r1 * r1 + r1 + 1) % n != 0 or (r2 * r2 + r2 + 1) % n != 0:
        return None, None
    lam_min = min(r1, r2)
    lam_max = max(r1, r2)
    return lam_min, lam_max

# ---------------------------------------------------------------------------
# Curve search: populate bins by λ/n ratio
# ---------------------------------------------------------------------------

def find_curves_by_lambda_ratio(target_bins, min_bits=11, max_bits=14):
    """
    Search for j=0 curves whose λ_min/n falls into each bin.
    target_bins: list of (lo, hi) float pairs.
    Returns dict: bin_idx -> list of (p, b, n, lam, G) tuples.
    """
    results = {i: [] for i in range(len(target_bins))}
    needed = {i: 2 for i in range(len(target_bins))}  # want 2 per bin

    p = 2 ** min_bits + 1
    if p % 2 == 0: p += 1

    checked = 0
    while any(n > 0 for n in needed.values()):
        # Advance to next odd number in range
        if p >= 2 ** max_bits:
            # Wrap: we've exhausted the range
            print(f"  [search] exhausted [{2**min_bits}, {2**max_bits}]; "
                  f"still need: "
                  f"{[(i, needed[i]) for i in range(len(target_bins)) if needed[i] > 0]}")
            break

        if p % 2 == 0:
            p += 1
            continue

        if not is_prime(p) or p % 3 != 1:
            p += 2
            continue

        eis = eisenstein_decompose(p)
        if eis is None:
            p += 2
            continue

        a_e, b_e = eis
        traces = j0_traces(a_e, b_e)
        for t in traces:
            n_cand = p + 1 - t
            if n_cand < 2 ** (min_bits - 1):
                continue
            if not is_prime(n_cand) or n_cand % 3 != 1:
                continue

            lam_min, lam_max = glv_eigenvalue_min(n_cand)
            if lam_min is None:
                continue

            ratio_min = lam_min / n_cand
            # Check each bin
            for i, (lo, hi) in enumerate(target_bins):
                if needed[i] <= 0:
                    continue
                if lo <= ratio_min < hi:
                    # Find b parameter
                    found_b = None
                    for b_try in range(1, 50):
                        # Quick order check via random point
                        rng_tmp = random.Random(b_try * 999 + 7)
                        for _ in range(20):
                            x = rng_tmp.randint(0, p - 1)
                            rhs = (pow(x, 3, p) + b_try) % p
                            y = tonelli_shanks(rhs, p)
                            if y is not None and y != 0:
                                Q = ec_mul((x, y), n_cand, p)
                                if Q is None:
                                    found_b = b_try
                                    break
                        if found_b is not None:
                            break

                    if found_b is None:
                        continue

                    G = find_generator(p, found_b, n_cand)
                    if G is None:
                        continue

                    results[i].append((p, found_b, n_cand, lam_min, G))
                    needed[i] -= 1
                    break  # one bin per (p, trace) pair

        checked += 1
        p += 2

    return results

# ---------------------------------------------------------------------------
# Lattice construction and LLL attack (from phase2_20bit.py)
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

def run_lll_attack(curve_params, m, d_secret, k1_bound, seed=42):
    p, b, n, lam, G = curve_params
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

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

# Bin boundaries for λ_min/n
BINS = [
    (0.04, 0.10),   # bin 0: [0.04, 0.10) — contains known fail at 0.07
    (0.10, 0.16),   # bin 1
    (0.16, 0.22),   # bin 2
    (0.22, 0.28),   # bin 3
    (0.28, 0.34),   # bin 4
    (0.34, 0.40),   # bin 5 — contains known success at 0.34
    (0.40, 0.46),   # bin 6
    (0.46, 0.50),   # bin 7
]

BIN_LABELS = [f"[{lo:.2f},{hi:.2f})" for (lo, hi) in BINS]

SEEDS = [42, 1234, 9999]
M_RANGE = list(range(5, 11))    # m ∈ {5,...,10}
EFF_TARGET = 0.10               # fixed efficiency

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection")
print("  Bins: 8 bins from [0.04, 0.50), width ~0.06")
print(f"  eff_target ≈ {EFF_TARGET:.2f}  |  m ∈ {M_RANGE}  |  seeds: {SEEDS}")
print("=" * 70)

print("\nSearching for representative curves (bit range 11–14)...")
bin_curves = find_curves_by_lambda_ratio(BINS, min_bits=11, max_bits=14)

# Report what we found
print("\nFound curves per bin:")
for i, (lo, hi) in enumerate(BINS):
    curves = bin_curves[i]
    if not curves:
        print(f"  bin {i} [{lo:.2f},{hi:.2f}): NONE FOUND")
    else:
        for (p, b, n, lam, G) in curves:
            print(f"  bin {i} [{lo:.2f},{hi:.2f}): p={p}, n={n} ({n.bit_length()}b), "
                  f"lam={lam}, lam/n={lam/n:.4f}")

# Attack sweep
print("\n" + "=" * 70)
print("LLL attack sweep (3/3 seeds threshold)")
print("=" * 70)

results_summary = []  # (bin_label, lam_n_ratio, first_m_3of3, curve_label)

for i, (lo, hi) in enumerate(BINS):
    curves = bin_curves[i]
    if not curves:
        print(f"\nbin {i} {BIN_LABELS[i]}: NO CURVE — skipping")
        results_summary.append((BIN_LABELS[i], (lo+hi)/2, None, "no curve"))
        continue

    # Take first (most representative) curve
    p, b, n, lam, G = curves[0]
    ratio = lam / n
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(EFF_TARGET * n / k2_bound))
    eff = k1_bound * k2_bound / n

    print(f"\nbin {i} {BIN_LABELS[i]}: p={p}, n={n} ({n.bit_length()}b), "
          f"lam={lam}, lam/n={ratio:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}")

    d_rng = random.Random(77777)
    d_secret = d_rng.randint(1, n - 1)

    first_m_3of3 = None
    for m in M_RANGE:
        wins = 0
        for seed in SEEDS:
            ok = run_lll_attack((p, b, n, lam, G), m, d_secret, k1_bound, seed)
            wins += ok
        status = f"{wins}/{len(SEEDS)}"
        marker = " ← 3/3 ✓" if wins == len(SEEDS) else ""
        print(f"  m={m}: {status}{marker}")
        if wins == len(SEEDS) and first_m_3of3 is None:
            first_m_3of3 = m
            break  # stop once we hit 3/3

    if first_m_3of3 is None:
        print(f"  → never 3/3 in m ∈ {M_RANGE}")

    results_summary.append((BIN_LABELS[i], ratio, first_m_3of3,
                             f"n={n}({n.bit_length()}b),lam={lam}"))

# Final summary table
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print(f"{'bin λ/n range':<18} {'lam/n':<8} {'1st m (3/3)':<13} {'curve'}")
print("-" * 70)
for (blabel, ratio, first_m, clabel) in results_summary:
    if first_m is None:
        m_str = "never (>10)"
    else:
        m_str = str(first_m)
    print(f"  {blabel:<16} {ratio:<8.4f} {m_str:<13} {clabel}")

# Identify threshold
print("\n" + "=" * 70)
print("THRESHOLD ANALYSIS")
passing = [(ratio, first_m) for (_, ratio, first_m, _) in results_summary if first_m is not None]
failing = [(ratio, label) for (label, ratio, first_m, _) in results_summary if first_m is None]
if passing and failing:
    max_fail_ratio = max(r for (r, _) in failing)
    min_pass_ratio = min(r for (r, _) in passing)
    print(f"  Highest λ/n that FAILS (never 3/3): {max_fail_ratio:.4f}")
    print(f"  Lowest  λ/n that PASSES (3/3):      {min_pass_ratio:.4f}")
    print(f"  → LLL threshold for 3/3 recovery: λ/n ∈ ({max_fail_ratio:.2f}, {min_pass_ratio:.2f})")
elif not failing:
    print("  All bins pass — consider testing smaller λ/n values.")
elif not passing:
    print("  All bins fail — consider testing larger λ/n values.")

print("\nDone.")
