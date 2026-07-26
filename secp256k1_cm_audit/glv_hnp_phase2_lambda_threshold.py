"""
GLV-HNP Phase 2 — Lambda/n threshold study.

Bisects the lambda/n threshold between LLL failure and success
for the GLV-aware Phase 2 attack on 12-bit j=0 prime-order curves.

Known data points:
  lambda/n = 0.07 (p=2677, n=2647): FAIL at all m <= 12 (even BKZ-40)
  lambda/n = 0.34 (p=524347, n=523969, 20-bit): SUCCESS at m=9
  lambda/n = 0.53 (p=211, n=199, 8-bit): SUCCESS at m=6
  lambda/n = 0.66 (p=2557, n=2659, 12-bit): SUCCESS at m=7

Strategy:
  1. Enumerate all 12-bit primes p with p ≡ 1 (mod 3).
  2. For each p, compute all 6 j=0 curve orders n via Eisenstein decomposition.
  3. Keep (p, n) where n is prime and n ≡ 1 (mod 3).
  4. Compute lambda = min root of x^2+x+1 = 0 mod n.
  5. Bin by lambda/n in intervals of 0.05 (bins: [0.05,0.10), ..., [0.45,0.50)).
  6. For each bin, pick one representative curve.
  7. Sweep m = 2..18 with 3 seeds each; record min m for 3/3 success.

Output: table of lambda/n vs min_m_for_3of3_success (or FAIL).

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass a=0)
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
    m_iter, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_iter - i - 1), p)
        m_iter, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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
    rng = random.Random(12345)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def find_twist_b(p, n):
    """Find smallest b such that #(y^2=x^3+b over F_p) = n."""
    rng = random.Random(42)
    for b in range(1, p):
        count = 0
        ok = False
        for _ in range(200):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    ok = True
                    break
        if ok:
            return b
    return None

# ---------------------------------------------------------------------------
# CM / Eisenstein theory
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
    if n % 3 != 1: return None, None
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

# ---------------------------------------------------------------------------
# Curve enumeration: find 12-bit j=0 prime-order GLV curves, bin by lambda/n
# ---------------------------------------------------------------------------

def enumerate_12bit_glv_curves():
    """
    Enumerate all 12-bit primes p ≡ 1 (mod 3), find all j=0 prime-order curves,
    compute lambda/n, return list of (p, n, lam, lam_over_n).
    """
    curves = []
    p = 2053  # first prime > 2048 with p ≡ 1 mod 3
    while p < 4096:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 3: continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam, lam2 = glv_eigenvalue(n_cand)
                        if lam is None: continue
                        ratio = lam / n_cand
                        curves.append((p, n_cand, lam, ratio))
        p = sympy.nextprime(p)
    return curves

def bin_curves_by_ratio(curves, bin_width=0.05):
    """Group curves by lambda/n in bins [k*bw, (k+1)*bw), k=1..9 (0.05..0.50)."""
    bins = {}
    for (p, n, lam, ratio) in curves:
        if ratio <= 0.0 or ratio >= 0.50: continue
        bin_k = int(ratio / bin_width)
        if bin_k not in bins:
            bins[bin_k] = []
        bins[bin_k].append((p, n, lam, ratio))
    return bins

# ---------------------------------------------------------------------------
# Attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return True
    return False

def run_attack(curve_params, m, k1_bound, seeds):
    """Run the attack on curve_params with m signatures and given seeds.
    Returns number of successes."""
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    wins = 0
    for seed in seeds:
        d_secret = random.Random(seed + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
        if len(sigs) < m:
            continue
        M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        A = IntegerMatrix.from_matrix(M)
        LLL.reduction(A)
        reduced = [[A[i][j] for j in range(2*m+2)] for i in range(2*m+2)]
        if recover_d(reduced, m, n, S_KANNAN, d_secret):
            wins += 1
    return wins

def sweep_m_until_success(curve_params, k1_bound, seeds, m_range, target_wins):
    """Sweep m values until target_wins/len(seeds) success, or return None."""
    for m in m_range:
        wins = run_attack(curve_params, m, k1_bound, seeds)
        if wins >= target_wins:
            return m, wins
    return None, 0

# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — Lambda/n threshold study")
print("=" * 70)
print()

# Step 1: Enumerate 12-bit curves
print("Step 1: Enumerating 12-bit j=0 prime-order GLV curves...")
curves = enumerate_12bit_glv_curves()
print(f"  Found {len(curves)} (p, n, lambda) triples.")
print()

# Step 2: Bin by lambda/n
BIN_WIDTH = 0.05
bins = bin_curves_by_ratio(curves, BIN_WIDTH)
print("Step 2: Binning by lambda/n:")
for k in sorted(bins):
    lo, hi = k * BIN_WIDTH, (k + 1) * BIN_WIDTH
    print(f"  [{lo:.2f}, {hi:.2f}): {len(bins[k])} curves  "
          f"(e.g. p={bins[k][0][0]}, n={bins[k][0][1]}, lam/n={bins[k][0][3]:.3f})")
print()

# Step 3: Pick representatives and run attack
print("Step 3: Attack sweep for each lambda/n bin")
print("-" * 70)

SEEDS = [42, 1234, 9999]
M_RANGE = range(2, 20)
TARGET_WINS = 3  # need 3/3 for "success"

# Fixed k1_bound for all 12-bit curves: eff ~ K1*K2/n ~ K1*sqrt(n)/n = K1/sqrt(n)
# For 12-bit n ≈ 2048-4096, sqrt(n) ≈ 50-64.
# K1=8 gives eff ≈ 8/55 ≈ 0.15 (same as prior work on p=2677 and p=2557)
K1_BOUND = 8

results_table = []
failed_to_build = []

for k in sorted(bins):
    lo, hi = k * BIN_WIDTH, (k + 1) * BIN_WIDTH
    # Try first a few candidates in the bin until we find one we can build
    found_curve = None
    for (p, n, lam, ratio) in bins[k][:5]:
        b = find_twist_b(p, n)
        if b is None:
            continue
        G = find_generator(p, b, n)
        if G is None:
            continue
        found_curve = (p, b, n, lam, G, ratio)
        break

    if found_curve is None:
        failed_to_build.append((lo, hi))
        print(f"  [{lo:.2f}, {hi:.2f}): SKIP — could not build curve")
        continue

    p, b, n, lam, G, ratio = found_curve
    curve_params = (p, b, n, lam, G)
    k2_bound = math.isqrt(n) + 1
    eff = K1_BOUND * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else 999

    print(f"  [{lo:.2f}, {hi:.2f}): p={p}, n={n}, lam={lam}, lam/n={ratio:.3f}")
    print(f"    K1={K1_BOUND}, K2={k2_bound}, eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")

    m_success, wins_success = sweep_m_until_success(
        curve_params, K1_BOUND, SEEDS, M_RANGE, TARGET_WINS
    )

    if m_success is not None:
        print(f"    LLL SUCCESS at m={m_success} ({wins_success}/{len(SEEDS)} seeds) ✓")
        results_table.append((lo, hi, ratio, p, n, lam, m_success, "SUCCESS"))
    else:
        # Count wins at max m
        m_max = max(M_RANGE)
        wins_at_max = run_attack(curve_params, m_max, K1_BOUND, SEEDS)
        print(f"    LLL FAIL: 0/3 at m={m_max}. wins_at_max={wins_at_max}")
        results_table.append((lo, hi, ratio, p, n, lam, None, f"FAIL (0/3 at m={m_max})"))
    print()

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print("=" * 70)
print("RESULTS: Lambda/n threshold for LLL success in GLV Phase 2 attack")
print(f"  (12-bit curves, K1={K1_BOUND}, K2=sqrt(n)+1, seeds={SEEDS})")
print("=" * 70)
print()
print(f"  {'lambda/n bin':<18} {'repr lam/n':>10} {'n':>6} {'min m (3/3)':>12} {'outcome'}")
print(f"  {'-'*18:<18} {'-'*10:>10} {'-'*6:>6} {'-'*12:>12} {'-'*20}")
for (lo, hi, ratio, p, n, lam, m_succ, outcome) in results_table:
    bin_str = f"[{lo:.2f}, {hi:.2f})"
    m_str = str(m_succ) if m_succ is not None else "—"
    print(f"  {bin_str:<18} {ratio:>10.3f} {n:>6} {m_str:>12}   {outcome}")
print()

# Report threshold
success_ratios = [(ratio, m_succ) for (lo, hi, ratio, p, n, lam, m_succ, outcome) in results_table
                  if outcome == "SUCCESS"]
fail_ratios    = [(ratio,) for (lo, hi, ratio, p, n, lam, m_succ, outcome) in results_table
                  if outcome != "SUCCESS"]

if success_ratios and fail_ratios:
    min_success_ratio = min(r for (r, _) in success_ratios)
    max_fail_ratio    = max(r for (r,) in fail_ratios)
    print(f"THRESHOLD ESTIMATE:")
    print(f"  Highest failure lam/n:  {max_fail_ratio:.3f}")
    print(f"  Lowest success lam/n:   {min_success_ratio:.3f}")
    print(f"  Threshold bracket:      ({max_fail_ratio:.3f}, {min_success_ratio:.3f})")
elif success_ratios:
    min_s = min(r for (r, _) in success_ratios)
    print(f"  All tested bins SUCCEED (min lam/n tested = {min_s:.3f})")
    print(f"  Threshold is below {min_s:.3f} — need smaller-lam/n curves to find it.")
else:
    print(f"  All tested bins FAIL — unexpected (check setup).")

# Known data points for comparison
print()
print("Known data points from prior runs:")
print("  lam/n=0.07 (p=2677, 12-bit): FAIL for all m <= 12 (incl. BKZ-40)")
print("  lam/n=0.34 (p=524347, 20-bit): SUCCESS at m=9")
print("  lam/n=0.53 (p=211, 8-bit): SUCCESS at m=6")
print("  lam/n=0.66 (p=2557, 12-bit): SUCCESS at m=7")
print()
print("Done.")
