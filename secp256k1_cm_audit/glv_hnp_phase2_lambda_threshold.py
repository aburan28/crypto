"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Goal: bisect the LLL-success threshold between λ/n=0.07 (known fail)
and λ/n=0.34 (known pass).

Method:
  1. Search for j=0 prime-order curves in the 12-15 bit range with
     λ/n values spanning [0.05, 0.50] in ~0.05 increments.
  2. For each λ/n bin, run the GLV-HNP Phase 2 lattice attack
     (column-scaled, Kannan embedding) with m=8, 5 seeds.
  3. Record 0-5 success counts per bin.
  4. Report the empirical threshold.

From the 2026-07-26 session:
  Known FAIL: p=2677, n=2647, λ/n=0.07  (LLL 0/5, BKZ-40 0/5)
  Known PASS: p=524347, n=523969, λ/n=0.34  (LLL 3/3 at m=9)

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0, y^2 = x^3 + b over F_p)
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
    mc, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mc - i - 1), p)
        mc, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition for j=0 curves over F_p
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
                b_val = num // 2
                if b_val >= 0 and a * a - a * b_val + b_val * b_val == p:
                    return (a, b_val)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

def find_b_for_order(p, n):
    """Find curve parameter b such that #(y^2=x^3+b/F_p) == n.
    Checks b=1..12 (covers all 6 twist classes twice for safety)."""
    for b in range(1, 13):
        rng = random.Random(b * 7)
        for _ in range(200):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b
                break
    return None

# ---------------------------------------------------------------------------
# Curve database: find j=0 prime-order curves for each λ/n bin
# ---------------------------------------------------------------------------

TARGET_BINS = [
    (0.04, 0.08),   # includes known-fail λ/n=0.07
    (0.08, 0.12),
    (0.12, 0.16),
    (0.16, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.36),   # includes known-pass λ/n=0.34
    (0.36, 0.42),
    (0.42, 0.50),   # secp256k1: λ/n ≈ 0.44
]

# Fine bisection bins: narrow the [0.04,0.12] range
FINE_BINS = [
    (0.040, 0.060),
    (0.060, 0.072),  # includes known-fail 0.07
    (0.072, 0.082),
    (0.082, 0.092),
    (0.092, 0.100),  # includes new-pass 0.094
    (0.100, 0.115),
    (0.115, 0.130),
]

def find_curves_by_lambda_ratio(bit_lo=12, bit_hi=14, max_per_bin=1):
    """
    Search for j=0 prime-order curves in [2^bit_lo, 2^bit_hi] and
    categorize by λ/n ratio. Returns a dict: bin_key -> list of
    (p, b, n, lam, lam_ratio).
    """
    bins = {b: [] for b in range(len(TARGET_BINS))}
    found_count = {b: 0 for b in range(len(TARGET_BINS))}

    p_lo, p_hi = 2**bit_lo, 2**bit_hi
    p = sympy.nextprime(p_lo - 1)
    checked = 0

    print(f"Searching j=0 curves in [{p_lo}, {p_hi}] (bits {bit_lo}-{bit_hi})...")

    while p < p_hi:
        checked += 1
        if p % 6 != 1:
            p = sympy.nextprime(p)
            continue

        eis = eisenstein_decompose(p)
        if eis is None:
            p = sympy.nextprime(p)
            continue
        a_eis, b_eis = eis

        traces = j0_traces(a_eis, b_eis)
        for t in traces:
            n = p + 1 - t
            if n <= 1:
                continue
            if not sympy.isprime(n):
                continue
            if n % 3 != 1:  # need GLV eigenvalue
                continue

            lam, lam2 = glv_eigenvalue(n)
            if lam is None:
                continue

            lam_ratio = lam / n
            # Also consider the other eigenvalue ratio
            lam2_ratio = lam2 / n

            # Use the min ratio (smaller eigenvalue) for threshold study
            r = min(lam_ratio, lam2_ratio)
            lam_use = lam if lam_ratio == r else lam2

            for bin_idx, (r_lo, r_hi) in enumerate(TARGET_BINS):
                if found_count[bin_idx] >= max_per_bin:
                    continue
                if r_lo <= r < r_hi:
                    # Find the actual curve coefficient b
                    b_curve = find_b_for_order(p, n)
                    if b_curve is None:
                        continue
                    bins[bin_idx].append((p, b_curve, n, lam_use, r))
                    found_count[bin_idx] += 1
                    print(f"  Bin [{r_lo:.2f},{r_hi:.2f}): p={p}, n={n}, λ={lam_use}, λ/n={r:.4f}")
                    break

        if all(found_count[b] >= max_per_bin for b in range(len(TARGET_BINS))):
            break
        p = sympy.nextprime(p)

    print(f"  Searched {checked} primes.")
    return bins

# ---------------------------------------------------------------------------
# GLV-HNP Phase 2 lattice attack (column-scaled)
# ---------------------------------------------------------------------------

def gen_signatures_curve(p, b_coeff, n, lam, G, d_secret, m, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 50000:
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
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def build_glv_lattice_scaled(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2

    S_K1     = n // k1_bound
    S_D      = 1
    S_K2     = n // k2_bound
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]

    # Rows 0..m-1: mod-n constraints (k1_i columns scaled by S_K1)
    for i in range(m):
        M[i][i] = n * S_K1

    # Row m: d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D

    # Rows m+1..2m: k2_i rows
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2

    # Row 2m+1: Kannan target
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN

    return M, S_K1, S_D, S_K2, S_KANNAN

def run_attack(p, b_coeff, n, lam, k1_bound, k2_bound, m_sigs, seed=42):
    """
    Run the GLV-HNP Phase 2 attack. Returns True if d is recovered.
    """
    rng = random.Random(seed + 99999)
    d_secret = rng.randint(1, n - 1)

    G = find_generator(p, b_coeff, n)
    if G is None:
        return False

    sigs = gen_signatures_curve(p, b_coeff, n, lam, G, d_secret, m_sigs,
                                k1_bound, k2_bound, seed=seed)
    if len(sigs) < m_sigs:
        return False

    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice_scaled(
        sigs, n, lam, k1_bound, k2_bound)

    dim = 2 * m_sigs + 2
    A_mat = IntegerMatrix.from_matrix(M)
    LLL.reduction(A_mat)

    # Scan reduced rows for d: look for rows where last slot == ±S_KANNAN
    for i in range(dim):
        row = [A_mat[i][j] for j in range(dim)]
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m_sigs]) % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return True
        # Also check negated (sign flip)
        d_cand2 = (-sign * row[m_sigs]) % n
        if d_cand2 == d_secret:
            return True

    return False

# ---------------------------------------------------------------------------
# Main: find curves, run attack for each bin, report threshold
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 7, 314159]
# m chosen to stay comfortably above the info-theoretic threshold for eff≈0.15:
# m_thresh = ceil(log(n)/log(1/eff)) ≈ 11-12 bits / log2(1/0.15) ≈ 12/2.74 ≈ 4.4 → m=5
# But LLL needs more slack: use m=12 to match the 12-bit pass experiments.
M_SIGS = 12   # number of signatures; above info-theoretic threshold for eff≈0.15

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection study")
print(f"Signatures per trial: m={M_SIGS}  Seeds: {len(SEEDS)}")
print("=" * 70)
print()

# Seed the curve search with a known-fail (λ/n=0.07) and known-pass (λ/n=0.34)
# These are fixed from the 2026-07-26 session; confirm they reproduce.
KNOWN = [
    {"p": 2677, "b": None, "n": 2647, "lam_target": 0.07,
     "label": "known-fail (2026-07-26)"},
    # 20-bit known-pass too large for quick run; use 12-bit proxy
]

bins = find_curves_by_lambda_ratio(bit_lo=11, bit_hi=16, max_per_bin=1)
print()

# Show curve inventory
print("Curve inventory by λ/n bin:")
print(f"{'Bin':20s}  {'p':>8}  {'n':>8}  {'λ':>8}  {'λ/n':>6}")
print("-" * 60)
for bin_idx, (r_lo, r_hi) in enumerate(TARGET_BINS):
    curves = bins[bin_idx]
    if curves:
        p, b, n, lam, r = curves[0]
        print(f"[{r_lo:.2f}, {r_hi:.2f})  {p:>8}  {n:>8}  {lam:>8}  {r:.4f}")
    else:
        print(f"[{r_lo:.2f}, {r_hi:.2f})  {'(not found)':>30}")
print()

# Choose k bounds calibrated to eff ≈ 0.15, matching the 12-bit passing experiments.
# k1_bound = max(2, int(0.15 * sqrt(n))), k2_bound = isqrt(n) + 1
# This gives eff = k1*k2/n ≈ 0.15*sqrt(n)*sqrt(n)/n = 0.15 independent of n.
def make_bounds(n):
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(0.15 * math.sqrt(n)))
    return k1, k2

# Run attack for each bin
print(f"Attack sweep (m={M_SIGS}, {len(SEEDS)} seeds per curve):")
print(f"{'Bin':20s}  {'p':>8}  {'λ/n':>6}  {'LLL success':>12}")
print("-" * 60)

results = []
for bin_idx, (r_lo, r_hi) in enumerate(TARGET_BINS):
    curves = bins[bin_idx]
    if not curves:
        print(f"[{r_lo:.2f}, {r_hi:.2f})  {'--':>8}  {'--':>6}  {'N/A (no curve)':>12}")
        results.append((r_lo, r_hi, None, None, None))
        continue

    p, b, n, lam, r = curves[0]
    k1_bound, k2_bound = make_bounds(n)

    wins = 0
    for seed in SEEDS:
        ok = run_attack(p, b, n, lam, k1_bound, k2_bound, M_SIGS, seed=seed)
        if ok:
            wins += 1

    outcome = f"{wins}/{len(SEEDS)}"
    flag = "PASS" if wins >= 3 else ("MARGINAL" if wins >= 1 else "FAIL")
    print(f"[{r_lo:.2f}, {r_hi:.2f})  {p:>8}  {r:.4f}  {outcome:>5}  ({flag})")
    results.append((r_lo, r_hi, p, r, wins))

print()

# Bisection summary
print("=" * 70)
print("Threshold summary:")
fail_max = None
pass_min = None
for (r_lo, r_hi, p, r, wins) in results:
    if wins is None:
        continue
    if wins < 3:
        if fail_max is None or r > fail_max:
            fail_max = r
    else:
        if pass_min is None or r < pass_min:
            pass_min = r

if fail_max is not None and pass_min is not None:
    midpoint = (fail_max + pass_min) / 2
    print(f"  Highest confirmed FAIL:  λ/n = {fail_max:.4f}")
    print(f"  Lowest confirmed PASS:   λ/n = {pass_min:.4f}")
    print(f"  Threshold estimate:      λ/n ≈ {midpoint:.4f} (midpoint)")
    print(f"  Threshold range:         ({fail_max:.4f}, {pass_min:.4f})")
elif fail_max is None:
    print("  All bins passed — threshold below tested range")
elif pass_min is None:
    print("  All bins failed — threshold above tested range")
print()

# k-bounds reference
print("k-bounds used (per curve):")
for bin_idx, (r_lo, r_hi) in enumerate(TARGET_BINS):
    curves = bins[bin_idx]
    if not curves:
        continue
    p, b, n, lam, r = curves[0]
    k1, k2 = make_bounds(n)
    eff = k1 * k2 / n
    print(f"  [λ/n={r:.3f}] n={n}: K1={k1}, K2={k2}, eff={eff:.4f}, "
          f"m_thresh≈{math.ceil(math.log(n)/math.log(1/eff)) if eff < 1 else '∞'}")

print()

# ---------------------------------------------------------------------------
# Phase 2: Fine bisection in the [0.04, 0.13] region where threshold lies
# ---------------------------------------------------------------------------

print("=" * 70)
print("PHASE 2: Fine bisection in [0.04, 0.13]")
print("=" * 70)

# Add the known-fail anchor: p=2677, n=2647, lam=185, lam/n=0.0699
# from the 2026-07-26 session (confirmed fail with BKZ-40 too).
ANCHOR_FAIL = (2677, 2, 2647, 185)   # (p, b, n, lam)
print("\nAnchor (known-fail from 2026-07-26): p=2677, n=2647, lam=185, λ/n=0.0699")
G_anchor = find_generator(2677, 2, 2647)
if G_anchor is not None:
    anchor_k1, anchor_k2 = make_bounds(2647)
    anchor_wins = sum(
        1 for seed in SEEDS
        if run_attack(2677, 2, 2647, 185, anchor_k1, anchor_k2, M_SIGS, seed=seed)
    )
    print(f"  LLL m={M_SIGS}: {anchor_wins}/{len(SEEDS)}  (expected: 0/5)")

fine_bins_db = {b: [] for b in range(len(FINE_BINS))}
fine_count = {b: 0 for b in range(len(FINE_BINS))}
FINE_SEEDS = SEEDS + [271828, 31415, 12321, 77777, 100003]  # 10 seeds total

# Search in a wider prime range for fine bins
print("\nSearching for fine-grain curves in [2048, 65536]...")
p_cur = sympy.nextprime(2048)
fp_checked = 0
while p_cur < 65536:
    fp_checked += 1
    if p_cur % 6 == 1:
        eis = eisenstein_decompose(p_cur)
        if eis is not None:
            a_e, b_e = eis
            for t in j0_traces(a_e, b_e):
                n_cand = p_cur + 1 - t
                if n_cand <= 1: continue
                if not sympy.isprime(n_cand): continue
                if n_cand % 3 != 1: continue
                lam_c, lam2_c = glv_eigenvalue(n_cand)
                if lam_c is None: continue
                r_min = min(lam_c, lam2_c) / n_cand
                # Use the SMALLER eigenvalue (may be the bottleneck)
                lam_use_c = lam_c if lam_c / n_cand == r_min else lam2_c
                for bin_idx, (r_lo, r_hi) in enumerate(FINE_BINS):
                    if fine_count[bin_idx] >= 1: continue
                    if r_lo <= r_min < r_hi:
                        b_cur = find_b_for_order(p_cur, n_cand)
                        if b_cur is None: continue
                        fine_bins_db[bin_idx].append((p_cur, b_cur, n_cand, lam_use_c, r_min))
                        fine_count[bin_idx] += 1
                        print(f"  Fine [{r_lo:.3f},{r_hi:.3f}): p={p_cur}, n={n_cand}, "
                              f"λ={lam_use_c}, λ/n={r_min:.5f}")
                        break
    if all(fine_count[b] >= 1 for b in range(len(FINE_BINS))):
        break
    p_cur = sympy.nextprime(p_cur)

print(f"  Searched {fp_checked} primes for fine bins.\n")

print(f"Fine sweep (m={M_SIGS}, {len(FINE_SEEDS)} seeds per curve):")
print(f"{'Fine bin':22s}  {'p':>8}  {'λ/n':>7}  {'LLL':>8}")
print("-" * 55)

fine_results = []
for bin_idx, (r_lo, r_hi) in enumerate(FINE_BINS):
    if not fine_bins_db[bin_idx]:
        print(f"[{r_lo:.3f}, {r_hi:.3f})  {'--':>8}  {'--':>7}  N/A")
        fine_results.append(None)
        continue
    p_f, b_f, n_f, lam_f, r_f = fine_bins_db[bin_idx][0]
    k1_f, k2_f = make_bounds(n_f)
    wins_f = sum(
        1 for seed in FINE_SEEDS
        if run_attack(p_f, b_f, n_f, lam_f, k1_f, k2_f, M_SIGS, seed=seed)
    )
    flag_f = "PASS" if wins_f >= 6 else ("MARGINAL" if wins_f >= 2 else "FAIL")
    print(f"[{r_lo:.3f}, {r_hi:.3f})  {p_f:>8}  {r_f:.5f}  "
          f"{wins_f}/{len(FINE_SEEDS)}  ({flag_f})")
    fine_results.append((r_lo, r_hi, p_f, r_f, wins_f))

print()

# Fine bisection summary
print("Fine bisection threshold:")
fine_fail = None
fine_pass = None
for item in fine_results:
    if item is None: continue
    r_lo, r_hi, p_f, r_f, wins_f = item
    if wins_f < 6:
        if fine_fail is None or r_f > fine_fail:
            fine_fail = r_f
    else:
        if fine_pass is None or r_f < fine_pass:
            fine_pass = r_f

if fine_fail is not None and fine_pass is not None:
    print(f"  Highest fine FAIL:  λ/n = {fine_fail:.5f}")
    print(f"  Lowest fine PASS:   λ/n = {fine_pass:.5f}")
    print(f"  Threshold window:   ({fine_fail:.5f}, {fine_pass:.5f})")
elif fine_pass is None:
    print("  All fine bins failed — threshold above fine range")
elif fine_fail is None:
    print("  All fine bins passed — threshold below fine range")

# ---------------------------------------------------------------------------
# Phase 3: Eigenvalue choice check — does using max(lam) fix marginal cases?
# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("PHASE 3: Eigenvalue choice — check if lam_max fixes marginal coarse bins")
print("=" * 70)

MARGINAL_BINS = [(r_lo, r_hi, p, r, wins)
                 for (r_lo, r_hi, p, r, wins) in results
                 if wins is not None and 0 < wins < 3]

if not MARGINAL_BINS:
    print("  No marginal coarse bins to investigate.")
else:
    print(f"  Marginal bins (1-2/5): testing with lam_max")
    # For each marginal bin, find the curve and try the OTHER eigenvalue
    for (r_lo, r_hi, p_m, r_m, wins_m) in MARGINAL_BINS:
        bin_idx = next(i for i, (rl, rh) in enumerate(TARGET_BINS)
                       if rl == r_lo and rh == r_hi)
        if not bins[bin_idx]:
            continue
        p_c, b_c, n_c, lam_c, r_c = bins[bin_idx][0]
        # Other eigenvalue
        lam_other = (n_c - 1 - lam_c) % n_c
        r_other = lam_other / n_c

        # Regenerate signatures using lam_other (signer uses lam_other)
        print(f"\n  Bin [{r_lo:.2f},{r_hi:.2f}): p={p_c}, n={n_c}")
        print(f"    lam_min={lam_c} (ratio={r_c:.4f}): {wins_m}/{len(SEEDS)} (marginal)")
        wins_other = sum(
            1 for seed in SEEDS
            if run_attack(p_c, b_c, n_c, lam_other, make_bounds(n_c)[0],
                          make_bounds(n_c)[1], M_SIGS, seed=seed)
        )
        flag_o = "PASS" if wins_other >= 3 else ("MARGINAL" if wins_other >= 1 else "FAIL")
        print(f"    lam_max={lam_other} (ratio={r_other:.4f}): "
              f"{wins_other}/{len(SEEDS)} ({flag_o})")

print()
print("Done.")

