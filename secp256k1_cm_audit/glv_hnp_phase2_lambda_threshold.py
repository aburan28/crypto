"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Bisects the threshold between:
  - λ/n = 0.07 (p=2677, n=2647):   FAILS (LLL, BKZ-20, BKZ-40 all fail)
  - λ/n = 0.34 (p=524347, n=523969): SUCCEEDS at m=9

Method: find one 12-14 bit j=0 prime-order curve per 0.05-wide λ/n bin,
then sweep LLL at m=4..14 (3 seeds each) to find the 3/3-recovery threshold.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass, a=0)
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

def find_point(p, b, rng=None):
    rng = rng or random.Random(42)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n, rng=None):
    rng = rng or random.Random(12345)
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
# CM theory: Eisenstein decomposition + GLV eigenvalue
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
    """Return the smaller GLV eigenvalue lambda < n/2."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

# ---------------------------------------------------------------------------
# Find one curve per λ/n bin
# ---------------------------------------------------------------------------

def find_curves_by_ratio_bins(bins, min_p=2**11, max_p=2**15):
    """
    Find one j=0 prime-order curve per bin in `bins`.
    bins: list of (lo, hi) for λ/n range.
    Returns: list of (bin_lo, bin_hi, p, b, n, lam, G) or None per bin.
    """
    results = {(lo, hi): None for lo, hi in bins}
    remaining = set((lo, hi) for lo, hi in bins)

    p = sympy.nextprime(min_p - 1)
    searched = 0
    while p <= max_p and remaining:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 4 or n_cand >= 2 * max_p:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    for (lo, hi) in list(remaining):
                        if lo <= ratio < hi:
                            # Find b such that #E(F_p) = n_cand
                            found_b = None
                            rng_b = random.Random(p + n_cand)
                            for b_try in range(1, min(p, 500)):
                                pt = find_point(p, b_try, random.Random(42 + b_try))
                                if pt is None:
                                    continue
                                if ec_mul(pt, n_cand, p) is None:
                                    found_b = b_try
                                    break
                            if found_b is None:
                                continue
                            G = find_generator(p, found_b, n_cand, random.Random(p))
                            if G is None:
                                continue
                            results[(lo, hi)] = (p, found_b, n_cand, lam, G)
                            remaining.discard((lo, hi))
                            break
        p = sympy.nextprime(p)
        searched += 1

    return [(lo, hi, results[(lo, hi)]) for lo, hi in bins]

# ---------------------------------------------------------------------------
# Lattice build + LLL recovery
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
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

def run_single(p, b_curve, n, lam, G, m, d_secret, k1_bound, k2_bound, seed, use_bkz=False, bkz_beta=20):
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b_curve, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, _, _, _, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep_lll(p, b_curve, n, lam, G, m_range, seeds, use_bkz=False, bkz_beta=20):
    """Returns dict m -> (wins, total)."""
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 9999).randint(1, n - 1)
            ok = run_single(p, b_curve, n, lam, G, m, d, k1_bound, k2_bound, seed,
                            use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
    return results, k1_bound, k2_bound

def first_3of3(results):
    for m, (w, t) in sorted(results.items()):
        if w == t:
            return m
    return None

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 72)
print("GLV-HNP Phase 2: λ/n threshold bisection study")
print("Known: λ/n=0.07 FAILS (BKZ-40 too); λ/n=0.34 SUCCEEDS; λ/n=0.53 SUCCEEDS")
print("Goal: find threshold in [0.07, 0.34] using intermediate-ratio curves")
print("=" * 72)

# λ/n bins to probe (lo, hi)
BINS = [
    (0.05, 0.10),
    (0.10, 0.15),
    (0.15, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.35),
    (0.35, 0.50),
]

SEEDS = [42, 1234, 9999]
M_RANGE = list(range(4, 16))    # m=4..15

# ---- Phase 1: Find curves --------------------------------------------------
print("\nPhase 1: Searching for one j=0 curve per λ/n bin (12-15 bit)...")
bin_curves = find_curves_by_ratio_bins(BINS, min_p=2**11, max_p=2**15)

print(f"\n{'Bin':<12} {'p':>7} {'n':>7} {'bits':>5} {'lam/n':>8} {'lam':>8}")
print("-" * 56)
found_curves = []
for lo, hi, curve in bin_curves:
    if curve is None:
        print(f"[{lo:.2f},{hi:.2f})   NOT FOUND")
    else:
        p, b_c, n, lam, G = curve
        ratio = lam / n
        bits = n.bit_length()
        print(f"[{lo:.2f},{hi:.2f})  {p:>7} {n:>7} {bits:>5}  {ratio:>8.4f} {lam:>8}")
        found_curves.append((lo, hi, p, b_c, n, lam, G))

# ---- Phase 2: LLL sweep for each found curve -------------------------------
print("\n\nPhase 2: LLL sweep (m=4..15, 3 seeds each)")
print("K1_BOUND = max(2, 0.05·√n), K2_BOUND = ⌊√n⌋+1")
print("=" * 72)

SUMMARY = []

for lo, hi, p, b_c, n, lam, G in found_curves:
    ratio = lam / n
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    eff = k1_bound * k2_bound / n
    m_thresh = (math.ceil(math.log(n) / math.log(1.0 / eff))
                if eff < 1 else float('inf'))

    print(f"\n{'='*72}")
    print(f"Bin [{lo:.2f},{hi:.2f})  |  p={p}, n={n} ({n.bit_length()}b), "
          f"lam={lam}, lam/n={ratio:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")

    results, k1b, k2b = sweep_lll(p, b_c, n, lam, G, M_RANGE, SEEDS)
    first = first_3of3(results)

    for m, (w, t) in sorted(results.items()):
        marker = " ← 3/3" if w == t and (first is None or m == first) else ""
        above = "(≥thresh)" if m >= m_thresh else "(<thresh)"
        print(f"  m={m:2d}: {w}/{t}  {above}{marker}")

    status = f"3/3 at m={first}" if first else "NEVER 3/3 (m≤15)"
    print(f"  Summary: {status}")
    SUMMARY.append((lo, hi, ratio, first, p, n, lam))

# ---- Phase 3: Known anchor points (pre-confirmed) --------------------------
print("\n\nPhase 3: Anchor points (from prior sessions)")
print("=" * 72)

ANCHORS = [
    # (label, lo, hi, lam/n, m_thresh_3of3_or_None)
    ("p=211,  n=199,    lam=106",  0.53, "3/3 at m=4  (8-bit, Session 2026-07-26)"),
    ("p=2557, n=2659,   lam=1755", 0.66, "3/3 at m=7  (12-bit, Session 2026-07-26)"),
    ("p=524347,n=523969,lam=177902", 0.34, "3/3 at m=9  (20-bit, Session 2026-07-26)"),
    ("p=2677, n=2647,   lam=185",   0.07, "NEVER 3/3   (LLL+BKZ-40, Session 2026-07-26)"),
]
for label, ratio, result in ANCHORS:
    print(f"  lam/n={ratio:.2f}  {label:40s}  {result}")

# ---- Phase 4: Threshold analysis -------------------------------------------
print("\n\nPhase 4: Threshold analysis")
print("=" * 72)

# Combine all results
all_pts = []
for lo, hi, ratio, first, p, n, lam in SUMMARY:
    all_pts.append((ratio, first))
# Add anchors
all_pts += [(0.53, 4), (0.66, 7), (0.34, 9), (0.07, None)]
all_pts.sort(key=lambda x: x[0])

print(f"\n{'lam/n':>8} {'3/3 at m':>10}  status")
print("-" * 40)
for ratio, first in all_pts:
    status = f"3/3 at m={first}" if first else "FAIL (m≤15)"
    marker = " ← known anchor" if ratio in (0.53, 0.66, 0.34, 0.07) else ""
    print(f"  {ratio:6.4f}  {str(first) if first else 'None':>8}   {status}{marker}")

# Find the threshold interval
successes = [(r, m) for r, m in all_pts if m is not None]
failures  = [(r, m) for r, m in all_pts if m is None]

if successes and failures:
    lo_fail = max(r for r, _ in failures)
    hi_ok   = min(r for r, _ in successes)
    print(f"\nConclusion:")
    print(f"  Last failure: lam/n = {lo_fail:.4f}")
    print(f"  First success: lam/n = {hi_ok:.4f}")
    print(f"  Threshold in ({lo_fail:.4f}, {hi_ok:.4f})")
    mid = (lo_fail + hi_ok) / 2
    print(f"  Midpoint estimate: lam/n ≈ {mid:.4f}")
elif failures:
    print("\nAll tested ratios FAIL — need smaller step or larger m.")
elif successes:
    print("\nAll tested ratios SUCCEED — threshold is below the smallest tested.")
else:
    print("\nNo results to analyze.")

# ---- Phase 5: BKZ rescue on threshold region --------------------------------
# If there's a failure at some ratio, try BKZ-20 and BKZ-40 on it.
print("\n\nPhase 5: BKZ rescue on first failure above λ/n=0.07")
print("=" * 72)

# Find lowest-ratio failure from new results
new_fails = [(ratio, p, b_c, n, lam, G)
             for (lo, hi, ratio, first, p, n, lam) in SUMMARY
             if first is None]
new_successes = [(ratio, p, b_c, n, lam, G)
                 for (lo, hi, ratio, first, p, n, lam) in SUMMARY
                 if first is not None]

# Rebuild curve lookup for BKZ
curve_by_ratio = {ratio: (p, b_c, n, lam, G)
                  for lo, hi, ratio, first, p, n, lam in SUMMARY}
_G_map = {}
for lo, hi, p, b_c, n, lam, G in found_curves:
    ratio = lam / n
    _G_map[ratio] = (p, b_c, n, lam, G)

if new_fails:
    highest_fail_ratio = max(r for r, *_ in new_fails)
    # get the curve
    for lo, hi, p, b_c, n, lam, G in found_curves:
        if abs(lam / n - highest_fail_ratio) < 1e-9:
            target_curve = (p, b_c, n, lam, G)
            break
    else:
        target_curve = None

    if target_curve:
        p_t, b_t, n_t, lam_t, G_t = target_curve
        print(f"Testing BKZ rescue on failure curve: p={p_t}, n={n_t}, "
              f"lam={lam_t}, lam/n={lam_t/n_t:.4f}")
        for bkz_b in [20, 40]:
            res_bkz, k1b, k2b = sweep_lll(p_t, b_t, n_t, lam_t, G_t,
                                           range(6, 16), SEEDS,
                                           use_bkz=True, bkz_beta=bkz_b)
            first_bkz = first_3of3(res_bkz)
            status_b = f"3/3 at m={first_bkz}" if first_bkz else "NEVER 3/3"
            print(f"  BKZ(beta={bkz_b}): {status_b}")
    else:
        print("  (Could not reconstruct failure curve for BKZ test)")
else:
    print("  No new failures found — all tested λ/n bins succeed with LLL.")

print("\n\nDone.")
