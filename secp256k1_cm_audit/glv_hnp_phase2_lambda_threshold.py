"""
GLV-HNP Phase 2: lambda/n threshold study.

Systematically searches for j=0 GLV curves at various lambda/n ratios
in [0.05, 0.50] using small (10-14 bit) primes, tests LLL recovery,
and identifies the critical threshold.

Known data points from prior runs (2026-07-26):
  lam/n = 0.53  (8-bit,  n=199,    p=211):  LLL 3/3 at m=4
  lam/n = 0.66  (12-bit, n=2659,   p=2557): LLL 3/3 at m=7
  lam/n = 0.34  (20-bit, n=523969, p=...):  LLL 3/3 at m=9
  lam/n = 0.07  (12-bit, n=2647,   p=2677): LLL NEVER 3/3 (up to m=12, BKZ-40 also fails)

Goal: bisect [0.07, 0.34] to identify the critical lambda/n below which
LLL (and BKZ-20) fail, regardless of m.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass a=0, y^2 = x^3 + b)
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
    m_ts, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_ts - i - 1), p)
        m_ts, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(12345 + p)
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
# CM theory
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
# Curve search: find j=0 curves with lam/n in a target band
# ---------------------------------------------------------------------------

def find_curves_in_band(lam_n_target, band_half=0.02, bit_lo=10, bit_hi=14,
                        max_curves=3):
    """
    Search for j=0 GLV curves with lam/n in
    [lam_n_target - band_half, lam_n_target + band_half].

    Returns list of (p, b, n, lam, G, lam/n).
    Stops when max_curves found.
    """
    lo = lam_n_target - band_half
    hi = lam_n_target + band_half
    candidates = []

    p = sympy.nextprime(2**bit_lo)
    p_limit = 2**bit_hi
    while p < p_limit and len(candidates) < max_curves:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 4: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    if lo <= ratio <= hi:
                        # Find b parameter (twist class)
                        found_b = None
                        for b_try in range(1, min(p, 500)):
                            rhs_test = (pow(b_try, 3, p) + b_try) % p  # dummy check
                            # Actually find a point on y^2 = x^3 + b_try
                            rng_tmp = random.Random(77 + b_try + p)
                            for _ in range(200):
                                x = rng_tmp.randint(0, p - 1)
                                rhs = (pow(x, 3, p) + b_try) % p
                                y = tonelli_shanks(rhs, p)
                                if y is not None and y != 0:
                                    P = (x, y)
                                    if ec_mul(P, n_cand, p) is None:
                                        found_b = b_try
                                        break
                            if found_b is not None:
                                break
                        if found_b is None:
                            continue
                        G = find_generator(p, found_b, n_cand)
                        if G is None:
                            continue
                        candidates.append((p, found_b, n_cand, lam, G, ratio))
                        if len(candidates) >= max_curves:
                            break
        p = sympy.nextprime(p)
    return candidates

# ---------------------------------------------------------------------------
# Lattice attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs = []
    for _ in range(500000):
        if len(sigs) >= m: break
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
    M[m][m] = 1  # S_D = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def recover_d(M_reduced, m, n, d_secret):
    dim = 2 * m + 2
    S_KANNAN = n
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KANNAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_attack(curve_params, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G, ratio = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, d_secret)

def sweep_m(curve_params, k1_bound, seeds, m_range, use_bkz=False, bkz_beta=20):
    """
    For a given curve, find the minimum m at which LLL achieves 3/3.
    Returns (min_m_for_3of3, results_dict) where results_dict[m] = (wins, total).
    """
    p, b, n, lam, G, ratio = curve_params
    results = {}
    min_m = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 3333 + n).randint(1, n - 1)
            ok = run_attack(curve_params, m, d_trial, k1_bound, seed,
                            use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
        if wins == len(seeds) and min_m is None:
            min_m = m
    return min_m, results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: lambda/n threshold study")
print("Locates the critical lambda/n below which LLL fails regardless of m")
print("=" * 70)
print()

SEEDS = [42, 1234, 9999]
# For each target band, we test 1-2 representative curves.
# m range: 4..16 (generous upper bound)
M_RANGE = range(4, 17)

# K1_BOUND choice: eff = K1 * K2 / n, want eff ~ 0.02..0.08
# For 12-bit curves (n ~3000, sqrt(n)~55): K2~55, K1 ~ 0.04 * 3000/55 ~ 2..4
# We use K1_BOUND = 4 uniformly for 12-bit curves for comparability.
K1_BOUND = 4

# Target lambda/n bands to probe.
# We have: 0.07 fails, 0.34 succeeds. Bisect [0.07, 0.34].
TARGET_RATIOS = [0.10, 0.14, 0.18, 0.22, 0.26, 0.30]
# Also add known-bad (0.07) and known-good (0.34) as reference points
REFERENCE_RATIOS = [0.07, 0.34]

ALL_RATIOS = [(r, False) for r in REFERENCE_RATIOS] + \
             [(r, True)  for r in TARGET_RATIOS]  # (ratio, is_new)

results_summary = []

for target_ratio, is_new in ALL_RATIOS:
    label = "NEW " if is_new else "REF "
    print(f"\n{'='*70}")
    print(f"{label}Target lambda/n ≈ {target_ratio:.2f}")
    print(f"{'='*70}")

    # Use a wider band for reference points to ensure we find them
    band = 0.025 if is_new else 0.035

    curves = find_curves_in_band(target_ratio, band_half=band,
                                 bit_lo=11, bit_hi=14, max_curves=2)

    if not curves:
        print(f"  No curve found in band [{target_ratio-band:.3f}, {target_ratio+band:.3f}]")
        # Try wider search
        curves = find_curves_in_band(target_ratio, band_half=0.05,
                                     bit_lo=10, bit_hi=15, max_curves=1)
    if not curves:
        print(f"  SKIP: No curve found even with band ±0.05")
        results_summary.append((target_ratio, None, None, None))
        continue

    # Take the first curve found
    curve = curves[0]
    p, b, n, lam, G, actual_ratio = curve
    k2 = math.isqrt(n) + 1
    eff = K1_BOUND * k2 / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else 999

    print(f"  Curve: p={p}, b={b}, n={n} ({n.bit_length()}b)")
    print(f"  lam={lam}, lam/n={actual_ratio:.4f}")
    print(f"  K1={K1_BOUND}, K2={k2}, eff={eff:.4f}, m_thresh≈{m_thresh:.1f}")

    # LLL sweep
    print(f"  Running LLL sweep (m={M_RANGE.start}..{M_RANGE.stop-1}, {len(SEEDS)} seeds)...")
    min_m_lll, res_lll = sweep_m(curve, K1_BOUND, SEEDS, M_RANGE, use_bkz=False)

    lll_summary = []
    for m, (wins, total) in sorted(res_lll.items()):
        marker = " ←3/3" if wins == total else ""
        lll_summary.append(f"m={m}:{wins}/{total}{marker}")
    print(f"  LLL:    {', '.join(lll_summary)}")

    if min_m_lll is None:
        # LLL failed: try BKZ-20
        print(f"  LLL FAILED at all m. Trying BKZ(beta=20)...")
        min_m_bkz, res_bkz = sweep_m(curve, K1_BOUND, SEEDS, M_RANGE,
                                      use_bkz=True, bkz_beta=20)
        bkz_summary = []
        for m, (wins, total) in sorted(res_bkz.items()):
            marker = " ←3/3" if wins == total else ""
            bkz_summary.append(f"m={m}:{wins}/{total}{marker}")
        print(f"  BKZ20:  {', '.join(bkz_summary)}")
        outcome = "FAIL"
        min_m_winner = None
    else:
        min_m_bkz = None
        outcome = "SUCCESS"
        min_m_winner = min_m_lll

    results_summary.append((actual_ratio, outcome, min_m_winner, m_thresh))
    print(f"  => {outcome} | min_m_3of3={min_m_winner}")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'lam/n':>8}  {'outcome':>10}  {'min_m_3of3':>12}  {'m_thresh':>10}")
print("-" * 50)
for ratio, outcome, min_m, m_thresh in results_summary:
    mt = f"{m_thresh:.0f}" if m_thresh is not None else "N/A"
    mm = str(min_m) if min_m is not None else "—"
    out = outcome if outcome else "NO CURVE"
    print(f"{ratio:>8.3f}  {out:>10}  {mm:>12}  {mt:>10}")

print()
print("Threshold estimate: the critical lam/n is between the highest FAIL")
print("and lowest SUCCESS in the table above.")
print()
print("Done.")
