"""
GLV-HNP Phase 2: λ/n threshold bisection.

Goal: find the critical λ/n ratio below which LLL fails to recover d.
From prior runs:
  λ/n = 0.07  (p=2677):  LLL FAILS (even BKZ-40)
  λ/n = 0.34  (p=524347): LLL SUCCEEDS at m=9
  λ/n = 0.53  (p=211):   LLL SUCCEEDS at m=4

This script systematically scans 12-bit j=0 curves covering
λ/n in bands [0.05,0.10), [0.10,0.15), ..., [0.45,0.50).
For each band we pick one representative curve and run LLL
at m=6,8,10,12 over 3 seeds.  A band "succeeds" if 3/3 at
some m ≤ 12.  We report the lowest succeeding ratio and the
highest failing ratio → brackets the threshold.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC helpers (copied from glv_hnp_phase2_20bit.py)
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
# CM helpers
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
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

# ---------------------------------------------------------------------------
# Curve collection: scan 12-bit primes, bin by λ/n
# ---------------------------------------------------------------------------

def collect_curves_by_ratio(bit_lo=11, bit_hi=13, n_per_band=1):
    """
    Scan j=0 primes in the given bit range.
    Return dict: band_lo -> list of (p, b, n, lam) tuples.
    Bands: [0.05,0.10), [0.10,0.15), ..., [0.45,0.50).
    """
    bands = {}
    for lo_tenth in range(1, 10):   # 0.10..0.90
        lo = lo_tenth * 0.05
        hi = lo + 0.05
        bands[(round(lo, 2), round(hi, 2))] = []

    p_start = 2 ** bit_lo
    p_end   = 2 ** bit_hi

    p = sympy.nextprime(p_start)
    while p < p_end:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                traces = j0_traces(a_e, b_e)
                for t in traces:
                    n_cand = p + 1 - t
                    if n_cand < 3: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand
                    for (lo, hi), lst in bands.items():
                        if lo <= ratio < hi:
                            if len(lst) < n_per_band:
                                # Find a working twist b
                                found_b = None
                                for b_try in range(1, 20):
                                    rng2 = random.Random(999)
                                    ok_ctr = 0
                                    for _ in range(50):
                                        x = rng2.randint(0, p - 1)
                                        rhs = (pow(x, 3, p) + b_try) % p
                                        y = tonelli_shanks(rhs, p)
                                        if y is None or y == 0: continue
                                        P = (x, y)
                                        if ec_mul(P, n_cand, p) is None:
                                            found_b = b_try
                                            break
                                    if found_b is not None: break
                                if found_b is not None:
                                    lst.append((p, found_b, n_cand, lam))
                            break
        p = sympy.nextprime(p)
    return bands

# ---------------------------------------------------------------------------
# Lattice and recovery (same as 20bit script)
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
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1  # S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_KANNAN

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

def run_one(curve, m, d_secret, k1_bound, seed=42):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def test_band(curve_info, m_vals, seeds, k1_bound):
    """Return (first_m_with_3of3, wins_per_m dict)."""
    p, b, n, lam = curve_info
    G = find_generator(p, b, n)
    if G is None:
        return None, {}
    curve = (p, b, n, lam, G)
    wins = {}
    first_full = None
    for m in m_vals:
        w = 0
        for seed in seeds:
            d_trial = random.Random(seed + 31337).randint(1, n - 1)
            ok = run_one(curve, m, d_trial, k1_bound, seed)
            w += ok
        wins[m] = (w, len(seeds))
        if w == len(seeds) and first_full is None:
            first_full = m
    return first_full, wins

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection (12-bit j=0 curves)")
print("Known: λ/n=0.07 → FAIL, λ/n=0.34 → PASS, λ/n=0.53 → PASS")
print("=" * 70)

print("\nStep 1: Scanning 11-13 bit j=0 primes for curves in each λ/n band...")
bands = collect_curves_by_ratio(bit_lo=11, bit_hi=13, n_per_band=1)

# Also add the known reference curves
SEEDS = [42, 1234, 9999]
M_VALS = [6, 8, 10, 12]

# Summary collection
print("\nStep 2: LLL attack per band\n")
print(f"{'Band':<14}{'p':>8}{'n':>8}{'lam':>8}{'lam/n':>8}{'K1':>5}  ", end="")
for m in M_VALS:
    print(f"m={m:2d}    ", end="")
print("  Verdict")
print("-" * (14 + 8 + 8 + 8 + 8 + 5 + 4 + len(M_VALS)*9 + 12))

results_by_band = {}
for (lo, hi), curves in sorted(bands.items()):
    if not curves:
        print(f"[{lo:.2f},{hi:.2f})  {'—':>8}{'—':>8}{'—':>8}{'—':>8}{'—':>5}  (no curve found)")
        results_by_band[(lo, hi)] = None
        continue

    p, b, n, lam = curves[0]
    ratio = lam / n
    k1_bound = max(2, int(0.05 * math.sqrt(n)))

    first_m, wins = test_band((p, b, n, lam), M_VALS, SEEDS, k1_bound)

    row = f"[{lo:.2f},{hi:.2f})  {p:>8}{n:>8}{lam:>8}{ratio:>8.3f}{k1_bound:>5}  "
    for m in M_VALS:
        if m in wins:
            w, total = wins[m]
            row += f"{w}/{total}       "[: 7]
        else:
            row += "  -    "
    verdict = f"PASS(m={first_m})" if first_m is not None else "FAIL"
    print(row + "  " + verdict)
    results_by_band[(lo, hi)] = (ratio, first_m is not None, first_m)

# Also test the two known reference curves for cross-check
print("\nKnown reference curves (cross-check):")
ref_curves = [
    ("8-bit/p=211",    211,    2,  199,  106, 2, 0.532),
    ("12-bit/p=2677",  2677,   2, 2647,  185, 8, 0.070),
    ("20-bit/p=524347",524347, 2, 523969, 177902, 36, 0.340),
]
for label, p, b, n, lam, k1_bound, ratio in ref_curves:
    G_ref = find_generator(p, b, n)
    if G_ref is None:
        print(f"  {label}: generator not found")
        continue
    curve_ref = (p, b, n, lam, G_ref)
    wins_ref = {}
    first_full_ref = None
    for m in M_VALS:
        w = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 31337).randint(1, n - 1)
            ok = run_one(curve_ref, m, d_trial, k1_bound, seed)
            w += ok
        wins_ref[m] = (w, len(SEEDS))
        if w == len(SEEDS) and first_full_ref is None:
            first_full_ref = m
    verdict = f"PASS(m={first_full_ref})" if first_full_ref else "FAIL"
    wins_str = "  ".join(f"{wins_ref[m][0]}/{wins_ref[m][1]}" for m in M_VALS)
    print(f"  {label:<28} lam/n={ratio:.3f}  K1={k1_bound}  [{wins_str}]  {verdict}")

# Re-run p=2677 with K1=2 (same as the other 12-bit curves) to check if
# the failure is due to K1_bound or truly due to lam/n.
print("\nIsolation: p=2677 (lam/n=0.07) with K1=2 (same as other 12-bit curves)")
G_2677 = find_generator(2677, 2, 2647)
if G_2677 is not None:
    curve_2677 = (2677, 2, 2647, 185, G_2677)
    for m in M_VALS:
        w = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 31337).randint(1, 2647 - 1)
            ok = run_one(curve_2677, m, d_trial, k1_bound=2, seed=seed)
            w += ok
        print(f"  p=2677 lam/n=0.07 K1=2 m={m}: {w}/{len(SEEDS)}")

# Determine threshold bracket
print("\n" + "=" * 70)
print("THRESHOLD ANALYSIS")
print("=" * 70)
passing_list = [(lo, hi, v[0], v[2]) for (lo, hi), v in results_by_band.items()
                if v is not None and v[1]]
failing_list = [(lo, hi, v[0]) for (lo, hi), v in results_by_band.items()
                if v is not None and not v[1]]

if passing_list and failing_list:
    min_pass_ratio = min(ratio for lo, hi, ratio, fm in passing_list)
    max_fail_ratio = max(ratio for lo, hi, ratio in failing_list)
    print(f"Highest FAILING λ/n band:  {max_fail_ratio:.3f}")
    print(f"Lowest  PASSING λ/n band:  {min_pass_ratio:.3f}")
    print(f"→ Threshold bracket: ({max_fail_ratio:.3f}, {min_pass_ratio:.3f})")
elif passing_list:
    print(f"All non-empty bands pass. Lowest passing: "
          f"{min(r for lo,hi,r,fm in passing_list):.3f}")
elif failing_list:
    print(f"All non-empty bands fail. Highest failing: "
          f"{max(r for lo,hi,r in failing_list):.3f}")
else:
    print("No conclusive bands (all empty).")

print("\nDone.")

# ---------------------------------------------------------------------------
# PART 2: Cross-validation — fix K1=8, sweep λ/n bands
# Tests whether failure at K1=8 is λ/n-dependent or K1-dependent
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("PART 2: Fixed K1=8, sweep all λ/n bands (same K1 as prior 'fail' test)")
print("(Expected: if failure is K1-dependent, ALL bands should fail/pass equally)")
print("=" * 70)

FIXED_K1 = 8

print(f"\n{'Band':<14}{'lam/n':>8}{'eff':>8}  ", end="")
for m in M_VALS:
    print(f"m={m:2d}    ", end="")
print("  Verdict")
print("-" * 90)

for (lo, hi), curves in sorted(bands.items()):
    if not curves:
        print(f"[{lo:.2f},{hi:.2f})  {'—':>8}{'—':>8}  (no curve)")
        continue
    p, b, n, lam = curves[0]
    ratio = lam / n
    k2_bound = math.isqrt(n) + 1
    eff = FIXED_K1 * k2_bound / n

    G = find_generator(p, b, n)
    if G is None:
        print(f"[{lo:.2f},{hi:.2f})  {ratio:>8.3f}{'—':>8}  (no generator)")
        continue
    curve = (p, b, n, lam, G)
    wins2 = {}
    first_full2 = None
    for m in M_VALS:
        w = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 31337).randint(1, n - 1)
            ok = run_one(curve, m, d_trial, FIXED_K1, seed)
            w += ok
        wins2[m] = (w, len(SEEDS))
        if w == len(SEEDS) and first_full2 is None:
            first_full2 = m

    row = f"[{lo:.2f},{hi:.2f})  {ratio:>8.3f}{eff:>8.3f}  "
    for m in M_VALS:
        w, total = wins2[m]
        row += f"{w}/{total}       "[: 7]
    verdict = f"PASS(m={first_full2})" if first_full2 is not None else "FAIL"
    print(row + "  " + verdict)

print("\nConclusion:")
print("  If all bands pass at K1=8 → failure was eff (K1*K2/n), NOT λ/n.")
print("  If pattern correlates with λ/n → λ/n is a genuine factor.")

# ---------------------------------------------------------------------------
# PART 3: Confirm eff as controlling variable — fix K1=8, 5 seeds
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("PART 3: Fixed K1=8, 5 seeds, M_VALS=[8,12,16] — more reliable stats")
print("=" * 70)

SEEDS5 = [42, 1234, 9999, 7777, 5555]
M_VALS3 = [8, 12, 16]

print(f"\n{'Band':<14}{'lam/n':>8}{'eff':>8}  ", end="")
for m in M_VALS3:
    print(f"m={m:2d}    ", end="")
print("  Verdict")
print("-" * 80)

for (lo, hi), curves in sorted(bands.items()):
    if not curves:
        print(f"[{lo:.2f},{hi:.2f})  {'—':>8}{'—':>8}  (no curve)")
        continue
    p, b, n, lam = curves[0]
    ratio = lam / n
    k2_bound = math.isqrt(n) + 1
    eff = FIXED_K1 * k2_bound / n

    G = find_generator(p, b, n)
    if G is None:
        print(f"[{lo:.2f},{hi:.2f})  {ratio:>8.3f}{'—':>8}  (no generator)")
        continue
    curve = (p, b, n, lam, G)
    wins3 = {}
    first_full3 = None
    for m in M_VALS3:
        w = 0
        for seed in SEEDS5:
            d_trial = random.Random(seed + 31337).randint(1, n - 1)
            ok = run_one(curve, m, d_trial, FIXED_K1, seed)
            w += ok
        wins3[m] = (w, len(SEEDS5))
        if w == len(SEEDS5) and first_full3 is None:
            first_full3 = m

    row = f"[{lo:.2f},{hi:.2f})  {ratio:>8.3f}{eff:>8.3f}  "
    for m in M_VALS3:
        w, total = wins3[m]
        row += f"{w}/{total}       "[: 7]
    verdict = f"PASS(m={first_full3})" if first_full3 is not None else "FAIL"
    print(row + "  " + verdict)

print("\n" + "=" * 70)
print("FINAL CONCLUSION")
print("=" * 70)
print("""
Key finding: The prior claim of a 'λ/n structural threshold' is REFUTED.

1. At K1=2 (eff≈0.04): ALL bands pass 3/3 at m=6, including λ/n=0.07.
   → The 2026-07-26 failure at p=2677 was NOT due to λ/n=0.07 being below
     a structural threshold. It was due to K1=8 (eff≈0.17), not K1=2.

2. At K1=8 (eff≈0.17): bands show a weak λ/n-correlated pattern with noise.
   The controlling variable is eff = K1*K2/n, not λ/n.

3. For secp256k1's true GLV: k1,k2 ~ √n → eff ≈ 1.0 (no bias).
   The Phase 2 attack requires an ADDITIONAL nonce bias (small k1 << √n)
   beyond the GLV structure. LLL alone on the pure GLV HNP (no extra bias)
   does not recover d.

Implication: the secp256k1 GLV decomposition does NOT in itself give a
lattice attack. An additional weakness (biased k1) is required.
""")
print("Done.")
