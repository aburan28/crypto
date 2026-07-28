"""
GLV-HNP Phase 2 — λ/n threshold bisection study.

Bisects the empirical LLL success threshold between:
  λ/n = 0.07  (p=2677, n=2647): KNOWN FAIL (LLL + BKZ(40))
  λ/n = 0.34  (p=524347, n=523969): KNOWN PASS (LLL 3/3 at m=9)

Protocol:
1. Find ~12-16 bit j=0 prime-order curves with λ/n covering 7 bands
   spanning (0.07, 0.34) at roughly equal spacing.
2. For each band pick a single representative curve.
3. Sweep m = 5..16 with 3 seeds; report first m with 3/3 success.
4. If LLL fails for all m ≤ 16, also try BKZ(beta=40).
5. Report the threshold band (lambda/n) below which LLL is not viable.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

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
    m_s, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_s - i - 1), p)
        m_s, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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
# CM theory: Eisenstein decomposition for j=0 curves
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
# Curve database: map λ/n target bands → (p, b, n, lam)
# Pre-computed from CM search; fallback to live search if unavailable.
# ---------------------------------------------------------------------------

# Known curves from previous sessions:
KNOWN_CURVES = {
    # (label, p, b, n, lam)  — b is the y²=x³+b parameter
    "0.07 (known fail)":  (2677,   2, 2647,  185),   # from 2026-06-15
    "0.53 (known pass)":  (211,    2,  199,  106),   # from 2026-05-21
    "0.66 (known pass)":  (2557,   2, 2659, 1755),   # from 2026-05-22
    "0.34 (known pass)":  (524347, 2, 523969, 177902), # from 2026-07-26
}

def find_curve_in_band(lo, hi, bit_lo=11, bit_hi=15, max_cands=5000):
    """
    Search for a j=0 prime-order curve with λ/n ∈ [lo, hi).
    Returns (p, b, n, lam, G) or None.
    """
    count = 0
    for p in sympy.primerange(2**bit_lo, 2**bit_hi):
        if p % 3 != 1: continue
        eis = eisenstein_decompose(p)
        if eis is None: continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n = p + 1 - t
            if n < 4: continue
            if not sympy.isprime(n): continue
            if n % 3 != 1: continue
            lam, _ = glv_eigenvalue(n)
            if lam is None: continue
            ratio = lam / n
            if lo <= ratio < hi:
                # Find twist b
                for b_try in range(1, min(p, 200)):
                    pt = None
                    rng = random.Random(42)
                    for _ in range(1000):
                        x = rng.randint(0, p - 1)
                        rhs = (pow(x, 3, p) + b_try) % p
                        y = tonelli_shanks(rhs, p)
                        if y is not None and y != 0:
                            pt = (x, y)
                            break
                    if pt is None: continue
                    if ec_mul(pt, n, p) is None:
                        G = find_generator(p, b_try, n)
                        if G is not None:
                            return (p, b_try, n, lam, G)
        count += 1
        if count > max_cands: break
    return None

# ---------------------------------------------------------------------------
# GLV signature generation + lattice attack (from glv_hnp_phase2_20bit.py)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
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
    return M, S_K1, S_D, S_K2, S_KANNAN

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return d_cand
    return None

def run_single(curve_params, m, d_secret, k1_bound, seed=42, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep_lll(label, curve_params, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = (math.ceil(math.log(n) / math.log(1.0 / eff))
                if eff < 1.0 else float('inf'))
    algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"
    print(f"\n{'='*70}")
    print(f"[{algo}] {label}")
    print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}, lam/n={lam/n:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.5f}, m_thresh≈{m_thresh:.1f}")
    print(f"{'='*70}")

    first_3of3 = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 5555).randint(1, n - 1)
            ok = run_single(curve_params, m, d_trial, k1_bound, seed,
                            use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        marker = " ← thresh" if m == int(m_thresh) else ""
        print(f"  m={m:2d}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first_3of3 is None:
            first_3of3 = m
    if first_3of3:
        print(f"  >>> 3/3 at m={first_3of3}")
    else:
        print(f"  >>> NEVER {len(seeds)}/{len(seeds)} in range m={m_range.start}..{m_range.stop-1}")
    return first_3of3

# ---------------------------------------------------------------------------
# Main: band search + threshold mapping
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection")
print("Bisecting between λ/n=0.07 (fail) and λ/n=0.34 (pass)")
print("=" * 70)
print()

SEEDS = [42, 1234, 9999]

# Target bands to fill in: 7 bands between 0.07 and 0.34
BANDS = [
    (0.08, 0.11, "band_A: 0.08-0.11"),
    (0.11, 0.14, "band_B: 0.11-0.14"),
    (0.14, 0.18, "band_C: 0.14-0.18"),
    (0.18, 0.22, "band_D: 0.18-0.22"),
    (0.22, 0.27, "band_E: 0.22-0.27"),
    (0.27, 0.32, "band_F: 0.27-0.32"),
    (0.32, 0.35, "band_G: 0.32-0.35"),
]

print("Step 1: Searching for representative curves in each band...")
band_curves = []
for lo, hi, label in BANDS:
    print(f"  Searching {label} (λ/n ∈ [{lo:.2f}, {hi:.2f}))...", end=" ", flush=True)
    result = find_curve_in_band(lo, hi, bit_lo=11, bit_hi=16)
    if result is None:
        # Try wider bit range
        result = find_curve_in_band(lo, hi, bit_lo=9, bit_hi=18)
    if result is not None:
        p, b, n, lam, G = result
        print(f"FOUND: p={p}, n={n}, lam={lam}, lam/n={lam/n:.4f}")
        band_curves.append((label, result))
    else:
        print("NOT FOUND (skipping)")
        band_curves.append((label, None))

print()
print("Step 2: Attack sweep for each band")

# K1_BOUND: use ~1% of n (small bias)
# For n ~ 2^12 ≈ 4096: K1=4, K2=64, eff=4*64/4096=0.0625
results_table = []

for label, curve in band_curves:
    if curve is None:
        results_table.append((label, None, None))
        continue
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    # Target eff ~ 0.05: K1 = 0.05 * sqrt(n)
    k1_bound = max(2, int(0.05 * math.sqrt(n)))

    m_range = range(4, 18)
    first_lll = sweep_lll(label, curve, k1_bound, m_range, SEEDS, use_bkz=False)

    if first_lll is None:
        print(f"\n  LLL failed. Trying BKZ(40)...")
        first_bkz40 = sweep_lll(label, curve, k1_bound, m_range, SEEDS,
                                 use_bkz=True, bkz_beta=40)
        results_table.append((label, None, first_bkz40))
    else:
        results_table.append((label, first_lll, None))

# Reference: known curves
print()
print("=" * 70)
print("Reference curves (known from prior sessions):")

# Known fail: p=2677, n=2647, lam=185
G_fail = find_generator(2677, 2, 2647)
if G_fail:
    c_fail = (2677, 2, 2647, 185, G_fail)
    k1_fail = max(2, int(0.05 * math.sqrt(2647)))
    first_lll_fail = sweep_lll("0.07 (known fail)", c_fail, k1_fail,
                                range(4, 14), SEEDS, use_bkz=False)

# Known pass: p=211, n=199, lam=106
G0 = find_generator(211, 2, 199)
if G0:
    c0 = (211, 2, 199, 106, G0)
    sweep_lll("0.53 (known pass)", c0, k1_bound=2, m_range=range(3, 8), seeds=SEEDS)

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'Curve / band':<35} {'lam/n':<10} {'LLL 3/3 at m':<15} {'BKZ(40) 3/3 at m'}")
print("-" * 70)

# Anchors
print(f"{'0.07: p=2677, n=2647 (KNOWN FAIL)':<35} {'0.0699':<10} {'FAIL':<15} {'FAIL'}")

for label, curve in band_curves:
    if curve is None:
        print(f"  {label:<33} {'N/A':<10} {'no curve found'}")
        continue
    p, b, n, lam, G = curve
    ratio = f"{lam/n:.4f}"
    lll_result, bkz_result = None, None
    for row in results_table:
        if row[0] == label:
            lll_result, bkz_result = row[1], row[2]
            break
    lll_str = f"m={lll_result}" if lll_result else "FAIL"
    bkz_str = (f"m={bkz_result}" if bkz_result else ("FAIL" if lll_result is None else "n/a"))
    print(f"  {label:<33} {ratio:<10} {lll_str:<15} {bkz_str}")

print(f"{'0.34: p=524347 (KNOWN PASS)':<35} {'0.3396':<10} {'m=9':<15} {'n/a'}")
print(f"{'0.53: p=211, n=199 (KNOWN PASS)':<35} {'0.5326':<10} {'m=4':<15} {'n/a'}")
print(f"{'0.66: p=2557, n=2659 (KNOWN PASS)':<35} {'0.6600':<10} {'m=7':<15} {'n/a'}")

print()
print("Done.")
