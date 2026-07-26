"""
GLV-HNP Phase 2 — lambda/n threshold study.

Goal: determine the critical lam/n ratio below which LLL fails to recover d.

Prior data:
  lam/n = 0.07 (n=2647, lam=185):   LLL FAILS, BKZ(40) FAILS
  lam/n = 0.34 (n=523969):           LLL success at m=9
  lam/n = 0.53 (n=199, lam=106):     LLL success at m=4

Bisection target: what is the threshold in [0.07, 0.34]?

Protocol:
  Part 1 — bin sweep using 12-bit curves with lam_small/n in bins
           [0.07..0.10), [0.10..0.14), [0.14..0.18), [0.18..0.22),
           [0.22..0.28), [0.28..0.35), [0.35..0.45).
           For each bin, run LLL sweep m=3..16 (5 seeds). Report first 3/3 m.
  Part 2 — root-flip test: for the known-failure curve (p=2677, n=2647, lam=185,
           lam/n=0.07), swap to lam_large=2461 (ratio=0.93). Does LLL succeed?
  Part 3 — secp256k1 reference ratio: lambda_secp ~ 0.44*n; verify it's above threshold.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys

try:
    import sympy
    from fpylll import IntegerMatrix, LLL, BKZ
except ImportError as e:
    sys.exit(f"Missing dependency: {e}\n  pip install sympy fpylll")

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass a=0, y^2 = x^3 + b)
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
    m_var, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_var - i - 1), p)
        m_var, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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
# CM/Eisenstein utilities
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b_val = num // 2
                if b_val >= 0 and a * a - a * b_val + b_val * b_val == p:
                    return (a, b_val)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue_both(n):
    """Return (lam_small, lam_large) — both roots of x^2+x+1=0 mod n."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    # verify
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0
    assert (r2 * r2 + r2 + 1) % n == 0
    lam_small = min(r1, r2)
    lam_large = max(r1, r2)
    return lam_small, lam_large

# ---------------------------------------------------------------------------
# Curve search: find 12-bit j=0 curves with lam_small/n in a target range
# ---------------------------------------------------------------------------

def find_curve_with_ratio(lo, hi, bit_target=12, max_search=3000, seed_offset=0):
    """
    Find a j=0 prime-order curve with lam_small/n in [lo, hi).
    Returns (p, b_coeff, n, lam_small, lam_large, G) or None.
    """
    p_lo = max(2**(bit_target - 1), 200)
    p_hi = 2**bit_target
    p = sympy.nextprime(p_lo + seed_offset * 17)
    tried = 0
    while p < p_hi and tried < max_search:
        tried += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 8 or n_cand >= 2 * p:
                        continue
                    if not sympy.isprime(n_cand) or n_cand % 3 != 1:
                        continue
                    lam_small, lam_large = glv_eigenvalue_both(n_cand)
                    if lam_small is None:
                        continue
                    ratio = lam_small / n_cand
                    if lo <= ratio < hi:
                        # find b_coeff
                        for b_try in range(1, 100):
                            rng_tmp = random.Random(42 + b_try)
                            found_b = False
                            for _ in range(500):
                                x = rng_tmp.randint(0, p - 1)
                                rhs = (pow(x, 3, p) + b_try) % p
                                y = tonelli_shanks(rhs, p)
                                if y is not None and y != 0:
                                    P_test = (x, y)
                                    if ec_mul(P_test, n_cand, p) is None:
                                        found_b = True
                                        break
                            if found_b:
                                G = find_generator(p, b_try, n_cand)
                                if G is not None:
                                    return (p, b_try, n_cand, lam_small, lam_large, G)
                                break
        p = sympy.nextprime(p)
    return None

# ---------------------------------------------------------------------------
# Lattice construction and recovery
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

def recover_d_from_reduced(M_reduced, m, n, S_KANNAN, d_secret):
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

def run_single(curve_params_ext, lam_override, m, d_secret, k1_bound, seed=42,
               use_bkz=False, bkz_beta=20):
    """Run one LLL/BKZ experiment. lam_override allows testing non-default root."""
    p, b_coeff, n, lam_small, lam_large, G = curve_params_ext
    lam = lam_override
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b_coeff, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d_from_reduced(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep_curve_ext(label, curve_params_ext, lam_override, k1_bound,
                    m_range, seeds, use_bkz=False, bkz_beta=20):
    p, b_coeff, n, lam_small, lam_large, G = curve_params_ext
    lam = lam_override
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"
    ratio = lam / n
    print(f"  [{algo}] {label}  lam/n={ratio:.4f}  K1={k1_bound}  "
          f"eff={eff:.4f}  m_thresh≈{m_thresh:.1f}")
    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 91919).randint(1, n - 1)
            ok = run_single(curve_params_ext, lam, m, d_trial, k1_bound, seed,
                            use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
    # find first m with all seeds success
    first_all = None
    for m in sorted(results):
        wins, total = results[m]
        if wins == total:
            first_all = m
            break
    line = "  m: " + " ".join(f"{m}={w}/{t}" for m, (w,t) in sorted(results.items()))
    print(f"    {line}")
    if first_all:
        print(f"    -> 3/3 first at m={first_all}")
    else:
        print(f"    -> NEVER 3/3 in m={list(m_range)}")
    return results, first_all

# ---------------------------------------------------------------------------
# Ratio bins to study (using smaller root convention)
# ---------------------------------------------------------------------------

BINS = [
    (0.05, 0.10, "bin_A [0.05,0.10)"),
    (0.10, 0.14, "bin_B [0.10,0.14)"),
    (0.14, 0.18, "bin_C [0.14,0.18)"),
    (0.18, 0.23, "bin_D [0.18,0.23)"),
    (0.23, 0.28, "bin_E [0.23,0.28)"),
    (0.28, 0.34, "bin_F [0.28,0.34)"),
    (0.34, 0.42, "bin_G [0.34,0.42)"),
    (0.42, 0.50, "bin_H [0.42,0.50)"),
]

SEEDS = [42, 1234, 9999]
M_RANGE = range(3, 17)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — lambda/n threshold study")
print("=" * 70)
print(f"Searching for 12-bit j=0 curves, one per lam/n bin.")
print()

# ---- Part 1: bin sweep ------------------------------------------------------
print("=" * 70)
print("PART 1: bin sweep — LLL success vs lam/n (smaller root)")
print("=" * 70)

bin_results = []
for (lo, hi, bin_label) in BINS:
    print(f"\n--- {bin_label} ---")
    curve = find_curve_with_ratio(lo, hi, bit_target=12)
    if curve is None:
        # Try 11-bit
        curve = find_curve_with_ratio(lo, hi, bit_target=11)
    if curve is None:
        # Try 10-bit
        curve = find_curve_with_ratio(lo, hi, bit_target=10)
    if curve is None:
        print(f"  SKIP: no curve found in this bin at 10-12 bit range")
        bin_results.append((bin_label, lo, hi, None, None, None))
        continue

    p, b_coeff, n, lam_small, lam_large, G = curve
    k2_bound = math.isqrt(n) + 1
    # K1_BOUND: target eff ~ 0.06 (want ~4-6 signatures minimum)
    k1_bound = max(2, int(0.06 * math.sqrt(n)))
    eff = k1_bound * k2_bound / n

    print(f"  Curve: p={p}, b={b_coeff}, n={n} ({n.bit_length()}b), "
          f"lam_small={lam_small}, lam_large={lam_large}")
    print(f"  lam_small/n={lam_small/n:.4f}, lam_large/n={lam_large/n:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}")

    # LLL sweep with smaller root
    res_lll, first_lll = sweep_curve_ext(
        f"n={n} lam_small={lam_small}", curve, lam_small,
        k1_bound, M_RANGE, SEEDS, use_bkz=False
    )

    # If LLL fails, try BKZ(40)
    first_bkz = None
    if first_lll is None:
        res_bkz, first_bkz = sweep_curve_ext(
            f"n={n} lam_small={lam_small}", curve, lam_small,
            k1_bound, M_RANGE, SEEDS, use_bkz=True, bkz_beta=40
        )

    bin_results.append((bin_label, lo, hi, lam_small/n, first_lll, first_bkz))

# ---- Part 2: root-flip test -------------------------------------------------
print("\n" + "=" * 70)
print("PART 2: root-flip test on known-failure curve (p=2677, n=2647, lam_small=185)")
print("=" * 70)

p_fail, b_fail, n_fail = 2677, 2, 2647
lam_fail_small = 185   # ratio = 0.07 (KNOWN TO FAIL)
# Compute lam_large
lam_fail_large = (n_fail - 1 - lam_fail_small) % n_fail
print(f"  lam_small={lam_fail_small} (ratio={lam_fail_small/n_fail:.4f})  "
      f"lam_large={lam_fail_large} (ratio={lam_fail_large/n_fail:.4f})")

G_fail = find_generator(p_fail, b_fail, n_fail)
if G_fail is None:
    print("  ERROR: could not find generator")
else:
    curve_fail = (p_fail, b_fail, n_fail, lam_fail_small, lam_fail_large, G_fail)
    k2_fail = math.isqrt(n_fail) + 1
    k1_fail = 8  # same as prior experiments

    print(f"\n  Test A: lam_small={lam_fail_small} (ratio=0.07 — known failure)")
    res_a, first_a = sweep_curve_ext(
        f"p=2677 lam_small=185", curve_fail, lam_fail_small,
        k1_fail, range(5, 15), SEEDS, use_bkz=False
    )

    print(f"\n  Test B: lam_large={lam_fail_large} (ratio={lam_fail_large/n_fail:.4f} — ROOT FLIP)")
    res_b, first_b = sweep_curve_ext(
        f"p=2677 lam_large={lam_fail_large}", curve_fail, lam_fail_large,
        k1_fail, range(5, 15), SEEDS, use_bkz=False
    )

    print(f"\n  Comparison: small-root LLL -> {'fail' if first_a is None else f'ok m={first_a}'} | "
          f"large-root LLL -> {'fail' if first_b is None else f'ok m={first_b}'}")

# ---- Part 3: secp256k1 reference ratio --------------------------------------
print("\n" + "=" * 70)
print("PART 3: secp256k1 reference lambda/n ratio")
print("=" * 70)
# secp256k1: n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
# lambda_secp = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
n_secp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
lam_secp = 0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
ratio_secp = lam_secp / n_secp
print(f"  secp256k1 n (256b): {hex(n_secp)}")
print(f"  secp256k1 lambda:   {hex(lam_secp)}")
print(f"  lam/n = {ratio_secp:.6f}")
# The OTHER root
lam_secp_large = n_secp - 1 - lam_secp
ratio_secp_large = lam_secp_large / n_secp
print(f"  other root lam'/n = {ratio_secp_large:.6f}")
print()

# ---- Summary ----------------------------------------------------------------
print("=" * 70)
print("SUMMARY TABLE — LLL threshold vs lam/n (smaller root)")
print("=" * 70)
print(f"{'Bin':<22} {'lam/n':>7} {'LLL first-m':>12} {'BKZ(40) first-m':>16}")
print("-" * 60)
for (bin_label, lo, hi, ratio, first_lll, first_bkz) in bin_results:
    if ratio is None:
        print(f"  {bin_label:<20}   -       SKIP             SKIP")
        continue
    lll_str = str(first_lll) if first_lll else "FAIL"
    bkz_str = str(first_bkz) if first_bkz else ("-" if first_lll is not None else "FAIL")
    print(f"  {bin_label:<20} {ratio:7.4f}   {lll_str:^12} {bkz_str:^16}")

print()
print(f"secp256k1 reference: lam/n = {ratio_secp:.4f}  (other root: {ratio_secp_large:.4f})")
print()
print("Interpretation:")
for (bin_label, lo, hi, ratio, first_lll, first_bkz) in bin_results:
    if ratio is None: continue
    if first_lll is not None:
        print(f"  {bin_label}: lam/n={ratio:.4f} -> LLL SUCCESS (m={first_lll})")
    elif first_bkz is not None:
        print(f"  {bin_label}: lam/n={ratio:.4f} -> LLL fail, BKZ(40) SUCCESS (m={first_bkz})")
    else:
        print(f"  {bin_label}: lam/n={ratio:.4f} -> BOTH FAIL")

print("\nDone.")
