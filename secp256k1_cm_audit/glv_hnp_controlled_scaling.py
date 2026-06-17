"""
GLV-HNP Phase 2: Controlled scaling experiment (2026-06-17).

Goal: isolate the effect of BIT SIZE on m/m_thresh, by fixing eff ≈ 0.15
across 8-bit, 12-bit, and 20-bit j=0 GLV curves.

2026-06-16 finding: 20-bit at eff=0.05 required m/m_thresh=1.80, vs 1.33/1.40
for 8/12-bit at eff=0.15.  Confound: eff was different.

This run: all three curves at eff ≈ 0.15, so any remaining m/m_thresh growth
is attributable purely to increasing lattice dimension.

Run: python3 glv_hnp_controlled_scaling.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic
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
    m_val, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_val - i - 1), p)
        m_val, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(100000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory for j=0 curves
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
    assert (r1 * r1 + r1 + 1) % n == 0
    lam = min(r1, r2)
    return lam, n - 1 - lam

def find_curve_for_order(p_start, target_min_bits, lam_ratio_lo=0.25, lam_ratio_hi=0.75):
    """Find a j=0 curve with prime order near p_start."""
    p = sympy.nextprime(p_start - 1)
    p_max = p_start * 2
    while p < p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam, _ = glv_eigenvalue(n_cand)
                        if lam is None:
                            continue
                        ratio = lam / n_cand
                        if lam_ratio_lo <= ratio <= lam_ratio_hi:
                            for b_try in range(1, 300):
                                pt = None
                                rng_tmp = random.Random(42)
                                for _ in range(500):
                                    x = rng_tmp.randint(0, p - 1)
                                    rhs = (pow(x, 3, p) + b_try) % p
                                    y = tonelli_shanks(rhs, p)
                                    if y is not None and y != 0:
                                        if ec_mul((x, y), n_cand, p) is None:
                                            pt = (x, y)
                                            break
                                if pt is not None:
                                    G = find_generator(p, b_try, n_cand)
                                    if G:
                                        return (p, b_try, n_cand, lam, G)
        p = sympy.nextprime(p)
    return None

# ---------------------------------------------------------------------------
# Signature generation (nonce bias: k = k1 + lam*k2 with k1<K1, k2<K2)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Build GLV lattice (column-diagonal scaling, from glv_hnp_phase2_20bit.py)
# ---------------------------------------------------------------------------

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
    M[m][m] = 1  # S_D=1 (d is full-size, no scaling)
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
        if abs(row[dim - 1]) != S_KANNAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Single experiment (one (m, seed) pair)
# ---------------------------------------------------------------------------

def run_experiment(curve_params, m, d_secret, k1_bound, seed=42):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

# ---------------------------------------------------------------------------
# Sweep (fixed eff, sweep m)
# ---------------------------------------------------------------------------

def sweep(label, curve_params, k1_bound, m_range, seeds):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / (-math.log(eff)))

    print(f"\n{'='*60}")
    print(f"Curve: {label}")
    print(f"  n={n} ({n.bit_length()}b), lam/n={lam/n:.4f}")
    print(f"  K1={k1_bound}, K2={k2_bound}, eff={eff:.4f}, m_thresh={m_thresh}")
    print(f"{'='*60}")

    first_full = None
    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 9001).randint(1, n - 1)
            ok = run_experiment(curve_params, m, d, k1_bound, seed)
            wins += ok
        results[m] = (wins, len(seeds))
        marker = " ← m_thresh" if m == m_thresh else ""
        print(f"  m={m:2d}: {wins}/{len(seeds)}{marker}")
        if wins == len(seeds) and first_full is None:
            first_full = m

    ratio = first_full / m_thresh if first_full else None
    print(f"  First 3/3: m={first_full}  →  ratio m/m_thresh = "
          f"{ratio:.2f}" if ratio else "  First 3/3: NOT FOUND in sweep range")
    return results, first_full, m_thresh, eff

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 60)
print("GLV-HNP Controlled Scaling: fix eff≈0.15, vary bit size")
print("=" * 60)

SEEDS = [42, 1234, 9999]

# ---- Curve 1: 8-bit (p=211, n=199, lam=106) --------------------------------
print("\n[1] Building 8-bit curve generator...")
G_8 = find_generator(211, 2, 199)
curve_8 = (211, 2, 199, 106, G_8)
# K1=2: eff = 2*(isqrt(199)+1)/199 = 2*15/199 = 0.1508
res_8, m1_8, mt_8, eff_8 = sweep(
    "8-bit: y²=x³+2 / F_211, n=199, lam=106",
    curve_8, k1_bound=2, m_range=range(3, 9), seeds=SEEDS
)

# ---- Curve 2: 12-bit (p=2557, n=2659, lam=1755) ----------------------------
print("\n[2] Building 12-bit curve generator...")
G_12 = find_generator(2557, 2, 2659)
curve_12 = (2557, 2, 2659, 1755, G_12)
# K1=8: eff = 8*(isqrt(2659)+1)/2659 = 8*52/2659 = 0.1564
res_12, m1_12, mt_12, eff_12 = sweep(
    "12-bit: y²=x³+2 / F_2557, n=2659, lam=1755",
    curve_12, k1_bound=8, m_range=range(4, 13), seeds=SEEDS
)

# ---- Curve 3: 20-bit (find fresh) ------------------------------------------
print("\n[3] Finding 20-bit curve (targeting eff≈0.15)...")
curve_20 = find_curve_for_order(2**19, 19)
if curve_20 is None:
    print("ERROR: No 20-bit curve found — aborting.")
else:
    p20, b20, n20, lam20, G20 = curve_20
    k2_20 = math.isqrt(n20) + 1
    # K1 = round(0.15 * n / K2) = round(0.15 * sqrt(n))
    k1_20 = max(2, round(0.15 * n20 / k2_20))
    eff_20_actual = k1_20 * k2_20 / n20
    print(f"  Found: p={p20}, b={b20}, n={n20} ({n20.bit_length()}b), lam={lam20}")
    print(f"  K1={k1_20}, K2={k2_20}, eff={eff_20_actual:.4f}")

    # m_thresh for 20-bit at eff≈0.15
    mt_20_est = math.ceil(math.log(n20) / (-math.log(eff_20_actual)))
    m_lo = max(3, mt_20_est - 2)
    m_hi = mt_20_est + 10

    res_20, m1_20, mt_20, eff_20 = sweep(
        f"20-bit: y²=x³+{b20} / F_{p20}, n={n20}, lam={lam20}",
        curve_20, k1_bound=k1_20, m_range=range(m_lo, m_hi), seeds=SEEDS
    )

# ---- Summary ----------------------------------------------------------------
print("\n" + "=" * 60)
print("CONTROLLED SCALING SUMMARY (eff ≈ 0.15 for all)")
print("=" * 60)
print(f"{'Curve':<28} {'n_bits':>6} {'eff':>7} {'m_thresh':>8} {'first_3/3':>9} {'ratio':>6}")
print("-" * 60)

for label, n_bits, eff_v, mt, m1 in [
    ("8-bit/199 (lam/n=0.53)", 8, eff_8, mt_8, m1_8),
    ("12-bit/2659 (lam/n=0.66)", 12, eff_12, mt_12, m1_12),
]:
    ratio_str = f"{m1/mt:.2f}" if m1 else "N/A"
    m1_str = str(m1) if m1 else "—"
    print(f"  {label:<26} {n_bits:>6} {eff_v:>7.4f} {mt:>8} {m1_str:>9} {ratio_str:>6}")

if curve_20 is not None:
    lam_ratio_20 = lam20 / n20
    ratio_str = f"{m1_20/mt_20:.2f}" if m1_20 else "N/A"
    m1_str = str(m1_20) if m1_20 else "—"
    label_20 = f"20-bit/{n20} (lam/n={lam_ratio_20:.2f})"
    print(f"  {label_20:<26} {n20.bit_length():>6} {eff_20:>7.4f} {mt_20:>8} {m1_str:>9} {ratio_str:>6}")

print()
print("If ratio m/m_thresh is stable (all ~1.3–1.5), LLL scales well with n.")
print("If ratio grows (>1.8 for 20-bit), LLL degrades as lattice dimension grows.")
print("Done.")
