"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Sweeps λ/n ∈ {0.08, 0.10, 0.13, 0.17, 0.20, 0.25, 0.30} on 10–12 bit
j=0 prime-order curves with fixed K1_BOUND=8, K2=sqrt(n)+1.

Goal: find the critical ratio below which LLL fails, giving a concrete
threshold separating feasible from infeasible Phase 2 attacks.

Background (from 2026-07-26 run):
  lam/n = 0.53  (8-bit,  n=199):   LLL 3/3 at m=4   ✓
  lam/n = 0.66  (12-bit, n=2659):  LLL 3/3 at m=7   ✓
  lam/n = 0.34  (20-bit, n=524k):  LLL 3/3 at m=9   ✓
  lam/n = 0.07  (12-bit, n=2647):  LLL never 3/3    ✗ (BKZ(40) also fails)

For secp256k1: lam/n ≈ 0.44  (in the "safe" regime).

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0, y^2 = x^3 + b)
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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b2 = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b2 * b2 % p, t * b2 * b2 % p, r * b2 % p

def find_generator(p, b, n):
    """Find a generator of order n on y^2=x^3+b over F_p (n prime)."""
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
# CM theory: Eisenstein decomposition for j=0
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - a*b + b^2 = p, a,b > 0."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b > 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    """Six Frobenius traces for j=0 sextic twists given Eisenstein (a,b)."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Smaller GLV eigenvalue λ with λ²+λ+1≡0 (mod n). Returns (lam, lam/n)."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    lam = min(r1, n - 1 - r1)
    return lam, lam / n

# ---------------------------------------------------------------------------
# Curve finder: search for a j=0 GLV curve with lam/n in [lo, hi]
# ---------------------------------------------------------------------------

def find_curve_with_ratio(lo, hi, p_lo=500, p_hi=8000):
    """
    Search for a prime p in [p_lo, p_hi] with p ≡ 1 (mod 6), prime order n,
    and GLV eigenvalue lam/n ∈ [lo, hi]. Returns (p, b, n, lam, G) or None.
    """
    for p in sympy.primerange(p_lo, p_hi):
        if p % 6 != 1: continue
        eis = eisenstein_decompose(p)
        if eis is None: continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n_cand = p + 1 - t
            if n_cand < 100: continue
            if not sympy.isprime(n_cand): continue
            if n_cand % 3 != 1: continue
            lam, ratio = glv_eigenvalue(n_cand)
            if lam is None: continue
            if not (lo <= ratio <= hi): continue
            # Find b parameter (twist) for this n
            found_b = None
            for b_try in range(1, p):
                # Quick: check one point
                rng_tmp = random.Random(b_try)
                x = rng_tmp.randint(0, p - 1)
                rhs = (pow(x, 3, p) + b_try) % p
                y = tonelli_shanks(rhs, p)
                if y is None: continue
                P = (x, y)
                if ec_mul(P, n_cand, p) is None:
                    found_b = b_try
                    break
            if found_b is None: continue
            G = find_generator(p, found_b, n_cand)
            if G is None: continue
            return (p, found_b, n_cand, lam, G)
    return None

# ---------------------------------------------------------------------------
# Signature generation with k1 bias
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Build GLV Phase 2 lattice (column-diagonal scaling)
# ---------------------------------------------------------------------------

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

    return M, S_KANNAN

def recover_d(A_mat, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for i in range(dim):
        last = A_mat[i, dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A_mat[i, m]) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Single experiment: run LLL on one (curve, m, seed) triple
# ---------------------------------------------------------------------------

def run_lll(curve_params, m, d_secret, k1_bound, seed=42, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    return recover_d(A, m, n, S_KANNAN, d_secret)

# ---------------------------------------------------------------------------
# Sweep: m_range x seeds for one curve
# ---------------------------------------------------------------------------

def sweep(label, curve_params, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    if eff < 1.0:
        m_thresh_theory = math.ceil(math.log(n) / math.log(1.0 / eff))
    else:
        m_thresh_theory = 1
    algo = f"BKZ(β={bkz_beta})" if use_bkz else "LLL"
    print(f"\n  [{algo}] {label}")
    print(f"   p={p}, n={n} ({n.bit_length()}b), λ={lam}, λ/n={lam/n:.4f}")
    print(f"   K1={k1_bound}, K2={k2_bound}, eff={eff:.5f}, m_theory≈{m_thresh_theory}")
    first_success = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed * 31337 + 7).randint(1, n - 1)
            if run_lll(curve_params, m, d_trial, k1_bound, seed,
                       use_bkz=use_bkz, bkz_beta=bkz_beta):
                wins += 1
        total = len(seeds)
        marker = " ← 3/3!" if wins == total else ""
        print(f"   m={m}: {wins}/{total}{marker}")
        if wins == total and first_success is None:
            first_success = m
    return first_success  # None = never 3/3

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 68)
print("GLV-HNP Phase 2 — λ/n threshold bisection")
print("=" * 68)
print()
print("Known: λ/n=0.07 → LLL fails; λ/n=0.34 → LLL 3/3 at m=9.")
print("Goal:  find the critical threshold between these.")
print()

# K1 scales with n to keep eff = K1*K2/n ≈ sqrt(n)/300 (small bias fraction).
# Fixed K1=8 is only valid for 12-bit curves (n~2647); for smaller n it makes eff>0.5.
# Formula: K1 = max(2, n // 300).  For 9-bit n≈350: K1=2, eff≈0.11.
#                                   For 12-bit n≈2650: K1=8, eff≈0.16.
def k1_for_n(n): return max(2, n // 300)

SEEDS = [42, 1234, 9999]
M_RANGE = range(4, 15)

# Reference anchor curves (known from 2026-07-26 run)
ANCHORS = [
    # (label, p, b, n, lam) — already verified
    ("REF-low  [FAIL] lam/n≈0.07", 2677, 2, 2647, 185),
    ("REF-high [PASS] lam/n≈0.53", 211,  2, 199,  106),
]

# Target λ/n bands for bisection
# Each tuple: (center, lo, hi)
BANDS = [
    (0.08,  0.07,  0.10),
    (0.10,  0.09,  0.12),
    (0.13,  0.12,  0.15),
    (0.17,  0.15,  0.19),
    (0.20,  0.19,  0.22),
    (0.25,  0.23,  0.28),
    (0.30,  0.28,  0.33),
]

# Step 1: rebuild anchor curves
print("─" * 68)
print("Step 1: Anchor verification (known results)")
print("─" * 68)
anchor_params = []
for label, p_a, b_a, n_a, lam_a in ANCHORS:
    G_a = find_generator(p_a, b_a, n_a)
    if G_a is None:
        print(f"  {label}: ERROR — generator not found")
        anchor_params.append(None)
    else:
        anchor_params.append((p_a, b_a, n_a, lam_a, G_a))

for (label, *_), params in zip(ANCHORS, anchor_params):
    if params is None: continue
    p_a, b_a, n_a, lam_a, G_a = params
    K1_a = k1_for_n(n_a)
    sweep(label, params, K1_a, range(4, 12), SEEDS[:2])

# Step 2: find curves at each target band and sweep them
print()
print("─" * 68)
print("Step 2: Threshold bisection curves")
print("─" * 68)

summary = []  # (center_ratio, lam/n, first_success_m or None)

for (center, lo, hi) in BANDS:
    print(f"\n>>> Searching for curve with λ/n ∈ [{lo:.2f}, {hi:.2f}] (target {center:.2f})")
    # Prefer 12-bit curves (p in [1000, 6000]) so K1=max(2,n//300) ≈ 7-20
    curve = find_curve_with_ratio(lo, hi, p_lo=1000, p_hi=6000)
    if curve is None:
        print(f"  NOT FOUND at 12-bit — trying wider range")
        curve = find_curve_with_ratio(lo * 0.9, hi * 1.1, p_lo=300, p_hi=30000)
    if curve is None:
        print(f"  SKIP: no qualifying curve found.")
        summary.append((center, None, None))
        continue
    p_c, b_c, n_c, lam_c, G_c = curve
    ratio = lam_c / n_c
    K1_c = k1_for_n(n_c)
    m_first = sweep(f"λ/n={ratio:.4f} (K1={K1_c})", curve, K1_c, M_RANGE, SEEDS)
    summary.append((center, ratio, m_first))

# Step 3: BKZ rescue for near-threshold curves (those where LLL fails)
print()
print("─" * 68)
print("Step 3: BKZ(β=20) rescue for LLL-failing curves")
print("─" * 68)

failing = [(center, r, mf) for (center, r, mf) in summary
           if r is not None and mf is None and r < 0.20]
if not failing:
    print("  No LLL-failing curves in bisection range (or none found).")
else:
    for (center, ratio, _) in failing:
        # Re-find the same curve
        lo, hi = ratio * 0.98, ratio * 1.02
        curve = find_curve_with_ratio(lo, hi, p_lo=200, p_hi=50000)
        if curve is None: continue
        K1_c = k1_for_n(curve[2])
        mf_bkz = sweep(f"BKZ(20) λ/n={ratio:.4f} (K1={K1_c})", curve,
                       K1_c, range(5, 15), SEEDS, use_bkz=True, bkz_beta=20)
        print(f"  BKZ(20): first_success_m={mf_bkz}")

# Step 4: Summary table
print()
print("─" * 68)
print("SUMMARY TABLE")
print("─" * 68)
print(f"{'Target λ/n':>12}  {'Actual λ/n':>12}  {'LLL 3/3 at m':>14}  {'Status':}")
print("-" * 60)

# Add anchor rows
print(f"  0.07 (ref) →    0.0699       never 3/3    FAIL (known)")
print(f"  0.34 (ref) →    0.3400       m=9          PASS (known, 20b)")
print(f"  0.53 (ref) →    0.5327       m=4          PASS (known, 8b)")

for (center, ratio, m_first) in summary:
    if ratio is None:
        status = "NOT FOUND"
        ratio_str = "—"
        m_str = "—"
    elif m_first is None:
        status = "FAIL"
        ratio_str = f"{ratio:.4f}"
        m_str = "never"
    else:
        status = "PASS"
        ratio_str = f"{ratio:.4f}"
        m_str = str(m_first)
    print(f"  {center:.2f}          {ratio_str:>10}  {m_str:>14}  {status}")

print()
print("Interpretation:")
print("  PASS = LLL recovers d (3/3 seeds) at some m in range [4,13].")
print("  FAIL = LLL never recovers d at m ≤ 13 for any seed.")
print()
passes = [(r, m) for (c, r, m) in summary if r is not None and m is not None]
fails = [(r, m) for (c, r, m) in summary if r is not None and m is None]
if passes and fails:
    max_fail = max(r for r, _ in fails)
    min_pass = min(r for r, _ in passes)
    print(f"  λ/n threshold: between {max_fail:.4f} (highest FAIL) and {min_pass:.4f} (lowest PASS)")
elif passes:
    print(f"  All tested ratios PASS. Threshold is below {min(r for r, _ in passes):.4f}.")
elif fails:
    print(f"  All tested ratios FAIL. Threshold is above {max(r for r, _ in fails):.4f}.")
else:
    print("  Insufficient data.")

print()
print("Done.")
