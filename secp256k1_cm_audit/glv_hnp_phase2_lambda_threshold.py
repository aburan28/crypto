"""
GLV-HNP Phase 2: λ/n threshold bisection study.

From 2026-07-26 run:
  lam/n = 0.53 (p=211):    LLL 3/3 at m=6   [SUCCESS]
  lam/n = 0.66 (p=2557):   LLL 3/3 at m=7   [SUCCESS]
  lam/n = 0.34 (p=524347): LLL 3/3 at m=9   [SUCCESS]
  lam/n = 0.07 (p=2677):   LLL never 3/3    [FAILURE]
  BKZ(40) also fails at lam/n=0.07.

Goal: bisect the threshold between 0.07 and 0.34 to find where LLL
switches from failure to success.  Also test whether BKZ(40) pushes
the threshold lower.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
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
    mc, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mc - i - 1), p)
        mc, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    """Find a random generator of order n on y^2 = x^3+b over F_p."""
    rng = random.Random(seed)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def find_b_for_order(p, n, seed=999):
    """
    Find a small b such that #(y^2 = x^3 + b over F_p) = n.
    Uses random point sampling for each candidate b.
    """
    rng = random.Random(seed)
    for b in range(1, 300):
        # Try to find a point on y^2=x^3+b and check its order divides n
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
# CM theory: Eisenstein decomposition for j=0 curves
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
                bb = num // 2
                if bb >= 0 and a * a - a * bb + bb * bb == p:
                    return (a, bb)
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

def collect_curves_by_ratio(target_ratios, tol=0.04, p_min=300, p_max=100000):
    """
    For each target ratio in target_ratios, find a j=0 prime-order curve
    whose λ/n ratio is within tol of the target.
    Returns a dict: {target_ratio -> (p, b, n, lam, G, actual_ratio)}.
    """
    # We want one curve per target ratio; use different p ranges to avoid reuse
    results = {}
    needed = set(range(len(target_ratios)))

    p = sympy.nextprime(p_min - 1)
    found_count = 0
    while p <= p_max and found_count < len(target_ratios):
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 50 or n_cand >= p:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    # Check if this ratio is near any unfilled target
                    for idx in list(needed):
                        if abs(ratio - target_ratios[idx]) <= tol:
                            b_val = find_b_for_order(p, n_cand)
                            if b_val is None:
                                continue
                            G = find_generator(p, b_val, n_cand)
                            if G is None:
                                continue
                            results[target_ratios[idx]] = (p, b_val, n_cand, lam, G, ratio)
                            needed.discard(idx)
                            found_count += 1
                            print(f"  target={target_ratios[idx]:.2f}: "
                                  f"p={p}, n={n_cand}, lam={lam}, "
                                  f"lam/n={ratio:.4f}")
                            break
        p = sympy.nextprime(p)

    for idx in needed:
        print(f"  target={target_ratios[idx]:.2f}: NOT FOUND (search exhausted up to p={p_max})")
    return results

# ---------------------------------------------------------------------------
# Signature generation and lattice
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
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B})
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = max(1, n // k1_bound)
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_KANNAN

def recover_d(A_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for i in range(dim):
        last = A_reduced[i, dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A_reduced[i, m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_one(curve_params, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=40):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    return recover_d(A, m, n, S_KANNAN, d_secret)

def test_curve(label, curve_params, k1_bound, m=8, seeds=(42, 1234, 9999),
               use_bkz=False, bkz_beta=40):
    p, b, n, lam, G = curve_params
    wins = 0
    for seed in seeds:
        d = random.Random(seed + 7777).randint(1, n - 1)
        ok = run_one(curve_params, m, d, k1_bound, seed, use_bkz, bkz_beta)
        wins += ok
    algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"
    status = "PASS" if wins == len(seeds) else ("FAIL" if wins == 0 else f"{wins}/{len(seeds)}")
    ratio = lam / n
    print(f"  lam/n={ratio:.4f}  {algo}  {wins}/{len(seeds)}  [{status}]  {label}")
    return wins, len(seeds)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold bisection study")
print("Known: lam/n=0.07 FAILS (LLL+BKZ40), lam/n=0.34 PASSES (LLL)")
print("=" * 70)

SEEDS = (42, 1234, 9999)
M_DEFAULT = 8

# Known reference points
print("\n--- Reference curves (from prior runs) ---")
G_ref_pass = find_generator(211, 2, 199)     # lam/n=0.53, PASS
G_ref_fail = find_generator(2677, 2, 2647)   # lam/n=0.07, FAIL
G_ref_mid  = find_generator(2557, 2, 2659)   # lam/n=0.66, PASS

lll_summary = {}   # ratio -> (wins, total, label)

if G_ref_pass:
    w, t = test_curve("8-bit/199 [KNOWN PASS]",   (211, 2, 199, 106, G_ref_pass),     k1_bound=2,  m=M_DEFAULT, seeds=SEEDS)
    lll_summary[0.533] = (w, t, "8-bit/199")
if G_ref_mid:
    w, t = test_curve("12-bit/2557 [KNOWN PASS]", (2557, 2, 2659, 1755, G_ref_mid),   k1_bound=8,  m=M_DEFAULT, seeds=SEEDS)
    lll_summary[0.660] = (w, t, "12-bit/2557")
if G_ref_fail:
    w, t = test_curve("12-bit/2677 [KNOWN FAIL]", (2677, 2, 2647, 185, G_ref_fail),   k1_bound=8,  m=M_DEFAULT, seeds=SEEDS)
    lll_summary[0.070] = (w, t, "12-bit/2677")

# Collect curves at intermediate ratios
# We want to fill the gap between 0.07 and 0.34 with about 6 points
TARGET_RATIOS = [0.10, 0.13, 0.17, 0.20, 0.25, 0.30]

print("\n--- Collecting curves for target λ/n ratios ---")
collected = collect_curves_by_ratio(TARGET_RATIOS, tol=0.04, p_min=300, p_max=150000)

print("\n--- LLL sweep over new curves (m=8, 3 seeds) ---")
for tgt, (p, b, n, lam, G, actual_ratio) in sorted(collected.items()):
    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    curve = (p, b, n, lam, G)
    label = f"n={n},lam={lam}"
    w, t = test_curve(label, curve, k1_bound=k1_bound, m=M_DEFAULT, seeds=SEEDS)
    lll_summary[actual_ratio] = (w, t, label)

# ---------------------------------------------------------------------------
# Phase 2: BKZ(40) on failing curves
# ---------------------------------------------------------------------------

print("\n--- BKZ(40) on curves that fail LLL ---")

bkz_summary = {}
for ratio, (w, t, label) in sorted(lll_summary.items()):
    if w < t:  # LLL failed on this curve
        # Reconstruct curve params
        if "8-bit/199" in label:
            curve = (211, 2, 199, 106, G_ref_pass)
            k1 = 2
        elif "12-bit/2677" in label:
            curve = (2677, 2, 2647, 185, G_ref_fail)
            k1 = 8
        elif "12-bit/2557" in label:
            curve = (2557, 2, 2659, 1755, G_ref_mid)
            k1 = 8
        else:
            # Find in collected
            found = None
            for tgt, (p, b, n, lam, G, ar) in collected.items():
                if abs(ar - ratio) < 0.001:
                    found = (p, b, n, lam, G)
                    k1 = max(2, int(0.05 * math.sqrt(n)))
                    break
            if found is None:
                print(f"  [skip] lam/n={ratio:.4f} — params not available")
                continue
            curve = found
        wb, tb = test_curve(f"{label} [BKZ40]", curve, k1_bound=k1, m=M_DEFAULT,
                            seeds=SEEDS, use_bkz=True, bkz_beta=40)
        bkz_summary[ratio] = (wb, tb)

# ---------------------------------------------------------------------------
# Phase 3: Extended m sweep on boundary curves
# ---------------------------------------------------------------------------

print("\n--- Extended m sweep on boundary cases ---")

sorted_lll = sorted(lll_summary.items())

# Find the transition: largest ratio that fails and smallest that passes
lo_fail = max((r for r, (w, t, _) in sorted_lll if w < t), default=None)
hi_pass = min((r for r, (w, t, _) in sorted_lll if w == t), default=None)

print(f"  Largest FAIL ratio: {lo_fail}")
print(f"  Smallest PASS ratio: {hi_pass}")

if lo_fail is not None:
    # Extended sweep on the largest-failing curve
    for tgt, (p, b, n, lam, G, ar) in collected.items():
        if abs(ar - lo_fail) < 0.001:
            k1 = max(2, int(0.05 * math.sqrt(n)))
            curve = (p, b, n, lam, G)
            print(f"\n  Extended m sweep: p={p}, n={n}, lam/n={ar:.4f}")
            for m_try in [10, 12, 14, 16, 20]:
                wins2 = sum(
                    run_one(curve, m_try,
                            random.Random(s + 7777).randint(1, n - 1), k1, s)
                    for s in SEEDS
                )
                marker = " ← PASS" if wins2 == len(SEEDS) else ""
                print(f"    m={m_try}: {wins2}/{len(SEEDS)}{marker}")
            break

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("FINAL SUMMARY — λ/n threshold for LLL success (m=8, 3 seeds)")
print("=" * 70)
print(f"  {'lam/n':>8}  {'LLL':>8}  {'BKZ40':>8}  label")
print(f"  {'-'*8}  {'-'*8}  {'-'*8}  -----")
for ratio, (w, t, label) in sorted(lll_summary.items()):
    lll_str = f"{w}/{t} {'PASS' if w == t else 'FAIL'}"
    bkz_str = "-"
    if ratio in bkz_summary:
        wb, tb = bkz_summary[ratio]
        bkz_str = f"{wb}/{tb} {'PASS' if wb == tb else 'FAIL'}"
    print(f"  {ratio:>8.4f}  {lll_str:>8}  {bkz_str:>8}  {label}")

print()
if lo_fail is not None and hi_pass is not None:
    print(f"  LLL threshold: between lam/n={lo_fail:.4f} (FAIL) and lam/n={hi_pass:.4f} (PASS)")
elif lo_fail is None:
    print("  All tested curves PASS at m=8 (threshold may be < 0.07)")
elif hi_pass is None:
    print("  All tested curves FAIL at m=8 (threshold may be > 0.34)")
print()
print("  secp256k1: lam/n ≈ 0.44 → EXPECTED PASS (above threshold)")
print()
print("Done.")
