"""
GLV-HNP Phase 2: Structural obstruction threshold mapping.

Objective: find the critical δ*(3λ,n)/n that separates "obstructed" (LLL fails
at K1=72) from "unobstructed" (LLL succeeds).  We know:
  - Curve A: δ/n ≈ 0.019 → obstructed (BKZ-40 fails, 0/3 at all m)
  - Curve C: δ/n ≈ 0.140 → unobstructed (LLL 3/3 at m=14)
  - secp256k1: δ/n ≈ 0.023 → expected obstructed

Method: scan 20-bit (n≈2^20) j=0 CM curves, bin by δ/n into five windows,
pick 2 representatives per window, run K1=72 LLL at m=10..20 (3 seeds each).

Run: python3 glv_hnp_delta_threshold.py
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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b_param, n):
    rng = random.Random(12345)
    for _ in range(50000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_param) % p
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
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b2 = num // 2
                if b2 >= 0 and a * a - a * b2 + b2 * b2 == p:
                    return (a, b2)
    return None

def j0_traces(a, b2):
    return [2*a - b2, -2*a + b2, -(a + b2), a + b2, 2*b2 - a, a - 2*b2]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def delta_ratio(lam, n):
    """δ(3λ,n)/n = min(3λ mod n, n - 3λ mod n) / n"""
    x = (3 * lam) % n
    return min(x, n - x) / n

def find_b_for_n(p, n):
    """Find curve parameter b such that #E(F_p): y²=x³+b = n."""
    for b_try in range(1, min(p, 500)):
        rng2 = random.Random(99)
        for _ in range(20):
            x = rng2.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b_try
                break
    return None

# ---------------------------------------------------------------------------
# Curve search by δ/n target range
# ---------------------------------------------------------------------------

DELTA_BINS = [
    ("bin0_tiny",   0.015, 0.035),   # ~secp256k1 region (δ/n≈0.023)
    ("bin1_low",    0.040, 0.065),   # bridge region
    ("bin2_midlow", 0.065, 0.095),   # mid-low
    ("bin3_mid",    0.095, 0.125),   # mid
    ("bin4_high",   0.130, 0.165),   # ~Curve C region (δ/n≈0.140)
]

def scan_curves(p_lo=2**19, p_hi=2**20, max_per_bin=2):
    """
    Scan 20-bit j=0 CM primes.  For each prime group order n, compute δ/n
    and slot into DELTA_BINS (up to max_per_bin per bin).
    Returns dict: bin_label -> list of (p, b, n, lam, delta_ratio).
    """
    bins = {label: [] for label, _, _ in DELTA_BINS}
    p = sympy.nextprime(p_lo - 1)
    scanned = 0
    while p < p_hi:
        scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    dr = delta_ratio(lam, n_cand)
                    for label, lo, hi in DELTA_BINS:
                        if lo <= dr < hi and len(bins[label]) < max_per_bin:
                            b_param = find_b_for_n(p, n_cand)
                            if b_param is None:
                                continue
                            G = find_generator(p, b_param, n_cand)
                            if G is None:
                                continue
                            bins[label].append((p, b_param, n_cand, lam, dr))
                            print(f"  [{label}] p={p}, n={n_cand} ({n_cand.bit_length()}b), "
                                  f"lam={lam}, δ/n={dr:.4f}")
                            break
        p = sympy.nextprime(p)
        # Stop early if all bins are full
        if all(len(bins[label]) >= max_per_bin for label, _, _ in DELTA_BINS):
            break

    print(f"  (scanned {scanned} primes in [{p_lo}, {p_hi}])")
    return bins

# ---------------------------------------------------------------------------
# GLV signatures + lattice attack
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
        attempts += 1
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(0, K2 - 1)
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
        sigs.append({'A': A, 'B': B})
    return d_secret, sigs

def build_lattice(sigs, n, lam, K1):
    m = len(sigs)
    K2 = math.isqrt(n) + 1
    dim = 2 * m + 2
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1            # d row
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_lll(p, b_param, n, lam, G, m, K1, seed):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, K1)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep(p, b_param, n, lam, G, K1, m_range, seeds):
    """Run LLL sweep; return dict m -> (wins, total)."""
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_lll(p, b_param, n, lam, G, m, K1, s))
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Reference curves from prior runs
# ---------------------------------------------------------------------------

# Curve A (from 2026-06-20): p=524347, n≈520000, lam=177902, δ/n=0.019
# Curve C (from 2026-06-20): p=624517, n=622957, lam=178615, δ/n=0.140
# These are used as anchors; we verify them here.

REFERENCE_CURVES = [
    # (label, p, b, n, lam, expected_delta, expected_outcome)
    ("CurveA_anchor", 524347, None, 522243, 177902, 0.019, "obstructed"),
    ("CurveC_anchor", 624517, None, 622957, 178615, 0.140, "unobstructed"),
]

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: Structural obstruction threshold mapping")
print("=" * 70)
print()

K1 = 72         # Same as Curve A / Curve C experiments
SEEDS = [42, 1234, 9999]
M_RANGE = range(10, 21)

# ---- Step 1: Verify reference curves ----------------------------------------
print("Step 1: Verify reference anchor curves")
print("-" * 40)

ref_results = {}
for label, p, b, n, lam, exp_delta, exp_outcome in REFERENCE_CURVES:
    K2 = math.isqrt(n) + 1
    eff = K1 * K2 / n
    dr = delta_ratio(lam, n)
    print(f"\n{label}: p={p}, n={n}, lam={lam}")
    print(f"  δ/n={dr:.4f} (expected ~{exp_delta}), eff≈{eff:.4f}")

    # Find b if not given
    if b is None:
        b = find_b_for_n(p, n)
        if b is None:
            print(f"  ERROR: could not find b for n={n}; skipping anchor")
            ref_results[label] = None
            continue
        print(f"  Found b={b}")

    G = find_generator(p, b, n)
    if G is None:
        print(f"  ERROR: could not find generator; skipping anchor")
        ref_results[label] = None
        continue

    print(f"  Running K1={K1}, m=10..20, seeds×3 ...")
    res = sweep(p, b, n, lam, G, K1, M_RANGE, SEEDS)
    first_win = next((m for m, (w, t) in sorted(res.items()) if w == t), None)
    if first_win:
        print(f"  RESULT: 3/3 first at m={first_win}  [expected: {exp_outcome}]")
    else:
        max_wins = max(w for w, t in res.values())
        print(f"  RESULT: never 3/3 (max {max_wins}/3)  [expected: {exp_outcome}]")
    ref_results[label] = res

# ---- Step 2: Scan for new curves in δ-bins ----------------------------------
print("\n" + "=" * 70)
print("Step 2: Scan 20-bit j=0 CM curves by δ/n bin")
print("-" * 40)

bins = scan_curves(p_lo=2**19, p_hi=2**20, max_per_bin=2)

# ---- Step 3: Run threshold experiment ---------------------------------------
print("\n" + "=" * 70)
print("Step 3: LLL sweep for each δ-bin representative")
print("-" * 40)

threshold_results = {}  # bin_label -> list of (delta_ratio, first_win_m_or_None)

for label, lo, hi in DELTA_BINS:
    curves = bins[label]
    print(f"\n[{label}]  δ/n ∈ [{lo:.3f}, {hi:.3f})")
    if not curves:
        print("  No curves found in this bin.")
        threshold_results[label] = []
        continue

    for (p, b, n, lam, dr) in curves:
        K2 = math.isqrt(n) + 1
        eff = K1 * K2 / n
        print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}, δ/n={dr:.4f}, eff={eff:.4f}")

        G = find_generator(p, b, n)
        if G is None:
            print(f"  ERROR: generator not found; skip")
            if label not in threshold_results:
                threshold_results[label] = []
            threshold_results[label].append((dr, None))
            continue

        res = sweep(p, b, n, lam, G, K1, M_RANGE, SEEDS)
        first_win = next((m for m, (w, t) in sorted(res.items()) if w == t), None)
        detail = " ".join(f"m{m}={w}/{t}" for m, (w, t) in sorted(res.items()))
        if first_win:
            print(f"  3/3 first at m={first_win}   detail: {detail}")
        else:
            max_wins = max(w for w, t in res.values())
            print(f"  never 3/3 (max {max_wins}/3)   detail: {detail}")

        if label not in threshold_results:
            threshold_results[label] = []
        threshold_results[label].append((dr, first_win))

# ---- Step 4: Summary table --------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print(f"{'δ/n range':<20} {'δ/n':<8} {'K1=72, m=10..20':<30} outcome")
print("-" * 70)

# Anchors first
for label, p, b, n, lam, exp_delta, exp_outcome in REFERENCE_CURVES:
    if ref_results.get(label) is not None:
        res = ref_results[label]
        dr = delta_ratio(lam, n)
        first_win = next((m for m, (w, t) in sorted(res.items()) if w == t), None)
        outcome = f"3/3 at m={first_win}" if first_win else "never 3/3"
        print(f"  {label:<20} {dr:.4f}   {outcome:<30} [anchor {exp_outcome}]")

for label, lo, hi in DELTA_BINS:
    for dr, fw in threshold_results.get(label, []):
        outcome = f"3/3 at m={fw}" if fw else "never 3/3"
        print(f"  [{label}]  {f'{lo:.3f}–{hi:.3f}':<14} {dr:.4f}   {outcome}")

# secp256k1 reference
print(f"\n  secp256k1 (256b):              δ/n≈0.023   [projected: obstructed, same as Curve A]")
print()
print("Done.")
