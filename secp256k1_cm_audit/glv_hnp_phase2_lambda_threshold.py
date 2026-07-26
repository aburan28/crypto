"""
GLV-HNP Phase 2: Lambda/n threshold bisection study.

Motivation (from 2026-07-26 log):
  - lambda/n = 0.07 (p=2677): LLL fails, BKZ(40) fails (structural)
  - lambda/n = 0.34 (p=524347): LLL succeeds at m=9
  - lambda/n = 0.53 (p=211): LLL succeeds at m=4
  - lambda/n = 0.66 (p=2557): LLL succeeds at m=7

Goal: bisect the threshold between 0.07 and 0.34 by finding j=0 curves
  at lambda/n ≈ {0.10, 0.15, 0.20, 0.25, 0.28} and running LLL+BKZ sweeps.
  Report: minimum lambda/n where LLL succeeds 3/3 seeds.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass a=0)
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
# CM theory helpers (Eisenstein, j=0 traces, GLV eigenvalue)
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

def find_b_for_curve(p, n):
    """Find curve parameter b such that #(y^2=x^3+b)(F_p) = n."""
    rng = random.Random(99)
    # Try a range of b values
    for b in list(range(1, 200)) + [rng.randint(200, p-1) for _ in range(50)]:
        rhs0 = (pow(0, 3, p) + b) % p
        P = None
        for x in range(p):
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                break
        if P is None: continue
        if ec_mul(P, n, p) is None:
            return b
    return None

# ---------------------------------------------------------------------------
# Find curves with specific lambda/n ranges
# ---------------------------------------------------------------------------

TARGET_BINS = [
    (0.08, 0.12, "bin_010"),
    (0.12, 0.17, "bin_015"),
    (0.17, 0.23, "bin_020"),
    (0.23, 0.27, "bin_025"),
    (0.27, 0.31, "bin_028"),
    (0.31, 0.37, "bin_033"),
]

def find_curves_by_lambda_ratio(p_min, p_max, max_per_bin=2):
    """
    Search j=0 curves over F_p for p in [p_min, p_max].
    Collect up to max_per_bin curves per lambda/n ratio bin.
    Returns dict: bin_label -> list of (p, b, n, lam, lam/n, G)
    """
    bins = {label: [] for (_, _, label) in TARGET_BINS}
    full_bins = set()

    p = sympy.nextprime(p_min - 1)
    checked = 0
    while p <= p_max and len(full_bins) < len(TARGET_BINS):
        checked += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand >= p: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam, _ = glv_eigenvalue(n_cand)
                    if lam is None: continue
                    ratio = lam / n_cand

                    for (lo, hi, label) in TARGET_BINS:
                        if label in full_bins: continue
                        if lo <= ratio < hi:
                            b_curve = find_b_for_curve(p, n_cand)
                            if b_curve is None: break
                            G = find_generator(p, b_curve, n_cand)
                            if G is None: break
                            bins[label].append((p, b_curve, n_cand, lam, ratio, G))
                            if len(bins[label]) >= max_per_bin:
                                full_bins.add(label)
                            break
        p = sympy.nextprime(p)

    print(f"  (Checked {checked} primes in [{p_min}, {p_max}])")
    return bins

# ---------------------------------------------------------------------------
# GLV signature generation
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

# ---------------------------------------------------------------------------
# GLV lattice and recovery
# ---------------------------------------------------------------------------

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
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def run_once(curve_params, m, k1_bound, seed=42, use_bkz=False, bkz_beta=20):
    p, b, n, lam, ratio, G = curve_params
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
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
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

def sweep(label, curve_params, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    """Run lattice attack for multiple m and seeds; return first m where all seeds succeed."""
    p, b, n, lam, ratio, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1.0 else float('inf')
    algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"

    print(f"  [{algo}] {label}  lam/n={ratio:.4f}  n={n}b={n.bit_length()}  "
          f"K1={k1_bound}  eff={eff:.4f}  m_thresh≈{m_thresh:.1f}")

    first_all = None
    for m in m_range:
        wins = sum(run_once(curve_params, m, k1_bound, s, use_bkz, bkz_beta) for s in seeds)
        marker = "*" if m >= m_thresh else " "
        print(f"    m={m}{marker}: {wins}/{len(seeds)}", end="")
        if wins == len(seeds) and first_all is None:
            first_all = m
            print("  <- first 3/3", end="")
        print()
    return first_all

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-aware HNP Phase 2 — lambda/n threshold bisection")
print("=" * 70)
print()
print("Prior data points:")
print("  p=2677,   n=2647,   lam=185,   lam/n=0.070: FAIL (LLL+BKZ all m)")
print("  p=524347, n=523969, lam=~0.34n:             SUCCESS (LLL m=9)")
print("  p=211,    n=199,    lam=106,   lam/n=0.53:  SUCCESS (LLL m=4)")
print("  p=2557,   n=2659,   lam=1755,  lam/n=0.66:  SUCCESS (LLL m=7)")
print()
print(f"Searching for curves at lam/n in {[f'{lo:.2f}-{hi:.2f}' for lo,hi,_ in TARGET_BINS]}...")
print()

SEEDS = [42, 1234, 9999]
# Search 12-16 bit primes (3000..65000) for targeted lambda/n ratios
bins = find_curves_by_lambda_ratio(3000, 65000, max_per_bin=2)

print("\nFound curves per bin:")
for lo, hi, label in TARGET_BINS:
    curves = bins[label]
    if curves:
        for (p, b, n, lam, ratio, G) in curves:
            print(f"  {label} [{lo:.2f}-{hi:.2f}]: p={p}, n={n}({n.bit_length()}b), "
                  f"lam={lam}, lam/n={ratio:.4f}")
    else:
        print(f"  {label} [{lo:.2f}-{hi:.2f}]: NOT FOUND in range")

print()
print("=" * 70)
print("Running LLL and BKZ sweeps per bin")
print("=" * 70)

# K1_BOUND selection: keep eff ~ 0.05 so m_thresh ≈ 5-8
# eff = K1 * sqrt(n) / n = K1 / sqrt(n)
# K1 ≈ 0.05 * sqrt(n)

# Store threshold results: ratio -> (lll_first_3of3, bkz20_first_3of3, bkz40_first_3of3)
threshold_data = []

M_RANGE = range(3, 20)  # sweep m from 3 to 19

for lo, hi, label in TARGET_BINS:
    curves = bins[label]
    if not curves:
        print(f"\n[{label}] No curve found — skipping")
        continue

    curve_params = curves[0]   # use first found curve
    p, b, n, lam, ratio, G = curve_params

    k1_bound = max(2, int(0.05 * math.sqrt(n)))
    print(f"\n[{label}] lam/n={ratio:.4f}  p={p}  n={n}  K1={k1_bound}")

    lll_thresh = sweep(f"y^2=x^3+{b}/F_{p}", curve_params, k1_bound,
                       M_RANGE, SEEDS, use_bkz=False)

    # Only run BKZ if LLL fails (or at a higher m)
    if lll_thresh is None:
        bkz20_thresh = sweep(f"y^2=x^3+{b}/F_{p}", curve_params, k1_bound,
                             M_RANGE, SEEDS, use_bkz=True, bkz_beta=20)
        if bkz20_thresh is None:
            bkz40_thresh = sweep(f"y^2=x^3+{b}/F_{p}", curve_params, k1_bound,
                                 range(3, 16), SEEDS, use_bkz=True, bkz_beta=40)
        else:
            bkz40_thresh = None  # BKZ(20) already succeeded; skip BKZ(40)
    else:
        bkz20_thresh = None
        bkz40_thresh = None

    threshold_data.append({
        'label': label, 'ratio': ratio,
        'n': n, 'p': p, 'lam': lam, 'k1_bound': k1_bound,
        'lll': lll_thresh, 'bkz20': bkz20_thresh, 'bkz40': bkz40_thresh,
    })

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"{'lam/n':>8} {'n_bits':>6} {'K1':>4} {'LLL 3/3 at':>12} {'BKZ20 3/3':>10} {'BKZ40 3/3':>10}")
print("-" * 70)

# Known reference points
known = [
    (0.070, 12, 8, None, None, None, "2677 (known FAIL)"),
]
for ratio, nb, k1, lll, bkz20, bkz40, note in known:
    lll_s = f"m={lll}" if lll else "FAIL"
    b20_s = f"m={bkz20}" if bkz20 else "FAIL"
    b40_s = f"m={bkz40}" if bkz40 else ""
    print(f"{ratio:8.3f} {nb:6} {k1:4} {lll_s:>12} {b20_s:>10} {b40_s:>10}  {note}")

for d in threshold_data:
    ratio = d['ratio']
    nb = d['n'].bit_length()
    k1 = d['k1_bound']
    lll_s = f"m={d['lll']}" if d['lll'] else "FAIL"
    b20_s = f"m={d['bkz20']}" if d['bkz20'] else ("skip" if d['lll'] else "FAIL")
    b40_s = f"m={d['bkz40']}" if d['bkz40'] else ("skip" if (d['lll'] or d['bkz20']) else "FAIL")
    print(f"{ratio:8.3f} {nb:6} {k1:4} {lll_s:>12} {b20_s:>10} {b40_s:>10}")

# Reference success points
for ratio, label in [(0.34, "524347 (known OK)"), (0.53, "211 (known OK)"), (0.66, "2557 (known OK)")]:
    print(f"{ratio:8.3f} {'?':>6} {'?':>4} {'OK':>12} {'skip':>10} {'skip':>10}  {label}")

print()
print("Interpretation:")
# Find threshold
fail_ratios = [d['ratio'] for d in threshold_data if d['lll'] is None and d['bkz20'] is None and d['bkz40'] is None]
pass_ratios = [d['ratio'] for d in threshold_data if d['lll'] is not None]
if fail_ratios and pass_ratios:
    top_fail = max(fail_ratios)
    bot_pass = min(pass_ratios)
    print(f"  LLL threshold: lam/n in ({top_fail:.3f}, {bot_pass:.3f})")
elif not fail_ratios:
    print(f"  LLL succeeds for all tested ratios >= {min(d['ratio'] for d in threshold_data):.3f}")
elif not pass_ratios:
    print(f"  LLL fails for all tested ratios <= {max(d['ratio'] for d in threshold_data):.3f}")
else:
    print("  Inconclusive (need more data points).")

bkz_rescues = [d['ratio'] for d in threshold_data if d['lll'] is None and (d['bkz20'] or d['bkz40'])]
if bkz_rescues:
    print(f"  BKZ rescues LLL at ratios: {[f'{r:.3f}' for r in bkz_rescues]}")
else:
    print("  BKZ does NOT rescue LLL failures in any bin.")

print()
print("Done.")
