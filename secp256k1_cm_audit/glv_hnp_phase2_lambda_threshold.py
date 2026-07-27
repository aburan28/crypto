"""
GLV-HNP Phase 2: λ/n threshold sweep (Thread 20).

Goal: bisect the λ/n threshold for the GLV-aware column-scaled lattice attack.
Known data points (from 2026-06-15 + 2026-07-26 runs with K1=8, K2≈√n):
  λ/n = 0.07 (p=2677, n=2647, λ=185):  LLL FAIL for all m in [5,12]
  λ/n ≥ 0.34 (various 8-12bit curves): LLL PASS

Structural hypothesis: when λ/n < threshold, the -λ·S_K1 column entry
in the k2-row is small relative to the n·S_K1 modular row. After LLL,
the k2 constraint gets absorbed into modular-row combinations and d is
not recoverable.

Correct comparison: fix K1=8, K2=√n, n∈[1000,6000] and vary only λ/n.
Use m=6 fixed (lattice dim=14) across all curves.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random

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
    m_r, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m_r - i - 1), p)
        m_r, c, t, r = i, b * b % p, t * b * b % p, r * b % p

# ---------------------------------------------------------------------------
# Primality (Miller-Rabin, deterministic for n < 3.2e9)
# ---------------------------------------------------------------------------

def is_prime(n):
    if n < 2: return False
    if n in (2, 3, 5, 7, 11, 13): return True
    if n % 2 == 0: return False
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

# ---------------------------------------------------------------------------
# Eisenstein CM for j=0 curves
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

def glv_min_eigenvalue(n):
    """Smaller root of x²+x+1 ≡ 0 (mod n); requires n ≡ 1 (mod 3)."""
    if n % 3 != 1: return None
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + n - sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return min(r1, r2)

def find_generator(p, b, n, tries=10000):
    rng = random.Random((p * 31337) ^ b)
    for _ in range(tries):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            if ec_mul((x, y), n, p) is None:
                return (x, y)
    return None

def find_twist_b(p, n, tries=1000):
    """Find b ∈ [1,tries) with #E_{y²=x³+b}(F_p) = n."""
    rng = random.Random(p ^ (n * 7))
    for b in range(1, tries):
        for _ in range(40):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                if ec_mul((x, y), n, p) is None:
                    return b
                break  # point not on order-n curve, next b
    return None

# ---------------------------------------------------------------------------
# Curve database by λ/n bin
# ---------------------------------------------------------------------------

TARGET_BINS = [
    (0.04, 0.09),   # A: prior FAIL at λ/n=0.07 (p=2677,n=2647)
    (0.09, 0.14),   # B
    (0.14, 0.20),   # C
    (0.20, 0.27),   # D
    (0.27, 0.35),   # E
    (0.35, 0.43),   # F
    (0.43, 0.50),   # G: prior PASS region
]
BIN_LABELS = 'ABCDEFG'

N_MIN = 1000   # n must be at least 1000 for non-trivial lattices
N_MAX = 6000   # keep runtime manageable

def find_curves_per_bin():
    """Search j=0 curves with n ∈ [N_MIN, N_MAX] and p ≤ N_MAX+200."""
    bins = [None] * len(TARGET_BINS)
    p = 3
    p_max = N_MAX + 300  # n ≈ p, so p slightly above N_MAX
    while p < p_max:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n = p + 1 - t
                    if n < N_MIN or n > N_MAX: continue
                    if not is_prime(n) or n % 3 != 1: continue
                    lam = glv_min_eigenvalue(n)
                    if lam is None: continue
                    ratio = lam / n
                    for i, (lo, hi) in enumerate(TARGET_BINS):
                        if bins[i] is None and lo <= ratio < hi:
                            b_twist = find_twist_b(p, n)
                            if b_twist is None: continue
                            G = find_generator(p, b_twist, n)
                            if G is None: continue
                            bins[i] = {'p': p, 'b': b_twist, 'n': n,
                                       'lam': lam, 'ratio': ratio, 'G': G}
                            break
        if all(c is not None for c in bins):
            break
        p += 2 if p > 2 else 1
        while p < p_max and not is_prime(p):
            p += 2
    return bins

# ---------------------------------------------------------------------------
# Attack (K1=8 to match prior failure benchmark)
# ---------------------------------------------------------------------------

K1_BOUND = 8   # match 2026-06-15 / 2026-07-26 sessions

def gen_signatures(curve, d_secret, m, seed):
    p, b, n, lam, G = curve['p'], curve['b'], curve['n'], curve['lam'], curve['G']
    k2_bound = math.isqrt(n) + 1
    rng = random.Random(seed)
    sigs, attempts = [], 0
    while len(sigs) < m and attempts < 2_000_000:
        attempts += 1
        k1 = rng.randint(0, K1_BOUND - 1)
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

def build_lattice(sigs, n, lam):
    m = len(sigs)
    dim = 2 * m + 2
    k2_bound = math.isqrt(n) + 1
    S_K1 = n // K1_BOUND
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
    return M, S_KANNAN

def try_recover(reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        if abs(row[dim - 1]) != S_KANNAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret: return True
    return False

def run_attack(curve, m, d_secret, seed, use_bkz=False, bkz_beta=20):
    n, lam = curve['n'], curve['lam']
    sigs = gen_signatures(curve, d_secret, m, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

# ---------------------------------------------------------------------------
# Sweep at fixed m and multiple m values
# ---------------------------------------------------------------------------

def sweep(curve, m_values, seeds, use_bkz=False, bkz_beta=20):
    results = {}
    for m in m_values:
        wins = 0
        for seed in seeds:
            d = random.Random(seed ^ 0xBEEF_CAFE).randint(1, curve['n'] - 1)
            wins += run_attack(curve, m, d, seed, use_bkz, bkz_beta)
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 7, 314159]
M_FIXED = 6  # fixed m for all curves (lattice dim = 14)

print("=" * 74)
print("GLV-HNP Phase 2 — λ/n threshold sweep (Thread 20)")
print(f"K1_BOUND={K1_BOUND} (matching prior failure benchmark), seeds={SEEDS}")
print(f"n ∈ [{N_MIN},{N_MAX}], lattice dim = {2*M_FIXED+2} (m={M_FIXED} fixed)")
print(f"Known: λ/n=0.07 (p=2677,n=2647) → LLL FAIL; λ/n≥0.34 → LLL PASS")
print("=" * 74)

# ---- Step 1: find curves ---------------------------------------------------
print(f"\n[Step 1] Searching for j=0 curves per λ/n bin ...")
curves = find_curves_per_bin()
for i, (lo, hi) in enumerate(TARGET_BINS):
    c = curves[i]
    if c:
        k2 = math.isqrt(c['n']) + 1
        eff = K1_BOUND * k2 / c['n']
        mt = math.ceil(math.log(c['n']) / math.log(1.0 / eff)) if eff < 1 else '∞'
        print(f"  Bin {BIN_LABELS[i]} [{lo:.2f},{hi:.2f}): p={c['p']}, n={c['n']} "
              f"({c['n'].bit_length()}b), λ={c['lam']}, λ/n={c['ratio']:.4f}, "
              f"K2={k2}, eff={eff:.3f}, m_thresh≈{mt}")
    else:
        print(f"  Bin {BIN_LABELS[i]} [{lo:.2f},{hi:.2f}): NOT FOUND")

# Inject the known-failure curve from prior run as Bin A validation
known_fail = {
    'p': 2677, 'b': 2, 'n': 2647, 'lam': 185,
    'ratio': 185/2647, 'G': None
}
known_fail_G = find_generator(known_fail['p'], known_fail['b'], known_fail['n'])
known_fail['G'] = known_fail_G
print(f"\n  [Reference] p=2677, n=2647, λ=185, λ/n={185/2647:.4f} — prior FAIL curve")

# ---- Step 2: LLL sweep at m=M_FIXED ----------------------------------------
print(f"\n[Step 2] LLL at fixed m={M_FIXED}, {len(SEEDS)} seeds")
lll_fixed = {}
for i, (lo, hi) in enumerate(TARGET_BINS):
    c = curves[i]
    if c is None: continue
    res = sweep(c, [M_FIXED], SEEDS)
    wins, total = res[M_FIXED]
    lll_fixed[i] = wins
    status = "PASS" if wins == total else ("PARTIAL" if wins > 0 else "FAIL")
    print(f"  Bin {BIN_LABELS[i]} λ/n={c['ratio']:.4f}: {wins}/{total} [{status}]")

# Reference curve
if known_fail['G']:
    res_ref = sweep(known_fail, [M_FIXED], SEEDS)
    wins_ref, total_ref = res_ref[M_FIXED]
    print(f"  [Ref]      λ/n={185/2647:.4f}: {wins_ref}/{total_ref} "
          f"[{'PASS' if wins_ref==total_ref else 'FAIL'}]  (expected FAIL)")

# ---- Step 3: LLL sweep over m range ----------------------------------------
print(f"\n[Step 3] LLL sweep over m values (m=4..9) for all bins")
M_VALUES = [4, 5, 6, 7, 8, 9]
lll_sweep = {}
for i, (lo, hi) in enumerate(TARGET_BINS):
    c = curves[i]
    if c is None: continue
    res = sweep(c, M_VALUES, SEEDS[:3])  # 3 seeds for efficiency
    lll_sweep[i] = res
    any_pass = any(w == t for w, t in res.values())
    best_m = min((m for m, (w, t) in res.items() if w == t), default=None)
    parts = [f"m={m}:{w}/{t}" for m, (w, t) in sorted(res.items())]
    tag = f"PASS@m={best_m}" if best_m else "FAIL"
    print(f"  Bin {BIN_LABELS[i]} λ/n={c['ratio']:.4f} [{tag}]: " + "  ".join(parts))

# Reference
if known_fail['G']:
    res_ref_s = sweep(known_fail, M_VALUES, SEEDS[:3])
    parts_r = [f"m={m}:{w}/{t}" for m, (w, t) in sorted(res_ref_s.items())]
    any_pass_r = any(w == t for w, t in res_ref_s.values())
    print(f"  [Ref]      λ/n=0.0699 [{'PASS' if any_pass_r else 'FAIL'}]: " + "  ".join(parts_r))

# ---- Step 4: BKZ(20) on failing bins ----------------------------------------
fail_bins = [i for i in range(len(TARGET_BINS))
             if i in lll_sweep and
             not any(w == t for w, t in lll_sweep[i].values())]
print(f"\n[Step 4] BKZ(20) on {len(fail_bins)} strictly-failing bins: "
      f"{[BIN_LABELS[i] for i in fail_bins]}")
bkz_sweep = {}
for i in fail_bins:
    c = curves[i]
    if c is None: continue
    res = sweep(c, M_VALUES, SEEDS[:3], use_bkz=True, bkz_beta=20)
    bkz_sweep[i] = res
    any_pass = any(w == t for w, t in res.values())
    parts = [f"m={m}:{w}/{t}" for m, (w, t) in sorted(res.items())]
    tag = "PASS" if any_pass else "FAIL"
    print(f"  Bin {BIN_LABELS[i]} BKZ(20) [{tag}]: " + "  ".join(parts))

# ---- Step 5: Planted-vector norm vs GH estimate ----------------------------
print(f"\n[Step 5] Planted-vector geometry (m={M_FIXED})")
print(f"{'Bin':<5} {'λ/n':<10} {'||v_plant||':>12} {'GH_est':>10} {'ratio':>8} {'k2-row resid':>14}")
print("-" * 65)
for i, (lo, hi) in enumerate(TARGET_BINS):
    c = curves[i]
    if c is None: continue
    n, lam = c['n'], c['lam']
    k2_bound = math.isqrt(n) + 1
    S_K1 = n // K1_BOUND
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    m = M_FIXED
    # Typical planted norm: k1 avg=K1/2, k2 avg=k2_bound/2, d avg=n/2
    k1_avg = K1_BOUND / 2
    k2_avg = k2_bound / 2
    d_avg = n / 2
    v_norm = math.sqrt(m * (k1_avg * S_K1)**2 + d_avg**2 +
                       m * (k2_avg * S_K2)**2 + S_KANNAN**2)
    # GH estimate for lambda_1
    dim = 2 * m + 2
    log_det = (m * math.log(n * S_K1) + math.log(1) +
               m * math.log(S_K2) + math.log(S_KANNAN))
    gh = math.exp(1.0) * math.exp(log_det / dim)
    # k2-row GS residual magnitude: λ * S_K1 (proj onto mod-row = λ/n * n*S_K1 = λ*S_K1 → residual=0)
    # But the ABSOLUTE magnitude before projection:
    k2_entry = lam * S_K1  # magnitude of -lam*S_K1
    mod_entry = n * S_K1
    k2_resid = k2_entry / mod_entry * mod_entry  # ≡ k2_entry (the projection cancels it exactly)
    # The key ratio: λ/n (how much of mod-row is "visible" in k2-row)
    # After projection: k2-row residual is purely [0,...,S_K2,...,0] regardless of λ
    # But the SIZE REDUCTION step sees -lam*S_K1, and μ = -lam*S_K1 / (n*S_K1) = -λ/n
    mu_k2_mod = lam / n  # |μ_{k2-row, mod-row}|
    print(f"  {BIN_LABELS[i]}    {c['ratio']:.4f}    {v_norm:>12.0f} {gh:>10.0f} "
          f"{v_norm/gh:>8.2f}  μ={mu_k2_mod:.4f}")

# ---- Summary ----------------------------------------------------------------
print("\n" + "=" * 74)
print("SUMMARY")
print(f"{'Bin':<5} {'λ/n range':<14} {'λ/n actual':<12} {'n':<7} {'LLL(m=4..9)':<14} {'BKZ(20)'}")
print("-" * 72)
for i, (lo, hi) in enumerate(TARGET_BINS):
    c = curves[i]
    if c is None:
        print(f"{BIN_LABELS[i]:<5} [{lo:.2f},{hi:.2f})  {'N/A':<12}")
        continue
    lll_ok = (i in lll_sweep and any(w == t for w, t in lll_sweep[i].values()))
    bkz_ok = (i in bkz_sweep and any(w == t for w, t in bkz_sweep[i].values()))
    lll_s = "PASS" if lll_ok else "FAIL"
    bkz_s = ("PASS" if bkz_ok else "FAIL") if i in bkz_sweep else "n/a"
    print(f"{BIN_LABELS[i]:<5} [{lo:.2f},{hi:.2f})  {c['ratio']:.4f}       "
          f"{c['n']:<7} {lll_s:<14} {bkz_s}")
if known_fail['G']:
    lll_ref = any(w==t for w,t in res_ref_s.values()) if 'res_ref_s' in dir() else False
    print(f"[Ref]  [0.04,0.09)  0.0699       2647    "
          f"{'PASS' if lll_ref else 'FAIL'}           (prior-run benchmark)")

print()
# Threshold
pass_bins = [i for i in range(len(TARGET_BINS))
             if i in lll_sweep and any(w==t for w,t in lll_sweep[i].values())]
fail_bins_f = [i for i in range(len(TARGET_BINS))
               if i in lll_sweep and not any(w==t for w,t in lll_sweep[i].values())]
if pass_bins and fail_bins_f:
    max_fail = max(fail_bins_f)
    min_pass = min(pass_bins)
    r_f = curves[max_fail]['ratio'] if curves[max_fail] else 'N/A'
    r_p = curves[min_pass]['ratio'] if curves[min_pass] else 'N/A'
    print(f"LLL threshold (K1={K1_BOUND}): λ/n ∈ ({r_f:.4f}, {r_p:.4f})")
else:
    print("Monotone: all bins passed or all failed — check intermediate bins.")

print("\nDone.")
