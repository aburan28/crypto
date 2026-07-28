"""
GLV-HNP Phase 2 — λ/n threshold bisection study (Thread 20).

Prior results (from glv_hnp_phase2_20bit.py, 2026-07-26):
  lam/n = 0.07  (p=2677, n=2647, K1=8):  LLL fails at all m, BKZ(40) also fails
  lam/n = 0.34  (p=524347, n=523969):     LLL succeeds at m=9

Goal: sweep λ/n ∈ [0.05, 0.45] using 12-bit j=0 CM curves.
Use K1=8 (fixed, matching prior p=2677 experiment) so lam/n is the only variable.

Method:
  - Scan primes p in [4096, 32768] with p ≡ 1 (mod 3)
  - For each valid (p, n, λ), bin by lam/n into 9 bins
  - Take the first 2 curves per bin
  - Run LLL at m=3..16 with 5 seeds; run BKZ(20) if LLL fails

Run: python3 glv_hnp_phase2_lambda_threshold.py
Expected runtime: 3-6 minutes.
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (a=0 short Weierstrass over F_p)
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
    rng = random.Random(12345 + p)
    for _ in range(30000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# Eisenstein CM: decompose p = a^2 - ab + b^2, list traces, get lambda
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

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

# ---------------------------------------------------------------------------
# Scan primes, bin by lam/n
# ---------------------------------------------------------------------------

def find_curves_by_lam_bins(p_lo, p_hi, bin_centers, width, max_per_bin=2):
    """
    Returns dict: bin_center (float) -> list of (p, b, n, lam, G)
    """
    results = {c: [] for c in bin_centers}
    full = set()

    p = sympy.nextprime(p_lo - 1)
    while p <= p_hi and len(full) < len(bin_centers):
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_c = p + 1 - t
                    if n_c < 200: continue
                    if not sympy.isprime(n_c): continue
                    if n_c % 3 != 1: continue
                    lam = glv_eigenvalue(n_c)
                    if lam is None: continue
                    ratio = lam / n_c
                    for c in bin_centers:
                        if c in full: continue
                        if abs(ratio - c) <= width:
                            # Find curve twist b and generator
                            found_b, found_G = None, None
                            for b_try in range(1, min(p, 200)):
                                G = find_generator(p, b_try, n_c)
                                if G is not None:
                                    found_b = b_try
                                    found_G = G
                                    break
                            if found_b is None: continue
                            results[c].append((p, found_b, n_c, lam, found_G))
                            if len(results[c]) >= max_per_bin:
                                full.add(c)
                            break  # don't put same curve in multiple bins
        p = sympy.nextprime(p)
    return results

# ---------------------------------------------------------------------------
# GLV lattice + LLL attack
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
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
        sigs.append({'A': A, 'B': B})
    return sigs

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

def run_trial(p, b, n, lam, G, m, d_secret, k1_bound, k2_bound, seed, use_bkz, bkz_beta):
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[i][m]) % n
        if d_cand == d_secret:
            return True
    return False

def sweep_m(p, b, n, lam, G, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    k2_bound = math.isqrt(n) + 1
    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 31337).randint(1, n - 1)
            ok = run_trial(p, b, n, lam, G, m, d, k1_bound, k2_bound, seed, use_bkz, bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
    return results

def first_win(res, need=3):
    for m in sorted(res):
        if res[m][0] >= need: return m
    return None

def fmt_res(res):
    return "  ".join(f"m={m}:{w}/{t}" + ("*" if w>=3 else "") for m,(w,t) in sorted(res.items()))

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 72)
print("GLV-HNP Phase 2 — λ/n threshold bisection (Thread 20)")
print("Fixed K1=8 across all curves (matches 2026-07-26 p=2677 experiment).")
print("=" * 72)

BIN_CENTERS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
WIDTH = 0.028   # ±0.028 window per bin
K1 = 8
SEEDS = [42, 1234, 9999, 7, 314159]
M_RANGE = range(3, 17)

print(f"\nSweep config: K1={K1}, seeds={SEEDS}, m ∈ {list(M_RANGE)}")
print(f"\nScanning p ∈ [4096, 32768] for λ/n bins {BIN_CENTERS}...")

curves_by_bin = find_curves_by_lam_bins(4096, 32768, BIN_CENTERS, WIDTH, max_per_bin=2)

print("\nCurves found per bin:")
for c in BIN_CENTERS:
    clist = curves_by_bin[c]
    print(f"  λ/n ≈ {c:.2f}: {len(clist)} curve(s)", end="")
    for p, b, n, lam, G in clist:
        print(f"  (p={p}, n={n}, lam={lam}, lam/n={lam/n:.3f})", end="")
    print()

# Add anchor curve (p=2677, known failure from 2026-07-26 at K1=8)
G_2677 = find_generator(2677, 2, 2647)
anchor = (2677, 2, 2647, 185, G_2677)

print("\n" + "=" * 72)
print("Running sweep...")
print("=" * 72)

all_results = []  # (ratio, label, lll_win, bkz_win)

# ---- Anchor: p=2677, lam/n=0.070 ----
print(f"\n[ANCHOR] p=2677, n=2647, lam=185, lam/n=0.070 (K1={K1})")
k2_anc = math.isqrt(2647) + 1
print(f"  eff = K1*K2/n = {K1}*{k2_anc}/{2647} = {K1*k2_anc/2647:.4f}")
res_lll = sweep_m(2677, 2, 2647, 185, G_2677, K1, M_RANGE, SEEDS)
win_lll = first_win(res_lll)
print(f"  LLL: {fmt_res(res_lll)}")
print(f"  → LLL wins at m={win_lll}")
bkz_win = None
if win_lll is None:
    res_bkz = sweep_m(2677, 2, 2647, 185, G_2677, K1, M_RANGE, SEEDS, use_bkz=True, bkz_beta=20)
    bkz_win = first_win(res_bkz)
    print(f"  BKZ20: {fmt_res(res_bkz)}")
    print(f"  → BKZ20 wins at m={bkz_win}")
all_results.append((0.070, "p=2677/n=2647", win_lll, bkz_win))

# ---- Per-bin sweep ----
for c in BIN_CENTERS:
    clist = curves_by_bin[c]
    if not clist:
        print(f"\n[λ/n≈{c:.2f}] NO CURVES FOUND — skipping bin")
        continue
    for idx, (p, b, n, lam, G) in enumerate(clist):
        ratio = lam / n
        k2 = math.isqrt(n) + 1
        eff = K1 * k2 / n
        print(f"\n[λ/n≈{c:.2f} #{idx+1}] p={p}, n={n}, lam={lam}, "
              f"lam/n={ratio:.3f}, K1={K1}, K2={k2}, eff={eff:.4f}")
        res_lll = sweep_m(p, b, n, lam, G, K1, M_RANGE, SEEDS)
        win_lll = first_win(res_lll)
        print(f"  LLL: {fmt_res(res_lll)}")
        print(f"  → LLL wins at m={win_lll}")
        bkz_win = None
        if win_lll is None:
            res_bkz = sweep_m(p, b, n, lam, G, K1, M_RANGE, SEEDS, use_bkz=True, bkz_beta=20)
            bkz_win = first_win(res_bkz)
            print(f"  BKZ20: {fmt_res(res_bkz)}")
            print(f"  → BKZ20 wins at m={bkz_win}")
        all_results.append((ratio, f"p={p}/n={n}", win_lll, bkz_win))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 72)
print("SUMMARY TABLE — λ/n vs LLL/BKZ threshold")
print("K1=8 fixed; K2=floor(sqrt(n))+1; 5 seeds; need 3/5 success")
print("=" * 72)
print(f"{'lam/n':>7}  {'curve':>18}  {'LLL m*':>8}  {'BKZ20 m*':>10}")
print("-" * 50)
for ratio, label, lll_m, bkz_m in sorted(all_results):
    lll_s = str(lll_m) if lll_m is not None else "FAIL"
    bkz_s = ("(not run)" if lll_m is not None else
              (str(bkz_m) if bkz_m is not None else "FAIL"))
    print(f"  {ratio:.3f}  {label:>18}  {lll_s:>8}  {bkz_s:>10}")

print()
lll_success = sorted([(r, m) for r, _, m, _ in all_results if m is not None])
lll_fail    = sorted([(r, _) for r, _, m, _ in all_results if m is None])
if lll_success and lll_fail:
    max_fail = max(r for r, _ in lll_fail)
    min_succ = min(r for r, _ in lll_success)
    print(f"Empirical LLL threshold: λ/n in ({max_fail:.3f}, {min_succ:.3f})")
    print(f"  (above {min_succ:.3f}: LLL succeeds; below {max_fail:.3f}: fails)")
elif lll_success:
    min_succ = min(r for r, _ in lll_success)
    print(f"LLL succeeds at all tested λ/n (min={min_succ:.3f})")
else:
    max_fail = max(r for r, _ in lll_fail)
    print(f"LLL fails at all tested λ/n (max={max_fail:.3f})")

print("\nDone.")
