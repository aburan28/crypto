"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Thread 20 from 2026-07-26 autolab run.

Goal: find the critical λ/n ratio below which LLL fails to recover d.
Known data points:
  - λ/n = 0.07  (p=2677,   n=2647):  LLL FAILS at all m, any algo
  - λ/n = 0.34  (p=524347, n=523969): LLL SUCCEEDS at m=9

This script finds 12-bit j=0 CM curves at λ/n targets between 0.07 and 0.34,
then sweeps m=m_thresh..m_thresh+8 with 5 seeds each to find the critical
threshold. K1=8 is fixed to match the 2677 failure-curve parameters, keeping
eff = K1 * sqrt(n) / n consistent across all curves.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sys
import sympy
from fpylll import IntegerMatrix, LLL

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
    mv, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mv - i - 1), p)
        mv, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(200000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: j=0 curves (p ≡ 1 mod 3), Eisenstein decomposition
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
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def find_b_for_n(p, n):
    """Find b s.t. #E(F_p): y²=x³+b equals n. Tries b=1..500."""
    rng = random.Random(99)
    for b_try in range(1, min(p, 600)):
        for _ in range(30):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b_try
                break
    return None

# ---------------------------------------------------------------------------
# Curve search: find 12-bit j=0 curves at target λ/n bins
# ---------------------------------------------------------------------------

# Target bins: bisect between 0.07 (fail) and 0.34 (success)
# Use half-width ±0.018 so bins don't overlap and cover the range
LAM_N_TARGETS = [0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30, 0.33]
LAM_N_TOL = 0.018

# Also include the two anchor points from 2026-07-26 run
ANCHOR_FAIL  = dict(p=2677,   b=2, n=2647,   lam=185,    label="anchor-FAIL",  lamn=0.0699)
ANCHOR_SUCC  = dict(p=524347, b=2, n=523969, lam=177902, label="anchor-SUCC",  lamn=0.3397)


def find_curves_for_bins(targets=LAM_N_TARGETS, tol=LAM_N_TOL,
                          n_min=1200, n_max=6000):
    """
    Scan small primes p ≡ 1 mod 3 and collect j=0 CM curves at each λ/n target.
    Returns dict: target → (p, b, n, lam, lam_n) or None if not found.
    """
    found = {t: None for t in targets}
    p = 3  # start from small primes
    scanned = 0

    while p < 100000 and any(v is None for v in found.values()):
        p = sympy.nextprime(p)
        scanned += 1
        if p % 3 != 1:
            continue
        eis = eisenstein_decompose(p)
        if eis is None:
            continue
        a_e, b_e = eis
        for tr in j0_traces(a_e, b_e):
            n_cand = p + 1 - tr
            if n_cand < n_min or n_cand > n_max:
                continue
            if not sympy.isprime(n_cand):
                continue
            if n_cand % 3 != 1:
                continue
            lam = glv_eigenvalue(n_cand)
            if lam is None:
                continue
            ratio = lam / n_cand
            for t in targets:
                if found[t] is not None:
                    continue
                if abs(ratio - t) <= tol:
                    b_param = find_b_for_n(p, n_cand)
                    if b_param is None:
                        continue
                    found[t] = (p, b_param, n_cand, lam, ratio)
                    print(f"  [bin {t:.2f}] p={p}, b={b_param}, n={n_cand}, "
                          f"lam={lam}, lam/n={ratio:.4f}")
                    break

    missing = [t for t in targets if found[t] is None]
    if missing:
        print(f"  WARNING: no curves found for bins {missing} "
              f"(searched {scanned} primes up to p < 100000)")
    return found

# ---------------------------------------------------------------------------
# GLV-HNP attack infrastructure
# ---------------------------------------------------------------------------

K1 = 8  # fixed K1 to match 2677 failure-curve parameters

def gen_signatures(p, b, n, lam, G, m, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
        attempts += 1
        k1_v = rng.randint(0, k1_bound - 1)
        k2_v = rng.randint(0, k2_bound - 1)
        k_full = (k1_v + lam * k2_v) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        sigs.append({'A': h * modinv(s, n) % n, 'B': r * modinv(s, n) % n})
    return d_secret, sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
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
    M[m][m] = 1
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
        if abs(row[dim - 1]) != S_KANNAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_lll(p, b, n, lam, G, m, k1_bound, k2_bound, seed):
    d_secret, sigs = gen_signatures(p, b, n, lam, G, m, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

SEEDS = [42, 1234, 9999, 7, 314159]

def sweep_curve(p, b, n, lam, G, k1_bound, m_range, label=""):
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
    lam_n = lam / n
    delta = min((3 * lam) % n, n - (3 * lam) % n) / n   # |lam/n - 1/3| proxy

    print(f"  lam/n={lam_n:.4f}  delta/n={delta:.4f}  "
          f"K1={k1_bound} K2={k2_bound} eff={eff:.4f} m_thresh≈{m_thresh:.1f}")
    sys.stdout.flush()

    results = {}
    for m in m_range:
        wins = sum(
            run_lll(p, b, n, lam, G, m, k1_bound, k2_bound, seed)
            for seed in SEEDS
        )
        results[m] = (wins, len(SEEDS))
        bar = "*" * wins + "." * (len(SEEDS) - wins)
        print(f"    m={m:2d}: {wins}/{len(SEEDS)} [{bar}]", end="")
        if wins == len(SEEDS):
            print("  ← 5/5 ✓")
        elif wins >= 3:
            print("  ← 3/5+ ✓")
        else:
            print()
        sys.stdout.flush()

    first_win = next((m for m, (w, t) in sorted(results.items()) if w >= 3), None)
    return results, first_win, m_thresh, eff, lam_n, delta

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection (Thread 20)")
print("=" * 70)
print(f"K1_BOUND = {K1}  (fixed, matching 2677 failure-curve parameters)")
print(f"Targets: λ/n ∈ {LAM_N_TARGETS}")
print(f"Tolerance: ±{LAM_N_TOL}")
print()

# ---- Step 1: Find curves for each bin -------------------------------------
print("Searching for 12-bit j=0 CM curves in each λ/n bin...")
bin_curves = find_curves_for_bins()
print()

# ---- Step 2: Sweep each curve -------------------------------------------
print("=" * 70)
print("LLL sweeps (m_thresh..m_thresh+8, 5 seeds each)")
print("=" * 70)

summary_rows = []  # (lam_n, p, n, lam, first_win, m_thresh)

# Anchor: FAIL curve (p=2677, known from 2026-07-26)
print(f"\n[ANCHOR-FAIL: p=2677, n=2647, lam=185, lam/n=0.0699]")
p_af, b_af, n_af, lam_af = (ANCHOR_FAIL['p'], ANCHOR_FAIL['b'],
                              ANCHOR_FAIL['n'], ANCHOR_FAIL['lam'])
G_af = find_generator(p_af, b_af, n_af)
if G_af is not None:
    _, fw_af, mt_af, eff_af, ln_af, dl_af = sweep_curve(
        p_af, b_af, n_af, lam_af, G_af, K1, range(5, 13), label="ANCHOR-FAIL")
    summary_rows.append((ln_af, p_af, n_af, lam_af, fw_af, mt_af))
else:
    print("  ERROR: generator not found for anchor-fail curve")

# Intermediate curves
for t in LAM_N_TARGETS:
    entry = bin_curves[t]
    if entry is None:
        print(f"\n[SKIP bin {t:.2f}: no curve found]")
        summary_rows.append((t, None, None, None, None, None))
        continue
    p_c, b_c, n_c, lam_c, ratio_c = entry
    print(f"\n[bin {t:.2f}: p={p_c}, n={n_c}, lam={lam_c}, lam/n={ratio_c:.4f}]")
    G_c = find_generator(p_c, b_c, n_c)
    if G_c is None:
        print("  ERROR: generator not found; skipping")
        summary_rows.append((ratio_c, p_c, n_c, lam_c, None, None))
        continue
    k2_c = math.isqrt(n_c) + 1
    eff_c = K1 * k2_c / n_c
    m_thresh_c = math.ceil(math.log(n_c) / math.log(1.0 / eff_c)) if eff_c < 1 else 99
    m_lo = max(4, m_thresh_c)
    m_hi = m_lo + 9
    _, fw_c, mt_c, _, ln_c, dl_c = sweep_curve(
        p_c, b_c, n_c, lam_c, G_c, K1, range(m_lo, m_hi), label=f"bin{t:.2f}")
    summary_rows.append((ln_c, p_c, n_c, lam_c, fw_c, mt_c))

# Anchor: SUCCESS curve (p=524347, K1=8 used here for comparability)
print(f"\n[ANCHOR-SUCC: p=524347, n=523969, lam=177902, lam/n=0.3397]")
p_as = ANCHOR_SUCC['p']; b_as = ANCHOR_SUCC['b']
n_as = ANCHOR_SUCC['n']; lam_as = ANCHOR_SUCC['lam']
G_as = find_generator(p_as, b_as, n_as)
if G_as is not None:
    k2_as = math.isqrt(n_as) + 1
    eff_as = K1 * k2_as / n_as
    m_thresh_as = math.ceil(math.log(n_as) / math.log(1.0 / eff_as)) if eff_as < 1 else 99
    _, fw_as, mt_as, _, ln_as, dl_as = sweep_curve(
        p_as, b_as, n_as, lam_as, G_as, K1, range(m_thresh_as, m_thresh_as + 9),
        label="ANCHOR-SUCC")
    summary_rows.append((ln_as, p_as, n_as, lam_as, fw_as, mt_as))
else:
    print("  ERROR: generator not found for anchor-succ curve")

# ---- Step 3: Summary table -----------------------------------------------
print()
print("=" * 70)
print("THRESHOLD SUMMARY TABLE")
print("=" * 70)
print(f"{'lam/n':>8}  {'n':>8}  {'lam':>8}  {'m_thresh':>9}  {'first 3/5+':>11}  outcome")
print("-" * 65)

last_fail = None
first_succ = None

for (lam_n, p_r, n_r, lam_r, fw_r, mt_r) in sorted(summary_rows, key=lambda x: x[0]):
    if lam_n is None:
        continue
    fw_str  = str(fw_r) if fw_r is not None else "never"
    mt_str  = str(mt_r) if mt_r is not None else "?"
    n_str   = str(n_r) if n_r is not None else "?"
    lam_str = str(lam_r) if lam_r is not None else "?"
    if fw_r is None:
        outcome = "FAIL"
        last_fail = lam_n
    else:
        outcome = f"PASS (m={fw_r})"
        if first_succ is None:
            first_succ = lam_n
    print(f"  {lam_n:6.4f}  {n_str:>8}  {lam_str:>8}  {mt_str:>9}  "
          f"{fw_str:>11}  {outcome}")

print()
if last_fail is not None and first_succ is not None:
    print(f"Threshold bracket: λ/n ∈ ({last_fail:.4f}, {first_succ:.4f})")
    mid = (last_fail + first_succ) / 2
    print(f"Midpoint estimate: λ/n ≈ {mid:.4f}")
elif last_fail is not None:
    print(f"All tested curves FAIL. Need curves with λ/n > {last_fail:.4f} in new bins.")
elif first_succ is not None:
    print(f"All tested curves PASS. Threshold may be below {first_succ:.4f}.")
else:
    print("No conclusive bracket found.")

print()
print("Done.")
