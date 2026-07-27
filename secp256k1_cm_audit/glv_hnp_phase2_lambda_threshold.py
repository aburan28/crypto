"""
GLV-HNP Phase 2: λ/n threshold bisection study.

Bisects the success/failure boundary for the GLV-aware HNP LLL attack
in the range λ/n ∈ (0.07, 0.34).

Known data (from glv_hnp_phase2_20bit.py, 2026-07-26):
  λ/n = 0.07 (p=2677, n=2647): LLL FAILS at all m ≤ 12, BKZ-40 also FAILS
  λ/n = 0.34 (p=524347, n=523969): LLL SUCCEEDS at m=9

This script samples λ/n targets: 0.09, 0.12, 0.15, 0.18, 0.21, 0.27
For each target, finds a small (~12-16 bit) j=0 CM curve with that λ/n
and runs the GLV-aware LLL sweep (m=4..20, 3 seeds).

Reports the minimum λ/n for which 3/3 LLL succeeds, and the m needed.

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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b_param, n, seed=12345):
    rng = random.Random(seed)
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

def find_b_for_n(p, n):
    rng2 = random.Random(99)
    for b_try in range(1, min(p, 500)):
        for _ in range(30):
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
# Curve search by λ/n target range
# ---------------------------------------------------------------------------

def find_curve_by_lamn(lo, hi, p_lo=2**11, p_hi=2**17):
    """
    Find a small prime p with a j=0 curve having n prime, n≡1(mod 3),
    λ/n ∈ [lo, hi].  Returns (p, b, n, lam) or None.
    """
    p = sympy.nextprime(p_lo - 1)
    while p < p_hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for tr in j0_traces(a_e, b_e):
                    n_cand = p + 1 - tr
                    if n_cand < 7:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    if lo <= ratio <= hi:
                        b_param = find_b_for_n(p, n_cand)
                        if b_param is None:
                            continue
                        return (p, b_param, n_cand, lam)
        p = sympy.nextprime(p)
    return None

# ---------------------------------------------------------------------------
# GLV lattice attack
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
        attempts += 1
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(0, K2 - 1)
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
        sigs.append({'A': A, 'B': B})
    return d_secret, sigs

def build_lattice(sigs, n, lam, K1):
    m = len(sigs)
    K2 = math.isqrt(n) + 1
    dim = 2 * m + 2
    S_K1 = max(1, n // K1)
    S_K2 = max(1, n // K2)
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

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_lll(p, b_param, n, lam, G, m, K1, seed, use_bkz=False, bkz_beta=20):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, K1)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep_m(p, b_param, n, lam, G, K1, m_range, seeds, use_bkz=False, bkz_beta=20):
    """Sweep m; return (first_m_with_3of3, {m: (wins, total)})."""
    results = {}
    first_win = None
    for m in m_range:
        wins = sum(1 for s in seeds
                   if run_lll(p, b_param, n, lam, G, m, K1, s,
                              use_bkz=use_bkz, bkz_beta=bkz_beta))
        results[m] = (wins, len(seeds))
        if wins == len(seeds) and first_win is None:
            first_win = m
    return first_win, results

# ---------------------------------------------------------------------------
# Anchor curves (known results from prior sessions)
# ---------------------------------------------------------------------------

ANCHORS = [
    # (label, p, b, n, lam, expected_outcome)
    ("8-bit/fail  lam/n=0.07", 2677, 2, 2647, 185, "FAIL"),
    ("20-bit/pass lam/n=0.34", 524347, None, 523969, 177902, "PASS"),
    ("8-bit/pass  lam/n=0.53", 211,  2, 199,  106, "PASS"),
]

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection (0.07 → 0.34)")
print("=" * 70)
print()

SEEDS = [42, 1234, 9999]

# Targets to bisect: evenly distributed in (0.07, 0.34)
TARGETS = [
    (0.085, 0.095),   # ~0.09
    (0.110, 0.125),   # ~0.12
    (0.140, 0.160),   # ~0.15
    (0.170, 0.195),   # ~0.18
    (0.200, 0.225),   # ~0.21
    (0.255, 0.285),   # ~0.27
]

# ---- Step 1: Find curves for each target -----------------------------------
print("Step 1: Finding j=0 CM curves for each λ/n target window...")
print("-" * 60)

found_curves = []  # (lam_n, p, b, n, lam)
for (lo, hi) in TARGETS:
    result = find_curve_by_lamn(lo, hi)
    if result is None:
        print(f"  λ/n ∈ [{lo:.3f},{hi:.3f}]: NOT FOUND")
        found_curves.append((None, None, None, None, None))
    else:
        p, b, n, lam = result
        ratio = lam / n
        print(f"  λ/n ∈ [{lo:.3f},{hi:.3f}]: p={p}, n={n} ({n.bit_length()}b), "
              f"lam={lam}, lam/n={ratio:.4f}")
        found_curves.append((ratio, p, b, n, lam))

# Also add the known failure anchor (p=2677) and known success (p=524347)
print()
print("Anchor verification: p=2677 (lam/n=0.07, expected FAIL)")
G_anchor_fail = find_generator(2677, 2, 2647)
if G_anchor_fail:
    print(f"  Generator found: {G_anchor_fail}")

print()

# ---- Step 2: LLL sweep for each curve -------------------------------------
print("=" * 70)
print("Step 2: LLL sweep (m=4..20, 3 seeds each)")
print("-" * 70)

K1_BASE = 8   # K1 bound: for small n (~2^12), K1=8 gives eff ≈ K1/sqrt(n) ≈ 0.12

sweep_summary = []   # (lam_n, first_win_lll, first_win_bkz)

for (lam_n, p, b, n, lam) in found_curves:
    if p is None:
        sweep_summary.append((None, None, None))
        continue

    # Adaptive K1: target eff ≈ 0.08..0.12
    K2 = math.isqrt(n) + 1
    K1 = max(2, int(0.10 * math.sqrt(n)))

    G = find_generator(p, b, n)
    if G is None:
        print(f"  lam/n={lam_n:.4f} p={p}: generator not found; SKIP")
        sweep_summary.append((lam_n, None, None))
        continue

    eff = K1 * K2 / n
    print(f"\nlam/n={lam_n:.4f}  p={p}  n={n} ({n.bit_length()}b)  "
          f"lam={lam}  K1={K1}  K2={K2}  eff={eff:.4f}")

    # LLL sweep
    first_lll, res_lll = sweep_m(p, b, n, lam, G, K1,
                                  range(4, 21), SEEDS, use_bkz=False)
    detail = "  ".join(f"m{m}:{w}/{t}" for m, (w, t) in sorted(res_lll.items())
                       if w > 0 or m <= 8)
    if first_lll:
        print(f"  LLL:  3/3 first at m={first_lll}  [{detail}]")
    else:
        max_wins = max(w for w, t in res_lll.values())
        print(f"  LLL:  never 3/3 (best {max_wins}/3)  [{detail}]")

    # BKZ(20) if LLL failed
    first_bkz = None
    if first_lll is None:
        first_bkz, res_bkz = sweep_m(p, b, n, lam, G, K1,
                                      range(5, 16), SEEDS,
                                      use_bkz=True, bkz_beta=20)
        if first_bkz:
            print(f"  BKZ20: 3/3 at m={first_bkz}")
        else:
            max_bkz = max(w for w, t in res_bkz.values())
            print(f"  BKZ20: never 3/3 (best {max_bkz}/3)")

    sweep_summary.append((lam_n, first_lll, first_bkz))

# ---- Step 3: Anchor check (p=2677 failure) ---------------------------------
print()
print("=" * 70)
print("Step 3: Anchor check — p=2677 (lam/n=0.07, EXPECTED FAIL)")
print("-" * 70)
if G_anchor_fail:
    n_anc, lam_anc = 2647, 185
    K2_anc = math.isqrt(n_anc) + 1
    K1_anc = max(2, int(0.10 * math.sqrt(n_anc)))
    eff_anc = K1_anc * K2_anc / n_anc
    print(f"  K1={K1_anc}, K2={K2_anc}, eff={eff_anc:.4f}")
    first_anc, res_anc = sweep_m(2677, 2, n_anc, lam_anc, G_anchor_fail,
                                  K1_anc, range(5, 16), SEEDS)
    if first_anc:
        print(f"  LLL: 3/3 at m={first_anc}  [UNEXPECTED SUCCESS]")
    else:
        max_anc = max(w for w, t in res_anc.values())
        print(f"  LLL: never 3/3 (best {max_anc}/3)  [confirmed FAIL]")
    sweep_summary.append((185/2647, first_anc, None))
else:
    print("  Generator not found; anchor check skipped")
    sweep_summary.append((0.07, None, None))

# ---- Step 4: Summary table -------------------------------------------------
print()
print("=" * 70)
print("SUMMARY: λ/n threshold for GLV-HNP Phase 2 LLL success")
print("-" * 70)
print(f"  {'λ/n':>8}  {'LLL first 3/3':>14}  {'BKZ20 first 3/3':>16}  outcome")
print(f"  {'-'*8}  {'-'*14}  {'-'*16}  {'-'*20}")

# Add anchor known-pass rows to summary
all_rows = []
all_rows.append((0.07, None, None, "FAIL (anchor: p=2677)"))
all_rows.append((0.34, 9, None, "PASS (anchor: p=524347, from 2026-07-26 log)"))
all_rows.append((0.44, 4, None, "PASS (secp256k1 proxy)"))
all_rows.append((0.53, 4, None, "PASS (8-bit p=211)"))
all_rows.append((0.66, 7, None, "PASS (12-bit p=2557)"))

for (lam_n, first_lll, first_bkz) in sweep_summary:
    if lam_n is None:
        continue
    if first_lll is not None:
        outcome = f"PASS at m={first_lll}"
    elif first_bkz is not None:
        outcome = f"LLL-FAIL, BKZ20-PASS at m={first_bkz}"
    else:
        outcome = "FAIL (LLL and BKZ20)"
    all_rows.append((lam_n, first_lll, first_bkz, outcome))

all_rows.sort(key=lambda x: x[0])

for (lam_n, first_lll, first_bkz, outcome) in all_rows:
    lll_str = str(first_lll) if first_lll is not None else "—"
    bkz_str = str(first_bkz) if first_bkz is not None else ("—" if first_lll is not None else "tried")
    print(f"  {lam_n:>8.4f}  {lll_str:>14}  {bkz_str:>16}  {outcome}")

print()

# Identify threshold
pass_ratios = [lam_n for lam_n, fl, fb, _ in all_rows
               if fl is not None or fb is not None]
fail_ratios = [lam_n for lam_n, fl, fb, _ in all_rows
               if fl is None and fb is None]

if pass_ratios and fail_ratios:
    threshold_lo = max(fail_ratios)
    threshold_hi = min(pass_ratios)
    print(f"Threshold: λ/n ∈ ({threshold_lo:.4f}, {threshold_hi:.4f})")
    print(f"LLL succeeds iff λ/n > ~{(threshold_lo + threshold_hi)/2:.3f}")
else:
    print("Threshold not yet determined from this sweep.")

print()
print("Done.")
