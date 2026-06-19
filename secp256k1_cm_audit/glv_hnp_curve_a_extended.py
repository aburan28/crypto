"""
GLV-HNP Phase 2: Curve A extended sweep (2026-06-19).

Yesterday (2026-06-18) falsified the lam/n hypothesis:
  - Curve A (lam/n=0.3395, p=524347): fails at K1≥55 (eff≥0.076) in m=5..19
  - Curve C (lam/n=0.2867, p=624517): succeeds at all K1 up to 109 in m≤13
  - lam/n alone does NOT predict the eff ceiling.

Two candidate explanations for Curve A's K1=55 failure:
  (a) Insufficient m — Curve A just needs m>19 (ratio>3.17); should eventually succeed.
  (b) Structural lattice obstruction — lam≈n/3 creates a near-singular lattice
      structure that LLL cannot overcome regardless of m.

This script resolves (a) vs (b) by:
  1. Sweeping Curve A at K1=55 for m=5..35 with 10 seeds.
     If it succeeds at some m≤35, explanation (a) is confirmed.
     If it still fails at m=35, explanation (b) is more likely (or m>>35 needed).
  2. Checking K1=36 and K1=55 on Curve A with BKZ(β=20) at m=15..25 to see
     if stronger reduction helps where LLL fails.
  3. Running 10-seed check for K1=36 on Curve A (previously only 3 seeds)
     to verify robustness of the m=7 result.

Curve A parameters:
  p=524347, b=2, n=523969, lam=177902
  lam/n = 0.3395, 3*lam = 533706, n=523969, δ(3*lam, n) = 9737
  Note: lam ≈ n/3 to within 9737/523969 ≈ 1.9%

Run: python3 glv_hnp_curve_a_extended.py
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (same as lamn_isolation.py)
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
            pt = (x, y)
            if ec_mul(pt, n, p) is None:
                return pt
    return None

# ---------------------------------------------------------------------------
# GLV-HNP lattice construction (same as lamn_isolation.py)
# ---------------------------------------------------------------------------

def gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed):
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
        sigs.append({
            'A': h * modinv(s, n) % n,
            'B': r * modinv(s, n) % n,
            'k1': k1, 'k2': k2,
        })
    return sigs

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
        if abs(row[dim-1]) != S_KANNAN: continue
        sign = 1 if row[dim-1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def experiment_lll(curve, m, d_secret, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def experiment_bkz(curve, m, d_secret, k1_bound, seed, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_sigs(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    # Use empty strategies list to avoid missing strategies file
    par = BKZ.Param(block_size=bkz_beta, strategies=[], auto_abort=True)
    BKZ.reduction(A, par)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

# ---------------------------------------------------------------------------
# Setup: Curve A
# ---------------------------------------------------------------------------

P_A, B_A, N_A, LAM_A = 524347, 2, 523969, 177902
K2_A = math.isqrt(N_A) + 1

print("=" * 70)
print("GLV-HNP Phase 2: Curve A extended sweep (2026-06-19)")
print("=" * 70)
print(f"Curve A: p={P_A}, b={B_A}, n={N_A}, lam={LAM_A}")
print(f"  lam/n = {LAM_A/N_A:.4f}")
print(f"  3*lam mod n = {3*LAM_A % N_A}  (δ from 0: {3*LAM_A % N_A}, want ≠ 0)")
print(f"  K2 = {K2_A}")
print()

G_A = find_generator(P_A, B_A, N_A)
assert G_A is not None and ec_mul(G_A, N_A, P_A) is None
CURVE_A = (P_A, B_A, N_A, LAM_A, G_A)
print("Generator OK.")
print()

SEEDS_10 = [42, 1234, 9999, 7, 31337, 999999, 0xDEAD, 0xBEEF, 0xC0DE, 0xCAFE]

# ---------------------------------------------------------------------------
# Experiment 1: K1=36, 10 seeds, m=5..20
# (Sanity check: verify K1=36 is robust, previously m=7 with 3 seeds)
# ---------------------------------------------------------------------------

K1_36 = 36
eff_36 = K1_36 * K2_A / N_A
m_thresh_36 = math.ceil(math.log(N_A) / (-math.log(eff_36)))

print(f"--- Experiment 1: K1={K1_36}, eff={eff_36:.4f}, m_thresh={m_thresh_36} ---")
print(f"    10 seeds, m=5..20, LLL")
print()

first_full_36 = None
for m in range(5, 21):
    wins = sum(
        experiment_lll(CURVE_A, m, random.Random(s + 9001).randint(1, N_A - 1), K1_36, s)
        for s in SEEDS_10
    )
    marker = " ← m_thresh" if m == m_thresh_36 else ""
    print(f"  m={m:2d}: {wins:2d}/10{marker}")
    if wins == 10 and first_full_36 is None:
        first_full_36 = m
        print(f"  *** 10/10 first achieved at m={m} ***")

ratio_36 = f"{first_full_36/m_thresh_36:.2f}" if first_full_36 else "N/A"
print(f"  K1={K1_36}: first 10/10 at m={first_full_36}, ratio={ratio_36}")
print()

# ---------------------------------------------------------------------------
# Experiment 2: K1=55, 10 seeds, m=5..35 (the main question)
# ---------------------------------------------------------------------------

K1_55 = 55
eff_55 = K1_55 * K2_A / N_A
m_thresh_55 = math.ceil(math.log(N_A) / (-math.log(eff_55)))

print(f"--- Experiment 2: K1={K1_55}, eff={eff_55:.4f}, m_thresh={m_thresh_55} ---")
print(f"    10 seeds, m=5..35, LLL")
print(f"    Resolves: is Curve A failure at K1=55 structural or just needs m>19?")
print()

t0 = time.time()
first_full_55 = None
max_wins_55 = 0
for m in range(5, 36):
    wins = sum(
        experiment_lll(CURVE_A, m, random.Random(s + 9001).randint(1, N_A - 1), K1_55, s)
        for s in SEEDS_10
    )
    marker = " ← m_thresh" if m == m_thresh_55 else ""
    print(f"  m={m:2d}: {wins:2d}/10{marker}")
    if wins > max_wins_55:
        max_wins_55 = wins
    if wins == 10 and first_full_55 is None:
        first_full_55 = m
        print(f"  *** 10/10 first achieved at m={m} ***")
    elif wins == 10:
        pass

t1 = time.time()
ratio_55 = f"{first_full_55/m_thresh_55:.2f}" if first_full_55 else "N/A"
print(f"  K1={K1_55}: first 10/10 at m={first_full_55}, ratio={ratio_55}, max_wins={max_wins_55}/10")
print(f"  Wall time: {t1-t0:.1f}s")
print()

# ---------------------------------------------------------------------------
# Experiment 3: BKZ(β=20) on Curve A at K1=55, m=10..20
# Does BKZ rescue what LLL cannot?
# ---------------------------------------------------------------------------

print(f"--- Experiment 3: K1={K1_55}, BKZ(β=20), 3 seeds, m=10..22 ---")
print()

first_full_bkz = None
for m in range(10, 23):
    wins = sum(
        experiment_bkz(CURVE_A, m, random.Random(s + 9001).randint(1, N_A - 1), K1_55, s, bkz_beta=20)
        for s in SEEDS_10[:3]
    )
    print(f"  m={m:2d}: {wins}/3 (BKZ-20)")
    if wins == 3 and first_full_bkz is None:
        first_full_bkz = m

print(f"  BKZ-20 first 3/3 at m={first_full_bkz}")
print()

# ---------------------------------------------------------------------------
# Experiment 4: K1=36 BKZ rescue vs LLL comparison (sanity)
# ---------------------------------------------------------------------------

print(f"--- Experiment 4: K1={K1_36}, BKZ(β=20) vs LLL, 3 seeds, m=5..12 ---")
print()

for m in range(5, 13):
    wins_lll = sum(
        experiment_lll(CURVE_A, m, random.Random(s + 9001).randint(1, N_A - 1), K1_36, s)
        for s in SEEDS_10[:3]
    )
    wins_bkz = sum(
        experiment_bkz(CURVE_A, m, random.Random(s + 9001).randint(1, N_A - 1), K1_36, s, bkz_beta=20)
        for s in SEEDS_10[:3]
    )
    print(f"  m={m:2d}: LLL={wins_lll}/3  BKZ-20={wins_bkz}/3")

print()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Curve A: p={P_A}, n={N_A}, lam/n={LAM_A/N_A:.4f}")
print(f"  lam ≈ n/3? 3*lam={3*LAM_A}, n={N_A}, δ={3*LAM_A - N_A}")
print()
print(f"K1=36 (eff={eff_36:.4f}): first 10/10 at m={first_full_36} (ratio={ratio_36})")
print(f"K1=55 (eff={eff_55:.4f}): first 10/10 at m={first_full_55} (ratio={ratio_55})")
print()

if first_full_55 is not None:
    print(f"CONCLUSION: Curve A DOES succeed at K1=55 given enough m (m={first_full_55}).")
    print(f"  The K1=55 failure at m≤19 is an insufficient-m issue, NOT structural.")
    print(f"  Explanation (a) confirmed: just needs m>{first_full_55-1}.")
elif max_wins_55 >= 3:
    print(f"CONCLUSION: Curve A shows partial success (max {max_wins_55}/10 at K1=55).")
    print(f"  Suggests m>35 might succeed, but very high overhead required.")
else:
    print(f"CONCLUSION: Curve A NEVER succeeded at K1=55 in m=5..35 (max {max_wins_55}/10).")
    print(f"  Consistent with explanation (b): lam≈n/3 creates a structural obstruction.")
    print(f"  Near-singular sub-lattice hypothesis: lam≈n/3 means the GLV row")
    print(f"  M[m+1+i][i] = -lam*S_K1 ≈ -(n/3)*S_K1 and M[i][i] = n*S_K1,")
    print(f"  so the ratio is -1/3, creating repeated near-linear dependence across rows.")

print()
print(f"BKZ-20 rescue at K1=55: first 3/3 at m={first_full_bkz}")
if first_full_bkz is not None and (first_full_55 is None or first_full_bkz < first_full_55):
    print("  BKZ-20 rescues Curve A at K1=55 at lower m than LLL.")
print()
print("Done.")
