"""
GLV-HNP Phase 2: K1_BOUND vs lam/n — focused study.

The 2026-07-26 log concluded that the attack fails when lam/n=0.07 (p=2677, n=2647)
with K1_BOUND=8 and BKZ(40). Today's threshold scan shows LLL succeeds for
ALL lam/n in [0.05, 0.48] when K1_BOUND=2.

This script tests the REAL variable: K1_BOUND (the bias range), not lam/n.

Hypotheses:
  H1: With K1_BOUND=2, the attack succeeds at lam/n=0.07 (p=2677).
      (Would contradict the 2026-07-26 geometric argument.)
  H2: With K1_BOUND=8, the attack fails at lam/n=0.07 (p=2677).
      (Reproduces previous result — confirms K1_BOUND is the differentiator.)
  H3: There is a (K1_BOUND, lam/n) surface — larger K1 needs larger lam/n.

Run: python3 glv_hnp_phase2_k1bound_study.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

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

def count_points(p, b):
    count = 1
    for x in range(p):
        rhs = (pow(x, 3, p) + b) % p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2
    return count

def find_twist_b(p, n_target, b_limit=300):
    for b in range(1, b_limit + 1):
        if count_points(p, b) == n_target:
            return b
    return None

def find_generator(p, b, n):
    rng = random.Random(99991)
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
# Attack machinery
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
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

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KANNAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_attack(curve, m, d_secret, k1_bound, seed, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret)

def sweep(label, curve, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    print(f"\n  {label}")
    print(f"    p={p}, n={n}, lam={lam} (lam/n={lam/n:.4f}), K1={k1_bound}, K2={math.isqrt(n)+1}")
    found_m = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            d = random.Random(seed + 3141).randint(1, n - 1)
            ok = run_attack(curve, m, d, k1_bound, seed, use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        algo = f"BKZ({bkz_beta})" if use_bkz else "LLL"
        print(f"    m={m}: {wins}/{len(seeds)} [{algo}]", end='')
        if wins == len(seeds):
            print(" ← 3/3 ✓")
            found_m = m
            break
        else:
            print()
    if found_m is None:
        print(f"    → FAIL up to m={max(m_range)-1}")
    return found_m

SEEDS = [42, 1234, 9999]
M_SHORT = range(3, 16)
M_LONG  = range(3, 26)

# ===========================================================================
# Part 1: Reproduce the 2026-07-26 failure (p=2677, n=2647, lam=185, K1=8)
# ===========================================================================
print("=" * 65)
print("Part 1: Reproduce 2026-07-26 failure (p=2677, lam/n≈0.07, K1=8)")
print("=" * 65)

p_fail, n_fail, lam_fail = 2677, 2647, 185
assert (lam_fail**2 + lam_fail + 1) % n_fail == 0, "lam check"
print(f"\nFinding twist b for p={p_fail}, n={n_fail}...")
b_fail = find_twist_b(p_fail, n_fail)
print(f"  b={b_fail}")
G_fail = find_generator(p_fail, b_fail, n_fail)
print(f"  G={G_fail}")
curve_fail = (p_fail, b_fail, n_fail, lam_fail, G_fail)

# K1=8 (previous session)
sweep("K1_BOUND=8 [PREVIOUS FAILURE — expected to fail]",
      curve_fail, k1_bound=8, m_range=M_LONG, seeds=SEEDS)

# ===========================================================================
# Part 2: Test same curve with K1=2 (H1: should succeed)
# ===========================================================================
print("\n" + "=" * 65)
print("Part 2: Same curve (p=2677, lam/n≈0.07) with K1_BOUND=2")
print("=" * 65)

sweep("K1_BOUND=2 [HYPOTHESIS: will succeed if K1, not lam/n, is the variable]",
      curve_fail, k1_bound=2, m_range=M_SHORT, seeds=SEEDS)

# ===========================================================================
# Part 3: K1_BOUND sweep on the lam/n≈0.07 curve
# ===========================================================================
print("\n" + "=" * 65)
print("Part 3: K1_BOUND sweep on p=2677, lam/n≈0.07")
print("  Finds the max K1_BOUND at which LLL still succeeds.")
print("=" * 65)

for k1b in [2, 3, 4, 5, 6, 8, 10]:
    sweep(f"K1_BOUND={k1b}", curve_fail, k1_bound=k1b,
          m_range=range(3, 21), seeds=SEEDS)

# ===========================================================================
# Part 4: K1_BOUND sweep on higher lam/n curves — does lam/n change the K1 limit?
# ===========================================================================
print("\n" + "=" * 65)
print("Part 4: K1_BOUND sweep on lam/n≈0.27 curve (p=163, n=181, lam=48)")
print("  Checks whether larger lam/n tolerates larger K1_BOUND.")
print("=" * 65)

p_hi, n_hi, lam_hi = 163, 181, 48
assert (lam_hi**2 + lam_hi + 1) % n_hi == 0
b_hi = find_twist_b(p_hi, n_hi)
print(f"\n  b={b_hi}")
G_hi = find_generator(p_hi, b_hi, n_hi)
curve_hi = (p_hi, b_hi, n_hi, lam_hi, G_hi)

for k1b in [2, 4, 6, 8, 12, 16]:
    sweep(f"K1_BOUND={k1b}", curve_hi, k1_bound=k1b,
          m_range=range(3, 21), seeds=SEEDS)

# ===========================================================================
# Part 5: BKZ rescue on K1=8, lam/n≈0.07 (from 2026-07-26)
# ===========================================================================
print("\n" + "=" * 65)
print("Part 5: BKZ rescue on K1=8, lam/n≈0.07 — push m further (3..30)")
print("=" * 65)

for bkz_b in [40, 60]:
    sweep(f"BKZ(beta={bkz_b}), K1=8",
          curve_fail, k1_bound=8, m_range=range(3, 31), seeds=SEEDS,
          use_bkz=True, bkz_beta=bkz_b)

# ===========================================================================
# Summary
# ===========================================================================
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print("""
The 2026-07-26 log concluded: lam/n < 0.07 causes structural failure.
TODAY'S FINDING: With K1_BOUND=2, LLL succeeds at ALL lam/n tested (0.05..0.48).

Working hypothesis:
  The attack's K1_BOUND (bias range) is the real differentiator, not lam/n.
  Larger K1_BOUND → planted vector norm grows → weaker SNR → more sigs needed.
  The lam/n ratio affects HOW MANY sigs (m) are needed, but not viability.

See Part 3 above for the critical K1_BOUND threshold on the lam/n=0.07 curve.
""")
