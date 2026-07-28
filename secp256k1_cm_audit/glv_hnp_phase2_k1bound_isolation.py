"""
GLV-HNP Phase 2: K1_BOUND isolation test.

The 2026-07-26 session reported lam/n=0.07 (p=2677) as a "structural failure"
even with BKZ-40. But the lambda-threshold run today (2026-07-28) found success
down to lam/n=0.042 using K1_BOUND=2.

Hypothesis: the failure at p=2677 was NOT structural (not caused by small lam/n)
but was instead caused by LARGE K1_BOUND=8, which makes k_full harder to recover.

Test: sweep K1_BOUND ∈ {2,4,6,8,10,12} on two curves:
  - p=2677, n=2647, lam=185 (lam/n=0.070)  — the prior failure at K1_BOUND=8
  - p=547,  n=547,  lam=40  (lam/n=0.073)  — this session's success at K1_BOUND=2

For each (curve, K1_BOUND), find min m for 3/3 LLL success (or FAIL if >18).

Run: python3 glv_hnp_phase2_k1bound_isolation.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (copied from threshold script)
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
    if k < 0: return ec_mul((P[0], (-P[1]) % p), -k, p)
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

def find_generator(p, b, n):
    rng = random.Random(1234)
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
# Signature & lattice
# ---------------------------------------------------------------------------

def gen_sigs(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed):
    rng = random.Random(seed)
    sigs = []
    for _ in range(500000):
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})
        if len(sigs) == m: break
    return sigs

def run_lll(p, b, n, lam, G, m, k1_bound, d, seed):
    k2_bound = math.isqrt(n) + 1
    sigs = gen_sigs(G, d, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m: return False
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
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    for row_idx in range(dim):
        last = A[row_idx][dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[row_idx][m]) % n
        if d_cand == d: return True
    return False

SEEDS = [42, 1234, 9999]

def min_m(p, b, n, lam, G, k1_bound, m_max=22):
    for m in range(3, m_max + 1):
        wins = sum(
            run_lll(p, b, n, lam, G, m, k1_bound,
                    random.Random(s * 31337).randint(1, n - 1), s)
            for s in SEEDS
        )
        if wins == len(SEEDS):
            return m
    return None

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 68)
print("K1_BOUND isolation: does failure at lam/n=0.07 depend on K1_BOUND?")
print("=" * 68)

# Curve A: the prior "structural failure" (p=2677, K1_BOUND=8)
p_A, b_A, n_A, lam_A = 2677, 2, 2647, 185
G_A = find_generator(p_A, b_A, n_A)
print(f"\nCurve A: p={p_A}, b={b_A}, n={n_A}, lam={lam_A}, lam/n={lam_A/n_A:.4f}")
print(f"  (Prior result: FAIL at K1_BOUND=8, m≤12, BKZ-40 also fails)")
if G_A is None:
    print("  ERROR: could not find generator")
else:
    print(f"  G={G_A}")
    print(f"\n  K1_BOUND sweep (K2_BOUND={math.isqrt(n_A)+1}):")
    print(f"  {'K1_BOUND':>10}  {'eff=K1*K2/n':>12}  {'m*':>4}")
    print("  " + "-" * 35)
    for k1 in [2, 3, 4, 6, 8, 10, 12]:
        eff = k1 * (math.isqrt(n_A) + 1) / n_A
        m_star = min_m(p_A, b_A, n_A, lam_A, G_A, k1, m_max=22)
        s = f"{m_star}" if m_star else "FAIL(>22)"
        print(f"  {k1:>10d}  {eff:>12.5f}  {s:>4}")

# Curve B: today's success (p=547, n=547, lam=40)
p_B, b_B, n_B, lam_B = 547, 1, 547, 40
G_B = find_generator(p_B, b_B, n_B)
print(f"\nCurve B: p={p_B}, b={b_B}, n={n_B}, lam={lam_B}, lam/n={lam_B/n_B:.4f}")
print(f"  (This session: SUCCESS at K1_BOUND=2, m*=5)")
if G_B is None:
    # Try other b values
    for b_try in range(1, 100):
        G_B = find_generator(p_B, b_try, n_B)
        if G_B:
            b_B = b_try
            print(f"  (found generator at b={b_B})")
            break
if G_B is None:
    print("  ERROR: could not find generator")
else:
    print(f"  G={G_B}")
    print(f"\n  K1_BOUND sweep (K2_BOUND={math.isqrt(n_B)+1}):")
    print(f"  {'K1_BOUND':>10}  {'eff=K1*K2/n':>12}  {'m*':>4}")
    print("  " + "-" * 35)
    for k1 in [2, 3, 4, 6, 8, 10, 12]:
        eff = k1 * (math.isqrt(n_B) + 1) / n_B
        m_star = min_m(p_B, b_B, n_B, lam_B, G_B, k1, m_max=22)
        s = f"{m_star}" if m_star else "FAIL(>22)"
        print(f"  {k1:>10d}  {eff:>12.5f}  {s:>4}")

print("\nDone.")
