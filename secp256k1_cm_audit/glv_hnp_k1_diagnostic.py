"""
GLV-HNP Phase 2: K1 diagnostic — vary K1 at fixed 20-bit curve (2026-06-17).

Previous: controlled experiment showed 20-bit eff=0.15 (K1=109) fails at m=16.
Prior session: 20-bit eff=0.05 (K1=36) succeeds at m=9.

Question: what is the K1 threshold at 20-bit?  Is the failure structural or
does it just need larger m?

Diagnostic A: K1 scan (36, 55, 72, 90, 109) at 20-bit, sweep m=5..18.
Diagnostic B: extended sweep for K1=109 up to m=22.

Run: python3 glv_hnp_k1_diagnostic.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (same as other scripts)
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
# GLV-HNP experiment
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
        s_inv = modinv(s, n)
        sigs.append({
            'A': h * s_inv % n,
            'B': r * s_inv % n,
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

def experiment(curve, m, d_secret, k1_bound, seed):
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

def sweep(label, curve, k1_bound, m_range, seeds):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = math.ceil(math.log(n) / (-math.log(eff)))
    print(f"  [{label}]  K1={k1_bound}, eff={eff:.4f}, m_thresh={m_thresh}")
    first_full = None
    for m in m_range:
        wins = sum(experiment(curve, m, random.Random(s+9001).randint(1, n-1), k1_bound, s)
                   for s in seeds)
        marker = "←" if m == m_thresh else " "
        print(f"    m={m:2d}: {wins}/{len(seeds)} {marker}")
        if wins == len(seeds) and first_full is None:
            first_full = m
    ratio = f"{first_full/m_thresh:.2f}" if first_full else "N/A"
    print(f"    → first 3/3: m={first_full}, ratio={ratio}\n")
    return first_full, m_thresh, eff

# ---------------------------------------------------------------------------
# Curve: 20-bit y²=x³+2 / F_524347, n=523969, lam=177902  (from 2026-06-16)
# ---------------------------------------------------------------------------

P20 = 524347
B20 = 2
N20 = 523969
LAM20 = 177902
K2_20 = math.isqrt(N20) + 1  # 724

print("Building 20-bit generator...")
G20 = find_generator(P20, B20, N20)
assert G20 is not None, "Generator not found"
assert ec_mul(G20, N20, P20) is None, "Order check failed"
print(f"  G20 = {G20}")
CURVE20 = (P20, B20, N20, LAM20, G20)

SEEDS = [42, 1234, 9999]

print(f"\n20-bit curve: p={P20}, n={N20} ({N20.bit_length()}b), lam={LAM20}, lam/n={LAM20/N20:.4f}")
print(f"K2={K2_20}")
print("=" * 55)

# Diagnostic A: vary K1 (eff from 0.05 to 0.15)
print("\n=== DIAGNOSTIC A: K1 scan (m=5..18) ===")
k1_values = [36, 55, 72, 90, 109]
results_a = {}
for k1 in k1_values:
    eff_expected = k1 * K2_20 / N20
    f3, mt, eff = sweep(f"K1={k1:3d} eff={eff_expected:.4f}", CURVE20, k1, range(5, 19), SEEDS)
    results_a[k1] = (f3, mt, eff)

# Diagnostic B: extended sweep for K1=109 (the failed case)
print("\n=== DIAGNOSTIC B: K1=109 extended sweep (m=14..24) ===")
f3_ext, mt_ext, eff_ext = sweep("K1=109 extended", CURVE20, 109, range(14, 25), SEEDS)

# Summary table
print("=" * 55)
print("SUMMARY: 20-bit, vary K1")
print(f"{'K1':>5} {'eff':>7} {'m_thresh':>9} {'first_3/3':>10} {'ratio':>7}")
for k1, (f3, mt, eff) in results_a.items():
    ratio = f"{f3/mt:.2f}" if f3 else "N/A"
    m_str = str(f3) if f3 else "—"
    print(f"  {k1:3d}  {eff:7.4f}  {mt:9d}  {m_str:>10}  {ratio:>7}")
# K1=109 extended
ext_f3 = f3_ext if f3_ext else "—"
ext_ratio = f"{f3_ext/mt_ext:.2f}" if f3_ext else "N/A"
print(f"  109  {eff_ext:7.4f}  {mt_ext:9d}  {str(ext_f3):>10}  {ext_ratio:>7}  (ext sweep)")
print("\nDone.")
