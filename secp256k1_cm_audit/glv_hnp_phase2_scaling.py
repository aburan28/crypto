"""
GLV-aware HNP Phase 2 — scaling test on 12-bit and 14-bit toy curves.

Extends glv_hnp_phase2_attack.py to verify the attack generalises beyond
the original 8-bit toy (n=199).  Tests curves found by ellcard enumeration:
  - 8-bit:  y^2 = x^3 + 2 over F_211,  n=199,   lam=106
  - 12-bit: y^2 = x^3 + 2 over F_2557, n=2659,  lam=1755
  - 12-bit: y^2 = x^3 + 2 over F_2677, n=2647,  lam=185
  - 12-bit: y^2 = x^3 + 2 over F_3571, n=3571,  lam=3467

For each curve, sets K2_BOUND = ceil(sqrt(n)) and sweeps K1_BOUND
to demonstrate the threshold where LLL reliably recovers d.

Run: python3 glv_hnp_phase2_scaling.py
"""

import random
import math
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass a=0)
# ---------------------------------------------------------------------------
def modinv(a, m): return pow(a, -1, m)

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
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def find_generator(p, b, n):
    for x in range(p):
        rhs = (x**3 + b) % p
        for y in range(p):
            if y * y % p == rhs:
                G = (x, y)
                if ec_mul(G, n, p) is None:
                    return G
    return None

# ---------------------------------------------------------------------------
# GLV decomposition: brute force for small n
# ---------------------------------------------------------------------------
def glv_decompose(k_full, n, lam, k2_bound):
    for k2 in range(k2_bound):
        k1 = (k_full - lam * k2) % n
        if k1 < 50:  # small k1 signals decomposition found
            return (k1, k2)
    return None

# ---------------------------------------------------------------------------
# Signature generation
# ---------------------------------------------------------------------------
def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 100000:
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
# Build GLV lattice with column-diagonal scaling
# ---------------------------------------------------------------------------
def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D  = 1
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

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------
def recover_d(M_reduced, m, n, S_KANNAN, d_secret=None):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_secret is not None and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Single experiment
# ---------------------------------------------------------------------------
def run_experiment(curve_params, m, d_secret, k1_bound, seed=42):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1

    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False

    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2

    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]

    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
def sweep_curve(label, curve_params, k1_bound, m_range, seeds, verbose=True):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    if eff >= 1.0:
        m_thresh = float('inf')
    else:
        m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))

    if verbose:
        print(f"\n{'='*60}")
        print(f"Curve: {label}")
        print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}")
        print(f"  K1_BOUND={k1_bound}, K2_BOUND={k2_bound}")
        print(f"  eff = K1*K2/n = {eff:.4f}, m_thresh ≈ {m_thresh:.1f}")
        print(f"{'='*60}")

    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok = run_experiment(curve_params, m, d_trial, k1_bound, seed)
            wins += ok
        results[m] = (wins, len(seeds))
        marker = " (≥thresh)" if m >= m_thresh else " (below)"
        if verbose:
            print(f"  m={m}: {wins}/{len(seeds)} recovered{marker}")
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
print("=" * 60)
print("GLV-aware HNP Phase 2 — scaling test")
print("=" * 60)

SEEDS = [42, 1234, 9999, 7, 314159]

# Curve 1: original 8-bit toy (p=211, n=199, lam=106)
print("\nBuilding generator for 8-bit curve (p=211)...")
G0 = find_generator(211, 2, 199)
print(f"  G = {G0}")
curve0 = (211, 2, 199, 106, G0)
res0 = sweep_curve("8-bit: y²=x³+2/F_211, n=199, lam=106",
                   curve0, k1_bound=2, m_range=range(3, 8), seeds=SEEDS)

# Curve 2: 12-bit toy (p=2557, n=2659, lam=1755)
print("\nBuilding generator for 12-bit curve (p=2557)...")
G1 = find_generator(2557, 2, 2659)
print(f"  G = {G1}")
curve1 = (2557, 2, 2659, 1755, G1)
# K1_BOUND = 8: k1 in [0,8), about 3 bits of bias out of ~6 bits (sqrt(2659)≈51)
res1 = sweep_curve("12-bit: y²=x³+2/F_2557, n=2659, lam=1755",
                   curve1, k1_bound=8, m_range=range(3, 9), seeds=SEEDS)

# Curve 3: 12-bit (p=2677, n=2647, lam=185)
print("\nBuilding generator for 12-bit curve (p=2677)...")
G2 = find_generator(2677, 2, 2647)
print(f"  G = {G2}")
curve2 = (2677, 2, 2647, 185, G2)
res2 = sweep_curve("12-bit: y²=x³+2/F_2677, n=2647, lam=185",
                   curve2, k1_bound=8, m_range=range(3, 9), seeds=SEEDS)

# Summary
print("\n" + "=" * 60)
print("Summary: first m where 5/5 seeds recovered")
print("=" * 60)
for label, results in [("8-bit (K1=2)", res0), ("12-bit/2557 (K1=8)", res1), ("12-bit/2677 (K1=8)", res2)]:
    for m, (wins, total) in sorted(results.items()):
        if wins == total:
            print(f"  {label}: first 5/5 at m={m}")
            break
    else:
        print(f"  {label}: never reached 5/5 in sweep")
