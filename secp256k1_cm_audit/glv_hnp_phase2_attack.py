"""
GLV-aware HNP toy attack (Phase 2) — Python/fpylll implementation.

Threat model: GLV-domain k1-bias.
  k = k1 + lambda * k2 mod n  (GLV decomposition)
  k1 in [0, K1_BOUND)          (biased: attacker knows k1 is small)
  k2 in [0, K2_BOUND)          (GLV domain: |k2| <= sqrt(n))

Toy curve: y^2 = x^3 + 2 over F_211, n=199, lambda=106.
  (same as glv_hnp_phase2_lattice.gp)

Lattice: dimension 2m+2, with column-diagonal scaling so that
all entries of the expected short vector are O(n * K1 / K2).

Run: python3 glv_hnp_phase2_attack.py
"""

import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Curve parameters (from PARI: y^2 = x^3 + 2 over F_211, n=199, lambda=106)
# ---------------------------------------------------------------------------
P_FIELD = 211
B_COEFF = 2
N = 199
LAM = 106        # GLV eigenvalue: LAM^2 + LAM + 1 = 0 mod N
LAM_INV = pow(LAM, -1, N)   # lambda^{-1} mod n

# ---------------------------------------------------------------------------
# Minimal EC arithmetic over F_p (a=0, y^2 = x^3 + b)
# ---------------------------------------------------------------------------
def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P;  x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % P_FIELD == 0:
            return None  # point at infinity
        s = 3 * x1 * x1 * modinv(2 * y1, P_FIELD) % P_FIELD
    else:
        s = (y2 - y1) * modinv(x2 - x1, P_FIELD) % P_FIELD
    x3 = (s * s - x1 - x2) % P_FIELD
    y3 = (s * (x1 - x3) - y1) % P_FIELD
    return (x3, y3)

def ec_mul(P, k):
    R, Q = None, P
    while k > 0:
        if k & 1:
            R = ec_add(R, Q)
        Q = ec_add(Q, Q)
        k >>= 1
    return R

# Find a generator (first point of order N)
G = None
for _x in range(P_FIELD):
    rhs = (_x**3 + B_COEFF) % P_FIELD
    for _y in range(P_FIELD):
        if _y * _y % P_FIELD == rhs:
            _P = (_x, _y)
            if ec_mul(_P, N) is None:
                G = _P
                break
    if G:
        break
assert G is not None, "Generator not found"

# Verify GLV eigenvalue: G * (lambda) should equal phi(G) = (beta*x, y)
# phi: (x,y) -> (beta*x, y) where beta = cube root of unity in F_p
# We just verify lambda^2 + lambda + 1 == 0 mod n
assert (LAM * LAM + LAM + 1) % N == 0, "Bad lambda"

# ---------------------------------------------------------------------------
# GLV decomposition: given k_full in [1,n), return (k1, k2) with
#   k_full = k1 + LAM * k2  mod n,  0 <= k1 < K2_BOUND,  0 <= k2 < K2_BOUND
# For tiny n=199 we brute-force the search over k2 in [0, K2_BOUND).
# ---------------------------------------------------------------------------
K2_BOUND = 15     # ~= ceil(sqrt(199))
K1_BOUND = 2      # bias bound: k1 in {0, 1}

def glv_decompose(k_full, k2_bound=K2_BOUND):
    """Return (k1, k2) with k1 + LAM*k2 == k_full (mod N), k2 in [0, k2_bound).
    Returns None if no such (k1, k2) with k1 in [0, K1_BOUND) exists."""
    for k2 in range(k2_bound):
        k1 = (k_full - LAM * k2) % N
        if 0 <= k1 < K1_BOUND:
            return (k1, k2)
    return None

# ---------------------------------------------------------------------------
# Generate signatures with k in GLV domain and k1 biased
# ---------------------------------------------------------------------------
def gen_signatures(d_secret, m, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 100000:
        attempts += 1
        # Sample k in GLV domain: k1 in [0,K1_BOUND), k2 in [0,K2_BOUND)
        k1 = rng.randint(0, K1_BOUND - 1)
        k2 = rng.randint(0, K2_BOUND - 1)
        k_full = (k1 + LAM * k2) % N
        if k_full == 0:
            continue
        R = ec_mul(G, k_full)
        if R is None:
            continue
        r = R[0] % N
        if r == 0:
            continue
        h = rng.randint(0, N - 1)
        s = modinv(k_full, N) * (h + d_secret * r) % N
        if s == 0:
            continue
        s_inv = modinv(s, N)
        A = h * s_inv % N
        B = r * s_inv % N
        # Verify: A + B*d == k_full mod n
        assert (A + B * d_secret) % N == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Build the 2D-GLV lattice with COLUMN-DIAGONAL BALANCING.
#
# Standard rows (before column scaling):
#   Rows 0..m-1:  n * e_i               [mod-n constraint on k1_i slot]
#   Row m:        (B_1,...,B_m, 1, 0,...,0, 0)   [d-row; d-diagonal = 1]
#   Row m+1+i:    (-lam * e_i, 0, 1 in k2_i slot, 0)  [k2 rows; k2-diagonal = 1]
#   Row 2m+1:     (A_1,...,A_m, 0,...,0, 1)  [Kannan; Kannan-diagonal = 1]
#
# After column-scaling by diag(S_k1,...,S_k1, S_d, S_k2,...,S_k2, S_kannan):
#   k1_i slot j: multiply column j by S_k1
#   d slot m:    multiply column m by S_d
#   k2 slot m+1+i: multiply column m+1+i by S_k2
#   Kannan slot: multiply last column by S_kannan
#
# Target (planted) short vector entries become:
#   k1_i * S_k1,  d * S_d,  k2_i * S_k2,  S_kannan
#
# For BALANCE: all entries should have the same magnitude T.
#   k1_i * S_k1 ≈ T  →  S_k1 ≈ T / K1_BOUND
#   d    * S_d  ≈ T  →  S_d  ≈ T / N
#   k2_i * S_k2 ≈ T  →  S_k2 ≈ T / K2_BOUND
#   S_kannan        ≈ T
#
# Setting T = N (the modulus): S_k1 = N // K1_BOUND, S_d = 1, S_k2 = N // K2_BOUND, S_kannan = N.
#
# Spurious check: combination 2*(k2-row_i) + (mod-n-row_i) produces
#   k1_i slot: (2*(-lam) + n) * S_k1 = (-13) * S_k1
#   k2_i slot: 2 * S_k2
#   All other slots: 0
# Spurious norm = sqrt((-13*S_k1)^2 + (2*S_k2)^2)
#               = sqrt(169 * (N/K1)^2 + 4 * (N/K2)^2)
#   For N=199, K1=2, K2=15: sqrt(169*9900^2/... wait let me compute.
#
# With S_k1 = 199//2 = 99, S_k2 = 199//15 = 13, S_d = 1, S_kannan = 199:
# Spurious norm = sqrt(169*99^2 + 4*13^2) = sqrt(169*9801 + 4*169) = sqrt(169*(9801+4)) = 13*sqrt(9805) ≈ 13*99 = 1287.
# Planted norm = sqrt(m*(K1*S_k1)^2 + (N*S_d)^2 + m*(K2*S_k2)^2 + S_kannan^2)
#              = sqrt(4*(2*99)^2 + 199^2 + 4*(14*13)^2 + 199^2)
#              = sqrt(4*39204 + 39601 + 4*33124 + 39601)
#              = sqrt(156816 + 79202 + 132496) = sqrt(368514) ≈ 607.
# Here planted (607) < spurious (1287). GOOD!
# ---------------------------------------------------------------------------

# Column scale factors
S_K1     = N // K1_BOUND     # 199 // 2 = 99
S_D      = 1                 # d appears unscaled (d ≤ n ≈ 200)
S_K2     = N // K2_BOUND     # 199 // 15 = 13
S_KANNAN = N                 # Kannan scale = N

def build_glv_lattice(sigs):
    m = len(sigs)
    dim = 2 * m + 2
    M = [[0] * dim for _ in range(dim)]

    # Rows 0..m-1: mod-n constraints (scaled k1_i column by S_K1)
    for i in range(m):
        M[i][i] = N * S_K1

    # Row m: d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1   # B_i in k1_i column (scaled)
    M[m][m] = S_D                         # d-diagonal (S_D = 1)

    # Rows m+1..2m: k2_i rows
    for i in range(m):
        M[m + 1 + i][i] = -LAM * S_K1   # -lambda in k1_i column (scaled)
        M[m + 1 + i][m + 1 + i] = S_K2  # k2_i diagonal

    # Row 2m+1: Kannan target row
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1  # A_i in k1 column (scaled)
    M[2 * m + 1][dim - 1] = S_KANNAN

    return M

def matrix_to_fpylll(M):
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    return A

# ---------------------------------------------------------------------------
# Verify that the planted combination is in the lattice with CORRECT signs.
# Correct coefficients:
#   -q_i * row_i  +  d * row_m  +  k2_i * row_{m+1+i}  +  1 * row_{2m+1}
# gives the vector with k1_i slots = k1_i * S_K1.
#
# Expected short vector (after column scaling):
#   slot i (k1_i):     k1_i * S_K1
#   slot m (d):        d * S_D  = d
#   slot m+1+i (k2_i): k2_i * S_K2
#   slot 2m+1 (Kannan): S_KANNAN
# ---------------------------------------------------------------------------
def verify_planted_vector(M, sigs, d_secret):
    m = len(sigs)
    dim = 2 * m + 2

    combos = []
    for s in sigs:
        unreduced = s['A'] + s['B'] * d_secret - LAM * s['k2']
        q = (unreduced - s['k1']) // N
        assert unreduced - s['k1'] == q * N, f"q_i not integer"
        combos.append(q)

    v = [0] * dim
    for i in range(m):
        for j in range(dim):
            v[j] -= combos[i] * M[i][j]
    for j in range(dim):
        v[j] += d_secret * M[m][j]
    for i in range(m):
        for j in range(dim):
            v[j] += sigs[i]['k2'] * M[m + 1 + i][j]
    for j in range(dim):
        v[j] += M[2 * m + 1][j]

    ok = True
    for i in range(m):
        expected = sigs[i]['k1'] * S_K1
        if v[i] != expected:
            print(f"  MISMATCH slot {i}: got {v[i]}, expected k1*S_K1={expected}")
            ok = False
    if v[m] != d_secret * S_D:
        print(f"  MISMATCH d-slot: got {v[m]}, expected d*S_D={d_secret*S_D}")
        ok = False
    for i in range(m):
        expected = sigs[i]['k2'] * S_K2
        if v[m + 1 + i] != expected:
            print(f"  MISMATCH k2 slot {m+1+i}: got {v[m+1+i]}, expected k2*S_K2={expected}")
            ok = False
    if v[dim - 1] != S_KANNAN:
        print(f"  MISMATCH Kannan: got {v[dim-1]}, expected S_KANNAN={S_KANNAN}")
        ok = False
    return ok, v

# ---------------------------------------------------------------------------
# Recovery: after LLL, scan reduced vectors for d.
# With S_D=1: v[m] == ±d (direct d recovery; no division needed).
# With S_KANNAN=N: v[dim-1] == ±N marks the Kannan embedding.
# ---------------------------------------------------------------------------
def verify_d_against_sigs(d_cand, sigs):
    """Check that d_cand is consistent with all signatures' GLV structure."""
    for s in sigs:
        kf = (s['A'] + s['B'] * d_cand) % N
        # kf must decompose as k1 + LAM*k2 with k1 in [0,K1_BOUND) and k2 in [0,K2_BOUND)
        found = False
        for k2 in range(K2_BOUND):
            k1 = (kf - LAM * k2) % N
            if 0 <= k1 < K1_BOUND:
                found = True
                break
        if not found:
            return False
    return True

def recover_d(reduced_rows, m, d_secret=None):
    """Scan LLL-reduced rows for the correct d.

    For each row where |last| == S_KANNAN, extract d = (±row[m]) mod N
    and verify it against the signatures.  Returns (d, witness_row) or (None, None).
    """
    dim = 2 * m + 2
    # All reduced rows sorted by norm (already sorted after LLL)
    for row in reduced_rows:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_candidate = (sign * row[m]) % N
        if d_candidate == 0:
            continue
        if d_secret is not None:
            if d_candidate == d_secret:
                return d_candidate, row
        # Without ground truth: verify against signatures (requires sigs in scope)
    return None, None

# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------
def run_experiment(m, d_secret, seed=42, verbose=True):
    sigs = gen_signatures(d_secret, m, seed)
    assert len(sigs) == m, f"Only generated {len(sigs)}/{m} signatures"

    if verbose:
        print(f"m={m} signatures (seed={seed}), d={d_secret}")
        for i, s in enumerate(sigs):
            print(f"  sig {i+1}: k1={s['k1']}, k2={s['k2']}, k_full={s['k_full']}")

    M = build_glv_lattice(sigs)
    dim = len(M)

    if verbose:
        print(f"\nLattice dim: {dim}x{dim}")
        print("Basis rows:")
        for i, row in enumerate(M):
            print(f"  [{i}]: {row}")

    # Verify planted vector
    ok, v = verify_planted_vector(M, sigs, d_secret)
    norm_v = sum(x*x for x in v) ** 0.5
    if verbose:
        print(f"\nPlanted combination verification: {'PASS' if ok else 'FAIL'}")
        print(f"Planted vector: {v}")
        print(f"Planted vector norm: {norm_v:.1f}")
        # Compare to basis vector norms
        row_norms = [sum(x*x for x in row)**0.5 for row in M]
        print(f"Basis row norms: {[f'{n:.1f}' for n in row_norms]}")
        print(f"Planted norm / min basis norm: {norm_v / min(row_norms):.2f}")

    # LLL reduction
    A = matrix_to_fpylll(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]

    if verbose:
        print("\nTop 5 LLL-reduced rows (sorted by norm):")
        norms = [(sum(x*x for x in r)**0.5, r) for r in reduced]
        norms.sort()
        for norm, row in norms[:5]:
            print(f"  norm={norm:.1f}  last={row[dim-1]}  d-slot={row[m]}  row={row}")

    d_recovered, witness = recover_d(reduced, m, d_secret)

    # Fallback: if direct k1-based recovery failed, try all rows with sig verification
    if d_recovered is None:
        for row in reduced:
            last = row[dim - 1]
            if abs(last) != S_KANNAN:
                continue
            sign = 1 if last > 0 else -1
            d_cand = (sign * row[m]) % N
            if d_cand == 0:
                continue
            if verify_d_against_sigs(d_cand, sigs):
                d_recovered = d_cand
                witness = row
                break

    if verbose:
        if d_recovered is not None and d_recovered == d_secret:
            print(f"\nRECOVERED d = {d_recovered}  (planted = {d_secret}) ✓")
            print(f"Witness row: {witness}")
        else:
            print(f"\nNOT RECOVERED (planted d = {d_secret}, best candidate = {d_recovered})")

    return d_recovered is not None and d_recovered == d_secret

# ---------------------------------------------------------------------------
# Sweep: test multiple seeds and m values
# ---------------------------------------------------------------------------
print("=" * 60)
print("GLV-aware HNP toy attack — Phase 2 lattice recovery")
print(f"Curve: y^2 = x^3 + {B_COEFF} over F_{P_FIELD}, n={N}, lambda={LAM}")
print(f"K1_BOUND={K1_BOUND} (k1 bias), K2_BOUND={K2_BOUND} (GLV domain)")
print("=" * 60)
print()

rng_main = random.Random(99)
d_test = rng_main.randint(1, N - 1)
print(f"Test secret d = {d_test}")
print()

# Information-theoretic threshold: need (K1*K2/n)^m << 1
# For K1=2, K2=15, n=199: threshold m where (2*15/199)^m ~ 1/199
# m >= log(199) / log(199/(2*15)) = log(199)/log(6.63) ~ 7.6/1.89 ~ 4.0 -> m=5
import math
eff = K1_BOUND * K2_BOUND / N
m_thresh = math.ceil(math.log(N) / math.log(1.0 / eff)) if eff < 1 else float('inf')
print(f"Info-theoretic threshold: (K1*K2/n)^m < 1/n  =>  m >= {m_thresh:.1f}")
print()

print("-" * 60)
print(f"Running with m={m_thresh+1:.0f} signatures (one above threshold):")
m_run = int(m_thresh) + 1
success = run_experiment(m_run, d_test, seed=42, verbose=True)
print()

print("-" * 60)
print("Sweep: multiple seeds, m values:")
results = {}
for m in [m_thresh - 1, m_thresh, m_thresh + 1, m_thresh + 2, m_thresh + 3, m_thresh + 4]:
    m = max(2, int(m))
    wins = 0
    for seed in [42, 1234, 9999, 7, 314159]:
        d_trial = random.Random(seed + 1000).randint(1, N - 1)
        w = run_experiment(m, d_trial, seed=seed, verbose=False)
        wins += w
    results[m] = wins
    print(f"  m={m}: {wins}/5 seeds recovered d")
print()
print("Summary table:")
for m, w in sorted(results.items()):
    print(f"  m={m}: {w}/5  {'(above threshold)' if m >= m_thresh else '(below threshold)'}")
