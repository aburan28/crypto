"""
GLV-aware HNP attack — lambda/n threshold study (Thread 20).

Question: what is the minimum lambda/n ratio such that the column-scaled
(2m+2)×(2m+2) lattice admits LLL recovery of d?

Method: SYNTHETIC curves — no EC arithmetic needed.
  - Fix n (prime), K1_BOUND=2, K2_BOUND=ceil(sqrt(n)).
  - Choose lambda = round(ratio * n) for ratio in {0.05, 0.10, ..., 0.50}.
  - Generate artificial signatures: k_i = k1_i + lambda*k2_i mod n,
    choose random B_i ≠ 0, A_i = (k_i - B_i*d) mod n.
  - Build and LLL-reduce the (2m+2)×(2m+2) lattice (same column scaling as Phase 2).
  - Count successes over 5 seeds for each (ratio, m) pair.

Also tests: BKZ(β=20, 40) rescue for borderline ratios.

Expected output: a threshold ratio ~[0.07, 0.34] at which 3/5 success first appears.

Reference: autolab 2026-07-26, Thread 5 findings.
"""

import random
import math
from fpylll import IntegerMatrix, LLL, BKZ


# ---------------------------------------------------------------------------
# Primality helper
# ---------------------------------------------------------------------------
def isprime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for p in range(3, min(n, int(n**0.5) + 2), 2):
        if n % p == 0: return False
    return True

def nextprime(n):
    n = max(n, 2)
    while not isprime(n):
        n += 1
    return n


# ---------------------------------------------------------------------------
# Synthetic signature generation (no real curve needed)
# ---------------------------------------------------------------------------
def gen_synthetic_sigs(d, m, n, lam, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 10 * m:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0:
            continue
        # Artificial ECDSA-like: choose random B ≠ 0, set A = (k - B*d) mod n
        B = rng.randint(1, n - 1)
        A = (k_full - B * d) % n
        assert (A + B * d) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs


# ---------------------------------------------------------------------------
# Lattice construction (same as Phase 2 column-scaled design)
# ---------------------------------------------------------------------------
def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2

    S_K1     = n // k1_bound
    S_D      = 1
    S_K2     = n // k2_bound
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]

    # Rows 0..m-1: mod-n constraints
    for i in range(m):
        M[i][i] = n * S_K1

    # Row m: d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D

    # Rows m+1..2m: k2_i rows
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2

    # Row 2m+1: Kannan target row
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN

    return M, S_K1, S_D, S_K2, S_KANNAN


def recover_d(reduced, m, n, d_true):
    """Scan LLL-reduced rows for correct d via Kannan slot = ±S_KANNAN."""
    dim = 2 * m + 2
    S_KANNAN = n
    for row in reduced:
        if abs(row[dim - 1]) != S_KANNAN:
            continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_true:
            return True
    return False


def run_attack(m, n, lam, k1_bound, k2_bound, d, seed, reduction='lll', bkz_beta=20):
    sigs = gen_synthetic_sigs(d, m, n, lam, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False

    M, _, _, _, _ = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)

    if reduction == 'lll':
        LLL.reduction(A)
    elif reduction == 'bkz':
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)

    dim = 2 * m + 2
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, d)


# ---------------------------------------------------------------------------
# Compute planted vector norm and smallest spurious norm
# ---------------------------------------------------------------------------
def planted_info(n, lam, k1_bound, k2_bound, m):
    S_K1 = n // k1_bound
    S_K2 = n // k2_bound
    S_KANNAN = n
    # Worst-case planted norm (k1 ≈ k1_bound, k2 ≈ k2_bound, d ≈ n)
    planted_sq = m * ((k1_bound - 1) * S_K1)**2 + n**2 + m * ((k2_bound - 1) * S_K2)**2 + S_KANNAN**2
    # Simplest spurious: 2*(k2 row_i) + (mod-n row_i)  →  entries: (2*(-lam)+n)*S_K1 and 2*S_K2
    # For each i independently. Combined norm (1 pair):
    lam_effective = (n - lam) % n if lam > n // 2 else lam   # min(lam, n-lam)
    spurious_sq = (2 * lam_effective * S_K1)**2 + (2 * S_K2)**2
    return math.sqrt(planted_sq), math.sqrt(spurious_sq), lam_effective / n


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
N_BITS = 12
N = nextprime(2**N_BITS)
K2_BOUND = math.isqrt(N) + 1
# Use K1_BOUND=8 to match the 2026-07-26 failure condition (n=2647, K1_BOUND=8).
# K1_BOUND=2 trivially passes for all lambda/n (planted vector is n/2 times shorter
# than any basis vector — too easy). With K1_BOUND=8 the competition is realistic.
K1_BOUND = 8
SEEDS = [42, 1234, 9999, 7, 314159]
N_SEEDS = len(SEEDS)
SUCCESS_THRESHOLD = 3   # 3/5 seeds = "success"

print("=" * 70)
print("GLV-aware HNP — lambda/n threshold study")
print(f"n={N} ({N_BITS}-bit prime), K1_BOUND={K1_BOUND}, K2_BOUND={K2_BOUND}")
print(f"eff = K1*K2/n = {K1_BOUND*K2_BOUND/N:.4f}")
print(f"Seeds: {SEEDS}, success_threshold={SUCCESS_THRESHOLD}/{N_SEEDS}")
print("=" * 70)

# Ratios to test (SYMMETRIC: ratio and 1-ratio should both be tested)
RATIOS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
M_VALUES = [5, 7, 10, 15, 20]

D_FIXED = random.Random(0xABCD).randint(1, N - 1)
print(f"Fixed d = {D_FIXED}")
print()

# For each ratio, find the minimum m achieving SUCCESS_THRESHOLD/N_SEEDS
threshold_by_ratio = {}

print(f"{'ratio':>6}  {'lambda':>6}  {'lam/n':>6}  {'min(l,n-l)/n':>13} | ", end="")
for m in M_VALUES:
    print(f" m={m:2d} ", end="")
print()
print("-" * 70)

for ratio in RATIOS:
    lam = max(1, min(N - 1, round(ratio * N)))
    actual_ratio = lam / N
    planted_norm, spurious_norm, sym_ratio = planted_info(N, lam, K1_BOUND, K2_BOUND, M_VALUES[0])

    win_counts = []
    for m in M_VALUES:
        wins = sum(
            run_attack(m, N, lam, K1_BOUND, K2_BOUND, D_FIXED, seed, 'lll')
            for seed in SEEDS
        )
        win_counts.append(wins)

    # Find minimum m with >= SUCCESS_THRESHOLD wins
    min_m = None
    for m, w in zip(M_VALUES, win_counts):
        if w >= SUCCESS_THRESHOLD:
            min_m = m
            break

    threshold_by_ratio[ratio] = (lam, actual_ratio, sym_ratio, min_m, win_counts)

    print(f"{ratio:>6.2f}  {lam:>6d}  {actual_ratio:>6.3f}  {sym_ratio:>13.3f} | ", end="")
    for w in win_counts:
        marker = "*" if w >= SUCCESS_THRESHOLD else " "
        print(f" {w}/{N_SEEDS}{marker}", end="")
    print(f"  → min_m={'N/A' if min_m is None else min_m}")

print()

# Summary
print("=" * 70)
print("SUMMARY — LLL threshold")
print()
success_ratios = [(r, d[0], d[2], d[3]) for r, d in threshold_by_ratio.items() if d[3] is not None]
fail_ratios    = [(r, d[0], d[2])       for r, d in threshold_by_ratio.items() if d[3] is None]

if success_ratios:
    min_success = min(success_ratios, key=lambda x: x[2])
    print(f"Smallest PASSING ratio: {min_success[0]:.2f}  (lam={min_success[1]}, min(l,n-l)/n={min_success[2]:.3f}, min_m={min_success[3]})")
if fail_ratios:
    max_fail = max(fail_ratios, key=lambda x: x[2])
    print(f"Largest FAILING ratio:  {max_fail[0]:.2f}  (lam={max_fail[1]}, min(l,n-l)/n={max_fail[2]:.3f})")
print()

# --------------------------------------------------------------------------
# BKZ rescue: test borderline ratios (just below the LLL threshold) with BKZ
# --------------------------------------------------------------------------
print("=" * 70)
print("BKZ RESCUE — testing fail-ratios with BKZ(beta=20) and BKZ(beta=40)")
print()

# Find up to 3 ratios that LLL fails on
borderline = [r for r, d in threshold_by_ratio.items() if d[3] is None][:3]

for ratio in borderline:
    lam = threshold_by_ratio[ratio][0]
    actual_ratio = threshold_by_ratio[ratio][1]
    sym_ratio = threshold_by_ratio[ratio][2]
    print(f"ratio={ratio:.2f}  lam={lam}  min(l,n-l)/n={sym_ratio:.3f}")

    for beta in [20, 40]:
        wins_m = {}
        for m in [10, 15, 20]:
            wins = sum(
                run_attack(m, N, lam, K1_BOUND, K2_BOUND, D_FIXED, seed, 'bkz', beta)
                for seed in SEEDS
            )
            wins_m[m] = wins
        win_str = ", ".join(f"m={m}:{w}/{N_SEEDS}" for m, w in wins_m.items())
        status = "*PASS*" if any(w >= SUCCESS_THRESHOLD for w in wins_m.values()) else "FAIL"
        print(f"  BKZ(beta={beta:2d}): {win_str}  [{status}]")
    print()

# --------------------------------------------------------------------------
# Geometric analysis: compare effective GS length vs planted norm
# --------------------------------------------------------------------------
print("=" * 70)
print("GEOMETRIC ANALYSIS — planted norm vs spurious-1 norm vs GS-k2 component")
print()
print(f"{'ratio':>6}  {'lam':>6}  {'planted_n':>10}  {'spurious_n':>10}  {'gs_k2_len':>10}  {'ratio_sq':>9}")
print("-" * 65)

M_ANAL = 10  # fixed m for analysis
for ratio in RATIOS:
    lam = max(1, min(N - 1, round(ratio * N)))
    S_K1 = N // K1_BOUND
    S_K2 = N // K2_BOUND
    S_KANNAN = N

    # Planted norm (approximate, using expected values)
    planted_sq = M_ANAL * ((K1_BOUND - 1) * S_K1)**2 + N**2 + M_ANAL * ((K2_BOUND - 1) * S_K2)**2 + S_KANNAN**2
    planted_n = math.sqrt(planted_sq)

    # Smallest spurious: sum over 1 pair (2*(k2-row) + mod-n-row)
    lam_eff = min(lam, N - lam)
    spurious_sq = ((2 * lam_eff) * S_K1)**2 + (2 * S_K2)**2
    spurious_n = math.sqrt(spurious_sq)

    # GS k2-row component: the -lam*S_K1 entry reduced mod n*S_K1 → GS length = lam_eff * S_K1
    # (after projecting out mod-n rows)
    gs_k2 = lam_eff * S_K1

    print(f"{ratio:>6.2f}  {lam:>6d}  {planted_n:>10.1f}  {spurious_n:>10.1f}  {gs_k2:>10.1f}  {(gs_k2/S_K1)**2/N:>9.4f}")

print()
print("gs_k2_len = min(lam, n-lam) * S_K1 = effective GS component of k2 row in k1-column basis")
print("ratio_sq  = (gs_k2_len / (n*S_K1))^2 = fractional GS-squared contribution")
