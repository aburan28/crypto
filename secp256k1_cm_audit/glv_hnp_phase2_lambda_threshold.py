"""
GLV-HNP Phase 2 — λ/n threshold sweep.

Bisects the failure threshold between λ/n = 0.07 (known FAIL) and
λ/n = 0.34 (known PASS) for the Phase 2 lattice attack.

Strategy:
- Find toy prime-order j=0 curves covering λ/n in bins [0.05, 0.50).
- Fix K1_BOUND = 2 (k1 ∈ {0,1}) so the IT threshold is always met with m≥4.
- Run LLL (m=6 sigs, 3 seeds) for each bin and record win rate.
- Report the empirical λ/n threshold.

Known data points from 2026-07-26 run:
  λ/n ≈ 0.44  (p=211, n=199):   success 5/5
  λ/n ≈ 0.07  (p=2677, n=2647): failure 0/3 (also BKZ failed)

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Arithmetic helpers
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def tonelli_shanks(n, p):
    """Return r with r^2 ≡ n mod p (smallest), or None."""
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
        i, tmp = 1, t * t % p
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mc - i - 1), p)
        mc, c, t, r = i, b * b % p, t * c * c % p, r * b % p

def is_prime_mr(n):
    """Miller-Rabin (deterministic for n < 3.3e24)."""
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

# ---------------------------------------------------------------------------
# Elliptic curve arithmetic (short Weierstrass, a=0)
# ---------------------------------------------------------------------------

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

def curve_order_brute(p, b):
    """#E(F_p) for y^2 = x^3 + b by Legendre sum (O(p))."""
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x * x * x + b) % p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2
    return count

def find_glv_eigenvalue(n):
    """Find smallest root of x^2+x+1=0 mod n by brute force (n < ~10000)."""
    for lam in range(1, n):
        if (lam * lam + lam + 1) % n == 0:
            return lam
    return None

def find_generator(p, b, n):
    """Generator of order n on y^2 = x^3 + b / F_p."""
    for x in range(p):
        rhs = (x * x * x + b) % p
        sq = tonelli_shanks(rhs, p)
        if sq is None: continue
        for y in [sq, p - sq]:
            if ec_mul((x, y), n, p) is None:
                return (x, y)
    return None

# ---------------------------------------------------------------------------
# Curve catalogue: search for prime-order j=0 curves, grouped by λ/n bin
# ---------------------------------------------------------------------------

BIN_WIDTH = 0.03   # λ/n bins of width 0.03
BIN_COUNT = 16     # bins covering [0.03, 0.51)

def bin_idx(ratio):
    idx = int(ratio / BIN_WIDTH)
    return min(idx, BIN_COUNT - 1)

def find_curves(p_max=6000, b_vals=(1, 2, 3, 5, 6, 7, 10, 11)):
    """Return best representative curve per λ/n bin."""
    bins = {}
    for p in range(7, p_max + 1):
        if not is_prime_mr(p): continue
        if p % 6 != 1: continue
        for b in b_vals:
            if b % p == 0: continue
            n = curve_order_brute(p, b)
            if not is_prime_mr(n): continue
            if n % 3 != 1: continue  # GLV requires 3 | (n-1)
            lam = find_glv_eigenvalue(n)
            if lam is None: continue
            lam_small = min(lam, n - 1 - lam)
            ratio = lam_small / n
            bi = bin_idx(ratio)
            if bi not in bins:
                bins[bi] = (p, b, n, lam_small, ratio)
        if len(bins) >= BIN_COUNT:
            # Quick-exit: check if all interesting bins filled
            covered = sum(1 for i in range(1, 15) if i in bins)
            if covered >= 13:
                break
    return bins

# ---------------------------------------------------------------------------
# Phase 2 GLV lattice attack
# ---------------------------------------------------------------------------

K1_BOUND = 2       # k1 ∈ {0, 1}: very tight bias, IT threshold easily met
M_SIGS   = 6       # signatures per trial
SEEDS    = [42, 1234, 9999]

def gen_sigs(G, d, n, lam, p, k1_bound, k2_bound, m, rng):
    """Generate m signatures with k1 ∈ [0,k1_bound), k2 ∈ [0,k2_bound)."""
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200_000:
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
        s = modinv(k_full, n) * (h + d * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs if len(sigs) == m else None

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Build column-scaled (2m+2)×(2m+2) GLV Phase 2 lattice."""
    m = len(sigs)
    dim = 2 * m + 2

    S_K1 = n // k1_bound
    S_D  = 1
    S_K2 = n // k2_bound
    S_KAN = n

    M = [[0] * dim for _ in range(dim)]

    # Rows 0..m-1: mod-n on k1 slots
    for i in range(m):
        M[i][i] = n * S_K1

    # Row m: d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D

    # Rows m+1..2m: k2 rows
    for i in range(m):
        M[m + 1 + i][i] = (-lam * S_K1) % (n * S_K1)
        # Store the negative explicitly (LLL handles negative integers)
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2

    # Row 2m+1: Kannan embedding
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KAN

    return M, S_K1, S_D, S_K2, S_KAN

def recover_d(reduced, dim, m, n, lam, k2_bound, S_D, S_KAN):
    """Scan reduced rows for d (Kannan marker row[dim-1] == ±S_KAN)."""
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m] * modinv(S_D, n)) % n
        if d_cand == 0: continue
        yield d_cand

def verify_d(d_cand, sigs, n, lam, k1_bound, k2_bound):
    """Check d_cand is consistent with all sigs' GLV structure."""
    for s in sigs:
        kf = (s['A'] + s['B'] * d_cand) % n
        ok = False
        for k2 in range(k2_bound):
            k1 = (kf - lam * k2) % n
            if 0 <= k1 < k1_bound:
                ok = True
                break
        if not ok: return False
    return True

def run_attack(p, b, n, lam, seed, verbose=False):
    """Run Phase 2 attack for one seed. Returns True if d recovered."""
    k2_bound = math.isqrt(n) + 1
    rng = random.Random(seed)
    G = find_generator(p, b, n)
    if G is None: return False

    d_secret = rng.randint(1, n - 1)
    sigs = gen_sigs(G, d_secret, n, lam, p, K1_BOUND, k2_bound, M_SIGS, rng)
    if sigs is None:
        if verbose: print(f"  seed={seed}: could not generate {M_SIGS} sigs")
        return False

    M, S_K1, S_D, S_K2, S_KAN = build_lattice(sigs, n, lam, K1_BOUND, k2_bound)
    dim = len(M)

    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]

    for d_cand in recover_d(reduced, dim, M_SIGS, n, lam, k2_bound, S_D, S_KAN):
        if d_cand == d_secret:
            if verbose: print(f"  seed={seed}: RECOVERED d={d_secret}")
            return True
        if verify_d(d_cand, sigs, n, lam, K1_BOUND, k2_bound):
            if d_cand == d_secret:
                if verbose: print(f"  seed={seed}: RECOVERED d={d_secret}")
                return True
    if verbose: print(f"  seed={seed}: FAIL (planted d={d_secret})")
    return False

# ---------------------------------------------------------------------------
# Main: sweep
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2 — λ/n threshold sweep")
print(f"K1_BOUND={K1_BOUND}, m_sigs={M_SIGS}, seeds={SEEDS}")
print(f"Bins: width={BIN_WIDTH}, covering [0.03, {BIN_WIDTH*BIN_COUNT:.2f})")
print("=" * 70)
print()

# Include the two known data points as anchors
KNOWN = [
    # (p, b, n, lam, tag) — validated in prior runs
    (211, 2, 199, 106, "known-pass λ/n=0.44"),
    (2677, 7, 2647, 185, "known-fail λ/n=0.07"),
]

print("Step 1: Searching for prime-order j=0 curves by λ/n bin...")
bins = find_curves(p_max=6000)
print(f"Found representatives for {len(bins)} bins.\n")

# Show catalogue
print("Bin catalogue (λ/n → curve):")
for bi in sorted(bins):
    p, b, n, lam, ratio = bins[bi]
    print(f"  bin {bi:2d}  λ/n=[{bi*BIN_WIDTH:.2f},{(bi+1)*BIN_WIDTH:.2f})  "
          f"p={p}, b={b}, n={n}, λ={lam}  ratio={ratio:.4f}")
print()

# Add anchors to test set
test_cases = []
for bi in sorted(bins):
    p, b, n, lam, ratio = bins[bi]
    test_cases.append((p, b, n, lam, ratio, f"bin{bi}"))
# Append known anchors
test_cases.append((211, 2, 199, 106, 106/199, "known-pass"))
test_cases.append((2677, 7, 2647, 185, 185/2647, "known-fail"))

print("Step 2: Running Phase 2 attack (LLL, m=6, 3 seeds) per curve...")
print()

results = []
for p, b, n, lam, ratio, tag in test_cases:
    wins = 0
    for seed in SEEDS:
        if run_attack(p, b, n, lam, seed, verbose=False):
            wins += 1
    results.append((ratio, lam, n, wins, tag))
    bar = "✓" * wins + "✗" * (len(SEEDS) - wins)
    print(f"  λ/n={ratio:.4f}  n={n}  λ={lam}  wins={wins}/{len(SEEDS)}  {bar}  [{tag}]")

print()
print("=" * 70)
print("Summary: λ/n vs. LLL success (3-seed majority = PASS if wins ≥ 2)")
print("=" * 70)

results_sorted = sorted(results, key=lambda x: x[0])
for ratio, lam, n, wins, tag in results_sorted:
    verdict = "PASS" if wins >= 2 else "FAIL"
    print(f"  λ/n={ratio:.4f}  n={n}  {wins}/{len(SEEDS)}  {verdict}  [{tag}]")

# Identify threshold
passes = [(r, w) for r, _, _, w, _ in results_sorted if w >= 2]
fails  = [(r, w) for r, _, _, w, _ in results_sorted if w < 2]

if passes and fails:
    min_pass = min(r for r, _ in passes)
    max_fail = max(r for r, _ in fails)
    print()
    if max_fail < min_pass:
        print(f"Threshold: λ/n in ({max_fail:.4f}, {min_pass:.4f})")
        print(f"  All λ/n ≤ {max_fail:.4f} → FAIL")
        print(f"  All λ/n ≥ {min_pass:.4f} → PASS")
    else:
        print(f"Non-monotone: max_fail={max_fail:.4f}, min_pass={min_pass:.4f}")
        print("Threshold is not clean — need more granularity or larger m.")
elif not fails:
    print("\nAll bins PASS — threshold below 0.03 or attack is robust everywhere.")
elif not passes:
    print("\nAll bins FAIL — threshold above 0.50 or attack broken.")

print()
print("Done.")
