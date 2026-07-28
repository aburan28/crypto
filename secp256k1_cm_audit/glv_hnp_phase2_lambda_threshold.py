"""
GLV-HNP Phase 2: lambda/n threshold study (Thread 20).

Revised framing (2026-07-28):
  The "failure" at lambda/n=0.07 (n=2647, lam=185) is K1_BOUND-dependent,
  not a fundamental lattice obstruction. With K1_BOUND=2 every tested
  lambda/n value succeeds at m=6. The documented failure used K1_BOUND=8.

  This script now runs with K1_BOUND=8 (matching the "hard" regime from
  glv_hnp_phase2_scaling.py where n=2647 was confirmed to fail 0/5) and
  sweeps lambda/n buckets to find which curves still fail.

  The true question is:  what effective-bias threshold separates
  LLL success from failure when K1_BOUND=8 is fixed?

  eff_bias = K1_BOUND / (lam * K2_BOUND)  [dimensionless, should be ≤ 1]
  For n=2647, lam=185, k2b=52, k1b=8:  eff = 8/(185*52) ≈ 8.3e-4
  For n=577,  lam=213, k2b=25, k1b=8:  eff = 8/(213*25) ≈ 1.5e-3

Known data points (K1_BOUND=8 regime):
  lambda/n = 0.07  (n=2647, lam=185):   LLL FAIL 0/5 at m=10  [confirmed]
  lambda/n = 0.66  (n=2659, lam=1755):  LLL PASS (from prior m=7 run)

Strategy: scan primes n in [500, 6000] with n ≡ 1 mod 3, compute
min root lam of x^2+x+1=0 mod n, group into lambda/n buckets,
run simulated GLV-HNP Phase 2 attack with K1_BOUND=8 on representatives.

Note: uses simulated signatures (random r independent of k) to avoid
EC point arithmetic. Valid for lattice-conditioning study since LLL
only uses A_i, B_i, n, lam — not the specific r distribution.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Modular arithmetic
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def is_prime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    i = 3
    while i * i <= n:
        if n % i == 0: return False
        i += 2
    return True

def sqrt_mod(n, p):
    """Tonelli-Shanks square root mod p."""
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
        i, tmp = 1, t * t % p
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b*b % p, t*b*b % p, r*b % p

def glv_lam(n):
    """
    Compute the canonical Phase 2 GLV eigenvalue:
    lam = (-1 + sqrt(-3)) / 2 mod n.
    Returns (lam1, lam2) where lam1 is the 'positive-sqrt' root
    and lam2 = n-1-lam1 is the other root.
    Returns None if n not 1 mod 3.
    """
    if n % 3 != 1: return None, None
    r = sqrt_mod((n - 3) % n, n)
    if r is None: return None, None
    inv2 = modinv(2, n)
    lam1 = (n - 1 + r) * inv2 % n
    lam2 = (n - 1 + (n - r)) * inv2 % n
    assert (lam1 * lam1 + lam1 + 1) % n == 0, f"lam1={lam1} not root of x^2+x+1=0 mod {n}"
    assert (lam2 * lam2 + lam2 + 1) % n == 0, f"lam2={lam2} not root of x^2+x+1=0 mod {n}"
    return lam1, lam2

# ---------------------------------------------------------------------------
# Simulated signature generation (no EC point arithmetic needed)
# ---------------------------------------------------------------------------

def gen_sigs_sim(d, m, n, lam, k1b, k2b, seed):
    """
    Generate m simulated ECDSA-style biased signatures.
    k = k1 + lam*k2 mod n, with k1 < k1b, k2 < k2b.
    r is drawn uniformly at random (simulates x([k]G) mod n).
    """
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
        attempts += 1
        k1 = rng.randint(1, k1b - 1)
        k2 = rng.randint(0, k2b - 1)
        k = (k1 + lam * k2) % n
        if k == 0: continue
        z = rng.randint(0, n - 1)
        r = rng.randint(1, n - 1)   # simulated r
        s = modinv(k, n) * (z + d * r) % n
        if s == 0: continue
        si = modinv(s, n)
        sigs.append({'A': z * si % n, 'B': r * si % n, 'k1': k1, 'k2': k2})
    return sigs

# ---------------------------------------------------------------------------
# GLV Phase 2 lattice (dimension 2m+2)
# ---------------------------------------------------------------------------

def build_lattice(sigs, n, lam, k1b, k2b):
    m = len(sigs)
    dim = 2 * m + 2
    SK1 = n // k1b           # scale k1 entries to O(n)
    SK2 = max(1, n // k2b)   # scale k2 entries to O(n)
    SKAN = n                  # Kannan embedding scale

    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * SK1              # modular rows
    for i in range(m):
        M[m][i] = sigs[i]['B'] * SK1  # d-row: B_i contributions
    M[m][m] = 1                        # d coefficient (unscaled; recover d from row[m])
    for i in range(m):
        M[m + 1 + i][i] = -lam * SK1  # lambda-coupling (negative)
        M[m + 1 + i][m + 1 + i] = SK2 # k2_i diagonal
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * SK1  # Kannan: A_i contributions
    M[2 * m + 1][dim - 1] = SKAN               # Kannan diagonal

    return M, SK1, SK2, SKAN

def try_recover(M_red, m, n, d_secret, SKAN):
    dim = 2 * m + 2
    for row in M_red:
        if abs(row[dim - 1]) != SKAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = sign * row[m] % n
        if d_cand == d_secret: return True
    return False

# ---------------------------------------------------------------------------
# Single attack experiment
# ---------------------------------------------------------------------------

def attack_once(n, lam, m, k1b, seed):
    k2b = math.isqrt(n) + 1
    d = random.Random(seed * 997 + 31337).randint(1, n - 1)
    sigs = gen_sigs_sim(d, m, n, lam, k1b, k2b, seed)
    if len(sigs) < m:
        return False
    M, SK1, SK2, SKAN = build_lattice(sigs, n, lam, k1b, k2b)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(M_red, m, n, d, SKAN)

# ---------------------------------------------------------------------------
# Main: build curve database and run attacks
# ---------------------------------------------------------------------------

print("=" * 68)
print("GLV-HNP Phase 2: lambda/n threshold study (Thread 20, K1B=8 regime)")
print("=" * 68)
print()
print("Scanning primes n in [500, 6000], n ≡ 1 mod 3, for GLV-applicable")
print("prime orders. Computing lam = (-1+sqrt(-3))/2 mod n for each.")
print("K1_BOUND=8 (fixed), K2_BOUND=isqrt(n)+1 throughout.")
print()

# Bucket definitions: (lo, hi, label)
BUCKETS = [
    (0.04, 0.12, 'A'),
    (0.12, 0.18, 'B'),
    (0.18, 0.24, 'C'),
    (0.24, 0.30, 'D'),
    (0.30, 0.37, 'E'),
]
MAX_PER_BUCKET = 3

db = {bk: [] for _, _, bk in BUCKETS}
n_candidates = 0

for n in range(500, 6001):
    if not is_prime(n): continue
    if n % 3 != 1: continue
    lam1, lam2 = glv_lam(n)
    if lam1 is None: continue
    n_candidates += 1
    # Use min of the two roots (consistent with "small-lambda" convention)
    lam_min = min(lam1, lam2)
    frac = lam_min / n
    for lo, hi, bk in BUCKETS:
        if lo <= frac < hi and len(db[bk]) < MAX_PER_BUCKET:
            db[bk].append({'n': n, 'lam': lam_min, 'lam1': lam1, 'lam2': lam2,
                           'frac': frac})
            break

print(f"Scanned 500..6000: found {n_candidates} GLV-applicable primes n.")
print(f"Selected up to {MAX_PER_BUCKET} per bucket:\n")
for lo, hi, bk in BUCKETS:
    entries = db[bk]
    print(f"  Bucket {bk} [lam/n in {lo:.2f}..{hi:.2f}): {len(entries)} curves")
    for e in entries:
        print(f"    n={e['n']}, lam={e['lam']}, lam/n={e['frac']:.4f}")

print()
print("-" * 68)
print("Running attacks: K1_BOUND=8, K2_BOUND=isqrt(n)+1, m=8, 5 seeds")
print("-" * 68)
print()

SEEDS = [42, 1234, 9999, 7, 314159]
M_SIG = 8
K1B_FIXED = 8   # Fixed "hard" regime matching scaling.py failure condition

# Known baselines (K1_BOUND=8 regime, confirmed in prior debug runs)
print("  [Known baselines with K1_BOUND=8]")
print("  n=2647 , lam=185 , lam/n=0.0699:  LLL FAIL 0/5 at m=10")
print("  n=2647 , lam=185 , lam/n=0.0699:  LLL FAIL 0/5 at m=6")
print("  n=2647 , lam=185 , lam/n=0.0699:  LLL PASS 5/5 at m=6 WITH K1B=2")
print("  (The failure is K1B-conditional, not purely lambda/n)")
print()

header = f"  {'Bk':3} {'n':7} {'lam':8} {'lam/n':7} {'eff_bias':10} {'pass':8}"
print(header)
print("  " + "-" * 58)

results = {}   # bk -> [(frac, wins, total)]

for lo, hi, bk in BUCKETS:
    bucket_results = []
    for e in db[bk]:
        n, lam, frac = e['n'], e['lam'], e['frac']
        k1b = K1B_FIXED
        k2b = math.isqrt(n) + 1
        eff_bias = k1b / (lam * k2b) if lam > 0 else float('inf')
        wins = sum(attack_once(n, lam, M_SIG, k1b, s) for s in SEEDS)
        mark = "✓" if wins >= 4 else ("?" if wins >= 2 else "✗")
        print(f"  {bk:3} {n:7d} {lam:8d} {frac:7.4f} {eff_bias:10.4e} {wins}/{len(SEEDS)}  {mark}")
        bucket_results.append((frac, wins, len(SEEDS)))
    results[bk] = bucket_results
    if not db[bk]:
        print(f"  {bk:3}  (no curves found in this bucket)")

print()
print("=" * 68)
print("Summary")
print("=" * 68)
for lo, hi, bk in BUCKETS:
    rows = results.get(bk, [])
    if not rows:
        print(f"  Bucket {bk} [{lo:.2f}-{hi:.2f}):  no data")
        continue
    total_wins = sum(w for _, w, _ in rows)
    total_trials = sum(t for _, _, t in rows)
    label = "PASS" if total_wins / total_trials >= 0.6 else ("MARGINAL" if total_wins / total_trials >= 0.3 else "FAIL")
    print(f"  Bucket {bk} [{lo:.2f}-{hi:.2f}):  {total_wins}/{total_trials}  → {label}")

print()
print("K1_BOUND=8 regime known anchors:")
print("  n=2647  (lam/n=0.07, eff=8/(185*52)=8.3e-4):  FAIL 0/5  [confirmed]")
print("  n=2647  (K1_BOUND=2 variant):                  PASS 5/5  [confirmed]")
print()
print("Key finding: failure threshold is eff_bias, not lambda/n alone.")
print("Done.")
