"""
GLV-HNP Phase 2: λ/n threshold bisection (Thread 20).

Objective: Pinpoint the lower threshold below which the Phase 2 lattice attack
fails on j=0 GLV-aware curves.  From 2026-07-26:
  - λ/n = 0.07 (p=2677, 12-bit):  LLL + BKZ(40) FAIL
  - λ/n = 0.34 (p=524347, 20-bit): LLL succeeds at m=9 (3/3)
  - λ/n ≈ 0.33 (secp256k1):       OBSTRUCTED (Effect A, near-1/3 resonance)

This script sweeps λ/n ∈ {0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28}
on 20-bit j=0 CM curves to locate the Effect B failure boundary.

Effect A (resonance near 1/3): λ/n ≈ 0.333 — separate obstruction, not studied here.
Effect B (small λ): λ/n small → -λ·S_K1 entry in λ-rows is small relative to
    n·S_K1 diagonal; λ-rows become near-zero, leading to ill-conditioned lattice.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (j=0, a=0)
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

def find_generator(p, b_param, n):
    rng = random.Random(12345)
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
# CM / GLV theory
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
    for b_try in range(1, min(p, 500)):
        rng2 = random.Random(99)
        for _ in range(20):
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
# Curve discovery: scan 20-bit j=0 CM primes, bucket by λ/n target
# ---------------------------------------------------------------------------

LAM_TARGETS = [0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28]
LAM_TOL = 0.015      # ±1.5% of n
MAX_PER_TARGET = 2   # 2 curves per target for reliability

def find_curves_by_lamn(p_lo=2**18, p_hi=2**21, targets=LAM_TARGETS,
                         tol=LAM_TOL, max_per=MAX_PER_TARGET):
    found = {t: [] for t in targets}
    p = sympy.nextprime(p_lo - 1)
    scanned = 0
    while p < p_hi:
        scanned += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for tr in j0_traces(a_e, b_e):
                    n_cand = p + 1 - tr
                    if n_cand < 2:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    lam_n = lam / n_cand
                    for t in targets:
                        if abs(lam_n - t) <= tol and len(found[t]) < max_per:
                            b_param = find_b_for_n(p, n_cand)
                            if b_param is None:
                                continue
                            found[t].append((p, b_param, n_cand, lam, lam_n))
                            print(f"  [lam/n≈{t:.2f}] p={p}, n={n_cand}, "
                                  f"lam={lam}, lam/n={lam_n:.5f}")
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    print(f"  (scanned {scanned} primes in [{p_lo}, {p_hi}])")
    return found

# ---------------------------------------------------------------------------
# LLL attack: Phase 2 lattice
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m, K1, K2, seed):
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
        attempts += 1
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(0, K2 - 1)
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
        sigs.append({'A': A, 'B': B})
    return d_secret, sigs

def build_lattice(sigs, n, lam, K1, K2):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // K1
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
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_lll(p, b_param, n, lam, G, m, K1, K2, seed):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, K2, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, K1, K2)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep(p, b_param, n, lam, G, K1, K2, m_range, seeds):
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_lll(p, b_param, n, lam, G, m, K1, K2, s))
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

K1 = 72
SEEDS = [42, 1234, 9999]
M_RANGE = range(10, 23)   # m = 10..22

print("=" * 72)
print("Thread 20: λ/n threshold bisection — Effect B lower boundary")
print(f"K1={K1}, m=10..22, seeds={SEEDS}")
print(f"Targets: {[f'{t:.2f}' for t in LAM_TARGETS]}")
print("=" * 72)
print()

# ---- Scan for curves -------------------------------------------------------
print("Scanning 20-bit j=0 CM primes for λ/n targets...")
bucket = find_curves_by_lamn()
print()

# ---- LLL sweep per curve ---------------------------------------------------
print("=" * 72)
print("LLL sweep results")
print("-" * 72)

all_results = []  # (lam_n_actual, p, n, lam, first_win_or_none, outcome_str)

for t in LAM_TARGETS:
    curves = bucket[t]
    print(f"\n[λ/n≈{t:.2f}]  ({len(curves)} curve(s) found)")
    if not curves:
        print("  NO CURVES FOUND — skipping target")
        all_results.append((t, None, None, None, None, "NO CURVE"))
        continue
    for (p, b_param, n, lam, lam_n) in curves:
        K2 = math.isqrt(n) + 1
        G = find_generator(p, b_param, n)
        if G is None:
            print(f"  p={p}, n={n}: generator not found, skip")
            continue
        res = sweep(p, b_param, n, lam, G, K1, K2, M_RANGE, SEEDS)
        first_win = next((m for m, (w, _) in sorted(res.items()) if w == len(SEEDS)), None)
        detail = "  ".join(f"m{m}:{w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if first_win:
            outcome = f"SUCCESS m={first_win}"
        else:
            max_w = max(w for w, _ in res.values())
            outcome = f"OBSTRUCTED (max {max_w}/3)"
        print(f"  p={p}, n={n}, lam={lam}, lam/n={lam_n:.5f}")
        print(f"  {detail}")
        print(f"  → {outcome}")
        all_results.append((lam_n, p, n, lam, first_win, outcome))

# ---- Summary table --------------------------------------------------------
print()
print("=" * 72)
print("SUMMARY TABLE: λ/n vs LLL Phase 2 attack outcome")
print(f"{'λ/n':<9} {'p':<10} {'n':<10} {'first_3/3':<12} {'outcome'}")
print("-" * 72)
for (lam_n, p, n, lam, fw, outcome) in sorted(all_results):
    if p is None:
        print(f"  {lam_n:.5f}  (no curve found)                          {outcome}")
    else:
        fw_str = str(fw) if fw else "never"
        print(f"  {lam_n:.5f}  {p:<10} {n:<10} {fw_str:<12} {outcome}")

# Known reference points
print()
print("Reference points (from prior runs):")
print(f"  0.07000   p=2677,   n=2647,   never        OBSTRUCTED [2026-07-26 log]")
print(f"  0.34000   p=524347, n=523969, m=9           SUCCESS    [2026-07-26 log]")
print(f"  0.32600   (secp256k1 analogue)               OBSTRUCTED [Effect A, ≈1/3]")
print()

# Threshold summary
successes = [(ln, fw) for (ln, _, _, _, fw, _) in all_results if fw is not None]
failures  = [(ln,) for (ln, _, _, _, fw, out) in all_results
             if fw is None and out not in ("NO CURVE",)]
if successes and failures:
    min_success = min(successes, key=lambda x: x[0])
    max_failure  = max(failures, key=lambda x: x[0])
    print(f"Lowest success : λ/n = {min_success[0]:.5f}  (first 3/3 at m={min_success[1]})")
    print(f"Highest failure: λ/n = {max_failure[0]:.5f}")
    print(f"Threshold bracket: [{max_failure[0]:.5f}, {min_success[0]:.5f}]")

print()
print("Done.")
