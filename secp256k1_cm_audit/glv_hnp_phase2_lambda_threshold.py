"""
GLV-HNP Phase 2: λ/n lower-bound threshold experiment (Thread 20, 2026-07-27).

Context:
  From 2026-06-21 (glv_hnp_delta_threshold.py), using 20-bit j=0 CM curves,
  K1=72, m=10..20, seeds=3:
    λ/n ≈ 0.016 → SUCCESS at m=10   (Effect B: small-λ near-orthogonal rows)
    λ/n ≈ 0.039 → SUCCESS at m=20
    λ/n ≈ 0.281-0.387 → OBSTRUCTED (Effect A: near-linear dependence at λ/n≈1/3)

  From 2026-06-15/16 (glv_hnp_phase2_20bit.py), using 12-bit curve p=2677:
    λ/n = 0.07 → FAIL at m≤12, LLL and BKZ(40) both fail.
    Root cause: spurious Kannan vectors shorter than planted vector.

Gap to fill: λ/n ∈ (0.04, 0.20) at 20-bit scale.

Goal: find the lower threshold τ_lo such that:
  λ/n < τ_lo  → SUCCESS (Effect B, small λ)
  τ_lo ≤ λ/n  → FAIL (spurious Kannan vector dominates)
  λ/n ≥ ~0.20 → SUCCESS again (Effect B recovers, or Effect A hasn't kicked in)

Method: scan 20-bit j=0 CM primes (Eisenstein decomp), find 2 curves per
λ/n target ∈ {0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21},
run K1=72 LLL at m=10..25 (3 seeds), and report first 3/3 or FAIL.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass a=0, y^2 = x^3 + b)
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

def find_generator(p, b, n):
    rng = random.Random(12345)
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
# CM theory: Eisenstein decomposition + j=0 traces
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

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

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
# Curve search by λ/n target
# ---------------------------------------------------------------------------

LAM_TARGETS = [0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21]
LAM_TOL = 0.012        # ±1.2% of n
MAX_PER_TARGET = 2

def find_curves_by_lamn(p_lo=2**19, p_hi=2**21, targets=LAM_TARGETS,
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
                            print(f"  [lam/n≈{t:.2f}] p={p}, n={n_cand}, lam={lam}, "
                                  f"lam/n={lam_n:.5f}")
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    print(f"  (scanned {scanned} primes up to {p})")
    return found

# ---------------------------------------------------------------------------
# Anchor curves from June 2026 experiments (verified)
# ---------------------------------------------------------------------------

# From 2026-06-21: λ/n=0.016 SUCCESS at m=10, λ/n=0.039 SUCCESS at m=20
ANCHORS = [
    dict(p=524983, n=526429, lam=8367,   lam_n=0.0159, label="anchor-016"),
    dict(p=525127, n=524149, lam=20387,  lam_n=0.0389, label="anchor-039"),
    # Control: known FAIL at 12-bit (lam/n=0.07); tested here at 20-bit for comparison
    # (no 20-bit equivalent found in prior sessions)
]

# ---------------------------------------------------------------------------
# GLV-HNP lattice attack
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
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
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})
    return d_secret, sigs

def build_lattice(sigs, n, lam, K1):
    m = len(sigs)
    K2 = math.isqrt(n) + 1
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

def run_lll(p, b_param, n, lam, G, m, K1, seed):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, K1)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep_curve(p, b_param, n, lam, G, K1, m_range, seeds):
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_lll(p, b_param, n, lam, G, m, K1, s))
        results[m] = (wins, len(seeds))
    first_win = next((m for m, (w, t) in sorted(results.items()) if w == t), None)
    return results, first_win

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

K1 = 72
SEEDS = [42, 1234, 9999]
M_RANGE = range(10, 26)   # m=10..25 (same as lamn_boundary.py)

print("=" * 72)
print("GLV-HNP Phase 2: λ/n lower-bound threshold (Thread 20, 2026-07-27)")
print(f"K1={K1}, m={M_RANGE.start}..{M_RANGE.stop-1}, seeds={SEEDS}")
print("=" * 72)

# Step 1: Find curves at target λ/n values
print("\n[Step 1] Scanning 20-bit j=0 CM primes for λ/n targets...")
bucket = find_curves_by_lamn()
print()

# Step 2: Check anchor curves
print("[Step 2] Verifying anchor curves from 2026-06-21...")
anchor_curves = []
for anc in ANCHORS:
    p, n, lam = anc['p'], anc['n'], anc['lam']
    b = find_b_for_n(p, n)
    if b is None:
        print(f"  {anc['label']}: b not found, skip")
        continue
    G = find_generator(p, b, n)
    if G is None:
        print(f"  {anc['label']}: generator not found, skip")
        continue
    anchor_curves.append((p, b, n, lam, anc['lam_n'], G, anc['label']))
    print(f"  {anc['label']}: p={p}, n={n}, lam={lam}, lam/n={anc['lam_n']:.4f} OK")
print()

# Step 3: Run LLL sweeps
print("[Step 3] Running LLL sweeps...")
sweep_results = {}   # label -> (lam_n, first_win_or_None, results)

# Run anchors first
for (p, b, n, lam, lam_n, G, label) in anchor_curves:
    print(f"\n  {label}: p={p}, n={n}, lam/n={lam_n:.4f}")
    res, fw = sweep_curve(p, b, n, lam, G, K1, M_RANGE, SEEDS)
    detail = " ".join(f"m{m}={w}/{t}" for m, (w, t) in sorted(res.items()))
    if fw:
        print(f"    → 3/3 first at m={fw}   [{detail}]")
    else:
        max_w = max(w for w, t in res.values())
        print(f"    → never 3/3 (max {max_w}/3)   [{detail}]")
    sweep_results[label] = (lam_n, fw, res)

# Run target curves
for t in LAM_TARGETS:
    curves = bucket[t]
    print(f"\n  [λ/n≈{t:.2f}]  ({len(curves)} curves found)")
    if not curves:
        print("    No curves found; skipping.")
        sweep_results[f"lam{t:.2f}_0"] = (t, None, {})
        continue
    for idx, (p, b, n, lam, lam_n) in enumerate(curves):
        label = f"lam{t:.2f}_{idx}"
        G = find_generator(p, b, n)
        if G is None:
            print(f"    {label}: generator not found, skip")
            continue
        K2 = math.isqrt(n) + 1
        eff = K1 * K2 / n
        print(f"    {label}: p={p}, n={n}, lam={lam}, lam/n={lam_n:.5f}, eff={eff:.4f}")
        res, fw = sweep_curve(p, b, n, lam, G, K1, M_RANGE, SEEDS)
        detail = " ".join(f"m{m}={w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if fw:
            print(f"    → 3/3 first at m={fw}   [{detail}]")
        else:
            max_w = max(w for w, t2 in res.values())
            print(f"    → never 3/3 (max {max_w}/3)   [{detail}]")
        sweep_results[label] = (lam_n, fw, res)

# Step 4: Summary table
print()
print("=" * 72)
print("SUMMARY TABLE (K1=72, m=10..25, seeds=3, 20-bit j=0 CM curves)")
print(f"{'label':<18} {'lam/n':>7} {'lam/n-0.333':>12} {'first 3/3':>10}  outcome")
print("-" * 65)

all_entries = []
for label, (lam_n, fw, res) in sorted(sweep_results.items(), key=lambda x: x[1][0]):
    dist = abs(lam_n - 1/3)
    outcome = f"SUCCESS m={fw}" if fw else "FAIL (never 3/3)"
    print(f"  {label:<16} {lam_n:>7.4f} {dist:>12.4f} {str(fw) if fw else 'never':>10}  {outcome}")
    all_entries.append((lam_n, fw))

# Find the transition
successes = [(lam_n, fw) for lam_n, fw in all_entries if fw is not None]
failures = [(lam_n, fw) for lam_n, fw in all_entries if fw is None]

print()
print("INTERPRETATION:")
if successes:
    min_success = min(lam_n for lam_n, _ in successes)
    print(f"  Lowest λ/n with 3/3 success: {min_success:.4f}")
if failures:
    max_failure = max(lam_n for lam_n, _ in failures)
    min_failure = min(lam_n for lam_n, _ in failures)
    print(f"  Failure range: λ/n ∈ [{min_failure:.4f}, {max_failure:.4f}]")
if successes and failures:
    s_below = max((lam_n for lam_n, _ in successes if lam_n < min(lam_n for lam_n, _ in failures)), default=None)
    f_above = min((lam_n for lam_n, _ in failures if lam_n > min_success), default=None) if successes else None
    if s_below:
        print(f"  Lower transition bracket: SUCCESS at {s_below:.4f}, FAIL at {min_failure:.4f}")
    if f_above:
        print(f"  Upper transition bracket: FAIL at {max_failure:.4f}, SUCCESS at {min_success:.4f}")

print()
print("June 2026 anchors for reference:")
print("  λ/n=0.016 → SUCCESS at m=10 (2026-06-21)")
print("  λ/n=0.039 → SUCCESS at m=20 (2026-06-21)")
print("  λ/n=0.07  → FAIL all m≤12 (2026-06-15, 12-bit curve p=2677)")
print("  λ/n=0.281 → FAIL (Effect A borderline)")
print("  λ/n=0.287 → SUCCESS at m=15 (Curve C, 2026-06-21)")
print()
print("Done.")
