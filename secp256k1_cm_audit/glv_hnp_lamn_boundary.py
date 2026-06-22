"""
GLV-HNP Phase 2: Effect A boundary — λ/n proximity to 1/3 obstruction map.

Objective: characterise the precise λ/n range that triggers the near-linear
dependence obstruction.  From 2026-06-21 we know:
  - λ/n ≈ 0.327-0.340 → fully obstructed (0/3 at all m=10..20)
  - λ/n ≈ 0.016-0.039 → trivially unobstructed (Effect B, small λ)
  - λ/n ≈ 0.281-0.387 borderline (some fail at m≤20)
  - secp256k1: λ/n ≈ 0.326 → obstructed

Today: find 20-bit j=0 CM curves at λ/n targets
    {0.20, 0.23, 0.27, 0.30, 0.33, 0.35, 0.37, 0.40, 0.43, 0.47}
covering the range 0.20..0.50 (note 1/3≈0.333 is in the middle).
Run K1=72 LLL at m=10..25, 3 seeds each.

The hypothesis (Effect A mechanism): when λ/n ≈ k/3 for k=1 (or k=2, k=3...),
the GLV rows M[m+1+i][i] = -λ·S_K1 create a near-linear dependence with the
primary rows M[i][i] = n·S_K1 (since -λ·S_K1 ≈ -(1/3)·n·S_K1).
If the hypothesis is correct, we should see failure whenever:
    |3λ mod n| / n < threshold  (equivalently: λ/n within ε of 1/3)
and success outside that band.

Run: python3 glv_hnp_lamn_boundary.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (shared with delta_threshold.py)
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
# CM theory
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

def delta_ratio(lam, n):
    x = (3 * lam) % n
    return min(x, n - x) / n

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

# We want curves at these λ/n ratios (tolerance ±0.015).
# Note: by convention, lam = min(r1,r2) ≤ n/2, so λ/n ∈ (0, 0.5].
LAM_TARGETS = [0.20, 0.23, 0.27, 0.30, 0.33, 0.35, 0.37, 0.40, 0.43, 0.47]
LAM_TOL = 0.015          # ±1.5% of n
MAX_PER_TARGET = 2

def find_curves_by_lamn(p_lo=2**19, p_hi=2**20, targets=LAM_TARGETS,
                         tol=LAM_TOL, max_per=MAX_PER_TARGET):
    """
    Scan 20-bit j=0 CM primes and bucket by λ/n proximity to each target.
    Returns dict: target -> list of (p, b, n, lam, lam_n, delta_r).
    """
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
                    dr = delta_ratio(lam, n_cand)
                    for t in targets:
                        if abs(lam_n - t) <= tol and len(found[t]) < max_per:
                            b_param = find_b_for_n(p, n_cand)
                            if b_param is None:
                                continue
                            found[t].append((p, b_param, n_cand, lam, lam_n, dr))
                            print(f"  [lam/n≈{t:.2f}] p={p}, n={n_cand}, lam={lam}, "
                                  f"lam/n={lam_n:.4f}, δ/n={dr:.4f}")
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    print(f"  (scanned {scanned} primes in [{p_lo}, {p_hi}])")
    return found

# ---------------------------------------------------------------------------
# GLV signatures + lattice attack
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m, K1, seed):
    K2 = math.isqrt(n) + 1
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

def sweep(p, b_param, n, lam, G, K1, m_range, seeds):
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_lll(p, b_param, n, lam, G, m, K1, s))
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Curve C extended test (λ/n=0.287, m=10..30)
# ---------------------------------------------------------------------------

CURVE_C = dict(p=624517, b=None, n=622957, lam=178615, lam_n=0.2866, delta_r=0.140)

def run_curve_c_extended(K1=72, seeds=(42, 1234, 9999), m_hi=30):
    print("=" * 70)
    print(f"Curve C extended: m=10..{m_hi}, K1={K1}")
    print("=" * 70)
    p, n, lam = CURVE_C['p'], CURVE_C['n'], CURVE_C['lam']
    b = find_b_for_n(p, n)
    if b is None:
        print("  ERROR: no b found for Curve C")
        return None
    G = find_generator(p, b, n)
    if G is None:
        print("  ERROR: no generator for Curve C")
        return None
    dr = delta_ratio(lam, n)
    lam_n = lam / n
    print(f"  p={p}, n={n}, lam={lam}, lam/n={lam_n:.4f}, δ/n={dr:.4f}")
    results = sweep(p, b, n, lam, G, K1, range(10, m_hi + 1), seeds)
    first_win = next((m for m, (w, t) in sorted(results.items()) if w == t), None)
    # Print per-m results
    for m, (w, t) in sorted(results.items()):
        marker = " ← first 3/3" if m == first_win else ""
        print(f"  m={m:2d}: {w}/{t}{marker}")
    if first_win:
        print(f"  → 3/3 first at m={first_win}")
    else:
        max_w = max(w for w, t in results.values())
        print(f"  → never 3/3 (max {max_w}/3)")
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: Effect A boundary (λ/n proximity to 1/3)")
print("=" * 70)
print()

K1 = 72
SEEDS = [42, 1234, 9999]
M_RANGE = range(10, 26)     # m=10..25

# ---- Part 1: Curve C extended test ----------------------------------------
curve_c_results = run_curve_c_extended(K1=K1, seeds=SEEDS, m_hi=30)
print()

# ---- Part 2: Scan curves by λ/n -------------------------------------------
print("=" * 70)
print("Scanning 20-bit j=0 CM curves by λ/n target")
print("-" * 40)
bucket = find_curves_by_lamn()
print()

# ---- Part 3: LLL sweep -------------------------------------------------------
print("=" * 70)
print(f"LLL sweep: K1={K1}, m=10..25, seeds×3")
print("-" * 40)

sweep_results = {}    # (lam_n_actual, p, n) -> (first_win_m_or_None, res_dict)

for t in LAM_TARGETS:
    curves = bucket[t]
    dist_to_third = abs(t - 1/3)
    print(f"\n[λ/n≈{t:.2f}, |λ/n - 1/3|={dist_to_third:.4f}]  ({len(curves)} curves found)")
    if not curves:
        print("  No curves found; skipping.")
        continue
    for (p, b, n, lam, lam_n, dr) in curves:
        K2 = math.isqrt(n) + 1
        eff = K1 * K2 / n
        print(f"  p={p}, n={n}, lam={lam}, lam/n={lam_n:.4f}, δ/n={dr:.4f}, eff={eff:.4f}")
        G = find_generator(p, b, n)
        if G is None:
            print(f"  ERROR: generator not found; skip")
            continue
        res = sweep(p, b, n, lam, G, K1, M_RANGE, SEEDS)
        first_win = next((m for m, (w, t2) in sorted(res.items()) if w == t2), None)
        detail = " ".join(f"m{m}={w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if first_win:
            print(f"  → 3/3 first at m={first_win}   [{detail}]")
        else:
            max_wins = max(w for w, t2 in res.values())
            print(f"  → never 3/3 (max {max_wins}/3)   [{detail}]")
        sweep_results[(lam_n, p, n)] = (first_win, res)

# ---- Part 4: Summary table ---------------------------------------------------
print()
print("=" * 70)
print("SUMMARY: λ/n vs LLL outcome at K1=72")
print(f"{'λ/n':<8} {'|λ/n-1/3|':<12} {'δ/n':<8} {'first 3/3':<12} outcome")
print("-" * 60)

for t in LAM_TARGETS:
    dist = abs(t - 1/3)
    for (p, b, n, lam, lam_n, dr) in bucket[t]:
        key = (lam_n, p, n)
        if key not in sweep_results:
            continue
        fw, _ = sweep_results[key]
        if fw:
            outcome = f"SUCCESS at m={fw}"
        else:
            outcome = "OBSTRUCTED"
        print(f"  {lam_n:.4f}  {dist:.4f}       {dr:.4f}  {str(fw) if fw else 'never':<12} {outcome}")

# Curve C reference
if curve_c_results:
    dr = delta_ratio(CURVE_C['lam'], CURVE_C['n'])
    lam_n = CURVE_C['lam_n']
    dist = abs(lam_n - 1/3)
    fw = next((m for m, (w, t) in sorted(curve_c_results.items()) if w == t), None)
    outcome = f"SUCCESS at m={fw}" if fw else "OBSTRUCTED"
    print(f"  {lam_n:.4f}  {dist:.4f}       {dr:.4f}  {str(fw) if fw else 'never':<12} {outcome}  [Curve C anchor]")

# secp256k1
print(f"  0.3257  {'|0.3257-0.333|':<12}0.023    never        OBSTRUCTED [secp256k1, 256-bit]")

print()
print("Done.")
