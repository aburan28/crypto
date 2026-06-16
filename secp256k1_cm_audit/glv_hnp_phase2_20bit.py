"""
GLV-HNP Phase 2: 20-bit scaling test using CM curve finding.

Extends the scaling study to 20-bit prime-order j=0 curves.
Uses Eisenstein CM theory to find candidate group orders without
brute-force point counting (O(sqrt(p)) instead of O(p) per prime).

Also runs BKZ (beta=20) on the small-lambda failure curve (p=2677, lam=185)
from the 2026-06-15 session to check if BKZ rescues the attack.

Run: python3 glv_hnp_phase2_20bit.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
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
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_point(p, b):
    """Find a random point on y^2 = x^3 + b over F_p."""
    rng = random.Random(42)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            assert (y * y - pow(x, 3, p) - b) % p == 0
            return (x, y)
    return None

def find_generator(p, b, n):
    """Find a generator of E(F_p): y^2=x^3+b of order n (n prime)."""
    rng = random.Random(12345)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition for j=0 curves over F_p
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """
    Find (a,b) with a^2 - a*b + b^2 = p (Eisenstein norm), a,b >= 0.
    Uses the quadratic-formula trick: for each candidate a, solve
    b^2 - a*b + (a^2 - p) = 0 => b = (a +/- sqrt(4p - 3a^2)) / 2.
    Runs in O(sqrt(p)).
    """
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    """
    The 6 Frobenius traces for the 6 sextic twists of j=0 over F_p,
    given Eisenstein decomposition a^2 - a*b + b^2 = p.

    The six associates of pi = a + b*omega (omega = primitive 3rd root of unity)
    have real parts {a - b/2, -a + b/2, -(a+b)/2, (a+b)/2, b - a/2, a/2 - b}.
    Traces = 2 * (real part).
    """
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """
    Compute the GLV eigenvalue lambda: the smaller root of x^2 + x + 1 = 0 mod n.
    Requires n ≡ 1 (mod 3).
    Returns (lam, n-1-lam) sorted so lam <= n//2.
    """
    # sqrt(-3 mod n)
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    # lam = (-1 + sq) / 2 mod n, or (-1 - sq) / 2 mod n
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    # Verify
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0, f"GLV eigenvalue check failed for n={n}"
    lam = min(r1, r2)
    return lam, n - 1 - lam

def identify_b_for_n(p, n, a_eis, b_eis):
    """
    Find the curve parameter b (twist) such that #(y^2=x^3+b over F_p) = n.
    Tries b=1,2,...,6 (covering all 6 twist classes), verifying by checking
    a random point.
    """
    # The 6 sextic twist classes are b * c^6 for c in F_p^*.
    # Representatives: b=1, g, g^2, g^3, g^4, g^5 where g is a primitive root mod p.
    # But for small checks, trying b=1..6 usually suffices (they're often in different classes).
    for b in range(1, p):
        # Quick check: find one point and verify order
        P = find_point(p, b)
        if P is None:
            continue
        Q = ec_mul(P, n, p)
        if Q is None:
            return b
    return None

def find_20bit_curve():
    """
    Search for a 20-bit prime p with:
    - p ≡ 1 (mod 3)   (j=0 GLV structure)
    - n = #E(F_p) prime, n ≡ 1 (mod 3)   (GLV eigenvalue exists)
    - lambda/n ∈ [0.25, 0.75]   (good bias regime)
    Returns (p, b, n, lam, G).
    """
    print("Searching for a 20-bit j=0 GLV curve...")
    p = sympy.nextprime(2**19)
    count = 0
    while p < 2**20:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                traces = j0_traces(a_e, b_e)
                for t in traces:
                    n_cand = p + 1 - t
                    if n_cand < 2:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam, lam2 = glv_eigenvalue(n_cand)
                        if lam is None:
                            continue
                        ratio = lam / n_cand
                        if 0.25 <= ratio <= 0.75:
                            # Found a candidate. Now find the twist b.
                            count += 1
                            if count <= 3:
                                print(f"  Candidate: p={p}, n={n_cand}, lam={lam}, "
                                      f"lam/n={ratio:.3f}")
                            # Find b: try small values
                            found_b = None
                            for b_try in range(1, 200):
                                P = find_point(p, b_try)
                                if P is None:
                                    continue
                                Q = ec_mul(P, n_cand, p)
                                if Q is None:
                                    found_b = b_try
                                    break
                            if found_b is None:
                                continue
                            G = find_generator(p, found_b, n_cand)
                            if G is None:
                                continue
                            return (p, found_b, n_cand, lam, G)
        p = sympy.nextprime(p)
    return None

# ---------------------------------------------------------------------------
# GLV decomposition (brute force for small n)
# ---------------------------------------------------------------------------

def glv_decompose_bf(k_full, n, lam, k2_bound):
    """Find (k1,k2) with k_full ≡ k1 + lam*k2 (mod n), k2 in [0,k2_bound)."""
    for k2 in range(k2_bound):
        k1 = (k_full - lam * k2) % n
        if k1 < k2_bound * 3:   # k1 small relative to k2_bound
            return (k1, k2)
    return None

# ---------------------------------------------------------------------------
# Signature generation
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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
# Build GLV lattice (column-diagonal scaling from 2026-06-15)
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
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
# Single LLL experiment
# ---------------------------------------------------------------------------

def run_experiment(curve_params, m, d_secret, k1_bound, seed=42, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def sweep_curve(label, curve_params, k1_bound, m_range, seeds, verbose=True,
                use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    if eff >= 1.0:
        m_thresh = float('inf')
    else:
        m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))

    algo = f"BKZ(beta={bkz_beta})" if use_bkz else "LLL"
    if verbose:
        print(f"\n{'='*65}")
        print(f"Curve: {label}  [{algo}]")
        print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}")
        print(f"  lam/n={lam/n:.4f}, K1={k1_bound}, K2={k2_bound}")
        print(f"  eff={eff:.5f}, m_thresh≈{m_thresh:.1f}")
        print(f"{'='*65}")

    results = {}
    for m in m_range:
        wins = 0
        for seed in seeds:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok = run_experiment(curve_params, m, d_trial, k1_bound, seed,
                                use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
        marker = " (≥thresh)" if m >= m_thresh else ""
        if verbose:
            print(f"  m={m}: {wins}/{len(seeds)} recovered{marker}")
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 65)
print("GLV-aware HNP Phase 2 — 20-bit scaling + BKZ rescue")
print("=" * 65)

SEEDS = [42, 1234, 9999]

# ---- Part 1: Find a 20-bit curve ------------------------------------------
params_20 = find_20bit_curve()
if params_20 is None:
    print("ERROR: No suitable 20-bit curve found. Aborting Part 1.")
else:
    p20, b20, n20, lam20, G20 = params_20
    print(f"\nSelected 20-bit curve:")
    print(f"  p={p20}, b={b20}, n={n20} ({n20.bit_length()}b), lam={lam20}")
    print(f"  lam/n={lam20/n20:.4f}")
    print(f"  G={G20}")

    # Compute K1_BOUND, K2_BOUND for 20-bit
    k2_20 = math.isqrt(n20) + 1
    # K1_BOUND: want eff ~ 0.03..0.10 so that m_thresh ~ 4-6
    # eff = K1 * K2 / n ~ K1 * sqrt(n) / n = K1 / sqrt(n)
    # Target eff=0.05: K1 = 0.05 * sqrt(n) ≈ 0.05 * 1024 ≈ 51
    k1_20 = max(2, int(0.05 * math.sqrt(n20)))
    print(f"  K1_BOUND={k1_20}, K2_BOUND={k2_20}")

    res_20 = sweep_curve(
        f"20-bit: y²=x³+{b20}/F_{p20}, n={n20}, lam={lam20}",
        params_20, k1_bound=k1_20, m_range=range(3, 10), seeds=SEEDS
    )

# ---- Part 2: BKZ rescue on the small-lambda failure (p=2677) ---------------
print("\n" + "=" * 65)
print("Part 2: BKZ rescue on small-lambda failure curve (p=2677, lam=185)")
print("=" * 65)

# Reconstruct curve from 2026-06-15 (p=2677, n=2647, lam=185)
p_fail = 2677
b_fail = 2
n_fail = 2647
lam_fail = 185
print(f"  p={p_fail}, n={n_fail}, lam={lam_fail}, lam/n={lam_fail/n_fail:.4f}")

# Rebuild generator (same as glv_hnp_phase2_scaling.py)
G_fail = None
rng_tmp = random.Random(12345)
for _ in range(10000):
    x = rng_tmp.randint(0, p_fail - 1)
    rhs = (pow(x, 3, p_fail) + b_fail) % p_fail
    y = tonelli_shanks(rhs, p_fail)
    if y is not None and y != 0:
        P = (x, y)
        if ec_mul(P, n_fail, p_fail) is None:
            G_fail = P
            break

if G_fail is None:
    print("ERROR: Could not find generator for p=2677 curve")
else:
    params_fail = (p_fail, b_fail, n_fail, lam_fail, G_fail)
    k2_fail = math.isqrt(n_fail) + 1  # ≈52

    # LLL sweep (baseline from 2026-06-15: always fails at m≤14)
    print("\n--- LLL baseline (known to fail at all m for lam/n=0.07) ---")
    res_lll_fail = sweep_curve(
        f"12-bit/2677 LLL (lam/n=0.07, EXPECTED FAIL)",
        params_fail, k1_bound=8, m_range=range(5, 13), seeds=SEEDS,
        use_bkz=False
    )

    # BKZ with beta=20
    print("\n--- BKZ(beta=20) on the same small-lambda curve ---")
    res_bkz20 = sweep_curve(
        f"12-bit/2677 BKZ(beta=20) (lam/n=0.07)",
        params_fail, k1_bound=8, m_range=range(5, 13), seeds=SEEDS,
        use_bkz=True, bkz_beta=20
    )

    # BKZ with beta=40
    print("\n--- BKZ(beta=40) on the same small-lambda curve ---")
    res_bkz40 = sweep_curve(
        f"12-bit/2677 BKZ(beta=40) (lam/n=0.07)",
        params_fail, k1_bound=8, m_range=range(5, 13), seeds=SEEDS,
        use_bkz=True, bkz_beta=40
    )

# ---- Part 3: Reference 8-bit and 12-bit (from 2026-06-15) -----------------
print("\n" + "=" * 65)
print("Part 3: Reference curves (8-bit and 12-bit/2557) — LLL")
print("=" * 65)

# 8-bit (p=211)
G0 = find_generator(211, 2, 199)
curve0 = (211, 2, 199, 106, G0)
res0 = sweep_curve("8-bit: y²=x³+2/F_211, n=199, lam=106",
                   curve0, k1_bound=2, m_range=range(3, 8), seeds=SEEDS)

# 12-bit (p=2557)
G1 = find_generator(2557, 2, 2659)
curve1 = (2557, 2, 2659, 1755, G1)
res1 = sweep_curve("12-bit/2557: y²=x³+2/F_2557, n=2659, lam=1755",
                   curve1, k1_bound=8, m_range=range(3, 9), seeds=SEEDS)

# ---- Summary -----------------------------------------------------------------
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)

results_list = [
    ("8-bit/199 (lam/n=0.53)", res0, SEEDS),
    ("12-bit/2557 (lam/n=0.66)", res1, SEEDS),
]
if params_20 is not None:
    results_list.append((f"20-bit/{n20} (lam/n={lam20/n20:.2f})", res_20, SEEDS))

for label, results, seeds in results_list:
    first5 = None
    for m, (wins, total) in sorted(results.items()):
        if wins == total:
            first5 = m
            break
    if first5:
        print(f"  {label}: LLL 3/3 at m={first5}")
    else:
        print(f"  {label}: LLL never 3/3 in sweep")

if G_fail is not None:
    print(f"\n  12-bit/2677 (lam/n=0.07, small-lambda FAILURE):")
    for algo, results in [("LLL", res_lll_fail), ("BKZ(20)", res_bkz20), ("BKZ(40)", res_bkz40)]:
        first3 = None
        for m, (wins, total) in sorted(results.items()):
            if wins == total:
                first3 = m
                break
        status = f"3/3 at m={first3}" if first3 else "never 3/3"
        print(f"    {algo}: {status}")

print("\nDone.")
