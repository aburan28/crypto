"""
GLV-aware HNP Phase 2 — 32-bit toy curve verification.

Scales the GLV-HNP lattice attack to a 32-bit prime-order j=0 curve.
Uses Eisenstein CM theory (O(sqrt(p))) to find candidate group orders.

Attack model: k = k1 + lam*k2  mod n, with k1 small (biased).
              k2 in [0, K2_BOUND) where K2_BOUND = ceil(sqrt(n)) ~ 2^16.
              k1 in [0, K1_BOUND), bias chosen so eff = K1*K2/n ~ 0.05.
              m_thresh ≈ log(n) / log(1/eff) ≈ 32 / 4.3 ≈ 7.5 -> test m = 8..14.

Lattice: column-diagonal balanced (2m+2) x (2m+2) integer matrix.
Scales: S_K1 = n // K1, S_D = 1, S_K2 = n // K2, S_KANNAN = n.

Run: python3 glv_hnp_phase2_32bit.py
"""

import math
import random
import sys

try:
    import sympy
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

try:
    from fpylll import IntegerMatrix, LLL, BKZ
    HAS_FPYLLL = True
except ImportError:
    print("ERROR: fpylll not installed.  Run: pip3 install fpylll cysignals")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass a=0, y^2 = x^3 + b over F_p)
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

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(100000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition and trace enumeration
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) >= 0 with a^2 - a*b + b^2 = p.  O(sqrt(p))."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    """Six Frobenius traces for j=0 sextic twists, given Eisenstein (a,b)."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Compute the GLV eigenvalue: smaller root of x^2+x+1=0 mod n."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

def is_prime_miller_rabin(n, k=20):
    """Deterministic for n < 3.3e24 using known witnesses; probabilistic otherwise."""
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    if HAS_SYMPY:
        return sympy.isprime(n)
    # Fallback: Miller-Rabin with fixed witnesses (deterministic for n < 3.3e24)
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in witnesses:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def find_32bit_curve():
    """
    Search for a 32-bit prime p with:
      p ≡ 1 mod 3,  n = #E(F_p) prime,  n ≡ 1 mod 3,  lam/n ∈ [0.25, 0.75].
    Uses Eisenstein CM to enumerate n candidates in O(sqrt(p)) per prime.
    """
    print("Searching for a 32-bit j=0 GLV curve (this may take ~30s)...")
    # Start slightly above 2^31 to get genuine 32-bit candidates
    target_bits = 32
    start = 2 ** (target_bits - 1) + 1
    if HAS_SYMPY:
        p = sympy.nextprime(start)
    else:
        p = start | 1
        while not is_prime_miller_rabin(p): p += 2

    checked = 0
    while p < 2 ** target_bits:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                traces = j0_traces(a_e, b_e)
                for t in traces:
                    n_cand = p + 1 - t
                    if n_cand < 2 ** (target_bits - 2):
                        continue
                    if not is_prime_miller_rabin(n_cand):
                        continue
                    if n_cand % 3 != 1:
                        continue
                    lam, lam2 = glv_eigenvalue(n_cand)
                    if lam is None:
                        continue
                    ratio = lam / n_cand
                    if not (0.25 <= ratio <= 0.75):
                        continue
                    # Find curve parameter b
                    found_b = None
                    for b_try in range(1, 200):
                        rng_tmp = random.Random(42)
                        ok = False
                        for _ in range(50):
                            x = rng_tmp.randint(0, p - 1)
                            rhs = (pow(x, 3, p) + b_try) % p
                            y = tonelli_shanks(rhs, p)
                            if y is not None and y != 0:
                                P = (x, y)
                                if ec_mul(P, n_cand, p) is None:
                                    found_b = b_try
                                    ok = True
                                    break
                        if ok:
                            break
                    if found_b is None:
                        continue
                    G = find_generator(p, found_b, n_cand)
                    if G is None:
                        continue
                    print(f"  Found: p={p}, b={found_b}, n={n_cand} ({n_cand.bit_length()}b), "
                          f"lam={lam}, lam/n={ratio:.4f}")
                    return (p, found_b, n_cand, lam, G)
        checked += 1
        if HAS_SYMPY:
            p = sympy.nextprime(p)
        else:
            p += 2
            while not is_prime_miller_rabin(p): p += 2

    return None

# ---------------------------------------------------------------------------
# Signature generation (k1 biased small, k2 in GLV range)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

# ---------------------------------------------------------------------------
# Lattice construction with column-diagonal balancing
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound          # scale k1 columns up
    S_D = 1                        # d is unscaled (d < n ~ 2^32)
    S_K2 = max(1, n // k2_bound)  # scale k2 columns up
    S_KANNAN = n                   # Kannan marker

    M = [[0] * dim for _ in range(dim)]
    # Rows 0..m-1: mod-n constraint
    for i in range(m):
        M[i][i] = n * S_K1
    # Row m: d-row
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    # Rows m+1..2m: k2-rows
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    # Row 2m+1: Kannan
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN

    return M, S_K1, S_D, S_K2, S_KANNAN

# ---------------------------------------------------------------------------
# Recovery: scan LLL-reduced rows for correct d
# ---------------------------------------------------------------------------

def recover_d(M_reduced, m, n, S_KANNAN, d_secret=None):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_secret is not None:
            if d_cand == d_secret:
                return d_cand
        else:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Single experiment
# ---------------------------------------------------------------------------

def run_experiment(curve_params, m, d_secret, k1_bound, seed=42,
                   use_bkz=False, bkz_beta=20):
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
# Sweep over m and seeds
# ---------------------------------------------------------------------------

def sweep(label, curve_params, k1_bound, m_range, seeds, use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    algo = f"BKZ(beta={bkz_beta})" if use_bkz else "LLL"
    if eff < 1.0:
        m_thresh = math.ceil(math.log(n) / math.log(1.0 / eff))
    else:
        m_thresh = float('inf')

    print(f"\n{'='*65}")
    print(f"Curve: {label}  [{algo}]")
    print(f"  p={p}, n={n} ({n.bit_length()}b), lam={lam}")
    print(f"  lam/n={lam/n:.4f}, K1={k1_bound}, K2={k2_bound}")
    print(f"  eff={eff:.5f}, m_thresh≈{m_thresh:.1f}")
    print(f"{'='*65}")

    results = {}
    first_perfect = None
    for m in m_range:
        wins = 0
        for seed in seeds:
            rng_d = random.Random(seed + 7777)
            d_trial = rng_d.randint(1, n - 1)
            ok = run_experiment(curve_params, m, d_trial, k1_bound,
                                seed=seed, use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
        results[m] = (wins, len(seeds))
        thresh_marker = " (≥thresh)" if m >= m_thresh else ""
        print(f"  m={m}: {wins}/{len(seeds)} recovered{thresh_marker}")
        if wins == len(seeds) and first_perfect is None:
            first_perfect = m
    return results, first_perfect

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 65)
print("GLV-aware HNP Phase 2 — 32-bit toy curve")
print("=" * 65)

SEEDS_LIGHT = [42, 1234, 9999]
SEEDS_FULL  = [42, 1234, 9999, 7, 314159]

# --- Find 32-bit curve -------------------------------------------------
params_32 = find_32bit_curve()
if params_32 is None:
    print("ERROR: No suitable 32-bit curve found.  Aborting.")
    sys.exit(1)

p32, b32, n32, lam32, G32 = params_32
k2_32 = math.isqrt(n32) + 1

# K1_BOUND: target eff = K1*K2/n ~ 0.05 → K1 ≈ 0.05 * n / K2 = 0.05 * sqrt(n)
k1_32 = max(2, int(0.05 * math.sqrt(n32)))
eff_32 = k1_32 * k2_32 / n32

print(f"\n32-bit curve found:")
print(f"  p={p32}, b={b32}, n={n32} ({n32.bit_length()}b)")
print(f"  lam={lam32}, lam/n={lam32/n32:.4f}")
print(f"  K1_BOUND={k1_32}, K2_BOUND={k2_32}")
print(f"  eff={eff_32:.5f}")
if eff_32 < 1.0:
    m_thresh_32 = math.ceil(math.log(n32) / math.log(1.0 / eff_32))
    print(f"  m_thresh≈{m_thresh_32:.1f}")
print()

# --- LLL sweep (light: 3 seeds, m = thresh-2..thresh+6) ---------------
if eff_32 < 1.0:
    m_lo = max(3, m_thresh_32 - 2)
    m_hi = m_thresh_32 + 7
else:
    m_lo, m_hi = 6, 14

res_lll, fp_lll = sweep(
    f"32-bit: y²=x³+{b32}/F_{p32}, n={n32}, lam={lam32}",
    params_32, k1_bound=k1_32,
    m_range=range(m_lo, m_hi + 1),
    seeds=SEEDS_LIGHT,
    use_bkz=False
)

# --- BKZ(20) sweep on same curve if LLL did not achieve 3/3 ----------
lll_3of3 = fp_lll is not None
if not lll_3of3:
    print("\nLLL did not reach 3/3; trying BKZ(beta=20) ...")
    res_bkz, fp_bkz = sweep(
        f"32-bit: y²=x³+{b32}/F_{p32}, n={n32}, lam={lam32}",
        params_32, k1_bound=k1_32,
        m_range=range(m_lo, m_hi + 1),
        seeds=SEEDS_LIGHT,
        use_bkz=True, bkz_beta=20
    )
else:
    fp_bkz = None

# --- Full 5-seed sweep at the first-3/3 m (if found) -----------------
if fp_lll is not None:
    print(f"\n--- Full 5-seed confirmation at m={fp_lll} (first 3/3 LLL) ---")
    wins5 = 0
    for seed in SEEDS_FULL:
        rng_d = random.Random(seed + 7777)
        d_trial = rng_d.randint(1, n32 - 1)
        ok = run_experiment(params_32, fp_lll, d_trial, k1_32, seed=seed)
        wins5 += ok
    print(f"  m={fp_lll}: {wins5}/5 seeds recovered d")
elif fp_bkz is not None:
    print(f"\n--- Full 5-seed confirmation at m={fp_bkz} (first 3/3 BKZ-20) ---")
    wins5 = 0
    for seed in SEEDS_FULL:
        rng_d = random.Random(seed + 7777)
        d_trial = rng_d.randint(1, n32 - 1)
        ok = run_experiment(params_32, fp_bkz, d_trial, k1_32,
                            seed=seed, use_bkz=True, bkz_beta=20)
        wins5 += ok
    print(f"  m={fp_bkz}: {wins5}/5 seeds recovered d")

# --- Summary ----------------------------------------------------------
print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"32-bit curve: p={p32}, n={n32} ({n32.bit_length()}b), lam/n={lam32/n32:.4f}")
print(f"K1={k1_32}, K2={k2_32}, eff={eff_32:.5f}")
if fp_lll is not None:
    print(f"LLL:     3/3 first at m={fp_lll}")
    if fp_lll is not None and wins5:
        print(f"         5/5 confirmed at m={fp_lll}: {'YES' if wins5 == 5 else wins5 + '/5'}")
elif fp_bkz is not None:
    print(f"LLL:     never 3/3 in sweep")
    print(f"BKZ(20): 3/3 first at m={fp_bkz}")
else:
    print("LLL:     never 3/3 in sweep")
    print("BKZ(20): never 3/3 in sweep  [may need larger m or BKZ-40]")
print("\nDone.")
