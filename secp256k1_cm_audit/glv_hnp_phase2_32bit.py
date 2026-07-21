"""
GLV-aware HNP Phase 2 — 32-bit scaling test.

Extends glv_hnp_phase2_20bit.py to 32-bit j=0 prime-order curves.
Uses the same 2(2m+2)-dimensional lattice with column-diagonal balancing.

Threat model (BOTH-bounded GLV):
  k = k1 + lambda * k2  (mod n)
  k1 in [0, K1_BOUND)    (biased: small)
  k2 in [0, K2_BOUND)    (GLV domain: |k2| <= sqrt(n))

Research questions:
1. Does LLL succeed at 32-bit? (expected: yes, similar m_thresh to 20-bit)
2. How does wall-clock scale from 20-bit -> 32-bit?
3. Does BKZ rescue help if LLL fails?

Run: python3 glv_hnp_phase2_32bit.py
"""

import math
import random
import time
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic over F_p (short Weierstrass, a=0)
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

def find_generator(p, b, n, attempts=20000, seed=12345):
    rng = random.Random(seed)
    for _ in range(attempts):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: Eisenstein decomposition + sextic trace enumeration
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - ab + b^2 = p, a,b > 0. O(sqrt(p))."""
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
    """Six Frobenius traces for the j=0 sextic twists given Eisenstein (a,b)."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Smaller root of x^2 + x + 1 = 0 mod n (requires n ≡ 1 mod 3)."""
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    lam = min(r1, r2)
    return lam, n - 1 - lam

def identify_twist(p, n, max_b=500):
    """Find smallest b in [1, max_b] such that #(y^2=x^3+b over F_p)=n."""
    rng = random.Random(99)
    for b in range(1, max_b + 1):
        # Quick point check: find a point and test its order
        for _ in range(30):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b
                break  # this b doesn't have order n
    return None

# ---------------------------------------------------------------------------
# Find a suitable 32-bit curve
# ---------------------------------------------------------------------------

def find_32bit_curve(lam_ratio_lo=0.25, lam_ratio_hi=0.75, max_p=None):
    """
    Search for a 32-bit prime p with:
    - p ≡ 1 (mod 6)
    - n = #E prime, n ≡ 1 (mod 3)
    - lam_ratio_lo <= lambda/n <= lam_ratio_hi  (balanced GLV)
    Returns (p, b, n, lam, G) or None.
    """
    if max_p is None:
        max_p = 2**32
    print(f"Searching for 32-bit j=0 GLV curve (p in [2^31, 2^32))...")
    p = sympy.nextprime(2**31)
    found = 0
    t0 = time.time()
    while p < max_p:
        if p % 6 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2:
                        continue
                    if sympy.isprime(n_cand) and n_cand % 3 == 1:
                        lam, lam2 = glv_eigenvalue(n_cand)
                        if lam is None:
                            continue
                        ratio = lam / n_cand
                        if lam_ratio_lo <= ratio <= lam_ratio_hi:
                            found += 1
                            if found <= 3:
                                print(f"  Candidate {found}: p={p}, n={n_cand}, "
                                      f"lam={lam}, lam/n={ratio:.4f} "
                                      f"(elapsed {time.time()-t0:.1f}s)")
                            b = identify_twist(p, n_cand)
                            if b is None:
                                continue
                            G = find_generator(p, b, n_cand)
                            if G is None:
                                continue
                            print(f"  Selected: p={p}, n={n_cand}, lam={lam}, "
                                  f"lam/n={ratio:.4f}, b={b}")
                            return (p, b, n_cand, lam, G)
        p = sympy.nextprime(p)
    return None

# ---------------------------------------------------------------------------
# Signature generation (GLV k1+k2 both bounded)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# GLV lattice (2m+2 dimensional, column-diagonal balancing)
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_D = 1
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

    M = [[0] * dim for _ in range(dim)]
    # Rows 0..m-1: mod-n constraints (k1_i columns scaled by S_K1)
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

# ---------------------------------------------------------------------------
# Recovery: scan LLL-reduced rows for d
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

def run_experiment(curve_params, m, d_secret, k1_bound, seed=42,
                   use_bkz=False, bkz_beta=20, verbose=False):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    t0 = time.time()
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False, 0.0
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    d_rec = recover_d(reduced, m, n, S_KANNAN, d_secret)
    elapsed = time.time() - t0
    if verbose:
        print(f"    seed={seed}, m={m}: {'OK' if d_rec else 'FAIL'} ({elapsed:.2f}s)")
    return d_rec is not None, elapsed

# ---------------------------------------------------------------------------
# Sweep across m values and seeds
# ---------------------------------------------------------------------------

def sweep_curve(label, curve_params, k1_bound, m_range, seeds,
                use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve_params
    k2_bound = math.isqrt(n) + 1
    eff = k1_bound * k2_bound / n
    m_thresh = (math.ceil(math.log(n) / math.log(1.0 / eff))
                if eff < 1 else float('inf'))
    algo = f"BKZ(beta={bkz_beta})" if use_bkz else "LLL"
    print(f"\n{'='*70}")
    print(f"Curve: {label}  [{algo}]")
    print(f"  p={p} ({p.bit_length()}b), n={n} ({n.bit_length()}b), lam={lam}")
    print(f"  lam/n={lam/n:.4f}, K1={k1_bound}, K2={k2_bound} (~sqrt(n))")
    print(f"  eff={eff:.5f}, m_thresh≈{m_thresh:.1f}")
    print(f"{'='*70}")

    results = {}
    first_full = None
    for m in m_range:
        wins = 0
        total_t = 0.0
        for seed in seeds:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok, t = run_experiment(curve_params, m, d_trial, k1_bound, seed,
                                   use_bkz=use_bkz, bkz_beta=bkz_beta)
            wins += ok
            total_t += t
        results[m] = (wins, len(seeds))
        marker = " ← FULL" if wins == len(seeds) else ""
        marker2 = " (≥thresh)" if m >= m_thresh else ""
        avg_t = total_t / len(seeds)
        print(f"  m={m}: {wins}/{len(seeds)} recovered{marker2}  avg {avg_t:.2f}s/trial{marker}")
        if first_full is None and wins == len(seeds):
            first_full = m
    return results, first_full

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-aware HNP Phase 2 — 32-bit toy curve scaling test")
print("=" * 70)
print()

SEEDS = [42, 1234, 9999]

# ------ Part 1: Find a 32-bit curve -----------------------------------------
t_start = time.time()
params_32 = find_32bit_curve()
if params_32 is None:
    print("ERROR: No 32-bit curve found. Aborting.")
    raise SystemExit(1)

p32, b32, n32, lam32, G32 = params_32
k2_32 = math.isqrt(n32) + 1
# K1_BOUND: target eff ≈ 0.05 → K1 ≈ 0.05 * sqrt(n)
k1_32 = max(2, int(0.05 * math.sqrt(n32)))
print(f"\n32-bit curve selected (search took {time.time()-t_start:.1f}s):")
print(f"  p={p32}, b={b32}, n={n32} ({n32.bit_length()}b)")
print(f"  lam={lam32}, lam/n={lam32/n32:.4f}")
print(f"  K1_BOUND={k1_32}, K2_BOUND={k2_32}")
eff32 = k1_32 * k2_32 / n32
m_thresh32 = math.ceil(math.log(n32) / math.log(1.0 / eff32)) if eff32 < 1 else 99
print(f"  eff={eff32:.5f}, m_thresh≈{m_thresh32:.1f}")
print()

# ------ Part 2: LLL sweep ---------------------------------------------------
m_lo = max(3, m_thresh32 - 2)
m_hi = m_thresh32 + 5
res_lll, first_lll = sweep_curve(
    f"32-bit: y²=x³+{b32}/F_{p32}, n={n32}, lam={lam32}",
    params_32, k1_bound=k1_32, m_range=range(m_lo, m_hi + 1), seeds=SEEDS,
    use_bkz=False
)

# ------ Part 3: BKZ rescue if LLL fails at threshold -----------------------
need_bkz = (first_lll is None or first_lll > m_thresh32 + 2)
if need_bkz:
    print("\nLLL struggled — trying BKZ(beta=20):")
    res_bkz, first_bkz = sweep_curve(
        f"32-bit BKZ(20): y²=x³+{b32}/F_{p32}, n={n32}",
        params_32, k1_bound=k1_32, m_range=range(m_lo, m_hi + 1), seeds=SEEDS,
        use_bkz=True, bkz_beta=20
    )
else:
    res_bkz, first_bkz = None, None

# ------ Part 4: Reference comparison (reproduce 8-bit & 12-bit) ------------
print("\n" + "=" * 70)
print("Reference: 8-bit and 12-bit curves (from 2026-06-15 baseline)")
print("=" * 70)

from fpylll import IntegerMatrix, LLL

def find_ref_generator(p, b, n):
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

G_8bit  = find_ref_generator(211, 2, 199)
G_12bit = find_ref_generator(2557, 2, 2659)

if G_8bit and G_12bit:
    curve_8  = (211,  2, 199,  106,  G_8bit)
    curve_12 = (2557, 2, 2659, 1755, G_12bit)
    res_8,  first_8  = sweep_curve("8-bit:  y²=x³+2/F_211,  n=199,  lam=106",
                                   curve_8,  k1_bound=2,  m_range=range(3, 9),  seeds=SEEDS)
    res_12, first_12 = sweep_curve("12-bit: y²=x³+2/F_2557, n=2659, lam=1755",
                                   curve_12, k1_bound=8,  m_range=range(3, 10), seeds=SEEDS)
else:
    first_8 = first_12 = None
    print("  (could not reproduce reference generators)")

# ------ Summary -------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY — GLV-HNP Phase 2 scaling")
print("=" * 70)

def status(first_full, seeds):
    if first_full is None:
        return f"never {len(seeds)}/{len(seeds)} in sweep"
    return f"LLL {len(seeds)}/{len(seeds)} at m={first_full}"

if first_8  is not None: print(f"  8-bit  / n=199  (lam/n=0.53):  {status(first_8, SEEDS)}")
if first_12 is not None: print(f"  12-bit / n=2659 (lam/n=0.66):  {status(first_12, SEEDS)}")

ratio32 = lam32 / n32
print(f"  32-bit / n={n32} (lam/n={ratio32:.2f}):  {status(first_lll, SEEDS)}")
if res_bkz is not None:
    bstatus = f"BKZ(20) {len(SEEDS)}/{len(SEEDS)} at m={first_bkz}" if first_bkz else f"BKZ(20) never {len(SEEDS)}/{len(SEEDS)}"
    print(f"  32-bit / BKZ rescue:  {bstatus}")

print()
if first_lll is not None:
    print(f"CONCLUSION: GLV-HNP Phase 2 lattice attack WORKS at 32-bit.")
    print(f"  First full recovery at m={first_lll} (m_thresh={m_thresh32}).")
    print(f"  Ratio m/m_thresh = {first_lll/m_thresh32:.2f}x (expected ~1.5-2x from 20-bit experience).")
else:
    print(f"CONCLUSION: LLL failed at 32-bit up to m={m_hi}.")
    print(f"  m_thresh={m_thresh32}. May need higher m or BKZ.")
    if first_bkz is not None:
        print(f"  BKZ(20) rescues: {len(SEEDS)}/{len(SEEDS)} at m={first_bkz}.")
