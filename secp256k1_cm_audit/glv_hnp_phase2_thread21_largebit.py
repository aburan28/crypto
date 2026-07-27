"""
Thread 21: λ/n threshold study at 16-20 bit scale.

Motivation: the 11-13 bit study (Thread 20, 2026-07-27) showed non-monotone
LLL success vs λ/n. Hypothesis: at larger bit scale, law-of-large-numbers
averaging removes curve-specific noise and a monotone threshold emerges.

Method: Same as glv_hnp_phase2_lambda_threshold.py but bit_lo=16, bit_hi=20.

Run: python3 glv_hnp_phase2_thread21_largebit.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# Reuse helpers from the threshold script (inline for self-containment)

def modinv(a, m): return pow(a, -1, m)

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
    mc, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mc - i - 1), p)
        mc, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b_val = num // 2
                if b_val >= 0 and a * a - a * b_val + b_val * b_val == p:
                    return (a, b_val)
    return None

def j0_traces(a, b): return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
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

def find_b_for_order(p, n):
    for b in range(1, 13):
        rng = random.Random(b * 7)
        for _ in range(200):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    return b
                break
    return None

def gen_signatures_curve(p, b_coeff, n, lam, G, d_secret, m, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 50000:
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
    return M, S_KANNAN

def run_attack(p, b_coeff, n, lam, k1_bound, k2_bound, m_sigs, seed=42):
    rng = random.Random(seed + 99999)
    d_secret = rng.randint(1, n - 1)
    G = find_generator(p, b_coeff, n)
    if G is None: return False
    sigs = gen_signatures_curve(p, b_coeff, n, lam, G, d_secret, m_sigs, k1_bound, k2_bound, seed)
    if len(sigs) < m_sigs: return False
    M, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m_sigs + 2
    A_mat = IntegerMatrix.from_matrix(M)
    LLL.reduction(A_mat)
    for i in range(dim):
        row = [A_mat[i][j] for j in range(dim)]
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m_sigs]) % n
        if d_cand == d_secret: return True
        if ((-sign * row[m_sigs]) % n) == d_secret: return True
    return False

def make_bounds(n):
    k2 = math.isqrt(n) + 1
    k1 = max(2, int(0.15 * math.sqrt(n)))
    return k1, k2

# ---------------------------------------------------------------------------
# Main: 16-20 bit sweep
# ---------------------------------------------------------------------------

BINS_16_20 = [
    (0.04, 0.08),
    (0.08, 0.12),
    (0.12, 0.16),
    (0.16, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.36),
    (0.36, 0.42),
    (0.42, 0.50),
]

SEEDS = [42, 1234, 9999, 7, 314159]
M_SIGS = 14  # slightly more signatures than the 12-bit run (larger n → more needed)

print("=" * 70)
print("Thread 21: GLV-HNP λ/n threshold at 16-20 bit scale")
print(f"m={M_SIGS}, {len(SEEDS)} seeds per curve")
print("=" * 70)
print()

BIT_LO, BIT_HI = 16, 20
p_lo, p_hi = 2**BIT_LO, 2**BIT_HI

print(f"Searching j=0 prime-order curves in [{p_lo}, {p_hi}]...")
bins_db = {b: [] for b in range(len(BINS_16_20))}
found_count = {b: 0 for b in range(len(BINS_16_20))}

p_cur = sympy.nextprime(p_lo - 1)
fp_checked = 0
while p_cur < p_hi:
    fp_checked += 1
    if p_cur % 6 == 1:
        eis = eisenstein_decompose(p_cur)
        if eis is not None:
            a_e, b_e = eis
            for t in j0_traces(a_e, b_e):
                n_cand = p_cur + 1 - t
                if n_cand <= 1: continue
                if not sympy.isprime(n_cand): continue
                if n_cand % 3 != 1: continue
                lam_c, lam2_c = glv_eigenvalue(n_cand)
                if lam_c is None: continue
                r_min = min(lam_c, lam2_c) / n_cand
                lam_use = lam_c if lam_c / n_cand == r_min else lam2_c
                for bin_idx, (r_lo, r_hi) in enumerate(BINS_16_20):
                    if found_count[bin_idx] >= 1: continue
                    if r_lo <= r_min < r_hi:
                        b_cur = find_b_for_order(p_cur, n_cand)
                        if b_cur is None: continue
                        bins_db[bin_idx].append((p_cur, b_cur, n_cand, lam_use, r_min))
                        found_count[bin_idx] += 1
                        print(f"  Bin [{r_lo:.2f},{r_hi:.2f}): p={p_cur}, n={n_cand}, "
                              f"λ={lam_use}, λ/n={r_min:.5f}")
                        break
    if all(found_count[b] >= 1 for b in range(len(BINS_16_20))):
        break
    p_cur = sympy.nextprime(p_cur)

print(f"  Searched {fp_checked} primes.\n")

print(f"Attack sweep:")
print(f"{'Bin':20s}  {'n-bits':>6}  {'λ/n':>7}  {'LLL':>6}  flag")
print("-" * 55)

results = []
for bin_idx, (r_lo, r_hi) in enumerate(BINS_16_20):
    if not bins_db[bin_idx]:
        print(f"[{r_lo:.2f},{r_hi:.2f})  --  --  N/A")
        results.append(None)
        continue
    p_c, b_c, n_c, lam_c, r_c = bins_db[bin_idx][0]
    k1_c, k2_c = make_bounds(n_c)
    eff = k1_c * k2_c / n_c
    wins = sum(1 for seed in SEEDS
               if run_attack(p_c, b_c, n_c, lam_c, k1_c, k2_c, M_SIGS, seed=seed))
    flag = "PASS" if wins >= 3 else ("MARGINAL" if wins >= 1 else "FAIL")
    print(f"[{r_lo:.2f},{r_hi:.2f})  {n_c.bit_length():>6}  {r_c:.5f}  "
          f"{wins}/{len(SEEDS)}  ({flag})  [K1={k1_c},K2={k2_c},eff={eff:.3f}]")
    results.append((r_lo, r_hi, p_c, r_c, wins))

print()
print("=" * 70)
print("Threshold analysis:")
fail_max = None
pass_min = None
# Only consider consecutive monotone region for threshold
for item in results:
    if item is None: continue
    r_lo, r_hi, p_c, r_c, wins = item
    if wins < 3:
        if fail_max is None or r_c > fail_max:
            fail_max = r_c
    else:
        if pass_min is None or r_c < pass_min:
            pass_min = r_c

if fail_max is not None and pass_min is not None:
    print(f"  Highest FAIL:  λ/n = {fail_max:.5f}")
    print(f"  Lowest PASS:   λ/n = {pass_min:.5f}")
    if fail_max < pass_min:
        print(f"  Monotone threshold estimate: λ/n ≈ {(fail_max+pass_min)/2:.5f}")
    else:
        print(f"  NON-MONOTONE behavior detected (fail_max > pass_min)")
        print(f"  → Same non-monotone pattern as 11-13 bit study (Thread 20)")
elif fail_max is None:
    print("  All passed — threshold below tested range")
elif pass_min is None:
    print("  All failed — threshold above tested range")

print("\nDone.")
