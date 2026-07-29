"""
GLV-HNP Phase 2: λ/n threshold bisection.

Objective: bisect the λ/n boundary where LLL transitions from success to failure.

Known data from prior sessions (autolab 2026-07-26):
  - p=524347, n=523969, λ=177902, λ/n=0.34 → LLL 3/3 at m=9 (SUCCESS)
  - p=2677,   n=2647,   λ=185,    λ/n=0.07 → LLL 0/3 at all m (FAIL)

Hypothesis: there is a sharp threshold λ/n = τ* in (0.07, 0.34) below which LLL
fails. This script locates τ* by testing curves at λ/n ∈ {0.07, 0.10, 0.13, 0.16,
0.19, 0.22, 0.25, 0.28, 0.31, 0.34} using 13-bit j=0 CM curves (fast LLL).

Also tests: does BKZ(beta=20) push below τ*?

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic (short Weierstrass, a=0)
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

def find_generator(p, b_param, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(100000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_param) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# CM theory: j=0 curves
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

def find_b_for_n(p, n, max_b=500):
    rng2 = random.Random(7777)
    for b_try in range(1, min(p, max_b)):
        for _ in range(10):
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
# Curve search: bucket by λ/n target
# ---------------------------------------------------------------------------

LAM_TARGETS = [0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34]
LAM_TOL = 0.015
MAX_PER_TARGET = 2

def find_curves_by_lamn(p_lo, p_hi, targets=LAM_TARGETS, tol=LAM_TOL,
                        max_per=MAX_PER_TARGET):
    """
    Scan j=0 CM primes in [p_lo, p_hi] and bucket by λ/n proximity to each target.
    Returns dict: target -> list of (p, b, n, lam, lam_n).
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
                    if n_cand < 4:
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
                            print(f"  [{t:.2f}] p={p}, n={n_cand}, λ={lam}, "
                                  f"λ/n={lam_n:.4f}")
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    print(f"  Scanned {scanned} primes in [{p_lo}, {p_hi}]")
    return found

# ---------------------------------------------------------------------------
# Attack: GLV signatures + lattice
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m, K1, K2, seed):
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 500000:
        attempts += 1
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(0, K2 - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0:
            continue
        R = ec_mul(G, k_full, p)
        if R is None:
            continue
        r = R[0] % n
        if r == 0:
            continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0:
            continue
        s_inv = modinv(s, n)
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})
    return d_secret, sigs

def build_lattice(sigs, n, lam, K1, K2):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = max(1, n // K1)
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
    M[2*m+1][dim - 1] = S_KANNAN
    return M, S_KANNAN

def try_recover(M_red, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

def run_attack(p, b_param, n, lam, G, m, K1, K2, seed,
               use_bkz=False, bkz_beta=20):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, K2, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, K1, K2)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep(p, b_param, n, lam, G, K1, K2, m_range, seeds,
          use_bkz=False, bkz_beta=20):
    results = {}
    for m in m_range:
        wins = sum(run_attack(p, b_param, n, lam, G, m, K1, K2, s,
                              use_bkz=use_bkz, bkz_beta=bkz_beta)
                   for s in seeds)
        results[m] = (wins, len(seeds))
    return results

def first_full_success(results):
    for m, (w, t) in sorted(results.items()):
        if w == t:
            return m
    return None

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: λ/n threshold bisection")
print("Targets: λ/n in {0.07, 0.10, ..., 0.34} — where does LLL cross fail→success?")
print("=" * 70)

SEEDS = [42, 1234, 9999]
M_RANGE = range(3, 15)

# Use 13-bit curves (p ~ 4096-8192) for speed; fallback to 16-bit if sparse
P_LO, P_HI = 2**12, 2**14
print(f"\nScanning {P_LO}–{P_HI} (13-bit) j=0 CM primes:")
bucket = find_curves_by_lamn(P_LO, P_HI)

# Check coverage; extend to 16-bit for targets with no curves
missing = [t for t in LAM_TARGETS if not bucket[t]]
if missing:
    print(f"\nMissing targets {missing}, extending scan to 16-bit...")
    P_LO2, P_HI2 = 2**14, 2**16
    bucket2 = find_curves_by_lamn(P_LO2, P_HI2, targets=missing)
    for t in missing:
        bucket[t].extend(bucket2[t])

# K1 calibration: target eff = K1*K2/n ≈ 0.12 for each curve
# With K2 = √n, eff = K1/√n → K1 = 0.12 * √n
def get_k1_k2(n):
    K2 = math.isqrt(n) + 1
    K1 = max(4, int(0.12 * math.sqrt(n)))
    return K1, K2

print("\n" + "=" * 70)
print("LLL sweep: m=3..14, 3 seeds, K1 ~ 0.12·√n (eff ≈ 0.12)")
print("=" * 70)

lll_results = {}   # (lam_n, p, n) -> (first_success_m_or_None, detail_dict)

for t in LAM_TARGETS:
    curves = bucket[t]
    print(f"\n[λ/n ≈ {t:.2f}]  ({len(curves)} curve(s))")
    if not curves:
        print("  NO CURVES FOUND — gap in search; skip.")
        lll_results[(t, 0, 0)] = (None, {})
        continue
    for (p, b_param, n, lam, lam_n) in curves:
        K1, K2 = get_k1_k2(n)
        eff = K1 * K2 / n
        print(f"  p={p}, n={n} ({n.bit_length()}b), λ={lam}, λ/n={lam_n:.4f}, "
              f"K1={K1}, K2={K2}, eff={eff:.3f}")
        G = find_generator(p, b_param, n)
        if G is None:
            print("  ERROR: generator not found; skip")
            continue
        res = sweep(p, b_param, n, lam, G, K1, K2, M_RANGE, SEEDS)
        fw = first_full_success(res)
        detail = "  ".join(f"m{m}:{w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if fw:
            print(f"  → 3/3 at m={fw}  [{detail}]")
        else:
            max_w = max(w for w, t2 in res.values())
            print(f"  → FAIL (max {max_w}/3)  [{detail}]")
        lll_results[(lam_n, p, n)] = (fw, res)

# BKZ rescue: test targets where LLL failed
print("\n" + "=" * 70)
print("BKZ(beta=20) rescue: re-test targets where LLL failed")
print("=" * 70)

bkz_results = {}

for t in LAM_TARGETS:
    curves = bucket[t]
    for entry in curves:
        p, b_param, n, lam, lam_n = entry
        if (lam_n, p, n) not in lll_results:
            continue
        fw_lll, _ = lll_results[(lam_n, p, n)]
        if fw_lll is not None:
            bkz_results[(lam_n, p, n)] = (fw_lll, "LLL OK — skip BKZ")
            continue
        K1, K2 = get_k1_k2(n)
        print(f"\n  BKZ rescue: λ/n={lam_n:.4f}, p={p}, n={n}")
        G = find_generator(p, b_param, n)
        if G is None:
            bkz_results[(lam_n, p, n)] = (None, "ERROR: no generator")
            continue
        res = sweep(p, b_param, n, lam, G, K1, K2, M_RANGE, SEEDS,
                    use_bkz=True, bkz_beta=20)
        fw_bkz = first_full_success(res)
        detail = "  ".join(f"m{m}:{w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if fw_bkz:
            print(f"  → BKZ 3/3 at m={fw_bkz}  [{detail}]")
        else:
            max_w = max(w for w, t2 in res.values())
            print(f"  → BKZ FAIL (max {max_w}/3)  [{detail}]")
        bkz_results[(lam_n, p, n)] = (fw_bkz, res)

# Summary table
print("\n" + "=" * 70)
print("SUMMARY TABLE: λ/n vs LLL outcome")
print(f"{'λ/n':<8} {'curve':<22} {'K1':>4} {'K2':>5} {'eff':>6}  {'LLL':>10}  {'BKZ(20)':>10}")
print("-" * 70)

for t in LAM_TARGETS:
    curves = bucket[t]
    if not curves:
        print(f"  {t:.2f}   (no curves found)")
        continue
    for (p, b_param, n, lam, lam_n) in curves:
        K1, K2 = get_k1_k2(n)
        eff = K1 * K2 / n
        key = (lam_n, p, n)
        fw_lll, _ = lll_results.get(key, (None, {}))
        lll_str = f"m={fw_lll}" if fw_lll else "FAIL"
        bkz_entry = bkz_results.get(key, (None, ""))
        bkz_str = f"m={bkz_entry[0]}" if bkz_entry[0] is not None else \
                  ("skip" if fw_lll is not None else "FAIL")
        label = f"p={p},n={n}"
        print(f"  {lam_n:.4f}  {label:<22} {K1:>4} {K2:>5} {eff:>6.3f}  {lll_str:>10}  {bkz_str:>10}")

# Locate transition
print("\n" + "=" * 70)
print("THRESHOLD ESTIMATE:")
last_fail = None
first_succ = None
for t in LAM_TARGETS:
    curves = bucket[t]
    for (p, b_param, n, lam, lam_n) in curves:
        key = (lam_n, p, n)
        fw, _ = lll_results.get(key, (None, {}))
        if fw is None:
            last_fail = lam_n
        else:
            if first_succ is None or lam_n < first_succ:
                first_succ = lam_n

if last_fail is not None and first_succ is not None:
    print(f"  LLL last failure:  λ/n = {last_fail:.4f}")
    print(f"  LLL first success: λ/n = {first_succ:.4f}")
    mid = (last_fail + first_succ) / 2
    print(f"  Threshold estimate: τ* ≈ {mid:.4f}")
elif first_succ is not None:
    print(f"  LLL succeeds from λ/n = {first_succ:.4f} (no failure detected in range)")
elif last_fail is not None:
    print(f"  LLL fails up to λ/n = {last_fail:.4f} (no success detected in range)")
else:
    print("  Insufficient data.")

print("\nDone.")
