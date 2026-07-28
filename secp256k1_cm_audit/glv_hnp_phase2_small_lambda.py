"""
GLV-HNP Phase 2: Small-λ threshold study (Thread 20).

Goal: find the minimum λ/n at which the Phase 2 GLV-aware HNP attack succeeds,
using 20-bit j=0 CM curves at a fixed eff≈0.05.

Context (from 2026-07-26 run):
  - 8-bit  (n=199,    λ/n=0.53): LLL 3/3 at m=6        [SUCCESS]
  - 12-bit (n=2659,   λ/n=0.66): LLL 3/3 at m=7        [SUCCESS]
  - 20-bit (n=523969, λ/n=0.34): LLL 3/3 at m=9        [SUCCESS]
  - 12-bit (n=2647,   λ/n=0.07): LLL+BKZ(40) NEVER 3/3 [FAIL]

The "small-λ failure" at λ/n=0.07 may be a 12-bit artefact or a genuine
structural obstruction. Thread 20 fixes the curve size at 20-bit and varies
λ/n ∈ {0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.26} to find the threshold.

Phase 2 attack parameters: K1=36 (≈0.05·√n for n≈500K), K2=√n+1, eff≈0.05.

Run: python3 glv_hnp_phase2_small_lambda.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# EC arithmetic
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

LAM_TARGETS = [0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.26]
LAM_TOL = 0.015
MAX_PER_TARGET = 2

def find_curves_by_lamn(p_lo=2**19, p_hi=2**20, targets=LAM_TARGETS,
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
                            print(f"  [λ/n≈{t:.2f}] p={p}, n={n_cand}, lam={lam}, "
                                  f"λ/n={lam_n:.4f}")
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    print(f"  (scanned {scanned} primes in [{p_lo}, {p_hi}])")
    return found

# ---------------------------------------------------------------------------
# Phase 2 lattice attack
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

def run_lll(p, b_param, n, lam, G, m, K1, seed, use_bkz=False, bkz_beta=30):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, K1)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep(p, b_param, n, lam, G, K1, m_range, seeds, use_bkz=False, bkz_beta=30):
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_lll(p, b_param, n, lam, G, m, K1, s,
                                                use_bkz=use_bkz, bkz_beta=bkz_beta))
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: Small-λ threshold study (Thread 20)")
print("Fixed eff≈0.05 (K1=36 for 20-bit), λ/n ∈ {0.05..0.26}")
print("=" * 70)
print()

K1_FIXED = 36   # ≈ 0.05 * sqrt(500K) — fixed for all 20-bit curves
SEEDS    = [42, 1234, 9999]
M_RANGE  = range(5, 21)   # m=5..20

# ---- Part 0: Reference point from prior work (20-bit, λ/n=0.34) --------
# This was the first 20-bit success from 2026-07-26; include as anchor.
REF_20 = dict(p=524347, b=2, n=523969, lam=177902, lam_n=0.34)
print("Reference anchor: 20-bit λ/n=0.34 (known to succeed at m=9)")
print(f"  p={REF_20['p']}, n={REF_20['n']}, lam={REF_20['lam']}, K1={K1_FIXED}")
G_ref = find_generator(REF_20['p'], REF_20['b'], REF_20['n'])
if G_ref:
    res_ref = sweep(REF_20['p'], REF_20['b'], REF_20['n'], REF_20['lam'], G_ref,
                    K1_FIXED, range(5, 14), SEEDS)
    fw_ref = next((m for m, (w, t) in sorted(res_ref.items()) if w == t), None)
    print(f"  LLL sweep m=5..13: first 3/3 at m={fw_ref}")
else:
    print("  ERROR: generator not found for reference curve")
print()

# ---- Part 1: Reference failure (12-bit, λ/n=0.07) ----------------------
# Reproduced from 2026-07-26: K1=8 (≈0.05*sqrt(2647)≈2.6→use 8 as in prior run).
print("Reference failure: 12-bit λ/n=0.07 (known to fail)")
p_f, b_f, n_f, lam_f = 2677, 2, 2647, 185
K1_12 = max(2, round(0.05 * math.sqrt(n_f)))   # K1=2 at eff=0.05
print(f"  p={p_f}, n={n_f}, lam={lam_f}, λ/n={lam_f/n_f:.4f}, K1={K1_12}")
G_f = find_generator(p_f, b_f, n_f)
if G_f:
    res_f = sweep(p_f, b_f, n_f, lam_f, G_f, K1_12, range(5, 16), SEEDS)
    fw_f = next((m for m, (w, t) in sorted(res_f.items()) if w == t), None)
    detail_f = " ".join(f"m{m}={w}/{t}" for m, (w, t) in sorted(res_f.items()))
    print(f"  K1={K1_12} (eff≈{K1_12*(math.isqrt(n_f)+1)/n_f:.3f}): {detail_f}")
    print(f"  → first 3/3: {fw_f}")
else:
    print("  ERROR: generator not found")
print()

# ---- Part 2: Scan 20-bit curves at small λ/n targets -------------------
print("=" * 70)
print("Scanning 20-bit j=0 CM curves at λ/n targets (below Effect A zone)")
print("-" * 40)
bucket = find_curves_by_lamn()
print()

# ---- Part 3: LLL sweep for each target ---------------------------------
print("=" * 70)
print(f"LLL sweep: K1={K1_FIXED} (eff≈{K1_FIXED*724/524e3:.4f}), m=5..20, 3 seeds")
print("-" * 40)

sweep_results = {}   # (lam_n, p, n) -> (first_win_m or None, detail_str)

for t in LAM_TARGETS:
    curves = bucket[t]
    print(f"\n[λ/n≈{t:.2f}]  ({len(curves)} curves found)")
    if not curves:
        print("  No curves found; skipping.")
        continue
    for (p, b, n, lam, lam_n) in curves:
        K2 = math.isqrt(n) + 1
        eff = K1_FIXED * K2 / n
        print(f"  p={p}, n={n}, lam={lam}, λ/n={lam_n:.4f}, eff={eff:.4f}")
        G = find_generator(p, b, n)
        if G is None:
            print("  ERROR: generator not found; skip")
            continue
        res = sweep(p, b, n, lam, G, K1_FIXED, M_RANGE, SEEDS)
        fw = next((m for m, (w, t2) in sorted(res.items()) if w == t2), None)
        detail = " ".join(f"m{m}={w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if fw:
            print(f"  → LLL 3/3 first at m={fw}   [{detail}]")
        else:
            max_wins = max(w for w, t2 in res.values())
            print(f"  → FAIL never 3/3 (max {max_wins}/3)   [{detail}]")
        sweep_results[(lam_n, p, n)] = (fw, detail)

# ---- Part 4: BKZ rescue on failures ------------------------------------
print()
print("=" * 70)
print(f"BKZ(β=30) rescue pass on LLL failures")
print("-" * 40)

for t in LAM_TARGETS:
    curves = bucket[t]
    for (p, b, n, lam, lam_n) in curves:
        key = (lam_n, p, n)
        if key not in sweep_results:
            continue
        fw, _ = sweep_results[key]
        if fw is not None:
            continue   # already succeeded with LLL
        print(f"\n  Retrying p={p}, n={n}, lam={lam}, λ/n={lam_n:.4f} with BKZ(30)")
        G = find_generator(p, b, n)
        if G is None:
            print("  ERROR: generator not found; skip")
            continue
        res_bkz = sweep(p, b, n, lam, G, K1_FIXED, M_RANGE, SEEDS,
                        use_bkz=True, bkz_beta=30)
        fw_bkz = next((m for m, (w, t2) in sorted(res_bkz.items()) if w == t2), None)
        detail_bkz = " ".join(f"m{m}={w}/{t2}" for m, (w, t2) in sorted(res_bkz.items()))
        if fw_bkz:
            print(f"  → BKZ(30) 3/3 first at m={fw_bkz}   [{detail_bkz}]")
        else:
            max_wins = max(w for w, t2 in res_bkz.values())
            print(f"  → BKZ(30) also fails (max {max_wins}/3)   [{detail_bkz}]")
        sweep_results[key] = (fw_bkz, f"BKZ30:{detail_bkz}")

# ---- Part 5: Summary table ---------------------------------------------
print()
print("=" * 70)
print("SUMMARY: Phase 2 attack (K1=36, eff≈0.05) vs λ/n")
print(f"{'λ/n':<8} {'n':>10} {'lam':>10} {'first 3/3 m':>13} {'outcome'}")
print("-" * 65)

# Reference anchor
if G_ref:
    fw_r = next((m for m, (w, t) in sorted(res_ref.items()) if w == t), None)
    print(f"  0.340   {REF_20['n']:>10}  {REF_20['lam']:>10}  {str(fw_r) if fw_r else 'never':>13}  "
          f"{'SUCCESS' if fw_r else 'FAIL'}  [20-bit anchor]")

# Main sweep
for t in LAM_TARGETS:
    for (p, b, n, lam, lam_n) in bucket[t]:
        key = (lam_n, p, n)
        if key not in sweep_results:
            continue
        fw, _ = sweep_results[key]
        outcome = f"SUCCESS m={fw}" if fw else "FAIL"
        print(f"  {lam_n:.4f}  {n:>10}  {lam:>10}  {str(fw) if fw else 'never':>13}  {outcome}")

# 12-bit reference failure
if G_f:
    fw_f2 = next((m for m, (w, t) in sorted(res_f.items()) if w == t), None)
    print(f"  {lam_f/n_f:.4f}  {n_f:>10}  {lam_f:>10}  {str(fw_f2) if fw_f2 else 'never':>13}  "
          f"{'SUCCESS' if fw_f2 else 'FAIL (12-bit ref)'}  [12-bit K1={K1_12}]")

print()
print("Done.")
