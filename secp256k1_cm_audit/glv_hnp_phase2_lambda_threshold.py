"""
GLV-HNP Phase 2: small-λ threshold bisection (Thread 20, 2026-07-27).

From the 2026-07-26 Phase 2 run:
  - λ/n ≥ 0.34 (8–20 bit curves): LLL succeeds (3/3 seeds).
  - λ/n = 0.07 (p=2677, 12-bit): LLL AND BKZ(40) both fail (structural).

Hypothesis: there is a sharp threshold T* such that
  - λ/n < T*  →  LLL always fails (small-λ structural obstruction)
  - λ/n > T*  →  LLL succeeds at some finite m

This script bisects T* by finding 12-bit j=0 CM curves at λ/n targets
{0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34} and running
the Phase 2 attack with m=4..16, K1=8, seeds×3.

Note: glv_hnp_lamn_boundary.py covers the 1/3 obstruction (Effect A) at 20-bit.
This script covers the SMALL-λ regime at 12-bit — a distinct phenomenon.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

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
    mv, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mv - i - 1), p)
        mv, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b_param, n, rng_seed=12345):
    rng = random.Random(rng_seed)
    for _ in range(100000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b_param) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def find_b(p, n):
    for b in range(1, min(p, 200)):
        rng = random.Random(99)
        for _ in range(20):
            x = rng.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                if ec_mul((x, y), n, p) is None:
                    return b
                break
    return None

# ---------------------------------------------------------------------------
# CM theory for j=0 curves
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

# ---------------------------------------------------------------------------
# Curve search: 12-bit j=0 CM curves, bucketed by λ/n target
# ---------------------------------------------------------------------------

LAM_TARGETS = [0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34]
LAM_TOL = 0.015

def find_12bit_curves(targets=LAM_TARGETS, tol=LAM_TOL, max_per=3,
                      p_lo=2**11, p_hi=2**13):
    """Find 12-bit j=0 CM curves bucketed by λ/n proximity."""
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
                    if n_cand < 3:
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
                            b_param = find_b(p, n_cand)
                            if b_param is None:
                                continue
                            found[t].append((p, b_param, n_cand, lam, lam_n))
                            print(f"  [λ/n≈{t:.2f}] p={p}, n={n_cand}, "
                                  f"λ={lam}, λ/n={lam_n:.4f}")
                            break
        p = sympy.nextprime(p)
        if all(len(found[t]) >= max_per for t in targets):
            break
    print(f"  (scanned {scanned} primes in [{p_lo}, {p_hi}])")
    return found

# ---------------------------------------------------------------------------
# Phase 2 lattice attack
# ---------------------------------------------------------------------------

K1_BOUND = 8     # small k1 bias: k1 in [0, K1_BOUND)

def gen_signatures(p, b_param, n, lam, G, m, k1_bound, seed):
    k2_bound = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
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
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})
    return d_secret, sigs

def build_lattice(sigs, n, lam, k1_bound):
    m = len(sigs)
    k2_bound = math.isqrt(n) + 1
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
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
        if abs(row[dim-1]) != S_KANNAN: continue
        sign = 1 if row[dim-1] > 0 else -1
        if (sign * row[m]) % n == d_secret:
            return True
    return False

def run_lll(p, b_param, n, lam, G, m, k1_bound, seed):
    d_secret, sigs = gen_signatures(p, b_param, n, lam, G, m, k1_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KANNAN = build_lattice(sigs, n, lam, k1_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return try_recover(reduced, m, n, S_KANNAN, d_secret)

def sweep_curve(p, b_param, n, lam, G, k1_bound, m_range, seeds):
    results = {}
    for m in m_range:
        wins = sum(1 for s in seeds if run_lll(p, b_param, n, lam, G, m, k1_bound, s))
        results[m] = (wins, len(seeds))
    return results

# ---------------------------------------------------------------------------
# Reference: reproduce p=2677 failure
# ---------------------------------------------------------------------------

def probe_p2677():
    """Verify the p=2677 failure from 2026-07-26 log."""
    print("=" * 70)
    print("Anchor: p=2677 failure (λ/n=0.07) from 2026-07-26 log")
    print("=" * 70)
    p = 2677
    # Find n by Eisenstein decompose
    eis = eisenstein_decompose(p)
    if eis is None:
        print("  ERROR: p=2677 not an Eisenstein prime")
        return
    a_e, b_e = eis
    traces = j0_traces(a_e, b_e)
    n_found = None
    for tr in traces:
        n_cand = p + 1 - tr
        if n_cand >= 2 and sympy.isprime(n_cand) and n_cand % 3 == 1:
            lam = glv_eigenvalue(n_cand)
            if lam is not None:
                lam_n = lam / n_cand
                if abs(lam_n - 0.07) < 0.03:
                    n_found = n_cand
                    lam_found = lam
                    lam_n_found = lam_n
                    break
    if n_found is None:
        # Fall back: pick the n with smallest λ/n
        candidates = []
        for tr in traces:
            n_cand = p + 1 - tr
            if n_cand >= 2 and sympy.isprime(n_cand) and n_cand % 3 == 1:
                lam = glv_eigenvalue(n_cand)
                if lam is not None:
                    candidates.append((lam / n_cand, n_cand, lam))
        if not candidates:
            print("  ERROR: no valid j=0 prime-order curve found for p=2677")
            return
        candidates.sort()
        lam_n_found, n_found, lam_found = candidates[0]
    print(f"  p={p}, n={n_found}, λ={lam_found}, λ/n={lam_n_found:.4f}")
    b_param = find_b(p, n_found)
    if b_param is None:
        print("  ERROR: no b found for p=2677")
        return
    G = find_generator(p, b_param, n_found)
    if G is None:
        print("  ERROR: no generator for p=2677")
        return
    SEEDS = [42, 1234, 9999]
    M_RANGE = range(4, 17)
    res = sweep_curve(p, b_param, n_found, lam_found, G, K1_BOUND, M_RANGE, SEEDS)
    first_win = next((m for m, (w, t) in sorted(res.items()) if w == t), None)
    for m, (w, t) in sorted(res.items()):
        marker = " ← first 3/3" if m == first_win else ""
        print(f"    m={m}: {w}/{t}{marker}")
    if first_win:
        print(f"  Anchor: 3/3 first at m={first_win}  [UNEXPECTED — log said failure]")
    else:
        max_w = max(w for w, t in res.values())
        print(f"  Anchor: FAIL (max {max_w}/3). Reproduces 2026-07-26 failure. ✓")
    return res

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print("=" * 70)
print("GLV-HNP Phase 2: small-λ threshold bisection (Thread 20)")
print(f"K1_BOUND={K1_BOUND}, m=4..16, seeds×3, 12-bit curves")
print("=" * 70)
print()

# Step 1: confirm p=2677 failure
p2677_res = probe_p2677()
print()

# Step 2: find 12-bit curves at each λ/n target
print("=" * 70)
print("Finding 12-bit j=0 CM curves by λ/n target")
print("-" * 40)
bucket = find_12bit_curves()
print()

# Step 3: run LLL sweep on each bucket
SEEDS = [42, 1234, 9999]
M_RANGE = range(4, 17)

print("=" * 70)
print(f"LLL sweep: K1={K1_BOUND}, m=4..16, seeds×3")
print("-" * 40)

sweep_results = {}   # lam_n -> list of (first_win, n)

for t in LAM_TARGETS:
    curves = bucket[t]
    print(f"\n[λ/n≈{t:.2f}]  ({len(curves)} curves found)")
    if not curves:
        print("  No curves found; skipping.")
        sweep_results[t] = []
        continue
    for (p, b_param, n, lam, lam_n) in curves:
        k2_bound = math.isqrt(n) + 1
        eff = K1_BOUND * k2_bound / n
        G = find_generator(p, b_param, n)
        if G is None:
            print(f"  p={p}: no generator; skip")
            continue
        res = sweep_curve(p, b_param, n, lam, G, K1_BOUND, M_RANGE, SEEDS)
        first_win = next((m for m, (w, t2) in sorted(res.items()) if w == t2), None)
        detail = " ".join(f"m{m}={w}/{t2}" for m, (w, t2) in sorted(res.items()))
        if first_win:
            print(f"  p={p}, n={n}, λ/n={lam_n:.4f}, eff={eff:.4f} → 3/3 at m={first_win}  [{detail}]")
        else:
            max_w = max(w for w, t2 in res.values())
            print(f"  p={p}, n={n}, λ/n={lam_n:.4f}, eff={eff:.4f} → FAIL (max {max_w}/3)  [{detail}]")
        if t not in sweep_results:
            sweep_results[t] = []
        sweep_results[t].append((first_win, lam_n, n))

# Step 4: summary table
print()
print("=" * 70)
print("SUMMARY: λ/n vs LLL outcome (small-λ bisection)")
print(f"{'λ/n':>6}  {'outcome':>20}  first_3of3_m  note")
print("-" * 60)
for t in LAM_TARGETS:
    entries = sweep_results.get(t, [])
    if not entries:
        print(f"  {t:.2f}   (no curves found)")
        continue
    for (fw, lam_n, n) in entries:
        if fw:
            outcome = f"SUCCESS at m={fw}"
        else:
            outcome = "OBSTRUCTED"
        note = "← anchor failure" if abs(lam_n - 0.07) < 0.01 else ""
        print(f"  {lam_n:.4f}  {outcome:<22}  {str(fw) if fw else 'never':<14} {note}")

# Secp256k1 reference
print(f"  0.4400  (reference: secp256k1, 256-bit; not tested at toy scale)")
print()
print("Done.")
