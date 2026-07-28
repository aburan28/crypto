"""
Thread 20: λ/n threshold bisection for GLV-HNP Phase 2.

Sweeps over abstract GLV-HNP instances with different λ/n ratios to
find the empirical threshold below which the Phase 2 LLL attack fails.

Known endpoints (from glv_hnp_phase2_20bit.py, 2026-07-26 log):
  λ/n ≈ 0.07 (n=2647, λ=185, K1=8):  FAILS (LLL + BKZ-40) for all m≤12
  λ/n ≈ 0.53 (n=199,  λ=106, K1=2):  SUCCEEDS at m=5

KEY PARAMETER: K1_BOUND is proportional to n (matching prior experiments):
  K1_BOUND = max(2, round(K1_FRAC * n))  where K1_FRAC = 0.003 (≈ 0.3% of n)

This script:
  1. Enumerates primes n ≡ 1 (mod 3) in [100, 50000]
  2. Computes λ_min = smaller root of x²+x+1 ≡ 0 (mod n)
  3. Groups by λ/n ratio bins
  4. Uses EC-based signature generation (REAL GLV sigs from j=0 curves)
  5. Reports LLL success/failure table → threshold identified

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Bias fraction (K1 = max(2, round(K1_FRAC * n))), matching 20bit script.
# For n=199: K1=2.  For n=2647: K1=8.  Scales as O(n).
# ---------------------------------------------------------------------------
K1_FRAC = 0.003   # 0.3% of n

SEEDS = [42, 1234, 9999]

# ---------------------------------------------------------------------------
# Primality
# ---------------------------------------------------------------------------

def is_prime(n):
    if n < 2: return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
        if n == p: return True
        if n % p == 0: return False
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1: break
        else:
            return False
    return True

# ---------------------------------------------------------------------------
# Tonelli-Shanks
# ---------------------------------------------------------------------------

def sqrt_mod(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, rr = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while True:
        if t == 1: return rr
        i, tmp = 1, t * t % p
        while tmp != 1:
            tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, rr = i, b * b % p, t * b * b % p, rr * b % p

# ---------------------------------------------------------------------------
# GLV lambda: smaller root of x^2+x+1 ≡ 0 (mod n)
# ---------------------------------------------------------------------------

def find_glv_lambda(n):
    """Return (lam_small, lam_large) or (None, None) if n ≢ 1 (mod 3)."""
    if n % 3 != 1: return None, None
    sq = sqrt_mod(n - 3, n)  # sqrt(-3) mod n
    if sq is None: return None, None
    inv2 = pow(2, -1, n)
    r1 = (-1 + sq) * inv2 % n
    r2 = (-1 - sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0
    assert (r2 * r2 + r2 + 1) % n == 0
    lam = min(r1, r2)
    return lam, max(r1, r2)

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass a=0, y^2 = x^3 + b)
# ---------------------------------------------------------------------------

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        s = (y2 - y1) * pow(x2 - x1, -1, p) % p
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

def find_generator(p, b, n, max_tries=5000):
    """Find a point of order n on y^2 = x^3 + b over F_p."""
    rng = random.Random(12345)
    for _ in range(max_tries):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = sqrt_mod(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# Eisenstein CM decomposition (find j=0 curve with given trace)
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """Find (a,b) with a^2 - ab + b^2 = p."""
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in [a + s, a - s]:
            if num % 2 == 0:
                bv = num // 2
                if bv >= 0 and a * a - a * bv + bv * bv == p:
                    return (a, bv)
    return None

def j0_traces(a, b):
    """6 possible traces for j=0 sextic twists given Eisenstein decomp (a,b)."""
    return [2*a - b, -(2*a - b), a + b, -(a + b), 2*b - a, -(2*b - a)]

def find_j0_curve(p_min, p_max, lam_ratio_target, lam_ratio_tol=0.04,
                  require_prime_n=True):
    """Find a j=0 curve over F_p (p_min ≤ p ≤ p_max) with λ/n near target.

    Returns (p, b, n, lam) or None.
    """
    for p in range(p_min | 1, p_max + 1, 2):
        if not is_prime(p): continue
        if p % 3 != 1: continue
        eis = eisenstein_decompose(p)
        if eis is None: continue
        a_e, b_e = eis
        for t in j0_traces(a_e, b_e):
            n_cand = p + 1 - t
            if n_cand < 7: continue
            if require_prime_n and not is_prime(n_cand): continue
            if n_cand % 3 != 1: continue
            ls, ll = find_glv_lambda(n_cand)
            if ls is None: continue
            # Test BOTH roots (attack might use either)
            for lam in [ls, ll]:
                ratio = lam / n_cand
                if abs(ratio - lam_ratio_target) <= lam_ratio_tol:
                    # Find curve parameter b
                    for b_try in range(1, min(p, 200)):
                        G = find_generator(p, b_try, n_cand, max_tries=200)
                        if G is not None:
                            return (p, b_try, n_cand, lam)
    return None

# ---------------------------------------------------------------------------
# Signature generation (REAL GLV sigs from EC arithmetic)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 50 * m:
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
        s = pow(k_full, -1, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = pow(s, -1, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

# ---------------------------------------------------------------------------
# Lattice construction
# ---------------------------------------------------------------------------

def build_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_D = 1
    S_KAN = n
    dim = 2 * m + 2
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
    M[2 * m + 1][dim - 1] = S_KAN
    return M, S_KAN

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def try_recover_d(M_lll, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in M_lll:
        if abs(row[dim - 1]) != S_KAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Single experiment
# ---------------------------------------------------------------------------

def run_experiment(p, b, n, lam, G, d_secret, m, k1_bound, k2_bound,
                   seed=42, use_bkz=False, bkz_beta=20):
    """Run Phase 2 attack. Returns True iff d recovered."""
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M, S_KAN = build_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(2 * m + 2)] for i in range(2 * m + 2)]
    return try_recover_d(reduced, m, n, S_KAN, d_secret)

# ---------------------------------------------------------------------------
# Sweep one ratio bin: pick representative curve, run LLL at multiple m
# ---------------------------------------------------------------------------

def sweep_ratio(ratio_target, ratio_tol, n_range, label=""):
    """Find curve with given λ/n and run LLL/BKZ sweep. Returns dict."""
    p_lo, p_hi = n_range
    result = find_j0_curve(p_lo, p_hi, ratio_target, ratio_tol)
    if result is None:
        return None

    p, b, n, lam = result
    k2_bound = math.isqrt(n) + 1
    k1_bound = max(2, round(K1_FRAC * n))
    eff = k1_bound * k2_bound / n
    if eff >= 1.0:
        m_thresh = 2
    else:
        m_thresh = max(2, math.ceil(math.log(n) / math.log(1.0 / eff)))

    G = find_generator(p, b, n)
    if G is None:
        return None

    m_range = range(m_thresh, m_thresh + 7)
    lll_wins = {}
    for m in m_range:
        wins = 0
        for seed in SEEDS:
            d = random.Random(seed * 7 + 13).randint(1, n - 1)
            if run_experiment(p, b, n, lam, G, d, m, k1_bound, k2_bound, seed):
                wins += 1
        lll_wins[m] = wins
        if wins == len(SEEDS):
            break  # found m that works; don't need more

    # Report the result at m_thresh and at max m tried
    m_success = next((m for m, w in sorted(lll_wins.items()) if w == len(SEEDS)), None)
    max_m_tried = max(lll_wins)
    wins_at_max = lll_wins[max_m_tried]

    return {
        'p': p, 'b': b, 'n': n, 'lam': lam,
        'ratio': lam / n,
        'k1': k1_bound, 'k2': k2_bound,
        'm_thresh': m_thresh, 'eff': eff,
        'm_success': m_success,        # None = never 3/3 in sweep
        'max_m_tried': max_m_tried,
        'wins_at_max': wins_at_max,
        'wins_table': lll_wins,
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("Thread 20: GLV-HNP Phase 2 — λ/n threshold bisection")
    print(f"K1_FRAC={K1_FRAC:.3f} (proportional to n), seeds={SEEDS}")
    print("=" * 72)
    print()

    # ---- Step 1: Reproduce known endpoints to validate framework ----
    print("STEP 1: Reproduce known endpoints")
    print("-" * 72)

    # Known FAILURE: p=2677, n=2647, lam=185 (lam/n≈0.07), K1=8
    p_f, b_f, n_f, lam_f = 2677, 2, 2647, 185
    k1_f = max(2, round(K1_FRAC * n_f))  # should be 8
    k2_f = math.isqrt(n_f) + 1
    G_f = find_generator(p_f, b_f, n_f)
    if G_f:
        print(f"Failure case: p={p_f}, n={n_f}, λ={lam_f}, λ/n={lam_f/n_f:.4f}, K1={k1_f}, K2={k2_f}")
        for m in [6, 8, 10, 12]:
            wins = sum(
                run_experiment(p_f, b_f, n_f, lam_f, G_f,
                               random.Random(s * 7 + 13).randint(1, n_f - 1),
                               m, k1_f, k2_f, seed=s)
                for s in SEEDS
            )
            print(f"  m={m}: {wins}/3 seeds  {'(FAIL expected)' if wins < 3 else '(UNEXPECTED PASS)'}")

    print()

    # Known SUCCESS: p=211, n=199, lam=106 (lam/n≈0.53), K1=2
    p_s, b_s, n_s, lam_s = 211, 2, 199, 106
    k1_s = max(2, round(K1_FRAC * n_s))  # should be 2 (or 1 → max 2)
    k2_s = math.isqrt(n_s) + 1
    G_s = find_generator(p_s, b_s, n_s)
    if G_s:
        print(f"Success case: p={p_s}, n={n_s}, λ={lam_s}, λ/n={lam_s/n_s:.4f}, K1={k1_s}, K2={k2_s}")
        for m in [4, 5, 6]:
            wins = sum(
                run_experiment(p_s, b_s, n_s, lam_s, G_s,
                               random.Random(s * 7 + 13).randint(1, n_s - 1),
                               m, k1_s, k2_s, seed=s)
                for s in SEEDS
            )
            print(f"  m={m}: {wins}/3 seeds  {'(SUCCESS expected)' if wins == 3 else ''}")

    print()

    # ---- Step 2: Sweep ratio bins [0.05, 0.55] step=0.05 ----
    print("=" * 72)
    print("STEP 2: Coarse sweep λ/n ∈ [0.05, 0.55], step=0.05, K1∝n")
    print("=" * 72)
    print(f"{'λ/n target':>12}  {'actual λ/n':>10}  {'n':>7}  {'λ':>7}  "
          f"{'K1':>4}  {'K2':>4}  {'m*':>4}  {'max':>4}  result")
    print("-" * 80)

    coarse_results = {}
    for ratio_target_10 in range(5, 60, 5):
        ratio_target = ratio_target_10 / 100.0
        # Search in n ∈ [200, 30000]
        res = sweep_ratio(ratio_target, 0.03, (200, 30000))
        if res is None:
            print(f"  [{ratio_target:.2f}]  NO CURVE FOUND")
            continue
        m_star = res['m_success']
        w_max = res['wins_at_max']
        label = ("✓ SUCCEED" if m_star is not None else
                 ("△ PARTIAL" if w_max > 0 else "✗ FAIL"))
        m_str = f"m={m_star}" if m_star is not None else f"---"
        print(f"  [{ratio_target:.2f}]  {res['ratio']:.4f}     {res['n']:>7}  {res['lam']:>7}  "
              f"{res['k1']:>4}  {res['k2']:>4}  {m_str:>6}  {res['max_m_tried']:>4}  "
              f"{w_max}/{len(SEEDS)}  {label}")
        coarse_results[ratio_target] = res

    print()

    # ---- Step 3: Fine bisection in the transition region ----
    successes = [r for r, res in coarse_results.items() if res.get('m_success') is not None]
    failures  = [r for r, res in coarse_results.items() if res.get('m_success') is None]

    if successes and failures:
        max_fail  = max(failures)
        min_succ  = min(successes)
        print(f"Coarse transition: fails at λ/n ≤ {max_fail:.2f}, succeeds at λ/n ≥ {min_succ:.2f}")
        print()
        print("=" * 72)
        print(f"STEP 3: Fine bisection λ/n ∈ [{max_fail:.2f}, {min_succ:.2f}], step=0.01")
        print("=" * 72)
        print(f"{'λ/n target':>12}  {'actual λ/n':>10}  {'n':>7}  {'K1':>4}  "
              f"{'m_thresh':>8}  result")
        print("-" * 65)

        for ratio_target_100 in range(int(max_fail * 100), int(min_succ * 100) + 2, 1):
            ratio_target = ratio_target_100 / 100.0
            res = sweep_ratio(ratio_target, 0.015, (500, 50000))
            if res is None:
                print(f"  [{ratio_target:.2f}]  NO CURVE")
                continue
            m_star = res['m_success']
            w_max = res['wins_at_max']
            label = ("✓ SUCCEED" if m_star is not None else
                     ("△ PARTIAL" if w_max > 0 else "✗ FAIL"))
            m_str = f"m={m_star}" if m_star is not None else f"---"
            print(f"  [{ratio_target:.2f}]  {res['ratio']:.4f}     {res['n']:>7}  {res['k1']:>4}  "
                  f"m≥{res['m_thresh']}    {w_max}/{len(SEEDS)} @ max_m={res['max_m_tried']}  "
                  f"{m_str}  {label}")
            coarse_results[ratio_target] = res

        # Final threshold estimate
        all_fails  = [r for r, res in coarse_results.items() if res.get('m_success') is None]
        all_succs  = [r for r, res in coarse_results.items() if res.get('m_success') is not None]
        if all_fails and all_succs:
            thr_lo = max(all_fails)
            thr_hi = min(all_succs)
            print()
            print(f"THRESHOLD ESTIMATE: λ/n ∈ ({thr_lo:.2f}, {thr_hi:.2f})")
            print(f"  → LLL Phase-2 attack fails when λ/n < {thr_lo:.2f}")
            print(f"  → LLL Phase-2 attack succeeds when λ/n > {thr_hi:.2f}")
            print()
            print("SECP256K1 IMPLICATION:")
            lam_sec = 37718080363155996902926221483475592026636932068901688066867496373061895675762
            n_sec   = 115792089237316195423570985008687907852837564279074904382605163141518161494337
            sec_ratio = lam_sec / n_sec
            print(f"  secp256k1 λ/n ≈ {sec_ratio:.4f}")
            if sec_ratio < thr_lo:
                print(f"  BELOW threshold: Phase 2 attack would FAIL for secp256k1.")
            elif sec_ratio > thr_hi:
                print(f"  ABOVE threshold: Phase 2 attack would SUCCEED for secp256k1 (modulo precision).")
            else:
                print(f"  IN transition zone: inconclusive.")
    else:
        print("Transition region not found in coarse sweep.")

    print()
    print("Done.")

if __name__ == "__main__":
    main()
