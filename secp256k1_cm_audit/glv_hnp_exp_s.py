"""
GLV-HNP Phase 2: Exp S — test max_a (max partial quotient) as C1-predictor
on 50 fresh CM curves (2026-06-29).

Hypothesis from Exp R data:
  max_a = max{a_k : q_k ≤ 300} in CF expansion of lam/n.
  From 30-curve Exp R sample: C1 max_a ∈ [3,85], C2 max_a ∈ [4,17].
  PROPOSED RULE: max_a ≤ 17 → predict C2 / max_a ≥ 18 → predict C1.

  NOTE: ranges overlap at [4,17], so false positives expected.
  This experiment tests accuracy on a fresh 50-curve sample.

IMPLEMENTATION: Pure Python LLL (no fpylll dependency) sufficient for
26×26 matrices with ~34-bit entries (17-20 bit primes, m=12, K1=72).

Run: python3 glv_hnp_exp_s.py
"""

import math
import random
import time

# ---------------------------------------------------------------------------
# Primality (Miller-Rabin)
# ---------------------------------------------------------------------------

def is_prime(n, rounds=12):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    small = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small:
        if n == p: return True
        if n % p == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0: r += 1; d //= 2
    for a in small[:rounds]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

# ---------------------------------------------------------------------------
# EC arithmetic
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
    return (x3, (s * (x1 - x3) - y1) % p)

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p-1)//2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p+1)//4, p)
    q, s = p-1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    mv, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q+1)//2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mv-i-1), p)
        mv, c, t, r = i, b*b%p, t*b*b%p, r*b%p

def glv_eigenvalue(n):
    sq = tonelli((n-3) % n, n)
    if sq is None: return None
    inv2 = pow(2, -1, n)
    r1 = (n-1+sq) * inv2 % n
    r2 = (n-1+(n-sq)) * inv2 % n
    if (r1*r1+r1+1) % n != 0: r1, r2 = r2, r1
    if (r1*r1+r1+1) % n != 0: return None
    return min(r1, r2)

def find_generator(p, b, n):
    rng = random.Random(9999)
    for _ in range(100000):
        x = rng.randint(0, p-1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli(rhs, p)
        if y and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

# ---------------------------------------------------------------------------
# Cornacchia: a²+3b²=4p for j=0 CM curves
# ---------------------------------------------------------------------------

def cornacchia_j0(p):
    if p % 3 != 1: return None
    x0 = tonelli(p-3, p)
    if x0 is None: return None
    limit = math.isqrt(4*p)
    u, v = 2*p, x0
    while v > limit:
        u, v = v, u % v
    a = v
    b2 = 4*p - a*a
    if b2 <= 0 or b2 % 3 != 0:
        return None
    b = math.isqrt(b2 // 3)
    if b*b*3 + a*a != 4*p: return None
    return abs(a), b

# ---------------------------------------------------------------------------
# CM curve generation (Cornacchia-based)
# ---------------------------------------------------------------------------

KNOWN_N_SET = {523597, 345109, 326479, 337999, 321163, 949423, 904369, 811297, 319129}

def gen_cm_curve(seed, bit_lo=17, bit_hi=20):
    rng = random.Random(seed * 7919 + 31337)
    for _ in range(5000):
        bits = rng.randint(bit_lo, bit_hi)
        p = rng.randint(2**(bits-1), 2**bits - 1) | 1
        p += 0 if p % 3 == 1 else (1 if p % 3 == 0 else 2)
        if p % 2 == 0: p += 1
        if not is_prime(p): continue
        if p % 3 != 1: continue
        corn = cornacchia_j0(p)
        if corn is None: continue
        a_corn, _ = corn
        for sign in [1, -1]:
            n = p + 1 - sign * a_corn
            if n < 4: continue
            if not is_prime(n): continue
            if n in KNOWN_N_SET: continue
            lam = glv_eigenvalue(n)
            if lam is None: continue
            for b_try in range(1, 30):
                y2 = (pow(0, 3, p) + b_try) % p  # x=0
                rhs = b_try % p
                # try a few x values
                found_pt = None
                for x in range(min(50, p)):
                    rhs_x = (pow(x, 3, p) + b_try) % p
                    y = tonelli(rhs_x, p)
                    if y and y != 0:
                        P = (x, y)
                        if ec_mul(P, n, p) is None:
                            found_pt = P
                            break
                if found_pt:
                    return p, n, b_try, lam, found_pt
    return None

# ---------------------------------------------------------------------------
# CF analysis — max_a = max partial quotient with q_k ≤ 300
# ---------------------------------------------------------------------------

def cf_analysis(lam, n, q_cutoff=300):
    """Returns (max_q, n_convs, max_a) for the CF of lam/n up to q_cutoff."""
    mq, count, max_a = 1, 0, 1
    p0, p1 = 1, 0
    q0, q1 = 0, 1
    x, y = lam, n
    for _ in range(200):
        if y == 0: break
        a = x // y
        x, y = y, x - a * y
        p_k = a * p0 + p1
        q_k = a * q0 + q1
        if q_k > q_cutoff: break
        mq = q_k; count += 1; max_a = max(max_a, a)
        p1, p0 = p0, p_k
        q1, q0 = q0, q_k
    return mq, count, max_a

# ---------------------------------------------------------------------------
# Pure Python LLL reduction
# ---------------------------------------------------------------------------

def lll_reduce(M, delta=0.75):
    """
    LLL reduction with Gram-Schmidt coefficients computed in float64.
    M: list of lists (integer row vectors), modified in-place.
    Returns M (same object).
    """
    n = len(M)
    if n == 0: return M
    m_dim = len(M[0])

    # Compute GS coefficients and norms
    B_star = [list(M[i]) for i in range(n)]  # float copies
    mu = [[0.0] * n for _ in range(n)]
    B_norm = [0.0] * n  # ||b*_i||^2

    def gs_update(k):
        """Recompute B_star[k] and mu[k][j] for j<k, and B_norm[k]."""
        v = [float(x) for x in M[k]]
        for j in range(k):
            dot = sum(v[l] * B_star[j][l] for l in range(m_dim))
            if B_norm[j] > 0:
                mu[k][j] = dot / B_norm[j]
            else:
                mu[k][j] = 0.0
            for l in range(m_dim):
                v[l] -= mu[k][j] * B_star[j][l]
        B_star[k] = v
        B_norm[k] = sum(x*x for x in v)

    # Initial GS
    for i in range(n):
        gs_update(i)

    def size_reduce(k, j):
        """Size-reduce b_k with respect to b_j."""
        q = round(mu[k][j])
        if q == 0: return
        for l in range(m_dim):
            M[k][l] -= q * M[j][l]
        # Update mu[k][i] for i <= j
        for i in range(j + 1):
            mu[k][i] -= q * mu[j][i]
        mu[k][j] -= q  # correct for i == j case

    k = 1
    max_iter = 200 * n * n
    it = 0
    while k < n and it < max_iter:
        it += 1
        # Size reduce b_k against b_{k-1}
        for j in range(k - 1, -1, -1):
            size_reduce(k, j)
        # Lovász condition
        if B_norm[k] >= (delta - mu[k][k-1]**2) * B_norm[k-1]:
            k += 1
        else:
            # Swap b_k and b_{k-1}
            M[k], M[k-1] = M[k-1], M[k]
            # Recompute GS from k-1
            for i in range(max(0, k-1), n):
                gs_update(i)
            k = max(1, k - 1)
    return M

# ---------------------------------------------------------------------------
# GLV-HNP LLL attack
# ---------------------------------------------------------------------------

def lll_attack(p, n, b_param, lam, G, m, K1, seed):
    """Returns (recovered_bool, d_secret)."""
    K2 = math.isqrt(n) + 1
    S_K1 = max(1, n // K1)
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    dim = 2 * m + 2

    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs = []
    for _ in range(m * 10):
        if len(sigs) >= m: break
        k1 = rng.randint(0, K1 - 1)
        k2 = rng.randint(0, K2 - 1)
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
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})

    if len(sigs) < m:
        return False, d_secret

    # Build lattice matrix
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = 1
    for i in range(m):
        M[m+1+i][i] = (-lam * S_K1) % (n * S_K1)  # keep non-negative for GS stability
        M[m+1+i][m+1+i] = S_K2
    # Correct: use negative lam
    for i in range(m):
        M[m+1+i][i] = (-lam * S_K1 + n * S_K1 * (lam * S_K1 // (n * S_K1) + 1))
    # Re-do: just use negative integers directly
    for i in range(m):
        M[m+1+i][i] = -lam * S_K1
        M[m+1+i][m+1+i] = S_K2
    for i in range(m):
        M[2*m+1][i] = sigs[i]['A'] * S_K1
    M[2*m+1][dim-1] = S_KANNAN

    lll_reduce(M)

    for row in M:
        last = row[dim-1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_secret:
            return True, d_secret

    return False, d_secret

# ---------------------------------------------------------------------------
# Exp S main
# ---------------------------------------------------------------------------

def exp_s():
    print("=" * 72)
    print("Exp S: max_a (max partial quotient) as C1-predictor — 50 fresh curves")
    print("  RULE: max_a ≤ 17 → predict C2  (LLL succeeds at K1=72)")
    print("        max_a ≥ 18 → predict C1  (LLL fails at K1=72)")
    print("  Based on Exp R: C1 max_a ∈ [3,85], C2 max_a ∈ [4,17] (known overlap)")
    print("=" * 72)
    print()

    M_SIGS = 12
    K1 = 72
    N_SEEDS = 6
    THRESHOLD = 18  # max_a >= THRESHOLD → predict C1

    print("Generating 50 fresh CM curves (Cornacchia)...", flush=True)
    curves = []
    for seed in range(2000):
        if len(curves) >= 50: break
        r = gen_cm_curve(seed)
        if r is None: continue
        p, n, b, lam, G = r
        curves.append((p, n, b, lam, G))
    print(f"Generated {len(curves)} curves.")
    print()

    print(f"  {'#':>3}  {'n':>9}  {'λ/n':>7}  {'max_a':>6}  "
          f"{'predict':>8}  {'actual':>8}  {'ok?':>5}  rec/seeds",
          flush=True)
    print("  " + "-" * 65)

    results = []
    for idx, (p, n, b, lam, G) in enumerate(curves):
        _mq, _nc, max_a = cf_analysis(lam, n, q_cutoff=300)
        predict_c1 = (max_a >= THRESHOLD)

        recovered = 0
        for seed in range(N_SEEDS):
            ok, _ = lll_attack(p, n, b, lam, G, M_SIGS, K1, seed)
            if ok: recovered += 1
        actual_c1 = (recovered < N_SEEDS // 2)

        match = (predict_c1 == actual_c1)
        pred_str = "C1" if predict_c1 else "C2"
        act_str  = "C1" if actual_c1  else "C2"
        ok_str   = "OK" if match else "MISS"

        print(f"  {idx+1:>3}  {n:>9}  {lam/n:>7.4f}  {max_a:>6}  "
              f"{pred_str:>8}  {act_str:>8}  {ok_str:>5}  {recovered}/{N_SEEDS}",
              flush=True)

        results.append({
            'n': n, 'lam': lam, 'max_a': max_a,
            'predict_c1': predict_c1, 'actual_c1': actual_c1,
            'match': match, 'recovered': recovered
        })

    print()
    n_total = len(results)
    n_correct = sum(r['match'] for r in results)
    acc = n_correct / n_total * 100 if n_total else 0
    print(f"ACCURACY: {n_correct}/{n_total} = {acc:.1f}%")

    c1_actual = [r for r in results if r['actual_c1']]
    c2_actual = [r for r in results if not r['actual_c1']]
    c1_pred_ok = [r for r in c1_actual if r['predict_c1']]
    c2_pred_ok = [r for r in c2_actual if not r['predict_c1']]
    print(f"C1 sensitivity: {len(c1_pred_ok)}/{len(c1_actual)}")
    print(f"C2 specificity: {len(c2_pred_ok)}/{len(c2_actual)}")

    c1_ma = sorted(r['max_a'] for r in c1_actual)
    c2_ma = sorted(r['max_a'] for r in c2_actual)
    print(f"\nmax_a distributions:")
    print(f"  C1 (LLL fails): {c1_ma}")
    print(f"  C2 (LLL ok):    {c2_ma}")

    misses = [r for r in results if not r['match']]
    if misses:
        print(f"\nMisses ({len(misses)}):")
        for r in misses:
            print(f"  n={r['n']}, λ/n={r['lam']/r['n']:.4f}, max_a={r['max_a']}, "
                  f"pred={'C1' if r['predict_c1'] else 'C2'}, "
                  f"actual={'C1' if r['actual_c1'] else 'C2'}, "
                  f"rec={r['recovered']}/{N_SEEDS}")

    # K1_threshold regression attempt: is K1_threshold ≈ f(max_a)?
    # Proxy: recovered count at K1=72 as a surrogate for continuous threshold.
    # Correlation: max_a vs recovered count.
    print(f"\nmax_a vs recovery count at K1=72:")
    for r in sorted(results, key=lambda x: x['max_a']):
        flag = " ← C1" if r['actual_c1'] else ""
        print(f"  max_a={r['max_a']:>4}  rec={r['recovered']}/{N_SEEDS}{flag}")

    return results, c1_actual, c2_actual

if __name__ == "__main__":
    t0 = time.time()
    print()
    print("GLV-HNP Phase 2: Exp S — max_a predictor (2026-06-29)")
    print()

    results, c1, c2 = exp_s()

    print(f"\nTotal wall time: {time.time() - t0:.1f}s")
