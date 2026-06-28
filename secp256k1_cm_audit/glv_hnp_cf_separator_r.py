"""
GLV-HNP Phase 2: Exp R (standalone) — max_q_cf separator on 30 fresh curves (2026-06-28)

REVISED HYPOTHESIS from Exp P:
  The correct separator between C1-type (K1_threshold~40-56) and C2-type
  (K1_threshold~128-200) is not 'q_cf in [20,100]' (FALSIFIED), but:

  max_q_cf(lam/n, cutoff=300) = max { q_k : q_k ≤ 300 } where (p_k, q_k) range
  over the CF convergents of the rational lam/n.

  FINDING from Exp P (9 known curves):
    C1-type: max_q_cf ∈ {194, 147, 227, 135, 111, 203}  ← all ≥ 111
    C2-type: max_q_cf ∈ {35, 23, 25}                     ← all ≤ 35
    Gap: [35, 111] is empty.

  PREDICTION RULE (to test here):
    max_q_cf > 60  →  predict C1-type  (LLL fails at K1=72)
    max_q_cf ≤ 60  →  predict C2-type  (LLL succeeds at K1=72)

  MECHANISM HYPOTHESIS:
    C1 curves have small CF partial quotients → denominator grows slowly → many
    convergents ≤ 300. C2 curves have at least one large partial quotient → denominator
    jumps past 300 quickly. This large partial quotient corresponds to a very good
    rational approximation lam/n ≈ p/q at small q, which creates a "structural
    stability" in the GLV lattice (the target vector is protected from being
    displaced by spurious short vectors at moderate K1).

FAST CURVE GENERATION (Cornacchia-based):
  For j=0 curves over F_p (p≡1 mod 3): Cornacchia gives (a,b) with a²+3b²=4p.
  Candidate orders: n1=p+1-a, n2=p+1+a. Verify primality, then check via ec_mul.
  Cost: O(log p) per candidate vs O(p) for naive point counting.

Run: python3 glv_hnp_cf_separator_r.py
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Primality (Miller-Rabin)
# ---------------------------------------------------------------------------

def is_prime_miller_rabin(n, rounds=12):
    if n < 2: return False
    if n == 2 or n == 3: return True
    if n % 2 == 0: return False
    small = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small:
        if n == p: return True
        if n % p == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    witnesses = small[:rounds] if n > 37 else []
    for a in witnesses:
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

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0: r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return min(r1, r2)

def find_point_on_curve(p, b):
    """Find any point on y²=x³+b over F_p."""
    for x in range(p):
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

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
# Cornacchia: solve a² + 3b² = 4p (for j=0 CM curves)
# ---------------------------------------------------------------------------

def cornacchia_j0(p):
    """Find (a, b) with a > 0, b > 0, a²+3b²=4p, a≡2 mod 3.
       Returns None if no solution exists.
       Reference: ALGORITHM 1.5.3 in Cohen 'CANT'."""
    if p == 2: return None
    if p % 3 != 1: return None
    # Solve x²≡-3 mod p
    neg3 = p - 3
    if pow(neg3, (p-1)//2, p) != 1:
        return None
    x0 = tonelli_shanks(neg3, p)
    if x0 is None:
        return None
    # Normalize: we want x0 to be the sqrt closest to the right range
    # Run Euclidean algorithm until remainder <= sqrt(4p)
    limit = math.isqrt(4 * p)
    u, v = 2 * p, x0  # start from (2p, x0)
    while v > limit:
        u, v = v, u % v
    # Now v should satisfy v² < 4p and 4p ≡ v² mod p*something
    a = v
    b_sq_num = 4 * p - a * a
    if b_sq_num <= 0 or b_sq_num % 3 != 0:
        # Try a=2p - a case
        a = 2*p - v
        b_sq_num = 4 * p - a * a
        if b_sq_num <= 0 or b_sq_num % 3 != 0:
            return None
    b_sq = b_sq_num // 3
    b = math.isqrt(b_sq)
    if b * b != b_sq:
        return None
    return abs(a), b

# ---------------------------------------------------------------------------
# Fast CM curve generation
# ---------------------------------------------------------------------------

KNOWN_N_SET = {523597, 345109, 326479, 337999, 321163, 949423, 904369, 811297, 319129}

def gen_cm_curve_fast(seed, bit_lo=17, bit_hi=20):
    """Generate a random j=0 CM curve using Cornacchia.
       Returns (p, n, b_param, lam) or None.
    """
    rng = random.Random(seed * 7919 + 31337)
    for _ in range(5000):
        bits = rng.randint(bit_lo, bit_hi)
        p_cand = rng.randint(2**(bits-1), 2**bits - 1)
        p_cand |= 1
        if p_cand % 3 != 1:
            p_cand += 2 if p_cand % 3 == 2 else 0
        if p_cand % 2 == 0:
            p_cand += 1
        if not is_prime_miller_rabin(p_cand):
            continue
        if p_cand % 3 != 1:
            continue

        corn = cornacchia_j0(p_cand)
        if corn is None:
            continue
        a_corn, b_corn = corn

        # Candidate group orders: p+1±a
        for sign in [1, -1]:
            n_cand = p_cand + 1 - sign * a_corn
            if n_cand < 4:
                continue
            if not is_prime_miller_rabin(n_cand):
                continue
            if n_cand in KNOWN_N_SET:
                continue
            # Check GLV eigenvalue
            lam = glv_eigenvalue(n_cand)
            if lam is None:
                continue

            # Find b_param: try small values, verify order via ec_mul
            for b_try in range(1, 30):
                P = find_point_on_curve(p_cand, b_try)
                if P is None:
                    continue
                # Fast order check: does n_cand * P = O?
                if ec_mul(P, n_cand, p_cand) is None:
                    # Confirm full generator (not cofactor issue)
                    # n_cand is prime so any non-identity point is a generator
                    return p_cand, n_cand, b_try, lam
    return None

# ---------------------------------------------------------------------------
# CF analysis
# ---------------------------------------------------------------------------

def cf_convergents(num, den, max_terms=100):
    """Return list of (p_k, q_k) convergents of num/den."""
    convergents = []
    a_prev, a_curr = 1, 0
    b_prev, b_curr = 0, 1
    x, y = num, den
    for _ in range(max_terms):
        if y == 0:
            break
        a = x // y
        x, y = y, x - a * y
        p_k = a * a_prev + a_curr
        q_k = a * b_prev + b_curr
        convergents.append((p_k, q_k, a))
        a_curr, a_prev = a_prev, p_k
        b_curr, b_prev = b_prev, q_k
    return convergents

def max_q_cf(lam, n, cutoff=None):
    """Return max CF convergent denominator of lam/n (up to cutoff, or unlimited)."""
    if cutoff is None:
        cutoff = n
    convs = cf_convergents(lam, n)
    mq = 1
    for _, q_k, _ in convs:
        if q_k > cutoff:
            break
        mq = q_k
    return mq

def max_partial_quotient(lam, n, max_q_cutoff=300):
    """Return the maximum partial quotient a_k before q_k exceeds max_q_cutoff."""
    convs = cf_convergents(lam, n)
    max_a = 1
    for _, q_k, a_k in convs:
        if q_k > max_q_cutoff:
            break
        max_a = max(max_a, a_k)
    return max_a

def cf_info(lam, n, q_cutoff=300):
    """Return (max_q, n_convs_in_range, max_partial_quot, all_qs)."""
    convs = cf_convergents(lam, n)
    mq = 1
    count = 0
    max_a = 1
    all_qs = []
    for _, q_k, a_k in convs:
        if q_k > q_cutoff:
            break
        mq = q_k
        count += 1
        max_a = max(max_a, a_k)
        all_qs.append(q_k)
    return mq, count, max_a, all_qs

# ---------------------------------------------------------------------------
# LLL attack
# ---------------------------------------------------------------------------

def gen_signatures(p, b_param, n, lam, G, m_sigs, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs, k1_vals, k2_vals = [], [], []
    attempts = 0
    while len(sigs) < m_sigs and attempts < 500000:
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
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})
        k1_vals.append(k1); k2_vals.append(k2)
    return d_secret, sigs

def lll_attack(p, b_param, n, lam, G, m, K1, seed):
    d_sec, sigs = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False, d_sec

    K2 = math.isqrt(n) + 1
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    dim = 2 * m + 2

    M_int = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M_int[i][i] = n * S_K1
    for i in range(m):
        M_int[m][i] = sigs[i]['B'] * S_K1
    M_int[m][m] = 1
    for i in range(m):
        M_int[m + 1 + i][i] = -lam * S_K1
        M_int[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M_int[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M_int[2 * m + 1][dim - 1] = S_KANNAN

    A = IntegerMatrix.from_matrix(M_int)
    LLL.reduction(A)
    M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == d_sec:
            return True, d_sec
    return False, d_sec

# ---------------------------------------------------------------------------
# Exp R: test max_q_cf separator on 30 fresh curves
# ---------------------------------------------------------------------------

def exp_r():
    print("=" * 72)
    print("Exp R: max_q_cf separator validated on 30 fresh CM curves")
    print("  RULE: max_q_cf(lam/n, cutoff=300) > 60 → predict C1 (LLL fails at K1=72)")
    print("        max_q_cf ≤ 60                → predict C2 (LLL succeeds at K1=72)")
    print("=" * 72)
    print()

    m = 12
    K1_TEST = 72
    N_SEEDS = 6
    MAX_Q_CUTOFF = 300
    SEPARATOR = 60

    print("  Generating 30 fresh CM curves (Cornacchia method)...", flush=True)

    fresh_curves = []
    for seed in range(500):
        if len(fresh_curves) >= 30:
            break
        result = gen_cm_curve_fast(seed)
        if result is None:
            continue
        fresh_curves.append(result)

    print(f"  Generated {len(fresh_curves)} curves.")
    print()

    print(f"  {'#':>3}  {'n':>9}  {'λ/n':>7}  {'max_q':>7}  {'max_a':>7}  "
          f"{'predict':>8}  {'actual':>8}  {'ok?':>5}  {'rec'}",
          flush=True)
    print(f"  {'-'*3}  {'-'*9}  {'-'*7}  {'-'*7}  {'-'*7}  "
          f"{'-'*8}  {'-'*8}  {'-'*5}  {'-'*6}")

    n_total = 0
    n_correct = 0
    results = []

    for i, (p_c, n_c, b_c, lam_c) in enumerate(fresh_curves):
        G_c = find_generator(p_c, b_c, n_c)
        if G_c is None:
            continue

        # CF analysis
        mq, n_convs, max_a, all_qs = cf_info(lam_c, n_c, q_cutoff=MAX_Q_CUTOFF)

        # Prediction
        predict_c1 = (mq > SEPARATOR)

        # Actual LLL classification
        recovered = 0
        for seed in range(N_SEEDS):
            ok, _ = lll_attack(p_c, b_c, n_c, lam_c, G_c, m, K1_TEST, seed)
            if ok:
                recovered += 1
        actual_c1 = (recovered < N_SEEDS // 2)

        match = (predict_c1 == actual_c1)
        n_correct += int(match)
        n_total += 1

        pred_str = "C1" if predict_c1 else "C2"
        act_str = "C1" if actual_c1 else "C2"
        ok_str = "OK" if match else "MISS"

        print(f"  {i+1:>3}  {n_c:>9}  {lam_c/n_c:>7.4f}  {mq:>7}  {max_a:>7}  "
              f"{pred_str:>8}  {act_str:>8}  {ok_str:>5}  {recovered}/{N_SEEDS}",
              flush=True)

        results.append({
            'p': p_c, 'n': n_c, 'b': b_c, 'lam': lam_c,
            'mq': mq, 'max_a': max_a, 'all_qs': all_qs,
            'predict_c1': predict_c1, 'actual_c1': actual_c1, 'match': match,
            'recovered': recovered
        })

    print()
    if n_total > 0:
        acc = n_correct / n_total * 100
        print(f"  ACCURACY: {n_correct}/{n_total} = {acc:.1f}%")

        c1_actual = [r for r in results if r['actual_c1']]
        c2_actual = [r for r in results if not r['actual_c1']]
        c1_pred_ok = [r for r in c1_actual if r['predict_c1']]
        c2_pred_ok = [r for r in c2_actual if not r['predict_c1']]
        print(f"  C1 sensitivity: {len(c1_pred_ok)}/{len(c1_actual)}")
        print(f"  C2 specificity: {len(c2_pred_ok)}/{len(c2_actual)}")

        misses = [r for r in results if not r['match']]
        if misses:
            print(f"\n  MISSES:")
            for r in misses:
                print(f"    n={r['n']}, λ/n={r['lam']/r['n']:.4f}, "
                      f"max_q={r['mq']}, max_a={r['max_a']}, "
                      f"pred={'C1' if r['predict_c1'] else 'C2'}, "
                      f"actual={'C1' if r['actual_c1'] else 'C2'}, "
                      f"qs={r['all_qs']}")
        else:
            print("\n  No misses — separator rule is perfect on this sample.")

        # Distribution of max_q for C1 vs C2
        c1_mq = sorted([r['mq'] for r in c1_actual])
        c2_mq = sorted([r['mq'] for r in c2_actual])
        print(f"\n  max_q distribution:")
        print(f"    C1 (fail) max_q values: {c1_mq}")
        print(f"    C2 (ok)   max_q values: {c2_mq}")

        # Max partial quotient distribution
        c1_ma = sorted([r['max_a'] for r in c1_actual])
        c2_ma = sorted([r['max_a'] for r in c2_actual])
        print(f"\n  max_partial_quotient distribution:")
        print(f"    C1 (fail): {c1_ma}")
        print(f"    C2 (ok):   {c2_ma}")

    return results


# ---------------------------------------------------------------------------
# Exp R2: fine-scan of the separator boundary
# ---------------------------------------------------------------------------

def exp_r2(results):
    """Use the Exp R data to refine the max_q separator threshold."""
    if not results:
        return
    print()
    print("=" * 72)
    print("Exp R2: Refine separator boundary from Exp R data")
    print("=" * 72)
    print()

    c1 = sorted([r['mq'] for r in results if r['actual_c1']])
    c2 = sorted([r['mq'] for r in results if not r['actual_c1']])

    if c1 and c2:
        min_c1 = min(c1)
        max_c2 = max(c2)
        print(f"  C1 min max_q = {min_c1}")
        print(f"  C2 max max_q = {max_c2}")
        print(f"  Gap: ({max_c2}, {min_c1})")

        if max_c2 < min_c1:
            mid = (max_c2 + min_c1) // 2
            print(f"  Optimal separator: max_q ≤ {max_c2} → C2, max_q ≥ {min_c1} → C1")
            print(f"  Suggested threshold: {mid}")
        else:
            print("  OVERLAP — separator is not clean; see misses above.")

    print()
    # Check max_a as separator
    c1_ma = sorted([r['max_a'] for r in results if r['actual_c1']])
    c2_ma = sorted([r['max_a'] for r in results if not r['actual_c1']])
    if c1_ma and c2_ma:
        print(f"  max_partial_quotient as separator:")
        print(f"    C1 range: [{min(c1_ma)}, {max(c1_ma)}]")
        print(f"    C2 range: [{min(c2_ma)}, {max(c2_ma)}]")
        if max(c1_ma) < min(c2_ma):
            print(f"    CLEAN: max_a ≤ {max(c1_ma)} → C1, max_a ≥ {min(c2_ma)} → C2")
        elif max(c2_ma) < min(c1_ma):
            print(f"    CLEAN: max_a ≤ {max(c2_ma)} → C2, max_a ≥ {min(c1_ma)} → C1")
        else:
            print(f"    OVERLAP in max_partial_quotient — not a clean separator")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time.time()
    print()
    print("GLV-HNP Phase 2: Exp R — max_q_cf separator (revised hypothesis)")
    print("  Date: 2026-06-28")
    print()

    results = exp_r()
    exp_r2(results)

    elapsed = time.time() - t0
    print(f"Total wall time: {elapsed:.1f}s")
