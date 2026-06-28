"""
GLV-HNP Phase 2: Exp P/Q — CF-convergent hypothesis for K1 threshold (2026-06-28)

HYPOTHESIS (from 2026-06-27 findings):
  C1-type curves (K1_threshold ~40-56) may have a CF convergent p/q of lam/n
  with q ≈ K1_threshold and q in [32,100]. C2-type curves (K1_threshold ~128-200)
  would have no such small convergent.

  MECHANISM: if lam/n ≈ p/q (CF convergent), then the vector
      v_spurious = q*(row_m+1+i) + p*(row_i)
  in the GLV lattice has a near-zero first component (K1-part), since
  p*n*S_K1 - q*lam*S_K1 ≈ 0. This vector has short K1-part and K2-part ≈ q*S_K2.
  When K1 ≈ q, this spurious vector's norm ≈ solution norm → LLL can't separate them.

Exp P: For every known curve (C1 main, 5 Exp-N C1-type, 3 Exp-O C2-type, 8 Exp-L
  additional), compute the CF expansion of lam/n and find the first convergent q_k
  in [20,300]. Compare q_k with empirical K1_threshold.

Exp Q: For C1 main (p=524743, lam_thresh≈46), at K1=46 (7/20 recovery), extract
  the short spurious lattice row and decode its k2 components. Check if they ≈ q_cf
  (the CF convergent denominator), confirming the spurious-vector mechanism.

Exp R: Extend to 20 fresh C1-type and C2-type curves — compute q_cf and use it
  to PREDICT C1/C2 classification before running LLL. Accuracy of q_cf < 100
  as C1-predictor tested on held-out curves.

Run: python3 glv_hnp_cf_separator.py
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (same as prior scripts)
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
# CF convergents of lam/n
# ---------------------------------------------------------------------------

def cf_convergents(num, den, max_terms=50):
    """Return list of (p_k, q_k) convergents of num/den."""
    convergents = []
    a_prev, a_curr = 1, 0   # p_{k-1}, p_{k-2}
    b_prev, b_curr = 0, 1   # q_{k-1}, q_{k-2}
    x, y = num, den
    for _ in range(max_terms):
        if y == 0:
            break
        a = x // y
        x, y = y, x - a * y
        p_k = a * a_prev + a_curr
        q_k = a * b_prev + b_curr
        convergents.append((p_k, q_k))
        a_curr, a_prev = a_prev, p_k
        b_curr, b_prev = b_prev, q_k
    return convergents

def first_convergent_in_range(lam, n, lo=20, hi=300):
    """Find the first CF convergent q_k of lam/n with lo <= q_k <= hi.
       Returns (p_k, q_k, residual) where residual = |lam*q_k - p_k*n|.
    """
    convs = cf_convergents(lam, n)
    for p_k, q_k in convs:
        if lo <= q_k <= hi:
            residual = abs(lam * q_k - p_k * n)
            return p_k, q_k, residual
    return None, None, None

def cf_summary(lam, n, max_q=300):
    """Return all convergents of lam/n with q_k <= max_q."""
    convs = cf_convergents(lam, n)
    result = []
    for p_k, q_k in convs:
        if q_k > max_q:
            break
        residual = abs(lam * q_k - p_k * n)
        result.append((p_k, q_k, residual))
    return result

# ---------------------------------------------------------------------------
# LLL attack (minimally adapted from glv_hnp_k1_threshold.py)
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
    return d_secret, sigs, k1_vals, k2_vals

def build_lattice(n, lam, sigs, m, K1):
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
    return M_int, S_K1, S_K2, S_KANNAN, dim

def lll_attack_full(p, b_param, n, lam, G, m, K1, seed):
    """Returns (recovered, d_cand, nom_norm, found_norm, d_sec, best_row, reduced_mat)."""
    d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False, None, None, None, d_sec, None, None

    K2 = math.isqrt(n) + 1
    M_int, S_K1, S_K2, S_KANNAN, dim = build_lattice(n, lam, sigs, m, K1)

    nom_norm_sq = (sum(k * k * S_K1 * S_K1 for k in k1v) +
                   d_sec * d_sec +
                   sum(k * k * S_K2 * S_K2 for k in k2v) +
                   S_KANNAN * S_KANNAN)
    nom_norm = math.sqrt(nom_norm_sq)

    A = IntegerMatrix.from_matrix(M_int)
    LLL.reduction(A)
    M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

    best_d = None; best_norm = float('inf'); best_row = None
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        rn = math.sqrt(sum(x * x for x in row))
        if rn < best_norm:
            best_norm = rn
            best_d = d_cand
            best_row = row

    if best_d is None:
        return False, None, nom_norm, None, d_sec, None, M_red

    return best_d == d_sec, best_d, nom_norm, best_norm, d_sec, best_row, M_red

def lll_attack(p, b_param, n, lam, G, m, K1, seed):
    ok, d_cand, nom, found, d_sec, _, _ = lll_attack_full(p, b_param, n, lam, G, m, K1, seed)
    return ok, d_cand, nom, found, d_sec

# ---------------------------------------------------------------------------
# Exp P: CF analysis on all known curves
# ---------------------------------------------------------------------------

def exp_p_cf_analysis():
    print("=" * 70)
    print("Exp P: CF-convergent analysis on all known C1/C2 curves")
    print("  HYPOTHESIS: C1-type ↔ first CF convergent q_k of lam/n in [20,100]")
    print("  C2-type ↔ no convergent in [20,100], first one in [100,300]")
    print("=" * 70)
    print()

    # All known curves with empirical K1_threshold and known C1/C2 classification
    # (class: 1=C1-type/fails, 2=C2-type/succeeds, at K1=72 from Exp L)
    curves = [
        # label, p, n, b, K1_threshold (empirical), class
        ("C1 main", 524743, 523597, None, 46, 1),   # Exp M
        # Exp N (C1-type):
        ("#2",      343963, 345109,  2, 40, 1),
        ("#9",      327553, 326479, 10, 48, 1),
        ("#15",     339139, 337999,  2, 40, 1),
        ("#5",      321889, 321163, 19, 40, 1),
        ("#25",     951001, 949423, 11, 56, 1),
        # Exp O (C2-type):
        ("#8",      903367, 904369,  3, 200, 2),
        ("#14",     812233, 811297, 15, 200, 2),
        ("#11",     318001, 319129, 13, 128, 2),
    ]

    # Find b_param for C1 main
    for entry in curves:
        if entry[0] == "C1 main":
            p, n = entry[1], entry[2]
            for b_try in range(1, 500):
                rng2 = random.Random(99)
                found_b = False
                for _ in range(20):
                    x = rng2.randint(0, p - 1)
                    rhs = (pow(x, 3, p) + b_try) % p
                    y = tonelli_shanks(rhs, p)
                    if y is not None and y != 0:
                        P = (x, y)
                        if ec_mul(P, n, p) is None:
                            found_b = True
                            break
                if found_b:
                    curves[0] = ("C1 main", p, n, b_try, 46, 1)
                    break
            break

    print(f"  {'Curve':>8}  {'class':>5}  {'lam/n':>8}  {'q_cf_small':>12}  {'residual':>12}  "
          f"{'K1_thresh':>10}  {'q_cf≈thresh?':>14}")
    print(f"  {'-'*8}  {'-'*5}  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*14}")

    results = []
    for (label, p, n, b_param, k1_thresh, cls) in curves:
        lam = glv_eigenvalue(n)
        if lam is None:
            print(f"  {label:>8}  no GLV eigenvalue, skip")
            continue

        # CF expansion: find first convergent with q in [20, 100]
        p_cf, q_cf, residual = first_convergent_in_range(lam, n, lo=20, hi=100)
        all_convs = cf_summary(lam, n, max_q=300)

        approx_match = "YES" if (q_cf is not None and abs(q_cf - k1_thresh) <= 10) else "NO"
        q_str = str(q_cf) if q_cf is not None else "none[20-100]"
        res_str = str(residual) if residual is not None else "---"

        print(f"  {label:>8}  C{cls}    {lam/n:>8.4f}  {q_str:>12}  {res_str:>12}  "
              f"{str(k1_thresh):>10}  {approx_match:>14}")

        results.append({
            'label': label, 'p': p, 'n': n, 'lam': lam,
            'lamn': lam/n, 'cls': cls,
            'k1_thresh': k1_thresh,
            'q_cf': q_cf, 'p_cf': p_cf, 'residual': residual,
            'all_convs': all_convs,
        })

    print()
    print("  Full CF convergent tables:")
    for r in results:
        convs = r['all_convs']
        print(f"\n  Curve {r['label']} (C{r['cls']}, λ/n={r['lamn']:.5f}):")
        print(f"    Convergents: ", end="")
        parts = [f"q={q}" for (p2, q, res) in convs if q <= 300]
        print(", ".join(parts[:10]) if parts else "(none ≤300)")

    return results


# ---------------------------------------------------------------------------
# Exp Q: Spurious-row anatomy at K1=46 for C1 main
# ---------------------------------------------------------------------------

def exp_q_spurious_anatomy():
    print()
    print("=" * 70)
    print("Exp Q: Spurious-row anatomy at K1=46 (C1 main)")
    print("  At K1=46 (7/20 success), extract the 'shortest Kannan row' from")
    print("  failing seeds. Decode k1_i, k2_i components. Check if k2 ≈ q_cf.")
    print("=" * 70)
    print()

    p, n = 524743, 523597
    # Find b
    b_param = None
    for b_try in range(1, 500):
        rng2 = random.Random(99)
        for _ in range(20):
            x = rng2.randint(0, p - 1)
            rhs = (pow(x, 3, p) + b_try) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                P = (x, y)
                if ec_mul(P, n, p) is None:
                    b_param = b_try
                    break
        if b_param is not None:
            break

    lam = glv_eigenvalue(n)
    G = find_generator(p, b_param, n)
    K1 = 46
    m = 12
    K2 = math.isqrt(n) + 1

    print(f"  p={p}, n={n}, λ={lam}, λ/n={lam/n:.4f}, K1={K1}, K2={K2}")

    # CF convergent for reference
    p_cf, q_cf, residual = first_convergent_in_range(lam, n, lo=20, hi=100)
    print(f"  CF convergent: p/q = {p_cf}/{q_cf}, residual |λ*q - p*n| = {residual}")
    print()

    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n

    print(f"  Scales: S_K1={S_K1}, S_K2={S_K2}, S_KANNAN={S_KANNAN}")
    print()

    FAIL_SEEDS = []
    for seed in range(20):
        ok, d_cand, nom, found, d_sec, best_row, M_red = lll_attack_full(
            p, b_param, n, lam, G, m, K1, seed)
        if not ok and best_row is not None:
            FAIL_SEEDS.append((seed, d_sec, best_row, M_red))
        if len(FAIL_SEEDS) >= 4:
            break

    print(f"  Collected {len(FAIL_SEEDS)} failing seeds for anatomy.")
    print()

    for (seed, d_sec, srow, M_red) in FAIL_SEEDS[:3]:
        print(f"  --- Seed {seed} (d_secret={d_sec}) ---")
        dim = 2 * m + 2
        # Decode the spurious row
        sign = 1 if srow[dim-1] > 0 else -1
        k1_parts = [sign * srow[i] / S_K1 for i in range(m)]
        d_part = sign * srow[m]
        k2_parts = [sign * srow[m+1+i] / S_K2 for i in range(m)]
        kannan_part = srow[dim-1]

        k2_scaled = [sign * srow[m+1+i] for i in range(m)]
        k2_gcd = math.gcd(*[abs(x) for x in k2_scaled if x != 0])

        print(f"    k1 components (raw/S_K1): {[f'{v:.2f}' for v in k1_parts[:5]]}...")
        print(f"    k2 components (raw/S_K2): {[f'{v:.2f}' for v in k2_parts[:5]]}...")
        print(f"    d component: {d_part} (should NOT be {d_sec})")
        print(f"    k2 raw values: {k2_scaled[:5]}...")
        print(f"    k2 gcd: {k2_gcd} (≈ q_cf={q_cf}?)")
        print(f"    k2_i / S_K2 magnitudes: {[abs(v) for v in k2_parts[:5]]}")

        # Check if k2_parts look like multiples of q_cf
        if q_cf is not None:
            k2_approx = [abs(v) / q_cf for v in k2_parts if v != 0]
            print(f"    |k2_i / S_K2| / q_cf: {[f'{v:.3f}' for v in k2_approx[:5]]}")

        row_norm = math.sqrt(sum(x*x for x in srow))
        print(f"    Row norm: {row_norm:.1f}")
        print()


# ---------------------------------------------------------------------------
# Exp R: Prediction power of q_cf on 20 fresh curves
# ---------------------------------------------------------------------------

def exp_r_prediction_accuracy():
    print("=" * 70)
    print("Exp R: q_cf < 100 as predictor of C1-classification (20 fresh curves)")
    print("  Generate 20 random CM toy curves; compute q_cf; run LLL at K1=72;")
    print("  measure: does q_cf < 100 predict LLL failure (C1-type)?")
    print("=" * 70)
    print()

    # Generate fresh CM curves (j=0, y²=x³+b over F_p, p≡1 mod 3, n = prime)
    def is_prime(x):
        if x < 2: return False
        if x % 2 == 0: return x == 2
        for f in range(3, min(int(x**0.5)+1, 1000), 2):
            if x % f == 0: return False
        return True

    def gen_cm_curve(seed):
        rng = random.Random(seed * 1337 + 42)
        for _ in range(5000):
            # Try a random prime p ≡ 1 mod 3, 18-20 bits
            bits = rng.randint(18, 20)
            p_cand = rng.randint(2**(bits-1), 2**bits - 1)
            p_cand |= 1  # make odd
            if p_cand % 3 != 1: continue
            if not is_prime(p_cand): continue
            # Count curve orders for small b values
            for b_try in range(1, 20):
                # Heuristic: try to verify we get prime order subgroup
                # Use b and the standard construction
                # Actually, just try a few points and use naive check
                # For a small prime, count all points
                count = 1  # point at infinity
                for x in range(p_cand):
                    rhs = (pow(x, 3, p_cand) + b_try) % p_cand
                    if rhs == 0:
                        count += 1
                    elif pow(rhs, (p_cand-1)//2, p_cand) == 1:
                        count += 2
                n_cand = count
                if n_cand < 4: continue
                if not is_prime(n_cand): continue
                # Check GLV eigenvalue exists
                lam = glv_eigenvalue(n_cand)
                if lam is None: continue
                # Check b gives a valid generator
                G = find_generator(p_cand, b_try, n_cand)
                if G is None: continue
                return p_cand, n_cand, b_try, lam
        return None

    print("  Generating 20 fresh CM curves (may take ~30s)...")
    print()

    fresh_curves = []
    seed_idx = 0
    while len(fresh_curves) < 20 and seed_idx < 200:
        result = gen_cm_curve(seed_idx)
        seed_idx += 1
        if result is None:
            continue
        p_c, n_c, b_c, lam_c = result
        # Skip if already in known set
        if n_c in {523597, 345109, 326479, 337999, 321163, 949423, 904369, 811297, 319129}:
            continue
        fresh_curves.append((p_c, n_c, b_c, lam_c))

    print(f"  Generated {len(fresh_curves)} fresh curves.")
    print()

    m = 12
    K1_TEST = 72  # The Exp L discrimination point
    n_seeds_quick = 6

    print(f"  {'#':>3}  {'p':>8}  {'n':>8}  {'λ/n':>7}  {'q_cf_small':>11}  "
          f"{'predict':>8}  {'actual':>8}  {'correct?':>9}")
    print(f"  {'-'*3}  {'-'*8}  {'-'*8}  {'-'*7}  {'-'*11}  {'-'*8}  {'-'*8}  {'-'*9}")

    n_correct = 0
    n_total = 0
    results_r = []

    for i, (p_c, n_c, b_c, lam_c) in enumerate(fresh_curves):
        G_c = find_generator(p_c, b_c, n_c)
        if G_c is None:
            continue

        # q_cf prediction
        p_cf, q_cf, _ = first_convergent_in_range(lam_c, n_c, lo=20, hi=100)
        predict_c1 = (q_cf is not None)  # q_cf < 100 → predict C1 (will fail)

        # Run LLL at K1=72 to get actual class
        recovered_count = 0
        for seed in range(n_seeds_quick):
            ok, _, _, _, _ = lll_attack(p_c, b_c, n_c, lam_c, G_c, m, K1_TEST, seed)
            if ok:
                recovered_count += 1
        actual_c2 = (recovered_count >= n_seeds_quick // 2)  # C2: succeeds at K1=72
        actual_c1 = not actual_c2

        # Did prediction match?
        match = (predict_c1 == actual_c1)
        n_correct += int(match)
        n_total += 1

        q_str = str(q_cf) if q_cf is not None else "none"
        pred_str = "C1" if predict_c1 else "C2"
        act_str = "C1" if actual_c1 else "C2"
        ok_str = "OK" if match else "MISS"

        print(f"  {i+1:>3}  {p_c:>8}  {n_c:>8}  {lam_c/n_c:>7.4f}  {q_str:>11}  "
              f"{pred_str:>8}  {act_str:>8}  {ok_str:>9}  "
              f"({recovered_count}/{n_seeds_quick} rec)")

        results_r.append({
            'p': p_c, 'n': n_c, 'b': b_c, 'lam': lam_c, 'q_cf': q_cf,
            'predict_c1': predict_c1, 'actual_c1': actual_c1, 'match': match,
            'recovered': recovered_count
        })

    print()
    if n_total > 0:
        accuracy = n_correct / n_total * 100
        print(f"  PREDICTION ACCURACY: {n_correct}/{n_total} = {accuracy:.1f}%")

        # Breakdown
        c1_actual = [r for r in results_r if r['actual_c1']]
        c2_actual = [r for r in results_r if not r['actual_c1']]
        c1_pred_c1 = [r for r in c1_actual if r['predict_c1']]
        c2_pred_c2 = [r for r in c2_actual if not r['predict_c1']]
        print(f"  C1 sensitivity (C1→predicted C1): {len(c1_pred_c1)}/{len(c1_actual)}")
        print(f"  C2 specificity (C2→predicted C2): {len(c2_pred_c2)}/{len(c2_actual)}")

        # Find misses
        misses = [r for r in results_r if not r['match']]
        if misses:
            print(f"\n  MISSES:")
            for r in misses:
                pred = "C1" if r['predict_c1'] else "C2"
                act = "C1" if r['actual_c1'] else "C2"
                print(f"    n={r['n']}, λ/n={r['lam']/r['n']:.4f}, "
                      f"q_cf={r['q_cf']}, predicted={pred}, actual={act}")

    return results_r


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time.time()
    print()
    print("GLV-HNP Phase 2: Exp P/Q/R — CF-convergent separator analysis")
    print("  Date: 2026-06-28")
    print()

    results_p = exp_p_cf_analysis()

    # Check: how many C1 curves have q_cf in [20,100] vs C2?
    c1_res = [r for r in results_p if r['cls'] == 1]
    c2_res = [r for r in results_p if r['cls'] == 2]
    c1_has_qcf = [r for r in c1_res if r['q_cf'] is not None]
    c2_has_qcf = [r for r in c2_res if r['q_cf'] is not None]

    print()
    print("  CF HYPOTHESIS PRELIMINARY CHECK (known curves):")
    print(f"    C1-type with q_cf in [20,100]: {len(c1_has_qcf)}/{len(c1_res)}")
    print(f"    C2-type with q_cf in [20,100]: {len(c2_has_qcf)}/{len(c2_res)}")
    if c1_has_qcf:
        q_vals = [r['q_cf'] for r in c1_has_qcf]
        k1_vals = [r['k1_thresh'] for r in c1_has_qcf]
        diffs = [abs(q - k) for q, k in zip(q_vals, k1_vals)]
        print(f"    C1: q_cf values = {q_vals}, K1_thresh = {k1_vals}")
        print(f"    C1: |q_cf - K1_thresh| = {diffs}")
    print()

    exp_q_spurious_anatomy()
    results_r = exp_r_prediction_accuracy()

    elapsed = time.time() - t0
    print()
    print(f"Total wall time: {elapsed:.1f}s")
    print()
    print("=" * 70)
    print("SUMMARY (Exp P/Q/R)")
    print("=" * 70)
    print()
    print("  See findings above for full tables.")
    print()
    print("  Key question: does q_cf (first CF convergent of λ/n in [20,100])")
    print("  separate C1 from C2 curves, and does q_cf ≈ K1_threshold for C1?")
