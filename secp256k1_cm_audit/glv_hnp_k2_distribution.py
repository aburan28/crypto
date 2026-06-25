"""
GLV-HNP Phase 2: Exp G/H/I — spurious vector anatomy + m-scaling for C1.

2026-06-25:

  Exp G — Spurious vector anatomy (C1 at K1=72, m=12, 50 seeds):
    When LLL returns wrong d (spurious Kannan row), what are the actual
    k1_i and k2_i encoded in that row?
    Hypothesis: spurious rows have genuinely small k1_i (within [0,K1-1])
    but WRONG d, meaning they satisfy a different ECDSA-like relation.
    Test: for each seed, extract both the CORRECT row (if present) and the
    SHORTEST Kannan row. Compare their k1_i/k2_i distributions.

  Exp H — Vary m at K1=72 for C1:
    m ∈ {8, 10, 12, 16, 20, 24}, 10 seeds each.
    Does recovery rate improve with more equations?
    If yes: C1 is solvable with more signatures.
    If no: structural obstruction — LLL cannot compress C1's solution
           regardless of lattice height.

  Exp I — k1-block comparison (C1 vs C2):
    For both curves at m=12, K1=72, 100 seeds:
    - Fraction of k1_i that fall in [0, K1-1] = [0,71]
    - Average k1_i for the SHORTEST Kannan row vs correct row
    - Is C1's shortest row achieved with k1_i *outside* [0,71]?
    This directly tests whether spurious solutions exploit larger k1-block.
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (copied from glv_hnp_bkz_c1.py)
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
# Lattice construction
# ---------------------------------------------------------------------------

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
    return M, S_K1, S_K2, S_KANNAN

def gen_signatures(p, b_param, n, lam, G, m_sigs, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    k1_vals, k2_vals = [], []
    sigs = []
    attempts = 0
    while len(sigs) < m_sigs and attempts < 200000:
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
        k1_vals.append(k1)
        k2_vals.append(k2)
    return d_secret, sigs, k1_vals, k2_vals

def nominal_target_norm(k1_vals, k2_vals, d, S_K1, S_K2, S_KANNAN):
    norm_sq = sum(k1 * k1 * S_K1 * S_K1 for k1 in k1_vals)
    norm_sq += d * d
    norm_sq += sum(k2 * k2 * S_K2 * S_K2 for k2 in k2_vals)
    norm_sq += S_KANNAN * S_KANNAN
    return math.sqrt(norm_sq)

def extract_kannan_rows(M_red, m, n, d_secret, S_KANNAN, S_K1, S_K2, K1, K2):
    """
    For each Kannan row (last entry == ±S_KANNAN), decode:
      - k1_i = row[i] / S_K1  (should be in [0, K1-1] for correct soln)
      - d_cand = row[m] * sign
      - k2_i = row[m+1+i] / S_K2  (should be in [0, K2-1] for correct soln)
    Returns a list of dicts.
    """
    dim = 2 * m + 2
    rows = []
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n

        k1_decoded = [sign * row[i] / S_K1 for i in range(m)]
        k2_decoded = [sign * row[m + 1 + i] / S_K2 for i in range(m)]

        k1_valid = all(0 <= k for k in k1_decoded) and all(k < K1 for k in k1_decoded)
        k2_valid = all(0 <= k for k in k2_decoded) and all(k < K2 for k in k2_decoded)

        is_correct = (d_cand == d_secret)

        row_norm = math.sqrt(sum(x * x for x in row))

        rows.append({
            'is_correct': is_correct,
            'd_cand': d_cand,
            'k1_decoded': k1_decoded,
            'k2_decoded': k2_decoded,
            'k1_valid': k1_valid,
            'k2_valid': k2_valid,
            'row_norm': row_norm,
        })

    rows.sort(key=lambda r: r['row_norm'])
    return rows


# ---------------------------------------------------------------------------
# Exp G: Spurious vector anatomy
# ---------------------------------------------------------------------------

def exp_g_spurious_anatomy():
    print("=" * 70)
    print("Exp G: Spurious vector anatomy (C1, K1=72, m=12, 50 seeds)")
    print("  Q: When LLL fails for C1, what does the returned (shortest)")
    print("     Kannan row actually look like? Are k1_i/k2_i in-range?")
    print("=" * 70)

    K1 = 72
    m = 12
    N_SEEDS = 50
    p, n = 524743, 523597
    label = "C1-FAIL"

    lam = glv_eigenvalue(n)
    b_param = find_b_for_n(p, n)
    G = find_generator(p, b_param, n)
    K2 = math.isqrt(n) + 1
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    dim = 2 * m + 2

    print(f"  {label}: p={p}, n={n}, K1={K1}, K2={K2}, lam/n={lam/n:.4f}")
    print(f"  S_K1={S_K1}, S_K2={S_K2}, S_KANNAN={S_KANNAN}\n")

    correct_recovered = 0
    spuri_k1_valid_count = 0
    spuri_k2_valid_count = 0
    spuri_both_valid_count = 0
    total_spuri = 0
    correct_k1_avgs = []
    correct_k2_avgs = []
    spuri_k1_avgs = []
    spuri_k2_avgs = []
    spuri_norm_ratios = []
    correct_norm_ratios = []

    for seed_i in range(N_SEEDS):
        d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed_i)
        if len(sigs) < m:
            continue

        nom_norm = nominal_target_norm(k1v, k2v, d_sec, S_K1, S_K2, S_KANNAN)

        M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
        A = IntegerMatrix.from_matrix(M_int)
        LLL.reduction(A)
        M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

        rows = extract_kannan_rows(M_red, m, n, d_sec, S_KANNAN, S_K1, S_K2, K1, K2)
        if not rows:
            continue

        shortest = rows[0]
        for r in rows:
            if r['is_correct']:
                correct_recovered += 1
                correct_norm_ratios.append(r['row_norm'] / nom_norm)
                correct_k1_avgs.append(sum(r['k1_decoded']) / m)
                correct_k2_avgs.append(sum(r['k2_decoded']) / m)
                break

        if not shortest['is_correct']:
            total_spuri += 1
            spuri_norm_ratios.append(shortest['row_norm'] / nom_norm)
            if shortest['k1_valid']:
                spuri_k1_valid_count += 1
            if shortest['k2_valid']:
                spuri_k2_valid_count += 1
            if shortest['k1_valid'] and shortest['k2_valid']:
                spuri_both_valid_count += 1
            spuri_k1_avgs.append(sum(shortest['k1_decoded']) / m)
            spuri_k2_avgs.append(sum(shortest['k2_decoded']) / m)

    print(f"  Recovery: {correct_recovered}/{N_SEEDS}")
    print(f"  Seeds where shortest Kannan = SPURIOUS: {total_spuri}/{N_SEEDS}")
    print()

    if total_spuri > 0:
        print("  SPURIOUS shortest-row analysis (when LLL returns wrong d):")
        print(f"  {'Metric':<40} {'Value':>10}")
        print("  " + "-" * 52)
        print(f"  {'k1_valid fraction (all k1 in [0,K1-1])':<40} {spuri_k1_valid_count/total_spuri:>10.3f}")
        print(f"  {'k2_valid fraction (all k2 in [0,K2-1])':<40} {spuri_k2_valid_count/total_spuri:>10.3f}")
        print(f"  {'both_valid fraction':<40} {spuri_both_valid_count/total_spuri:>10.3f}")
        if spuri_k1_avgs:
            print(f"  {'avg k1_i (spurious row)':<40} {sum(spuri_k1_avgs)/len(spuri_k1_avgs):>10.2f}")
            print(f"  {'avg k2_i (spurious row)':<40} {sum(spuri_k2_avgs)/len(spuri_k2_avgs):>10.2f}")
        if spuri_norm_ratios:
            spuri_norm_ratios.sort()
            print(f"  {'norm ratio min (spuri)':<40} {min(spuri_norm_ratios):>10.4f}")
            print(f"  {'norm ratio median (spuri)':<40} {spuri_norm_ratios[len(spuri_norm_ratios)//2]:>10.4f}")
            print(f"  {'norm ratio max (spuri)':<40} {max(spuri_norm_ratios):>10.4f}")

    if correct_recovered > 0:
        print()
        print("  CORRECT row analysis (when found by LLL):")
        print(f"  {'avg k1_i (correct row)':<40} {sum(correct_k1_avgs)/len(correct_k1_avgs):>10.2f}")
        print(f"  {'avg k2_i (correct row)':<40} {sum(correct_k2_avgs)/len(correct_k2_avgs):>10.2f}")
        correct_norm_ratios.sort()
        print(f"  {'norm ratio median (correct)':<40} {correct_norm_ratios[len(correct_norm_ratios)//2]:>10.4f}")

    print(f"\n  For reference: K1={K1}, K2={K2}")
    print(f"  Expected avg k1_i ~ {(K1-1)/2:.1f} (uniform [0,K1-1])")
    print(f"  Expected avg k2_i ~ {(K2-1)/2:.1f} (uniform [0,K2-1])")


# ---------------------------------------------------------------------------
# Exp H: Vary m for C1 and C2
# ---------------------------------------------------------------------------

def exp_h_vary_m():
    print("\n" + "=" * 70)
    print("Exp H: Vary m ∈ {8,10,12,16,20,24} at K1=72 (C1 and C2, 10 seeds)")
    print("  Q: Does more equations help C1 recover?")
    print("=" * 70)

    K1 = 72
    M_VALUES = [8, 10, 12, 16, 20, 24]
    N_SEEDS = 10

    curves = [
        ("C1-FAIL",    524743, 523597),
        ("C2-SUCCEED", 525043, 524269),
    ]

    for (label, p, n) in curves:
        lam = glv_eigenvalue(n)
        b_param = find_b_for_n(p, n)
        G = find_generator(p, b_param, n)
        K2 = math.isqrt(n) + 1

        print(f"\n  {label}: p={p}, n={n}, λ/n={lam/n:.4f}")
        print(f"  {'m':>4} {'dim':>5} {'recovery':>10} {'med_ratio':>11} {'t(s)':>7}")
        print("  " + "-" * 40)

        for m in M_VALUES:
            S_K1 = n // K1
            S_K2 = max(1, n // K2)
            S_KANNAN = n
            dim = 2 * m + 2

            recover_count = 0
            ratios = []
            t0 = time.time()

            for seed_i in range(N_SEEDS):
                d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed_i)
                if len(sigs) < m:
                    continue

                nom_norm = nominal_target_norm(k1v, k2v, d_sec, S_K1, S_K2, S_KANNAN)
                M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
                A = IntegerMatrix.from_matrix(M_int)
                LLL.reduction(A)
                M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

                for row in M_red:
                    last = row[dim - 1]
                    if abs(last) == S_KANNAN:
                        sign = 1 if last > 0 else -1
                        d_cand = (sign * row[m]) % n
                        if d_cand == d_sec:
                            recover_count += 1
                            rn = math.sqrt(sum(x * x for x in row))
                            ratios.append(rn / nom_norm)
                        break

            t_elapsed = time.time() - t0
            med = sorted(ratios)[len(ratios)//2] if ratios else float('nan')
            print(f"  {m:>4} {dim:>5} {recover_count:>4}/{N_SEEDS}  {med:>11.4f} {t_elapsed:>7.1f}")


# ---------------------------------------------------------------------------
# Exp I: k1/k2 component breakdown for C1 vs C2 (100 seeds)
# ---------------------------------------------------------------------------

def exp_i_component_breakdown():
    print("\n" + "=" * 70)
    print("Exp I: Solution-vector component breakdown (C1 vs C2, m=12, K1=72, 100 seeds)")
    print("  Decompose nominal norm into: k1-block, d, k2-block, Kannan.")
    print("  Q: Is C1's k2-block fraction larger → harder to compress?")
    print("=" * 70)

    K1 = 72
    m = 12
    N_SEEDS = 100

    curves = [
        ("C1-FAIL",    524743, 523597),
        ("C2-SUCCEED", 525043, 524269),
    ]

    for (label, p, n) in curves:
        lam = glv_eigenvalue(n)
        b_param = find_b_for_n(p, n)
        G = find_generator(p, b_param, n)
        K2 = math.isqrt(n) + 1
        S_K1 = n // K1
        S_K2 = max(1, n // K2)
        S_KANNAN = n

        k1_frac_list = []
        k2_frac_list = []
        d_frac_list = []
        kannan_frac_list = []
        nom_norm_list = []

        for seed_i in range(N_SEEDS):
            d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed_i)
            if len(sigs) < m:
                continue
            k1_sq = sum(k * k * S_K1 * S_K1 for k in k1v)
            k2_sq = sum(k * k * S_K2 * S_K2 for k in k2v)
            d_sq = d_sec * d_sec
            kan_sq = S_KANNAN * S_KANNAN
            total_sq = k1_sq + k2_sq + d_sq + kan_sq
            if total_sq == 0:
                continue
            k1_frac_list.append(k1_sq / total_sq)
            k2_frac_list.append(k2_sq / total_sq)
            d_frac_list.append(d_sq / total_sq)
            kannan_frac_list.append(kan_sq / total_sq)
            nom_norm_list.append(math.sqrt(total_sq))

        def avg(lst): return sum(lst) / len(lst) if lst else float('nan')
        def med(lst):
            if not lst: return float('nan')
            s = sorted(lst)
            return s[len(s)//2]

        print(f"\n  {label}: p={p}, n={n}, K1={K1}, K2={K2}")
        print(f"  S_K1={S_K1}, S_K2={S_K2}, S_KANNAN={S_KANNAN}")
        print(f"  {'Component':<25} {'avg frac':>10} {'med frac':>10}")
        print("  " + "-" * 47)
        print(f"  {'k1-block (k1_i*S_K1)^2':<25} {avg(k1_frac_list):>10.4f} {med(k1_frac_list):>10.4f}")
        print(f"  {'d^2':<25} {avg(d_frac_list):>10.4f} {med(d_frac_list):>10.4f}")
        print(f"  {'k2-block (k2_i*S_K2)^2':<25} {avg(k2_frac_list):>10.4f} {med(k2_frac_list):>10.4f}")
        print(f"  {'S_KANNAN^2':<25} {avg(kannan_frac_list):>10.4f} {med(kannan_frac_list):>10.4f}")
        print(f"  {'nom_norm (avg)':<25} {avg(nom_norm_list):>10.0f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=== GLV-HNP Phase 2: Spurious anatomy + m-scaling + component breakdown ===")
    print("Date: 2026-06-25\n")

    exp_g_spurious_anatomy()
    exp_h_vary_m()
    exp_i_component_breakdown()

    print("\nDone.")
