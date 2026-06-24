"""
GLV-HNP Phase 2: BKZ-on-C1 and 50-seed statistical analysis.

2026-06-24: Yesterday (Exp A/B/C) showed:
  - LLL finds Kannan rows for BOTH C1 and C2 (volume obstruction ruled out)
  - C2's correct solution compresses to ~65-75% of nominal norm → LLL finds it
  - C1's correct solution does NOT compress → LLL returns spurious rows instead
  - Open Q: WHY does C2 compress but C1 doesn't?

Today's experiments:
  Exp D — BKZ escalation on C1:
    Run BKZ with β=15,20,25,30 on C1 (failing) at m=12.
    If BKZ finds the correct solution, it confirms:
      (a) the solution IS present in the lattice at some short vector
      (b) it's just too hard for LLL to find without extra sieve help
    If BKZ also fails, suggests the effective norm for C1 is genuinely above
    the GH threshold even in LLL-reduced basis.

  Exp E — 50-seed statistical distribution:
    For C1 and C2, over 50 seeds:
      - Measure "effective norm ratio" = (norm of correct solution in LLL basis) / (nominal target norm)
      - For C2 this should cluster below 1.0 (compression); for C1 above 1.0 (no compression)
      - This is the key discriminating statistic

Curves (hardcoded from 2026-06-23 search):
  C1 (FAIL): p=524743, n=523597, lam/n≈0.2114, δ/n≈0.366
  C2 (SUCCEED): p=525043, n=524269, lam/n≈0.2122, δ/n≈0.364
"""

import math
import random
import time
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

# ---------------------------------------------------------------------------
# Signature generation
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Key-recovery check on a reduced basis
# ---------------------------------------------------------------------------

def check_recovery(M_red, m, n, d_secret, S_KANNAN):
    dim = 2 * m + 2
    for row in M_red:
        last = row[dim - 1]
        if abs(last) == S_KANNAN:
            sign = 1 if last > 0 else -1
            d_cand = (sign * row[m]) % n
            if d_cand == d_secret:
                return True
    return False

def effective_norm_of_solution(M_red, m, n, d_secret, S_KANNAN, S_K1, S_K2):
    """
    Find the correct-solution row (if present) in the reduced basis.
    Return (found, row_norm, row_index, norm_of_shortest_kannan_row).
    We need to FIND the row by scanning the Kannan rows.
    """
    dim = 2 * m + 2
    correct_norm = None
    correct_idx = None
    shortest_kannan_norm = None
    for row_idx, row in enumerate(M_red):
        last = row[dim - 1]
        if abs(last) == S_KANNAN:
            rn = math.sqrt(sum(x * x for x in row))
            if shortest_kannan_norm is None or rn < shortest_kannan_norm:
                shortest_kannan_norm = rn
            sign = 1 if last > 0 else -1
            d_cand = (sign * row[m]) % n
            if d_cand == d_secret:
                correct_norm = rn
                correct_idx = row_idx
    found = correct_norm is not None
    return found, correct_norm, correct_idx, shortest_kannan_norm

def nominal_target_norm(k1_vals, k2_vals, d, S_K1, S_K2, S_KANNAN):
    norm_sq = sum(k1 * k1 * S_K1 * S_K1 for k1 in k1_vals)
    norm_sq += d * d
    norm_sq += sum(k2 * k2 * S_K2 * S_K2 for k2 in k2_vals)
    norm_sq += S_KANNAN * S_KANNAN
    return math.sqrt(norm_sq)

# ---------------------------------------------------------------------------
# Exp D: BKZ escalation on C1
# ---------------------------------------------------------------------------

def exp_d_bkz_escalation():
    print("=" * 70)
    print("Exp D: BKZ escalation on C1 (failing curve) vs C2 (succeeding)")
    print("  K1=72, m=12, seeds=[0xDEAD, 0xBEEF, 0xCAFE]")
    print("  BKZ blocks β ∈ {15, 20, 25, 30}")
    print("=" * 70)

    K1 = 72
    m = 12
    SEEDS = [0xDEAD, 0xBEEF, 0xCAFE]
    BKZ_BETAS = [15, 20, 25, 30]

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

        lam_n = lam / n
        dim = 2 * m + 2

        log2_vol = (m * math.log2(n * S_K1)
                    + m * math.log2(S_K2)
                    + math.log2(S_KANNAN))
        gh_est = math.sqrt(dim / (2 * math.pi * math.e)) * 2 ** (log2_vol / dim)

        print(f"\n--- {label}: p={p}, n={n}, λ/n={lam_n:.4f} ---")
        print(f"{'Seed':<10} {'NomNorm':>11} {'GH':>11} {'v/GH':>6} "
              f"{'LLL':>5} {'BKZ15':>6} {'BKZ20':>6} {'BKZ25':>6} {'BKZ30':>6}")
        print("-" * 72)

        for seed in SEEDS:
            d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
            if len(sigs) < m:
                print(f"  seed={hex(seed)}: only {len(sigs)} sigs, skipping")
                continue

            nom_norm = nominal_target_norm(k1v, k2v, d_sec, S_K1, S_K2, S_KANNAN)

            M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
            A = IntegerMatrix.from_matrix(M_int)

            # LLL baseline
            t0 = time.time()
            LLL.reduction(A)
            t_lll = time.time() - t0
            M_lll = [[A[i][j] for j in range(dim)] for i in range(dim)]
            lll_ok = check_recovery(M_lll, m, n, d_sec, S_KANNAN)

            results = {'LLL': lll_ok}

            # BKZ escalation
            for beta in BKZ_BETAS:
                # Re-run from scratch (LLL pre-processes, then BKZ)
                A2 = IntegerMatrix.from_matrix(M_int)
                LLL.reduction(A2)
                t0 = time.time()
                BKZ.reduction(A2, BKZ.Param(beta))
                t_bkz = time.time() - t0
                M_bkz = [[A2[i][j] for j in range(dim)] for i in range(dim)]
                bkz_ok = check_recovery(M_bkz, m, n, d_sec, S_KANNAN)
                results[f'BKZ{beta}'] = bkz_ok

            row = [
                f"{hex(seed):<10}",
                f"{nom_norm:11.0f}",
                f"{gh_est:11.0f}",
                f"{nom_norm/gh_est:6.3f}",
                f"{'Y' if results['LLL'] else 'N':>5}",
            ]
            for beta in BKZ_BETAS:
                row.append(f"{'Y' if results[f'BKZ{beta}'] else 'N':>6}")
            print(" ".join(row))

# ---------------------------------------------------------------------------
# Exp E: 50-seed statistical analysis
# ---------------------------------------------------------------------------

def exp_e_statistical():
    print("\n" + "=" * 70)
    print("Exp E: 50-seed statistical analysis — effective norm ratio distribution")
    print("  K1=72, m=12, 50 seeds per curve")
    print("  Metric: (norm of shortest Kannan row in LLL-reduced basis) / (nominal target norm)")
    print("=" * 70)

    K1 = 72
    m = 12
    N_SEEDS = 50

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
        dim = 2 * m + 2

        print(f"\n--- {label}: p={p}, n={n} ---")

        ratios = []
        correct_ratios = []  # norm of CORRECT row / nominal (for seeds where correct row found)
        recover_count = 0
        no_kannan_count = 0

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

            found, correct_norm, correct_idx, shortest_kannan = effective_norm_of_solution(
                M_red, m, n, d_sec, S_KANNAN, S_K1, S_K2)

            if shortest_kannan is None:
                no_kannan_count += 1
                continue

            ratio = shortest_kannan / nom_norm
            ratios.append(ratio)

            if found:
                recover_count += 1
                correct_ratios.append(correct_norm / nom_norm)

        t_elapsed = time.time() - t0

        if not ratios:
            print("  No valid seeds!")
            continue

        ratios.sort()
        n_valid = len(ratios)

        print(f"  Valid seeds: {n_valid}/{N_SEEDS}  (time: {t_elapsed:.1f}s)")
        print(f"  Recovery rate: {recover_count}/{n_valid} ({100*recover_count/n_valid:.0f}%)")
        print(f"  No-Kannan seeds: {no_kannan_count}")
        print(f"\n  Distribution of (shortest Kannan norm) / (nominal target norm):")
        print(f"  {'Metric':<20} {'Value':>10}")
        print(f"  {'min':<20} {min(ratios):>10.4f}")
        print(f"  {'p10':<20} {ratios[max(0,n_valid//10)]:>10.4f}")
        print(f"  {'p25':<20} {ratios[n_valid//4]:>10.4f}")
        print(f"  {'median':<20} {ratios[n_valid//2]:>10.4f}")
        print(f"  {'p75':<20} {ratios[3*n_valid//4]:>10.4f}")
        print(f"  {'p90':<20} {ratios[min(9*n_valid//10, n_valid-1)]:>10.4f}")
        print(f"  {'max':<20} {max(ratios):>10.4f}")
        print(f"  {'mean':<20} {sum(ratios)/len(ratios):>10.4f}")
        print(f"  {'frac < 0.9':<20} {sum(1 for r in ratios if r < 0.9)/n_valid:>10.4f}")
        print(f"  {'frac < 1.0':<20} {sum(1 for r in ratios if r < 1.0)/n_valid:>10.4f}")
        print(f"  {'frac < 1.1':<20} {sum(1 for r in ratios if r < 1.1)/n_valid:>10.4f}")

        if correct_ratios:
            correct_ratios.sort()
            print(f"\n  Distribution of (CORRECT row norm) / (nominal) [n={len(correct_ratios)}]:")
            print(f"  {'min':<20} {min(correct_ratios):>10.4f}")
            print(f"  {'median':<20} {correct_ratios[len(correct_ratios)//2]:>10.4f}")
            print(f"  {'max':<20} {max(correct_ratios):>10.4f}")
            print(f"  {'mean':<20} {sum(correct_ratios)/len(correct_ratios):>10.4f}")

# ---------------------------------------------------------------------------
# Exp F: Effective norm vs δ/n correlation across 5 curves
# ---------------------------------------------------------------------------

def exp_f_delta_correlation():
    """
    Scan a small set of j=0 CM curves with varying δ/n and measure:
    - LLL recovery rate (3 seeds each)
    - Shortest Kannan row norm / nominal ratio (3 seeds each)
    Goal: does frac(shortest_row < 0.9 * nominal) correlate with δ/n?
    """
    print("\n" + "=" * 70)
    print("Exp F: δ/n correlation — does effective norm ratio correlate with δ/n?")
    print("  Scan j=0 CM curves (18–20 bit) at K1=72, m=12, 3 seeds each")
    print("=" * 70)

    K1 = 72
    m = 12
    SEEDS = [0xDEAD, 0xBEEF, 0xCAFE]
    TARGET_SEEDS = 3

    import sympy

    def eisenstein_decompose(p):
        for a in range(1, 2 * math.isqrt(p // 3) + 3):
            disc = 4 * p - 3 * a * a
            if disc < 0: break
            s = math.isqrt(disc)
            if s * s != disc: continue
            for num in [a + s, a - s]:
                if num % 2 == 0:
                    b2 = num // 2
                    if b2 >= 0 and a * a - a * b2 + b2 * b2 == p:
                        return (a, b2)
        return None

    def j0_traces(a, b2):
        return [2*a - b2, -2*a + b2, -(a + b2), a + b2, 2*b2 - a, a - 2*b2]

    def delta_ratio(lam, n):
        x = (3 * lam) % n
        return min(x, n - x) / n

    # collect at most 8 curves spanning different δ/n values
    found_curves = []
    p = sympy.nextprime(2**18 - 1)
    while p < 2**20 and len(found_curves) < 8:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for tr in j0_traces(a_e, b_e):
                    n_cand = p + 1 - tr
                    if n_cand < 2: continue
                    if not sympy.isprime(n_cand): continue
                    if n_cand % 3 != 1: continue
                    lam_cand = glv_eigenvalue(n_cand)
                    if lam_cand is None: continue
                    dr = delta_ratio(lam_cand, n_cand)
                    lam_n = lam_cand / n_cand
                    if lam_n < 0.08 or lam_n > 0.45: continue
                    b_p = find_b_for_n(p, n_cand)
                    if b_p is None: continue
                    found_curves.append((p, n_cand, lam_cand, lam_n, dr, b_p))
                    break
        p = sympy.nextprime(p)

    found_curves.sort(key=lambda x: x[4])  # sort by δ/n

    print(f"  Found {len(found_curves)} curves. Running LLL analysis...\n")
    print(f"  {'p':>7} {'n':>7} {'λ/n':>6} {'δ/n':>6} {'recovery':>9} {'median_ratio':>13}")
    print("  " + "-" * 55)

    for (p, n, lam, lam_n, dr, b_param) in found_curves:
        G = find_generator(p, b_param, n)
        if G is None:
            print(f"  {p:>7} {n:>7}: no generator found")
            continue

        K2 = math.isqrt(n) + 1
        S_K1 = n // K1
        S_K2 = max(1, n // K2)
        S_KANNAN = n
        dim = 2 * m + 2

        recover_count = 0
        ratios = []
        for seed in SEEDS:
            d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
            if len(sigs) < m:
                continue
            nom_norm = nominal_target_norm(k1v, k2v, d_sec, S_K1, S_K2, S_KANNAN)
            M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
            A = IntegerMatrix.from_matrix(M_int)
            LLL.reduction(A)
            M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

            found_flag, correct_norm, _, shortest_k = effective_norm_of_solution(
                M_red, m, n, d_sec, S_KANNAN, S_K1, S_K2)
            if found_flag:
                recover_count += 1
            if shortest_k is not None:
                ratios.append(shortest_k / nom_norm)

        med = sorted(ratios)[len(ratios)//2] if ratios else float('nan')
        print(f"  {p:>7} {n:>7} {lam_n:>6.3f} {dr:>6.3f}  {recover_count}/{TARGET_SEEDS}     {med:>13.4f}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=== GLV-HNP Phase 2: BKZ escalation + Statistical analysis ===")
    print("Date: 2026-06-24\n")

    exp_d_bkz_escalation()
    exp_e_statistical()
    exp_f_delta_correlation()

    print("\nDone.")
