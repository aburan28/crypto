"""
GLV-HNP Phase 2: Exp M/N/O — K1 threshold mapping (2026-06-27)

  Exp M — Fine K1 scan on C1 (p=524743, n=523597):
    From Exp J (2026-06-26): K1=32 gives 10/10, K1=48 gives 4/10, K1=64 gives 1/10.
    Fine scan: K1 in {33,35,37,40,42,44,46,50,55,60,70} × 20 seeds
    Track: recovery rate, spurious-vec norm ratio (found_norm/nom_norm).
    Goal: pin down K1_threshold(C1) = K1 where recovery crosses 50%.

  Exp N — K1 threshold sweep on 5 Exp-L C1-type curves (diverse lam/n):
    Selected curves (all 0/10 at K1=72 in Exp L):
      #2:  p=343963, n=345109, b=2,  lam/n≈0.034  (very low)
      #9:  p=327553, n=326479, b=10, lam/n≈0.097  (low)
      #15: p=339139, n=337999, b=2,  lam/n≈0.237  (medium)
      #5:  p=321889, n=321163, b=19, lam/n≈0.288  (medium-high)
      #25: p=951001, n=949423, b=11, lam/n≈0.404  (high but C1-type!)
    K1 sweep: {8,12,16,24,32,40,48,56,64,72,96} × 6 seeds
    Goal: find K1_threshold for each curve, correlate with lam/n.

  Exp O — C2-type upper threshold (3 curves, K1 escalation):
    Selected curves (all 10/10 at K1=72 in Exp L):
      #8:  p=903367,  n=904369,  lam/n≈0.086
      #14: p=812233,  n=811297,  lam/n≈0.131
      #11: p=318001,  n=319129,  lam/n≈0.480
    K1 sweep: {72, 96, 128, 160, 200, 256} × 5 seeds
    Goal: find where C2-type curves start to fail.

Run: python3 glv_hnp_k1_threshold.py
"""

import math
import random
import time
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
# Lattice construction + attack (same scaling as glv_hnp_structural_source.py)
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

def lll_attack(p, b_param, n, lam, G, m, K1, seed):
    """Returns (recovered, d_cand, nom_norm, found_norm, d_sec)."""
    d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
    if len(sigs) < m:
        return False, None, None, None, d_sec

    K2 = math.isqrt(n) + 1
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    dim = 2 * m + 2

    nom_norm_sq = (sum(k * k * S_K1 * S_K1 for k in k1v) +
                   d_sec * d_sec +
                   sum(k * k * S_K2 * S_K2 for k in k2v) +
                   S_KANNAN * S_KANNAN)
    nom_norm = math.sqrt(nom_norm_sq)

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

    best_d = None; best_norm = float('inf')
    for row in M_red:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        rn = math.sqrt(sum(x * x for x in row))
        if rn < best_norm:
            best_norm = rn
            best_d = d_cand

    if best_d is None:
        return False, None, nom_norm, None, d_sec

    return best_d == d_sec, best_d, nom_norm, best_norm, d_sec


# ---------------------------------------------------------------------------
# Exp M: Fine K1 scan on C1
# ---------------------------------------------------------------------------

def exp_m_c1_fine_scan():
    print("=" * 70)
    print("Exp M: Fine K1 scan on C1 (p=524743, n=523597)")
    print("  From Exp J: K1=32→10/10, K1=48→4/10, K1=64→1/10")
    print("  Fine scan to pin down 50% crossing and norm-ratio behaviour")
    print("=" * 70)

    p, n = 524743, 523597
    b_param = None
    # Find b
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
    print(f"  b={b_param}, λ={lam}, λ/n={lam/n:.4f}")
    print()

    m = 12
    n_seeds = 20
    K1_VALUES = [33, 35, 37, 40, 42, 44, 46, 50, 55, 60, 70]

    print(f"  {'K1':>4}  {'K1/n':>8}  {'recov':>8}  {'avg_ratio':>10}  {'min_ratio':>10}")
    print(f"  {'-'*4}  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}")

    results_m = []
    for K1 in K1_VALUES:
        recovered = 0
        ratios = []
        for seed in range(n_seeds):
            ok, d_cand, nom, found, d_sec = lll_attack(p, b_param, n, lam, G, m, K1, seed)
            if ok:
                recovered += 1
            if nom is not None and found is not None:
                ratios.append(found / nom)
        avg_ratio = sum(ratios) / len(ratios) if ratios else float('nan')
        min_ratio = min(ratios) if ratios else float('nan')
        print(f"  {K1:>4}  {K1/n:>8.5f}  {recovered:>3}/{n_seeds}  {avg_ratio:>10.4f}  {min_ratio:>10.4f}")
        results_m.append((K1, recovered, avg_ratio, min_ratio))

    # Find 50% crossing
    threshold_50 = None
    for i, (K1, rec, _, _) in enumerate(results_m):
        if rec < n_seeds / 2:
            threshold_50 = K1
            break
    print()
    print(f"  K1 threshold (first K1 with <50% recovery): {threshold_50}")
    print()
    return results_m, threshold_50


# ---------------------------------------------------------------------------
# Exp N: K1 threshold sweep on 5 Exp-L C1-type curves
# ---------------------------------------------------------------------------

def exp_n_expl_curves_k1_sweep():
    print("=" * 70)
    print("Exp N: K1 threshold sweep on 5 Exp-L C1-type curves")
    print("  Curves: diverse lam/n (0.034 to 0.404), all 0/10 at K1=72 in Exp L")
    print("=" * 70)
    print()

    curves = [
        # (label, p, n, b, expected_lam_n)
        ("#2",  343963, 345109,  2, 0.034),
        ("#9",  327553, 326479, 10, 0.097),
        ("#15", 339139, 337999,  2, 0.237),
        ("#5",  321889, 321163, 19, 0.288),
        ("#25", 951001, 949423, 11, 0.404),
    ]

    K1_VALUES = [8, 12, 16, 24, 32, 40, 48, 56, 64, 72, 96]
    m = 12
    n_seeds = 6

    all_thresholds = []

    for (label, p, n, b_param, _) in curves:
        lam = glv_eigenvalue(n)
        if lam is None:
            print(f"  {label}: no GLV eigenvalue, skip")
            continue
        G = find_generator(p, b_param, n)
        if G is None:
            print(f"  {label}: generator not found, skip")
            continue

        print(f"  Curve {label}: p={p}, n={n}, b={b_param}, λ/n={lam/n:.4f}")
        print(f"  {'K1':>5}  {'K1/n':>8}  {'recov':>8}  {'avg_ratio':>10}")

        threshold_50 = None
        for K1 in K1_VALUES:
            recovered = 0
            ratios = []
            for seed in range(n_seeds):
                ok, _, nom, found, _ = lll_attack(p, b_param, n, lam, G, m, K1, seed)
                if ok:
                    recovered += 1
                if nom is not None and found is not None:
                    ratios.append(found / nom)
            avg_ratio = sum(ratios) / len(ratios) if ratios else float('nan')
            print(f"  {K1:>5}  {K1/n:>8.5f}  {recovered:>3}/{n_seeds}  {avg_ratio:>10.4f}")
            if threshold_50 is None and recovered < n_seeds / 2:
                threshold_50 = K1

        print(f"  → K1_threshold(Curve {label}) ≈ {threshold_50}")
        print()
        all_thresholds.append((label, lam / n, threshold_50))

    print("Summary: lam/n vs K1 threshold")
    print(f"  {'Curve':>6}  {'lam/n':>8}  {'K1_thresh':>10}")
    for (label, lamn, thr) in all_thresholds:
        print(f"  {label:>6}  {lamn:>8.4f}  {str(thr):>10}")
    print()
    return all_thresholds


# ---------------------------------------------------------------------------
# Exp O: C2-type upper threshold
# ---------------------------------------------------------------------------

def exp_o_c2_upper_threshold():
    print("=" * 70)
    print("Exp O: C2-type upper threshold (K1 escalation)")
    print("  Curves: all 10/10 at K1=72 in Exp L")
    print("=" * 70)
    print()

    curves = [
        ("#8",  903367, 904369,  3, 0.086),
        ("#14", 812233, 811297, 15, 0.131),
        ("#11", 318001, 319129, 13, 0.480),
    ]

    K1_VALUES = [72, 96, 128, 160, 200, 256, 320]
    m = 12
    n_seeds = 6

    for (label, p, n, b_param, _) in curves:
        lam = glv_eigenvalue(n)
        if lam is None:
            print(f"  {label}: no GLV eigenvalue, skip")
            continue
        G = find_generator(p, b_param, n)
        if G is None:
            print(f"  {label}: generator not found, skip")
            continue

        print(f"  Curve {label}: p={p}, n={n}, b={b_param}, λ/n={lam/n:.4f}")
        print(f"  {'K1':>5}  {'K1/n':>8}  {'recov':>8}  {'avg_ratio':>10}")

        threshold_fail = None
        for K1 in K1_VALUES:
            recovered = 0
            ratios = []
            for seed in range(n_seeds):
                ok, _, nom, found, _ = lll_attack(p, b_param, n, lam, G, m, K1, seed)
                if ok:
                    recovered += 1
                if nom is not None and found is not None:
                    ratios.append(found / nom)
            avg_ratio = sum(ratios) / len(ratios) if ratios else float('nan')
            print(f"  {K1:>5}  {K1/n:>8.5f}  {recovered:>3}/{n_seeds}  {avg_ratio:>10.4f}")
            if threshold_fail is None and recovered < n_seeds / 2:
                threshold_fail = K1

        thr_str = str(threshold_fail) if threshold_fail else ">320"
        print(f"  → C2 starts failing at K1 ≈ {thr_str}")
        print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time.time()
    print()
    print("GLV-HNP Phase 2: Exp M/N/O — K1 threshold mapping")
    print(f"  Date: 2026-06-27")
    print()

    results_m, thr_c1 = exp_m_c1_fine_scan()
    thresholds_n = exp_n_expl_curves_k1_sweep()
    exp_o_c2_upper_threshold()

    elapsed = time.time() - t0
    print(f"Total wall time: {elapsed:.1f}s")

    print()
    print("=" * 70)
    print("KEY FINDINGS (Exp M/N/O)")
    print("=" * 70)
    print(f"  C1 (p=524743, λ/n=0.211): K1_threshold = {thr_c1}")
    print("  Exp-L C1-type thresholds vs lam/n:")
    for (label, lamn, thr) in thresholds_n:
        print(f"    Curve {label}: λ/n={lamn:.4f}  K1_thresh={thr}")
    print()
    print("  Interpretation:")
    print("  If threshold(C1-type) ~ constant regardless of lam/n → lam/n is NOT the predictor")
    print("  If threshold increases with lam/n → lam/n IS the predictor")
