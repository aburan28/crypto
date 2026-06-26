"""
GLV-HNP Phase 2: Exp J/K/L — K1 variation + CM-artifact check + cross-curve screening

2026-06-26:

  Exp J — K1 variation for C1 (K1 ∈ {8,16,32,48,64,96,128}, m=12, 10 seeds):
    Hypothesis: 0% recovery persists at ALL K1 → structural obstruction is
    not a K1-scaling artifact. Alternatively, small K1 might succeed.

  Exp K — Spurious-vector CM-artifact check (C1, K1=72, m=12, 20 seeds):
    For each seed: extract d_cand from the spurious Kannan row.
    If d_cand is the SAME small value across all seeds (independent of A_i,B_i),
    the spurious vector is a CM ring artifact.
    If d_cand varies, it depends on signature data.
    Also: extract the full spurious row for seed=0 and compute its norm
    fraction vs the CM-only sublattice (rows 0..m-1 and m+1..2m).

  Exp L — Cross-curve screening (j=0 CM curves, 18-20 bit range, target 30):
    Generate random j=0 CM curves with prime group order.
    For each: K1=72, m=12, 10 seeds → recovery rate.
    Classify C1-type (≤20%) vs C2-type (≥80%). Record algebraic invariants:
    p, n, b, lam, lam/n, n mod 9, (n-1)/3 mod 3, disc, conductor proxy.
    Goal: find algebraic discriminant between the two types.
"""

import math
import random
import time
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Primality (trial division up to sqrt, works for n ≤ 2^22)
# ---------------------------------------------------------------------------

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
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
    m2, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m2 - i - 1), p)
        m2, c, t, r = i, b * b % p, t * b * b % p, r * b % p

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

def random_point_on_curve(p, b, rng):
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_prime_order_curve(p, b_range=50, n_check_points=3):
    """
    Find b in [1, b_range] such that y^2 = x^3 + b over F_p has prime group order.
    Returns (b, n) or (None, None).
    """
    hasse_r = 2 * math.isqrt(p) + 2
    hasse_lo = p + 1 - hasse_r
    hasse_hi = p + 1 + hasse_r

    # Precompute primes in Hasse interval
    prime_candidates = [n for n in range(max(3, hasse_lo), hasse_hi + 1) if is_prime(n)]

    for b in range(1, b_range + 1):
        rng = random.Random(b * p + 7)
        pts = []
        for _ in range(n_check_points * 5):
            pt = random_point_on_curve(p, b, rng)
            if pt is not None:
                pts.append(pt)
            if len(pts) >= n_check_points:
                break
        if len(pts) < n_check_points:
            continue

        for n_cand in prime_candidates:
            lam_cand = glv_eigenvalue(n_cand)
            if lam_cand is None:
                continue
            if all(ec_mul(pt, n_cand, p) is None for pt in pts):
                return b, n_cand

    return None, None

# ---------------------------------------------------------------------------
# Lattice construction + attack
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
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, n // K1, max(1, n // (math.isqrt(n) + 1)), n

def gen_signatures(p, b_param, n, lam, G, m_sigs, K1, seed):
    K2 = math.isqrt(n) + 1
    rng = random.Random(seed)
    d_secret = rng.randint(1, n - 1)
    sigs, k1_vals, k2_vals = [], [], []
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
        sigs.append({'A': h * s_inv % n, 'B': r * s_inv % n})
        k1_vals.append(k1); k2_vals.append(k2)
    return d_secret, sigs, k1_vals, k2_vals

def lll_attack(p, b_param, n, lam, G, m, K1, seed):
    """Returns (recovered, d_cand, nom_norm, found_norm) for the shortest Kannan row."""
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

    M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
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
# Exp J: K1 variation for C1
# ---------------------------------------------------------------------------

def exp_j_k1_variation():
    print("=" * 70)
    print("Exp J: K1 variation for C1 (m=12, 10 seeds per K1)")
    print("  Q: Does 0% recovery persist at ALL K1?")
    print("=" * 70)

    K1_VALUES = [8, 16, 32, 48, 64, 96, 128]
    m = 12
    N_SEEDS = 10
    p, n = 524743, 523597

    lam = glv_eigenvalue(n)
    b_param = find_b_for_n(p, n)
    G = find_generator(p, b_param, n)

    print(f"  C1-FAIL: p={p}, n={n}, lam/n={lam/n:.4f}, b={b_param}\n")
    print(f"  {'K1':>5} {'K1/n':>8} {'recovery':>10} {'med_ratio':>11} {'t(s)':>7}")
    print("  " + "-" * 45)

    for K1 in K1_VALUES:
        recover = 0
        ratios = []
        t0 = time.time()
        for seed in range(N_SEEDS):
            ok, d_cand, nom, found, d_sec = lll_attack(p, b_param, n, lam, G, m, K1, seed)
            if ok:
                recover += 1
                if nom and found:
                    ratios.append(found / nom)
        t_el = time.time() - t0
        med = sorted(ratios)[len(ratios)//2] if ratios else float('nan')
        print(f"  {K1:>5} {K1/n:>8.5f} {recover:>4}/{N_SEEDS}  {med:>11.4f} {t_el:>7.2f}")

    print()

# ---------------------------------------------------------------------------
# Exp K: Spurious-vector CM-artifact check
# ---------------------------------------------------------------------------

def exp_k_spurious_cm_check():
    print("=" * 70)
    print("Exp K: Spurious-vector CM-artifact check (C1, K1=72, m=12, 20 seeds)")
    print("  Q: Is d_cand in the spurious row the same across seeds?")
    print("     If yes → CM artifact. If varies → signature-dependent.")
    print("=" * 70)

    K1 = 72
    m = 12
    N_SEEDS = 20
    p, n = 524743, 523597

    lam = glv_eigenvalue(n)
    b_param = find_b_for_n(p, n)
    G = find_generator(p, b_param, n)
    K2 = math.isqrt(n) + 1
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    dim = 2 * m + 2

    print(f"  C1: p={p}, n={n}, K1={K1}, K2={K2}, lam/n={lam/n:.4f}\n")

    d_cand_list = []
    d_cand_small_list = []   # min(d_cand, n - d_cand) — centered representation
    k1_avg_list = []
    k2_avg_list = []

    seed0_row = None  # store seed=0 spurious row for anatomy

    for seed in range(N_SEEDS):
        d_sec, sigs, k1v, k2v = gen_signatures(p, b_param, n, lam, G, m, K1, seed)
        if len(sigs) < m:
            continue

        M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
        A = IntegerMatrix.from_matrix(M_int)
        LLL.reduction(A)
        M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

        # Find shortest Kannan row
        best_row = None; best_norm = float('inf')
        for row in M_red:
            if abs(row[dim - 1]) != S_KANNAN: continue
            rn = math.sqrt(sum(x * x for x in row))
            if rn < best_norm:
                best_norm = rn
                best_row = list(row)

        if best_row is None:
            continue

        sign = 1 if best_row[dim - 1] > 0 else -1
        d_cand = (sign * best_row[m]) % n
        k1_decoded = [sign * best_row[i] / S_K1 for i in range(m)]
        k2_decoded = [sign * best_row[m + 1 + i] / S_K2 for i in range(m)]

        d_cand_list.append(d_cand)
        d_cand_small_list.append(min(d_cand, n - d_cand))
        k1_avg_list.append(sum(k1_decoded) / m)
        k2_avg_list.append(sum(k2_decoded) / m)

        if seed == 0:
            seed0_row = (best_row, sign, d_cand, k1_decoded, k2_decoded)

    print(f"  d_cand distribution across {len(d_cand_list)} seeds (where spurious row found):")
    print(f"  {'min(d_cand, n-d_cand)':>30}: min={min(d_cand_small_list)}, max={max(d_cand_small_list)}, mean={sum(d_cand_small_list)/len(d_cand_small_list):.1f}")
    d_set = set(d_cand_list)
    print(f"  Distinct d_cand values: {len(d_set)}")
    if len(d_set) <= 5:
        print(f"  Values: {sorted(d_set)}")
    else:
        sorted_d = sorted(d_cand_list)
        print(f"  Sample (first 5, last 5): {sorted_d[:5]} ... {sorted_d[-5:]}")

    print(f"\n  k1_i avg (across seeds): mean={sum(k1_avg_list)/len(k1_avg_list):.3f}, "
          f"std={math.sqrt(sum((x - sum(k1_avg_list)/len(k1_avg_list))**2 for x in k1_avg_list)/len(k1_avg_list)):.3f}")
    print(f"  k2_i avg (across seeds): mean={sum(k2_avg_list)/len(k2_avg_list):.3f}, "
          f"std={math.sqrt(sum((x - sum(k2_avg_list)/len(k2_avg_list))**2 for x in k2_avg_list)/len(k2_avg_list)):.3f}")

    # Seed-0 anatomy
    if seed0_row is not None:
        row, sign, d_cand_0, k1d_0, k2d_0 = seed0_row
        print(f"\n  Seed=0 spurious row anatomy:")
        print(f"  d_cand = {d_cand_0} (min-repr={min(d_cand_0, n-d_cand_0)})")
        print(f"  k1_i decoded (first 6): {[round(x,3) for x in k1d_0[:6]]}")
        print(f"  k2_i decoded (first 6): {[round(x,3) for x in k2d_0[:6]]}")
        # Are k1_i integers?
        k1_int = all(abs(x - round(x)) < 0.5 for x in k1d_0)
        k2_int = all(abs(x - round(x)) < 0.5 for x in k2d_0)
        print(f"  k1_i are near-integers: {k1_int}  (max frac: {max(abs(x-round(x)) for x in k1d_0):.4f})")
        print(f"  k2_i are near-integers: {k2_int}  (max frac: {max(abs(x-round(x)) for x in k2d_0):.4f})")
        # Check if d_cand satisfies a CM relation: d_cand * lam ≡ d_cand (mod n)?
        # This would mean d_cand is in the eigenspace
        lam_check = (d_cand_0 * lam) % n
        print(f"  CM check: d_cand * lam mod n = {lam_check} (vs d_cand={d_cand_0})")
        print(f"  d_cand mod K1 = {d_cand_0 % K1}")
        print(f"  d_cand mod K2 = {d_cand_0 % (math.isqrt(n)+1)}")

    # Consistency test: for seeds 1-5, print d_cand individually
    print(f"\n  Per-seed d_cand (first 10 seeds):")
    print(f"  {'seed':>5} {'d_cand':>10} {'min(d,n-d)':>12} {'k1_avg':>8} {'k2_avg':>8}")
    for i, (d, d_sm, k1a, k2a) in enumerate(zip(d_cand_list[:10], d_cand_small_list[:10],
                                                   k1_avg_list[:10], k2_avg_list[:10])):
        print(f"  {i:>5} {d:>10} {d_sm:>12} {k1a:>8.3f} {k2a:>8.3f}")

    print()

# ---------------------------------------------------------------------------
# Exp L: Cross-curve screening
# ---------------------------------------------------------------------------

def exp_l_cross_curve():
    print("=" * 70)
    print("Exp L: Cross-curve screening — 30 random j=0 CM curves (18-20 bit)")
    print("  Q: What algebraic invariant distinguishes C1-type (0%) from C2-type (>80%)?")
    print("=" * 70)

    K1 = 72
    m = 12
    N_SEEDS = 10
    TARGET_CURVES = 30

    # Generate primes p ≡ 1 (mod 3) in [2^18, 2^20]
    p_lo, p_hi = 1 << 18, 1 << 20

    curves = []
    rng_p = random.Random(20260626)
    tried = 0
    t0_total = time.time()

    while len(curves) < TARGET_CURVES and tried < 5000:
        p_cand = rng_p.randrange(p_lo | 1, p_hi, 2)  # odd candidates
        if not is_prime(p_cand):
            tried += 1
            continue
        if p_cand % 3 != 1:  # need p ≡ 1 (mod 3) for GLV endomorphism
            tried += 1
            continue

        b_cand, n_cand = find_prime_order_curve(p_cand, b_range=30, n_check_points=3)
        tried += 1
        if b_cand is None:
            continue

        lam_cand = glv_eigenvalue(n_cand)
        if lam_cand is None:
            continue

        G_cand = find_generator(p_cand, b_cand, n_cand)
        if G_cand is None:
            continue

        curves.append({
            'p': p_cand, 'n': n_cand, 'b': b_cand,
            'lam': lam_cand, 'G': G_cand,
            'lam_over_n': lam_cand / n_cand,
            'n_mod9': n_cand % 9,
            'n_mod6': n_cand % 6,
        })
        if len(curves) % 5 == 0:
            print(f"  Found {len(curves)}/{TARGET_CURVES} curves so far (tried {tried} primes)...")

    print(f"\n  Generated {len(curves)} curves from {tried} prime candidates in "
          f"{time.time()-t0_total:.1f}s\n")

    print(f"  Running LLL attack (K1={K1}, m={m}, {N_SEEDS} seeds per curve)...")
    print()
    print(f"  {'#':>3} {'p':>8} {'n':>8} {'b':>3} {'lam/n':>7} {'n%9':>4} {'recov':>7} {'type':>8}")
    print("  " + "-" * 60)

    c1_type = []  # recovery <= 20%
    c2_type = []  # recovery >= 80%
    mid_type = []  # 20% < recovery < 80%

    for idx, cv in enumerate(curves):
        p, n, b, lam, G = cv['p'], cv['n'], cv['b'], cv['lam'], cv['G']
        recover = 0
        for seed in range(N_SEEDS):
            ok, _, _, _, _ = lll_attack(p, b, n, lam, G, m, K1, seed)
            if ok:
                recover += 1
        cv['recovery'] = recover
        cv['recovery_rate'] = recover / N_SEEDS

        typ = "C2-type" if recover >= 8 else ("C1-type" if recover <= 2 else "MID")
        cv['type'] = typ

        print(f"  {idx+1:>3} {p:>8} {n:>8} {b:>3} {lam/n:>7.4f} {n%9:>4} {recover:>4}/{N_SEEDS}  {typ:>8}")

        if recover <= 2:
            c1_type.append(cv)
        elif recover >= 8:
            c2_type.append(cv)
        else:
            mid_type.append(cv)

    print(f"\n  Summary:")
    print(f"  C1-type (≤2/10): {len(c1_type)}")
    print(f"  C2-type (≥8/10): {len(c2_type)}")
    print(f"  MID-type (3-7/10): {len(mid_type)}")

    # Correlation analysis
    print(f"\n  --- Invariant analysis ---")
    if c1_type and c2_type:
        c1_lam = [cv['lam_over_n'] for cv in c1_type]
        c2_lam = [cv['lam_over_n'] for cv in c2_type]
        print(f"  lam/n: C1-type mean={sum(c1_lam)/len(c1_lam):.4f} "
              f"[{min(c1_lam):.4f},{max(c1_lam):.4f}]  "
              f"C2-type mean={sum(c2_lam)/len(c2_lam):.4f} "
              f"[{min(c2_lam):.4f},{max(c2_lam):.4f}]")

        c1_n9 = [cv['n_mod9'] for cv in c1_type]
        c2_n9 = [cv['n_mod9'] for cv in c2_type]
        from collections import Counter
        print(f"  n mod 9: C1-type counts={dict(Counter(c1_n9))}  "
              f"C2-type counts={dict(Counter(c2_n9))}")

        c1_n6 = [cv['n_mod6'] for cv in c1_type]
        c2_n6 = [cv['n_mod6'] for cv in c2_type]
        print(f"  n mod 6: C1-type counts={dict(Counter(c1_n6))}  "
              f"C2-type counts={dict(Counter(c2_n6))}")

        # lam ratio: does lam fall near 1/3 or 2/3 of n?
        c1_lam_frac = [min(cv['lam_over_n'], 1 - cv['lam_over_n']) for cv in c1_type]
        c2_lam_frac = [min(cv['lam_over_n'], 1 - cv['lam_over_n']) for cv in c2_type]
        print(f"  min(lam/n, 1-lam/n): C1-type mean={sum(c1_lam_frac)/len(c1_lam_frac):.4f} "
              f"C2-type mean={sum(c2_lam_frac)/len(c2_lam_frac):.4f}")

        # Check 3*lam mod n: is it small?
        c1_3lam = [(3 * cv['lam']) % cv['n'] for cv in c1_type]
        c2_3lam = [(3 * cv['lam']) % cv['n'] for cv in c2_type]
        c1_3lam_small = [min(x, cv['n'] - x) for x, cv in zip(c1_3lam, c1_type)]
        c2_3lam_small = [min(x, cv['n'] - x) for x, cv in zip(c2_3lam, c2_type)]
        print(f"  min(3*lam mod n, n - 3*lam mod n): "
              f"C1-type mean={sum(c1_3lam_small)/len(c1_3lam_small):.1f} "
              f"C2-type mean={sum(c2_3lam_small)/len(c2_3lam_small):.1f}")

        # Check lam mod K1
        K2 = math.isqrt(curves[0]['n']) + 1  # approximate
        c1_lam_K1 = [cv['lam'] % K1 for cv in c1_type]
        c2_lam_K1 = [cv['lam'] % K1 for cv in c2_type]
        print(f"  lam mod K1={K1}: C1-type mean={sum(c1_lam_K1)/len(c1_lam_K1):.1f} "
              f"C2-type mean={sum(c2_lam_K1)/len(c2_lam_K1):.1f}")

        # delta = lam - isqrt(lam) i.e. fractional part of lam in [0,n) : already lam/n
        # Check floor(lam/sqrt(n)) — the "effective GLV ratio"
        c1_ratio = [cv['lam'] // math.isqrt(cv['n']) for cv in c1_type]
        c2_ratio = [cv['lam'] // math.isqrt(cv['n']) for cv in c2_type]
        print(f"  floor(lam/sqrt(n)): C1-type mean={sum(c1_ratio)/len(c1_ratio):.1f} "
              f"[{min(c1_ratio)},{max(c1_ratio)}]  "
              f"C2-type mean={sum(c2_ratio)/len(c2_ratio):.1f} "
              f"[{min(c2_ratio)},{max(c2_ratio)}]")

        # lam mod (sqrt(n)): how lam decomposes in the K2-grid
        c1_lam_K2 = [cv['lam'] % (math.isqrt(cv['n']) + 1) for cv in c1_type]
        c2_lam_K2 = [cv['lam'] % (math.isqrt(cv['n']) + 1) for cv in c2_type]
        print(f"  lam mod (isqrt(n)+1)=K2: C1-type mean={sum(c1_lam_K2)/len(c1_lam_K2):.1f} "
              f"C2-type mean={sum(c2_lam_K2)/len(c2_lam_K2):.1f}")

    print()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("=== GLV-HNP Phase 2: Exp J/K/L — Structural source of C1 obstruction ===")
    print("Date: 2026-06-26\n")

    exp_j_k1_variation()
    exp_k_spurious_cm_check()
    exp_l_cross_curve()

    print("Done.")
