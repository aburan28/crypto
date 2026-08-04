"""
GLV-HNP Phase 2: Target-vector length + spurious-vector analysis for PAIR 1.

From 2026-06-23: κ(M) hypothesis falsified. The discriminating property between
discordant PAIR 1 curves (same λ/n≈0.20, opposite LLL outcomes) is UNKNOWN.

Today's two sub-experiments for PAIR 1:
  FAIL  curve: p=524743, n=523597 (lam/n=0.2114, δ/n=0.366)
  SUCCEED curve: p=525043, n=524269 (lam/n=0.2122, δ/n=0.364)

Exp A — Target vector norm:
  Compute ‖v_target‖ for K1=72, m=12, seeds=[0xDEAD, 0xBEEF, 0xCAFE].
  If ‖v_fail‖ >> ‖v_succeed‖, obstruction is a size mismatch.

Exp B — Spurious-vector check:
  Run LLL at m=12, inspect ALL rows of reduced basis.
  Count rows with |last_col| == S_KANNAN (Kannan-embedded rows).
  If failing curve has 0 such rows, the solution was never found.
  If it has rows but none match d, there are spurious Kannan vectors.
"""

import math
import random
import sympy
import numpy as np
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# EC arithmetic (from glv_hnp_conditioning.py)
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
# Experiment A: Target vector norm
# ---------------------------------------------------------------------------

def target_vector_norm(k1_vals, k2_vals, d, S_K1, S_K2, S_KANNAN):
    norm_sq = sum(k1 * k1 * S_K1 * S_K1 for k1 in k1_vals)
    norm_sq += d * d
    norm_sq += sum(k2 * k2 * S_K2 * S_K2 for k2 in k2_vals)
    norm_sq += S_KANNAN * S_KANNAN
    return math.sqrt(norm_sq)

# ---------------------------------------------------------------------------
# Experiment B: Spurious-vector check
# ---------------------------------------------------------------------------

def spurious_vector_check(M_red, m, n, d_secret, S_KANNAN):
    dim = 2 * m + 2
    kannan_rows = 0
    correct_rows = 0
    row_details = []
    for row_idx, row in enumerate(M_red):
        last = row[dim - 1]
        if abs(last) == S_KANNAN:
            kannan_rows += 1
            sign = 1 if last > 0 else -1
            d_cand = (sign * row[m]) % n
            is_correct = (d_cand == d_secret)
            if is_correct:
                correct_rows += 1
            row_details.append((row_idx, last, d_cand, is_correct))
    return kannan_rows, correct_rows, kannan_rows - correct_rows, row_details

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_pair1_curve(label, p, n, K1, m, seeds):
    print(f"\n{'='*60}")
    print(f"  Curve: {label}  p={p}  n={n}")

    lam = glv_eigenvalue(n)
    if lam is None:
        print("  ERROR: no GLV eigenvalue found")
        return

    b_param = find_b_for_n(p, n)
    if b_param is None:
        print("  ERROR: no b_param found")
        return

    lam_n = lam / n
    K2 = math.isqrt(n) + 1
    S_K1 = n // K1
    S_K2 = max(1, n // K2)
    S_KANNAN = n
    print(f"  lam/n={lam_n:.4f}  K1={K1}  K2={K2}  S_K1={S_K1}  S_K2={S_K2}")

    G = find_generator(p, b_param, n)
    if G is None:
        print("  ERROR: no generator found")
        return

    print(f"\n  {'Seed':<12} {'‖v_target‖':>14} {'GH_estimate':>13} {'v/GH':>8} {'Recovered':>10} {'Kannan#':>8} {'Correct':>8} {'Spurious':>9}")
    print("  " + "-" * 90)

    dim = 2 * m + 2
    log2_vol = m * math.log2(n * S_K1) + m * math.log2(S_K2) + math.log2(S_KANNAN)
    gh_est = math.sqrt(dim / (2 * math.pi * math.e)) * 2 ** (log2_vol / dim)

    for seed in seeds:
        d_secret, sigs, k1_vals, k2_vals = gen_signatures(
            p, b_param, n, lam, G, m, K1, seed)

        if len(sigs) < m:
            print(f"  seed={hex(seed):<10}  FAIL: only {len(sigs)} sigs generated")
            continue

        v_norm = target_vector_norm(k1_vals, k2_vals, d_secret, S_K1, S_K2, S_KANNAN)

        M_int, _, _, _ = build_lattice(sigs, n, lam, K1)
        A = IntegerMatrix.from_matrix(M_int)
        LLL.reduction(A)
        M_red = [[A[i][j] for j in range(dim)] for i in range(dim)]

        recovered = False
        for row in M_red:
            last = row[dim - 1]
            if abs(last) == S_KANNAN:
                sign = 1 if last > 0 else -1
                d_cand = (sign * row[m]) % n
                if d_cand == d_secret:
                    recovered = True
                    break

        k_rows, c_rows, s_rows, details = spurious_vector_check(
            M_red, m, n, d_secret, S_KANNAN)

        print(f"  seed={hex(seed):<10} "
              f"{v_norm:14.0f} {gh_est:13.0f} {v_norm/gh_est:8.3f} "
              f"{'YES':>10}" if recovered else
              f"  seed={hex(seed):<10} "
              f"{v_norm:14.0f} {gh_est:13.0f} {v_norm/gh_est:8.3f} "
              f"{'NO':>10}", end="")
        print(f" {k_rows:8d} {c_rows:8d} {s_rows:9d}")

        if details:
            for (ridx, last_col, d_cand, ok) in details[:5]:
                status = "CORRECT" if ok else "spurious"
                print(f"    row[{ridx}]: last_col={last_col}  d_cand={d_cand}  [{status}]")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    K1 = 72
    m = 12
    SEEDS = [0xDEAD, 0xBEEF, 0xCAFE]

    print("=== GLV-HNP Target-Vector Length + Spurious-Vector Analysis ===")
    print("PAIR 1: λ/n≈0.20 discordant pair")
    print(f"K1={K1}, m={m}, seeds={[hex(s) for s in SEEDS]}\n")

    print("Searching for PAIR 1 curves (scanning 20-bit CM primes)...")

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

    TARGET = 0.20
    TOL = 0.015
    found = []
    p = sympy.nextprime(2**19 - 1)
    while p < 2**20 and len(found) < 2:
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
                    lam_n_c = lam_cand / n_cand
                    if abs(lam_n_c - TARGET) <= TOL:
                        b_p = find_b_for_n(p, n_cand)
                        if b_p is None: continue
                        dr = delta_ratio(lam_cand, n_cand)
                        found.append((p, n_cand, lam_cand, lam_n_c, dr))
                        print(f"  Found: p={p}, n={n_cand}, lam/n={lam_n_c:.4f}, δ/n={dr:.4f}")
                        break
        p = sympy.nextprime(p)

    if len(found) < 2:
        print(f"Only found {len(found)} curve(s) near lam/n=0.20. Need 2 for pair analysis.")

    for i, (p_c, n_c, lam_c, lam_n_c, dr_c) in enumerate(found):
        label = f"PAIR1_C{i+1} (lam/n={lam_n_c:.4f}, δ/n={dr_c:.4f})"
        analyze_pair1_curve(label, p_c, n_c, K1, m, SEEDS)

    print("\n\n=== Extended LLL sweep m=10..22 for PAIR 1 ===")
    for i, (p_c, n_c, lam_c, lam_n_c, dr_c) in enumerate(found):
        b_param = find_b_for_n(p_c, n_c)
        if b_param is None:
            print(f"C{i+1}: no b_param")
            continue
        G = find_generator(p_c, b_param, n_c)
        if G is None:
            print(f"C{i+1}: no generator")
            continue
        K2 = math.isqrt(n_c) + 1
        S_K1 = n_c // K1
        S_K2 = max(1, n_c // K2)
        S_KANNAN = n_c
        lam_c_val = glv_eigenvalue(n_c)
        print(f"\n  C{i+1}: p={p_c}, n={n_c}, lam/n={lam_n_c:.4f}", end="")
        sweep_str = ""
        first_3of3 = None
        for m_sweep in range(10, 23):
            wins = 0
            for seed in SEEDS:
                d_sec, sigs, k1v, k2v = gen_signatures(
                    p_c, b_param, n_c, lam_c_val, G, m_sweep, K1, seed)
                if len(sigs) < m_sweep:
                    continue
                M_int, _, _, _ = build_lattice(sigs, n_c, lam_c_val, K1)
                dim = 2 * m_sweep + 2
                A = IntegerMatrix.from_matrix(M_int)
                LLL.reduction(A)
                M_red = [[A[j2][j3] for j3 in range(dim)] for j2 in range(dim)]
                for row in M_red:
                    last = row[dim - 1]
                    if abs(last) == S_KANNAN:
                        sign = 1 if last > 0 else -1
                        d_cand = (sign * row[m_sweep]) % n_c
                        if d_cand == d_sec:
                            wins += 1
                            break
            sweep_str += f" m{m_sweep}:{wins}"
            if wins == 3 and first_3of3 is None:
                first_3of3 = m_sweep
        print(f"\n  {sweep_str}")
        f33 = f"m={first_3of3}" if first_3of3 else "NEVER"
        print(f"  First 3/3: {f33}")

    print("\nDone.")
