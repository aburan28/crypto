#!/usr/bin/env python3
"""
Thread 23 — GLV-HNP Phase 2 with the secret d ELIMINATED from the lattice.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5)
--------------------------------------------------------
The Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` (build_glv_lattice) carries the
secret d in its own column at scale S_D = 1.  Because d is a FULL-SIZE unknown
(d ~ U[0,n)) while its column scale is 1, the lattice contains the vector

    n * S_D * e_m  =  n * row_m  -  sum_i B_i * row_i          (norm exactly n)

whereas  ||v_planted|| ~ n * sqrt(2m/3 + 4/3) > n  for every m >= 1.  So the planted
vector is NEVER lambda_1, and the Kannan embedding is unsound as an SVP instance:
recovery is a BDD/coset condition, not an SVP condition.  T5 verified this empirically
(|sv[m]|/n = 1.0000 exactly on every curve tested) and noted that no choice of S_D
removes it, since both vectors scale linearly in S_D.

This script implements the fix proposed as Thread 23: eliminate d algebraically, so the
d-column -- and with it the trivial vector -- does not exist.

Construction
------------
Signatures i = 0..m-1 satisfy   A_i + B_i*d = k_i = k1_i + lam*k2_i  (mod n),
with 0 <= k1_i < K1 (small) and 0 <= k2_i < K2 ~ sqrt(n).

Pivot on i = 0:   d = B_0^{-1} (k1_0 + lam*k2_0 - A_0)  (mod n).
Substituting into equation i >= 1, with

    t_i = B_i * B_0^{-1} mod n,      u_i = A_i - t_i * A_0 mod n,

gives the d-free congruences

    k1_i  =  u_i + t_i*k1_0 + t_i*lam*k2_0 - lam*k2_i   (mod n),      i = 1..m-1.   (*)

Unknowns: k1_0 and k1_1..k1_{m-1} (each < K1), k2_0..k2_{m-1} (each < K2).

Lattice (dim = 2m+1), columns:
    0 .. m-2      congruence columns, one per i=1..m-1, scale S_K1;  hold k1_i * S_K1
    m-1           the k1_0 column,                      scale S_K1;  holds k1_0 * S_K1
    m .. 2m-1     the k2_j columns, j=0..m-1,           scale S_K2;  hold k2_j * S_K2
    2m            Kannan column,                        scale S_KAN = n

Rows (generators):
    0 .. m-2      n * S_K1 * e_col(i)                              (modular reduction)
    m-1           k1_0 generator:  t_i*S_K1 in col(i);  S_K1 in the k1_0 column
    m             k2_0 generator:  t_i*lam*S_K1 in col(i);  S_K2 in the k2_0 column
    m+1 .. 2m-1   k2_j generator (j>=1):  -lam*S_K1 in col(j-1);  S_K2 in the k2_j column
    2m            Kannan row:      u_i*S_K1 in col(i);  S_KAN in the last column

By (*), the combination  1*(Kannan) + k1_0*(k1_0 row) + k2_0*(k2_0 row)
+ sum_{j>=1} k2_j*(k2_j row)  -  (modular corrections)  equals

    v_planted = ( k1_1*S_K1, .., k1_{m-1}*S_K1 | k1_0*S_K1 | k2_0*S_K2, .., k2_{m-1}*S_K2 | n )

of expected norm  n * sqrt(2m/3 + 1)  -- the d-coordinate's n^2/3 is gone, and no vector
of norm ~n survives (the analogous "trivial" combinations now have norm n*S_K1 = n^2/K1
and n*S_K2, both >> n).

    det = n^m * S_K1^m * S_K2^m ~ n^{3m} / (K1*K2)^m   in dimension 2m+1,
so the Gaussian-heuristic viability condition is still K1*K2 < n, i.e. unchanged eff.

Recovery: any reduced row with |last| = S_KAN yields k1_0 = row[m-1]/S_K1 and
k2_0 = row[m]/S_K2, hence d = B_0^{-1}(k1_0 + lam*k2_0 - A_0) mod n.

Falsifier (stated in the 2026-07-29 log entry)
----------------------------------------------
  (i)  does sv/pv rise to 1 (planted vector becomes lambda_1)?
  (ii) does the K1 wall on the lam*=0.07 curve (12-bit/2677, currently K1 ~ 4-6)
       move OUTWARD?
If (ii) fails, the wall is information-theoretic and Phase 2 is at its ceiling.

Arithmetic helpers are reused verbatim from glv_hnp_phase2_lambda_threshold.py so the
comparison against 2026-07-26 / 2026-07-29 is exact.
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Arithmetic helpers (verbatim from glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = 3 * x1 * x1 % p * modinv(2 * y1 % p, p) % p
    else:
        lam = (y2 - y1) % p * modinv((x2 - x1) % p, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    R = None
    Q = P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(a, p):
    a %= p
    if a == 0: return 0
    if p == 2: return a
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m_, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m_ - i - 1), p)
        m_, c = i, b * b % p
        t = t * c % p
        r = r * b % p
    return r

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(2000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Signatures (verbatim) and the ORIGINAL Phase-2 lattice (verbatim, for A/B)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
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
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    """ORIGINAL (d-in-lattice) construction, verbatim from glv_hnp_phase2_20bit.py:262."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
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
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector_orig(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_orig(M_reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# THREAD 23 — the d-ELIMINATED lattice
# ---------------------------------------------------------------------------

def eliminate_d(sigs, n):
    """Pivot on signature 0.  Returns (t, u) with t[i], u[i] defined for i=1..m-1
    (index 0 unused), such that  k1_i = u_i + t_i*k1_0 + t_i*lam*k2_0 - lam*k2_i (mod n)."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'] % n, n)
    A0 = sigs[0]['A'] % n
    t = [0] * m
    u = [0] * m
    for i in range(1, m):
        t[i] = sigs[i]['B'] % n * B0inv % n
        u[i] = (sigs[i]['A'] - t[i] * A0) % n
    return t, u

def build_elim_lattice(sigs, n, lam, k1_bound, k2_bound):
    """d-eliminated lattice, dim = 2m+1.  Column layout:
         0..m-2 : congruence cols (i=1..m-1), scale S_K1
         m-1    : k1_0 col,                   scale S_K1
         m..2m-1: k2_j cols (j=0..m-1),       scale S_K2
         2m     : Kannan col,                 scale S_KAN
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    t, u = eliminate_d(sigs, n)
    lm = lam % n

    M = [[0] * dim for _ in range(dim)]

    # modular reduction rows for the m-1 congruence columns
    for j in range(m - 1):
        M[j][j] = n * S_K1

    # k1_0 generator
    for j in range(m - 1):
        M[m - 1][j] = t[j + 1] * S_K1
    M[m - 1][m - 1] = S_K1

    # k2_0 generator
    for j in range(m - 1):
        M[m][j] = t[j + 1] * lm % n * S_K1
    M[m][m] = S_K2

    # k2_j generators, j = 1..m-1
    for j in range(1, m):
        M[m + j][j - 1] = -lm * S_K1
        M[m + j][m + j] = S_K2

    # Kannan row
    for j in range(m - 1):
        M[2 * m][j] = u[j + 1] * S_K1
    M[2 * m][2 * m] = S_KAN

    return M

def planted_vector_elim(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for j in range(m - 1):
        v[j] = sigs[j + 1]['k1'] * S_K1
    v[m - 1] = sigs[0]['k1'] * S_K1
    for j in range(m):
        v[m + j] = sigs[j]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def verify_membership(sigs, n, lam, k1_bound, k2_bound):
    """Prove v_planted is in the eliminated lattice by exhibiting the integer
    combination of rows and checking it equals planted_vector_elim()."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M = build_elim_lattice(sigs, n, lam, k1_bound, k2_bound)

    # coefficients: 1*Kannan + k1_0*(k1_0 row) + sum_j k2_j*(k2_j row)
    coef = [0] * dim
    coef[2 * m] = 1
    coef[m - 1] = sigs[0]['k1']
    for j in range(m):
        coef[m + j] = sigs[j]['k2']

    w = [0] * dim
    for r in range(dim):
        if coef[r] == 0: continue
        for c in range(dim):
            if M[r][c]:
                w[c] += coef[r] * M[r][c]

    # subtract the modular corrections on the congruence columns
    mod = n * S_K1
    for j in range(m - 1):
        q = w[j] // mod
        w[j] -= q * mod
        # (w[j] now in [0, n*S_K1); k1_{j+1}*S_K1 lies in that range)

    return w == planted_vector_elim(sigs, n, k1_bound, k2_bound), w

def recover_d_elim(M_reduced, m, n, sigs, lam, k1_bound, k2_bound, d_secret):
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    A0 = sigs[0]['A'] % n
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sgn = 1 if last > 0 else -1
        r = [sgn * x for x in row]
        if r[m - 1] % S_K1 or r[m] % S_K2:
            continue
        k1_0 = r[m - 1] // S_K1
        k2_0 = r[m] // S_K2
        d_cand = B0inv * (k1_0 + lam * k2_0 - A0) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def reduce_matrix(M, use_bkz=False, beta=20):
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def run_orig(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    red = reduce_matrix(M, use_bkz, beta)
    S_KAN = scales(n, k1_bound, k2_bound)[3]
    ok = recover_d_orig(red, m, n, S_KAN, d_secret) is not None
    pv = norm(planted_vector_orig(sigs, d_secret, n, k1_bound, k2_bound))
    sv = min(norm(r) for r in red if any(r))
    return {'ok': ok, 'pv': pv, 'sv': sv, 'ratio': sv / pv}

def run_elim(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    M = build_elim_lattice(sigs, n, lam, k1_bound, k2_bound)
    red = reduce_matrix(M, use_bkz, beta)
    ok = recover_d_elim(red, m, n, sigs, lam, k1_bound, k2_bound, d_secret) is not None
    pv = norm(planted_vector_elim(sigs, n, k1_bound, k2_bound))
    sv = min(norm(r) for r in red if any(r))
    return {'ok': ok, 'pv': pv, 'sv': sv, 'ratio': sv / pv}

def rate(fn, curve, m, k1_bound, seeds, use_bkz=False, beta=20):
    n = curve[2]
    good = 0
    ratios = []
    for s in seeds:
        d = random.Random(s * 7919 + 13).randint(1, n - 1)
        r = fn(curve, m, d, k1_bound, seed=s, use_bkz=use_bkz, beta=beta)
        if r is None: continue
        good += 1 if r['ok'] else 0
        ratios.append(r['ratio'])
    return good, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================
def main():
    print("=" * 78)
    print("Thread 23 — Phase-2 GLV-HNP lattice with d ELIMINATED")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,             p,    b, n,    lam,  K1, m
        ("8-bit/199",        211,  2, 199,  106,  2,  6),
        ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
    ]

    curves = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        curves.append((label, (p, b, n, lam, G), k1, m))

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E1: correctness — planted vector IS in the eliminated lattice")
    print("-" * 78)
    print(f"{'curve':<20} {'m':>3} {'K1':>4} {'dim_orig':>9} {'dim_elim':>9} {'member':>8}")
    for label, curve, k1, m in curves:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d = random.Random(42 * 7919 + 13).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, 42)
        ok, w = verify_membership(sigs, n, lam, k1, k2b)
        print(f"{label:<20} {m:>3} {k1:>4} {2*m+2:>9} {2*m+1:>9} {str(ok):>8}")
        assert ok, f"membership FAILED for {label}"
    print("All planted vectors verified as lattice members (exact integer check).")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E2: is the planted vector lambda_1 now?  (sv/pv, mean over 5 seeds)")
    print("-" * 78)
    print("sv = shortest reduced row, pv = ||v_planted||.  sv/pv = 1 means the")
    print("planted vector IS the shortest vector found.  T5 measured 0.34-0.60 for orig.")
    print()
    print(f"{'curve':<20} {'K1':>4} {'m':>3} {'sv/pv orig':>11} {'sv/pv elim':>11} "
          f"{'ok orig':>9} {'ok elim':>9}")
    for label, curve, k1, m in curves:
        go, tot, ro = rate(run_orig, curve, m, k1, SEEDS)
        ge, _,   re = rate(run_elim, curve, m, k1, SEEDS)
        print(f"{label:<20} {k1:>4} {m:>3} {ro:>11.3f} {re:>11.3f} "
              f"{str(go)+'/'+str(tot):>9} {str(ge)+'/'+str(tot):>9}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E3: the K1 wall — does elimination move it OUTWARD?  (the falsifier)")
    print("-" * 78)
    print("Grid replicates 2026-07-29 exp T4 exactly (same curves, K2=52, 5 seeds).")
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
    for label, curve, _k1, m in curves[1:]:
        p, b, n, lam, G = curve
        print(f"\n{label}   n={n}  lam*={lam_star(lam,n):.3f}  m={m}  K2={math.isqrt(n)+1}")
        hdr = "  ".join(f"{k:>5}" for k in K1_GRID)
        print(f"{'K1':<14} {hdr}")
        for name, fn in (("orig (d in L)", run_orig), ("elim (Thr23)", run_elim)):
            cells = []
            for k1 in K1_GRID:
                g, tot, _ = rate(fn, curve, m, k1, SEEDS)
                cells.append(f"{g}/{tot}")
            print(f"{name:<14} " + "  ".join(f"{c:>5}" for c in cells))

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E4: BKZ-20 at the wall — reduction strength or information-theoretic?")
    print("-" * 78)
    label, curve, _k1, m = curves[2]
    p, b, n, lam, G = curve
    print(f"{label}  n={n}  lam*={lam_star(lam,n):.3f}  m={m}")
    print(f"{'K1':>4} {'orig LLL':>9} {'orig BKZ20':>11} {'elim LLL':>9} {'elim BKZ20':>11}")
    for k1 in [6, 8, 12, 16]:
        a, tot, _ = rate(run_orig, curve, m, k1, SEEDS)
        bb, _, _  = rate(run_orig, curve, m, k1, SEEDS, use_bkz=True, beta=20)
        c, _, _   = rate(run_elim, curve, m, k1, SEEDS)
        dd, _, _  = rate(run_elim, curve, m, k1, SEEDS, use_bkz=True, beta=20)
        print(f"{k1:>4} {str(a)+'/'+str(tot):>9} {str(bb)+'/'+str(tot):>11} "
              f"{str(c)+'/'+str(tot):>9} {str(dd)+'/'+str(tot):>11}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E5: m-scaling at the wall (K1=8, lam*=0.07 curve)")
    print("-" * 78)
    print("T4b showed orig does not recover at K1=8 for any m in {8,12,16,24,32}.")
    print(f"{'m':>4} {'dim_orig':>9} {'orig':>7} {'dim_elim':>9} {'elim':>7}")
    for mm in [8, 12, 16, 24, 32]:
        a, tot, _ = rate(run_orig, curve, mm, 8, SEEDS)
        c, _, _   = rate(run_elim, curve, mm, 8, SEEDS)
        print(f"{mm:>4} {2*mm+2:>9} {str(a)+'/'+str(tot):>7} {2*mm+1:>9} "
              f"{str(c)+'/'+str(tot):>7}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
