"""
GLV-HNP Phase 2, Thread 23: remove the trivial short vector from the lattice.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  In the Phase-2 lattice of glv_hnp_phase2_20bit.py:263 (dim 2m+2, columns
  [k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan]) the shortest vector found by
  LLL is ALWAYS the trivial vector  n*S_D*e_m  (100% of its energy in the
  d-column, |sv[m]|/n = 1.0000 exactly).  It carries no information: d is only
  defined mod n.  Since ||v_planted||^2 ~ n^2*(2m/3 + 4/3) and ||n*S_D*e_m|| =
  n*S_D, the planted vector is never lambda_1 for any m >= 1.

  That entry claimed "no choice of S_D removes it - both vectors scale
  linearly in S_D", and proposed Thread 23: reformulate so the target IS
  lambda_1, with the falsifier "if sv/pv rises above 1 AND the K1 wall on the
  lam*=0.07 curve (currently K1 ~ 4-6) moves outward, the reformulation is a
  real improvement; if the wall stays, the wall is information-theoretic".

This script does two things.

  (A) CORRECTS the S_D claim.  ||v_planted||^2 = R^2 + d^2*S_D^2 where
      R^2 = m*(K1*S_K1)^2/3 + m*(K2*S_K2)^2/3 + S_KANNAN^2 does NOT depend on
      S_D, while ||n*S_D*e_m|| = n*S_D does.  So the trivial vector overtakes
      the planted one as soon as
            S_D > R / sqrt(n^2 - d^2),
      and the achievable norm ratio saturates at n/d (<= 2 for typical d).
      EXP S sweeps S_D and measures both the crossover and the recovery rate.

  (B) Builds the PIVOT lattice, which removes the trivial vector structurally
      by eliminating d instead of scaling it.

      Eliminate d using signature 0 as pivot.  With k_i = A_i + B_i*d (mod n)
      and k_i = k1_i + lam*k2_i, put gam_i = B_i/B_0 mod n and
      alpha_i = A_i - gam_i*A_0 mod n.  Then for i = 1..m-1

          k1_i  ==  alpha_i + gam_i*k1_0 + gam_i*lam*k2_0 - lam*k2_i   (mod n)

      Unknowns: k1_0, k2_0, and (k1_i, k2_i) for i >= 1.  There is no d
      coordinate, hence no free column, hence no n*e_d vector.  Dimension is
      2m+1 (one less than the original) and the determinant is det_orig / n.

      Columns: [k1_1..k1_{m-1} | k1_0 | k2_0 | k2_1..k2_{m-1} | kannan]
      Planted: (k1_i*S_K1 ..., k1_0*S_K1, k2_0*S_K2, k2_i*S_K2 ..., S_KANNAN)
      so ||v_pivot||^2 = ||v_orig||^2 - (d*S_D)^2 exactly.

      d is read back from the recovered k1_0, k2_0 as d = (k1_0 + lam*k2_0
      - A_0)/B_0 mod n.

Experiments
  E0  membership + norm identity check for the pivot lattice
  E1  sv/pv (shortest-reduced / planted) in both lattices, 3 historical curves
  E2  K1-wall grid, orig vs pivot, LLL and BKZ(20), 5 seeds  <- the falsifier
  E3  m-sweep at the wall
  S   S_D sweep in the ORIGINAL lattice (claim (A) above)

Helper functions (modinv .. lam_star, scales, gen_signatures,
build_glv_lattice, planted_vector, recover_d) are copied verbatim from
glv_hnp_phase2_lambda_threshold.py so the comparison to the 2026-07-29 numbers
is exact.

Run: python3 glv_hnp_phase2_pivot_lattice.py
Deps: fpylll (LLL/BKZ).  sympy is NOT needed (no curve search here).
"""

import math
import random
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Helpers -- verbatim from glv_hnp_phase2_lambda_threshold.py
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
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def norm(v):
    return math.sqrt(sum(x * x for x in v))

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

# ---------------------------------------------------------------------------
# ORIGINAL Phase-2 lattice (dim 2m+2), S_D made an explicit parameter
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, S_D=None):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D_def, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    if S_D is None:
        S_D = S_D_def
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

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound, S_D=None):
    m = len(sigs)
    S_K1, S_D_def, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if S_D is None:
        S_D = S_D_def
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d(M_reduced, m, n, S_KANNAN, d_secret, S_D=1):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % S_D != 0: continue
        d_cand = (val // S_D) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# PIVOT lattice (dim 2m+1) -- d eliminated, no trivial n*e_d vector
# ---------------------------------------------------------------------------

def pivot_params(sigs, n):
    """gam_i = B_i/B_0 mod n, alpha_i = A_i - gam_i*A_0 mod n, for i=1..m-1."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'], n)
    gam = [0] * m
    alp = [0] * m
    for i in range(1, m):
        gam[i] = sigs[i]['B'] * B0inv % n
        alp[i] = (sigs[i]['A'] - gam[i] * sigs[0]['A']) % n
    return gam, alp

def build_pivot_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Columns: [res_1..res_{m-1} | k1_0 | k2_0 | k2_1..k2_{m-1} | kannan]
       res_j (col j-1) carries k1_j*S_K1 for j = 1..m-1."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    gam, alp = pivot_params(sigs, n)
    L = lam % n

    M = [[0] * dim for _ in range(dim)]
    # rows 0..m-2 : modular reduction n in each residue column
    for j in range(1, m):
        M[j - 1][j - 1] = n * S_K1
    # row m-1 : the k1_0 generator
    for j in range(1, m):
        M[m - 1][j - 1] = gam[j] % n * S_K1
    M[m - 1][m - 1] = S_K1
    # row m : the k2_0 generator
    for j in range(1, m):
        M[m][j - 1] = (gam[j] * L) % n * S_K1
    M[m][m] = S_K2
    # rows m+1..2m-1 : the k2_j generators, j = 1..m-1
    for j in range(1, m):
        M[m + j][j - 1] = (-L) % n * S_K1
        M[m + j][m + j] = S_K2
    # row 2m : Kannan embedding
    for j in range(1, m):
        M[2 * m][j - 1] = alp[j] * S_K1
    M[2 * m][dim - 1] = S_KAN
    return M

def pivot_planted_vector(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * dim
    for j in range(1, m):
        v[j - 1] = sigs[j]['k1'] * S_K1
        v[m + j] = sigs[j]['k2'] * S_K2
    v[m - 1] = sigs[0]['k1'] * S_K1
    v[m] = sigs[0]['k2'] * S_K2
    v[dim - 1] = S_KAN
    return v

def pivot_recover_d(M_reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read k1_0,k2_0 off any row with |last| = S_KAN, then d = (k_0-A_0)/B_0."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        a, b = sign * row[m - 1], sign * row[m]
        if a % S_K1 != 0 or b % S_K2 != 0: continue
        k1_0, k2_0 = a // S_K1, b // S_K2
        k0 = (k1_0 + lam * k2_0) % n
        if k0 == 0: continue
        d_cand = (k0 - sigs[0]['A']) * B0inv % n
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def reduce_matrix(M, use_bkz=False, beta=20):
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def shortest_norm(rows):
    best = None
    for r in rows:
        nr = norm(r)
        if nr > 0 and (best is None or nr < best):
            best = nr
    return best

def run_orig(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20,
             S_D=None):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D_def, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if S_D is None:
        S_D = S_D_def
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, S_D=S_D)
    pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound, S_D=S_D)
    red = reduce_matrix(M, use_bkz, beta)
    ok = recover_d(red, m, n, S_KAN, d_secret, S_D=S_D) is not None
    sv = shortest_norm(red)
    triv = float(n) * S_D
    return {'ok': ok, 'pv': norm(pv), 'sv': sv, 'triv': triv,
            'sv_over_pv': sv / norm(pv), 'sigs': sigs}

def run_pivot(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_pivot_lattice(sigs, n, lam, k1_bound, k2_bound)
    pv = pivot_planted_vector(sigs, n, k1_bound, k2_bound)
    red = reduce_matrix(M, use_bkz, beta)
    ok = pivot_recover_d(red, sigs, n, lam, k1_bound, k2_bound, d_secret) is not None
    sv = shortest_norm(red)
    return {'ok': ok, 'pv': norm(pv), 'sv': sv,
            'sv_over_pv': sv / norm(pv), 'sigs': sigs}

def rate(fn, seeds):
    hits = 0
    ratios = []
    for s in seeds:
        r = fn(s)
        if r is None:
            continue
        hits += 1 if r['ok'] else 0
        ratios.append(r['sv_over_pv'])
    return hits, (sum(ratios) / len(ratios) if ratios else float('nan'))


# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]
D_SECRET_FRAC = 0.4213   # d = round(frac*n); fixed so runs are reproducible

HIST = [
    # label,          p,    b, n,    lam,  K1, m
    ("8-bit/199",     211,  2, 199,  106,  2,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",   2677, 2, 2647, 185,  8,  10),
]


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 23 - remove the trivial short vector (GLV-HNP Phase 2)")
    print("=" * 78)

    curves = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        d_sec = max(1, int(round(D_SECRET_FRAC * n)))
        curves.append((label, (p, b, n, lam, G), k1, m, d_sec))

    print(f"\nseeds = {SEEDS}   d_secret = round({D_SECRET_FRAC}*n)")
    for label, cv, k1, m, d in curves:
        p, b, n, lam, G = cv
        K2 = math.isqrt(n) + 1
        print(f"  {label:14s} p={p:6d} n={n:6d} lam={lam:6d} lam*={lam_star(lam,n):.4f} "
              f"K1={k1} K2={K2} m={m} eff={k1*K2/n:.4f} d={d}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E0: pivot-lattice sanity -- membership and the norm identity")
    print("-" * 78)
    print("Checks, for each curve: (i) the congruence k1_j == alpha_j + gam_j*k1_0")
    print("+ gam_j*lam*k2_0 - lam*k2_j (mod n) holds for the true nonces; (ii) the")
    print("planted vector is an integer combination of the pivot basis rows;")
    print("(iii) ||v_pivot||^2 == ||v_orig||^2 - (d*S_D)^2 exactly.")
    print()
    print(f"{'curve':<14} {'congr':>6} {'in-lat':>7} {'norm-id':>8} {'||v_orig||/n':>13} {'||v_piv||/n':>12}")
    e0_all_ok = True
    for label, cv, k1, m, d_sec in curves:
        p, b, n, lam, G = cv
        K2 = math.isqrt(n) + 1
        sigs = gen_signatures(G, d_sec, m, n, lam, p, k1, K2, seed=42)
        gam, alp = pivot_params(sigs, n)
        congr = all(
            (alp[j] + gam[j] * sigs[0]['k1'] + gam[j] * lam * sigs[0]['k2']
             - lam * sigs[j]['k2'] - sigs[j]['k1']) % n == 0
            for j in range(1, m))
        # (ii) rebuild the planted vector as an explicit integer combination
        M = build_pivot_lattice(sigs, n, lam, k1, K2)
        dim = 2 * m + 1
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, K2)
        coef = [0] * dim
        coef[2 * m] = 1
        coef[m - 1] = sigs[0]['k1']
        coef[m] = sigs[0]['k2']
        for j in range(1, m):
            coef[m + j] = sigs[j]['k2']
        comb = [sum(coef[r] * M[r][c] for r in range(dim)) for c in range(dim)]
        # subtract the modular rows to land in [0,n)
        for j in range(1, m):
            q = comb[j - 1] // (n * S_K1)
            coef[j - 1] -= q
            comb[j - 1] -= q * n * S_K1
        pv = pivot_planted_vector(sigs, n, k1, K2)
        in_lat = (comb == pv)
        ov = planted_vector(sigs, d_sec, n, k1, K2)
        norm_id = abs(norm(pv) ** 2 - (norm(ov) ** 2 - (d_sec * S_D) ** 2)) < 1e-6
        e0_all_ok &= congr and in_lat and norm_id
        print(f"{label:<14} {str(congr):>6} {str(in_lat):>7} {str(norm_id):>8} "
              f"{norm(ov)/n:>13.3f} {norm(pv)/n:>12.3f}")
    print(f"\nE0 verdict: {'ALL PASS' if e0_all_ok else 'FAIL'}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E1: is the planted vector lambda_1 now?  sv/pv in both lattices")
    print("-" * 78)
    print("sv = shortest row after LLL, pv = ||v_planted||.  sv/pv >= 1 means the")
    print("planted vector is (at least as short as) the shortest reduced vector.")
    print("'triv/pv' is the trivial vector n*S_D*e_d, absent from the pivot lattice.")
    print()
    print(f"{'curve':<14} {'K1':>3} {'m':>3} | {'orig sv/pv':>10} {'triv/pv':>8} {'orig ok':>8} | "
          f"{'piv sv/pv':>10} {'piv ok':>7}")
    for label, cv, k1, m, d_sec in curves:
        o = run_orig(cv, m, d_sec, k1, seed=42)
        pvt = run_pivot(cv, m, d_sec, k1, seed=42)
        print(f"{label:<14} {k1:>3} {m:>3} | {o['sv_over_pv']:>10.4f} "
              f"{o['triv']/o['pv']:>8.4f} {str(o['ok']):>8} | "
              f"{pvt['sv_over_pv']:>10.4f} {str(pvt['ok']):>7}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E2: the K1 wall -- orig vs pivot  (THE FALSIFIER)")
    print("-" * 78)
    print("2026-07-29 T4 measured, with LLL on the original lattice:")
    print("  12-bit/2557 (lam*=0.340): 5/5 up to K1=8, 4/5 at 12, 1/5 at 16, 0/5 at 24")
    print("  12-bit/2677 (lam*=0.070): 5/5 up to K1=4, 2/5 at 6,  0/5 from 8 on")
    print("Falsifier: if the pivot lattice moves the 2677 wall outward, the")
    print("reformulation is a real improvement; if not, the wall is not about the")
    print("trivial vector at all.")
    print()
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    for label, cv, k1_hist, m, d_sec in curves:
        if label == "8-bit/199":
            continue
        print(f"\n  {label}  (m={m}, 5 seeds; entries = #recovered/5)")
        hdr = "    " + f"{'lattice':<16}" + "".join(f"{k:>6}" for k in K1_GRID)
        print(hdr)
        for name, fn in (
            ("orig LLL",    lambda c, s, K: run_orig(cv, m, d_sec, K, seed=s)),
            ("pivot LLL",   lambda c, s, K: run_pivot(cv, m, d_sec, K, seed=s)),
            ("pivot BKZ20", lambda c, s, K: run_pivot(cv, m, d_sec, K, seed=s,
                                                      use_bkz=True, beta=20)),
        ):
            cells = []
            for K in K1_GRID:
                h, _ = rate(lambda s, K=K: fn(cv, s, K), SEEDS)
                cells.append(f"{h}/5")
            print("    " + f"{name:<16}" + "".join(f"{c:>6}" for c in cells))

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E3: does more data move the pivot wall?  (2677, K1=8)")
    print("-" * 78)
    print("2026-07-29 T4b: on the original lattice, m = 8/12/16/24/32 gave")
    print("0,0,1,0,1 out of 5 -- more signatures did not rescue K1=8.")
    print()
    print(f"{'m':>4} {'orig LLL':>10} {'pivot LLL':>10} {'pivot BKZ20':>12} {'piv sv/pv':>10}")
    cv2677 = curves[2][1]
    d2677 = curves[2][4]
    for m_try in (8, 12, 16, 24, 32):
        ho, _ = rate(lambda s: run_orig(cv2677, m_try, d2677, 8, seed=s), SEEDS)
        hp, rp = rate(lambda s: run_pivot(cv2677, m_try, d2677, 8, seed=s), SEEDS)
        hb, _ = rate(lambda s: run_pivot(cv2677, m_try, d2677, 8, seed=s,
                                         use_bkz=True, beta=20), SEEDS)
        print(f"{m_try:>4} {ho:>8}/5 {hp:>8}/5 {hb:>10}/5 {rp:>10.4f}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP S: S_D sweep in the ORIGINAL lattice")
    print("-" * 78)
    print("2026-07-29 T5 claimed 'no choice of S_D removes the trivial vector -")
    print("both vectors scale linearly in S_D'.  Only the trivial one does:")
    print("  ||n*S_D*e_d|| = n*S_D          (linear in S_D)")
    print("  ||v_planted||^2 = R^2 + d^2*S_D^2   with R independent of S_D")
    print("so trivial > planted iff S_D > R/sqrt(n^2-d^2); the ratio saturates at n/d.")
    print()
    for label, cv, k1, m, d_sec in curves:
        p, b, n, lam, G = cv
        K2 = math.isqrt(n) + 1
        S_K1, _, S_K2, S_KAN = scales(n, k1, K2)
        R2 = m * (k1 * S_K1) ** 2 / 3.0 + m * (K2 * S_K2) ** 2 / 3.0 + S_KAN ** 2
        sd_pred = math.sqrt(R2) / math.sqrt(float(n) ** 2 - float(d_sec) ** 2)
        print(f"\n  {label}  predicted crossover S_D > {sd_pred:.3f}   "
              f"saturation n/d = {n/d_sec:.3f}")
        print(f"    {'S_D':>6} {'triv/pv':>9} {'sv/pv':>8} {'sv is triv':>11} {'recovered':>10}")
        for S_D in (1, 2, 3, 4, 6, 8, 16, 64, 256):
            hits = 0
            acc_t, acc_s, triv_short = 0.0, 0.0, 0
            for s in SEEDS:
                r = run_orig(cv, m, d_sec, k1, seed=s, S_D=S_D)
                hits += 1 if r['ok'] else 0
                acc_t += r['triv'] / r['pv']
                acc_s += r['sv_over_pv']
                if abs(r['sv'] - r['triv']) < 1e-6:
                    triv_short += 1
            print(f"    {S_D:>6} {acc_t/len(SEEDS):>9.4f} {acc_s/len(SEEDS):>8.4f} "
                  f"{str(triv_short)+'/5':>11} {str(hits)+'/5':>10}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
