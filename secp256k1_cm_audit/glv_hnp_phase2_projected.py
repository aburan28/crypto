"""
GLV-HNP Phase 2, Thread 23: remove the trivial vector n*S_D*e_m from the
Phase-2 lattice and test whether the K1 wall moves.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, finding T5):
  The Phase-2 lattice L (dim 2m+2, built in glv_hnp_phase2_20bit.py:262)
  always has lambda_1(L) = n*S_D*e_m -- the "d is only defined mod n"
  vector.  It is 2-3x shorter than the planted vector on every curve
  tested, success and failure alike, so the planted vector is NEVER
  lambda_1 and recovery is a BDD/coset condition, not SVP.

  The 2026-07-29 next-step proposal (Thread 23) was: quotient out that
  direction and re-run the T4 K1 grid.  Falsifier as stated there:
    - if sv/pv rises above 1 AND the K1 wall on the lam*=0.07 curve
      (currently K1 ~ 4-6) moves outward, the reformulation is a real
      improvement;
    - if the wall stays put, the wall is information-theoretic and
      Phase 2 is at its ceiling.

Construction.
  L has coordinates  (k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan).
  Let pi drop the d coordinate.  ker(pi|_L) = Z*(n*S_D*e_m) exactly, so
  L' = pi(L) has rank 2m+1 and the trivial vector is gone by construction.

  Explicit basis of L' (cols: 0..m-1 = k1, m..2m-1 = k2, 2m = kannan).
  With C_i = B_i * B_0^{-1} mod n  (C_0 = 1):

    u_0     = S_K1 * (1, C_1, ..., C_{m-1} | 0 | 0)
    u_i     = n*S_K1 * e_i                         i = 1..m-1
    w_j     = -lam*S_K1 * e_j  +  S_K2 * e_{m+j}   j = 0..m-1
    z       = S_K1 * (A_0, ..., A_{m-1} | 0 | 0) + S_KANNAN * e_{2m}

  det(L') = S_K1^m * n^(m-1) * S_K2^m * S_KANNAN = det(L) / (n*S_D),
  which is the covolume identity for a projection along a primitive
  lattice direction of length n*S_D -- an independent check that ker is
  exactly what we claimed.

  This is the same lattice one gets by eliminating d algebraically before
  building the lattice (substitute d = (k_0 - A_0)/B_0 into the other m-1
  congruences); the determinants and dimensions coincide.

  Planted vector in L':  (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN), reached with
  the u_0 coefficient t = d*B_0 mod n.  Recovery inverts that:
  d = (k1_i + lam*k2_i - A_i) * B_i^{-1} mod n for any index i.

Experiments:
  E1  structural checks: planted vector is in L'; det identity; the trivial
      vector is present in L and dies under pi; Gaussian-heuristic ratios.
  E2  T4 grid replication on the historical anchors, K1 in
      {2,3,4,6,8,12,16,24}, lattice A (baseline) vs B (projected).
  E3  fresh 17-bit curves, m-sweep at fixed eff, A vs B, first-all-seeds m.
  E4  per-trial agreement of A and B.
  E5  separator head-to-head: sv/pv(A), sv/pv(B), nu_hat, lam*.
  E6  mechanism: is lambda_1(L') the lambda-block vector nu_hat measures?
  E7  scaling law for the nu_hat threshold (this script's own hypothesis).

RESULTS (2026-08-03 run; full output in
glv_hnp_phase2_projected_output.txt):

  E1  The projection is exact: det(L')*n*S_D == det(L) as exact integers on
      all three anchors, and the planted vector lies in L'.  sv/pv does rise
      (0.466 -> 0.591 and 0.468 -> 0.813 at m=6, K1=4) but not above 1 in
      general, so the falsifier's first clause already fails.

  E2/E3/E4  THREAD 23 FALSIFIED, second branch.  The K1 wall does not move
      at all: the lam*=0.070 anchor still walls at K1 ~ 4-6 and the
      lam*=0.340 anchor at K1 ~ 16, in BOTH lattices.  Every cell of the
      8-point K1 grid is identical; E3 totals 104/400 vs 104/400 with the
      same first-all m on all 10 curves; E4 finds 336/336 per-trial
      agreement and ZERO disagreements.  The trivial vector n*S_D*e_m was
      never an obstruction -- it is a spectator that LLL parks in one basis
      slot while recover_* scans all rows.  This is exactly what the
      2026-07-29 T5 analysis already said (recovery lives "in the
      (2m+1)-dimensional projection along e_m"); L' IS that projection, so
      A and B are literally the same problem.  The wall is
      information-theoretic and Phase 2 is at its ceiling.

  E5  What the projection DOES buy is a diagnostic.  On 40 fresh 17-bit
      curves (eff=0.10, m=8) sv/pv(B) separates C1 from C2 with AUC 0.987,
      while sv/pv(A) is a coin flip at 0.513 -- in L the statistic is pinned
      by the trivial vector and carries no signal.  nu_hat matches
      sv/pv(B) exactly (AUC 0.987) at O(log n) instead of a full LLL, so
      sv/pv(B) is strictly dominated as a practical predictor.  lam* 0.714.
      This also replicates the 2026-07-29 nu_hat result on fresh curves at a
      new bit size.

  E6  MECHANISM for nu_hat.  lambda_1(L') is the lambda-block vector mu
      exactly when the attack succeeds -- sv/pv(B) equals mu/||pv|| to 4
      decimals and the vector occupies 1 of the m lambda-blocks.  When the
      attack fails, mu is no longer shortest and lambda_1(L') is a diffuse
      combination spanning all m blocks.  So nu_hat is not merely correlated
      with success; it measures which of the two regimes the lattice is in.

  E7  This script's own scaling law is FALSIFIED.  The E6 mechanism implies
      nu_hat_crit = c*sqrt(eff)*sqrt(2m/3+1) with c universal.  Fitted
      independently over a 4x4 (eff, m) grid on 30 curves, c spreads 110.8%
      -- no better than nu_hat_crit's own 116.6%.  In a like-for-like
      comparison (one global free parameter each) the scaling rule scores
      0.9481 vs 0.9407 mean accuracy on the 9 well-separated settings and
      LOSES on all 14 (0.8381 vs 0.8571).  At fixed m, nu_hat_crit is
      roughly constant as eff varies 0.068 -> 0.129 (m=8: 0.482, 0.508,
      0.468, 0.474) where the law predicts 0.444 -> 0.613.
      Practical takeaway: a single fixed cut nu_hat <= 0.62 transfers across
      eff in [0.068, 0.129] and m in {6,8,10,12} without recalibration.

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sys

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Arithmetic helpers (verbatim from glv_hnp_phase2_lambda_threshold.py)
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a % m, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        if y1 == 0: return None
        lam = 3 * x1 * x1 % p * modinv(2 * y1, p) % p
    else:
        lam = (y2 - y1) % p * modinv((x2 - x1) % p, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    R = None
    while k > 0:
        if k & 1: R = ec_add(R, P, p)
        P = ec_add(P, P, p)
        k >>= 1
    return R

def tonelli_shanks(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
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

# ---------------------------------------------------------------------------
# j=0 curve construction
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0: break
        s = math.isqrt(disc)
        if s * s != disc: continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

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
    return min(lam, n - lam) / n

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                lam = glv_roots(n)[0]
                return (p, b, n, lam, G)
    return None

def find_curve_for_p(p, seed=12345):
    """All j=0 GLV curves over F_p with prime n = 1 mod 3."""
    if p % 3 != 1: return []
    eis = eisenstein_decompose(p)
    if eis is None: return []
    out = []
    for t in j0_traces(*eis):
        n = p + 1 - t
        if n < 5 or n % 3 != 1 or not sympy.isprime(n): continue
        if glv_roots(n) is None: continue
        cur = build_curve(p, n, seed)
        if cur is not None:
            out.append(cur)
    return out

def search_curves(lo, hi, want, seed=12345):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        for cur in find_curve_for_p(p, seed):
            out.append(cur)
            if len(out) >= want: break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Signatures + scales (verbatim)
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

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Lattice A: baseline, dim 2m+2 (verbatim from glv_hnp_phase2_20bit.py:262)
# ---------------------------------------------------------------------------

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound):
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

def planted_A(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_A(reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Lattice B: projected along e_m, dim 2m+1 (Thread 23)
# ---------------------------------------------------------------------------

def build_lattice_B(sigs, n, lam, k1_bound, k2_bound):
    """pi(L) with pi dropping the d coordinate.  Cols: 0..m-1 k1,
    m..2m-1 k2, 2m kannan."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    C = [sigs[i]['B'] * B0inv % n for i in range(m)]   # C[0] == 1
    assert C[0] == 1

    M = [[0] * dim for _ in range(dim)]
    # u_0
    for i in range(m):
        M[0][i] = C[i] * S_K1
    # u_1 .. u_{m-1}
    for i in range(1, m):
        M[i][i] = n * S_K1
    # w_0 .. w_{m-1}
    for j in range(m):
        M[m + j][j] = -lam * S_K1
        M[m + j][m + j] = S_K2
    # z
    for i in range(m):
        M[2 * m][i] = sigs[i]['A'] * S_K1
    M[2 * m][2 * m] = S_KANNAN
    return M

def planted_B(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_B(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_i, k2_i) off any row with |last| = S_KANNAN and invert the
    signature equation:  d = (k1_i + lam*k2_i - A_i) * B_i^{-1} mod n."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 or b % S_K2: continue
            k1i, k2i = a // S_K1, b // S_K2
            d_cand = (k1i + lam * k2i - sigs[i]['A']) % n * \
                     modinv(sigs[i]['B'], n) % n
            if d_cand and d_cand == d_secret:
                return True
    return False

# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_once(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False,
             bkz_beta=20):
    """Returns (recovered, planted_norm, shortest_reduced_norm) or None."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m: return None
    _S_K1, _S_D, _S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)

    if variant == 'A':
        M = build_lattice_A(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_A(sigs, d_secret, n, k1_bound, k2_bound)
    else:
        M = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_B(sigs, n, k1_bound, k2_bound)
    dim = len(M)
    X = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(X, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(X)
    reduced = [[X[i][j] for j in range(dim)] for i in range(dim)]

    if variant == 'A':
        ok = recover_A(reduced, m, n, S_KANNAN, d_secret)
    else:
        ok = recover_B(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret)
    nz = [r for r in reduced if any(r)]
    sn = min(norm(r) for r in nz)
    return (ok, norm(pv), sn)

def success_rate(curve, m, k1_bound, seeds, variant, use_bkz=False,
                 bkz_beta=20):
    n = curve[2]
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = run_once(curve, m, d_trial, k1_bound, seed, variant,
                       use_bkz, bkz_beta)
        if res is None: continue
        ok, pn, sn = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))

def gh_ratio(m, n, k1_bound, k2_bound, variant):
    """planted_norm / gaussian_heuristic(det, dim), from the exact scales."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if variant == 'A':
        dim = 2 * m + 2
        logdet = m * math.log(n * S_K1) + math.log(S_D) + \
                 m * math.log(S_K2) + math.log(S_KAN)
        pn2 = m * (k1_bound * S_K1) ** 2 / 3 + (n * S_D) ** 2 / 3 + \
              m * (k2_bound * S_K2) ** 2 / 3 + S_KAN ** 2
    else:
        dim = 2 * m + 1
        logdet = m * math.log(S_K1) + (m - 1) * math.log(n) + \
                 m * math.log(S_K2) + math.log(S_KAN)
        pn2 = m * (k1_bound * S_K1) ** 2 / 3 + \
              m * (k2_bound * S_K2) ** 2 / 3 + S_KAN ** 2
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)
    return math.sqrt(pn2) / gh

# ---------------------------------------------------------------------------
# E1  structural checks
# ---------------------------------------------------------------------------

def exp1(anchors):
    print("=" * 78)
    print("E1  Structural checks on the projected lattice L'")
    print("=" * 78)
    for label, curve in anchors:
        p, b, n, lam, G = curve
        m, k1_bound = 6, 4
        k2_bound = math.isqrt(n) + 1
        d = random.Random(1).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, 11)
        if len(sigs) < m:
            print(f"  {label}: not enough signatures"); continue
        S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

        MA = build_lattice_A(sigs, n, lam, k1_bound, k2_bound)
        MB = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
        pA = planted_A(sigs, d, n, k1_bound, k2_bound)
        pB = planted_B(sigs, n, k1_bound, k2_bound)

        # (a) both planted vectors lie in their lattices (exact integer solve)
        inA = in_lattice(MA, pA)
        inB = in_lattice(MB, pB)

        # (b) determinant identity det(L') == det(L)/(n*S_D)
        dA = abs(sympy.Matrix(MA).det())
        dB = abs(sympy.Matrix(MB).det())
        ident = (dA == dB * n * S_D)

        # (c) trivial vector: present in A, and its image is 0 (absent in B)
        triv = [0] * (2 * m + 2); triv[m] = n * S_D
        triv_in_A = in_lattice(MA, triv)

        # (d) lambda_1 after LLL
        okA, pnA, snA = run_once(curve, m, d, k1_bound, 11, 'A')
        okB, pnB, snB = run_once(curve, m, d, k1_bound, 11, 'B')

        print(f"\n  {label}: p={p} n={n} lam={lam} lam*={lam_star(lam,n):.4f} "
              f"m={m} K1={k1_bound}")
        print(f"    planted in L  : {inA}          planted in L' : {inB}")
        print(f"    det(L)        : {dA}")
        print(f"    det(L')*n*S_D : {dB * n * S_D}   identity holds: {ident}")
        print(f"    n*S_D*e_m in L: {triv_in_A}  (image under pi is 0 by "
              f"construction)")
        print(f"    LLL  sv/pv    : A {snA/pnA:.4f}   B {snB/pnB:.4f}"
              f"     recovered: A {okA}  B {okB}")
        print(f"    planted/GH    : A {gh_ratio(m,n,k1_bound,k2_bound,'A'):.4f}"
              f"   B {gh_ratio(m,n,k1_bound,k2_bound,'B'):.4f}")
        print(f"    lambda_1(L') energy: {sv_profile_B(curve, m, d, k1_bound, 11)}")

def sv_profile_B(curve, m, d, k1_bound, seed):
    """Where does the energy of lambda_1(L') sit?  Returns the fraction of
    ||sv||^2 in the k1 block / k2 block / kannan coord, plus how many of the
    m 2-D lambda-blocks (col j, col m+j) it occupies."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed)
    M = build_lattice_B(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    X = IntegerMatrix.from_matrix(M)
    LLL.reduction(X)
    rows = [[X[i][j] for j in range(dim)] for i in range(dim) if any(X[i])]
    sv = min(rows, key=norm)
    tot = sum(x * x for x in sv)
    k1e = sum(sv[i] ** 2 for i in range(m)) / tot
    k2e = sum(sv[m + i] ** 2 for i in range(m)) / tot
    kan = sv[2 * m] ** 2 / tot
    occ = sum(1 for j in range(m) if sv[j] or sv[m + j])
    return (f"k1 {k1e:.3f}  k2 {k2e:.3f}  kannan {kan:.3f}  "
            f"lambda-blocks occupied {occ}/{m}")

def in_lattice(M, v):
    """Is v in the row lattice of M?  Exact rational solve + integrality."""
    A = sympy.Matrix(M).T          # columns are the basis vectors
    b = sympy.Matrix(v)
    try:
        x = A.solve_least_squares(b) if A.rows != A.cols else A.solve(b)
    except Exception:
        try:
            x = A.gauss_jordan_solve(b)[0]
        except Exception:
            return None
    if A * x != b:
        return False
    return all(c.is_integer for c in x)

# ---------------------------------------------------------------------------
# E2  T4 grid replication, A vs B
# ---------------------------------------------------------------------------

def exp2(anchors, k1_grid, m, seeds):
    print()
    print("=" * 78)
    print(f"E2  T4 K1 grid replication (m={m}, {len(seeds)} seeds)  A = baseline"
          f", B = projected")
    print("=" * 78)
    rows = []
    for label, curve in anchors:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        print(f"\n  {label}: p={p} n={n} ({n.bit_length()}b) lam={lam} "
              f"lam*={lam_star(lam,n):.4f}")
        hdr = "    K1      eff  |  A wins  B wins  |  A sv/pv  B sv/pv  |  " \
              "A pl/GH  B pl/GH"
        print(hdr)
        for K1 in k1_grid:
            eff = K1 * k2_bound / n
            wA, tA, rA = success_rate(curve, m, K1, seeds, 'A')
            wB, tB, rB = success_rate(curve, m, K1, seeds, 'B')
            gA = gh_ratio(m, n, K1, k2_bound, 'A')
            gB = gh_ratio(m, n, K1, k2_bound, 'B')
            print(f"    {K1:<4} {eff:7.4f}  |  {wA}/{tA}     {wB}/{tB}     |  "
                  f"{rA:7.4f}  {rB:7.4f}  |  {gA:7.4f}  {gB:7.4f}")
            rows.append((label, K1, wA, wB, tA))
    return rows

# ---------------------------------------------------------------------------
# E3  fresh curves, m-sweep at fixed eff
# ---------------------------------------------------------------------------

def exp3(curves, eff_target, m_range, seeds):
    print()
    print("=" * 78)
    print(f"E3  fresh curves, eff~{eff_target}, first m reaching "
          f"{len(seeds)}/{len(seeds)}   A vs B")
    print("=" * 78)
    print("    curve                lam*    K1   |  first-all m  |  total wins"
          "  |  sv/pv")
    print("                                      |    A     B    |    A     B "
          "  |   A      B")
    agg = {'A': 0, 'B': 0, 'tot': 0}
    firsts = []
    for (p, b, n, lam, G) in curves:
        k2_bound = math.isqrt(n) + 1
        K1 = max(2, int(eff_target * n / k2_bound))
        fA = fB = None
        totA = totB = tot = 0
        rAs, rBs = [], []
        for m in m_range:
            wA, tA, rA = success_rate((p, b, n, lam, G), m, K1, seeds, 'A')
            wB, tB, rB = success_rate((p, b, n, lam, G), m, K1, seeds, 'B')
            totA += wA; totB += wB; tot += tA
            rAs.append(rA); rBs.append(rB)
            if fA is None and wA == tA: fA = m
            if fB is None and wB == tB: fB = m
        agg['A'] += totA; agg['B'] += totB; agg['tot'] += tot
        firsts.append((fA, fB))
        mA = sum(rAs) / len(rAs); mB = sum(rBs) / len(rBs)
        print(f"    p={p:<8} n={n:<9} {lam_star(lam,n):.4f}  {K1:<4} |  "
              f"{str(fA):>4}  {str(fB):>4}  |  {totA:>3}   {totB:>3}   |  "
              f"{mA:.3f}  {mB:.3f}")
    print(f"\n    TOTAL wins:  A {agg['A']}/{agg['tot']}   "
          f"B {agg['B']}/{agg['tot']}")
    both = [(a, bb) for a, bb in firsts if a is not None and bb is not None]
    if both:
        print(f"    first-all m (curves where both succeed, n={len(both)}): "
              f"mean A {sum(a for a,_ in both)/len(both):.2f}  "
              f"B {sum(b for _,b in both)/len(both):.2f}")
    print(f"    B better / same / worse on first-all m: "
          f"{sum(1 for a,bb in both if bb<a)} / "
          f"{sum(1 for a,bb in both if bb==a)} / "
          f"{sum(1 for a,bb in both if bb>a)}")
    return agg, firsts

# ---------------------------------------------------------------------------
# E4  per-trial agreement of A and B
# ---------------------------------------------------------------------------

def exp4(curves, k1_grid, m_range, seeds):
    print()
    print("=" * 78)
    print("E4  per-trial agreement between lattice A and lattice B")
    print("=" * 78)
    both = agree = disagree = 0
    dis_examples = []
    for (p, b, n, lam, G) in curves:
        k2_bound = math.isqrt(n) + 1
        for K1 in k1_grid:
            for m in m_range:
                for seed in seeds:
                    d = random.Random(seed + 7777).randint(1, n - 1)
                    rA = run_once((p, b, n, lam, G), m, d, K1, seed, 'A')
                    rB = run_once((p, b, n, lam, G), m, d, K1, seed, 'B')
                    if rA is None or rB is None: continue
                    both += 1
                    if rA[0] == rB[0]:
                        agree += 1
                    else:
                        disagree += 1
                        if len(dis_examples) < 10:
                            dis_examples.append((p, n, K1, m, seed, rA[0], rB[0]))
    print(f"    trials compared : {both}")
    print(f"    identical outcome: {agree}  ({100.0*agree/both:.2f}%)")
    print(f"    disagreements    : {disagree}")
    for e in dis_examples:
        print(f"      p={e[0]} n={e[1]} K1={e[2]} m={e[3]} seed={e[4]} "
              f"A={e[5]} B={e[6]}")
    return agree, disagree, both

# ---------------------------------------------------------------------------
# E5  separator head-to-head: sv/pv(A), sv/pv(B), nu_hat, lam*
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Lagrange-Gauss reduction; returns the shorter basis vector's norm."""
    def dot(a, c): return a[0] * c[0] + a[1] * c[1]
    while True:
        if dot(v, v) < dot(u, u): u, v = v, u
        mu = dot(u, v) / dot(u, u)
        r = round(mu)
        v = (v[0] - r * u[0], v[1] - r * u[1])
        if dot(v, v) >= dot(u, u):
            return math.sqrt(dot(u, u))

def nu_hat(n, lam, S_K1, S_K2):
    """lambda_1(L2)/sqrt(det L2) for the 2-D lambda block
    L2 = <(n*S_K1,0), (-lam*S_K1, S_K2)>.  O(log n), no LLL."""
    mu = gauss_reduce_2d((n * S_K1, 0), (-lam * S_K1, S_K2))
    return mu / math.sqrt(n * S_K1 * S_K2)

def auc(pos, neg):
    """P(score(pos) > score(neg)), ties at 0.5."""
    if not pos or not neg: return float('nan')
    s = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in pos for b in neg)
    return s / (len(pos) * len(neg))

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v); i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]: j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1): r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    nn = len(xs)
    mx, my = sum(rx) / nn, sum(ry) / nn
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) *
                    sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')

def exp5(curves, eff_target, m, seeds):
    print()
    print("=" * 78)
    print(f"E5  separator head-to-head at eff~{eff_target}, m={m}, "
          f"{len(seeds)} seeds/curve")
    print("=" * 78)
    print("    curve                   lam*   nu_hat  svpv_A  svpv_B  rate")
    rows = []
    for (p, b, n, lam, G) in curves:
        k2_bound = math.isqrt(n) + 1
        K1 = max(2, int(eff_target * n / k2_bound))
        S_K1, _S_D, S_K2, _S = scales(n, K1, k2_bound)
        wA, tA, rA = success_rate((p, b, n, lam, G), m, K1, seeds, 'A')
        wB, tB, rB = success_rate((p, b, n, lam, G), m, K1, seeds, 'B')
        assert wA == wB, (p, n, wA, wB)
        nh = nu_hat(n, lam, S_K1, S_K2)
        ls = lam_star(lam, n)
        rate = wA / tA
        rows.append(dict(p=p, n=n, lam_star=ls, nu=nh, svA=rA, svB=rB,
                         rate=rate))
        print(f"    p={p:<8} n={n:<9} {ls:.4f}  {nh:.4f}  {rA:.4f}  {rB:.4f}"
              f"  {wA}/{tA}")

    print(f"\n    curves: {len(rows)}   "
          f"C1 (rate>=0.75): {sum(1 for r in rows if r['rate']>=0.75)}   "
          f"C2 (rate<=0.25): {sum(1 for r in rows if r['rate']<=0.25)}   "
          f"mid: {sum(1 for r in rows if 0.25<r['rate']<0.75)}")
    print("\n    predictor   AUC(C1 vs C2)   spearman(pred, rate)   "
          "best single cut (acc)")
    for key, name in (('svB', 'sv/pv (B)'), ('svA', 'sv/pv (A)'),
                      ('nu', 'nu_hat   '), ('lam_star', 'lam*     ')):
        pos = [r[key] for r in rows if r['rate'] >= 0.75]
        neg = [r[key] for r in rows if r['rate'] <= 0.25]
        a = auc(pos, neg)
        # report the orientation-free AUC (max of a, 1-a) plus direction
        direction = 'low=easy' if a < 0.5 else 'high=easy'
        a_or = max(a, 1 - a)
        sp = spearman([r[key] for r in rows], [r['rate'] for r in rows])
        acc = best_cut([r[key] for r in rows],
                       [1 if r['rate'] >= 0.75 else 0 for r in rows])
        print(f"    {name}   {a_or:.3f}  ({direction})      {sp:+.3f}"
              f"               {acc:.3f}")
    return rows

def exp6(curves, eff_target, m, seeds):
    """Is lambda_1(L') the lambda-block vector that nu_hat measures?

    L2_j = <(n*S_K1) e_j, -lam*S_K1 e_j + S_K2 e_{m+j}> sits inside L' for
    every j, so mu = lambda_1(L2) is an upper bound on lambda_1(L').  If
    lambda_1(L') == mu whenever the attack fails, and == ||planted|| when it
    succeeds, then nu_hat is not just correlated with success -- it is the
    mechanism.
    """
    print()
    print("=" * 78)
    print(f"E6  mechanism: is lambda_1(L') the lambda-block vector?  "
          f"(eff~{eff_target}, m={m})")
    print("=" * 78)
    print("    curve                nu_hat  |  sv/pv_B  mu/pv  pv/pv  |  "
          "sv is  blocks  rate")
    tally = {}
    for (p, b, n, lam, G) in curves:
        k2_bound = math.isqrt(n) + 1
        K1 = max(2, int(eff_target * n / k2_bound))
        S_K1, _S_D, S_K2, S_KAN = scales(n, K1, k2_bound)
        nh = nu_hat(n, lam, S_K1, S_K2)
        mu = gauss_reduce_2d((n * S_K1, 0), (-lam * S_K1, S_K2))
        wins, svpv, mupv, kinds, blocks = 0, [], [], [], []
        for seed in seeds:
            d = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d, m, n, lam, p, K1, k2_bound, seed)
            if len(sigs) < m: continue
            M = build_lattice_B(sigs, n, lam, K1, k2_bound)
            dim = 2 * m + 1
            X = IntegerMatrix.from_matrix(M)
            LLL.reduction(X)
            rows = [[X[i][j] for j in range(dim)] for i in range(dim)
                    if any(X[i])]
            sv = min(rows, key=norm)
            pv = planted_B(sigs, n, K1, k2_bound)
            pn = norm(pv)
            svpv.append(norm(sv) / pn); mupv.append(mu / pn)
            occ = sum(1 for j in range(m) if sv[j] or sv[m + j])
            blocks.append(occ)
            if abs(norm(sv) - mu) < 1e-6 * mu:
                kinds.append('mu')
            elif sv[2 * m] and abs(abs(sv[2 * m]) - S_KAN) < 1e-9:
                kinds.append('planted')
            else:
                kinds.append('other')
            wins += recover_B(rows, sigs, n, lam, K1, k2_bound, d)
        if not svpv: continue
        kind = max(set(kinds), key=kinds.count)
        tally[kind] = tally.get(kind, 0) + len(kinds)
        for k in kinds: tally[k] = tally.get(k, 0)
        print(f"    p={p:<8} n={n:<9} {nh:.4f}  |  "
              f"{sum(svpv)/len(svpv):.4f}   {sum(mupv)/len(mupv):.4f}  "
              f"1.0000 |  {kind:<7} {sum(blocks)/len(blocks):4.1f}   "
              f"{wins}/{len(seeds)}")
    print()
    return tally

def exp7(curves, settings, seeds):
    """Does the nu_hat threshold obey the scaling law the E6 mechanism implies?

    E6 shows lambda_1(L') = mu (the lambda-block vector) exactly when the
    attack succeeds, and a diffuse all-block vector when it fails.  So the
    real decision variable is  mu / ||v_planted||, and

        mu            = nu_hat * sqrt(n*S_K1*S_K2) = nu_hat * n / sqrt(eff)
        ||v_planted|| ~ n * sqrt(2m/3 + 1)

        =>  mu/||pv||  =  nu_hat / ( sqrt(eff) * sqrt(2m/3 + 1) )

    PREDICTION: the critical mu/||pv|| is a universal constant c, so the
    critical nu_hat MOVES with (eff, m) as

        nu_hat_crit(eff, m) = c * sqrt(eff) * sqrt(2m/3 + 1).

    FALSIFIER: fit c independently at each (eff, m).  If the mechanism is
    right, c is constant across settings while nu_hat_crit is not.  If
    nu_hat_crit is instead the constant, the mechanism is wrong and nu_hat
    is a curve invariant after all.
    """
    print()
    print("=" * 78)
    print("E7  scaling law for the nu_hat threshold")
    print("=" * 78)
    print("     eff    m  | curves C1/C2 |  nu_hat_crit  acc  |  mu/pv_crit"
          "   c_fit  |  predicted nu_crit (c=const)")
    out = []
    cache = []
    for eff_target, m in settings:
        rows = []
        for (p, b, n, lam, G) in curves:
            k2_bound = math.isqrt(n) + 1
            K1 = max(2, int(eff_target * n / k2_bound))
            eff = K1 * k2_bound / n
            S_K1, _S_D, S_K2, _S = scales(n, K1, k2_bound)
            nh = nu_hat(n, lam, S_K1, S_K2)
            mu = gauss_reduce_2d((n * S_K1, 0), (-lam * S_K1, S_K2))
            wins = 0
            pvs = []
            for seed in seeds:
                d = random.Random(seed + 7777).randint(1, n - 1)
                r = run_once((p, b, n, lam, G), m, d, K1, seed, 'B')
                if r is None: continue
                wins += bool(r[0]); pvs.append(r[1])
            if not pvs: continue
            rows.append(dict(nu=nh, mupv=mu / (sum(pvs) / len(pvs)),
                             rate=wins / len(pvs), eff=eff))
        lab = [1 if r['rate'] >= 0.5 else 0 for r in rows]
        if min(sum(lab), len(lab) - sum(lab)) < 3:
            print(f"    {eff_target:5.3f} {m:3}  |  {sum(lab):3}/"
                  f"{len(lab)-sum(lab):<3}    | outside transition band, "
                  f"not fitted")
            continue
        nu_c, acc = best_cut_threshold([r['nu'] for r in rows], lab)
        mu_c, _ = best_cut_threshold([r['mupv'] for r in rows], lab)
        eff_mean = sum(r['eff'] for r in rows) / len(rows)
        denom = math.sqrt(eff_mean) * math.sqrt(2 * m / 3 + 1)
        c_fit = nu_c / denom
        out.append((eff_target, m, nu_c, mu_c, c_fit, denom, acc,
                    sum(lab), len(lab) - sum(lab)))
        cache.append(dict(eff=eff_mean, m=m, denom=denom, acc=acc,
                          nus=[r['nu'] for r in rows],
                          mupvs=[r['mupv'] for r in rows], lab=lab))
        print(f"    {eff_mean:5.3f} {m:3}  |  {sum(lab):3}/{len(lab)-sum(lab):<3}    "
              f"|    {nu_c:.4f}    {acc:.2f} |    {mu_c:.4f}   {c_fit:.4f}  |")
    if out:
        cs = [o[4] for o in out]
        mus = [o[3] for o in out]
        cbar = sum(cs) / len(cs)
        print(f"\n    c_fit across settings: {[f'{c:.4f}' for c in cs]}")
        print(f"    mean c = {cbar:.4f}   spread = "
              f"{(max(cs)-min(cs))/cbar*100:.1f}% of mean")
        print(f"    mu/pv_crit across settings: {[f'{v:.4f}' for v in mus]}"
              f"   spread = {(max(mus)-min(mus))/(sum(mus)/len(mus))*100:.1f}%")
        print(f"    nu_hat_crit across settings: "
              f"{[f'{o[2]:.4f}' for o in out]}   spread = "
              f"{(max(o[2] for o in out)-min(o[2] for o in out))/(sum(o[2] for o in out)/len(out))*100:.1f}%")
        print(f"\n     eff    m  | observed nu_crit | predicted (c={cbar:.3f}) "
              "| rel err")
        for (e, m, nu_c, mu_c, c_fit, denom, acc, np_, nn_) in out:
            pred = cbar * denom
            print(f"    {e:5.3f} {m:3}  |      {nu_c:.4f}      |     "
                  f"{pred:.4f}       | {100*(pred-nu_c)/nu_c:+6.1f}%")

        # -- fair head-to-head between the two competing rules ---------------
        # R1  constant cut:      nu_hat <= t          (t fitted globally)
        # R2  scaling-law cut:   nu_hat <= c*denom    (c fitted globally)
        # Both have exactly one global free parameter, fitted on the SAME
        # pooled data, so the comparison is like-for-like.
        clean = [s for s in cache if s['acc'] >= 0.95]
        print(f"\n    Rule comparison on the {len(clean)} well-separated "
              f"settings (cut accuracy >= 0.95):")
        for pool, tag in ((clean, 'well-separated'), (cache, 'all fitted')):
            if not pool: continue
            def mean_acc(score_fn, param):
                accs = []
                for s in pool:
                    pred = [1 if score_fn(s, i) <= param else 0
                            for i in range(len(s['lab']))]
                    accs.append(sum(1 for a, b in zip(pred, s['lab'])
                                    if a == b) / len(s['lab']))
                return sum(accs) / len(accs)
            r1 = max(((mean_acc(lambda s, i: s['nus'][i], t), t)
                      for t in [x / 200 for x in range(1, 201)]))
            r2 = max(((mean_acc(lambda s, i: s['nus'][i] / s['denom'], c), c)
                      for c in [x / 200 for x in range(1, 401)]))
            r3 = max(((mean_acc(lambda s, i: s['mupvs'][i], t), t)
                      for t in [x / 200 for x in range(1, 401)]))
            print(f"      [{tag}, {len(pool)} settings]")
            print(f"        R1 constant  nu_hat <= {r1[1]:.3f}          "
                  f"mean acc {r1[0]:.4f}")
            print(f"        R2 scaling   nu_hat <= {r2[1]:.3f}*denom    "
                  f"mean acc {r2[0]:.4f}")
            print(f"        R3 direct    mu/||pv|| <= {r3[1]:.3f}       "
                  f"mean acc {r3[0]:.4f}")
    return out

def best_cut_threshold(xs, labels):
    """Best single low-is-positive threshold, returns (threshold, accuracy)."""
    best_t, best_acc = None, -1.0
    for t in sorted(set(xs)):
        pred = [1 if x <= t else 0 for x in xs]
        acc = sum(1 for a, b in zip(pred, labels) if a == b) / len(labels)
        if acc > best_acc:
            best_acc, best_t = acc, t
    return best_t, best_acc

def best_cut(xs, labels):
    """Best single-threshold accuracy, either orientation."""
    best = 0.0
    for t in sorted(set(xs)):
        for sign in (1, -1):
            pred = [1 if sign * x <= sign * t else 0 for x in xs]
            acc = sum(1 for a, b in zip(pred, labels) if a == b) / len(labels)
            best = max(best, acc)
    return best

# ---------------------------------------------------------------------------

def main():
    random.seed(20260803)

    # historical anchors: the two 12-bit curves of the 2026-07-29 T4 K1 grid
    # (lam* = 0.340 and 0.070), plus the 8-bit curve of the T5 table.
    want = {(199, 211), (2557, 2659), (2677, 2647)}
    anchors = []
    for p_anchor in (199, 2557, 2677):
        for cur in find_curve_for_p(p_anchor):
            if (cur[0], cur[2]) in want:
                anchors.append((f"{p_anchor}/n={cur[2]}", cur))
    print("anchors: " + ", ".join(
        f"{l} lam*={lam_star(c[3], c[2]):.4f}" for l, c in anchors) + "\n")

    exp1(anchors)
    exp2(anchors, [2, 3, 4, 6, 8, 12, 16, 24], m=12, seeds=list(range(5)))

    curves17 = search_curves(2 ** 16, 2 ** 17, want=10)
    print(f"\n17-bit curves found: {len(curves17)}")
    exp3(curves17, eff_target=0.10, m_range=range(5, 13), seeds=list(range(5)))

    exp4([c for _, c in anchors] + curves17[:4], [3, 6, 12],
         range(6, 10), list(range(4)))

    curves17b = search_curves(2 ** 16, 2 ** 17, want=40)
    print(f"\n17-bit curves for E5: {len(curves17b)}")
    exp5(curves17b, eff_target=0.10, m=8, seeds=list(range(8)))

    exp6(curves17b[:14], eff_target=0.10, m=8, seeds=list(range(6)))

    exp7(curves17b[:30],
         [(e, m) for m in (6, 8, 10, 12) for e in (0.07, 0.09, 0.11, 0.13)],
         list(range(5)))

if __name__ == '__main__':
    main()
