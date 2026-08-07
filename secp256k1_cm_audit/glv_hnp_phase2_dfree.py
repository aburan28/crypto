"""
GLV-HNP Phase 2, Thread 23: remove the trivial vector n*S_D*e_m from the
Phase-2 lattice and re-measure the K1 wall.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, exp T4/T5):
  T5 found that the shortest vector of the Phase-2 lattice
  (glv_hnp_phase2_lambda_threshold.py:254, build_glv_lattice) is ALWAYS the
  trivial vector n*S_D*e_m -- 100% of its energy in the d column,
  |sv[m]|/n = 1.0000 -- and is 2-3x shorter than the planted vector on every
  curve, success and failure alike.  The planted vector is never lambda_1.
  T4 found the real predictor is the bias strength: both 12-bit curves recover
  for K1 <= 4, with the wall at K1 ~ 12-16 (lam*=0.34) / K1 ~ 4-6 (lam*=0.07).

  The 2026-07-29 next-step proposal (Thread 23) was:
    "project the lattice along e_m (quotient out the trivial n*e_m direction)
     and solve BDD in the projection, or replace the Kannan embedding with an
     explicit CVP call.  Falsifier: if sv/pv rises above 1 after the
     reformulation AND the K1 wall moves outward on the lam*=0.07 curve
     (currently K1 ~ 4-6), the reformulation is a real improvement; if the
     wall stays at K1 ~ 4-6, the wall is information-theoretic and Phase 2 is
     at its ceiling."

This script implements that reformulation three ways and runs the falsifier.

  The trivial vector exists because d gets its own column with S_D = 1 while
  d ~ n, so n*S_D*e_m has norm exactly n.  It carries no information: d is
  only defined mod n, so n*e_m is precisely the "d -> d+n" redundancy.  The
  fix is to DROP the d column and recover d algebraically afterwards.

  Dropping the column naively leaves 2m+1 generators spanning a rank-2m
  lattice (the d-row becomes a *rational* combination of the n*e_i rows,
  coefficients B_i/n), so LLL cannot be fed the generators directly.  An
  explicit basis is available without any HNF computation:

    tau_i = B_i * B_0^{-1} mod n      (i = 1..m-1;  B_0 is invertible mod n)

    columns:  0..m-1     k1 columns, scale S_K1
              m..2m-1    k2 columns, scale S_K2
              2m         Kannan column, scale S_KANNAN   (variant D only)

    R_i = (-lam*S_K1*e_i | S_K2*e_i | 0)          i = 0..m-1     [k2 rows]
    Q_0 = (S_K1, tau_1*S_K1, ..., tau_{m-1}*S_K1 | 0 | 0)        [d row]
    Q_i = (n*S_K1*e_i | 0 | 0)                    i = 1..m-1     [mod-n rows]
    K   = (A_i*S_K1 | 0 | S_KANNAN)                              [Kannan row]

  {R_i} u {Q_0, Q_1..Q_{m-1}} is a genuine basis of the rank-2m lattice:
  the vectors of L with all k2 columns zero are exactly
  S_K1 * {x in Z^m : x_i = tau_i * x_0 mod n}, whose determinant is n^(m-1),
  and det L = S_K1^m * S_K2^m * n^(m-1) = det(full-column L)/n, i.e. exactly
  the index-n quotient by the trivial direction.

  Planted vector (variant D):  v = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN)
    v = K + sum_i k2_i R_i + a Q_0 + sum_{i>=1} c_i Q_i
    with a = k1_0 - A_0 + lam*k2_0  and  c_i = (k1_i - A_i + lam*k2_i - a*tau_i)/n
    (integral because B_0*tau_i = B_i mod n).  Verified numerically in E1.

  d is recovered from the k1/k2 entries of the found row:
    d = B_0^{-1} * (k1_0 + lam*k2_0 - A_0)  mod n

Variants compared:
  A  baseline   Kannan lattice with the d column, dim 2m+2   (unchanged)
  D  d-free     Kannan lattice without the d column, dim 2m+1
  E  CVP        d-free lattice without the Kannan row, dim 2m, explicit CVP
                (fplll enumeration, Babai nearest-plane fallback)

Experiments:
  E1  correctness: planted vector really is in the variant-D lattice, and the
      d-free recovery returns the right d on the historical easy curves
  E2  sv/pv for A vs D vs E on the three historical curves + energy profile
  E3  THE FALSIFIER: K1 wall grid, A vs D vs E, 5 seeds, both 12-bit curves
  E4  BKZ(20) on variant D at the wall
  E5  control: does the *classical* Nguyen-Shparlinski HNP lattice also have a
      trivial vector shorter than its planted vector?  (i.e. is T5's finding
      pathological or generic?)
  E6  is the observed K1 wall information-theoretic?  (counting argument)
  E7  does over-determination (large m) rescue D or E at the wall?
  E8  Bleichenbacher nonce bias |B| of the GLV nonce distribution
  E9  falsifier for E8: does |B| predict LLL failure at fixed eff?

RESULTS (2026-08-07 run; full output in glv_hnp_phase2_dfree_output.txt)

  E1  The d-free basis is CORRECT.  The planted vector has integral
      coordinates in it and reconstructs exactly; det(L_A)/S_KANNAN divided by
      det(L_D-no-Kannan) equals n exactly on all three curves (199, 2659,
      2647), confirming L_D is precisely the index-n quotient by the trivial
      direction.  No information is lost: d is recovered algebraically from
      the k1_0/k2_0 entries and matches on every instance where A succeeds.

  E2  Dropping the d column does raise sv/pv (0.603 -> 0.843, 0.517 -> 0.532,
      0.422 -> 0.813) and moves the shortest vector's energy out of the d
      column into the k1/k2 blocks, as intended.  On the lam*=0.07 curve at
      K1 <= 4, sv/pv = 1.000 exactly: the planted vector IS lambda_1 of L_D.
      So the reformulation achieves what Thread 23 asked for.

  E3  THE FALSIFIER RESOLVES AGAINST THE REFORMULATION.  Making the planted
      vector lambda_1 does NOT move the wall.  On lam*=0.34 all three variants
      are bit-identical across K1 = 2..24.  On lam*=0.07 the wall moves only
      from K1=4 to K1=6 (A 0/5, D 2/5, E 4/5 at K1=6) and all three are 0/5 at
      K1=8.  Per the stated falsifier, the reformulation is not a real
      improvement.

  E5  T5's "planted vector is never lambda_1" is GENERIC, NOT A DEFECT.  The
      classical Nguyen-Shparlinski HNP lattice has n*e_d as its exact shortest
      vector too, with sv/pv = 0.11-0.18 -- far WORSE than Phase-2's
      0.34-0.61 -- and it still recovers d (3 of 4 configurations).  This
      corrects the 2026-07-29 reading of T5 and explains E3: the trivial
      vector was never the obstruction, so removing it could not help.

  E6  THE WALL IS NOT INFORMATION-THEORETIC.  N_spur = n*eff^m.  On the
      lam*=0.07 curve at K1=8 the surplus is +15.3 bits (N_spur = 2.4e-5) and
      every variant scores 0/5; the information-theoretic wall is at K1=24
      (surplus -0.5 bits).  The observed wall sits ~19 bits inside the
      feasible region.

  E7  OVER-DETERMINATION DOES NOT HELP, AND THE MECHANISM IS NOW EXPLICIT.  At
      K1=8, m=32 the surplus is +74 bits and D/E still score 0/5.  sv/pv for
      variant D falls MONOTONICALLY with m (0.796, 0.704, 0.607, 0.515, 0.451
      for m = 8,12,16,24,32).  Reason: ||v_planted|| ~ n*sqrt(2m/3) grows with
      m, while the competing short vector -- the 2-D lambda-block vector
      mu ~ n/sqrt(eff) (T2's mu, in the (col i, col m+i) plane) -- is
      m-INDEPENDENT, because each new signature merely adds another copy of
      the same 2-D block.  So sv/pv ~ 1/sqrt(2*m*eff/3) decreases as
      1/sqrt(m): each extra signature adds log2(1/eff) bits of information but
      makes the lattice geometry strictly worse.  This is the mechanism behind
      T4b, and it is a property of the formulation, not of the instance.

  E8/E9  The leak is only log2(sqrt(n)/K1) ~ 2.7 bits per signature, the
      regime where Bleichenbacher's FFT method is classically the right tool.
      The GLV nonce bias factors exactly as
        |B| = |D(1,K1)| * |D(lam,K2)|,  |D(c,K)| = |sin(pi cK/n)|/(K|sin(pi c/n)|)
      (closed form agrees with brute force to 4 decimals: 0.0021 / 0.0186 /
      0.0805 on the three historical curves).  |B| is LARGEST on the curve
      where LLL fails, suggesting complementarity -- but E9 FALSIFIES |B| as a
      predictor: at fixed n, m, K1, K2 over 40 lam values, AUC(|B| predicts
      total failure) = 0.445 (A) / 0.417 (D), i.e. below chance and far below
      the 0.725 majority baseline; Spearman(|B|, wins) = +0.05.  Recorded as a
      DEAD hypothesis.  E9 also independently re-confirms the 2026-07-29
      falsification of lam*: AUC(lam*) = 0.357 at fixed eff.  29 of 40 lam
      values recover at eff=0.157, so the historical lam=185 curve is not
      special for having small lam* -- it is simply an unlucky lam.

Run: python3 glv_hnp_phase2_dfree.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ, GSO, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic -- same as glv_hnp_phase2_lambda_threshold.py
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

def find_point(p, b, seed=42):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            return (x, y)
    return None

def find_generator(p, b, n):
    rng = random.Random(12345)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    return r1, r2

def lam_star(lam, n):
    lr = lam % n
    return min(lr, n - lr) / n

# ---------------------------------------------------------------------------
# Signatures (identical generation to the 2026-07-29 run, so the A-column of
# every grid below is directly comparable to that run's T4 table)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- verbatim from the 2026-07-29 script."""
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
# Variant A -- baseline lattice WITH the d column (dim 2m+2)
# ---------------------------------------------------------------------------

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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
    M[2 * m + 1][dim - 1] = S_KAN
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

def recover_A(reduced, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Variant D -- d-free Kannan lattice (dim 2m+1)
# ---------------------------------------------------------------------------

def dfree_taus(sigs, n):
    """tau_i = B_i * B_0^{-1} mod n.  Returns (taus, B0inv) or None."""
    B0 = sigs[0]['B'] % n
    if B0 == 0 or math.gcd(B0, n) != 1:
        return None
    B0inv = modinv(B0, n)
    taus = [sigs[i]['B'] * B0inv % n for i in range(len(sigs))]
    assert taus[0] == 1
    return taus, B0inv

def build_lattice_D(sigs, n, lam, k1_bound, k2_bound, with_kannan=True):
    """d-free lattice.  Columns 0..m-1 = k1 (S_K1), m..2m-1 = k2 (S_K2),
    and (if with_kannan) column 2m = Kannan (S_KANNAN)."""
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ncol = 2 * m + (1 if with_kannan else 0)
    tt = dfree_taus(sigs, n)
    if tt is None:
        return None
    taus, _ = tt
    rows = []
    # k2 rows R_i
    for i in range(m):
        r = [0] * ncol
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    # d row Q_0 = S_K1 * (1, tau_1, ..., tau_{m-1})
    q0 = [0] * ncol
    for i in range(m):
        q0[i] = taus[i] * S_K1
    rows.append(q0)
    # mod-n rows Q_i, i = 1..m-1
    for i in range(1, m):
        r = [0] * ncol
        r[i] = n * S_K1
        rows.append(r)
    if with_kannan:
        k = [0] * ncol
        for i in range(m):
            k[i] = sigs[i]['A'] * S_K1
        k[2 * m] = S_KAN
        rows.append(k)
    return rows

def planted_D(sigs, n, k1_bound, k2_bound, with_kannan=True):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ncol = 2 * m + (1 if with_kannan else 0)
    v = [0] * ncol
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    if with_kannan:
        v[2 * m] = S_KAN
    return v

def dfree_membership_coeffs(sigs, d_secret, n, lam, k1_bound, k2_bound):
    """Explicit integer coefficients of the planted vector in the variant-D
    basis.  Returns (k2coeffs, a, c_list) or None if any c_i is non-integral
    (which would mean the basis derivation is wrong)."""
    m = len(sigs)
    taus, _ = dfree_taus(sigs, n)
    a = sigs[0]['k1'] - sigs[0]['A'] + lam * sigs[0]['k2']
    cs = []
    for i in range(1, m):
        num = sigs[i]['k1'] - sigs[i]['A'] + lam * sigs[i]['k2'] - a * taus[i]
        if num % n != 0:
            return None
        cs.append(num // n)
    return [s['k2'] for s in sigs], a, cs

def recover_D(reduced, m, n, lam, sigs, k1_bound, k2_bound, d_secret):
    """Scan reduced rows for |Kannan| = S_KANNAN, read (k1_i, k2_i) off the
    scaled columns, solve for d.  Requires exact divisibility by the scales."""
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    _, B0inv = dfree_taus(sigs, n)
    A0 = sigs[0]['A']
    for row in reduced:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        sg = 1 if last > 0 else -1
        if (sg * row[0]) % S_K1 != 0 or (sg * row[m]) % S_K2 != 0:
            continue
        k1_0 = (sg * row[0]) // S_K1
        k2_0 = (sg * row[m]) // S_K2
        d_cand = B0inv * (k1_0 + lam * k2_0 - A0) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Variant E -- explicit CVP on the d-free lattice (no Kannan row), dim 2m
# ---------------------------------------------------------------------------

def babai_nearest_plane(basis_rows, target):
    """Babai nearest-plane on an already-LLL-reduced integer basis."""
    dim = len(basis_rows)
    ncol = len(basis_rows[0])
    # Gram-Schmidt in floats (dims here are <= 24, entries <= ~n^2, fine)
    gs = []
    mu = [[0.0] * dim for _ in range(dim)]
    for i in range(dim):
        v = [float(x) for x in basis_rows[i]]
        for j in range(i):
            nj = sum(x * x for x in gs[j])
            if nj == 0:
                mu[i][j] = 0.0
                continue
            mu[i][j] = sum(basis_rows[i][k] * gs[j][k] for k in range(ncol)) / nj
            for k in range(ncol):
                v[k] -= mu[i][j] * gs[j][k]
        gs.append(v)
    w = [float(x) for x in target]
    coeffs = [0] * dim
    for i in range(dim - 1, -1, -1):
        nj = sum(x * x for x in gs[i])
        if nj == 0:
            continue
        c = sum(w[k] * gs[i][k] for k in range(ncol)) / nj
        ci = int(round(c))
        coeffs[i] = ci
        for k in range(ncol):
            w[k] -= ci * basis_rows[i][k]
    out = [0] * ncol
    for i in range(dim):
        if coeffs[i]:
            for k in range(ncol):
                out[k] += coeffs[i] * basis_rows[i][k]
    return out

def run_variant_E(sigs, n, lam, k1_bound, k2_bound, d_secret, use_exact_cvp=True):
    """CVP formulation: lattice L (d-free, no Kannan), target t = -(A_i*S_K1|0).
    The error u - t of the closest vector u is the planted (k1|k2) vector."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    rows = build_lattice_D(sigs, n, lam, k1_bound, k2_bound, with_kannan=False)
    if rows is None:
        return None, None
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(2 * m)] for i in range(2 * m)]
    u = None
    if use_exact_cvp:
        try:
            u = list(CVP.closest_vector(A, tuple(t)))
        except Exception:
            u = None
    if u is None:
        u = babai_nearest_plane(red, t)
    err = [u[k] - t[k] for k in range(2 * m)]
    _, B0inv = dfree_taus(sigs, n)
    A0 = sigs[0]['A']
    for cand in (err, [-x for x in err]):
        if cand[0] % S_K1 != 0 or cand[m] % S_K2 != 0:
            continue
        k1_0 = cand[0] // S_K1
        k2_0 = cand[m] // S_K2
        d_cand = B0inv * (k1_0 + lam * k2_0 - A0) % n
        if d_cand == d_secret:
            return d_cand, norm(err)
    return None, norm(err)

# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def run_A(sigs, n, lam, k1_bound, k2_bound, d_secret, use_bkz=False, beta=20):
    m = len(sigs)
    dim = 2 * m + 2
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_lattice_A(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz: BKZ.reduction(A, BKZ.Param(beta))
    else:       LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_A(red, m, n, S_KAN, d_secret) is not None
    pv = norm(planted_A(sigs, d_secret, n, k1_bound, k2_bound))
    sv = min(norm(r) for r in red)
    svrow = min(red, key=norm)
    return ok, pv, sv, svrow

def run_D(sigs, n, lam, k1_bound, k2_bound, d_secret, use_bkz=False, beta=20):
    m = len(sigs)
    dim = 2 * m + 1
    rows = build_lattice_D(sigs, n, lam, k1_bound, k2_bound, with_kannan=True)
    if rows is None:
        return None
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz: BKZ.reduction(A, BKZ.Param(beta))
    else:       LLL.reduction(A)
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_D(red, m, n, lam, sigs, k1_bound, k2_bound, d_secret) is not None
    pv = norm(planted_D(sigs, n, k1_bound, k2_bound, with_kannan=True))
    sv = min(norm(r) for r in red)
    svrow = min(red, key=norm)
    return ok, pv, sv, svrow

def instance(curve, m, k1_bound, seed):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    return sigs, d, k2_bound


# ===========================================================================
print("=" * 78)
print("Thread 23 -- d-free reformulation of the GLV-HNP Phase-2 lattice")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist.append((label, (p, b, n, lam, G), k1, m))


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E1: correctness of the d-free basis and of d-recovery")
print("-" * 78)
print("(a) planted vector has INTEGRAL coordinates in the variant-D basis")
print("(b) reconstructing from those coordinates reproduces it exactly")
print("(c) det(L_D) == det(L_A)/n  -- the quotient is exactly index n\n")

print(f"{'curve':<18} {'m':>3} {'K1':>3} {'(a) integral':>13} {'(b) exact':>10} "
      f"{'(c) det ratio':>14}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    inst = instance(curve, m, k1, 42)
    sigs, d, k2b = inst
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    co = dfree_membership_coeffs(sigs, d, n, lam, k1, k2b)
    integral = co is not None
    exact = False
    if integral:
        k2c, a, cs = co
        rows = build_lattice_D(sigs, n, lam, k1, k2b, with_kannan=True)
        # rows: R_0..R_{m-1}, Q_0, Q_1..Q_{m-1}, K
        acc = [0] * (2 * m + 1)
        def addrow(r, coef):
            for j in range(2 * m + 1):
                acc[j] += coef * r[j]
        for i in range(m):
            addrow(rows[i], k2c[i])
        addrow(rows[m], a)
        for i in range(1, m):
            addrow(rows[m + i], cs[i - 1])
        addrow(rows[2 * m], 1)
        exact = (acc == planted_D(sigs, n, k1, k2b, with_kannan=True))
    # determinants (Gram determinant via sympy for exactness)
    MA = build_lattice_A(sigs, n, lam, k1, k2b)
    detA = abs(sympy.Matrix(MA).det())
    rows_nk = build_lattice_D(sigs, n, lam, k1, k2b, with_kannan=False)
    detD = abs(sympy.Matrix(rows_nk).det())
    # det(L_A) includes the Kannan column (factor S_KAN) and the d column
    # (factor S_D=1); strip the Kannan factor from both sides for comparison
    detA_nokan = detA // S_KAN
    ratio = sympy.Rational(detA_nokan, detD)
    print(f"{label:<18} {m:>3} {k1:>3} {str(integral):>13} {str(exact):>10} "
          f"{str(ratio):>14}  (n={n})")

print("\n(d) end-to-end: does variant D recover d where variant A does?")
print(f"{'curve':<18} {'K1':>3} {'A':>6} {'D':>6} {'E(CVP)':>8}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    wa = wd = we = 0
    for s in SEEDS:
        sigs, d, k2b = instance(curve, m, k1, s)
        ra = run_A(sigs, n, lam, k1, k2b, d)
        rd = run_D(sigs, n, lam, k1, k2b, d)
        de, _ = run_variant_E(sigs, n, lam, k1, k2b, d)
        wa += 1 if ra[0] else 0
        wd += 1 if (rd and rd[0]) else 0
        we += 1 if de is not None else 0
    print(f"{label:<18} {k1:>3} {str(wa)+'/5':>6} {str(wd)+'/5':>6} {str(we)+'/5':>8}")


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E2: sv/pv and shortest-vector energy profile, A vs D")
print("-" * 78)
print("Variant A's sv is always n*S_D*e_m (T5).  Does dropping the d column")
print("raise sv/pv, and where does the new sv put its energy?\n")

def energy(row, m, has_d):
    """Fraction of squared norm in each block."""
    tot = sum(x * x for x in row)
    if tot == 0: return (0, 0, 0, 0)
    if has_d:
        k1 = sum(row[i] ** 2 for i in range(m))
        dd = row[m] ** 2
        k2 = sum(row[m + 1 + i] ** 2 for i in range(m))
        kan = row[2 * m + 1] ** 2
    else:
        k1 = sum(row[i] ** 2 for i in range(m))
        dd = 0
        k2 = sum(row[m + i] ** 2 for i in range(m))
        kan = row[2 * m] ** 2
    return (k1 / tot, dd / tot, k2 / tot, kan / tot)

print(f"{'curve':<18} {'K1':>3} | {'A sv/pv':>8} {'A:k1':>6} {'A:d':>6} {'A:k2':>6} "
      f"{'A:kan':>6} | {'D sv/pv':>8} {'D:k1':>6} {'D:k2':>6} {'D:kan':>6}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    sigs, d, k2b = instance(curve, m, k1, 42)
    oka, pva, sva, rowa = run_A(sigs, n, lam, k1, k2b, d)
    rd = run_D(sigs, n, lam, k1, k2b, d)
    okd, pvd, svd, rowd = rd
    ea = energy(rowa, m, True)
    ed = energy(rowd, m, False)
    print(f"{label:<18} {k1:>3} | {sva/pva:>8.3f} {ea[0]:>6.3f} {ea[1]:>6.3f} "
          f"{ea[2]:>6.3f} {ea[3]:>6.3f} | {svd/pvd:>8.3f} {ed[0]:>6.3f} "
          f"{ed[2]:>6.3f} {ed[3]:>6.3f}")

print("\nAlso: |sv[d-col]|/n for variant A (T5 reported exactly 1.0000):")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    sigs, d, k2b = instance(curve, m, k1, 42)
    _, _, _, rowa = run_A(sigs, n, lam, k1, k2b, d)
    print(f"  {label:<18} |sv[m]|/n = {abs(rowa[m])/n:.4f}")


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E3: THE FALSIFIER -- K1 wall grid, A vs D vs E")
print("-" * 78)
print("2026-07-29 T4 baseline (variant A), for reference:")
print("  12-bit/2557 (lam*=0.340): 5,5,5,5,5,4,1,0  for K1=2,3,4,6,8,12,16,24")
print("  12-bit/2677 (lam*=0.070): 5,5,5,2,0,0,0,0  for K1=2,3,4,6,8,12,16,24")
print("If D/E move the 2677 wall outward from K1~4-6, the reformulation helps.")
print("If not, the wall is information-theoretic and Phase 2 is at its ceiling.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
GRID_CURVES = [(label, curve, m) for label, curve, k1, m in hist if label.startswith("12-bit")]

for label, curve, m in GRID_CURVES:
    p, b, n, lam, G = curve
    print(f"{label}  (n={n}, lam={lam}, lam*={lam_star(lam,n):.4f}, m={m})")
    hdr = "  " + f"{'variant':<10}" + "".join(f"{('K1='+str(k)):>8}" for k in K1_GRID)
    print(hdr)
    for vname in ("A", "D", "E(CVP)"):
        cells = []
        for k1 in K1_GRID:
            w = 0
            for s in SEEDS:
                inst = instance(curve, m, k1, s)
                if inst is None:
                    continue
                sigs, d, k2b = inst
                if vname == "A":
                    w += 1 if run_A(sigs, n, lam, k1, k2b, d)[0] else 0
                elif vname == "D":
                    r = run_D(sigs, n, lam, k1, k2b, d)
                    w += 1 if (r and r[0]) else 0
                else:
                    de, _ = run_variant_E(sigs, n, lam, k1, k2b, d)
                    w += 1 if de is not None else 0
            cells.append(f"{w}/5")
        print("  " + f"{vname:<10}" + "".join(f"{c:>8}" for c in cells))
    # mean sv/pv per K1 for A and D
    for vname in ("A", "D"):
        cells = []
        for k1 in K1_GRID:
            rs = []
            for s in SEEDS:
                inst = instance(curve, m, k1, s)
                if inst is None: continue
                sigs, d, k2b = inst
                if vname == "A":
                    _, pv, sv, _ = run_A(sigs, n, lam, k1, k2b, d)
                else:
                    r = run_D(sigs, n, lam, k1, k2b, d)
                    if r is None: continue
                    _, pv, sv, _ = r
                rs.append(sv / pv)
            cells.append(f"{sum(rs)/len(rs):.3f}" if rs else "-")
        print("  " + f"{'sv/pv '+vname:<10}" + "".join(f"{c:>8}" for c in cells))
    print()


# ===========================================================================
print("-" * 78)
print("EXP E4: BKZ(20) on variant D at the wall")
print("-" * 78)
print(f"{'curve':<18} {'K1':>3} {'D-LLL':>7} {'D-BKZ20':>9} {'A-BKZ20':>9}")
for label, curve, m in GRID_CURVES:
    p, b, n, lam, G = curve
    for k1 in (6, 8, 12):
        wl = wb = wab = 0
        for s in SEEDS:
            inst = instance(curve, m, k1, s)
            if inst is None: continue
            sigs, d, k2b = inst
            r = run_D(sigs, n, lam, k1, k2b, d)
            wl += 1 if (r and r[0]) else 0
            rb = run_D(sigs, n, lam, k1, k2b, d, use_bkz=True, beta=20)
            wb += 1 if (rb and rb[0]) else 0
            wab += 1 if run_A(sigs, n, lam, k1, k2b, d, use_bkz=True, beta=20)[0] else 0
        print(f"{label:<18} {k1:>3} {str(wl)+'/5':>7} {str(wb)+'/5':>9} {str(wab)+'/5':>9}")


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E5: control -- is 'planted vector is not lambda_1' pathological?")
print("-" * 78)
print("Build the CLASSICAL Nguyen-Shparlinski HNP lattice (no GLV, k_i < n/2^l,")
print("rows n*e_i and (B_1..B_m, 1/2^{l+1}) scaled to integers) and measure")
print("sv/pv.  If sv/pv < 1 there too, T5's finding is generic to HNP, not a")
print("defect of the Phase-2 construction.\n")

def ns_hnp_experiment(n, m, ell, seed):
    """Classical HNP: A_i + B_i*d = k_i mod n, 0 <= k_i < n/2^ell."""
    rng = random.Random(seed)
    d = rng.randint(1, n - 1)
    kb = n >> ell
    sigs = []
    for _ in range(m):
        B = rng.randint(1, n - 1)
        k = rng.randint(0, kb - 1)
        A = (k - B * d) % n
        sigs.append((A, B, k))
    # integer-scaled NS lattice: multiply the classical basis by 2^{ell+1}
    S = 1 << (ell + 1)
    dim = m + 2
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S
    for i in range(m):
        M[m][i] = sigs[i][1] * S
    M[m][m] = 1
    for i in range(m):
        M[m + 1][i] = sigs[i][0] * S
    M[m + 1][dim - 1] = n
    A_ = IntegerMatrix.from_matrix(M)
    LLL.reduction(A_)
    red = [[A_[i][j] for j in range(dim)] for i in range(dim)]
    pv = [0] * dim
    for i in range(m):
        pv[i] = sigs[i][2] * S
    pv[m] = d
    pv[dim - 1] = n
    sv = min(norm(r) for r in red)
    svrow = min(red, key=norm)
    # recovery
    ok = False
    for row in red:
        if abs(row[dim - 1]) != n: continue
        sg = 1 if row[dim - 1] > 0 else -1
        if (sg * row[m]) % n == d % n:
            ok = True
    return ok, norm(pv), sv, abs(svrow[m]) / 1.0, svrow

print(f"{'n bits':>7} {'m':>3} {'ell':>4} {'recovered':>10} {'sv/pv':>8} "
      f"{'|sv[d]|':>10} {'sv==n*e_d':>10}")
for bits, m, ell in ((64, 20, 8), (64, 40, 4), (128, 40, 8), (128, 60, 4)):
    n = int(sympy.nextprime(2 ** bits))
    ok, pv, sv, svd, svrow = ns_hnp_experiment(n, m, ell, 4242)
    trivial = (abs(svrow[m]) == n and all(x == 0 for j, x in enumerate(svrow) if j != m))
    print(f"{bits:>7} {m:>3} {ell:>4} {str(ok):>10} {sv/pv:>8.3f} {svd:>10.0f} "
          f"{str(trivial):>10}")


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E6: is the K1 wall information-theoretic?")
print("-" * 78)
print("Counting argument.  Candidates for (d, k2_1..k2_m) number n * K2^m; each")
print("forces k1_i = A_i + B_i d - lam k2_i mod n to land in [0,K1), which for a")
print("random candidate happens with probability (K1/n)^m.  Expected number of")
print("spurious solutions is therefore")
print("    N_spur = n * (K1*K2/n)^m = n * eff^m,")
print("so the solution is information-theoretically unique iff")
print("    m * log2(1/eff) > log2(n)      [surplus = LHS - RHS, in bits].")
print("If the observed wall sits at a K1 where the surplus is still large and")
print("positive, the wall is NOT information-theoretic.\n")

print(f"{'curve':<18} {'m':>3} {'K1':>3} {'eff':>7} {'bits needed':>12} "
      f"{'bits avail':>11} {'surplus':>9} {'N_spur':>11} {'A':>5} {'D':>5} {'E':>5}")
for label, curve, m in GRID_CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for k1 in K1_GRID:
        eff = k1 * k2b / n
        need = math.log2(n)
        avail = m * math.log2(1.0 / eff) if eff < 1 else 0.0
        nspur = n * (eff ** m)
        wa = wd = we = 0
        for s in SEEDS:
            inst = instance(curve, m, k1, s)
            if inst is None: continue
            sigs, d, _kb = inst
            wa += 1 if run_A(sigs, n, lam, k1, _kb, d)[0] else 0
            r = run_D(sigs, n, lam, k1, _kb, d)
            wd += 1 if (r and r[0]) else 0
            de, _ = run_variant_E(sigs, n, lam, k1, _kb, d)
            we += 1 if de is not None else 0
        print(f"{label:<18} {m:>3} {k1:>3} {eff:>7.3f} {need:>12.2f} "
              f"{avail:>11.2f} {avail-need:>+9.2f} {nspur:>11.2e} "
              f"{str(wa)+'/5':>5} {str(wd)+'/5':>5} {str(we)+'/5':>5}")
    print()


# ===========================================================================
print("-" * 78)
print("EXP E7: does over-determination rescue the d-free variants at the wall?")
print("-" * 78)
print("T4b (2026-07-29) found variant A at K1=8 is not rescued by more data")
print("(m = 8/12/16/24/32 -> 0,0,1,0,1 of 5).  E6 shows the surplus at K1=8 is")
print("already huge, so if D or E scale with m the loss is in the lattice, not")
print("in the information.  Exact CVP is capped at dim 40; Babai beyond.\n")

label, curve, _m0 = GRID_CURVES[1]   # the lam*=0.07 failure curve
p, b, n, lam, G = curve
k2b_fix = math.isqrt(n) + 1
print(f"{label} (n={n}, lam*={lam_star(lam,n):.4f}), K1=8, eff={8*k2b_fix/n:.3f}")
print(f"{'m':>4} {'surplus(bits)':>14} {'A':>6} {'D':>6} {'E':>6} {'sv/pv D':>9}")
for m in (8, 12, 16, 24, 32):
    eff = 8 * k2b_fix / n
    surplus = m * math.log2(1.0 / eff) - math.log2(n)
    wa = wd = we = 0
    ratios = []
    for s in SEEDS:
        inst = instance(curve, m, 8, s)
        if inst is None: continue
        sigs, d, _kb = inst
        wa += 1 if run_A(sigs, n, lam, 8, _kb, d)[0] else 0
        r = run_D(sigs, n, lam, 8, _kb, d)
        if r:
            wd += 1 if r[0] else 0
            ratios.append(r[2] / r[1])
        de, _ = run_variant_E(sigs, n, lam, 8, _kb, d,
                              use_exact_cvp=(2 * m <= 40))
        we += 1 if de is not None else 0
    rr = f"{sum(ratios)/len(ratios):.3f}" if ratios else "-"
    print(f"{m:>4} {surplus:>14.1f} {str(wa)+'/5':>6} {str(wd)+'/5':>6} "
          f"{str(we)+'/5':>6} {rr:>9}")


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E8: Bleichenbacher bias of the GLV nonce -- the complementary method")
print("-" * 78)
print("E6/E7 show the leak is only log2(1/eff) = log2(sqrt(n)/K1) ~ 2.7 bits per")
print("signature and that the lattice loses ground as 1/sqrt(m).  The classical")
print("remedy in the few-bits regime is Bleichenbacher's FFT/bias method, which")
print("needs no lattice at all.  Its figure of merit is the nonce bias")
print("    B(K) = E[exp(2*pi*i*k/n)],  k = k1 + lam*k2 mod n.")
print("Because k1 and k2 are independent this factors exactly:")
print("    B = D(1, K1) * D(lam, K2),   D(c,K) = (1/K) * sum_{j<K} e^{2pi i c j/n}")
print("             |D(c,K)| = |sin(pi*c*K/n)| / (K*|sin(pi*c/n)|).")
print("The ENVELOPE of |D(lam,K2)| is 1/(K2*sin(pi*lam*)), which grows as")
print("lam* -> 0 -- exactly the regime where the lattice attack dies.  That")
print("suggests an anti-correlation; E9 tests it and refutes it.\n")

import cmath

def bias_factor(c, K, n):
    """|D(c,K)| = |(1/K) sum_{j<K} e^{2 pi i c j / n}|, exact closed form."""
    c %= n
    if c == 0:
        return 1.0
    num = abs(math.sin(math.pi * c * K / n))
    den = K * abs(math.sin(math.pi * c / n))
    if den == 0:
        return 1.0
    return num / den

def bias_direct(lam, K1, K2, n):
    """Brute-force check of B = E[e^{2 pi i k/n}] over all (k1,k2) pairs."""
    acc = 0j
    for a in range(K1):
        for bb in range(K2):
            acc += cmath.exp(2j * math.pi * ((a + lam * bb) % n) / n)
    return abs(acc) / (K1 * K2)

print(f"{'curve':<18} {'lam*':>7} {'K1':>3} {'K2':>4} {'|D(1,K1)|':>10} "
      f"{'|D(lam,K2)|':>12} {'|B| closed':>11} {'|B| direct':>11} "
      f"{'M~1/|B|^2':>11} {'LLL A':>6}")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    f1 = bias_factor(1, k1, n)
    f2 = bias_factor(lam, k2b, n)
    Bc = f1 * f2
    Bd = bias_direct(lam, k1, k2b, n)
    wa = 0
    for s in SEEDS:
        sigs, d, _kb = instance(curve, m, k1, s)
        wa += 1 if run_A(sigs, n, lam, k1, _kb, d)[0] else 0
    print(f"{label:<18} {lam_star(lam,n):>7.4f} {k1:>3} {k2b:>4} {f1:>10.4f} "
          f"{f2:>12.4f} {Bc:>11.4f} {Bd:>11.4f} {1.0/Bc**2:>11.1f} "
          f"{str(wa)+'/5':>6}")

print("\nBias vs lam* at fixed K1=8, K2=sqrt(n), for the two 12-bit moduli:")
print(f"{'n':>6} {'lam*':>7} {'|B|':>9} {'M~1/|B|^2':>11}")
for label, curve, m in GRID_CURVES:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    for ls in (0.005, 0.01, 0.02, 0.05, 0.07, 0.15, 0.25, 0.34, 0.45):
        lam_f = max(1, int(round(ls * n)))
        Bc = bias_factor(1, 8, n) * bias_factor(lam_f, k2b, n)
        print(f"{n:>6} {ls:>7.3f} {Bc:>9.4f} {1.0/Bc**2:>11.1f}")
    break

print("\nCAUTION: |B| is NOT monotone in lam*.  Only its ENVELOPE")
print("1/(K2*sin(pi*lam*)) decreases; the numerator |sin(pi*lam**K2)| oscillates")
print("with period ~1/K2 in lam*, so |B| has near-zeros (see lam*=0.25 above,")
print("|B|=0.0004).  The three historical curves happen to line up with the")
print("anti-correlation (|B| = 0.0021 / 0.0186 / 0.0805 against LLL 5/5, 5/5,")
print("0/5), but that is 3 points against an oscillating function -- exactly the")
print("shape of the lam/n error made on 2026-07-26.  E9 tests it properly.")


# ===========================================================================
print("\n" + "-" * 78)
print("EXP E9: falsifier for E8 -- does |B| predict LLL failure at fixed eff?")
print("-" * 78)
print("Everything held fixed (n=2647, m=10, K1=8, K2=52, so eff and the")
print("information surplus are constant); only lam varies, over 40 values.")
print("For the HNP/lattice structure lam may be ANY residue -- the GLV")
print("eigenvalue condition lam^2+lam+1=0 matters only for realisability of")
print("such nonces by an implementation, not for the lattice geometry.")
print("If |B| separates LLL successes from failures here, E8's anti-correlation")
print("is real; if not, it is another 3-point coincidence.\n")

label9, curve9, m9 = GRID_CURVES[1]
p9, b9, n9, lam9, G9 = curve9
k2b9 = math.isqrt(n9) + 1
K1_9, M_9 = 8, 10

rows9 = []
for idx in range(40):
    lam_t = 1 + (idx * n9) // 80          # spread over (0, n/2)
    Bc = bias_factor(1, K1_9, n9) * bias_factor(lam_t, k2b9, n9)
    wa = wd = 0
    for s in SEEDS:
        d = random.Random(s + 7777).randint(1, n9 - 1)
        sigs = gen_signatures(G9, d, M_9, n9, lam_t, p9, K1_9, k2b9, s)
        if len(sigs) < M_9:
            continue
        wa += 1 if run_A(sigs, n9, lam_t, K1_9, k2b9, d)[0] else 0
        r = run_D(sigs, n9, lam_t, K1_9, k2b9, d)
        wd += 1 if (r and r[0]) else 0
    rows9.append((lam_t, lam_star(lam_t, n9), Bc, wa, wd))

print(f"{'lam':>6} {'lam*':>7} {'|B|':>9} {'A':>5} {'D':>5}   "
      f"{'lam':>6} {'lam*':>7} {'|B|':>9} {'A':>5} {'D':>5}")
half = (len(rows9) + 1) // 2
for i in range(half):
    l = rows9[i]
    left = f"{l[0]:>6} {l[1]:>7.4f} {l[2]:>9.5f} {str(l[3])+'/5':>5} {str(l[4])+'/5':>5}"
    if i + half < len(rows9):
        r_ = rows9[i + half]
        right = (f"{r_[0]:>6} {r_[1]:>7.4f} {r_[2]:>9.5f} "
                 f"{str(r_[3])+'/5':>5} {str(r_[4])+'/5':>5}")
    else:
        right = ""
    print(left + "   " + right)

def auc(scores, labels):
    """AUC of `scores` predicting label==1 (ties counted at 0.5)."""
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return float('nan')
    tot = 0.0
    for a in pos:
        for bb in neg:
            tot += 1.0 if a > bb else (0.5 if a == bb else 0.0)
    return tot / (len(pos) * len(neg))

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx = sum(rx) / len(rx); my = sum(ry) / len(ry)
    num = sum((a - mx) * (bb - my) for a, bb in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((bb - my) ** 2 for bb in ry))
    return num / (dx * dy) if dx and dy else float('nan')

for vi, vname in ((3, "A"), (4, "D")):
    Bs = [r[2] for r in rows9]
    ls_ = [r[1] for r in rows9]
    fails = [1 if r[vi] == 0 else 0 for r in rows9]
    wins = [r[vi] for r in rows9]
    nfail = sum(fails)
    base = max(nfail, len(rows9) - nfail) / len(rows9)
    print(f"\nvariant {vname}: {len(rows9)-nfail}/{len(rows9)} lam values with "
          f">=1/5 recovery; majority baseline {base:.3f}")
    print(f"  AUC(|B|  predicts total failure) = {auc(Bs, fails):.3f}")
    print(f"  AUC(lam* predicts total failure) = {auc([-x for x in ls_], fails):.3f}"
          "   (sign: small lam* -> failure)")
    print(f"  Spearman(|B|,  wins) = {spearman(Bs, wins):+.3f}")
    print(f"  Spearman(lam*, wins) = {spearman(ls_, wins):+.3f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
