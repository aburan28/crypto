"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
lambda_1 (or at least so the trivial vector is gone).

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5):
  The Phase-2 lattice built by `build_glv_lattice()` in
  glv_hnp_phase2_lambda_threshold.py:254 (dim 2m+2) always contains the
  trivial vector  n*S_D*e_m  (the "d column" direction).  Its norm is
  n*S_D, while ||v_planted||^2 ~ n^2*(2m/3 + 4/3), so the trivial vector is
  2-3x SHORTER than the planted vector on every curve, success and failure
  alike (sv/pv in [0.34, 0.61]).  The planted vector is therefore never
  lambda_1 and recovery is a BDD/coset condition, not an SVP condition.

  The trivial vector is pure noise: d is only defined mod n, and n*e_m
  carries no information.  It cannot be scaled away (both it and the planted
  vector are linear in S_D).

This script implements the two reformulations proposed in that entry.

  REFORMULATION A ("projected", dim 2m+1).
    Quotient out the e_m direction: drop the d column entirely.  The kernel
    of the projection, restricted to the lattice, is exactly Z*(n*S_D*e_m),
    so the image has rank 2m+1 and determinant det(L)/(n*S_D).  Concretely
    the k1-block becomes S_K1 * L_0 where
        L_0 = { (t*B_1 + c_1*n, ..., t*B_m + c_m*n) : t, c_i in Z }
    is the classical Boneh-Venkatesan HNP lattice of rank m and determinant
    n^(m-1).  Basis (using B_1 invertible mod n, n prime):
        r_0 = S_K1 * (1, b_2, ..., b_m)      with b_i = B_1^{-1} B_i mod n
        r_i = n * S_K1 * e_i                 i = 1 .. m-1
    plus the m unchanged k2 rows and the unchanged Kannan row.
    d is no longer read off a coordinate; it is reconstructed from the
    recovered nonce split:  d = B_1^{-1} (k1_1 + lam*k2_1 - A_1)  mod n.

  REFORMULATION B ("Babai CVP", dim 2m, no Kannan embedding).
    Same k1/k2 blocks, no Kannan row.  Centered target
        T = ((K1/2 - A_i)*S_K1)_i  ++  ((K2/2)*S_K2)_i
    and solve BDD with Babai nearest-plane on the LLL-reduced basis.
    This removes the embedding entirely, so there is no coset condition at
    all -- only a covering-radius condition.

  REFORMULATION C ("projc" = A + centering).  Added after the first run of
    U3 showed B beating A: reformulation A is a no-op for recovery while B
    moves the K1 wall, and B differs from A in TWO ways (no embedding, and
    a centered target).  C isolates the second: it is exactly A with the
    Kannan row embedding the CENTERED target
        (A_i - K1/2)*S_K1  in the k1 block,  -(K2/2)*S_K2  in the k2 block.
    The planted vector's nonce coordinates become (k1_i - K1/2) and
    (k2_i - K2/2), so
        ||v_planted||^2 : n^2*(2m/3 + 1)  ->  n^2*(m/6 + 1),
    a ~2x shorter BDD target for free, using only the public bounds K1, K2.

Falsifier (stated verbatim in the 2026-07-29 entry):
    "if sv/pv rises above 1 after the reformulation and the K1 wall in T4
     moves outward on the lam*=0.07 curve (currently K1~4-6), the
     reformulation is a real improvement; if the wall stays at K1~4-6, then
     the wall is information-theoretic and Phase 2 is at its ceiling."

Experiments:
  U1  sanity: planted vector is in the projected lattice; recovery works.
  U2  sv/pv in the projected lattice vs the original (the T5 table, redone).
  U3  the T4 K1 grid, original vs projected vs Babai, on both 12-bit curves.
  U4  T4b control: does more data (m sweep) rescue the lam*=0.07 curve at
      the K1 where the original wall sits?
  U5  17-bit sweep at eff = 0.15 and 0.25 (where the original scored 3/20
      and 0/20) -- does the projection move the population-level wall?
  U6  2x2 factorial {projection} x {centering} -- attributes the effect.
  U7  where the centered variant's population wall sits (eff sweep).
  (U3/U4/U5 all compare orig / proj / projc / cvp side by side.)

RESULTS (2026-08-04 run; full output in glv_hnp_phase2_projected_output.txt):

  U1  Sanity passes on all three historical curves: the planted vector is an
      integer combination of the projected basis, and det = S_K1^m * n^(m-1)
      * S_K2^m * S_KANNAN exactly (i.e. det(L)/(n*S_D), as predicted).

  U2  Projection does raise sv/pv as intended (the trivial vector is gone):
      0.543 -> 0.790, 0.445 -> 0.475, 0.387 -> 0.767.  Still < 1: the planted
      vector is STILL not lambda_1.  Centering pushes it to 0.906/0.629/0.998.

  U3/U6  REFORMULATION A IS A NO-OP.  Removing the trivial vector changes the
      recovery rate in NOT ONE of the 16 (curve, K1) cells, nor in any of the
      10 (curve, m) cells of U4, nor on any of the 40 curve-instances of U5.
      orig and proj agree cell-for-cell everywhere.  The trivial vector
      n*S_D*e_m identified in T5 was a red herring: it is short, but it never
      obstructed recovery.

  U6  CENTERING IS THE WHOLE EFFECT.  The 2x2 factorial {projection} x
      {centering} shows centering carries the gain with or without the
      projection (origc ~ projc), and projection carries none with or
      without centering (orig == proj).

  U3  The K1 wall MOVES OUTWARD under centering.  On the lam*=0.07 curve
      (p=2677, n=2647, m=10) the wall goes from K1 ~ 4-6 to K1 ~ 12-16:
        orig/proj   5/5 5/5 5/5 0/5 0/5 0/5 0/5 0/5   (K1 = 2,3,4,6,8,12,16,24)
        origc/projc 5/5 5/5 5/5 5/5 5/5 4/5 1/5 1/5
      On the lam*=0.34 curve the wall goes from K1 ~ 12 to K1 ~ 24.

  U4  This also overturns the T4b control: at K1=8 the lam*=0.07 curve was
      NOT rescued by more data (0,0,1,0,1 of 5 for m=8,12,16,24,32).  With
      centering it is 4/5 at m=8 and 5/5 for every m >= 12.

  U5  Population effect at 17 bits, m=12, 20 curves: full recovery on
        eff=0.15   orig 3/20  ->  projc 14/20
        eff=0.25   orig 0/20  ->  projc  4/20

  U7  The new population wall sits at eff ~ 0.35 (mean win rate 0.22) and is
      hard by eff ~ 0.45 (0.02).  The old wall was at eff ~ 0.15.  Centering
      therefore roughly DOUBLES the tolerable bias strength eff = K1*K2/n,
      from ~0.10 to ~0.25.

  Why: ||v_planted||^2 falls from n^2*(2m/3 + 1) to n^2*(m/6 + 1), i.e. a
  factor ~1.6-1.7x shorter BDD target at m=10, obtained for free from the
  public GLV bounds K1, K2.  Nothing about lam is involved.

  CONSEQUENCE for the 2026-07-29 falsifier: its two clauses split.  sv/pv did
  NOT rise above 1, yet the K1 wall DID move outward -- so "Phase 2 is at its
  ceiling" is refuted, but not by the mechanism T5 proposed.  Every Phase-2
  number logged between 2026-06-15 and 2026-07-29 was measured on an
  UNCENTERED target and understates the attack by ~2x in eff.

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, GSO

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + curve helpers
# (verbatim from glv_hnp_phase2_lambda_threshold.py so results are comparable)
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
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (Eisenstein norm form)."""
    if p % 3 != 1:
        return None
    lim = int(2 * math.isqrt(p)) + 2
    for a in range(0, lim + 1):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        r = math.isqrt(disc)
        if r * r != disc:
            continue
        if (a + r) % 2 == 0:
            b = (a + r) // 2
            if a * a - a * b + b * b == p:
                return (a, b)
        if (a - r) % 2 == 0:
            b = (a - r) // 2
            if a * a - a * b + b * b == p:
                return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None
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

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1:
                        continue
                    if not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

# ---------------------------------------------------------------------------
# Signatures and scales (verbatim)
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
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# ORIGINAL lattice (dim 2m+2), verbatim from Thread 20
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, center=False):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    c1 = k1_bound // 2 if center else 0
    c2 = k2_bound // 2 if center else 0
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
        M[2 * m + 1][i] = (sigs[i]['A'] - c1) * S_K1
        M[2 * m + 1][m + 1 + i] = -c2 * S_K2
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector_orig(sigs, d_secret, n, k1_bound, k2_bound, center=False):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1 = k1_bound // 2 if center else 0
    c2 = k2_bound // 2 if center else 0
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - c2) * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_d_orig(reduced, m, n, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# REFORMULATION A -- projected lattice (dim 2m+1), d column removed
# ---------------------------------------------------------------------------
# Columns:  0 .. m-1     k1 block (scale S_K1)
#           m .. 2m-1    k2 block (scale S_K2)
#           2m           Kannan   (scale S_KANNAN)

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound, center=False):
    """Basis of pi(L), where pi drops the d column.  Rank 2m+1.

    center=True embeds the *centered* CVP target instead of the raw one:
    the Kannan row becomes (A_i - K1/2)*S_K1 in the k1 block and
    -(K2/2)*S_K2 in the k2 block, so the planted vector's nonce coordinates
    become (k1_i - K1/2) and (k2_i - K2/2), i.e. symmetric about 0.
    This is a legitimate reformulation: K1, K2 are public GLV bounds.
    ||v_planted||^2 drops from ~n^2(2m/3 + 1) to ~n^2(m/6 + 1).
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B = [s['B'] % n for s in sigs]
    # pick a pivot j with B[j] invertible mod n (n prime -> any B[j] != 0)
    j = next((i for i in range(m) if B[i] % n != 0), None)
    if j is None:
        return None, None
    Binv = modinv(B[j], n)
    bb = [(Binv * B[i]) % n for i in range(m)]   # bb[j] == 1

    rows = []
    # rank-m basis of S_K1 * L_0 :  the "t" row plus n*e_i for i != j
    r = [0] * dim
    for i in range(m):
        r[i] = bb[i] * S_K1
    rows.append(r)
    for i in range(m):
        if i == j:
            continue
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    # k2 rows
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        rows.append(r)
    # Kannan row (raw or centered target)
    c1 = k1_bound // 2 if center else 0
    c2 = k2_bound // 2 if center else 0
    r = [0] * dim
    for i in range(m):
        r[i] = (sigs[i]['A'] - c1) * S_K1
        r[m + i] = -c2 * S_K2
    r[2 * m] = S_KANNAN
    rows.append(r)
    assert len(rows) == dim
    return rows, j

def planted_vector_proj(sigs, n, k1_bound, k2_bound, center=False):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1 = k1_bound // 2 if center else 0
    c2 = k2_bound // 2 if center else 0
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + i] = (sigs[i]['k2'] - c2) * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_proj(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret,
                   center=False):
    """Scan for a row with |last| = S_KANNAN, read off (k1_i, k2_i),
    reconstruct d = B_i^{-1}(k1_i + lam*k2_i - A_i) mod n and verify it
    against ALL signatures (self-check, no oracle needed) before comparing
    to d_secret."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1 = k1_bound // 2 if center else 0
    c2 = k2_bound // 2 if center else 0
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        # candidate nonce split
        ok = True
        k1s, k2s = [], []
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 != 0 or b % S_K2 != 0:
                ok = False; break
            k1s.append(a // S_K1 + c1); k2s.append(b // S_K2 + c2)
        if not ok:
            continue
        for i in range(m):
            if sigs[i]['B'] % n == 0:
                continue
            d_cand = (modinv(sigs[i]['B'], n)
                      * (k1s[i] + lam * k2s[i] - sigs[i]['A'])) % n
            if d_cand == 0:
                continue
            # self-consistency: same d must explain every signature
            good = all((sigs[t]['A'] + sigs[t]['B'] * d_cand
                        - k1s[t] - lam * k2s[t]) % n == 0 for t in range(m))
            if good and d_cand == d_secret:
                return d_cand
            break
    return None

# ---------------------------------------------------------------------------
# REFORMULATION B -- Babai nearest-plane CVP (dim 2m, no embedding)
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Same as the projected lattice minus the Kannan row/column.  Rank 2m."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    B = [s['B'] % n for s in sigs]
    j = next((i for i in range(m) if B[i] % n != 0), None)
    if j is None:
        return None
    Binv = modinv(B[j], n)
    bb = [(Binv * B[i]) % n for i in range(m)]
    rows = []
    r = [0] * dim
    for i in range(m):
        r[i] = bb[i] * S_K1
    rows.append(r)
    for i in range(m):
        if i == j:
            continue
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        rows.append(r)
    assert len(rows) == dim
    return rows

def babai_nearest_plane(basis_rows, target):
    """Babai nearest-plane using fpylll's exact-ish GSO on an LLL-reduced
    basis.  Returns the lattice vector (list of ints)."""
    dim = len(basis_rows)
    A = IntegerMatrix.from_matrix(basis_rows)
    LLL.reduction(A)
    Bmat = [[A[i][j] for j in range(dim)] for i in range(dim)]
    M = GSO.Mat(A, float_type='mpfr' if False else 'd')
    M.update_gso()
    # mu[i][j] and r[i] = ||b*_i||^2
    t = [float(x) for x in target]
    coeffs = [0] * dim
    # nearest plane, top down
    w = list(t)
    for i in range(dim - 1, -1, -1):
        # <w, b*_i> / ||b*_i||^2  -- reconstruct b*_i from mu and B
        bstar = _gram_schmidt_row(M, Bmat, i, dim)
        rr = sum(x * x for x in bstar)
        if rr == 0:
            continue
        c = sum(w[k] * bstar[k] for k in range(dim)) / rr
        ci = int(math.floor(c + 0.5))
        coeffs[i] = ci
        for k in range(dim):
            w[k] -= ci * Bmat[i][k]
    out = [0] * dim
    for i in range(dim):
        if coeffs[i]:
            for k in range(dim):
                out[k] += coeffs[i] * Bmat[i][k]
    return out

_GS_CACHE = {}

def _gram_schmidt_row(M, Bmat, i, dim):
    key = (id(M), i)
    if key in _GS_CACHE:
        return _GS_CACHE[key]
    # b*_i = b_i - sum_{j<i} mu[i][j] b*_j
    bs = [float(x) for x in Bmat[i]]
    for j in range(i):
        bj = _gram_schmidt_row(M, Bmat, j, dim)
        mu = M.get_mu(i, j)
        for k in range(dim):
            bs[k] -= mu * bj[k]
    _GS_CACHE[key] = bs
    return bs

def run_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    if rows is None:
        return False
    # lattice vector u = (S_K1*(k1_i - A_i), S_K2*k2_i); want it near T
    target = [(-sigs[i]['A'] + k1_bound // 2) * S_K1 for i in range(m)] + \
             [int(k2_bound // 2) * S_K2 for _ in range(m)]
    _GS_CACHE.clear()
    u = babai_nearest_plane(rows, target)
    _GS_CACHE.clear()
    k1s, k2s = [], []
    for i in range(m):
        a = u[i] + sigs[i]['A'] * S_K1
        b = u[m + i]
        if a % S_K1 != 0 or b % S_K2 != 0:
            return False
        k1s.append(a // S_K1); k2s.append(b // S_K2)
    for i in range(m):
        if sigs[i]['B'] % n == 0:
            continue
        d_cand = (modinv(sigs[i]['B'], n)
                  * (k1s[i] + lam * k2s[i] - sigs[i]['A'])) % n
        return d_cand == d_secret
    return False

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_once(curve, m, d_secret, k1_bound, seed, mode):
    """mode in {'orig','proj','projc','cvp'}. Returns (ok, planted_norm, sv_norm).
       orig  = Thread-20 lattice (dim 2m+2, d column present)
       origc = Thread-20 lattice (dim 2m+2) with the CENTERED target
       proj  = d column quotiented out (dim 2m+1), raw target
       projc = d column quotiented out (dim 2m+1), CENTERED target
       cvp   = no embedding, Babai nearest-plane on the centered target"""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

    if mode == 'cvp':
        ok = run_cvp(sigs, n, lam, k1_bound, k2_bound, d_secret)
        return (ok, float('nan'), float('nan'))

    if mode in ('orig', 'origc'):
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound,
                              center=(mode == 'origc'))
        dim = 2 * m + 2
        pv = norm(planted_vector_orig(sigs, d_secret, n, k1_bound, k2_bound,
                                      center=(mode == 'origc')))
    else:
        ctr = (mode == 'projc')
        M, _j = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound,
                                        center=ctr)
        if M is None:
            return None
        dim = 2 * m + 1
        pv = norm(planted_vector_proj(sigs, n, k1_bound, k2_bound, center=ctr))

    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    if mode in ('orig', 'origc'):
        ok = recover_d_orig(reduced, m, n, S_KAN, d_secret) is not None
    else:
        ok = recover_d_proj(reduced, sigs, n, lam, k1_bound, k2_bound,
                            d_secret, center=(mode == 'projc')) is not None
    nz = [norm(r) for r in reduced if any(r)]
    sv = min(nz) if nz else float('nan')
    return (ok, pv, sv)

def rate(curve, m, k1_bound, seeds, mode):
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = run_once(curve, m, d_trial, k1_bound, seed, mode)
        if res is None:
            continue
        ok, pv, sv = res
        wins += bool(ok)
        if pv == pv and pv:
            ratios.append(sv / pv)
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr


# ===========================================================================
if __name__ == '__main__':
    print("=" * 78)
    print("Thread 23 - projected / CVP reformulation of the Phase-2 GLV lattice")
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

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U1: sanity - planted vector lies in the projected lattice")
    print("-" * 78)
    from fractions import Fraction
    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d = random.Random(11).randint(1, n - 1)
        sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, 42)
        rows, j = build_projected_lattice(sigs, n, lam, k1, k2b)
        pv = planted_vector_proj(sigs, n, k1, k2b)
        # membership test: solve rows^T x = pv over Q, check x integral
        Mx = sympy.Matrix(rows).T
        sol = Mx.solve(sympy.Matrix(pv))
        integral = all(sympy.Rational(s).q == 1 for s in sol)
        detA = abs(sympy.Matrix(rows).det())
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
        pred = (S_K1 ** m) * (n ** (m - 1)) * (S_K2 ** m) * S_KAN
        print(f"  {label:18s} m={m:2d} dim={2*m+1:3d} pivot={j} "
              f"planted-in-lattice={integral}  det matches n^(m-1) form={detA == pred}")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U2: sv/pv  (T5 redone) - original vs projected")
    print("-" * 78)
    print(f"  {'curve':18s} {'K1':>3s} {'m':>3s} | {'orig sv/pv':>10s} "
          f"{'win':>5s} | {'proj sv/pv':>10s} {'win':>5s} | "
          f"{'projc sv/pv':>11s} {'win':>5s}")
    for label, curve, k1, m in hist:
        wo, to, ro = rate(curve, m, k1, SEEDS, 'orig')
        wp, tp, rp = rate(curve, m, k1, SEEDS, 'proj')
        wc, tc, rc = rate(curve, m, k1, SEEDS, 'projc')
        print(f"  {label:18s} {k1:3d} {m:3d} | {ro:10.3f} {wo:2d}/{to:<2d} | "
              f"{rp:10.3f} {wp:2d}/{tp:<2d} | {rc:11.3f} {wc:2d}/{tc:<2d}")

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U3: the T4 K1 grid - orig vs proj vs Babai-CVP")
    print("-" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    for label, curve, _k1, m in hist[1:]:
        p, b, n, lam, G = curve
        print(f"\n  {label}   lam*={lam_star(lam, n):.3f}  m={m}  "
              f"K2={math.isqrt(n)+1}")
        hdr = "    " + f"{'mode':6s}" + "".join(f"{k:>7d}" for k in K1_GRID)
        print(hdr)
        for mode in ('orig', 'proj', 'projc', 'cvp'):
            cells = []
            for k1 in K1_GRID:
                w, t, _ = rate(curve, m, k1, SEEDS, mode)
                cells.append(f"{w}/{t}")
            print("    " + f"{mode:6s}" + "".join(f"{c:>7s}" for c in cells))

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U4: T4b control - does more data rescue lam*=0.07 at K1=8?")
    print("-" * 78)
    label, curve, _k1, _m = hist[2]
    print(f"  {label}, K1=8")
    print("    " + f"{'mode':6s}" + "".join(f"{mm:>7d}" for mm in
                                            (8, 12, 16, 24, 32)))
    for mode in ('orig', 'proj', 'projc', 'cvp'):
        cells = []
        for mm in (8, 12, 16, 24, 32):
            w, t, _ = rate(curve, mm, 8, SEEDS, mode)
            cells.append(f"{w}/{t}")
        print("    " + f"{mode:6s}" + "".join(f"{c:>7s}" for c in cells))

    # -----------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP U5: 17-bit population sweep at eff = 0.15 and 0.25")
    print("-" * 78)
    curves17 = search_curves(2**16, 2**17, per_bin=2, nbins=10)
    print(f"  collected {len(curves17)} 17-bit j=0 GLV curves")
    M_SIGS = 12
    for eff in (0.15, 0.25):
        n_orig = n_proj = n_projc = tot = 0
        detail = []
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            cur = (p, b, n, lam, G)
            wo, t, _ = rate(cur, M_SIGS, k1b, SEEDS, 'orig')
            wp, _, _ = rate(cur, M_SIGS, k1b, SEEDS, 'proj')
            wc, _, _ = rate(cur, M_SIGS, k1b, SEEDS, 'projc')
            tot += 1
            n_orig += (wo == t)
            n_proj += (wp == t)
            n_projc += (wc == t)
            detail.append((lam_star(lam, n), wo, wp, wc, t))
        print(f"\n  eff={eff}:  full-recovery curves  orig={n_orig}/{tot}  "
              f"proj={n_proj}/{tot}  projc={n_projc}/{tot}")
        print(f"    {'lam*':>7s} {'orig':>7s} {'proj':>7s} {'projc':>7s}")
        for ls, wo, wp, wc, t in sorted(detail):
            print(f"    {ls:7.4f} {str(wo)+'/'+str(t):>7s} "
                  f"{str(wp)+'/'+str(t):>7s} {str(wc)+'/'+str(t):>7s}")

    print("\n" + "-" * 78)
    print("EXP U6: 2x2 factorial - {projection} x {centering}")
    print("-" * 78)
    print("  Isolates which of the two changes carries the effect.")
    for label, curve, _k1, m in hist[1:]:
        p, b, n, lam, G = curve
        print(f"\n  {label}   lam*={lam_star(lam, n):.3f}  m={m}")
        print("    " + f"{'mode':7s}" + "".join(f"{k:>7d}" for k in K1_GRID))
        for mode, desc in (('orig', 'proj- ctr-'), ('proj', 'proj+ ctr-'),
                           ('origc', 'proj- ctr+'), ('projc', 'proj+ ctr+')):
            cells = []
            for k1 in K1_GRID:
                w, t, _ = rate(curve, m, k1, SEEDS, mode)
                cells.append(f"{w}/{t}")
            print("    " + f"{mode:7s}" + "".join(f"{c:>7s}" for c in cells)
                  + f"   [{desc}]")

    print("\n" + "-" * 78)
    print("EXP U7: where is the NEW population wall?  (projc, m=12, 17-bit)")
    print("-" * 78)
    print(f"  {'eff':>5s} {'projc full-recovery':>22s} {'mean win rate':>15s}")
    for eff in (0.15, 0.25, 0.35, 0.45, 0.60):
        full = tot = wins = trials = 0
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            w, t, _ = rate((p, b, n, lam, G), 12, k1b, SEEDS, 'projc')
            full += (w == t); tot += 1; wins += w; trials += t
        print(f"  {eff:5.2f} {str(full)+'/'+str(tot):>22s} {wins/trials:15.2f}")

    print("\n" + "=" * 78)
    print("done")
