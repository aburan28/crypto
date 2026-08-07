"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the target is findable.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 Kannan-embedding lattice L_K (dim 2m+2) always contains the
  trivial vector  n*S_D*e_m , of norm n*S_D.  Since the planted vector has
  norm ~ n*sqrt(2m/3 + 4/3) >> n*S_D, the planted vector is NEVER lambda_1.
  Concretely lambda_1(L_K) = n*S_D for every curve tested, and sv/pv sat in
  [0.337, 0.368] for successes and failures alike.  So Kannan-SVP is the
  wrong tool: it is structurally unable to return the planted vector, and
  any success it had came from the planted vector merely appearing somewhere
  in the LLL-reduced basis rather than at position 1.

  Why the trivial vector exists: n*(row_d) - sum_i B_i*(row_i) = n*S_D*e_m.
  It is unavoidable because d is only defined mod n.  It is also harmless
  for the *decision* problem: w and w + n*e_m give the same d mod n.  It is
  fatal only for the SVP/Kannan formulation.

This script implements the 2026-07-29 "next step proposal" verbatim:
  drop the Kannan row/column and solve BDD directly by CVP on the
  (2m+1)-dimensional lattice.

  Lattice L (dim 2m+1), columns [k1_0..k1_{m-1} | d | k2_0..k2_{m-1}]:
      row i       (i<m):  n*S_K1 * e_i
      row m            :  (B_i*S_K1)_{i<m}  ,  S_D at col m
      row m+1+i (i<m)  :  -lam*S_K1 at col i,  S_K2 at col m+1+i

  From  A_i + B_i*d = k1_i + lam*k2_i  (mod n)  the lattice vector
      w = d*row_m - sum_i c_i*row_i + sum_i k2_i*row_{m+1+i}
  has coordinates  ( (k1_i - A_i)*S_K1 , d*S_D , k2_i*S_K2 ), hence with
      t = ( -A_i*S_K1 , 0 , 0 )                     [UNCENTERED target]
  we get  w - t = ( k1_i*S_K1 , d*S_D , k2_i*S_K2 ),  a BDD instance.
  Recovery: d = w[m] / S_D  mod n.

  CENTERED target: every unknown is uniform on a known box, so shifting the
  target to the box centre,
      t' = t + ( (K1/2)*S_K1 , n/2 , (K2/2)*S_K2 ),
  makes the error  ( (k1_i-K1/2)*S_K1 , d-n/2 , (k2_i-K2/2)*S_K2 ), whose
  expected norm is smaller by exactly sqrt(4) = 2  (uniform on [-c/2,c/2] has
  E[x^2] = c^2/12 versus c^2/3 for [0,c]).  Recovery is unchanged.

  PREDICTION (the falsifier for this run).  det(L) = (n*S_K1)^m * S_D * S_K2^m
  = n^{3m} / (K1*K2)^m, so the BDD radius scales as
      lambda_1 ~ det^{1/(2m+1)} ~ n * (K1*K2/n)^{-m/(2m+1)} ~ n*sqrt(n/(K1*K2)).
  Halving the error norm therefore buys a factor ~4 in the tolerable K1.
  So if the reformulation is real, the K1 wall measured on 2026-07-29
  (K1 ~ 12-16 for lam*=0.340, K1 ~ 4-6 for lam*=0.070) should move outward
  by ~4x under the centered CVP, and NOT move at all if the 2026-07-29 wall
  was information-theoretic.

Experiments:
  C1  sanity: the planted vector really is the CVP solution (exact error
      norms, and Babai's decoding radius (1/2)*min_i ||b*_i|| for comparison)
  C2  head-to-head K1 sweep, Kannan-SVP vs CVP-uncentered vs CVP-centered vs
      BKZ(20)+CVP-centered, on the two historical 12-bit curves (the exact
      T4 grid of 2026-07-29, so the numbers are directly comparable)
  C3  does more data rescue the wall?  m sweep at the first failing K1
      (the CVP analogue of 2026-07-29 T4b)
  C4  17-bit generalisation: fresh curves, K1 wall under SVP vs centered CVP

Run: python3 glv_hnp_phase2_cvp.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (a=0 short Weierstrass) -- identical to
# glv_hnp_phase2_lambda_threshold.py so results are directly comparable.
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
    mm, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        bb = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, bb * bb % p, t * bb * bb % p, r * bb % p

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

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    """Smaller root of x^2+x+1 = 0 mod n (so lam <= n//2), or None."""
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, r2)

def lam_star(lam, n):
    return min(lam % n, (-lam) % n) / n

# ---------------------------------------------------------------------------
# Scales, signatures  (verbatim from glv_hnp_phase2_lambda_threshold.py:227)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- unchanged from prior work."""
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
# V0: the historical Kannan-embedding lattice (dim 2m+2)
# ---------------------------------------------------------------------------

def build_kannan_lattice(sigs, n, lam, k1_bound, k2_bound):
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

def recover_d_kannan(reduced, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# V1/V2: the CVP reformulation (dim 2m+1, no Kannan row/column)
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """dim 2m+1; columns [k1_0..k1_{m-1} | d | k2_0..k2_{m-1}]."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    return M

def cvp_target(sigs, n, k1_bound, k2_bound, centered):
    """t (uncentered) or t' (centered).  See module docstring."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    t = [0] * dim
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    if centered:
        c1, cd, c2 = centres(n, k1_bound, k2_bound)
        for i in range(m):
            t[i] += c1
            t[m + 1 + i] = c2
        t[m] = cd
    return t

def centres(n, k1_bound, k2_bound):
    """Box centres in SCALED coordinates.  Centring in unscaled units
    ((K1//2)*S_K1) under-centres badly for small K1 -- e.g. K1=2 sends
    k1 in {0,1} to {-1,0}, which is no gain at all."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    return ((k1_bound * S_K1) // 2, (n * S_D) // 2, (k2_bound * S_K2) // 2)

def cvp_error_true(sigs, d_secret, n, k1_bound, k2_bound, centered):
    """The exact error vector w - t that the CVP solver should find."""
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    e = [0] * (2 * m + 1)
    c1, cd, c2 = centres(n, k1_bound, k2_bound) if centered else (0, 0, 0)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1 - c1
        e[m + 1 + i] = sigs[i]['k2'] * S_K2 - c2
    e[m] = d_secret * S_D - cd
    return e

# ---------------------------------------------------------------------------
# Gram-Schmidt + Babai nearest plane (exact integer basis, float GSO)
# ---------------------------------------------------------------------------

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def gram_schmidt(B):
    """Returns (Bstar, sqnorms) as floats.  B is a list of integer rows."""
    Bstar, sq = [], []
    for row in B:
        v = [float(x) for x in row]
        for bs, s in zip(Bstar, sq):
            if s == 0.0:
                continue
            mu = dot([float(x) for x in row], bs) / s
            v = [vi - mu * bi for vi, bi in zip(v, bs)]
        Bstar.append(v)
        sq.append(dot(v, v))
    return Bstar, sq

def babai_nearest_plane(B, t):
    """Return (w, e) with w in L(B) and e = t - w.  B should be LLL-reduced."""
    Bstar, sq = gram_schmidt(B)
    b = [float(x) for x in t]
    coeffs = [0] * len(B)
    for i in range(len(B) - 1, -1, -1):
        if sq[i] == 0.0:
            continue
        c = int(round(dot(b, Bstar[i]) / sq[i]))
        coeffs[i] = c
        if c:
            b = [bj - c * Bij for bj, Bij in zip(b, B[i])]
    w = [0] * len(t)
    for c, row in zip(coeffs, B):
        if c:
            for j, x in enumerate(row):
                w[j] += c * x
    e = [ti - wi for ti, wi in zip(t, w)]
    return w, e

def babai_round(B, t):
    """Babai rounding (solve t = x*B over R, round x).  Weaker than NP."""
    dim = len(B)
    # Gaussian elimination on the transpose, in floats.
    A = [[float(B[i][j]) for i in range(dim)] for j in range(dim)]  # A x = t
    rhs = [float(x) for x in t]
    for col in range(dim):
        piv = max(range(col, dim), key=lambda r: abs(A[r][col]))
        if abs(A[piv][col]) < 1e-9:
            return None, None
        A[col], A[piv] = A[piv], A[col]
        rhs[col], rhs[piv] = rhs[piv], rhs[col]
        pv = A[col][col]
        for r in range(col + 1, dim):
            f = A[r][col] / pv
            if f:
                for c2 in range(col, dim):
                    A[r][c2] -= f * A[col][c2]
                rhs[r] -= f * rhs[col]
    x = [0.0] * dim
    for r in range(dim - 1, -1, -1):
        s = rhs[r] - sum(A[r][c] * x[c] for c in range(r + 1, dim))
        x[r] = s / A[r][r]
    coeffs = [int(round(v)) for v in x]
    w = [0] * dim
    for c, row in zip(coeffs, B):
        if c:
            for j, v in enumerate(row):
                w[j] += c * v
    return w, [ti - wi for ti, wi in zip(t, w)]

def min_gs_norm(B):
    _, sq = gram_schmidt(B)
    return math.sqrt(min(s for s in sq if s > 0))

# ---------------------------------------------------------------------------
# Single trials
# ---------------------------------------------------------------------------

def lll_rows(M, beta=None):
    dim = len(M)
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def trial_kannan(curve, m, d_secret, k1_bound, seed, beta=None):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m:
        return None
    M = build_kannan_lattice(sigs, n, lam, k1_bound, k2)
    red = lll_rows(M, beta)
    _, _, _, S_KAN = scales(n, k1_bound, k2)
    return recover_d_kannan(red, m, n, S_KAN, d_secret)

def trial_cvp(curve, m, d_secret, k1_bound, seed, centered, beta=None,
              method='np', diag=False):
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2, seed)
    if len(sigs) < m:
        return None
    M = build_cvp_lattice(sigs, n, lam, k1_bound, k2)
    red = lll_rows(M, beta)
    t = cvp_target(sigs, n, k1_bound, k2, centered)
    if method == 'round':
        w, e = babai_round(red, t)
        if w is None:
            return None
    else:
        w, e = babai_nearest_plane(red, t)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2)
    # Centering shifts the TARGET, not the lattice: the sought lattice vector w
    # is the same in both cases and always has w[m] = d*S_D.  (An earlier
    # version added the centre offset back here, which is wrong.)
    d_cand = (w[m] // S_D) % n
    ok = (d_cand == d_secret)
    if not diag:
        return ok
    e_true = cvp_error_true(sigs, d_secret, n, k1_bound, k2, centered)
    return {
        'ok': ok,
        'e_found': norm(e),
        'e_true': norm(e_true),
        'babai_radius': 0.5 * min_gs_norm(red),
        'dim': 2 * m + 1,
    }

def rate(fn, curve, m, k1_bound, seeds, **kw):
    wins = 0
    for s in seeds:
        d = random.Random(s + 7777).randint(1, curve[2] - 1)
        r = fn(curve, m, d, k1_bound, s, **kw)
        if r:
            wins += 1
    return wins

# ---------------------------------------------------------------------------
# Curve construction
# ---------------------------------------------------------------------------

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        G = find_generator(p, b, n, seed=seed)
        if G is not None:
            return (p, b, n, glv_eigenvalue(n), G)
    return None

def search_curves(lo, hi, want, lam_lo=0.0, lam_hi=0.5):
    """Fresh j=0 GLV curves with prime n in [lo,hi]."""
    out = []
    p = sympy.nextprime(lo)
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 2 or not sympy.isprime(nc) or nc % 3 != 1:
                        continue
                    lam = glv_eigenvalue(nc)
                    if lam is None:
                        continue
                    ls = lam_star(lam, nc)
                    if not (lam_lo <= ls <= lam_hi):
                        continue
                    c = build_curve(p, nc)
                    if c and c[3] is not None:
                        out.append(c)
                        break
        p = sympy.nextprime(p)
    return out

# ===========================================================================
# MAIN
# ===========================================================================

print("=" * 78)
print("Thread 23: Phase-2 GLV lattice -- Kannan-SVP vs explicit CVP (BDD)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

# Exact historical curves from RESEARCH_AUTOLAB_LOG.md (2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,             p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",      2677, 2, 2647, 185,  8,  10),
]
curves = []
for label, p, b, n, lam, k1h, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    curves.append((label, (p, b, n, lam, G), k1h, m))

print("\nCurves (lam* = min(lam, n-lam)/n, the representation-invariant form):")
for label, c, k1h, m in curves:
    p, b, n, lam, G = c
    print(f"  {label:<14} p={p:<6} n={n:<6} lam={lam:<6} lam*={lam_star(lam,n):.4f}  m={m}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C1: is the planted vector the CVP solution?  (m and K1 as in T4)")
print("-" * 78)
print("e_true  = ||exact error||;  e_found = ||Babai error||;  radius = (1/2)min||b*_i||")
print("Babai nearest-plane is GUARANTEED correct when e_true < radius; it usually")
print("succeeds well beyond that, so radius is a sufficient-condition yardstick only.\n")
print(f"{'curve':<14} {'K1':>3} {'ctr':>4} {'dim':>4} {'e_true':>11} {'e_found':>11} "
      f"{'radius':>11} {'e_true/rad':>10} {'ok':>4}")
for label, c, k1h, m in curves:
    for k1 in (2, k1h):
        for centered in (False, True):
            d = random.Random(42 + 7777).randint(1, c[2] - 1)
            r = trial_cvp(c, m, d, k1, 42, centered, diag=True)
            print(f"{label:<14} {k1:>3} {str(centered):>5} {r['dim']:>4} "
                  f"{r['e_true']:>11.1f} {r['e_found']:>11.1f} {r['babai_radius']:>11.1f} "
                  f"{r['e_true']/r['babai_radius']:>10.3f} {str(r['ok']):>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C2: head-to-head K1 sweep -- the exact T4 grid of 2026-07-29")
print("-" * 78)
print("V0 = Kannan-SVP (historical);  V1 = CVP nearest-plane, uncentered;")
print("V2 = CVP nearest-plane, centered;  V3 = BKZ(20) + CVP centered;")
print("V2r = Babai *rounding*, centered (weaker rule, to isolate NP's contribution).")
print(f"Each cell = wins out of {len(SEEDS)} seeds.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

c2_summary = {}
for label, c, k1h, m in curves:
    n = c[2]
    print(f"\n  {label}   (lam*={lam_star(c[3],n):.4f}, m={m}, "
          f"K2={math.isqrt(n)+1}, dim_cvp={2*m+1})")
    hdr = "    " + " ".join(f"{k:>4}" for k in K1_GRID)
    print(f"    {'K1':<5}" + hdr[4:])
    rows = {}
    for name, fn, kw in (
        ("V0",  trial_kannan, {}),
        ("V1",  trial_cvp,    {'centered': False}),
        ("V2",  trial_cvp,    {'centered': True}),
        ("V3",  trial_cvp,    {'centered': True, 'beta': 20}),
        ("V2r", trial_cvp,    {'centered': True, 'method': 'round'}),
    ):
        res = [rate(fn, c, m, k1, SEEDS, **kw) for k1 in K1_GRID]
        rows[name] = res
        print(f"    {name:<5}" + " ".join(f"{v:>4}" for v in res))
    # wall = largest K1 with >= 3/5 (majority) recovery
    def wall(res):
        w = 0
        for k1, v in zip(K1_GRID, res):
            if v >= 3:
                w = k1
        return w
    walls = {k: wall(v) for k, v in rows.items()}
    c2_summary[label] = walls
    print("    walls (largest K1 with >=3/5): " +
          ", ".join(f"{k}={v}" for k, v in walls.items()))

print("\n  K1-wall summary (>=3/5 majority recovery):")
print(f"    {'curve':<14} {'V0 SVP':>8} {'V1 CVP':>8} {'V2 ctr':>8} {'V3 BKZ':>8} {'V2/V0':>8}")
for label, w in c2_summary.items():
    ratio = (w['V2'] / w['V0']) if w['V0'] else float('inf')
    print(f"    {label:<14} {w['V0']:>8} {w['V1']:>8} {w['V2']:>8} {w['V3']:>8} "
          f"{ratio:>8.2f}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C3: does more data rescue the CVP wall?  (analogue of T4b)")
print("-" * 78)
print("At the first K1 beyond the V2 wall, sweep m.  If recovery stays at 0 the")
print("wall is information-theoretic; if it climbs, it is a data wall.\n")

for label, c, k1h, m in curves:
    w2 = c2_summary[label]['V2']
    k1_over = next((k for k in K1_GRID if k > w2), None)
    if k1_over is None:
        print(f"  {label}: V2 wall beyond the grid; skipped.")
        continue
    ms = [4, 6, 8, 12, 16, 24, 32]
    res = [rate(trial_cvp, c, mm, k1_over, SEEDS, centered=True) for mm in ms]
    print(f"  {label}  K1={k1_over} (first beyond V2 wall {w2}):")
    print("    m      " + " ".join(f"{v:>4}" for v in ms))
    print("    V2     " + " ".join(f"{v:>4}" for v in res))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C4: 17-bit generalisation -- fresh curves, V0 vs V2 walls")
print("-" * 78)

fresh = search_curves(2**16, 2**17, 6)
print(f"  found {len(fresh)} fresh 17-bit j=0 GLV curves\n")
# 17-bit curves have K2 ~ 256, so eff = K1*K2/n reaches 1 (no bias at all) only
# near K1 ~ 256.  The K1_GRID of C2 tops out at 64 = eff 0.25, which CENSORS the
# V2 wall for these curves; extend it here so the wall is observed, not clipped.
K1_GRID_C4 = [2, 4, 8, 16, 32, 48, 64, 96, 128, 192, 254]
print(f"  K1 grid: {K1_GRID_C4}  (eff = K1*K2/n; '>=254' means wall not reached)\n")
print(f"  {'p':>7} {'n':>7} {'lam*':>7} {'m':>3} | {'V0 wall':>8} {'V2 wall':>8} {'ratio':>7}")
ratios = []
for c in fresh:
    p, b, n, lam, G = c
    m = 12
    def wall_of(fn, **kw):
        w = 0
        for k1 in K1_GRID_C4:
            if rate(fn, c, m, k1, SEEDS, **kw) >= 3:
                w = k1
        return w
    w0 = wall_of(trial_kannan)
    w2 = wall_of(trial_cvp, centered=True)
    r = (w2 / w0) if w0 else float('inf')
    ratios.append(r)
    print(f"  {p:>7} {n:>7} {lam_star(lam,n):>7.4f} {m:>3} | {w0:>8} {w2:>8} {r:>7.2f}")

fin = [r for r in ratios if math.isfinite(r)]
if fin:
    print(f"\n  median V2/V0 wall ratio: {sorted(fin)[len(fin)//2]:.2f}  "
          f"(prediction from det(L): ~4)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
