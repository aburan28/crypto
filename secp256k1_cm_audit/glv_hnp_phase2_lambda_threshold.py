"""
GLV-HNP Phase 2, Thread 20: is lambda/n really the predictor of LLL success?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-26 run #1):
  The Phase-2 GLV lattice attack recovered d for curves with lam/n in
  {0.53, 0.66, 0.34} but failed (LLL and BKZ-40 alike) for p=2677,
  n=2647, lam=185 (lam/n=0.07).  That run proposed "lam/n threshold".
  But an earlier run (2026-06-27, Exp M/N/O, log line ~3571) had already
  FALSIFIED lam/n as a continuous predictor for the Phase-1 lattice.

This script reconciles the two by testing a sharper hypothesis.

  H20:  LLL success is governed NOT by lam/n but by the shortest vector of
        the 2-dimensional "lambda block" sublattice
            B_i = < (n*S_K1, 0),  (-lam*S_K1, S_K2) >
        which lives in the (col i, col m+1+i) coordinate plane of the
        Phase-2 lattice.  Write mu = |shortest vector of B_i| (exact, via
        Gauss/Lagrange reduction).  The planted vector is only findable if
        mu is not much shorter than the planted vector's norm, i.e. if
            rho = mu / ||v_planted||   is >= ~1.

  Two immediate consequences, both testable:
   (T1) The lattice depends on lam only mod n, so replacing lam by lam - n
        (i.e. "lam/n = 0.07" -> "lam/n = -0.93") yields the SAME lattice and
        the SAME planted vector.  If lam/n were causal, this substitution
        would change the outcome.  It cannot.  => lam/n is not causal.
   (T2) mu is a continued-fraction-like quantity of (lam, n, S_K1, S_K2),
        not a monotone function of lam/n.  This explains why lam/n
        "worked" on 4 data points and failed on 50.

Experiments:
  T1  representation invariance (lam vs lam-n vs the other root n-1-lam)
  T2  mu / rho for the four historical Phase-2 curves
  T3  17-bit sweep, ~20 curves with lam* = min(lam,n-lam)/n spread over
      (0, 0.5); compare predictive accuracy of a lam* threshold vs a rho
      threshold.  BKZ(20) retry on every LLL failure.

RESULTS (2026-07-29 run; full output reproducible with the command below):

  T1  CONFIRMED.  HNF(lattice with lam) == HNF(lattice with lam-n) for all
      three historical curves, and the LLL success rates are identical
      (5/5 vs 5/5, 5/5 vs 5/5, 0/5 vs 0/5).  lam/n is therefore NOT a causal
      variable — only lam mod n is, and lam* = min(lam,n-lam)/n is the same
      for both roots of x^2+x+1.

  T2/T3  H20 FALSIFIED.  rho does not separate success from failure; the
      2026-07-26 failure curve has the LARGEST rho of the three (0.686 vs
      0.610 / 0.399).  In the 17-bit sweep the best rho threshold never beats
      the majority baseline.

  T3  The "lam* threshold" claim of 2026-07-26 is FALSIFIED.  At eff=0.05 all
      20 curves recover 5/5 (or 4/5), including lam*=0.0068 and lam*=0.027 —
      far below the 0.07 that was reported as a hard failure.  lam* is not a
      viability threshold at all.

  T4  The real variable is the bias strength eff = K1*K2/n, not lam.  Both
      historical 12-bit curves recover 5/5 for K1 <= 4.  lam* shifts the K1
      wall (lam*=0.34 -> wall at K1 ~ 12-16; lam*=0.07 -> wall at K1 ~ 4-6),
      a factor of ~3, but creates no hard structural obstruction.
      T4b: at K1=8 the lam*=0.07 curve is not rescued by more data
      (m = 8/12/16/24/32 -> 0,0,1,0,1 out of 5), so the K1 wall is genuine;
      it is simply a K1 wall and not a lam wall.

  T5  STRUCTURAL: the shortest vector of the Phase-2 lattice is always the
      trivial vector n*S_D*e_m (100% of its energy in the d column, |sv[m]|/n
      = 1.0000 exactly, last coordinate 0).  It is ~2-3x shorter than the
      planted vector on every curve tested, success and failure alike
      (sv/pv in [0.34, 0.61]).  The planted vector is therefore NEVER lambda_1
      of this lattice; recovery is a BDD/coset condition (shortest vector
      among those with last coordinate +-S_KANNAN), not an SVP condition.
      This explains the long run of failed closed-form separators logged
      between 2026-06-21 and 2026-06-29: no curve-level invariant can predict
      a BDD success on a (2m+1)-dimensional projected lattice.

Run: python3 glv_hnp_phase2_lambda_threshold.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)  -- same as glv_hnp_phase2_20bit.py
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

# ---------------------------------------------------------------------------
# CM theory for j=0 curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p)."""
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
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n.  Requires n = 1 mod 3."""
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
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    """Symmetric size of lam: min(lam, n-lam)/n in [0, 0.5]."""
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# API used by Threads 20b/20c/20d and 23.
#
# 2026-08-02: the 20b/20c/20d scripts committed in e845207 import these four
# names plus a 9-argument run_experiment from this module, and none of them
# existed here -- all three died with AttributeError on import, so the
# 2026-07-29 results were not reproducible from the committed tree.  Restored
# below; see RESEARCH_AUTOLAB_LOG.md 2026-08-02 "Repair" for the diagnosis.
# ---------------------------------------------------------------------------

def glv_eigenvalues(n):
    """Both roots of x^2+x+1 = 0 mod n, ascending.  Alias of glv_roots."""
    return glv_roots(n)

def mu_of(lam, n):
    """mu = min(lam, n-lam)/n in (0, 1/2].  Alias of lam_star."""
    return lam_star(lam, n)

def identify_twist(p, n, seed=12345):
    """The sextic-twist coefficient b with #E(F_p): y^2=x^3+b equal to n.
    n must be prime, so a single point of order n identifies the twist."""
    rng = random.Random((seed * 1000003) ^ (p << 1) ^ n)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            return b
    return None

def rival_sublattice_nu(n, lam, k1_bound, k2_bound):
    """(lambda_1(L2), sqrt(det L2), nu_hat) for the non-planted 2D sublattice

        L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >,   det L2 = n*S_K1*S_K2.

    nu_hat = lambda_1(L2)/sqrt(det L2) is scale-free and, since det L2 does not
    depend on lam, isolates the lam-dependence of the block geometry.
    Ratios go through math.log so 256-bit inputs do not overflow float."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    nsq = w[0] * w[0] + w[1] * w[1]
    det = n * S_K1 * S_K2
    lam1 = math.exp(0.5 * math.log(nsq)) if nsq else 0.0
    root_det = math.exp(0.5 * math.log(det))
    nu_hat = math.exp(0.5 * (math.log(nsq) - math.log(det))) if nsq else 0.0
    return lam1, root_det, nu_hat

# ---------------------------------------------------------------------------
# The lambda-block 2D sublattice and its exact shortest vector
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Lagrange/Gauss reduction of a 2D integer lattice.  Returns the exact
    shortest nonzero vector (u, v are tuples of ints)."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        # mu = round(<v,u>/<u,u>), exact integer rounding (round-half-away)
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def lambda_block_mu(n, lam, S_K1, S_K2):
    """Exact |shortest vector| of the 2D block < (n*S_K1,0), (-lam*S_K1,S_K2) >."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w

def planted_norm_expected(m, n, k1_bound, k2_bound, S_K1, S_D, S_K2, S_KANNAN):
    """E[||v_planted||] with k1~U[0,K1), k2~U[0,K2), d~U[0,n).  E[x^2]=X^2/3."""
    return math.sqrt(
        m * (k1_bound * S_K1) ** 2 / 3.0
        + (n * S_D) ** 2 / 3.0
        + m * (k2_bound * S_K2) ** 2 / 3.0
        + S_KANNAN ** 2
    )

# ---------------------------------------------------------------------------
# Phase-2 lattice (verbatim construction from glv_hnp_phase2_20bit.py)
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
    """lam is used verbatim (NOT reduced mod n) so that representation
    invariance can be tested."""
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

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound):
    """The exact planted lattice vector (k1_i*S_K1, d*S_D, k2_i*S_K2, S_KANNAN)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(x * x for x in v))


def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
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

def run_experiment(*a, **kw):
    """Two calling conventions, both in use across Threads 20a-20d and 23:

      run_experiment(curve, m, d, k1_bound, seed=..., lam_override=...)
          curve = (p, b, n, lam, G); returns (ok, planted_norm, shortest_norm).
      run_experiment(p, n, lam, G, m, d, k1_bound, k2_bound, seed)
          the flat form used by 20b/20c/20d; returns bool.  lam is used for
          BOTH signature generation and the lattice, which is what the Arm-5
          synthetic-lambda design requires.
    """
    if len(a) >= 9 and isinstance(a[0], int):
        p, n, lam, G, m, d, k1_bound, k2_bound, seed = a[:9]
        res = _run_experiment((p, None, n, lam, G), m, d, k1_bound, seed,
                              k2_override=k2_bound, **kw)
        return bool(res is not None and res[0])
    return _run_experiment(*a, **kw)

def _run_experiment(curve, m, d_secret, k1_bound, seed=42, lam_override=None,
                    use_bkz=False, bkz_beta=20, k2_override=None):
    """Returns (recovered: bool, planted_norm, shortest_reduced_norm)."""
    p, b, n, lam, G = curve
    lam_lat = lam if lam_override is None else lam_override
    k2_bound = math.isqrt(n) + 1 if k2_override is None else k2_override
    # signatures always use the canonical lam (the true GLV eigenvalue);
    # only the lattice representation varies.
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam_lat, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KANNAN = scales(n, k1_bound, k2_bound)
    ok = recover_d(reduced, m, n, S_KANNAN, d_secret) is not None
    pn = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
    sn = min(norm(r) for r in reduced)
    return (ok, pn, sn)

def success_rate(curve, m, k1_bound, seeds, lam_override=None,
                 use_bkz=False, bkz_beta=20):
    """Returns (wins, total, mean(shortest/planted) norm ratio)."""
    p, b, n, lam, G = curve
    wins = 0
    ratios = []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = run_experiment(curve, m, d_trial, k1_bound, seed,
                             lam_override=lam_override,
                             use_bkz=use_bkz, bkz_beta=bkz_beta)
        if res is None:
            continue
        ok, pn, sn = res
        wins += bool(ok)
        ratios.append(sn / pn if pn else float('nan'))
    mean_ratio = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mean_ratio

def hnf_fingerprint(M):
    """Row-style Hermite normal form via sympy, as a lattice fingerprint."""
    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form
    return hermite_normal_form(Matrix(M).T).T

# ---------------------------------------------------------------------------
# Curve construction helpers
# ---------------------------------------------------------------------------

def build_curve(p, n, seed=12345):
    """Find the sextic twist b with #E = n, plus a generator."""
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    """Collect j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3,
    bucketed by lam* = min(lam,n-lam)/n into `nbins` bins over [0,0.5]."""
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
                    ls = lam_star(lam, n_cand)
                    idx = min(nbins - 1, int(ls / (0.5 / nbins)))
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


def _driver():
    """The T1-T5 experiment suite.  Runs only on direct execution.

    2026-08-02: this block used to sit at module level, so Threads 20b/20c/20d
    (which import this file as a library) re-ran the entire suite, curve search
    included, on every import.  Wrapped without other changes.
    """
    # ===========================================================================
    # EXPERIMENT T1 — representation invariance
    # ===========================================================================

    print("=" * 78)
    print("Thread 20 — is lam/n the predictor?  (GLV-HNP Phase 2)")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    # Historical Phase-2 curves (from RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26)
    HIST = [
        # label,             p,    b, n,    lam,  K1, m
        ("8-bit/199",        211,  2, 199,  106,  2,  6),
        ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
        ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
    ]

    hist_curves = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist_curves.append((label, (p, b, n, lam, G), k1, m))

    print("\n" + "-" * 78)
    print("EXP T1: representation invariance of the lattice under lam -> lam - n")
    print("-" * 78)
    print("The Phase-2 lattice contains the rows n*S_K1*e_i, so replacing lam by")
    print("lam - n is a unimodular row operation: SAME lattice, SAME planted vector.")
    print("If lam/n were causal, 0.07 vs -0.93 would behave differently.\n")

    print(f"{'curve':<18} {'lam/n':>8} {'(lam-n)/n':>10} {'LLL(lam)':>9} "
          f"{'LLL(lam-n)':>11} {'HNF equal':>10}")
    for label, curve, k1, m in hist_curves:
        p, b, n, lam, G = curve
        w1, t1, _ = success_rate(curve, m, k1, SEEDS)
        w2, t2, _ = success_rate(curve, m, k1, SEEDS, lam_override=lam - n)
        # lattice identity check on one instance
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        k2b = math.isqrt(n) + 1
        sg = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
        same = (hnf_fingerprint(build_glv_lattice(sg, n, lam, k1, k2b))
                == hnf_fingerprint(build_glv_lattice(sg, n, lam - n, k1, k2b)))
        print(f"{label:<18} {lam/n:>8.4f} {(lam-n)/n:>10.4f} "
              f"{str(w1)+'/'+str(t1):>9} {str(w2)+'/'+str(t2):>11} {str(same):>10}")

    print("\nAlso: the two roots of x^2+x+1 mod n are lam and n-1-lam, so")
    print("lam* = min(lam, n-lam)/n is the SAME for both roots — any genuine")
    print("predictor must be a function of lam*, not of lam/n.")
    print(f"{'curve':<18} {'root1/n':>9} {'root2/n':>9} {'lam*':>8}")
    for label, curve, k1, m in hist_curves:
        p, b, n, lam, G = curve
        r1, r2 = glv_roots(n)
        print(f"{label:<18} {r1/n:>9.4f} {r2/n:>9.4f} {lam_star(lam, n):>8.4f}")


    # ===========================================================================
    # EXPERIMENT T2 — mu and rho for the historical curves
    # ===========================================================================

    print("\n" + "-" * 78)
    print("EXP T2: lambda-block shortest vector mu, and rho = mu / ||v_planted||")
    print("-" * 78)

    print("sv/pv = ||shortest vector after LLL|| / ||planted vector||;")
    print("       sv/pv < 1 means a NON-planted vector is shorter (structural failure).\n")
    print(f"{'curve':<18} {'lam*':>7} {'K1':>4} {'K2':>5} {'eff':>6} "
          f"{'mu':>12} {'||v||':>12} {'rho':>7} {'sv/pv':>7} {'LLL':>6}")
    t2_rows = []
    for label, curve, k1, m in hist_curves:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
        mu, wvec = lambda_block_mu(n, lam, S_K1, S_K2)
        vn = planted_norm_expected(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN)
        rho = mu / vn
        w, t, svpv = success_rate(curve, m, k1, SEEDS)
        eff = k1 * k2 / n
        t2_rows.append((label, lam_star(lam, n), rho, w == t))
        print(f"{label:<18} {lam_star(lam,n):>7.4f} {k1:>4} {k2:>5} {eff:>6.3f} "
              f"{mu:>12.1f} {vn:>12.1f} {rho:>7.3f} {svpv:>7.3f} {str(w)+'/'+str(t):>6}")


    # ===========================================================================
    # EXPERIMENT T3 — 17-bit sweep with lam* spread
    # ===========================================================================

    print("\n" + "-" * 78)
    print("EXP T3: 17-bit sweep, lam* spread over (0, 0.5)")
    print("-" * 78)

    LO, HI = 2**16, 2**17
    M_SIGS = 12

    print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
    curves = search_curves(LO, HI, per_bin=2, nbins=10)
    print(f"Found {len(curves)} curves.\n")

    def sweep_at_eff(curves, eff_target, m_sigs, seeds):
        print(f"\n=== eff target = {eff_target}  (m = {m_sigs}, {len(seeds)} seeds) ===")
        print(f"{'p':>7} {'n':>7} {'lam':>7} {'lam*':>6} {'K1':>4} {'eff':>6} "
              f"{'m_thr':>6} {'rho':>7} {'sv/pv':>7} {'LLL':>6} {'BKZ20':>6}")
        rows = []
        for (p, b, n, lam, G) in curves:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_target * n / k2))
            eff = k1 * k2 / n
            m_thr = math.log(n) / math.log(1.0 / eff) if eff < 1 else float('inf')
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
            mu, _ = lambda_block_mu(n, lam, S_K1, S_K2)
            vn = planted_norm_expected(m_sigs, n, k1, k2, S_K1, S_D, S_K2, S_KAN)
            rho = mu / vn
            curve = (p, b, n, lam, G)
            w, t, svpv = success_rate(curve, m_sigs, k1, seeds)
            if w == t:
                bkz_str = "-"
            else:
                wb, tb, _ = success_rate(curve, m_sigs, k1, seeds,
                                         use_bkz=True, bkz_beta=20)
                bkz_str = f"{wb}/{tb}"
            ls = lam_star(lam, n)
            rows.append({'p': p, 'n': n, 'lam': lam, 'lam_star': ls, 'rho': rho,
                         'svpv': svpv, 'wins': w, 'total': t, 'ok': w == t,
                         'any': w > 0})
            print(f"{p:>7} {n:>7} {lam:>7} {ls:>6.3f} {k1:>4} {eff:>6.3f} "
                  f"{m_thr:>6.1f} {rho:>7.3f} {svpv:>7.3f} "
                  f"{str(w)+'/'+str(t):>6} {bkz_str:>6}")
        return rows

    sweeps = {}
    for eff_t in (0.05, 0.15, 0.25):
        sweeps[eff_t] = sweep_at_eff(curves, eff_t, M_SIGS, SEEDS)


    # ===========================================================================
    # Predictor comparison
    # ===========================================================================

    print("\n" + "-" * 78)
    print("Predictor comparison on the T3 sweeps")
    print("-" * 78)

    def best_threshold(rows, key):
        """Best single-threshold classifier (predict success iff value >= thr)."""
        vals = sorted({r[key] for r in rows})
        best = (0, None)
        for i in range(len(vals) + 1):
            thr = vals[0] - 1 if i == 0 else vals[i - 1]
            acc = sum(1 for r in rows if (r[key] >= thr) == r['ok'])
            if acc > best[0]:
                best = (acc, thr)
        return best

    def rng_str(v):
        return f"[{min(v):.4f}, {max(v):.4f}]" if v else "(none)"

    for eff_t, rows in sweeps.items():
        if not rows:
            continue
        print(f"\n  --- eff = {eff_t} ---")
        for lab in ('ok', 'any'):
            nrows = len(rows)
            nsucc = sum(1 for r in rows if r[lab])
            if nsucc == 0 or nsucc == nrows:
                print(f"  [{lab}] degenerate: {nsucc}/{nrows} positives — "
                      f"no classifier is informative.")
                continue
            baseline = max(nsucc, nrows - nsucc)
            marked = [dict(r, ok=r[lab]) for r in rows]
            acc_ls, thr_ls = best_threshold(marked, 'lam_star')
            acc_rho, thr_rho = best_threshold(marked, 'rho')
            acc_sv, thr_sv = best_threshold(marked, 'svpv')
            print(f"  [{lab}] curves {nrows}, positives {nsucc}, "
                  f"majority baseline {baseline}/{nrows} = {baseline/nrows:.1%}")
            print(f"    lam* >= {thr_ls:.4f} -> {acc_ls}/{nrows} = {acc_ls/nrows:.1%}")
            print(f"    rho  >= {thr_rho:.4f} -> {acc_rho}/{nrows} = {acc_rho/nrows:.1%}")
            print(f"    sv/pv>= {thr_sv:.4f} -> {acc_sv}/{nrows} = {acc_sv/nrows:.1%}")
            for key in ('lam_star', 'rho', 'svpv'):
                s = [r[key] for r in rows if r[lab]]
                f = [r[key] for r in rows if not r[lab]]
                ov = bool(s and f and min(s) <= max(f) and min(f) <= max(s))
                print(f"    {key:<9} success {rng_str(s):<20} "
                      f"failure {rng_str(f):<20} overlap={ov}")


    # ===========================================================================
    # EXPERIMENT T4 — is the "small-lambda failure" really a bias-strength (eff)
    #                 effect?  Sweep K1 on the historical success/failure pair.
    # ===========================================================================

    print("\n" + "-" * 78)
    print("EXP T4: K1 sweep on the historical pair (same n-size, same K2, diff lam*)")
    print("-" * 78)
    print("2026-07-26 concluded lam*=0.07 fails 'structurally'.  Both curves were")
    print("run at K1=8 only.  Here we sweep K1 to separate 'lam effect' from")
    print("'bias-strength (eff) effect'.\n")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    print(f"{'curve':<18} {'lam*':>7} " + " ".join(f"K1={k:<5}" for k in K1_GRID))
    for label, curve, _k1, _m in hist_curves:
        p, b, n, lam, G = curve
        if n < 2000:
            continue          # 8-bit curve has K2=15, grid not comparable
        cells = []
        for k1 in K1_GRID:
            w, t, _ = success_rate(curve, 12, k1, SEEDS)
            cells.append(f"{w}/{t}  ")
        print(f"{label:<18} {lam_star(lam,n):>7.4f} " + " ".join(cells))

    print("\n(eff = K1*K2/n with K2=52, n~2650: "
          + ", ".join(f"K1={k}:{k*52/2647:.3f}" for k in K1_GRID) + ")")

    print("\nT4b: does MORE data rescue the lam*=0.07 curve at K1=8?")
    print("     (info-theoretic threshold at eff=0.157 is m ~ 4.3, so m=12 is")
    print("      already 3x over; if m=32 still fails the wall is not sample size)")
    fail_curve = [c for lbl, c, _, _ in hist_curves if c[2] == 2647][0]
    cells = []
    for m_try in (8, 12, 16, 24, 32):
        w, t, _ = success_rate(fail_curve, m_try, 8, SEEDS)
        cells.append(f"m={m_try}: {w}/{t}")
    print("     " + "   ".join(cells))


    # ===========================================================================
    # EXPERIMENT T5 — where does the (always shorter) non-planted vector live?
    # ===========================================================================

    print("\n" + "-" * 78)
    print("EXP T5: anatomy of the shortest LLL vector (sv/pv ~ 0.35 everywhere)")
    print("-" * 78)
    print("Energy split of the shortest reduced row across the coordinate blocks:")
    print("  k1-block = cols 0..m-1, d = col m, k2-block = cols m+1..2m, "
          "kan = col 2m+1\n")
    print(f"{'curve':<18} {'K1':>4} {'sv/pv':>7} {'k1-blk':>8} {'d':>8} "
          f"{'k2-blk':>8} {'kan':>8} {'kan=S?':>7} {'sv[m]/n':>9}")
    for label, curve, k1, m in hist_curves:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        sigs = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
        M = build_glv_lattice(sigs, n, lam, k1, k2b)
        dim = 2 * m + 2
        A = IntegerMatrix.from_matrix(M)
        LLL.reduction(A)
        rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
        sv = min(rows, key=norm)
        pv = planted_vector(sigs, d0, n, k1, k2b)
        tot = sum(x * x for x in sv) or 1
        e_k1 = sum(sv[i] ** 2 for i in range(m)) / tot
        e_d = sv[m] ** 2 / tot
        e_k2 = sum(sv[m + 1 + i] ** 2 for i in range(m)) / tot
        e_kan = sv[dim - 1] ** 2 / tot
        _, _, _, S_KAN = scales(n, k1, k2b)
        print(f"{label:<18} {k1:>4} {norm(sv)/norm(pv):>7.3f} {e_k1:>8.3f} "
              f"{e_d:>8.3f} {e_k2:>8.3f} {e_kan:>8.3f} "
              f"{str(abs(sv[dim-1]) == S_KAN):>7} {abs(sv[m])/n:>9.4f}")

    print("\nIf kan-energy is ~0 the shortest vector has last coordinate 0, i.e. it")
    print("lies in the non-Kannan sublattice and can never encode d — it is a")
    print("'parasitic' vector that LLL prefers over the planted one.")

    print("\nDone.")


if __name__ == '__main__':
    _driver()
