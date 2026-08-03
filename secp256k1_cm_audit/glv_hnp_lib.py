"""
GLV-HNP shared library — side-effect-free helpers for the Phase-2 lattice.

Why this file exists
--------------------
`glv_hnp_phase2_lambda_threshold.py` (Thread 20a, commit d525931) holds the
canonical Phase-2 helpers, but it executes experiments T1-T5 at module scope,
so importing it runs a multi-minute experiment suite.  The two Thread-20c/20d
scripts committed on 2026-07-29 (`glv_hnp_phase2_nuhat_control.py`,
`glv_hnp_nuhat_vs_c1c2.py`) import it and then reference four names that the
module does not define:

    rival_sublattice_nu, glv_eigenvalues, mu_of, identify_twist

so both are un-runnable as committed (AttributeError at import).  This module
carries the helpers verbatim, defines the four missing names, and adds the
Thread-23 lattice-geometry functions.  No code runs at import time.

Verbatim provenance: EC arithmetic / CM helpers / lattice builder are copied
unchanged from glv_hnp_phase2_lambda_threshold.py:87-350 so that every number
computed here is comparable to the 2026-07-26 .. 2026-07-29 log entries.
"""

import math
import random

import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0)
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

# --- names the 2026-07-29 Thread-20c/20d scripts expect -------------------

glv_eigenvalues = glv_roots      # alias: both roots of x^2+x+1 mod n
mu_of = lam_star                 # alias: mu = lam* = min(lam,n-lam)/n

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

def rival_sublattice_nu(n, lam, k1_bound, k2_bound):
    """nu_hat = lambda_1(L2) / sqrt(det L2)  for the non-planted 2D block L2.

    det(L2) = n*S_K1*S_K2 is independent of lam, so nu_hat isolates exactly the
    lam-dependence of the Phase-2 geometry.  nu_hat in (0, (4/3)^(1/4)].

    Takes the *bounds* K1, K2 (not the scales) and derives S_K1/S_K2 via
    scales(), matching the call sites in glv_hnp_phase2_nuhat_control.py:117,
    glv_hnp_nuhat_vs_c1c2.py:145 and glv_hnp_phase2_mu_response.py:196, all of
    which unpack the third component as nu_hat.

    Returns (w, mu, nu_hat) with w the exact shortest vector and mu = |w|.
    """
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    mu, w = lambda_block_mu(n, lam, S_K1, S_K2)
    det2 = n * S_K1 * S_K2
    return w, mu, mu / math.sqrt(det2)

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
    """The (2m+2)-dimensional Phase-2 lattice.  Column layout:
         0..m-1     k1 columns   (scale S_K1)
         m          d  column    (scale S_D = 1)
         m+1..2m    k2 columns   (scale S_K2)
         2m+1       Kannan column (scale S_KANNAN = n)
    """
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
    return math.sqrt(sum(float(x) * float(x) for x in v))

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

def run_experiment(curve, m, d_secret, k1_bound, seed=42, lam_override=None,
                   use_bkz=False, bkz_beta=20):
    """Returns (recovered: bool, planted_norm, shortest_reduced_norm)."""
    p, b, n, lam, G = curve
    lam_lat = lam if lam_override is None else lam_override
    k2_bound = math.isqrt(n) + 1
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

def run_experiment_flat(p, n, lam, G, m, d_secret, k1_bound, k2_bound, seed):
    """Flat-signature variant returning just the success bool.

    The three Thread-20c/20d scripts committed on 2026-07-29 call
    `run_experiment(p, n, lam, G, m, d, k1, k2, seed)` and use the result as an
    int (`wins += ...`); no such 9-argument entry point ever existed in the
    repo, which is the second reason those scripts could not run.  Unlike
    run_experiment() this honours the caller's k2_bound instead of recomputing
    isqrt(n)+1.
    """
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    _, _, _, S_KANNAN = scales(n, k1_bound, k2_bound)
    return recover_d(reduced, m, n, S_KANNAN, d_secret) is not None

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

def identify_twist(p, n, seed=12345):
    """The sextic twist parameter b with #E_b(F_p) = n, or None."""
    cur = build_curve(p, n, seed=seed)
    return None if cur is None else cur[1]

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

# ===========================================================================
# Thread 23 additions — exact lattice geometry
# ===========================================================================

def log_det_phase2(n, m, k1_bound, k2_bound):
    """log|det L| for the (2m+2)-dim Phase-2 lattice.

    The basis is block-triangular in the column order
    (k1 cols | d col | k2 cols | Kannan col), so
        det = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN
    exactly (no approximation).  Returned as a natural log because det itself
    overflows float64 for m >= 6.
    """
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    return (m * math.log(n * S_K1) + math.log(S_D)
            + m * math.log(S_K2) + math.log(S_KAN))

def gaussian_heuristic(log_det, dim):
    """GH(L) = sqrt(dim / (2*pi*e)) * det^(1/dim)."""
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def usvp_gap(n, m, k1_bound, k2_bound, planted_nrm):
    """nu_gh = ||v_planted|| / GH(L).

    nu_gh < 1  =>  v_planted is plausibly the unique shortest vector.
    nu_gh > 1  =>  the Gaussian heuristic predicts ~(nu_gh)^dim lattice vectors
                   shorter than the target, so uSVP is ill-posed no matter how
                   strong the reduction.
    """
    dim = 2 * m + 2
    ld = log_det_phase2(n, m, k1_bound, k2_bound)
    return planted_nrm / gaussian_heuristic(ld, dim)

def usvp_gap_closed_form(n, m, k1_bound, k2_bound):
    """Closed form of nu_gh under E[||v||^2] = V^2*(2m+4)/3 and det as above.

        nu_gh  ~=  sqrt(2*pi*e/3) * (n * eff^m)^(1/(2m+2)),    eff = K1*K2/n

    Derivation: with V = n (the common column scale), S_K1 = V/K1, S_D = V/n,
    S_K2 = V/K2, S_KANNAN = V, det = n^(m-1) V^(2m+2) / (K1 K2)^m, so
    det^(1/dim) = V * (1/(n*eff^m))^(1/dim), while ||v|| = V*sqrt((2m+4)/3).
    """
    eff = k1_bound * k2_bound / n
    dim = 2 * m + 2
    return math.sqrt(2 * math.pi * math.e / 3.0) * math.exp(
        (math.log(n) + m * math.log(eff)) / dim)

def gram_schmidt(rows):
    """Plain float64 Gram-Schmidt.  Returns (bstar, r) with r[i] = ||b*_i||.
    Used only on LLL-reduced bases (well conditioned, dim <= 40)."""
    bstar, r = [], []
    for row in rows:
        v = [float(x) for x in row]
        for u, nu in zip(bstar, r):
            if nu == 0.0:
                continue
            c = sum(a * b for a, b in zip(v, u)) / (nu * nu)
            v = [a - c * b for a, b in zip(v, u)]
        nv = math.sqrt(sum(a * a for a in v))
        bstar.append(v)
        r.append(nv)
    return bstar, r

def projected_target_profile(reduced_rows, v):
    """||pi_i(v)|| for i = 0..dim-1, where pi_i projects orthogonally to
    span(b_0..b_{i-1}) of the (reduced) basis.  pi_0(v) = v."""
    bstar, r = gram_schmidt(reduced_rows)
    w = [float(x) for x in v]
    out = []
    for i in range(len(reduced_rows)):
        out.append(math.sqrt(max(0.0, sum(a * a for a in w))))
        if r[i] == 0.0:
            continue
        c = sum(a * b for a, b in zip(w, bstar[i])) / (r[i] * r[i])
        w = [a - c * b for a, b in zip(w, bstar[i])]
    return out, r

def adps_margin(reduced_rows, v):
    """max_i ||b*_i|| / ||pi_i(v)||  (the Alkim-Ducas-Poeppelmann-Schwabe
    'target is visible at index i' condition, in ratio form).

    margin >= 1 at some index i  =>  the projection of the planted vector is
    no longer than the local Gram-Schmidt norm there, so a reduction that is
    strong enough at index i can lift it out.  margin < 1 everywhere => the
    planted vector is invisible to *any* index of this basis profile.
    Returns (margin, argmax_index).
    """
    pi, r = projected_target_profile(reduced_rows, v)
    # v is a lattice vector, so pi_i(v) collapses to 0 past the last nonzero
    # coefficient of v in this basis.  Those indices carry no information and
    # would give a spurious +inf ratio, so require pi_i(v) to be a nonneglible
    # fraction of ||v|| before scoring index i.
    floor = 1e-9 * (pi[0] if pi else 0.0)
    best, besti = 0.0, -1
    for i in range(len(r)):
        if pi[i] <= floor:
            break
        val = r[i] / pi[i]
        if val > best:
            best, besti = val, i
    return best, besti


def nuisance_sublattice(n, lam, m, S_K1, S_D, S_K2):
    """The rank-(m+1) 'nuisance' sublattice T of information-free short vectors:

        t_0 = n*S_D*e_m                              (Thread 20 T5 trivial vector)
        t_i = w[0]*e_i + w[1]*e_{m+1+i},  i < m      (Thread 20c rival blocks)

    with w = lambda_1 of the 2D block <(n*S_K1,0), (-lam*S_K1,S_K2)>.  The m+1
    generators have pairwise disjoint support, hence are pairwise ORTHOGONAL,
    so det(T) = (n*S_D) * mu^m exactly.  Every t in T has Kannan coordinate 0,
    so no vector of T can ever encode d.

    Returns (w, mu, log_det_T).
    """
    mu, w = lambda_block_mu(n, lam, S_K1, S_K2)
    return w, mu, math.log(n * S_D) + m * math.log(mu)


def usvp_gap_projected(n, lam, m, k1_bound, k2_bound, sigs, d_secret):
    """nu_proj: the uSVP gap in the projection orthogonal to the nuisance
    sublattice T (Thread 23).

        nu_proj = ||pi_T^perp(v_planted)|| / GH(pi_T^perp(L)),
        dim(pi_T^perp(L)) = (2m+2) - (m+1) = m+1,
        log det(pi_T^perp(L)) = log det(L) - log det(T).

    This unifies the two effects Thread 20 found separately.  A *small* nu_hat
    (skew rival block, small mu) shrinks det(T) and therefore INFLATES
    det(pi(L)) and GH, lowering nu_proj -- which is exactly the counter-
    intuitive sign observed in Thread 20c: short rivals make the attack easier
    because they get consumed by the first basis indices, leaving the planted
    vector to compete against a relatively bulkier projected lattice.
    """
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    w, mu, log_det_T = nuisance_sublattice(n, lam, m, S_K1, S_D, S_K2)
    v = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)

    nrm2 = sum(float(x) * float(x) for x in v)
    nrm2 -= float(v[m]) ** 2                                # component on t_0
    mu2 = float(mu) * float(mu)
    for i in range(m):                                      # components on t_i
        ip = float(v[i]) * w[0] + float(v[m + 1 + i]) * w[1]
        nrm2 -= ip * ip / mu2
    proj_norm = math.sqrt(max(0.0, nrm2))

    dim = m + 1
    log_det = log_det_phase2(n, m, k1_bound, k2_bound) - log_det_T
    return proj_norm / gaussian_heuristic(log_det, dim), proj_norm, mu

# --- Thread 23 reformulation: drop the (information-free) d column ---------

def build_glv_lattice_nod(sigs, n, lam, k1_bound, k2_bound):
    """Generating set for the rank-(2m+1) lattice pi(L), where pi deletes the
    d column (column m) of the Phase-2 lattice.

    Deleting the d column kills the trivial vector n*S_D*e_m (it maps to 0),
    which Thread 20 T5 showed is always lambda_1 of the (2m+2)-dim lattice.
    d is recovered afterwards from d = B_0^{-1}(k1_0 - A_0 + lam*k2_0) mod n,
    so no information is lost: d was only ever defined mod n.

    Returned as (2m+2) generators in dimension 2m+1 -- one more generator than
    the rank, because ker(pi|_L) = <n*S_D*e_m> is one dimensional.  LLL is run
    on the dependent generating set and the resulting zero row is dropped.

    Column layout: 0..m-1 k1 | m..2m-1 k2 | 2m Kannan.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(2 * m + 2)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1          # the d row, now with no d column
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_vector_nod(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_nod(M_reduced, sigs, n, k1_bound, k2_bound, d_secret):
    """Read (k1_0, k2_0) off any row with Kannan coordinate +-S_KANNAN and
    solve d = B_0^{-1}(k1_0 - A_0 + lam*k2_0) mod n.  lam is not needed
    explicitly: k_full_0 = k1_0 + lam*k2_0 is what the row encodes jointly, so
    we reconstruct it from the k1/k2 entries and the known scales."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        sgn = 1 if last > 0 else -1
        # k1_i and k2_i must be exact integer multiples of their scales
        ks = []
        ok = True
        for i in range(m):
            a, c = sgn * row[i], sgn * row[m + i]
            if a % S_K1 or c % S_K2:
                ok = False
                break
            ks.append((a // S_K1, c // S_K2))
        if not ok:
            continue
        k1_0, k2_0 = ks[0]
        B0 = sigs[0]['B']
        if B0 % n == 0:
            continue
        lam_roots = glv_roots(n)
        if lam_roots is None:
            continue
        for lam_c in lam_roots:
            k_full = (k1_0 + lam_c * k2_0) % n
            d_cand = (k_full - sigs[0]['A']) * modinv(B0, n) % n
            if d_cand == d_secret:
                return d_cand
    return None

def run_experiment_nod(curve, m, d_secret, k1_bound, seed=42,
                       use_bkz=False, bkz_beta=20):
    """Same as run_experiment but on the d-column-free (2m+1)-dim lattice."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice_nod(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)                       # handles the dependent generator
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    reduced = [r for r in rows if any(r)]  # drop the zero row from the kernel
    if use_bkz and len(reduced) == dim:
        A2 = IntegerMatrix.from_matrix(reduced)
        BKZ.reduction(A2, BKZ.Param(bkz_beta))
        reduced = [[A2[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_d_nod(reduced, sigs, n, k1_bound, k2_bound, d_secret) is not None
    pn = norm(planted_vector_nod(sigs, d_secret, n, k1_bound, k2_bound))
    sn = min(norm(r) for r in reduced)
    return (ok, pn, sn)

# ---------------------------------------------------------------------------
# Small statistics helpers
# ---------------------------------------------------------------------------

def auc(scores, labels):
    """Area under the ROC curve for `scores` predicting boolean `labels`.
    Rank-based (Mann-Whitney), handles ties by mid-rank."""
    pos = [s for s, y in zip(scores, labels) if y]
    neg = [s for s, y in zip(scores, labels) if not y]
    if not pos or not neg:
        return float('nan')
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        mid = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = mid
        i = j + 1
    rsum = sum(r for r, y in zip(ranks, labels) if y)
    return (rsum - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            mid = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = mid
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')
