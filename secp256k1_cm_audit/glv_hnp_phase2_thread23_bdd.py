"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the target is findable.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, commits d525931 / e845207):

  T5 of that run proved that the planted vector is NEVER lambda_1 of the
  Phase-2 Kannan lattice: the shortest vector is always the trivial vector
  n*S_D*e_m (|sv[m]|/n = 1.0000 exactly on every curve tested, success and
  failure alike, sv/pv in [0.34, 0.61]).  Recovery is therefore a BDD/coset
  condition, not an SVP condition.

  T4b of that run showed the K1 wall is NOT relieved by more data: at K1=8 on
  the lam*=0.07 curve, m = 8/12/16/24/32 gave 0,0,1,0,1 out of 5.

  The proposed next step (log line ~6089) was:
    "project the lattice along e_m (quotient out the trivial n*e_m direction)
     and solve BDD in the projection, or replace the Kannan embedding with an
     explicit CVP call (Babai nearest-plane on the reduced basis).
     Falsifier: if sv/pv rises above 1 after the reformulation and the K1 wall
     moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
     reformulation is a real improvement; if the wall stays at K1 ~ 4-6, then
     the wall is information-theoretic and Phase 2 is at its ceiling."

This script executes that falsifier.  Four experiments:

  T23-A  Gaussian-heuristic ratio gamma = ||v_planted|| / GH(L).
         Closed form (derived below, verified numerically): as m grows,
             gamma  ->  sqrt(2*pi*e/3) * sqrt(eff),   eff = K1*K2/n
         which is INDEPENDENT OF m.  If correct this explains T4b directly:
         adding signatures adds two unknown coordinates (k1_i, k2_i) per
         equation, so the planted vector's position relative to the Gaussian
         heuristic never improves.  Also gives a predicted wall at
         eff_crit = 3/(2*pi*e) = 0.1756.

  T23-B  S_D lift.  ||n*S_D*e_m|| = n*S_D scales linearly in S_D but
         ||v_planted||^2 = n^2(2m/3 + 1) + d^2*S_D^2 does not (only 1 of its
         2m+2 coordinates carries S_D).  So S_D > sqrt(m + 3/2) makes the
         planted vector shorter than the trivial vector.  The 2026-07-29 log
         asserted "no choice of S_D removes it -- both vectors scale linearly
         in S_D"; that assertion is tested here.

  T23-C  Babai nearest-plane CVP.  Drop the Kannan row/column entirely and
         solve the BDD instance directly on the (2m+1)-dimensional lattice
         L' = <rows 0..2m>, target t = (-A_i*S_K1)_i (+ zeros).  The error
         vector is then exactly e = (k1_i*S_K1, d*S_D, k2_i*S_K2).  This
         removes the trivial-vector obstruction by construction: nearest-plane
         does not care what lambda_1 is.

  T23-D  Exact per-instance Babai predicate.  Writing e = sum_j mu_j b*_j in
         the Gram-Schmidt basis of the reduced L', nearest-plane recovers the
         target IFF max_j |mu_j| < 1/2.  This is an exact, computable,
         per-instance success criterion -- the first one this project has had
         (the 2026-06-21..06-29 runs burned six curve-level invariants trying
         to predict this from curve data alone, which T5 explained is
         impossible).

Run: python3 glv_hnp_phase2_thread23_bdd.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (a=0 short Weierstrass) -- same as glv_hnp_phase2_20bit.py
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
# j=0 CM helpers
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
    """Both roots of x^2+x+1 = 0 mod n, or None."""
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
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=1, nbins=10):
    """j=0 GLV curves bucketed by lam* over [0, 0.5]."""
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1: continue
                    if not sympy.isprime(n_cand): continue
                    roots = glv_roots(n_cand)
                    if roots is None: continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_cand) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin: continue
                    cur = build_curve(p, n_cand)
                    if cur is None: continue
                    bins[idx].append((p, cur[1], n_cand, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()): break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

# ---------------------------------------------------------------------------
# Signatures and lattices
# ---------------------------------------------------------------------------

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

def scales(n, k1_bound, k2_bound, s_d=1):
    return (n // k1_bound, s_d, max(1, n // k2_bound), n)

def build_glv_lattice(sigs, n, lam, S_K1, S_D, S_K2, S_KANNAN):
    """Kannan-embedded Phase-2 lattice.  Identical to
    glv_hnp_phase2_20bit.py:262 except that S_D is a free parameter."""
    m = len(sigs)
    dim = 2 * m + 2
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

def build_bdd_lattice(sigs, n, lam, S_K1, S_D, S_K2):
    """The SAME lattice with the Kannan row and column removed: dim 2m+1.
    Rows 0..m-1: n*S_K1*e_i;  row m: (B_i*S_K1 | S_D | 0);
    rows m+1..2m: (-lam*S_K1*e_i | 0 | S_K2*e_i).
    Target for CVP is t = (-A_i*S_K1 | 0 | 0); the error vector is then
    exactly e = (k1_i*S_K1 | d*S_D | k2_i*S_K2)."""
    m = len(sigs)
    dim = 2 * m + 1
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    target = [0] * dim
    for i in range(m):
        target[i] = -sigs[i]['A'] * S_K1
    return M, target

# ---------------------------------------------------------------------------
# Recovery: Kannan embedding
# ---------------------------------------------------------------------------

def recover_kannan(reduced, m, n, S_D, S_KANNAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        v = sign * row[m]
        if v % S_D != 0: continue
        d_cand = (v // S_D) % n
        if d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Gram-Schmidt + Babai nearest-plane (float; entries here are < 2^40 so f64
# has ~13 significant decimal digits of headroom on the mu test)
# ---------------------------------------------------------------------------

def gram_schmidt(B):
    dim = len(B)
    bs = [[float(x) for x in row] for row in B]
    mu = [[0.0] * dim for _ in range(dim)]
    norms = [0.0] * dim
    for i in range(dim):
        for j in range(i):
            if norms[j] == 0.0: continue
            c = sum(bs[i][k] * bs[j][k] for k in range(dim)) / norms[j]
            mu[i][j] = c
            for k in range(dim):
                bs[i][k] -= c * bs[j][k]
        norms[i] = sum(x * x for x in bs[i])
    return bs, norms

def babai_nearest_plane(B, bstar, norms, target):
    """Return the lattice vector v closest to `target` under nearest-plane."""
    dim = len(B)
    t = [float(x) for x in target]
    coeffs = [0] * dim
    for i in range(dim - 1, -1, -1):
        if norms[i] == 0.0:
            c = 0
        else:
            c = round(sum(t[k] * bstar[i][k] for k in range(dim)) / norms[i])
        coeffs[i] = c
        if c:
            for k in range(dim):
                t[k] -= c * B[i][k]
    v = [0] * dim
    for i in range(dim):
        if coeffs[i]:
            for k in range(dim):
                v[k] += coeffs[i] * B[i][k]
    return v

def babai_mu_profile(bstar, norms, err):
    """mu_j = <e, b*_j> / ||b*_j||^2.  Nearest-plane succeeds IFF all |mu_j|<1/2."""
    dim = len(bstar)
    out = []
    for j in range(dim):
        if norms[j] == 0.0:
            out.append(float('inf')); continue
        out.append(sum(err[k] * bstar[j][k] for k in range(dim)) / norms[j])
    return out

# ---------------------------------------------------------------------------
# Predicted quantities
# ---------------------------------------------------------------------------

def planted_norm(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN):
    """E[||v_planted||] with k1_i ~ U[0,K1), k2_i ~ U[0,K2), d ~ U[0,n)."""
    return math.sqrt(m * (k1 * S_K1) ** 2 / 3.0
                     + (n * S_D) ** 2 / 3.0
                     + m * (k2 * S_K2) ** 2 / 3.0
                     + S_KAN ** 2)

def gaussian_heuristic(m, n, S_K1, S_D, S_K2, S_KAN):
    dim = 2 * m + 2
    log_det = (m * math.log(n * S_K1) + math.log(S_D)
               + m * math.log(S_K2) + math.log(S_KAN))
    return math.exp(0.5 * math.log(dim / (2 * math.pi * math.e))
                    + log_det / dim)

def gamma_ratio(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN):
    return (planted_norm(m, n, k1, k2, S_K1, S_D, S_K2, S_KAN)
            / gaussian_heuristic(m, n, S_K1, S_D, S_K2, S_KAN))

GAMMA_ASYMPTOTIC = math.sqrt(2 * math.pi * math.e / 3.0)   # 2.3862...

# ---------------------------------------------------------------------------
# One instance, all four methods
# ---------------------------------------------------------------------------

def run_instance(curve, m, k1_bound, seed, s_d=1, bkz_beta=0, want_mu=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, s_d)
    out = {'n': n, 'lam': lam, 'k1': k1_bound, 'k2': k2_bound, 'm': m,
           'eff': k1_bound * k2_bound / n, 's_d': s_d,
           'gamma': gamma_ratio(m, n, k1_bound, k2_bound, S_K1, S_D, S_K2, S_KAN)}

    # --- method 1/2: Kannan embedding + LLL (or BKZ) ---
    M = build_glv_lattice(sigs, n, lam, S_K1, S_D, S_K2, S_KAN)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if bkz_beta:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    out['kannan'] = recover_kannan(reduced, m, n, S_D, S_KAN, d_secret)
    # sv / pv  (the T5 diagnostic, recomputed at this S_D)
    svn = min(math.sqrt(sum(float(x) * x for x in row)) for row in reduced
              if any(row))
    pvn = math.sqrt(sum((sigs[i]['k1'] * S_K1) ** 2 for i in range(m))
                    + (d_secret * S_D) ** 2
                    + sum((sigs[i]['k2'] * S_K2) ** 2 for i in range(m))
                    + S_KAN ** 2)
    out['sv_pv'] = svn / pvn

    # --- method 3: Babai nearest-plane CVP on the Kannan-free lattice ---
    B, target = build_bdd_lattice(sigs, n, lam, S_K1, S_D, S_K2)
    dimb = 2 * m + 1
    Ab = IntegerMatrix.from_matrix(B)
    if bkz_beta:
        BKZ.reduction(Ab, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(Ab)
    Br = [[Ab[i][j] for j in range(dimb)] for i in range(dimb)]
    bstar, norms = gram_schmidt(Br)
    v = babai_nearest_plane(Br, bstar, norms, target)
    d_cand = (v[m] // S_D) % n if v[m] % S_D == 0 else None
    out['babai'] = (d_cand == d_secret)

    # --- method 4: exact Babai predicate on the true error vector ---
    if want_mu:
        # true closest vector v_true = target + e, with e the small part
        e = [0] * dimb
        for i in range(m):
            e[i] = sigs[i]['k1'] * S_K1
            e[m + 1 + i] = sigs[i]['k2'] * S_K2
        e[m] = d_secret * S_D
        mus = babai_mu_profile(bstar, norms, e)
        out['max_mu'] = max(abs(x) for x in mus)
        out['n_bad_mu'] = sum(1 for x in mus if abs(x) >= 0.5)
        out['err_norm'] = math.sqrt(sum(float(x) * x for x in e))
        out['last_gs'] = math.sqrt(norms[-1]) if norms[-1] > 0 else 0.0
    return out


def rate(curve, m, k1_bound, seeds, s_d=1, bkz_beta=0, want_mu=False):
    acc = {'kannan': 0, 'babai': 0, 'tot': 0,
           'sv_pv': [], 'gamma': None, 'max_mu': [], 'n_bad_mu': []}
    for s in seeds:
        r = run_instance(curve, m, k1_bound, s, s_d, bkz_beta, want_mu)
        if r is None: continue
        acc['tot'] += 1
        acc['kannan'] += bool(r['kannan'])
        acc['babai'] += bool(r['babai'])
        acc['sv_pv'].append(r['sv_pv'])
        acc['gamma'] = r['gamma']
        if want_mu:
            acc['max_mu'].append(r['max_mu'])
            acc['n_bad_mu'].append(r['n_bad_mu'])
    return acc


def mean(xs):
    return sum(xs) / len(xs) if xs else float('nan')


# ===========================================================================
def main():
    print("=" * 88)
    print("Thread 23 — reformulating the GLV-HNP Phase-2 lattice (BDD, not SVP)")
    print("=" * 88)

    SEEDS = [42, 1234, 9999, 555, 31337]

    # Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
    HIST = [
        # label,              p,    b, n,    lam,  lam*
        ("12-bit/2557",       2557, 2, 2659, 1755, 0.340),
        ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  0.070),
    ]
    hist = []
    for label, p, b, n, lam, ls in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), ls))
        print(f"  {label}: p={p} n={n} lam={lam} lam*={lam_star(lam,n):.3f} G={G}")


    # ===========================================================================
    # T23-A  gamma is m-independent  ->  explains T4b
    # ===========================================================================
    print("\n" + "-" * 88)
    print("T23-A: Gaussian-heuristic ratio gamma = ||v_planted|| / GH(L)")
    print("-" * 88)
    print("Claim: gamma -> sqrt(2*pi*e/3)*sqrt(eff) = %.4f*sqrt(eff), independent of m."
          % GAMMA_ASYMPTOTIC)
    print("If true, adding signatures cannot help: each new equation brings two new")
    print("unknown coordinates, so the planted vector never improves w.r.t. GH(L).\n")

    print(f"{'curve':>16} {'K1':>4} {'eff':>7} | " +
          " ".join(f"{'m='+str(mm):>7}" for mm in (4, 8, 12, 16, 24, 32, 64, 128)) +
          f" {'asympt':>8}")
    for label, curve, ls in hist:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in (4, 8, 16):
            eff = k1 * k2 / n
            S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
            row = [gamma_ratio(mm, n, k1, k2, S_K1, S_D, S_K2, S_KAN)
                   for mm in (4, 8, 12, 16, 24, 32, 64, 128)]
            print(f"{label:>16} {k1:>4} {eff:>7.4f} | " +
                  " ".join(f"{g:>7.4f}" for g in row) +
                  f" {GAMMA_ASYMPTOTIC*math.sqrt(eff):>8.4f}")

    print("\nPredicted unique-SVP wall gamma = 1  <=>  eff_crit = 3/(2*pi*e) = %.4f"
          % (3.0 / (2 * math.pi * math.e)))


    # ===========================================================================
    # T23-B  S_D lift: can the trivial vector n*S_D*e_m be demoted?
    # ===========================================================================
    print("\n" + "-" * 88)
    print("T23-B: S_D lift.  ||trivial|| = n*S_D (linear);")
    print("       ||v_planted||^2 = n^2(2m/3+1) + d^2*S_D^2  (only 1 of 2m+2 coords)")
    print("       => planted is shorter iff S_D > sqrt(m + 3/2).")
    print("-" * 88)

    M_SIGS = 12
    print(f"m = {M_SIGS}  =>  predicted crossover at S_D > sqrt({M_SIGS}+1.5) = "
          f"{math.sqrt(M_SIGS+1.5):.2f}\n")
    print(f"{'curve':>16} {'K1':>4} {'S_D':>5} {'gamma':>7} {'sv/pv':>7} "
          f"{'Kannan':>7} {'Babai':>7}")
    sd_rows = []
    for label, curve, ls in hist:
        for k1 in (4, 8):
            for s_d in (1, 2, 4, 8, 16, 64):
                a = rate(curve, M_SIGS, k1, SEEDS, s_d=s_d)
                print(f"{label:>16} {k1:>4} {s_d:>5} {a['gamma']:>7.3f} "
                      f"{mean(a['sv_pv']):>7.3f} "
                      f"{str(a['kannan'])+'/'+str(a['tot']):>7} "
                      f"{str(a['babai'])+'/'+str(a['tot']):>7}")
                sd_rows.append((label, k1, s_d, a))
        print()


    # ===========================================================================
    # T23-C  the falsifier: does the K1 wall move?
    # ===========================================================================
    print("-" * 88)
    print("T23-C: THE FALSIFIER — K1 wall, Kannan-LLL vs Babai-CVP vs BKZ(20)-Babai")
    print("       2026-07-29 T4 baseline: 2557 wall at K1~12-16, 2677 wall at K1~4-6")
    print("-" * 88)

    K1_GRID = (2, 3, 4, 6, 8, 12, 16, 24)
    print(f"{'curve':>16} {'K1':>4} {'eff':>7} {'gamma':>7} {'sv/pv':>7} "
          f"{'Kannan':>7} {'Babai':>7} {'BKZ20-B':>8} {'max|mu|':>8} {'bad_mu':>7}")
    wall_rows = []
    for label, curve, ls in hist:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        for k1 in K1_GRID:
            a = rate(curve, M_SIGS, k1, SEEDS, s_d=1, want_mu=True)
            ab = rate(curve, M_SIGS, k1, SEEDS, s_d=1, bkz_beta=20)
            print(f"{label:>16} {k1:>4} {k1*k2/n:>7.4f} {a['gamma']:>7.3f} "
                  f"{mean(a['sv_pv']):>7.3f} "
                  f"{str(a['kannan'])+'/'+str(a['tot']):>7} "
                  f"{str(a['babai'])+'/'+str(a['tot']):>7} "
                  f"{str(ab['babai'])+'/'+str(ab['tot']):>8} "
                  f"{mean(a['max_mu']):>8.3f} {mean(a['n_bad_mu']):>7.1f}")
            wall_rows.append({'label': label, 'k1': k1, 'eff': k1 * k2 / n,
                              'gamma': a['gamma'],
                              'kannan': a['kannan'], 'babai': a['babai'],
                              'bkz_babai': ab['babai'], 'tot': a['tot'],
                              'max_mu': mean(a['max_mu'])})
        print()


    # ===========================================================================
    # T23-D  exact Babai predicate vs observed outcome
    # ===========================================================================
    print("-" * 88)
    print("T23-D: exact per-instance predicate.  Nearest-plane recovers d IFF")
    print("       max_j |mu_j| < 1/2 where e = sum_j mu_j b*_j.  Check agreement.")
    print("-" * 88)

    agree = disagree = 0
    print(f"{'curve':>16} {'K1':>4} {'seed':>7} {'max|mu|':>9} {'pred':>5} "
          f"{'Babai':>6} {'Kannan':>7}")
    for label, curve, ls in hist:
        for k1 in (2, 4, 6, 8, 12, 16):
            for s in SEEDS[:3]:
                r = run_instance(curve, M_SIGS, k1, s, want_mu=True)
                if r is None: continue
                pred = r['max_mu'] < 0.5
                ok = (pred == r['babai'])
                agree += ok; disagree += (not ok)
                print(f"{label:>16} {k1:>4} {s:>7} {r['max_mu']:>9.4f} "
                      f"{str(pred):>5} {str(r['babai']):>6} {str(r['kannan']):>7}"
                      + ("" if ok else "   <-- MISMATCH"))
    print(f"\npredicate agreement: {agree}/{agree+disagree}")


    # ===========================================================================
    # T23-E  cross-curve confirmation on a fresh 17-bit sweep
    # ===========================================================================
    print("\n" + "-" * 88)
    print("T23-E: fresh 17-bit sweep — Kannan vs Babai at three eff levels")
    print("-" * 88)
    LO, HI = 2 ** 16, 2 ** 17
    print(f"searching j=0 GLV curves with p in [{LO},{HI}) ...")
    curves17 = search_curves(LO, HI, per_bin=1, nbins=10)
    print(f"found {len(curves17)} curves\n")

    for eff_t in (0.05, 0.15, 0.25):
        print(f"=== eff target {eff_t} (m={M_SIGS}, {len(SEEDS)} seeds) ===")
        print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>7} {'gamma':>7} "
              f"{'Kannan':>7} {'Babai':>7} {'max|mu|':>8}")
        tk = tb = tt = 0
        for (p, b, n, lam, G) in curves17:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            a = rate((p, b, n, lam, G), M_SIGS, k1, SEEDS, want_mu=True)
            tk += a['kannan']; tb += a['babai']; tt += a['tot']
            print(f"{p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {k1*k2/n:>7.4f} "
                  f"{a['gamma']:>7.3f} "
                  f"{str(a['kannan'])+'/'+str(a['tot']):>7} "
                  f"{str(a['babai'])+'/'+str(a['tot']):>7} "
                  f"{mean(a['max_mu']):>8.3f}")
        print(f"  TOTAL  Kannan {tk}/{tt}   Babai {tb}/{tt}\n")

    print("=" * 88)
    print("done")
    print("=" * 88)



if __name__ == "__main__":
    main()
