"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector can be
lambda_1, and test whether the K1 wall (T4, 2026-07-29) moves.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 entry, EXP T5):
  The Kannan-embedded Phase-2 lattice (dim 2m+2) always contains the trivial
  vector n*S_D*e_m (= n*row_m - sum_i B_i*row_i), of norm n*S_D, which is
  shorter than ||v_planted|| ~ n*sqrt(2m/3 + 4/3) for every m >= 1.  So the
  planted vector is NEVER lambda_1 and recovery is a BDD/coset condition, not
  an SVP condition.

This script builds three formulations and compares them on the same curves,
signatures and seeds:

  KAN   dim 2m+2  the historical Kannan embedding (verbatim from
                  glv_hnp_phase2_lambda_threshold.py:254)
  PROJ  dim 2m+1  KAN projected along e_m (the d-column is deleted).  Legal
                  because n*e_m is in KAN, so the coordinate deletion is an
                  orthogonal projection of the lattice.  d is recovered from
                  (k1_i,k2_i) afterwards via d = (k1_i + lam*k2_i - A_i)/B_i.
  BDD   dim 2m    PROJ with the Kannan column and the A-row removed, solved as
                  an explicit closest-vector problem against t = (A_i*S_K1, 0),
                  by Babai nearest-plane (exact rational GS) and by exact CVP.

Falsifier as stated on 2026-07-29:
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
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
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3."""
    if p % 3 != 1: return None
    lim = int(2 * math.isqrt(p)) + 2
    for a in range(0, lim):
        for b in range(0, lim):
            if a * a - a * b + b * b == p:
                return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None: return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0]*w[0] + w[1]*w[1]
    def dot(w, z): return w[0]*z[0] + w[1]*z[1]
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0: break
        q = (2*num + den)//(2*den) if num >= 0 else -((-2*num + den)//(2*den))
        v = (v[0] - q*u[0], v[1] - q*u[1])
        if nrm2(v) >= nrm2(u): break
        u, v = v, u
    return u

def nu_hat(n, lam, S_K1, S_K2):
    """lambda_1(L2)/sqrt(det L2) for L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)>
    (the 2026-07-29 separator).  det L2 = n*S_K1*S_K2, independent of lam."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    mu = math.sqrt(w[0]*w[0] + w[1]*w[1])
    return mu / math.sqrt(n * S_K1 * S_K2)

# ---------------------------------------------------------------------------
# Signatures + scales (verbatim)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 400000:
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
# Formulation KAN (dim 2m+2) -- historical
# ---------------------------------------------------------------------------

def build_kan(sigs, n, lam, k1_bound, k2_bound):
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

def planted_kan(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_kan(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Formulation PROJ (dim 2m+1) -- KAN projected along e_m
#   columns: [0,m)   k1 block (scale S_K1)
#            [m,2m)  k2 block (scale S_K2)
#            2m      Kannan  (scale S_KAN)
#   generators: 2m+2 of them for a rank-(2m+1) lattice (n*(B-row) is in the
#   span of the n*e_i rows), so LLL is fed a rank-deficient generating set and
#   returns one zero row.
# ---------------------------------------------------------------------------

def build_proj(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # n*S_K1*e_i
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                            # the d-elimination row
    for i in range(m): r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                       # lam rows
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * dim                            # Kannan / A row
    for i in range(m): r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KAN
    rows.append(r)
    return rows

def planted_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def d_from_kblocks(v, sigs, n, lam, S_K1, S_K2, m):
    """v = (k1_i*S_K1 | k2_i*S_K2 | ...).  Recover d, requiring consistency
    across every signature."""
    k1s, k2s = [], []
    for i in range(m):
        a, b = v[i], v[m + i]
        if a % S_K1 or b % S_K2:
            return None
        k1s.append(a // S_K1); k2s.append(b // S_K2)
    d = None
    for i in range(m):
        B = sigs[i]['B'] % n
        if B == 0: continue
        cand = (k1s[i] + lam * k2s[i] - sigs[i]['A']) * modinv(B, n) % n
        if d is None:
            d = cand
        elif d != cand:
            return None
    return d

def recover_proj(rows, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        v = [x if last > 0 else -x for x in row]
        d = d_from_kblocks(v, sigs, n, lam, S_K1, S_K2, m)
        if d is not None and d == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Formulation BDD (dim 2m) -- explicit CVP
# ---------------------------------------------------------------------------

def build_bdd(sigs, n, lam, k1_bound, k2_bound):
    """Lattice L0 (generating set) and target t with v = t + w, w in L0."""
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m): r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    t = [0] * dim
    for i in range(m): t[i] = sigs[i]['A'] * S_K1
    return rows, t

def lll_strip(rows, dim):
    """LLL-reduce a (possibly rank-deficient) generating set; drop zero rows."""
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    out = []
    for i in range(A.nrows):
        row = [A[i][j] for j in range(dim)]
        if any(row):
            out.append(row)
    return out

def gso_rational(basis):
    """Exact Gram-Schmidt over Q.  Returns (bstar, mu) with bstar rational."""
    bstar, mu = [], []
    for i, b in enumerate(basis):
        v = [Fraction(x) for x in b]
        mu_i = []
        for j in range(i):
            bj = bstar[j]
            den = sum(x * x for x in bj)
            num = sum(Fraction(b[k]) * bj[k] for k in range(len(b)))
            m_ij = num / den if den != 0 else Fraction(0)
            mu_i.append(m_ij)
            v = [v[k] - m_ij * bj[k] for k in range(len(v))]
        bstar.append(v); mu.append(mu_i)
    return bstar, mu

def babai_nearest_plane(basis, target):
    """Closest-ish lattice vector to `target` (exact rational GS)."""
    bstar, _ = gso_rational(basis)
    b = [Fraction(x) for x in target]
    coeffs = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        den = sum(x * x for x in bstar[i])
        if den == 0: continue
        num = sum(b[k] * bstar[i][k] for k in range(len(b)))
        c = num / den
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = ci
        b = [b[k] - ci * Fraction(basis[i][k]) for k in range(len(b))]
    w = [0] * len(target)
    for i, ci in enumerate(coeffs):
        if ci:
            for k in range(len(w)):
                w[k] += ci * basis[i][k]
    return w

def recover_bdd(sigs, n, lam, k1_bound, k2_bound, d_secret, exact_cvp=False):
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    gens, t = build_bdd(sigs, n, lam, k1_bound, k2_bound)
    basis = lll_strip(gens, dim)
    neg_t = [-x for x in t]
    if exact_cvp:
        B = IntegerMatrix.from_matrix(basis)
        try:
            w = list(CVP.closest_vector(B, tuple(neg_t)))
        except Exception:
            return None
    else:
        w = babai_nearest_plane(basis, neg_t)
    v = [t[k] + w[k] for k in range(dim)]
    for cand in (v, [-x for x in v]):
        d = d_from_kblocks(cand, sigs, n, lam, S_K1, S_K2, m)
        if d is not None and d == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# Fresh-curve search (mirrors glv_hnp_phase2_lambda_threshold.py:373)
# ---------------------------------------------------------------------------

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

def search_curves(lo, hi, per_bin=1, nbins=8):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by
    lam* = min(lam,n-lam)/n into `nbins` bins over [0,0.5]."""
    import sympy
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
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
    for i in range(nbins): out.extend(bins[i])
    return out

# ---------------------------------------------------------------------------
# Gaussian-heuristic margin
# ---------------------------------------------------------------------------

def gh_margin(n, m, k1_bound, k2_bound):
    """gamma = eff * n^(1/m) / (3/(2*pi*e)).

    Derivation (BDD form, dim 2m):
      det L0   = n^(m-1) * S_K1^m * S_K2^m  ~ n^(3m-1)/(K1*K2)^m
      GH(L0)   = sqrt(2m/(2*pi*e)) * det^(1/(2m))
      ||v||    ~ n*sqrt(2m/3)          (k1,k2 uniform, E[x^2]=X^2/3)
      ||v|| < GH  <=>  eff * n^(1/m) < 3/(2*pi*e) = 0.17566...
    gamma < 1 predicts recovery."""
    eff = k1_bound * k2_bound / n
    return eff * n ** (1.0 / m) / (3.0 / (2.0 * math.pi * math.e))

# ---------------------------------------------------------------------------
# Trial driver
# ---------------------------------------------------------------------------

METHODS = ("KAN", "PROJ", "BABAI", "CVP")

def trial(curve, m, d_secret, k1_bound, seed, methods=METHODS):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {}
    if "KAN" in methods:
        rows = lll_strip(build_kan(sigs, n, lam, k1_bound, k2_bound), 2 * m + 2)
        out["KAN"] = recover_kan(rows, m, n, S_KAN, d_secret)
        out["_kan_sv"] = min(norm(r) for r in rows)
        out["_kan_pv"] = norm(planted_kan(sigs, d_secret, n, k1_bound, k2_bound))
    if "PROJ" in methods:
        rows = lll_strip(build_proj(sigs, n, lam, k1_bound, k2_bound), 2 * m + 1)
        out["PROJ"] = recover_proj(rows, sigs, m, n, lam, k1_bound, k2_bound, d_secret)
        out["_proj_sv"] = min(norm(r) for r in rows)
        out["_proj_pv"] = norm(planted_proj(sigs, n, k1_bound, k2_bound))
    if "BABAI" in methods:
        out["BABAI"] = recover_bdd(sigs, n, lam, k1_bound, k2_bound, d_secret, False)
    if "CVP" in methods:
        r = recover_bdd(sigs, n, lam, k1_bound, k2_bound, d_secret, True)
        out["CVP"] = r
    return out

def sweep(curve, m, k1_bound, seeds, methods=METHODS):
    p, b, n, lam, G = curve
    agg = {k: 0 for k in methods}
    ratios = {"KAN": [], "PROJ": []}
    used = 0
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, methods)
        if res is None: continue
        used += 1
        for k in methods:
            agg[k] += bool(res.get(k))
        for k in ("KAN", "PROJ"):
            if f"_{k.lower()}_pv" in res:
                pv = res[f"_{k.lower()}_pv"]
                if pv:
                    ratios[k].append(res[f"_{k.lower()}_sv"] / pv)
    mean = {k: (sum(v) / len(v) if v else float('nan')) for k, v in ratios.items()}
    return agg, used, mean


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 23 — projected / BDD reformulation of the Phase-2 lattice")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    # Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
    HIST = [
        # label,             p,    b, n,    lam,  K1_hist, m
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
    print("EXP P1: the three formulations, dimensions and determinants")
    print("-" * 78)
    print(f"{'curve':<18}{'m':>3}{'K1':>4}{'lam*':>8}{'nu_hat':>8}"
          f"{'dimKAN':>8}{'dimPROJ':>8}{'dimBDD':>8}{'gamma':>8}")
    for label, cur, k1, m in curves:
        p, b, n, lam, G = cur
        k2 = math.isqrt(n) + 1
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
        print(f"{label:<18}{m:>3}{k1:>4}{lam_star(lam, n):>8.3f}"
              f"{nu_hat(n, lam, S_K1, S_K2):>8.3f}"
              f"{2*m+2:>8}{2*m+1:>8}{2*m:>8}{gh_margin(n, m, k1, k2):>8.2f}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP P2: is the planted vector lambda_1 after projection?")
    print("  sv/pv = ||shortest reduced row|| / ||planted vector||.")
    print("  sv/pv >= 1 means the planted vector IS the shortest (falsifier arm 1)")
    print("-" * 78)
    print(f"{'curve':<18}{'K1':>4}{'m':>4}{'sv/pv KAN':>12}{'sv/pv PROJ':>12}"
          f"{'KAN':>6}{'PROJ':>6}")
    for label, cur, k1_hist, m in curves:
        p, b, n, lam, G = cur
        for k1 in (2, k1_hist):
            agg, used, mean = sweep(cur, m, k1, SEEDS, methods=("KAN", "PROJ"))
            print(f"{label:<18}{k1:>4}{m:>4}{mean['KAN']:>12.3f}"
                  f"{mean['PROJ']:>12.3f}{agg['KAN']:>4}/{used:<2}"
                  f"{agg['PROJ']:>4}/{used:<2}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP P3: the K1 wall, all four formulations (T4 grid replicated)")
    print("-" * 78)
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    walls = {}
    for label, cur, k1_hist, m in curves:
        p, b, n, lam, G = cur
        k2 = math.isqrt(n) + 1
        print(f"\n  {label}  n={n}  lam*={lam_star(lam,n):.3f}  "
              f"nu_hat={nu_hat(n, lam, *[scales(n,2,k2)[i] for i in (0,2)]):.3f}  m={m}  K2={k2}")
        print(f"    {'K1':>4}{'eff':>8}{'gamma':>8}   "
              f"{'KAN':>7}{'PROJ':>7}{'BABAI':>7}{'CVP':>7}")
        walls[label] = {}
        for k1 in K1_GRID:
            agg, used, _ = sweep(cur, m, k1, SEEDS)
            eff = k1 * k2 / n
            g = gh_margin(n, m, k1, k2)
            cells = "".join(f"{agg[k]:>5}/{used:<2}" for k in METHODS)
            print(f"    {k1:>4}{eff:>8.4f}{g:>8.2f}   {cells}")
            for k in METHODS:
                walls[label].setdefault(k, []).append((k1, agg[k], used))

    print("\n  Wall = largest K1 with 5/5 recovery:")
    print(f"    {'curve':<18}" + "".join(f"{k:>8}" for k in METHODS))
    for label in walls:
        cells = ""
        for k in METHODS:
            w = 0
            for k1, wins, used in walls[label][k]:
                if used and wins == used: w = k1
            cells += f"{w:>8}"
        print(f"    {label:<18}{cells}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP P4: m-scaling prediction.  gamma = eff*n^(1/m)/0.17566 < 1")
    print("  predicts recovery, so at fixed K1 the predicted signature count is")
    print("  m* = ln(n) / ln(0.17566/eff).  T4b (2026-07-29) tested m<=32 at")
    print("  K1=8 on 12-bit/2677 and saw 0-1 of 5; the formula says m* ~ 70.")
    print("-" * 78)
    label, cur, _, _ = curves[2]           # 12-bit/2677
    p, b, n, lam, G = cur
    k2 = math.isqrt(n) + 1
    for k1 in (8, 12):
        eff = k1 * k2 / n
        if eff < 3.0 / (2 * math.pi * math.e):
            m_star = math.log(n) / math.log((3.0 / (2 * math.pi * math.e)) / eff)
        else:
            m_star = float('inf')
        print(f"\n  {label}  K1={k1}  eff={eff:.4f}  predicted m* = {m_star:.1f}")
        print(f"    {'m':>5}{'gamma':>8}   {'KAN':>7}{'PROJ':>7}{'BABAI':>7}")
        for m in (10, 20, 40, 60, 70, 80, 100, 140):
            # Babai here uses exact rational GS, O(dim^3) on Fractions with
            # n-bit numerators -- unusable past dim ~ 80, so drop it there.
            meth = ("KAN", "PROJ", "BABAI") if m <= 40 else ("KAN", "PROJ")
            agg, used, _ = sweep(cur, m, k1, SEEDS[:3], methods=meth)
            g = gh_margin(n, m, k1, k2)
            cells = "".join(f"{agg[k]:>5}/{used:<2}" if k in meth else f"{'-':>8}"
                            for k in ("KAN", "PROJ", "BABAI"))
            print(f"    {m:>5}{g:>8.2f}   {cells}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP P5: calibrate the wall against gamma on fresh curves.")
    print("  For each curve: m=10, find the largest K1 with 5/5 (method KAN),")
    print("  report gamma at that K1 ('gamma_wall').  gamma_wall == 1 would mean")
    print("  the Gaussian heuristic is exact; gamma_wall > 1 measures the excess.")
    print("-" * 78)
    fresh = search_curves(2 ** 13, 2 ** 14, per_bin=1, nbins=8)
    print(f"  found {len(fresh)} fresh j=0 GLV curves in [2^13, 2^14)")
    print(f"  {'p':>7}{'n':>7}{'lam*':>7}{'nu@2':>7}{'K1wall':>8}"
          f"{'eff_w':>8}{'g_wall':>8}{'nu@wall':>9}")
    rows5 = []
    for cur in fresh:
        p, b, n, lam, G = cur
        k2 = math.isqrt(n) + 1
        m = 10
        k1_wall = 0
        for k1 in range(2, 33):
            agg, used, _ = sweep(cur, m, k1, SEEDS, methods=("KAN",))
            if used and agg["KAN"] == used:
                k1_wall = k1
            elif k1_wall:
                break
        if not k1_wall:
            continue
        g_w = gh_margin(n, m, k1_wall, k2)
        nu2 = nu_hat(n, lam, *[scales(n, 2, k2)[i] for i in (0, 2)])
        nuw = nu_hat(n, lam, *[scales(n, k1_wall, k2)[i] for i in (0, 2)])
        rows5.append((p, n, lam_star(lam, n), nu2, k1_wall,
                      k1_wall * k2 / n, g_w, nuw))
        print(f"  {p:>7}{n:>7}{rows5[-1][2]:>7.3f}{nu2:>7.3f}{k1_wall:>8}"
              f"{rows5[-1][5]:>8.4f}{g_w:>8.2f}{nuw:>9.3f}")

    if len(rows5) >= 4:
        gws = [r[6] for r in rows5]
        nus = [r[3] for r in rows5]
        lss = [r[2] for r in rows5]
        def spearman(x, y):
            def rank(v):
                order = sorted(range(len(v)), key=lambda i: v[i])
                r = [0.0] * len(v)
                for pos, i in enumerate(order): r[i] = pos + 1.0
                return r
            rx, ry = rank(x), rank(y)
            mx, my = sum(rx)/len(rx), sum(ry)/len(ry)
            num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
            den = math.sqrt(sum((a-mx)**2 for a in rx) * sum((b-my)**2 for b in ry))
            return num/den if den else float('nan')
        print(f"\n  gamma_wall: min={min(gws):.2f} median="
              f"{sorted(gws)[len(gws)//2]:.2f} max={max(gws):.2f}  (n={len(gws)})")
        print(f"  spearman(nu_hat@K1=2, gamma_wall) = {spearman(nus, gws):+.3f}")
        print(f"  spearman(lam*,         gamma_wall) = {spearman(lss, gws):+.3f}")

    print("\nDone.")
