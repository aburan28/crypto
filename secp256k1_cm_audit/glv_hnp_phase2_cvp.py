"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the target is reachable.

Motivation (RESEARCH_AUTOLAB_LOG.md 2026-07-29, EXP T5): in the Kannan-embedded
Phase-2 lattice the shortest vector is ALWAYS the trivial vector n*S_D*e_m --
100% of its energy in the d-column, Kannan coordinate 0.  Algebraically
    ||v_planted||^2 ~ n^2 (2m/3 + 4/3)   vs   ||n*S_D*e_m||^2 = n^2,
so the planted vector is never lambda_1 for any m >= 1 and no choice of S_D
changes that (both scale linearly in S_D).  Recovery is therefore not an SVP
condition but a BDD/coset condition.

Thread 23 removes the Kannan embedding and states the problem as what it
actually is:

    L0  = < n*S_K1*e_i ,  (B_i*S_K1 | S_D*e_m) ,  (-lam*S_K1*e_i | S_K2*e_{m+i}) >
    a   = (A_i*S_K1, 0, ..., 0)
    find w in L0 minimising ||a + w||;  then a + w = (k1_i*S_K1, d*S_D, k2_i*S_K2).

dim L0 = 2m+1 (one less than the Kannan lattice).  The trivial vector n*e_m is
still in L0, but in CVP it is no longer a rival: it is exactly the "reduce d
mod n" freedom, and it moves the d-coordinate by n, which INCREASES the
distance once d is taken in the symmetric range (-n/2, n/2].

Three arms are compared on the historical T4 grid:
    K  Kannan embedding + LLL          (2026-07-29 baseline, dim 2m+2)
    B  Babai nearest plane on LLL(L0)  (dim 2m+1)
    X  exact CVP by enumeration on L0  (dim 2m+1)

Arm X is the point of the exercise.  It separates the two branches of the
2026-07-29 falsifier:
    * X recovers where K fails            -> the wall was algorithmic;
      the Kannan/SVP formulation was leaving recoverable instances on the table.
    * X fails too, with dist(found) < dist(planted)
                                          -> the planted vector is not the
      closest vector at all: the wall is information-theoretic and Phase 2 is
      at its ceiling for that (K1, m).

The ratio r = dist(X) / dist(planted) is the CVP analogue of the 2026-07-29
sv/pv statistic, and it is the diagnostic that sv/pv could not be: r = 1 means
"planted IS optimal", r < 1 means "provably not enough information".

Lattice construction, signature generation and curve data are reused verbatim
from glv_hnp_phase2_lambda_threshold.py / glv_hnp_phase2_20bit.py so the
comparison to 2026-07-26 and 2026-07-29 is exact.

Run: python3 glv_hnp_phase2_cvp.py
"""

import math
import random
import time

from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim, glv_hnp_phase2_20bit.py)
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
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None:
            return P
    return None

def lam_star(lam, n):
    """Symmetric size of lam: min(lam, n-lam)/n in [0, 0.5]."""
    return min(lam % n, n - (lam % n)) / n

# --- CM curve search + the L2 block (verbatim, glv_hnp_phase2_lambda_threshold.py)

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (Eisenstein norm form)."""
    if p % 3 != 1:
        return None
    for a in range(1, math.isqrt(4 * p // 3) + 2):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            continue
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for b in ((a + s) // 2, (a - s) // 2):
            if a * a - a * b + b * b == p:
                return (a, b)
    return None

def j0_traces(a, b):
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n, or None."""
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

def gauss_reduce_2d(u, v):
    """Lagrange/Gauss reduction: exact shortest nonzero vector of a 2D lattice."""
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
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def search_curves(lo, hi, per_bin=2, nbins=10, max_primes=100000):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, bucketed by lam*."""
    import sympy
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    seen = 0
    while p < hi and seen < max_primes:
        seen += 1
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
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
# Phase-2 lattice (verbatim, glv_hnp_phase2_lambda_threshold.py:227-286)
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
    """Kannan-embedded Phase-2 lattice, dim 2m+2 (the 2026-07-29 baseline)."""
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

def recover_d_kannan(M_reduced, m, n, S_KANNAN, d_secret):
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
# Thread 23: the CVP reformulation
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """L0: the Kannan lattice with the embedding row AND column removed.

    Rows (dim 2m+1, columns 0..2m):
      i < m        : n*S_K1 at col i                     (mod-n freedom c_i)
      i = m        : B_j*S_K1 at col j, S_D at col m     (the d direction)
      m < i <= 2m  : -lam*S_K1 at col i-m-1, S_K2 at col i   (the k2 directions)

    det L0 = (n*S_K1)^m * S_D * S_K2^m.  Full rank 2m+1.
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for j in range(m):
        M[m][j] = sigs[j]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    return M

def cvp_offset(sigs, n, k1_bound, k2_bound):
    """a = (A_i*S_K1, 0, ..., 0).  The planted vector is a + w, w in L0."""
    m = len(sigs)
    S_K1, _, _, _ = scales(n, k1_bound, k2_bound)
    a = [0] * (2 * m + 1)
    for i in range(m):
        a[i] = sigs[i]['A'] * S_K1
    return a

def cvp_planted(sigs, d_secret, n, k1_bound, k2_bound):
    """The planted point a + w_planted = (k1_i*S_K1, d_sym*S_D, k2_i*S_K2).

    d is taken in the SYMMETRIC range (-n/2, n/2]: the coset of L0 containing
    the planted point also contains every d + t*n, and the closest one is the
    symmetric representative.  Using d_secret raw would compare the CVP answer
    against a non-optimal member of its own coset.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    d_sym = d_secret if d_secret <= n // 2 else d_secret - n
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_sym * S_D
    return v

def norm(v):
    return math.sqrt(sum(int(x) * int(x) for x in v))

def read_d(point, m, n, S_D):
    """Read the secret out of a CVP answer point a + w."""
    dc = int(point[m])
    if S_D != 1:
        if dc % S_D:
            return None
        dc //= S_D
    return dc % n

# ---------------------------------------------------------------------------
# The three arms
# ---------------------------------------------------------------------------

def arm_kannan(sigs, n, lam, k1_bound, k2_bound, d_secret, use_bkz=False, beta=20):
    m = len(sigs)
    dim = 2 * m + 2
    _, _, _, S_KANNAN = scales(n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(build_glv_lattice(sigs, n, lam, k1_bound, k2_bound))
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d_kannan(reduced, m, n, S_KANNAN, d_secret) is not None

def _cvp_common(sigs, n, lam, k1_bound, k2_bound):
    """LLL-reduce L0 once and build the CVP target t = -a."""
    A = IntegerMatrix.from_matrix(build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound))
    LLL.reduction(A)
    a = cvp_offset(sigs, n, k1_bound, k2_bound)
    t = [-x for x in a]
    return A, a, t

def arm_babai(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Babai nearest plane on LLL(L0).  Returns (recovered, distance)."""
    m = len(sigs)
    _, S_D, _, _ = scales(n, k1_bound, k2_bound)
    A, a, t = _cvp_common(sigs, n, lam, k1_bound, k2_bound)
    w = CVP.babai(A, t)
    point = [a[j] + int(w[j]) for j in range(2 * m + 1)]
    return read_d(point, m, n, S_D) == d_secret, norm(point)

def arm_exact(sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Exact CVP by enumeration on L0.  Returns (recovered, distance)."""
    m = len(sigs)
    _, S_D, _, _ = scales(n, k1_bound, k2_bound)
    A, a, t = _cvp_common(sigs, n, lam, k1_bound, k2_bound)
    w = CVP.closest_vector(A, t)
    point = [a[j] + int(w[j]) for j in range(2 * m + 1)]
    return read_d(point, m, n, S_D) == d_secret, norm(point)

# ---------------------------------------------------------------------------
# Trial driver
# ---------------------------------------------------------------------------

def trial(curve, m, k1_bound, seed, arms=("K", "B", "X")):
    """One (curve, m, K1, seed) instance through the selected arms.

    d_secret is drawn exactly as in glv_hnp_phase2_lambda_threshold.py's
    success_rate(), so the K arm reproduces the 2026-07-29 T4 numbers.
    """
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {'planted': norm(cvp_planted(sigs, d_secret, n, k1_bound, k2_bound))}
    if "K" in arms:
        out['K'] = arm_kannan(sigs, n, lam, k1_bound, k2_bound, d_secret)
    if "B" in arms:
        out['B'], out['B_dist'] = arm_babai(sigs, n, lam, k1_bound, k2_bound, d_secret)
    if "X" in arms:
        out['X'], out['X_dist'] = arm_exact(sigs, n, lam, k1_bound, k2_bound, d_secret)
    return out

def cell(curve, m, k1_bound, seeds, arms=("K", "B", "X")):
    """Aggregate over seeds.  Returns dict arm -> wins, plus mean r = X/planted."""
    res = {a: 0 for a in arms}
    ratios = []
    tot = 0
    for s in seeds:
        r = trial(curve, m, k1_bound, s, arms)
        if r is None:
            continue
        tot += 1
        for a in arms:
            res[a] += bool(r[a])
        if "X" in arms:
            ratios.append(r['X_dist'] / r['planted'] if r['planted'] else float('nan'))
    res['_total'] = tot
    res['_r'] = sum(ratios) / len(ratios) if ratios else float('nan')
    # how many trials had the planted point EXACTLY optimal (r == 1)
    res['_r1'] = sum(1 for x in ratios if x == 1.0)
    return res

def min_m_for_full_recovery(curve, k1_bound, seeds, m_grid, arm="X"):
    """Smallest m in m_grid at which `arm` recovers on every seed."""
    for m in m_grid:
        c = cell(curve, m, k1_bound, seeds, arms=(arm,))
        if c['_total'] and c[arm] == c['_total']:
            return m, c
    return None, None


# ===========================================================================
print("=" * 78)
print("Thread 23 — CVP reformulation of the GLV-HNP Phase-2 lattice")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
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


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R1: sanity — is the planted point really in the CVP coset a + L0?")
print("-" * 78)
print("Checks that a + w_planted is a lattice point of L0 shifted by a, i.e.")
print("that solving CVP(L0, -a) is the same problem the Kannan row encoded.\n")

ok_all = True
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1, k2_bound, 42)
    L0 = build_cvp_lattice(sigs, n, lam, k1, k2_bound)
    a = cvp_offset(sigs, n, k1, k2_bound)
    v = cvp_planted(sigs, d_secret, n, k1, k2_bound)
    # w = v - a must be an integer combination of the rows of L0.  Solve by
    # back-substitution: the k2-rows fix their own columns, then the d-row,
    # then the mod-n rows.
    S_K1, S_D, S_K2, _ = scales(n, k1, k2_bound)
    w = [v[j] - a[j] for j in range(2 * m + 1)]
    coef = [0] * (2 * m + 1)
    for i in range(m):                      # k2 rows: column m+1+i is theirs alone
        assert w[m + 1 + i] % S_K2 == 0
        coef[m + 1 + i] = w[m + 1 + i] // S_K2
    assert w[m] % S_D == 0
    coef[m] = w[m] // S_D                   # d row: column m is its alone
    for i in range(m):                      # mod-n rows absorb the rest
        rem = w[i] - coef[m] * L0[m][i] - sum(coef[m + 1 + j] * L0[m + 1 + j][i]
                                              for j in range(m))
        if rem % (n * S_K1) != 0:
            ok_all = False
            print(f"  {label:<18} FAIL: residue not divisible at col {i}")
            break
        coef[i] = rem // (n * S_K1)
    else:
        recon = [sum(coef[r] * L0[r][j] for r in range(2 * m + 1)) + a[j]
                 for j in range(2 * m + 1)]
        match = recon == v
        ok_all &= match
        print(f"  {label:<18} m={m:<3} dim(L0)={2*m+1:<3} coset check: "
              f"{'OK' if match else 'FAIL'}   "
              f"||planted||={norm(v):.3e}  ||n*e_m||={n*S_D:.3e}  "
              f"ratio={n*S_D/norm(v):.3f}")
print(f"\n  All coset checks passed: {ok_all}")
print("  (ratio < 1 is exactly the 2026-07-29 T5 obstruction: n*e_m is shorter")
print("   than the planted vector.  In CVP it is no longer a rival — it moves")
print("   the d-coordinate by n and therefore increases the distance.)")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R2: the T4 falsifier grid — does the reformulation move the K1 wall?")
print("-" * 78)
print("2026-07-29 T4, Kannan+LLL, m=12, 5 seeds:")
print("    12-bit/2557 (lam*=0.340): 5/5 up to K1=8, 4/5 at 12, 1/5 at 16, 0/5 at 24")
print("    12-bit/2677 (lam*=0.070): 5/5 up to K1=4, 2/5 at 6,  0/5 at K1>=8")
print("Falsifier: if B/X push those walls outward the Kannan formulation was")
print("leaving recoverable instances on the table; if not, the wall is")
print("information-theoretic.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_R2 = 12
grid_rows = []
hdr = f"{'curve':<18} {'lam*':>6} {'arm':>4} " + " ".join(f"K1={k:<5}" for k in K1_GRID)
print(hdr)
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    if n < 2000:
        continue                       # 8-bit curve has K2=15, grid not comparable
    cells = {}
    t0 = time.time()
    for k1 in K1_GRID:
        cells[k1] = cell(curve, M_R2, k1, SEEDS)
    for arm in ("K", "B", "X"):
        line = " ".join(f"{cells[k][arm]}/{cells[k]['_total']}  " for k in K1_GRID)
        print(f"{label:<18} {lam_star(lam,n):>6.4f} {arm:>4} " + line)
    print(f"{'':<18} {'':>6} {'r':>4} "
          + " ".join(f"{cells[k]['_r']:<6.4f}" for k in K1_GRID))
    grid_rows.append((label, lam_star(lam, n), cells, time.time() - t0))
    print()

print("legend: K = Kannan+LLL (dim 2m+2, SVP)   B = Babai on LLL(L0) (dim 2m+1)")
print("        X = exact CVP on L0 (enumeration)")
print("        r = mean ||exact CVP answer|| / ||planted point||")
print("            r = 1.0000  -> planted IS the closest vector (recovery is")
print("                           then purely an algorithmic question)")
print("            r < 1       -> a strictly closer point exists: the planted")
print("                           vector is not determined by the instance,")
print("                           i.e. an information-theoretic wall")
print(f"\n(eff = K1*K2/n with K2=52, n~2650: "
      + ", ".join(f"K1={k}:{k*52/2647:.3f}" for k in K1_GRID) + ")")
for label, ls, cells, dt in grid_rows:
    print(f"  timing {label}: {dt:.1f}s for {len(K1_GRID)}x{len(SEEDS)} trials x 3 arms")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("R3: where is the crossover?  r as a function of K1, per seed")
print("-" * 78)
print("If the wall is information-theoretic there should be a K1 at which r")
print("drops below 1 and stays there, coinciding with the X-arm failure.\n")

fail_curve = [c for lbl, c, _, _ in hist if c[2] == 2647][0]
good_curve = [c for lbl, c, _, _ in hist if c[2] == 2659][0]
for name, cv in (("2677 (lam*=0.070)", fail_curve), ("2557 (lam*=0.340)", good_curve)):
    print(f"  curve {name}, m={M_R2}")
    print(f"    {'K1':>4} " + " ".join(f"{'s'+str(i):>8}" for i in range(len(SEEDS)))
          + f" {'X wins':>8}")
    for k1 in K1_GRID:
        rs, wins = [], 0
        for s in SEEDS:
            r = trial(cv, M_R2, k1, s, arms=("X",))
            if r is None:
                rs.append(float('nan')); continue
            rs.append(r['X_dist'] / r['planted'] if r['planted'] else float('nan'))
            wins += bool(r['X'])
        print(f"    {k1:>4} " + " ".join(f"{x:>8.4f}" for x in rs) + f" {wins:>8}")
    print()


# ---------------------------------------------------------------------------
print("-" * 78)
print("R4: does more data rescue the CVP arms at the wall?")
print("-" * 78)
print("2026-07-29 T4b: at K1=8 the lam*=0.070 curve stayed dead under Kannan+LLL")
print("for m = 8,12,16,24,32 (0,0,1,0,1 of 5).  Same grid, CVP arms.\n")

M_GRID_R4 = (6, 8, 10, 12, 14, 16, 20, 24, 32)
for k1_try in (6, 8, 12, 16):
    eff = k1_try * (math.isqrt(2647) + 1) / 2647
    m_star = math.ceil(math.log(2647) / math.log(1.0 / eff))
    print(f"  --- K1={k1_try}  (eff={eff:.3f}, counting bound m* = {m_star}) ---")
    print(f"  {'m':>4} {'K':>8} {'B':>8} {'X':>8} {'mean r':>9} {'#(r=1)':>8}")
    for m_try in M_GRID_R4:
        c = cell(fail_curve, m_try, k1_try, SEEDS)
        print(f"  {m_try:>4} {c['K']}/{c['_total']:<6} {c['B']}/{c['_total']:<6} "
              f"{c['X']}/{c['_total']:<6} {c['_r']:>9.4f} {c['_r1']:>6}/{c['_total']}")
    print()


# ---------------------------------------------------------------------------
print("-" * 78)
print("R5: the K1 wall re-measured — minimum m for 5/5 recovery, per arm")
print("-" * 78)
print("If the 'K1 wall' is really a fixed-m artefact of the Kannan formulation,")
print("then under exact CVP every K1 with eff < 1 should be recoverable once m")
print("exceeds the counting bound m* = ceil(log n / log(1/eff)).  Any K1 where")
print("no m in the grid works is a genuine wall.\n")

M_SEARCH = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48)
K1_R5 = [2, 3, 4, 6, 8, 12, 16, 24]   # eff > 0.5 costs minutes/cell and is far
                                      # past the wall on both curves; omitted.
for name, cv in (("2677 (lam*=0.070)", fail_curve), ("2557 (lam*=0.340)", good_curve)):
    p, b, n, lam, G = cv
    K2 = math.isqrt(n) + 1
    print(f"  curve {name}, n={n}, K2={K2}")
    print(f"    {'K1':>4} {'eff':>7} {'m*':>4} {'min m (K)':>10} {'min m (B)':>10} "
          f"{'min m (X)':>10} {'t(s)':>7}")
    for k1 in K1_R5:
        eff = k1 * K2 / n
        if eff >= 1.0:
            print(f"    {k1:>4} {eff:>7.3f} {'inf':>4} "
                  f"{'--':>10} {'--':>10} {'--':>10} {'':>7}   (eff>=1: no info)")
            continue
        m_star = math.ceil(math.log(n) / math.log(1.0 / eff))
        t0 = time.time()
        mk, _ = min_m_for_full_recovery(cv, k1, SEEDS, M_SEARCH, arm="K")
        mb, _ = min_m_for_full_recovery(cv, k1, SEEDS, M_SEARCH, arm="B")
        mx, _ = min_m_for_full_recovery(cv, k1, SEEDS, M_SEARCH, arm="X")
        f = lambda v: str(v) if v is not None else f">{M_SEARCH[-1]}"
        print(f"    {k1:>4} {eff:>7.3f} {m_star:>4} {f(mk):>10} {f(mb):>10} "
              f"{f(mx):>10} {time.time()-t0:>7.1f}")
    print()


# ---------------------------------------------------------------------------
print("-" * 78)
print("R6: why does the wall sit where it does?  Gaussian-heuristic prediction")
print("-" * 78)
print("""L0 has det = (n*S_K1)^m * S_D * S_K2^m = n^(3m) / (K1*K2)^m and dim 2m+1,
so as m grows  det^(1/dim) -> n^1.5 / (K1*K2)^0.5 = n / sqrt(eff),  and

    gh(L0) = sqrt(dim / (2*pi*e)) * det^(1/dim) ~ sqrt((2m+1)/(2*pi*e)) * n/sqrt(eff).

The planted point has  E||v||^2 = m*n^2/3 + m*n^2/3 + n^2/12,  so ||v|| ~ n*sqrt(2m/3).
Their ratio is asymptotically INDEPENDENT of m:

    ||v_planted|| / gh(L0)  ->  sqrt(2*pi*e/3) * sqrt(eff)  =  2.387 * sqrt(eff).

The planted point is the unique closest vector only while that ratio is < 1, i.e.

    eff  <  eff_c = 3/(2*pi*e) = 0.17549...

This is a *bias-strength* wall, not a data wall: no amount of extra signatures
moves it, which is exactly what R5 shows (>48 well above the counting bound m*).
""")

EFF_C = 3.0 / (2.0 * math.pi * math.e)
print(f"  eff_c = 3/(2*pi*e) = {EFF_C:.5f}\n")

def gh_ratio(n, k1_bound, k2_bound, m):
    """Finite-m  E||v_planted|| / gh(L0), using the actual integer scales."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    log_det = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2)
    gh = math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)
    pn = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3.0
                   + (n * S_D) ** 2 / 12.0
                   + m * (k2_bound * S_K2) ** 2 / 3.0)
    return pn / gh

def gh_min_m(n, k1_bound, k2_bound, m_grid):
    for m in m_grid:
        if gh_ratio(n, k1_bound, k2_bound, m) < 1.0:
            return m
    return None

print("  Finite-m form (uses the actual integer scales, not the asymptote).")
print("  'pred m' = smallest m in the R5 grid with E||planted||/gh(L0) < 1.\n")
print(f"  {'curve':<20} {'K1':>4} {'eff':>7} {'asym':>6} {'gh@m=32':>8} "
      f"{'pred m':>7} {'obs m (X)':>10} {'verdict':>9}")
R5_OBS = {                       # (curve label, K1) -> min m under exact CVP, from R5
    ("2677", 2): 6,  ("2677", 3): 6,  ("2677", 4): 8,  ("2677", 6): 16,
    ("2677", 8): 24, ("2677", 12): None, ("2677", 16): None, ("2677", 24): None,
    ("2557", 2): 6,  ("2557", 3): 6,  ("2557", 4): 6,  ("2557", 6): 8,
    ("2557", 8): 8,  ("2557", 12): 14, ("2557", 16): 28, ("2557", 24): None,
}
M_SEARCH_R6 = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48)
for name, cv in (("2677 (lam*=0.070)", fail_curve), ("2557 (lam*=0.340)", good_curve)):
    p, b, n, lam, G = cv
    K2 = math.isqrt(n) + 1
    key = "2677" if n == 2647 else "2557"
    for k1 in (2, 3, 4, 6, 8, 12, 16, 24):
        eff = k1 * K2 / n
        asym = math.sqrt(2 * math.pi * math.e / 3.0) * math.sqrt(eff)
        pm = gh_min_m(n, k1, K2, M_SEARCH_R6)
        obs_m = R5_OBS[(key, k1)]
        agree = (pm is None) == (obs_m is None)
        print(f"  {name:<20} {k1:>4} {eff:>7.3f} {asym:>6.3f} "
              f"{gh_ratio(n, k1, K2, 32):>8.3f} {str(pm) if pm else '>48':>7} "
              f"{str(obs_m) if obs_m else '>48':>10} "
              f"{'agree' if agree else 'MISMATCH':>9}")
    print()

print("  nu_hat (2026-07-29 separator) for the same instances:")
print("  L2 = <(n*S_K1,0), (-lam*S_K1, S_K2)>, nu_hat = lambda_1(L2)/sqrt(det L2)")
print(f"  {'curve':<20} {'K1':>4} {'nu_hat':>8}")
for name, cv in (("2677 (lam*=0.070)", fail_curve), ("2557 (lam*=0.340)", good_curve)):
    p, b, n, lam, G = cv
    K2 = math.isqrt(n) + 1
    for k1 in (8, 12, 16):
        S_K1, _, S_K2, _ = scales(n, k1, K2)
        w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
        mu = math.sqrt(w[0] ** 2 + w[1] ** 2)
        det2 = n * S_K1 * S_K2
        print(f"  {name:<20} {k1:>4} {mu / math.sqrt(det2):>8.4f}")

print("""
  Reading: the GH threshold predicts the 2677 (low lam*) wall exactly -- it
  falls between eff=0.157 (recovers at m=24) and eff=0.236 (dead at m=48).
  The 2557 curve BEATS the prediction, recovering out to eff=0.313.  A curve
  can only beat the Gaussian heuristic if L0 is non-random in a way that helps,
  and the known non-randomness is the L2 block: nu_hat below.
""")


# ---------------------------------------------------------------------------
print("-" * 78)
print("R7: does the eff wall track nu_hat across curves?")
print("-" * 78)
print("Cross-curve test of the R6 reading.  For each curve, measure the largest")
print("eff at which exact CVP still recovers 3/3 (m searched up to 24), and")
print("compare with nu_hat.  Prediction from 2026-07-29 (low nu_hat = easier):")
print("negative correlation between nu_hat and the eff wall.\n")

R7_SEEDS = [42, 1234, 9999]
R7_M = (6, 8, 12, 16, 24)
R7_EFF = [0.08, 0.12, 0.16, 0.24, 0.31, 0.39]

print("  building test curves (12-14 bit, lam* spread)...")
curves_r7 = search_curves(2**11, 2**14, per_bin=1, nbins=8)
print(f"  got {len(curves_r7)} curves\n")

print(f"  {'p':>7} {'n':>7} {'lam*':>7} {'nu_hat':>8} {'eff wall':>9} {'m at wall':>10}")
r7_rows = []
for (p, b, n, lam, G) in curves_r7:
    cv = (p, b, n, lam, G)
    K2 = math.isqrt(n) + 1
    # nu_hat is evaluated at the K1 nearest eff=0.16 (the GH threshold region)
    k1_ref = max(2, round(0.16 * n / K2))
    S_K1, _, S_K2, _ = scales(n, k1_ref, K2)
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    nu_hat = math.sqrt(w[0] ** 2 + w[1] ** 2) / math.sqrt(n * S_K1 * S_K2)
    wall_eff, wall_m = 0.0, None
    for eff_t in R7_EFF:
        k1 = max(2, round(eff_t * n / K2))
        if k1 * K2 / n >= 1.0:
            break
        mx, _ = min_m_for_full_recovery(cv, k1, R7_SEEDS, R7_M, arm="X")
        if mx is None:
            break
        wall_eff, wall_m = k1 * K2 / n, mx
    print(f"  {p:>7} {n:>7} {lam_star(lam,n):>7.4f} {nu_hat:>8.4f} "
          f"{wall_eff:>9.3f} {str(wall_m) if wall_m else '--':>10}")
    r7_rows.append((nu_hat, wall_eff, lam_star(lam, n)))

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
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')

if len(r7_rows) >= 4:
    nus = [r[0] for r in r7_rows]
    walls = [r[1] for r in r7_rows]
    lss = [r[2] for r in r7_rows]
    print(f"\n  spearman(nu_hat, eff wall)  = {spearman(nus, walls):+.3f}   "
          f"(2026-07-29 predicts negative)")
    print(f"  spearman(lam*,   eff wall)  = {spearman(lss, walls):+.3f}   "
          f"(lam* was falsified 2026-07-29; expect ~0)")
    print(f"  eff wall range observed: {min(walls):.3f} .. {max(walls):.3f}   "
          f"(GH prediction: {EFF_C:.3f} for a random lattice)")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
