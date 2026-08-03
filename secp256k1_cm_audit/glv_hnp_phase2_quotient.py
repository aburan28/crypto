"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector
can be lambda_1 -- quotient out the trivial direction, and measure the
information-theoretic ceiling with exact CVP.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 lattice L (dim 2m+2, cols: m k1-slots | d | m k2-slots | Kannan)
  always has shortest vector  n*S_D*e_m, which carries no information (d is
  only defined mod n).  The planted vector is therefore NEVER lambda_1, and
  recovery is a BDD/coset condition, not an SVP condition.  That run proposed
  reformulating the lattice; this script does it.

THE REFORMULATION (exact, not a heuristic):

  Claim.  L intersect Z*e_m  =  Z * (n*S_D*e_m).
  Proof.  Write v = sum_i a_i*row_i.  Vanishing on the k2 columns forces the
  coefficients of rows m+1..2m to 0 (those rows are the only ones with a
  nonzero k2 entry, and it is S_K2 on a distinct column each).  Vanishing on
  the Kannan column forces a_{2m+1} = 0.  On k1 column i what is left is
  a_i*n*S_K1 + a_m*B_i*S_K1 = 0, i.e. a_i = -a_m*B_i/n, which needs n | a_m*B_i
  for every i; gcd(B_i,n)=1 (B_i = r_i/s_i mod n, n prime), so n | a_m.  The
  generator is n*row_m - sum_i B_i*row_i = n*S_D*e_m.  QED

  Hence the coordinate-forgetting map  pi: Z^{2m+2} -> Z^{2m+1}  that DELETES
  column m has kernel exactly Z*(n*S_D*e_m) on L, so

        V1 := pi(L)  ==  L / (Z * n*S_D*e_m)

  is an honest integer lattice of rank 2m+1 realising the quotient.  The
  trivial vector maps to 0; the planted vector maps to
  (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN), and d is still recoverable, because any
  single signature gives  d = B_i^{-1}*(k1_i + lam*k2_i - A_i)  mod n.
  det(V1) = det(L)/(n*S_D) = (n*S_K1)^m * S_K2^m * S_KANNAN / n.

  V2 goes one step further: drop the Kannan row/column too and solve the
  underlying CVP directly (exact, by enumeration).  Target
  t = (A_i*S_K1 | 0), lattice spanned by {n*S_K1*e_i}, (B_i*S_K1), and
  {-lam*S_K1*e_i + S_K2*e_{m+i}}; then t - closest = (k1_i*S_K1 | k2_i*S_K2).
  Exact CVP is a strictly stronger oracle than LLL, so V2 measures the
  INFORMATION-THEORETIC ceiling: if V2 fails, no amount of reduction helps.

Variants compared:
  V0  baseline (2026-07-26 / 07-29 construction), LLL, Kannan, d column
  V1  quotient (col m deleted), LLL, Kannan
  V2  quotient + no Kannan, exact CVP  <-- information-theoretic ceiling

Falsifier stated by the 2026-07-29 entry:
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_quotient.py
"""

import math
import random
import sys
import time

import sympy
from fpylll import IntegerMatrix, LLL, SVP, CVP, GSO, Enumeration

# ---------------------------------------------------------------------------
# Minimal EC arithmetic -- verbatim from glv_hnp_phase2_lambda_threshold.py
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
        if ec_mul(P, n, p) is None:
            return P
    return None

# ---------------------------------------------------------------------------
# CM / GLV helpers -- verbatim
# ---------------------------------------------------------------------------

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

def glv_roots(n):
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
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Lattice construction (V0 verbatim; V1/V2 derived)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- verbatim from the 2026-07-26 attack."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    """Verbatim from glv_hnp_phase2_lambda_threshold.py:230 so that a given
    (curve, seed) produces byte-identical signatures across runs."""
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

def build_v0(sigs, n, lam, k1_bound, k2_bound):
    """Baseline lattice: dim 2m+2.
    cols  0..m-1   k1 slots (scale S_K1)
    col   m        d slot   (scale S_D)
    cols  m+1..2m  k2 slots (scale S_K2)
    col   2m+1     Kannan   (scale S_KANNAN)"""
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

def build_v1(sigs, n, lam, k1_bound, k2_bound):
    """Quotient lattice pi(L) = L / (Z * n*S_D*e_m): V0 with column m deleted.
    2m+2 generating rows of rank 2m+1 (one exact dependency:
    n*row_m - sum_i B_i*row_i = 0).
    cols  0..m-1     k1 slots
    cols  m..2m-1    k2 slots
    col   2m         Kannan"""
    m = len(sigs)
    M0 = build_v0(sigs, n, lam, k1_bound, k2_bound)
    return [row[:m] + row[m + 1:] for row in M0]

def build_v2(sigs, n, lam, k1_bound, k2_bound):
    """Quotient lattice with the Kannan row/column dropped: the underlying CVP.
    Returns (generators, target).
    cols 0..m-1 k1 slots, cols m..2m-1 k2 slots."""
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    gens = []
    for i in range(m):
        row = [0] * dim
        row[i] = n * S_K1
        gens.append(row)
    row = [0] * dim
    for i in range(m):
        row[i] = sigs[i]['B'] * S_K1
    gens.append(row)
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        gens.append(row)
    target = [0] * dim
    for i in range(m):
        target[i] = sigs[i]['A'] * S_K1
    return gens, target

# ---------------------------------------------------------------------------
# Planted vectors
# ---------------------------------------------------------------------------

def planted_v0(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def planted_v1(sigs, d_secret, n, k1_bound, k2_bound):
    v0 = planted_v0(sigs, d_secret, n, k1_bound, k2_bound)
    m = len(sigs)
    return v0[:m] + v0[m + 1:]

def planted_v2(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    return v

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_v0(rows, m, n, S_KAN, d_secret):
    for row in rows:
        last = row[2 * m + 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

def d_from_k(w, m, sigs, n, lam, S_K1, S_K2):
    """Given a candidate quotient-vector w = (k1_i*S_K1 | k2_i*S_K2 | ...),
    read off (k1_i,k2_i) and solve d = B_i^{-1}(k1_i + lam*k2_i - A_i) mod n.
    Returns the set of candidates over all i (an attacker checks each against
    the public key; here the caller compares to d_secret, the same oracle the
    baseline recover_v0 uses)."""
    cands = set()
    for i in range(m):
        k1 = int(round(w[i] / S_K1))
        k2 = int(round(w[m + i] / S_K2))
        try:
            d = (k1 + lam * k2 - sigs[i]['A']) * modinv(sigs[i]['B'], n) % n
        except ValueError:
            continue
        if d:
            cands.add(d)
    return cands

def recover_v1(rows, m, sigs, n, lam, S_K1, S_K2, S_KAN, d_secret):
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        w = [x if last > 0 else -x for x in row]
        if d_secret in d_from_k(w, m, sigs, n, lam, S_K1, S_K2):
            return d_secret
    return None

# ---------------------------------------------------------------------------
# Reduction drivers
# ---------------------------------------------------------------------------

def lll_rows(M):
    dim_r, dim_c = len(M), len(M[0])
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(dim_c)] for i in range(dim_r)]

def independent_basis(M):
    """LLL-reduce a (possibly rank-deficient) generating set and drop the
    zero rows fplll pushes to the top.  Returns an IntegerMatrix basis."""
    dim_c = len(M[0])
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim_c)] for i in range(A.nrows)]
    rows = [r for r in rows if any(r)]
    B = IntegerMatrix.from_matrix(rows)
    LLL.reduction(B)
    return B

def exact_svp(M, dim_cap=40):
    """Exact lambda_1 of the lattice generated by M (None if too big/failed).

    NOTE: SVP.shortest_vector() is unusable in this container -- fpylll's wheel
    does not ship fplll's pruning-strategy JSON, so it raises
    FileNotFoundError('.../strategies/default.json').  Enumeration takes no
    strategy file, so we enumerate directly inside a radius of ||b_0||, which
    is guaranteed to contain lambda_1."""
    if len(M[0]) > dim_cap:
        return None
    try:
        B = independent_basis(M)
        Mg = GSO.Mat(B)
        Mg.update_gso()
        radius = Mg.get_r(0, 0) * (1 + 1e-9)
        sol = Enumeration(Mg).enumerate(0, B.nrows, radius, 0)
        return math.sqrt(sol[0][0]) if sol else None
    except Exception:
        return None

def exact_cvp(gens, target, dim_cap=40):
    """Exact closest vector; returns the difference target - closest."""
    if len(target) > dim_cap:
        return None
    try:
        B = independent_basis(gens)
        c = CVP.closest_vector(B, tuple(int(x) for x in target))
        return [int(target[j]) - int(c[j]) for j in range(len(target))]
    except Exception:
        return None

# ---------------------------------------------------------------------------
# One trial, all three variants
# ---------------------------------------------------------------------------

CVP_AUDIT = {'n': 0, 'violations': 0}


def trial(curve, m, d_secret, k1_bound, seed, want=('v0', 'v1', 'v2')):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    out = {'m': m, 'n': n, 'K1': k1_bound, 'K2': k2_bound,
           'eff': k1_bound * k2_bound / n}

    if 'v0' in want:
        M0 = build_v0(sigs, n, lam, k1_bound, k2_bound)
        r0 = lll_rows(M0)
        pv0 = norm(planted_v0(sigs, d_secret, n, k1_bound, k2_bound))
        sv0 = min(norm(r) for r in r0 if any(r))
        out['v0'] = {'ok': recover_v0(r0, m, n, S_KAN, d_secret) is not None,
                     'pv': pv0, 'sv': sv0, 'ratio': sv0 / pv0}

    if 'v1' in want:
        M1 = build_v1(sigs, n, lam, k1_bound, k2_bound)
        r1 = lll_rows(M1)
        pv1 = norm(planted_v1(sigs, d_secret, n, k1_bound, k2_bound))
        nz = [r for r in r1 if any(r)]
        sv1 = min(norm(r) for r in nz)
        out['v1'] = {'ok': recover_v1(r1, m, sigs, n, lam, S_K1, S_K2, S_KAN,
                                      d_secret) is not None,
                     'pv': pv1, 'sv': sv1, 'ratio': sv1 / pv1,
                     'zero_rows': len(r1) - len(nz)}

    if 'v2' in want:
        gens, target = build_v2(sigs, n, lam, k1_bound, k2_bound)
        diff = exact_cvp(gens, target)
        pv2 = norm(planted_v2(sigs, n, k1_bound, k2_bound))
        if diff is None:
            out['v2'] = {'ok': None, 'pv': pv2, 'sv': float('nan'),
                         'ratio': float('nan')}
        else:
            ok = d_secret in d_from_k(diff, m, sigs, n, lam, S_K1, S_K2)
            dn = norm(diff)
            # Exactness self-check: the planted difference is a legal
            # candidate, so an EXACT CVP oracle can never return a longer one.
            CVP_AUDIT['n'] += 1
            if dn > pv2 * (1 + 1e-9):
                CVP_AUDIT['violations'] += 1
            out['v2'] = {'ok': ok, 'pv': pv2, 'sv': dn,
                         'ratio': dn / pv2 if pv2 else float('nan')}
    return out

def rates(curve, m, k1_bound, seeds, want=('v0', 'v1', 'v2')):
    acc = {k: {'win': 0, 'tot': 0, 'ratio': []} for k in want}
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve[2] - 1)
        res = trial(curve, m, d_trial, k1_bound, seed, want)
        if res is None:
            continue
        for k in want:
            acc[k]['tot'] += 1
            if res[k]['ok']:
                acc[k]['win'] += 1
            acc[k]['ratio'].append(res[k]['ratio'])
    for k in want:
        r = [x for x in acc[k]['ratio'] if x == x]
        acc[k]['mean_ratio'] = sum(r) / len(r) if r else float('nan')
    return acc

# ---------------------------------------------------------------------------
# Gaussian-heuristic predictor for V1
# ---------------------------------------------------------------------------

def gh_v1(m, n, k1_bound, k2_bound):
    """Gaussian heuristic for lambda_1 of the quotient lattice V1.
    det(V1) = det(V0)/(n*S_D) = (n*S_K1)^m * S_K2^m * S_KANNAN / n,
    rank 2m+1."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    log_det = m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_KAN) \
        - math.log(n * S_D)
    return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(log_det / dim)

def pv1_expected(m, n, k1_bound, k2_bound):
    """E[||planted||] in V1: k1_i*S_K1 and k2_i*S_K2 are ~U[0,n), Kannan = n."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    e2 = m * (k1_bound * S_K1) ** 2 / 3.0 \
        + m * (k2_bound * S_K2) ** 2 / 3.0 + float(S_KAN) ** 2
    return math.sqrt(e2)

def gamma_v1(m, n, k1_bound, k2_bound):
    """gamma = GH(V1) / E||planted||.  gamma > 1 predicts unique-SVP recovery."""
    return gh_v1(m, n, k1_bound, k2_bound) / pv1_expected(m, n, k1_bound, k2_bound)

# ---------------------------------------------------------------------------
# nu_hat -- the 2026-07-29 run #2 separator (recovered from commit e845207)
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Lagrange-Gauss reduction of a 2-D lattice; returns the shorter basis
    vector, which is exactly lambda_1."""
    def dot(a, b): return a[0] * b[0] + a[1] * b[1]
    while True:
        if dot(v, v) < dot(u, u):
            u, v = v, u
        q = round(dot(u, v) / dot(u, u))
        if q == 0:
            return u
        v = (v[0] - q * u[0], v[1] - q * u[1])

def nu_hat(n, lam, S_K1, S_K2):
    """nu_hat = lambda_1(L2)/sqrt(det L2) for the non-planted 2-D sublattice
    L2 = <(n*S_K1, 0), (-lam*S_K1, S_K2)>  spanned by rows i and m+1+i.
    det L2 = n*S_K1*S_K2 is independent of lam, so nu_hat isolates the
    lam-dependence.  Low nu_hat (skew L2) predicts EASIER recovery."""
    sv = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(sv[0] ** 2 + sv[1] ** 2) / math.sqrt(float(n) * S_K1 * S_K2)

def spearman(xs, ys):
    """Spearman rank correlation (average ranks for ties)."""
    def ranks(v):
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
    rx, ry = ranks(xs), ranks(ys)
    nn = len(xs)
    mx, my = sum(rx) / nn, sum(ry) / nn
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b - my) ** 2 for b in ry))
    return num / den if den else float('nan')

def auc(scores, labels):
    """Rank-based AUC of `scores` against boolean `labels` (ties get 0.5)."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    tot = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in pos for b in neg)
    return tot / (len(pos) * len(neg))

# ---------------------------------------------------------------------------
# Curve search (bucketed by lam*) -- verbatim behaviour
# ---------------------------------------------------------------------------

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


# ===========================================================================

print("=" * 78)
print("Thread 23 -- quotient reformulation of the GLV-HNP Phase-2 lattice")
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

# ===========================================================================
# EXP Q1 -- the quotient is exact: verify L cap Z*e_m = Z*n*S_D*e_m,
#           and that pi kills the trivial vector while keeping d recoverable
# ===========================================================================

print("\n" + "-" * 78)
print("EXP Q1: the kernel of the column-m deletion is exactly Z*(n*S_D*e_m)")
print("-" * 78)
print("Check per curve: (a) n*row_m - sum_i B_i*row_i == n*S_D*e_m exactly;")
print("(b) V1 = pi(V0) has exactly one zero row after LLL (rank 2m+1);")
print("(c) det(V1) == det(V0)/(n*S_D).\n")
print(f"{'curve':<18} {'m':>2} {'(a) generator':>14} {'(b) rank':>10} "
      f"{'(c) det ratio':>15}")

for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_trial, m, n, lam, p, k1, k2, 42)
    M0 = build_v0(sigs, n, lam, k1, k2)
    # (a) explicit generator of L cap Z*e_m
    dim = 2 * m + 2
    comb = [n * M0[m][j] - sum(sigs[i]['B'] * M0[i][j] for i in range(m))
            for j in range(dim)]
    expect = [0] * dim
    expect[m] = n * S_D
    gen_ok = (comb == expect)
    # (b) rank of V1
    M1 = build_v1(sigs, n, lam, k1, k2)
    r1 = lll_rows(M1)
    nzero = sum(1 for r in r1 if not any(r))
    rank_ok = (len(r1) - nzero == 2 * m + 1)
    # (c) determinant ratio (Gram determinant of the independent rows)
    B1 = independent_basis(M1)
    Mg = GSO.Mat(B1)
    Mg.update_gso()
    log_det_v1 = 0.5 * sum(math.log(Mg.get_r(i, i)) for i in range(B1.nrows))
    log_det_v0 = m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2) \
        + math.log(S_KAN)
    ratio = math.exp(log_det_v0 - log_det_v1) / (n * S_D)
    print(f"{label:<18} {m:>2} {str(gen_ok):>14} "
          f"{str(rank_ok) + ' (' + str(2*m+1) + ')':>10} {ratio:>15.6f}")

print("\n(a) True  => the trivial vector is exactly the kernel generator.")
print("(b) True  => pi(L) has rank 2m+1, i.e. exactly one dimension removed.")
print("(c) 1.000 => det(V1) = det(V0)/(n*S_D), as predicted.")

# ===========================================================================
# EXP Q2 -- norms: is the planted vector lambda_1 now?
# ===========================================================================

print("\n" + "-" * 78)
print("EXP Q2: shortest-vector / planted-vector ratio, V0 vs V1 vs V2")
print("-" * 78)
print("sv/pv < 1 means the planted vector is NOT the shortest (2026-07-29 T5).")
print("V2's 'sv' is ||target - closest||, and pv is the true (k1|k2) vector.")
print("lam_1 col = EXACT shortest vector of V1 by enumeration (not LLL).\n")
print(f"{'curve':<18} {'K1':>3} {'m':>3} {'V0 sv/pv':>9} {'V1 sv/pv':>9} "
      f"{'V2 sv/pv':>9} {'V1 lam1/pv':>11} {'V0':>4} {'V1':>4} {'V2':>4}")

q2_rows = []
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    res = trial(curve, m, d_trial, k1, 42)
    sigs = gen_signatures(G, d_trial, m, n, lam, p, k1, k2, 42)
    M1 = build_v1(sigs, n, lam, k1, k2)
    l1 = exact_svp(M1)
    pv1 = norm(planted_v1(sigs, d_trial, n, k1, k2))
    l1r = (l1 / pv1) if l1 else float('nan')
    print(f"{label:<18} {k1:>3} {m:>3} "
          f"{res['v0']['ratio']:>9.3f} {res['v1']['ratio']:>9.3f} "
          f"{res['v2']['ratio']:>9.3f} {l1r:>11.3f} "
          f"{str(res['v0']['ok']):>4} {str(res['v1']['ok']):>4} "
          f"{str(res['v2']['ok']):>4}")
    q2_rows.append((label, res, l1r))

# ===========================================================================
# EXP Q3 -- the falsifier: does the K1 wall move?
# ===========================================================================

print("\n" + "-" * 78)
print("EXP Q3 (THE FALSIFIER): K1 wall, V0 vs V1 vs V2, 5 seeds")
print("-" * 78)
print("2026-07-29 T4 baseline walls: 12-bit/2557 at K1 ~ 12-16,")
print("12-bit/2677 (lam*=0.07) at K1 ~ 4-6.  Does the reformulation move them?")
print("gamma = GH(V1)/E||planted|| is this run's closed-form predictor.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

for label, curve, _k1_hist, m in hist:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    print(f"\n{label}   n={n}  lam*={lam_star(lam, n):.4f}  m={m}  K2={k2}")
    hdr = f"  {'K1':>4} {'eff':>7} {'gamma':>7} " + \
        " ".join(f"{v:>7}" for v in ("V0", "V1", "V2"))
    print(hdr)
    for k1 in K1_GRID:
        acc = rates(curve, m, k1, SEEDS)
        g = gamma_v1(m, n, k1, k2)
        line = f"  {k1:>4} {k1*k2/n:>7.3f} {g:>7.3f} "
        line += " ".join(f"{acc[v]['win']}/{acc[v]['tot']:<5}" for v in
                         ("v0", "v1", "v2"))
        print(line)
        sys.stdout.flush()

# ===========================================================================
# EXP Q4 -- 17-bit sweep at the two regimes where the baseline degraded
# ===========================================================================

print("\n" + "-" * 78)
print("EXP Q4: 17-bit curve sweep at eff = 0.05 / 0.15 / 0.25")
print("-" * 78)
print("2026-07-29 T3 baseline: eff=0.05 -> 19/20 curves at 5/5;")
print("eff=0.15 -> 3/20; eff=0.25 -> 0/20.  Same curves, three variants.\n")

t0 = time.time()
curves17 = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
print(f"found {len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s\n")

M_SIGS = 12

def sweep(eff_target):
    print(f"eff = {eff_target}:  (per-curve wins out of {len(SEEDS)} seeds)")
    print(f"  {'p':>8} {'n':>8} {'lam*':>7} {'K1':>4} {'gamma':>7} "
          f"{'V0':>6} {'V1':>6} {'V2':>6}")
    tot = {'v0': 0, 'v1': 0, 'v2': 0}
    full = {'v0': 0, 'v1': 0, 'v2': 0}
    gam = []
    for (p, b, n, lam, G) in curves17:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_target * n / k2))
        acc = rates((p, b, n, lam, G), M_SIGS, k1, SEEDS)
        g = gamma_v1(M_SIGS, n, k1, k2)
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
        gam.append((g, acc['v1']['win']))
        Q6.append({'eff': eff_target, 'p': p, 'n': n, 'lam': lam,
                   'lam_star': lam_star(lam, n), 'K1': k1, 'K2': k2,
                   'gamma': g, 'nu_hat': nu_hat(n, lam, S_K1, S_K2),
                   'v0': acc['v0']['win'], 'v1': acc['v1']['win'],
                   'v2': acc['v2']['win']})
        for v in tot:
            tot[v] += acc[v]['win']
            full[v] += (acc[v]['win'] == len(SEEDS))
        print(f"  {p:>8} {n:>8} {lam_star(lam,n):>7.4f} {k1:>4} {g:>7.3f} "
              f"{acc['v0']['win']}/{len(SEEDS):<4} "
              f"{acc['v1']['win']}/{len(SEEDS):<4} "
              f"{acc['v2']['win']}/{len(SEEDS):<4}")
        sys.stdout.flush()
    nc = len(curves17)
    print(f"  TOTAL   curves at 5/5:  V0 {full['v0']}/{nc}  "
          f"V1 {full['v1']}/{nc}  V2 {full['v2']}/{nc}")
    print(f"          seed-level:     V0 {tot['v0']}/{nc*len(SEEDS)}  "
          f"V1 {tot['v1']}/{nc*len(SEEDS)}  V2 {tot['v2']}/{nc*len(SEEDS)}")
    return gam

gam_all = []
Q6 = []
for eff in (0.05, 0.15, 0.25):
    gam_all += sweep(eff)
    print()

# ===========================================================================
# EXP Q5 -- does gamma separate success from failure?
# ===========================================================================

print("-" * 78)
print("EXP Q5: is gamma = GH(V1)/E||planted|| a separator?")
print("-" * 78)
print("Six curve-level invariants failed between 2026-06-21 and 06-29.")
print("gamma is NOT curve-level: it depends on (m, K1, K2, n) jointly.\n")

succ = [g for g, w in gam_all if w == len(SEEDS)]
fail = [g for g, w in gam_all if w == 0]
mixed = [g for g, w in gam_all if 0 < w < len(SEEDS)]
print(f"  full-success instances (5/5): {len(succ):>3}  "
      f"gamma range [{min(succ):.3f}, {max(succ):.3f}]" if succ else
      "  full-success instances: none")
print(f"  full-failure instances (0/5): {len(fail):>3}  "
      f"gamma range [{min(fail):.3f}, {max(fail):.3f}]" if fail else
      "  full-failure instances: none")
print(f"  mixed instances:              {len(mixed):>3}  "
      f"gamma range [{min(mixed):.3f}, {max(mixed):.3f}]" if mixed else
      "  mixed instances: none")

best, best_acc = None, -1.0
allpts = [(g, w == len(SEEDS)) for g, w in gam_all]
for thr in sorted(set(g for g, _ in allpts)):
    acc = sum(1 for g, s in allpts if (g >= thr) == s) / len(allpts)
    if acc > best_acc:
        best, best_acc = thr, acc
base = max(sum(1 for _, s in allpts if s), sum(1 for _, s in allpts if not s)) \
    / len(allpts)
print(f"\n  best gamma threshold: {best:.4f}   accuracy {best_acc:.3f}   "
      f"majority baseline {base:.3f}")

print("\n" + "-" * 78)
print("EXP Q6: out-of-sample test of nu_hat, REGIME-CONTROLLED")
print("-" * 78)
print("nu_hat = lambda_1(L2)/sqrt(det L2) was validated on 2026-07-29 (run #2,")
print("commit e845207) at 20/24 bits and a single fixed regime (K1=72, m=12):")
print("AUC 0.935 against the June C1/C2 classes.  Prediction: LOW nu_hat =>")
print("EASIER recovery.  This sweep is independent data (17-bit, m=12, three")
print("eff regimes), and gamma is near-constant inside each regime, so testing")
print("nu_hat WITHIN a regime is a clean out-of-sample replication.\n")

print("AUC is on the binary 5/5 label and saturates when a regime is all-win or")
print("all-lose; spearman(nu_hat, wins) uses the graded 0..5 count instead.")
print("A NEGATIVE spearman is the predicted direction (low nu_hat = easier).\n")
print(f"  {'eff':>5} {'inst':>5} {'succ':>5} {'gamma range':>16} "
      f"{'nu_hat range':>16} {'AUC(nu_hat)':>12} {'dir':>5} {'spearman':>9}")
for eff in (0.05, 0.15, 0.25):
    rows = [r for r in Q6 if r['eff'] == eff]
    lab = [r['v1'] == len(SEEDS) for r in rows]
    nus = [r['nu_hat'] for r in rows]
    gs = [r['gamma'] for r in rows]
    # nu_hat predicts easier when LOW, so score with -nu_hat
    if not rows:
        continue
    a = auc([-x for x in nus], lab)
    if a != a:
        direction = "n/a"
    elif a > 0.5:
        direction = "ok"
    else:
        direction = "WRONG"
    grange = f"[{min(gs):.3f},{max(gs):.3f}]"
    nrange = f"[{min(nus):.3f},{max(nus):.3f}]"
    # graded score: AUC saturates when a regime is all-win or all-lose, so
    # also rank nu_hat against the raw win count (0..5).
    sp = spearman(nus, [r['v1'] for r in rows])
    print(f"  {eff:>5} {len(rows):>5} {sum(lab):>5} {grange:>16} "
          f"{nrange:>16} {a:>12.3f} {direction:>5} {sp:>9.3f}")
    sys.stdout.flush()

print("\n  Per-instance detail for the marginal regime eff=0.15")
print(f"  {'n':>8} {'lam*':>7} {'nu_hat':>7} {'gamma':>7} {'V0':>5} {'V1':>5} {'V2':>5}")
for r in sorted([r for r in Q6 if r['eff'] == 0.15], key=lambda r: r['nu_hat']):
    print(f"  {r['n']:>8} {r['lam_star']:>7.4f} {r['nu_hat']:>7.3f} "
          f"{r['gamma']:>7.3f} {r['v0']:>5} {r['v1']:>5} {r['v2']:>5}")

nus_all = [r['nu_hat'] for r in Q6]
lab_all = [r['v1'] == len(SEEDS) for r in Q6]
print(f"\n  pooled across all 3 regimes: AUC(-nu_hat) = "
      f"{auc([-x for x in nus_all], lab_all):.3f}   "
      f"AUC(gamma) = {auc([r['gamma'] for r in Q6], lab_all):.3f}")
print("  (pooled AUC(gamma) is inflated: gamma is essentially a monotone")
print("   function of eff, and the three regimes were chosen far apart.)")

print("\n" + "-" * 78)
print("CVP oracle exactness audit")
print("-" * 78)
print("The planted difference is always a legal candidate, so an exact CVP")
print("oracle can never return a strictly longer difference vector.")
print(f"  CVP calls: {CVP_AUDIT['n']}   "
      f"calls returning ||diff|| > ||planted||: {CVP_AUDIT['violations']}")
print("  0 violations => fplll's closest_vector behaved as an exact oracle "
      "on every instance.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
