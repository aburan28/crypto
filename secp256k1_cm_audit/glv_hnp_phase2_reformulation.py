"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector
is actually the target of the reduction.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, Exp T5):
  The Phase-2 lattice (dim 2m+2, `build_glv_lattice` in
  glv_hnp_phase2_20bit.py:263) always contains the "parasitic" vector

        n * S_D * e_m          (norm exactly n*S_D = n, since S_D = 1)

  obtained as n*row_d - sum_i B_i*row_i.  The planted vector has norm
  ~ n*sqrt(2m/3 + 4/3) > n for every m >= 1, so the planted vector is NEVER
  lambda_1.  T5 measured |sv|/|pv| in [0.34, 0.61] on every curve, success
  and failure alike.  The parasite carries no information (d is only defined
  mod n) and no choice of S_D removes it (both vectors scale linearly in S_D).

  The 2026-07-29 next-step proposal was therefore: quotient out the trivial
  e_m direction, or replace the Kannan embedding by an explicit CVP call, and
  see whether the K1 wall of Exp T4 moves outward.

Three formulations are compared here on identical signature sets:

  F0  BASELINE, dim 2m+2.  Verbatim `build_glv_lattice`.  Kannan embedding;
      d carried as its own lattice coordinate (scale S_D = 1).
      planted = (k1_i*S_K1 | d*S_D | k2_i*S_K2 | S_KANNAN)

  F1  NO-D + KANNAN, dim 2m+1.  The d column is deleted.  d survives only as
      the integer coefficient of the row (B_0*S_K1, ..., B_{m-1}*S_K1), and is
      reconstructed arithmetically after reduction from the recovered
      (k1_0, k2_0) via  d = (k1_0 - A_0 + lam*k2_0) * B_0^-1 mod n.
      planted = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN),  norm ~ n*sqrt(2m/3 + 1)
      The parasite n*S_D*e_m does not exist in this lattice.

  F2  NO-D + BABAI CVP, dim 2m.  Same as F1 with the Kannan row deleted too;
      instead LLL-reduce and run Babai nearest-plane against the target
            t = (-A_0*S_K1, ..., -A_{m-1}*S_K1 | 0, ..., 0)
      so that  v - t = (k1_i*S_K1 | k2_i*S_K2),  norm ~ n*sqrt(2m/3).
      This is the Nguyen-Shparlinski posture (CVP, not Kannan-SVP), which is
      immune to the parasite by construction: BDD does not care that the
      lattice contains short vectors in the quotiented direction.

Determinants (all three lattices, exact):
  det F0 = n^(3m+1) / (K1*K2)^m        dim 2m+2
  det F1 = n^(3m)   / (K1*K2)^m        dim 2m+1
  det F2 = n^(3m-1) / (K1*K2)^m        dim 2m
(using S_K1 = n//K1, S_K2 = n//K2, S_D = 1, S_KANNAN = n, and the fact that
 the sublattice n*Z^m + Z*B of Z^m has index n^(m-1) when some gcd(B_i,n)=1.)

Hence for F2 the Gaussian-heuristic ratio of the planted error to GH(L) is

    R = ||v|| / GH(L) = sqrt(2*pi*e/3) * sqrt(eff) * n^(1/(2m)),
        eff = K1*K2/n

which is the central quantitative prediction tested by Exp C/D below:
  * the n^(1/(2m)) factor means a toy curve needs LARGE m before the
    asymptotic eff-only criterion bites -- this retro-explains Exp T4b
    (2026-07-29), where m = 8..32 at eff = 0.157 did NOT rescue the attack;
  * asymptotically in m the criterion is eff < 3/(2*pi*e) = 0.1757,
    independent of the curve.

Experiments:
  A  harness reproduction: F0 on the Exp T4 K1-grid (must match 2026-07-29)
  B  F0 vs F1 vs F2 on the same grid; where does each wall sit?
  C  m-sweep at the Exp T4b point (K1 = 8, eff = 0.157), m up to 64.
     Falsifier for "the wall is information-theoretic".
  D  (eff, m) grid: is the empirical critical R* a constant?
  E  cross-size check on a 17-bit and a 20-bit curve.

Run: python3 glv_hnp_phase2_reformulation.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, GSO, FPLLL

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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p

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
                return (p, b, n, glv_roots(n)[0], G)
    return None

def find_curve_in_range(lo, hi, seed=12345):
    """First j=0 GLV curve with p in [lo,hi), n prime, n = 1 mod 3."""
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    if glv_roots(n_cand) is None:
                        continue
                    cur = build_curve(p, n_cand, seed=seed)
                    if cur is not None:
                        return cur
        p = int(sympy.nextprime(p))
    return None

# ---------------------------------------------------------------------------
# Signatures and scales (identical to the 2026-07-29 harness)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- verbatim from the Phase-2 harness."""
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

def inorm2(v):
    return sum(x * x for x in v)

# ---------------------------------------------------------------------------
# F0 -- baseline, dim 2m+2, Kannan, d as a lattice coordinate
# ---------------------------------------------------------------------------

def build_F0(sigs, n, lam, k1_bound, k2_bound):
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

def planted_F0(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_F0(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# F1 -- no d column, dim 2m+1, Kannan
# ---------------------------------------------------------------------------

def build_F1(sigs, n, lam, k1_bound, k2_bound):
    """Rows (2m+2 generators of a rank-(2m+1) lattice in 2m+1 columns).
    Columns: 0..m-1 = k1, m..2m-1 = k2, 2m = Kannan."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # n*S_K1*e_i
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                            # the d row (no d column now)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for j in range(m):                       # k2 rows
        r = [0] * dim
        r[j] = -lam * S_K1
        r[m + j] = S_K2
        rows.append(r)
    r = [0] * dim                            # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KAN
    rows.append(r)
    return rows

def planted_F1(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def d_from_k(k1_0, k2_0, sigs, n, lam):
    """d = (k1_0 - A_0 + lam*k2_0) * B_0^-1 mod n."""
    B0 = sigs[0]['B'] % n
    if math.gcd(B0, n) != 1:
        return None
    return (k1_0 - sigs[0]['A'] + lam * k2_0) * modinv(B0, n) % n

def recover_F1(rows, sigs, m, n, lam, k1_bound, k2_bound, d_secret):
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        r = row if last > 0 else [-x for x in row]
        if r[0] % S_K1 or r[m] % S_K2:
            continue
        d_cand = d_from_k(r[0] // S_K1, r[m] // S_K2, sigs, n, lam)
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# F2 -- no d column, no Kannan row, dim 2m, Babai nearest-plane CVP
# ---------------------------------------------------------------------------

def build_F2(sigs, n, lam, k1_bound, k2_bound):
    """Rows (2m+1 generators of a rank-2m lattice in 2m columns).
    Columns: 0..m-1 = k1, m..2m-1 = k2."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for j in range(m):
        r = [0] * dim
        r[j] = -lam * S_K1
        r[m + j] = S_K2
        rows.append(r)
    return rows

def target_F2(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, _S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m

def error_F2(sigs, n, k1_bound, k2_bound):
    """The planted CVP error v - t = (k1_i*S_K1 | k2_i*S_K2)."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    return ([sigs[i]['k1'] * S_K1 for i in range(m)]
            + [sigs[i]['k2'] * S_K2 for i in range(m)])

# ---------------------------------------------------------------------------
# Lattice plumbing
# ---------------------------------------------------------------------------

def lll_rows(rows, use_mpfr=False):
    """LLL-reduce, drop zero rows, return (list-of-rows, IntegerMatrix)."""
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    out = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    out = [r for r in out if any(r)]
    return out

def babai(rows, target, float_type='d'):
    """Babai nearest-plane against `target`; returns the closest lattice
    vector found.  `rows` must already be LLL-reduced and zero-row-free."""
    B = IntegerMatrix.from_matrix(rows)
    M = GSO.Mat(B, float_type=float_type)
    M.update_gso()
    c = M.babai([int(x) for x in target])
    v = [0] * B.ncols
    for i, ci in enumerate(c):
        if ci:
            for j in range(B.ncols):
                v[j] += ci * B[i][j]
    return v

def gaussian_heuristic_log(logdet, dim):
    """log GH = log(det)/dim + 0.5*log(dim/(2*pi*e))."""
    return logdet / dim + 0.5 * math.log(dim / (2 * math.pi * math.e))

def logdet_F(kind, m, n, k1_bound, k2_bound):
    """Exact log-determinants (see module docstring)."""
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    base = m * math.log(S_K1) + m * math.log(S_K2) + (m - 1) * math.log(n)
    if kind == 'F2':
        return base                              # dim 2m
    if kind == 'F1':
        return base + math.log(n)                # dim 2m+1  (Kannan col)
    if kind == 'F0':
        return base + 2 * math.log(n)            # dim 2m+2  (d col + Kannan)
    raise ValueError(kind)

# ---------------------------------------------------------------------------
# One trial, all three formulations, on the SAME signature set
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, kinds=('F0', 'F1', 'F2'),
          float_type='d'):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    res = {}

    if 'F0' in kinds:
        rows = lll_rows(build_F0(sigs, n, lam, k1_bound, k2_bound))
        pv = planted_F0(sigs, d_secret, n, k1_bound, k2_bound)
        res['F0'] = {
            'ok': recover_F0(rows, m, n, S_KAN, d_secret) is not None,
            'pv': norm(pv),
            'b1': min(norm(r) for r in rows),
        }

    if 'F1' in kinds:
        rows = lll_rows(build_F1(sigs, n, lam, k1_bound, k2_bound))
        pv = planted_F1(sigs, n, k1_bound, k2_bound)
        res['F1'] = {
            'ok': recover_F1(rows, sigs, m, n, lam,
                             k1_bound, k2_bound, d_secret) is not None,
            'pv': norm(pv),
            'b1': min(norm(r) for r in rows),
        }

    if 'F2' in kinds:
        rows = lll_rows(build_F2(sigs, n, lam, k1_bound, k2_bound))
        t = target_F2(sigs, n, k1_bound, k2_bound)
        v = babai(rows, t, float_type=float_type)
        e = [v[j] - t[j] for j in range(2 * m)]
        ok = False
        if e[0] % S_K1 == 0 and e[m] % S_K2 == 0:
            d_cand = d_from_k(e[0] // S_K1, e[m] // S_K2, sigs, n, lam)
            ok = (d_cand is not None and d_cand == d_secret)
        res['F2'] = {
            'ok': ok,
            'pv': norm(error_F2(sigs, n, k1_bound, k2_bound)),
            'b1': min(norm(r) for r in rows),
            'err': norm(e),
        }

    return res

def rate(curve, m, k1_bound, seeds, kinds=('F0', 'F1', 'F2'), float_type='d'):
    p, b, n, lam, G = curve
    wins = {k: 0 for k in kinds}
    acc = {k: [] for k in kinds}
    tot = 0
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = trial(curve, m, d_trial, k1_bound, seed, kinds, float_type)
        if r is None:
            continue
        tot += 1
        for k in kinds:
            wins[k] += bool(r[k]['ok'])
            acc[k].append(r[k]['pv'] / r[k]['b1'])
    mean = {k: (sum(acc[k]) / len(acc[k]) if acc[k] else float('nan'))
            for k in kinds}
    return wins, tot, mean

def R_pred(kind, m, n, k1_bound, k2_bound):
    """||planted|| / GH(L), the Gaussian-heuristic ratio."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if kind == 'F0':
        dim = 2 * m + 2
        pv2 = (m * (k1_bound * S_K1) ** 2 / 3.0 + (n * S_D) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0 + SKAN2(S_KAN))
    elif kind == 'F1':
        dim = 2 * m + 1
        pv2 = (m * (k1_bound * S_K1) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0 + SKAN2(S_KAN))
    else:
        dim = 2 * m
        pv2 = (m * (k1_bound * S_K1) ** 2 / 3.0
               + m * (k2_bound * S_K2) ** 2 / 3.0)
    lgh = gaussian_heuristic_log(logdet_F(kind, m, n, k1_bound, k2_bound), dim)
    return math.exp(0.5 * math.log(pv2) - lgh)

def SKAN2(S_KAN):
    return float(S_KAN) ** 2


# ===========================================================================
print("=" * 78)
print("Thread 23 -- reformulating the GLV-HNP Phase-2 lattice (F0 / F1 / F2)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),   # lam* = 0.340
    ("12-bit/2677", 2677, 2, 2647, 185),    # lam* = 0.070  (the "FAIL" curve)
]
curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G)))
    print(f"  {label}: n={n} ({n.bit_length()}b) lam={lam} "
          f"lam*={lam_star(lam,n):.3f} K2={math.isqrt(n)+1}")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_BASE = {"12-bit/2557": 8, "12-bit/2677": 10}   # the Exp T4 values

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP A -- harness reproduction: F0 on the Exp T4 K1 grid")
print("-" * 78)
print("Must match the 2026-07-29 T4 table exactly (5 seeds, same d-derivation).")
print("  expected 2557 (m=8):  5 5 5 5 5 4 1 0")
print("  expected 2677 (m=10): 5 5 5 2 0 0 0 0\n")
hdr = f"{'curve':<14} {'m':>3} " + " ".join(f"{'K1=' + str(k):>6}" for k in K1_GRID)
print(hdr)
expA = {}
for label, cur in curves:
    m = M_BASE[label]
    row = []
    for k1 in K1_GRID:
        w, t, _ = rate(cur, m, k1, SEEDS, kinds=('F0',))
        row.append(w['F0'])
    expA[label] = row
    print(f"{label:<14} {m:>3} " + " ".join(f"{v:>6}" for v in row))

print("\nmatch 2557:", expA["12-bit/2557"] == [5, 5, 5, 5, 5, 4, 1, 0])
print("match 2677:", expA["12-bit/2677"] == [5, 5, 5, 2, 0, 0, 0, 0])

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP B -- F0 vs F1 vs F2 on the same signature sets, same K1 grid")
print("-" * 78)
print("F1 = drop the d column (kills the parasite n*S_D*e_m).")
print("F2 = drop d column AND Kannan row; Babai nearest-plane CVP instead.\n")
for label, cur in curves:
    p, b, n, lam, G = cur
    m = M_BASE[label]
    k2b = math.isqrt(n) + 1
    print(f"{label}  m={m}  K2={k2b}")
    print(f"  {'K1':>4} {'eff':>7} | " + " ".join(f"{k:>5}" for k in ('F0', 'F1', 'F2'))
          + " | " + " ".join(f"{'R_' + k:>6}" for k in ('F0', 'F1', 'F2'))
          + " | " + " ".join(f"{'pv/b1_' + k:>10}" for k in ('F0', 'F1', 'F2')))
    for k1 in K1_GRID:
        w, t, mean = rate(cur, m, k1, SEEDS)
        eff = k1 * k2b / n
        Rs = [R_pred(k, m, n, k1, k2b) for k in ('F0', 'F1', 'F2')]
        print(f"  {k1:>4} {eff:>7.4f} | "
              + " ".join(f"{w[k]:>5}" for k in ('F0', 'F1', 'F2'))
              + " | " + " ".join(f"{r:>6.3f}" for r in Rs)
              + " | " + " ".join(f"{mean[k]:>10.3f}" for k in ('F0', 'F1', 'F2')))
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP C -- the Exp T4b point (2677, K1=8, eff=0.157): does MORE DATA")
print("         rescue F1/F2 where it did not rescue F0?")
print("-" * 78)
print("2026-07-29 T4b measured F0 at m = 8/12/16/24/32 -> 0,0,1,0,1 of 5.")
print("Prediction from R = sqrt(2*pi*e/3)*sqrt(eff)*n^(1/2m): R falls with m,")
print("so if a critical R* exists there is an m at which the wall breaks.\n")
label, cur = curves[1]
p, b, n, lam, G = cur
k2b = math.isqrt(n) + 1
print(f"{'m':>4} {'dim F2':>7} | {'F0':>4} {'F1':>4} {'F2':>4} | "
      f"{'R_F0':>6} {'R_F1':>6} {'R_F2':>6}")
for m in [10, 14, 18, 24, 32, 40, 48, 64]:
    w, t, mean = rate(cur, m, 8, SEEDS)
    Rs = [R_pred(k, m, n, 8, k2b) for k in ('F0', 'F1', 'F2')]
    print(f"{m:>4} {2*m:>7} | " + " ".join(f"{w[k]:>4}" for k in ('F0', 'F1', 'F2'))
          + " | " + " ".join(f"{r:>6.3f}" for r in Rs))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP D -- (eff, m) grid on 2677: is the critical R* a constant?")
print("-" * 78)
print("Each cell = F2 wins / 5.  R_F2 printed underneath.\n")
M_GRID = [8, 12, 16, 24, 32, 48]
K1_D = [4, 6, 8, 12, 16]
print(f"{'K1':>4} {'eff':>7} | " + " ".join(f"{'m=' + str(mm):>11}" for mm in M_GRID))
grid = {}
for k1 in K1_D:
    cells = []
    for mm in M_GRID:
        w, t, mean = rate(cur, mm, k1, SEEDS, kinds=('F2',))
        R = R_pred('F2', mm, n, k1, k2b)
        grid[(k1, mm)] = (w['F2'], R)
        cells.append(f"{w['F2']}/5 R={R:.2f}")
    print(f"{k1:>4} {k1*k2b/n:>7.4f} | " + " ".join(f"{c:>11}" for c in cells))

succ = [R for (k1, mm), (wv, R) in grid.items() if wv >= 4]
fail = [R for (k1, mm), (wv, R) in grid.items() if wv <= 1]
if succ and fail:
    print(f"\n  R over cells with >=4/5 wins : [{min(succ):.3f}, {max(succ):.3f}]")
    print(f"  R over cells with <=1/5 wins : [{min(fail):.3f}, {max(fail):.3f}]")
    print(f"  separating?  {max(succ) < min(fail)}   "
          f"(gap {min(fail) - max(succ):+.3f})")
else:
    print("\n  not enough contrast in the grid to bracket R*")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E -- cross-size check (17-bit and 20-bit j=0 GLV curves)")
print("-" * 78)
big = []
for lo, hi, tag in [(2**16, 2**17, "17-bit"), (2**19, 2**20, "20-bit")]:
    c = find_curve_in_range(lo, hi)
    if c is None:
        print(f"  {tag}: no curve found")
        continue
    p2, b2, n2, lam2, G2 = c
    big.append((tag, c))
    print(f"  {tag}: p={p2} n={n2} ({n2.bit_length()}b) lam={lam2} "
          f"lam*={lam_star(lam2,n2):.3f} K2={math.isqrt(n2)+1}")

FT = 'd'
for tag, c in big:
    p2, b2, n2, lam2, G2 = c
    k2b2 = math.isqrt(n2) + 1
    ft = 'd' if n2.bit_length() <= 18 else 'mpfr'
    if ft == 'mpfr':
        FPLLL.set_precision(212)
    print(f"\n{tag}  n={n2}  K2={k2b2}  (Babai float_type={ft})")
    print(f"  {'K1':>4} {'m':>4} {'eff':>7} | {'F0':>4} {'F1':>4} {'F2':>4} | "
          f"{'R_F2':>6}")
    for k1, mm in [(4, 12), (8, 12), (8, 24), (12, 24), (16, 32), (8, 40)]:
        w, t, mean = rate(c, mm, k1, SEEDS[:3], float_type=ft)
        R = R_pred('F2', mm, n2, k1, k2b2)
        print(f"  {k1:>4} {mm:>4} {k1*k2b2/n2:>7.4f} | "
              + " ".join(f"{w[k]:>4}" for k in ('F0', 'F1', 'F2'))
              + f" | {R:>6.3f}")

print("\nDone.")
