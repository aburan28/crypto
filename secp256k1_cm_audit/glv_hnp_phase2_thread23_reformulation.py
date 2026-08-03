"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector
is actually the thing the reduction is looking for.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, experiment T5):
  In the Phase-2 Kannan lattice of glv_hnp_phase2_20bit.py:263 the shortest
  vector is ALWAYS the trivial vector n*S_D*e_m (100% of its energy in the d
  column, |sv[m]|/n = 1.0000 exactly).  It is 1.6-3x shorter than the planted
  vector on every curve tested, success and failure alike, so the planted
  vector is never lambda_1 and recovery is a BDD/coset condition rather than
  an SVP condition.  The 2026-07-29 entry proposed:

     "project the lattice along e_m (quotient out the trivial n*e_m
      direction) and solve BDD in the projection, or replace the Kannan
      embedding with an explicit CVP call (Babai nearest-plane on the
      reduced basis)."

  Falsifier stated there: if sv/pv rises above 1 after the reformulation AND
  the K1 wall on the lam*=0.07 curve (currently K1 ~ 4-6) moves outward, the
  reformulation is a real improvement.  If the wall stays put, the wall is
  information-theoretic and Phase 2 is at its ceiling.

Four lattice variants are compared on the identical signature sets:

  V0  BASELINE.  Kannan embedding, dim 2m+2, LLL, scan reduced rows for
      |last| = S_KANNAN.  Verbatim from glv_hnp_phase2_20bit.py.

  V1  CVP.  Same lattice with the Kannan row deleted (dim 2m+1); recover d
      by Babai nearest-plane against the target t = (-A_i*S_K1, 0, ..., 0).
      The trivial vector n*S_D*e_m survives here, but it is harmless: it is
      exactly the "d is only defined mod n" ambiguity, and Babai simply picks
      the representative with |d| <= n/2, which still gives the right d mod n.

  V2  D-ELIMINATED + Kannan.  d is removed from the system algebraically
      before building the lattice.  With t_i = B_i * B_0^{-1} mod n and
      C_i = A_i - t_i*A_0 mod n, the m congruences
          k1_i + lam*k2_i - A_i = B_i*d   (mod n)
      become m-1 congruences with no d at all:
          k1_i + lam*k2_i - t_i*k1_0 - t_i*lam*k2_0 = C_i   (mod n),  i>=1.
      Unknowns are the 2m small values k1_*, k2_*; dim is 2m+1 (2m coords +
      Kannan).  The d column - and with it the trivial vector n*S_D*e_m - is
      gone.  The shortest trivial vectors that remain are n*S_K1*e_{x0} and
      n*S_K2*e_{yj}, of norm n^2/K1 and n^2/K2, i.e. a factor ~n/K1 LONGER
      than the planted vector instead of ~2x shorter.

  V3  D-ELIMINATED + CVP.  Same elimination, Kannan row deleted (dim 2m),
      Babai nearest-plane against t = (0, -C_i*S_K1, ..., 0).

Experiments:
  E1  correctness of every variant on the known-good 8-bit curve
  E2  sv/pv (shortest-reduced-vector / planted-vector norm) per variant:
      does the planted vector become lambda_1?
  E3  the T4 K1-wall grid, all four variants, both historical 12-bit curves
  E4  the 2026-07-29 T4b control: does more data rescue K1=8 under the new
      formulations?

Run: python3 glv_hnp_phase2_thread23_reformulation.py
"""

import math
import random
from fractions import Fraction

from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Arithmetic helpers (verbatim from glv_hnp_phase2_lambda_threshold.py:87)
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

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Signature generation (verbatim from glv_hnp_phase2_lambda_threshold.py:230)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) — the 2026-06-15 column-diagonal scaling."""
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
# V0 / V1 — the original (d-carrying) lattice
# ---------------------------------------------------------------------------
# Columns: [0 .. m-1] = k1_i,  [m] = d,  [m+1 .. 2m] = k2_i,  [2m+1] = Kannan.
# Planted vector: (k1_i*S_K1, d*S_D, k2_i*S_K2, S_KANNAN).

def build_lattice_v0(sigs, n, lam, k1_bound, k2_bound):
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

def build_lattice_v1(sigs, n, lam, k1_bound, k2_bound):
    """V0 with the Kannan row and column deleted: dim 2m+1, CVP formulation."""
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

def target_v1(sigs, n, k1_bound, k2_bound):
    """t with v - t = (k1_i*S_K1, d*S_D, k2_i*S_K2) for the correct v in L."""
    m = len(sigs)
    S_K1 = scales(n, k1_bound, k2_bound)[0]
    t = [0] * (2 * m + 1)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return t

# ---------------------------------------------------------------------------
# V2 / V3 — the d-eliminated lattice
# ---------------------------------------------------------------------------
# t_i = B_i/B_0 mod n,  C_i = A_i - t_i*A_0 mod n   (i = 1..m-1)
#   k1_i + lam*k2_i - t_i*k1_0 - t_i*lam*k2_0 = C_i  (mod n)
# Columns: [0 .. m-1] = k1_i (x), [m .. 2m-1] = k2_i (y), [2m] = Kannan.

def eliminate_d(sigs, n):
    """Returns (t_i, C_i) for i = 1..m-1 (index 0 of the lists is i=1)."""
    m = len(sigs)
    B0inv = modinv(sigs[0]['B'], n)
    ts, Cs = [], []
    for i in range(1, m):
        t_i = sigs[i]['B'] * B0inv % n
        C_i = (sigs[i]['A'] - t_i * sigs[0]['A']) % n
        ts.append(t_i)
        Cs.append(C_i)
    return ts, Cs

def build_lattice_v2(sigs, n, lam, k1_bound, k2_bound, kannan=True):
    """d-eliminated lattice.  kannan=True -> dim 2m+1 (V2); False -> 2m (V3)."""
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ts, Cs = eliminate_d(sigs, n)
    ncol = 2 * m + (1 if kannan else 0)
    rows = []
    # mod-n rows for the determined coordinates x_1 .. x_{m-1}
    for i in range(1, m):
        r = [0] * ncol
        r[i] = n * S_K1
        rows.append(r)
    # x_0 row
    r = [0] * ncol
    r[0] = S_K1
    for i in range(1, m):
        r[i] = ts[i - 1] * S_K1
    rows.append(r)
    # y_0 row  (k2_0 enters every determined coordinate through t_i*lam)
    r = [0] * ncol
    r[m] = S_K2
    for i in range(1, m):
        r[i] = (ts[i - 1] * lam % n) * S_K1
    rows.append(r)
    # y_j rows, j = 1 .. m-1
    for j in range(1, m):
        r = [0] * ncol
        r[m + j] = S_K2
        r[j] = -(lam % n) * S_K1
        rows.append(r)
    if kannan:
        r = [0] * ncol
        for i in range(1, m):
            r[i] = Cs[i - 1] * S_K1
        r[2 * m] = S_KAN
        rows.append(r)
    assert len(rows) == ncol, (len(rows), ncol)
    return rows

def target_v3(sigs, n, k1_bound, k2_bound):
    """t with v - t = (k1_i*S_K1, k2_j*S_K2) for the correct v in the V3 lattice."""
    m = len(sigs)
    S_K1 = scales(n, k1_bound, k2_bound)[0]
    _, Cs = eliminate_d(sigs, n)
    t = [0] * (2 * m)
    for i in range(1, m):
        t[i] = -Cs[i - 1] * S_K1
    return t

def planted_v2(sigs, n, k1_bound, k2_bound, kannan=True):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ncol = 2 * m + (1 if kannan else 0)
    v = [0] * ncol
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    if kannan:
        v[2 * m] = S_KAN
    return v

def d_from_k(k1_0, k2_0, sigs, n, lam):
    """Invert k1_0 + lam*k2_0 = A_0 + B_0*d (mod n) for d."""
    lhs = (k1_0 + lam * k2_0 - sigs[0]['A']) % n
    return lhs * modinv(sigs[0]['B'], n) % n

# ---------------------------------------------------------------------------
# Babai nearest-plane (exact rational GSO; dims here are <= 25)
# ---------------------------------------------------------------------------

def babai_nearest_plane(basis, target):
    """Exact-arithmetic nearest plane against an already-reduced `basis`."""
    dim = len(basis)
    B = [[Fraction(x) for x in row] for row in basis]
    # Gram-Schmidt
    Bs, mu = [], [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        v = B[i][:]
        for j in range(i):
            nj = sum(x * x for x in Bs[j])
            if nj == 0:
                mu[i][j] = Fraction(0)
                continue
            mu[i][j] = sum(B[i][k] * Bs[j][k] for k in range(len(v))) / nj
            v = [v[k] - mu[i][j] * Bs[j][k] for k in range(len(v))]
        Bs.append(v)
    w = [Fraction(x) for x in target]
    coeffs = [0] * dim
    for i in range(dim - 1, -1, -1):
        nj = sum(x * x for x in Bs[i])
        if nj == 0:
            continue
        c = sum(w[k] * Bs[i][k] for k in range(len(w))) / nj
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = ci
        w = [w[k] - ci * Fraction(basis[i][k]) for k in range(len(w))]
    v = [sum(coeffs[i] * basis[i][k] for i in range(dim)) for k in range(len(target))]
    return v

def lll_reduce(M):
    dim = len(M)
    ncol = len(M[0])
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(ncol)] for i in range(dim)]

# ---------------------------------------------------------------------------
# The four attacks
# ---------------------------------------------------------------------------

def attack_v0(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    red = lll_reduce(build_lattice_v0(sigs, n, lam, k1_bound, k2_bound))
    pv = norm(planted_v0(sigs, d_secret, n, k1_bound, k2_bound))
    sv = min(norm(r) for r in red)
    for row in red:
        if abs(row[2 * m + 1]) != S_KAN:
            continue
        sign = 1 if row[2 * m + 1] > 0 else -1
        if (sign * row[m]) % n == d_secret:
            return True, pv, sv
    return False, pv, sv

def attack_v1(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    red = lll_reduce(build_lattice_v1(sigs, n, lam, k1_bound, k2_bound))
    t = target_v1(sigs, n, k1_bound, k2_bound)
    v = babai_nearest_plane(red, t)
    w = [v[k] - t[k] for k in range(len(t))]
    pv = norm([s['k1'] * S_K1 for s in sigs] + [d_secret * S_D]
              + [s['k2'] * S_K2 for s in sigs])
    sv = min(norm(r) for r in red)
    return (w[m] // S_D) % n == d_secret, pv, sv, norm(w)

def attack_v2(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    red = lll_reduce(build_lattice_v2(sigs, n, lam, k1_bound, k2_bound, kannan=True))
    pv = norm(planted_v2(sigs, n, k1_bound, k2_bound, kannan=True))
    sv = min(norm(r) for r in red)
    for row in red:
        if abs(row[2 * m]) != S_KAN:
            continue
        sign = 1 if row[2 * m] > 0 else -1
        k1_0, k2_0 = sign * row[0], sign * row[m]
        if k1_0 % S_K1 or k2_0 % S_K2:
            continue
        if d_from_k(k1_0 // S_K1, k2_0 // S_K2, sigs, n, lam) == d_secret:
            return True, pv, sv
    return False, pv, sv

def attack_v3(sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    red = lll_reduce(build_lattice_v2(sigs, n, lam, k1_bound, k2_bound, kannan=False))
    t = target_v3(sigs, n, k1_bound, k2_bound)
    v = babai_nearest_plane(red, t)
    w = [v[k] - t[k] for k in range(len(t))]
    pv = norm(planted_v2(sigs, n, k1_bound, k2_bound, kannan=False))
    sv = min(norm(r) for r in red)
    if w[0] % S_K1 or w[m] % S_K2:
        return False, pv, sv, norm(w)
    ok = d_from_k(w[0] // S_K1, w[m] // S_K2, sigs, n, lam) == d_secret
    return ok, pv, sv, norm(w)

VARIANTS = ['V0', 'V1', 'V2', 'V3']

def run_one(variant, curve, m, k1_bound, seed):
    """Returns dict(ok, pv, sv, ...) or None if signature generation failed."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if variant == 'V0':
        ok, pv, sv = attack_v0(sigs, n, lam, k1_bound, k2_bound, d_secret)
        wn = None
    elif variant == 'V1':
        ok, pv, sv, wn = attack_v1(sigs, n, lam, k1_bound, k2_bound, d_secret)
    elif variant == 'V2':
        ok, pv, sv = attack_v2(sigs, n, lam, k1_bound, k2_bound, d_secret)
        wn = None
    elif variant == 'V3':
        ok, pv, sv, wn = attack_v3(sigs, n, lam, k1_bound, k2_bound, d_secret)
    else:
        raise ValueError(variant)
    return {'ok': ok, 'pv': pv, 'sv': sv, 'wn': wn, 'ratio': sv / pv if pv else float('nan')}

def rate(variant, curve, m, k1_bound, seeds):
    wins, ratios = 0, []
    for s in seeds:
        r = run_one(variant, curve, m, k1_bound, s)
        if r is None:
            continue
        wins += bool(r['ok'])
        ratios.append(r['ratio'])
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))

# ===========================================================================
# Driver
# ===========================================================================

def main():

    print("=" * 78)
    print("Thread 23 — reformulating the Phase-2 lattice (GLV-HNP)")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]

    HIST = [
        # label,             p,    b, n,    lam,  K1, m
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

    print("\nHistorical curves (lam* = min(lam, n-lam)/n, the representation-")
    print("invariant size of the eigenvalue; see 2026-07-29 T1):")
    print(f"  {'curve':<20} {'p':>6} {'n':>6} {'lam':>6} {'lam*':>7} {'K1':>4} {'m':>3}")
    for label, (p, b, n, lam, G), k1, m in curves:
        print(f"  {label:<20} {p:>6} {n:>6} {lam:>6} {lam_star(lam, n):>7.4f} {k1:>4} {m:>3}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E1: correctness — every variant must recover d in the easy regime")
    print("-" * 78)
    print("If a reformulation is wrong (bad target, bad elimination) it fails here.")
    print("Curve 8-bit/199, K1=2, m=6; and 12-bit/2557 at K1=2, m=8.\n")

    easy = [(curves[0][0], curves[0][1], 6, 2), (curves[1][0], curves[1][1], 8, 2)]
    print(f"  {'curve':<20} {'m':>3} {'K1':>3} " + " ".join(f"{v:>7}" for v in VARIANTS))
    e1_ok = True
    for label, curve, m, k1 in easy:
        cells = []
        for v in VARIANTS:
            w, tot, _ = rate(v, curve, m, k1, SEEDS)
            cells.append(f"{w}/{tot}")
            if w != tot:
                e1_ok = False
        print(f"  {label:<20} {m:>3} {k1:>3} " + " ".join(f"{c:>7}" for c in cells))
    print(f"\n  E1 verdict: {'PASS — all four formulations are correct' if e1_ok else 'FAIL — a variant is miswired'}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E2: sv/pv — is the planted vector lambda_1 now?")
    print("-" * 78)
    print("sv = norm of the shortest row after LLL, pv = norm of the planted vector.")
    print("2026-07-29 T5 measured sv/pv in [0.34, 0.61] for V0 (planted is NEVER")
    print("shortest).  The Thread 23 falsifier asks whether sv/pv rises above 1.\n")

    print(f"  {'curve':<20} {'K1':>3} " + " ".join(f"{v+' sv/pv':>10}" for v in VARIANTS))
    e2_rows = []
    for label, curve, k1, m in curves:
        cells = []
        for v in VARIANTS:
            _, _, ratio = rate(v, curve, m, k1, SEEDS)
            cells.append(ratio)
        e2_rows.append((label, k1, cells))
        print(f"  {label:<20} {k1:>3} " + " ".join(f"{c:>10.3f}" for c in cells))

    print("\n  Note: for V0/V2 the planted vector is a lattice vector and sv/pv < 1")
    print("  means it is not shortest.  For V1/V3 (CVP) pv is the target distance,")
    print("  so sv/pv > 1 is the BDD-friendly regime (dist < lambda_1).")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E3: the K1 wall (2026-07-29 T4 grid) under all four formulations")
    print("-" * 78)
    print("T4 found the wall at K1 ~ 12-16 for lam*=0.34 and K1 ~ 4-6 for lam*=0.07.")
    print("Falsifier: a real reformulation moves the wall outward (to larger K1).\n")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    e3 = {}
    for label, curve, _k1, m in curves[1:]:      # the two 12-bit curves
        print(f"  {label}  (m={m}, {len(SEEDS)} seeds)")
        hdr = "    " + f"{'variant':<8}" + " ".join(f"{'K1='+str(k):>7}" for k in K1_GRID)
        print(hdr)
        for v in VARIANTS:
            cells = []
            for k1 in K1_GRID:
                w, tot, _ = rate(v, curve, m, k1, SEEDS)
                cells.append(w)
            e3[(label, v)] = cells
            print(f"    {v:<8}" + " ".join(f"{str(c)+'/'+str(len(SEEDS)):>7}" for c in cells))
        # wall = largest K1 with a full sweep
        for v in VARIANTS:
            cells = e3[(label, v)]
            wall = max([K1_GRID[i] for i, c in enumerate(cells) if c == len(SEEDS)], default=0)
            print(f"    -> {v} wall (largest K1 with {len(SEEDS)}/{len(SEEDS)}): K1 = {wall}")
        print()

    # ---------------------------------------------------------------------------
    print("-" * 78)
    print("EXP E4: T4b control — does more data rescue K1=8 on the lam*=0.07 curve?")
    print("-" * 78)
    print("2026-07-29 T4b: under V0, m = 8/12/16/24/32 gave 0,0,1,0,1 out of 5.")
    print("If a reformulation is information-theoretically better, more data should")
    print("now help.\n")

    label, curve, _k1, _m = curves[2]
    M_GRID = [8, 12, 16, 24, 32]
    print(f"  {label}, K1=8")
    print("    " + f"{'variant':<8}" + " ".join(f"{'m='+str(mm):>7}" for mm in M_GRID))
    e4 = {}
    for v in VARIANTS:
        cells = []
        for mm in M_GRID:
            w, tot, _ = rate(v, curve, mm, 8, SEEDS)
            cells.append(w)
        e4[v] = cells
        print(f"    {v:<8}" + " ".join(f"{str(c)+'/'+str(len(SEEDS)):>7}" for c in cells))

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP E5: is sv/pv a PREDICTOR of the wall in the d-eliminated lattice?")
    print("-" * 78)
    print("With the trivial vector gone (V2), sv is a genuine lambda_1 proxy and pv")
    print("is the planted-vector norm, so BDD theory predicts success iff pv < sv,")
    print("i.e. sv/pv > 1.  Every curve-level separator tried between 2026-06-21 and")
    print("2026-06-29 failed; this one is instance-level, computable before running")
    print("the attack, and has a parameter-free threshold at 1.0.")
    print("GH = Gaussian-heuristic lambda_1 of the (kannan-free) d-eliminated lattice:")
    print("    det = n^(m-1) * S_K1^m * S_K2^m,  dim = 2m,  GH = sqrt(dim/2*pi*e)*det^(1/dim)\n")

    def gh_v3(n, m, k1_bound, k2_bound):
        S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
        dim = 2 * m
        logdet = (m - 1) * math.log(n) + m * math.log(S_K1) + m * math.log(S_K2)
        return math.sqrt(dim / (2 * math.pi * math.e)) * math.exp(logdet / dim)

    e5_rows = []
    for label, curve, _k1, m in curves:
        p, b, n, lam, G = curve
        k2_bound = math.isqrt(n) + 1
        print(f"  {label}  (m={m}, n={n}, K2={k2_bound}, lam*={lam_star(lam, n):.4f})")
        print(f"    {'K1':>4} {'eff':>7} {'V2 win':>8} {'sv/pv':>8} {'GH/pv':>8} {'pred(sv/pv>1)':>14}")
        for k1 in K1_GRID:
            w, tot, ratio = rate('V2', curve, m, k1, SEEDS)
            _, _, r3 = rate('V3', curve, m, k1, SEEDS)
            pvs = []
            for s in SEEDS:
                d_s = random.Random(s + 7777).randint(1, n - 1)
                sg = gen_signatures(G, d_s, m, n, lam, p, k1, k2_bound, s)
                if len(sg) == m:
                    pvs.append(norm(planted_v2(sg, n, k1, k2_bound, kannan=False)))
            pv_mean = sum(pvs) / len(pvs) if pvs else float('nan')
            ghpv = gh_v3(n, m, k1, k2_bound) / pv_mean
            eff = k1 * k2_bound / n
            pred = 'success' if ratio > 1.0 else 'failure'
            actual = 'success' if w == tot else ('partial' if w else 'failure')
            mark = 'ok' if pred == actual else ('--' if actual == 'partial' else 'MISS')
            e5_rows.append((label, k1, eff, w, tot, ratio, ghpv, pred, actual, mark))
            print(f"    {k1:>4} {eff:>7.3f} {str(w)+'/'+str(tot):>8} {ratio:>8.3f} "
                  f"{ghpv:>8.3f} {pred:>14}  [{actual} {mark}]")
        print()

    hits = sum(1 for r in e5_rows if r[9] == 'ok')
    miss = sum(1 for r in e5_rows if r[9] == 'MISS')
    part = sum(1 for r in e5_rows if r[9] == '--')
    print(f"  sv/pv>1 rule:  {hits} correct, {miss} wrong, {part} partial-outcome rows"
          f"  (of {len(e5_rows)})")
    ghhits = sum(1 for r in e5_rows
                 if ('success' if r[6] > 1.0 else 'failure') == r[8])
    print(f"  GH/pv>1 rule:  {ghhits} correct (of {len(e5_rows)}) — closed form, no LLL needed")

    # ---------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    for label, _c, _k1, _m in curves[1:]:
        parts = []
        for v in VARIANTS:
            cells = e3[(label, v)]
            wall = max([K1_GRID[i] for i, c in enumerate(cells) if c == len(SEEDS)], default=0)
            parts.append(f"{v}:K1<={wall}")
        print(f"  {label:<20} wall  " + "  ".join(parts))
    print("  T4b(K1=8) totals over m in " + str(M_GRID) + ":")
    for v in VARIANTS:
        print(f"    {v}: {sum(e4[v])}/{len(M_GRID) * len(SEEDS)}")
    print("=" * 78)



if __name__ == "__main__":
    main()
