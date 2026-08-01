"""
GLV-HNP Phase 2, Thread 23: make the planted vector lambda_1 by removing the
d-column (projection along e_m).

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, result T5):
  In the Phase-2 lattice of glv_hnp_phase2_20bit.py:262 (`build_glv_lattice`)
  the shortest vector is ALWAYS the trivial vector n*S_D*e_m -- 100% of its
  energy in the d column, |sv[m]|/n = 1.0000 exactly.  Reason: d is uniform in
  [0,n), so the planted vector's d-coordinate contributes (n*S_D)^2/3 while the
  lattice contains n*S_D*e_m of norm exactly n*S_D.  No choice of S_D removes
  this: both vectors scale linearly in S_D.  Recovery is therefore a BDD/coset
  condition (shortest among vectors with last coordinate +-S_KANNAN), not SVP.

Thread 23 (proposed 2026-07-29):
  Reformulate so the target IS lambda_1.  d does not need to be a lattice
  coordinate at all: once any single (k1_i, k2_i) pair is known,
      d = (k1_i + lam*k2_i - A_i) * B_i^{-1}  (mod n).
  So delete the d column.  The row `(B_i*S_K1)_i | S_D` becomes
  `(B_i*S_K1)_i` with no distinguished column, which is exactly the orthogonal
  projection of L along e_m, identified with a lattice of rank 2m+1.

  Geometry of the change (K = k1_bound, and eff = K1*K2/n):
    ORIGINAL   dim 2m+2   det = n^(3m+1)/(K1*K2)^m   ||v_planted||^2 has the
                                                     extra (n*S_D)^2/3 term
    PROJECTED  dim 2m+1   det = n^(2m)/eff^m         d-term gone
  the d row lowers the k1-block determinant from n^m to n^(m-1) (the k1
  columns become {(c_i) : c_i = d*B_i mod n}), so one dimension AND one factor
  of n leave together.  Predicted consequence: sv/pv rises above 1.

Falsifier (stated verbatim in the 2026-07-29 log):
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Experiments:
  E1  geometry of both lattices on the 3 historical curves: dim, det^(1/dim),
      Gaussian-heuristic lambda_1, ||v_planted||, observed sv, sv/pv, and
      whether the LLL shortest vector IS +-v_planted.
  E2  the T4 grid re-run on both lattices: curves 12-bit/2557 (lam*=0.340) and
      12-bit/2677 (lam*=0.070), m=12, K1 in {2,3,4,6,8,12,16,24}, 5 seeds.
      This is the falsifier proper -- does the K1 wall move?
  E3  m-sweep at the wall (K1=8, the 2026-07-26 failure point) on both
      lattices, m in {6,8,10,12,16,20,24}: does more data now help?
  E4  17-bit sweep at eff=0.15 (where the original recovered 3/20 on
      2026-07-29) on both lattices: does the operating regime widen?

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(4000):
        x = rng.randrange(p)
        rhs = (x * x % p * x + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0: continue
        P = (x, y)
        if ec_mul(P, n, p) is None and P is not None:
            return P
    return None

# ---------------------------------------------------------------------------
# Eisenstein CM / GLV eigenvalue
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p).  [verbatim]"""
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
    if sq is None: return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0: return None, None
    return min(r1, r2), max(r1, r2)

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures (verbatim construction)
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

# ---------------------------------------------------------------------------
# ORIGINAL lattice (2m+2), verbatim from glv_hnp_phase2_20bit.py:262
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
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
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

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

# ---------------------------------------------------------------------------
# PROJECTED lattice (2m+1): the d column is deleted.
#
#   columns:  0..m-1      k1_i   (scale S_K1)
#             m..2m-1     k2_i   (scale S_K2)
#             2m          Kannan (scale S_KANNAN)
#   rows:     i<m         n*S_K1*e_i                       (modular rows)
#             m           (B_i*S_K1)_i                     (the d direction --
#                                                           now a free shift)
#             m+1+i       -lam*S_K1 at col i, S_K2 at col m+i
#             2m+1        (A_i*S_K1)_i, S_KANNAN at col 2m (Kannan embedding)
#
#   2m+2 generators spanning rank 2m+1: LLL returns one zero row.
# ---------------------------------------------------------------------------

def build_glv_lattice_projected(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(2 * m + 2)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][2 * m] = S_KANNAN
    return M

def planted_vector_projected(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_d_projected(M_reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_i, k2_i) off any row with |last| = S_KANNAN, then solve for d.

    d is NOT a coordinate here; it is reconstructed as
        d = (k1_i + lam*k2_i - A_i) * B_i^{-1} mod n
    and cross-checked on all m signatures (a check that needs no knowledge of
    d_secret -- d_secret is used only to score the experiment).
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        k1s, k2s, bad = [], [], False
        for i in range(m):
            a, b = sign * row[i], sign * row[m + i]
            if a % S_K1 or b % S_K2 or not (0 <= a // S_K1 < k1_bound) \
               or not (0 <= b // S_K2 < k2_bound):
                bad = True; break
            k1s.append(a // S_K1); k2s.append(b // S_K2)
        if bad: continue
        B0 = sigs[0]['B'] % n
        if B0 == 0: continue
        d_cand = (k1s[0] + lam * k2s[0] - sigs[0]['A']) * modinv(B0, n) % n
        # self-consistency across all signatures (no oracle needed)
        if all((sigs[i]['A'] + sigs[i]['B'] * d_cand - k1s[i] - lam * k2s[i]) % n == 0
               for i in range(m)):
            if d_cand == d_secret:
                return d_cand
    return None

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(x * x for x in v))

def lattice_det(rows, dim):
    """|det| of the lattice spanned by `rows` (rank = dim), exactly."""
    from sympy import Matrix
    nz = [r for r in rows if any(r)]
    M = Matrix(nz)
    if M.rows == dim:
        return abs(M.det())
    G = M * M.T
    return sympy.sqrt(abs(G.det()))

def gaussian_heuristic(det, dim):
    return math.sqrt(dim / (2 * math.pi * math.e)) * float(det) ** (1.0 / dim)

def lll_rows(M, dim):
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(A.nrows)]

def bkz_rows(M, dim, beta=20):
    A = IntegerMatrix.from_matrix(M)
    BKZ.reduction(A, BKZ.Param(beta))
    return [[A[i][j] for j in range(dim)] for i in range(A.nrows)]

def same_up_to_sign(u, v):
    return u == v or all(a == -b for a, b in zip(u, v))

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def run_one(curve, m, d_secret, k1_bound, seed, projected, use_bkz=False, beta=20):
    """Returns dict with ok / planted-norm / shortest-norm / is-lambda1."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if projected:
        dim = 2 * m + 1
        M = build_glv_lattice_projected(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector_projected(sigs, n, k1_bound, k2_bound)
        rows = bkz_rows(M, dim, beta) if use_bkz else lll_rows(M, dim)
        ok = recover_d_projected(rows, sigs, n, lam, k1_bound, k2_bound,
                                 d_secret) is not None
    else:
        dim = 2 * m + 2
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
        rows = bkz_rows(M, dim, beta) if use_bkz else lll_rows(M, dim)
        _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
        ok = recover_d(rows, m, n, S_KAN, d_secret) is not None
    nz = [r for r in rows if any(r)]
    sv = min(nz, key=norm)
    return {'ok': ok, 'pn': norm(pv), 'sn': norm(sv),
            'is_l1': same_up_to_sign(sv, pv), 'rows': rows, 'pv': pv,
            'sigs': sigs, 'dim': dim}

def success_rate(curve, m, k1_bound, seeds, projected, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    wins, ratios, l1hits = 0, [], 0
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = run_one(curve, m, d_trial, k1_bound, seed, projected, use_bkz, beta)
        if r is None: continue
        wins += bool(r['ok']); l1hits += bool(r['is_l1'])
        ratios.append(r['sn'] / r['pn'])
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, len(seeds), mr, l1hits

# ---------------------------------------------------------------------------
# Curve construction
# ---------------------------------------------------------------------------

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        G = find_generator(p, b, n)
        if G is not None:
            return b, G
    return None, None

def search_curves(lo, hi, want, seed=7):
    """Prime-order j=0 GLV curves with n in [lo,hi]."""
    out = []
    p = sympy.nextprime(lo)
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_c = p + 1 - t
                    if n_c < 8 or not sympy.isprime(n_c) or n_c % 3 != 1:
                        continue
                    lam, _ = glv_roots(n_c)
                    if lam is None: continue
                    b, G = build_curve(p, n_c, seed)
                    if G is None: continue
                    out.append((f"p={p},n={n_c}", (p, b, n_c, lam, G)))
                    break
        p = sympy.nextprime(p)
    return out


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2, Thread 23 — projected lattice (d column deleted)")
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

    # ===========================================================================
    # E1 — geometry of the two lattices
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E1: lattice geometry, original (2m+2) vs projected (2m+1)")
    print("-" * 78)
    print("det, GH = Gaussian-heuristic lambda_1, pv = ||v_planted||,")
    print("sv = ||shortest LLL row||, l1 = is the LLL shortest vector +-v_planted?\n")
    print(f"{'curve':<18} {'lat':<5} {'dim':>4} {'det^(1/dim)':>12} {'GH':>11} "
          f"{'pv':>11} {'sv':>11} {'sv/pv':>7} {'pv/GH':>7} {'l1':>5} {'rec':>4}")

    for label, curve, k1, m in hist:
        p, b, n, lam, G = curve
        d0 = random.Random(42 + 7777).randint(1, n - 1)
        for projected in (False, True):
            r = run_one(curve, m, d0, k1, 42, projected)
            det = lattice_det(r['rows'], r['dim'])
            gh = gaussian_heuristic(det, r['dim'])
            print(f"{label:<18} {'proj' if projected else 'orig':<5} {r['dim']:>4} "
                  f"{float(det) ** (1.0 / r['dim']):>12.1f} {gh:>11.1f} "
                  f"{r['pn']:>11.1f} {r['sn']:>11.1f} {r['sn']/r['pn']:>7.3f} "
                  f"{r['pn']/gh:>7.3f} {str(r['is_l1']):>5} {str(r['ok']):>4}")

    # ===========================================================================
    # E2 — the falsifier: does the K1 wall move?
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E2: T4 grid re-run (m=12, 5 seeds) — ORIGINAL vs PROJECTED")
    print("-" * 78)
    print("2026-07-29 T4 baseline (original lattice, m=12):")
    print("  12-bit/2557 (lam*=0.340): 5,5,5,5,5,4,1,0 for K1=2,3,4,6,8,12,16,24")
    print("  12-bit/2677 (lam*=0.070): 5,5,5,2,0,0,0,0 for the same K1 grid\n")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
    M_T4 = 12
    T4 = [("12-bit/2557", hist[1][1]), ("12-bit/2677", hist[2][1])]

    for label, curve in T4:
        p, b, n, lam, G = curve
        print(f"{label}  n={n}  lam={lam}  lam*={lam_star(lam,n):.4f}")
        print(f"  {'lattice':<9} " + " ".join(f"K1={k:<5}" for k in K1_GRID))
        for projected in (False, True):
            cells = []
            for k1 in K1_GRID:
                w, t, mr, l1 = success_rate(curve, M_T4, k1, SEEDS, projected)
                cells.append(f"{w}/{t}".ljust(8))
            print(f"  {'proj' if projected else 'orig':<9} " + " ".join(cells))
        # sv/pv and lambda_1 hit-rate at the historical failure point K1=8
        for projected in (False, True):
            w, t, mr, l1 = success_rate(curve, M_T4, 8, SEEDS, projected)
            print(f"    K1=8 {'proj' if projected else 'orig'}: "
                  f"mean sv/pv = {mr:.3f}, planted-is-lambda_1 in {l1}/{t} runs")
        print()

    # ===========================================================================
    # E3 — does more data help now?  (the T4b question, K1=8)
    # ===========================================================================
    print("-" * 78)
    print("E3: m-sweep at the historical wall K1=8, curve 12-bit/2677 (lam*=0.070)")
    print("-" * 78)
    print("2026-07-29 T4b baseline (original): m=8/12/16/24/32 -> 0,0,1,0,1 of 5\n")
    M_GRID = [6, 8, 10, 12, 16, 20, 24]
    curve_2677 = hist[2][1]
    print(f"  {'lattice':<9} " + " ".join(f"m={m:<5}" for m in M_GRID))
    for projected in (False, True):
        cells = []
        for m in M_GRID:
            w, t, mr, l1 = success_rate(curve_2677, m, 8, SEEDS, projected)
            cells.append(f"{w}/{t}".ljust(7))
        print(f"  {'proj' if projected else 'orig':<9} " + " ".join(cells))

    # ===========================================================================
    # E4 — 17-bit sweep at eff = 0.15
    # ===========================================================================
    print("\n" + "-" * 78)
    print("E4: 17-bit curves at eff = K1*K2/n = 0.15, m=12, 5 seeds")
    print("-" * 78)
    print("2026-07-29 baseline (original lattice): 3/20 curves recovered 5/5\n")

    curves17 = search_curves(2 ** 16, 2 ** 17, want=12)
    print(f"found {len(curves17)} 17-bit j=0 GLV curves\n")
    print(f"{'curve':<24} {'lam*':>7} {'K1':>4} {'orig':>7} {'proj':>7} "
          f"{'sv/pv o':>8} {'sv/pv p':>8} {'l1 p':>5}")
    tot_o = tot_p = 0
    for label, curve in curves17:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(0.15 * n / k2))
        wo, t, mro, _ = success_rate(curve, 12, k1, SEEDS, False)
        wp, _, mrp, l1p = success_rate(curve, 12, k1, SEEDS, True)
        tot_o += (wo == t); tot_p += (wp == t)
        print(f"{label:<24} {lam_star(lam,n):>7.4f} {k1:>4} {str(wo)+'/'+str(t):>7} "
              f"{str(wp)+'/'+str(t):>7} {mro:>8.3f} {mrp:>8.3f} {l1p:>5}")
    print(f"\nfull recovery (5/5): original {tot_o}/{len(curves17)}, "
          f"projected {tot_p}/{len(curves17)}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
