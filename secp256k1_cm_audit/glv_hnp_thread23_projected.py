"""
GLV-HNP Phase 2, Thread 23: make the planted vector lambda_1.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5):
  The Phase-2 lattice L (dim 2m+2, built by build_glv_lattice in
  glv_hnp_phase2_20bit.py:262) always has

      lambda_1(L) = n * S_D * e_m        ("the trivial vector")

  because n*row_m - sum_i B_i*row_i = n*S_D*e_m.  Its norm is n*S_D while
  ||v_planted|| ~ n*S_D*sqrt(2m/3 + 4/3), so the planted vector is NEVER the
  shortest vector, for any m >= 1 and any choice of S_D (both scale linearly
  in S_D).  Recovery is therefore a BDD/coset condition, not SVP.

  The trivial vector carries no information: d is only defined mod n, and
  n*S_D*e_m is exactly "add n to d".

This script tests the fix proposed on 2026-07-29:

  P1  PROJECTION.  Quotient L by the trivial direction, i.e. drop column m
      entirely.  pi(L) = L' has rank 2m+1 and det L' = det L / (n*S_D)
      = n^m * S_K1^m * S_K2^m.  d is recovered *algebraically* afterwards
      from any single signature:  d = (k1_i + lam*k2_i - A_i) * B_i^-1 mod n.
      The generating set (2m+2 rows, rank 2m+1) is fed to fplll directly;
      fplll handles rank deficiency and emits one zero row.

  P2  CENTERING.  k1 ~ U[0,K1) and k2 ~ U[0,K2) are not centred, so each
      contributes X^2/3 instead of X^2/12 to ||v_planted||^2.  Replacing
      A_i by A_i - K1/2 - lam*(K2/2) mod n re-centres both nonce halves at
      zero cost.  This is a pure change of target, not of lattice.

  P3  = P1 + P2.

Falsifier stated on 2026-07-29:
  "if sv/pv rises above 1 after the reformulation and the K1 wall on the
   lam*=0.07 curve (currently K1 ~ 4-6) moves outward, the reformulation is
   a real improvement; if the wall stays at K1 ~ 4-6, the wall is
   information-theoretic and Phase 2 is at its ceiling."

Gaussian-heuristic prediction computed before running (m=10, n=2647, K2=52):
  GH(L')/n = sqrt(D/2*pi*e) * det^(1/D) / n  with D = 2m+1
  ||pv'||/n = sqrt(2m/3 + 1) = 2.77   (P1, uncentred)
  ||pv'||/n = sqrt(m/6   + 1) = 1.63   (P3, centred)
  => predicted wall  K1 ~ 3  for P1,  K1 ~ 12 for P3
     (versus the measured K1 ~ 4-6 wall for the original formulation)

Run: python3 glv_hnp_thread23_projected.py            (all sections)
     python3 glv_hnp_thread23_projected.py P3         (one section)
"""

import math
import random
import sys
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

ALL_SECTIONS = ('P0', 'P1', 'P2', 'P3', 'P4')
RUN = tuple(a for a in sys.argv[1:] if a in ALL_SECTIONS) or ALL_SECTIONS
# P4 reports the measured wall from P2, so it needs P2's numbers.
if 'P4' in RUN and 'P2' not in RUN:
    RUN = RUN + ('P2',)

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

def find_generator(p, b, n):
    rng = random.Random(12345)
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
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_eigenvalue(n):
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
    return min(lam, n - lam) / n

# ---------------------------------------------------------------------------
# Signatures and scales -- verbatim from glv_hnp_phase2_lambda_threshold.py
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN)"""
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
# The four lattice formulations
# ---------------------------------------------------------------------------
#
#   orig      (2m+2)-dim, uncentred  -- the 2026-07-29 baseline
#   cent      (2m+2)-dim, centred
#   proj      (2m+1)-dim, uncentred  -- d column dropped
#   projcent  (2m+1)-dim, centred
#
# "centred" replaces A_i by A_i - K1//2 - lam*(K2//2) mod n, which shifts the
# unknowns to c1 = k1 - K1//2 in [-K1/2, K1/2) and c2 = k2 - K2//2 likewise.

def shifted_A(sig, lam, n, k1_bound, k2_bound, centred):
    if not centred:
        return sig['A'] % n
    return (sig['A'] - k1_bound // 2 - lam * (k2_bound // 2)) % n

def build_lattice(sigs, n, lam, k1_bound, k2_bound, projected, centred):
    """Return (M, dim).  Rows are lattice generators.

    Un-projected layout (dim 2m+2):
        col i        : k1_i * S_K1        i = 0..m-1
        col m        : d * S_D
        col m+1+i    : k2_i * S_K2
        col 2m+1     : S_KANNAN
    Projected layout (dim 2m+1): column m removed, everything after shifts
    down by one.  The generating set then has 2m+2 rows of rank 2m+1;
    fplll reduces it and emits one zero row.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if projected:
        dim = 2 * m + 1
        c_d, c_k2, c_kan = None, m, 2 * m
    else:
        dim = 2 * m + 2
        c_d, c_k2, c_kan = m, m + 1, 2 * m + 1

    rows = []
    for i in range(m):                                   # n*S_K1 * e_i
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim                                        # the d row
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    if c_d is not None:
        r[c_d] = S_D
    rows.append(r)
    for i in range(m):                                   # the k2 rows
        r = [0] * dim
        r[i] = -lam * S_K1
        r[c_k2 + i] = S_K2
        rows.append(r)
    r = [0] * dim                                        # the target row
    for i in range(m):
        r[i] = shifted_A(sigs[i], lam, n, k1_bound, k2_bound, centred) * S_K1
    r[c_kan] = S_KAN
    rows.append(r)
    return rows, dim

def planted_vector(sigs, d_secret, n, lam, k1_bound, k2_bound, projected, centred):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if projected:
        dim, c_d, c_k2, c_kan = 2 * m + 1, None, m, 2 * m
    else:
        dim, c_d, c_k2, c_kan = 2 * m + 2, m, m + 1, 2 * m + 1
    off1 = k1_bound // 2 if centred else 0
    off2 = k2_bound // 2 if centred else 0
    v = [0] * dim
    for i in range(m):
        v[i] = (sigs[i]['k1'] - off1) * S_K1
        v[c_k2 + i] = (sigs[i]['k2'] - off2) * S_K2
    if c_d is not None:
        v[c_d] = d_secret * S_D
    v[c_kan] = S_KAN
    return v

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_unprojected(reduced, sigs, n, lam, k1_bound, k2_bound, centred, d_secret):
    """2026-07-29 recovery, verbatim in spirit: read d off the d column."""
    m = len(sigs)
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 2
    for row in reduced:
        if abs(row[dim - 1]) != S_KAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

def recover_projected(reduced, sigs, n, lam, k1_bound, k2_bound, centred, d_secret):
    """d is not in the lattice any more: read (k1_i, k2_i) off the vector and
    solve d = (k1_i + lam*k2_i - A_i) * B_i^-1 mod n from signature 0."""
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + 1
    off1 = k1_bound // 2 if centred else 0
    off2 = k2_bound // 2 if centred else 0
    for row in reduced:
        if abs(row[dim - 1]) != S_KAN: continue
        sign = 1 if row[dim - 1] > 0 else -1
        if sign * row[0] % S_K1 or sign * row[m] % S_K2:
            continue
        k1_0 = sign * row[0] // S_K1 + off1
        k2_0 = sign * row[m] // S_K2 + off2
        s0 = sigs[0]
        d_cand = (k1_0 + lam * k2_0 - s0['A']) * modinv(s0['B'], n) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

def run_one(curve, m, k1_bound, seed, projected, centred, use_bkz=False, beta=20,
            want_geometry=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d_secret = random.Random(seed * 7919 + 13).randint(1, n - 1)
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    rows, dim = build_lattice(sigs, n, lam, k1_bound, k2_bound, projected, centred)
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    rec = (recover_projected if projected else recover_unprojected)(
        reduced, sigs, n, lam, k1_bound, k2_bound, centred, d_secret)
    out = {'ok': rec is not None}
    if want_geometry:
        pv = planted_vector(sigs, d_secret, n, lam, k1_bound, k2_bound,
                            projected, centred)
        nz = [r for r in reduced if any(r)]
        sv = min(nz, key=norm)
        out.update({'sv': norm(sv), 'pv': norm(pv), 'ratio': norm(sv) / norm(pv),
                    'rank': len(nz), 'dim': dim, 'n': n})
    return out

def success(curve, m, k1_bound, seeds, projected, centred, use_bkz=False, beta=20):
    w = 0
    for s in seeds:
        r = run_one(curve, m, k1_bound, s, projected, centred, use_bkz, beta)
        if r and r['ok']:
            w += 1
    return w

def gh_over_n(n, m, k1_bound, k2_bound, projected):
    """Gaussian heuristic lambda_1 estimate, normalised by n."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if projected:
        dim = 2 * m + 1
        logdet = m * math.log(n * S_K1) + m * math.log(S_K2)
    else:
        dim = 2 * m + 2
        logdet = (m * math.log(n * S_K1) + math.log(S_D)
                  + m * math.log(S_K2) + math.log(S_KAN))
    log_gh = 0.5 * math.log(dim / (2 * math.pi * math.e)) + logdet / dim
    return math.exp(log_gh) / n

def pv_over_n_expected(m, centred, projected):
    """E[||v_planted||]/n under k1*S_K1, k2*S_K2 ~ U[0,n) (or centred)."""
    var = 1.0 / 12 if centred else 1.0 / 3
    t = 2 * m * var + 1.0                     # k1 block + k2 block + Kannan
    if not projected:
        t += 1.0 / 3                          # the d column, d ~ U[0,n)
    return math.sqrt(t)

# ===========================================================================

print("=" * 78)
print("Thread 23 — reformulate the Phase-2 GLV lattice so the planted vector")
print("            can be lambda_1 (projection + centring)")
print(f"sections: {' '.join(RUN)}")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]
CONFIGS = [('orig', False, False), ('cent', False, True),
           ('proj', True, False), ('projcent', True, True)]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 2026-07-26)
HIST = [
    # label,          p,    b, n,    lam,  m
    ("12-bit/2557",   2557, 2, 2659, 1755, 8),
    ("12-bit/2677",   2677, 2, 2647, 185,  10),
]
hist = []
for label, p, b, n, lam, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), m))

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32]
wall = {}

# ---------------------------------------------------------------------------
def exp_P0():
    print("\n" + "-" * 78)
    print("EXP P0: lattice sanity — the projection drops exactly one dimension")
    print("-" * 78)
    print("Check that pi(L) really is L / <n*S_D*e_m>: LLL returns rank 2m+1 from")
    print("the 2m+2 generators (one zero row), and log2 det drops by log2(n*S_D).\n")
    print(f"{'curve':<14} {'cfg':<6} {'dim':>4} {'rank':>5} {'log2 det':>10} "
          f"{'drop':>8} {'log2(n*S_D)':>12}")
    for label, curve, m in hist:
        p, b, n, lam, G = curve
        k2b, k1b = math.isqrt(n) + 1, 8
        S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
        base = (m * math.log2(n * S_K1) + m * math.log2(S_K2))
        dets = {}
        for cfg, proj in (('orig', False), ('proj', True)):
            r = run_one(curve, m, k1b, 42, proj, False, want_geometry=True)
            ld = base + (0.0 if proj else math.log2(S_D) + math.log2(S_KAN))
            dets[cfg] = ld
            drop = '—' if proj is False else f"{dets['orig'] - ld:.2f}"
            print(f"{label:<14} {cfg:<6} {r['dim']:>4} {r['rank']:>5} "
                  f"{ld:>10.2f} {drop:>8} {math.log2(n * S_D):>12.2f}")

# ---------------------------------------------------------------------------
def exp_P1():
    print("\n" + "-" * 78)
    print("EXP P1: does the planted vector become lambda_1?   (sv/pv ratio)")
    print("-" * 78)
    print("2026-07-29 measured sv/pv in [0.34, 0.61] for 'orig' — the planted")
    print("vector is 2-3x LONGER than lambda_1.  A ratio of 1.000 means the")
    print("shortest vector found IS the planted vector.\n")
    print(f"{'curve':<14} {'K1':>3} {'cfg':<9} {'sv/pv':>7} {'|pv|/n':>8} "
          f"{'|sv|/n':>8} {'GH/n':>7} {'E|pv|/n':>8} {'rec':>4}")
    for label, curve, m in hist:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for k1b in (2, 8):
            for cfg, proj, cen in CONFIGS:
                r = run_one(curve, m, k1b, 42, proj, cen, want_geometry=True)
                print(f"{label:<14} {k1b:>3} {cfg:<9} {r['ratio']:>7.3f} "
                      f"{r['pv']/n:>8.3f} {r['sv']/n:>8.3f} "
                      f"{gh_over_n(n, m, k1b, k2b, proj):>7.3f} "
                      f"{pv_over_n_expected(m, cen, proj):>8.3f} "
                      f"{('Y' if r['ok'] else 'n'):>4}")
            print()

# ---------------------------------------------------------------------------
def exp_P2():
    print("-" * 78)
    print("EXP P2: the K1 wall — does the reformulation move it outward?")
    print("-" * 78)
    print("2026-07-29 T4 grid, replayed for all four formulations, 5 seeds each.")
    print("Baseline ('orig'): 2557 walls at K1 ~ 12-16, 2677 walls at K1 ~ 4-6.\n")
    for label, curve, m in hist:
        p, b, n, lam, G = curve
        print(f"{label}  (n={n}, lam*={lam_star(lam, n):.4f}, m={m})")
        print(f"  {'cfg':<9}" + "".join(f"{k:>5}" for k in K1_GRID))
        for cfg, proj, cen in CONFIGS:
            cells, last_ok = [], None
            for k1b in K1_GRID:
                w = success(curve, m, k1b, SEEDS, proj, cen)
                cells.append(f"{w:>5}")
                if w >= 4:
                    last_ok = k1b
            wall[(label, cfg)] = last_ok
            print(f"  {cfg:<9}{''.join(cells)}")
        print()
    print("Wall = largest K1 with >=4/5 recoveries:")
    print(f"{'curve':<14}" + "".join(f"{c:>10}" for c, _, _ in CONFIGS))
    for label, curve, m in hist:
        print(f"{label:<14}" + "".join(f"{str(wall[(label, c)]):>10}"
                                       for c, _, _ in CONFIGS))

# ---------------------------------------------------------------------------
def find_j0_curves(lo, hi, want):
    """Fresh prime-order j=0 GLV curves with p in [lo, hi)."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 4 or not sympy.isprime(nc) or nc % 3 != 1:
                        continue
                    lam = glv_eigenvalue(nc)
                    if lam is None:
                        continue
                    for b_try in range(1, 60):
                        G0 = find_generator(p, b_try, nc)
                        if G0 is not None:
                            out.append((f"p={p}/n={nc}", (p, b_try, nc, lam, G0)))
                            break
                    break
        p = int(sympy.nextprime(p))
    return out

def exp_P3():
    print("\n" + "-" * 78)
    print("EXP P3: independent replication on fresh 17-bit curves")
    print("-" * 78)
    print("The two historical curves are 2 data points.  Sweep 8 fresh j=0 GLV")
    print("curves across the K1 range that brackets the 'orig' wall at 17 bits.")
    print("Wall position scales with eff = K1*K2/n, and K2 = isqrt(n)+1 ~ 256")
    print("here (vs 52 at 12 bits), so the discriminating K1 range is ~5x larger.\n")
    M17 = 12
    K1_FRESH = [16, 24, 32, 48, 64]
    fresh = find_j0_curves(2 ** 16, 2 ** 17, 8)
    ns = len(SEEDS)
    print(f"{len(fresh)} fresh 17-bit j=0 curves, m={M17}, {ns} seeds, "
          f"K1 in {K1_FRESH}\n")
    tot = {(c, k): 0 for c, _, _ in CONFIGS for k in K1_FRESH}
    print(f"{'curve':<20} {'lam*':>7} {'cfg':<9}"
          + "".join(f"{'K1=' + str(k):>7}" for k in K1_FRESH))
    for label, curve in fresh:
        p, b, n, lam, G = curve
        for cfg, proj, cen in CONFIGS:
            cells = []
            for k1b in K1_FRESH:
                w = success(curve, M17, k1b, SEEDS, proj, cen)
                tot[(cfg, k1b)] += w
                cells.append(f"{w:>7}")
            print(f"{label:<20} {lam_star(lam, n):>7.4f} {cfg:<9}"
                  + "".join(cells))
        print()
    denom = len(fresh) * ns
    print(f"{'TOTAL (/' + str(denom) + ')':<20} {'':>7} "
          + " " * 2 + "".join(f"{'K1=' + str(k):>7}" for k in K1_FRESH))
    for cfg, _, _ in CONFIGS:
        print(f"{'':<20} {'':>7} {cfg:<9}"
              + "".join(f"{tot[(cfg, k)]:>7}" for k in K1_FRESH))
    print("\nWall = largest K1 with >=80% of trials recovering:")
    for cfg, _, _ in CONFIGS:
        w = [k for k in K1_FRESH if tot[(cfg, k)] >= 0.8 * denom]
        print(f"  {cfg:<9} {max(w) if w else '<' + str(K1_FRESH[0])}")

# ---------------------------------------------------------------------------
def exp_P4():
    print("\n" + "-" * 78)
    print("EXP P4: GH crossing vs measured wall")
    print("-" * 78)
    print("Model: the planted vector is lambda_1 iff E||pv|| < GH(L).  Report the")
    print("largest K1 satisfying that, next to the measured wall from P2.\n")
    print(f"{'curve':<14} {'cfg':<9} {'GH-crossing K1':>15} {'measured wall':>14}")
    for label, curve, m in hist:
        p, b, n, lam, G = curve
        k2b = math.isqrt(n) + 1
        for cfg, proj, cen in CONFIGS:
            pv = pv_over_n_expected(m, cen, proj)
            cross = None
            for k1b in K1_GRID:
                if pv < gh_over_n(n, m, k1b, k2b, proj):
                    cross = k1b
            print(f"{label:<14} {cfg:<9} {str(cross):>15} "
                  f"{str(wall.get((label, cfg))):>14}")

# ---------------------------------------------------------------------------
for _name in ALL_SECTIONS:
    if _name in RUN:
        globals()['exp_' + _name]()

print("\n" + "=" * 78)
print("done")
print("=" * 78)
