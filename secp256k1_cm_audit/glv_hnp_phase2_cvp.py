"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the target IS the
shortest / closest vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, "T5"):
  The Phase-2 Kannan lattice (dim 2m+2, `glv_hnp_phase2_20bit.py:263`) always
  has  n*S_D*e_m  as its shortest vector: 100% of the reduced-basis shortest
  vector's energy sits in the d-column, |sv[m]|/n = 1.0000 on every curve
  tested, and sv/pv in [0.337, 0.368] for successes and failures alike.
  Recovery there is a BDD/coset condition, not an SVP condition.  That run
  proposed Thread 23: drop the d-column and the Kannan column, and solve the
  underlying CVP directly.

This script implements the reformulation and measures whether it moves the
empirical K1 wall, and separates TWO independent changes that the Kannan form
conflates:

  (R) REFORMULATION.  Replace the (2m+2)-dim Kannan embedding by the
      2m-dim lattice
          L = < n*S_K1*e_i ;  (B_i*S_K1)_i ;  (-lam*S_K1*e_i , S_K2*e_i) >
      and solve CVP against the target  t_i = -A_i*S_K1  (block 1), 0 (block 2).
      The unknown d is the *coefficient* of the (B_i*S_K1) generator, so it is
      read off the CVP coefficient vector rather than off a lattice coordinate.
      This removes both the trivial n*e_m vector and the d^2 ~ n^2/3 term from
      the error norm.

  (C) CENTERING.  The Kannan planted vector has coordinates k1_i*S_K1 and
      k2_i*S_K2 with k1_i ~ U[0,K1), k2_i ~ U[0,K2) — i.e. errors centered at
      n/2, not at 0, so E[x^2] = n^2/3 instead of n^2/12.  Offsetting the
      target by (K1/2)*S_K1 and (K2/2)*S_K2 halves ||error||.  Centering is
      independent of (R): it applies to the Kannan form too, by offsetting the
      A_i row.  Prior Phase-2 code never centered.

Four arms are measured so the two effects can be attributed separately:
      A0 = Kannan, uncentered  (the historical baseline, verbatim construction)
      A1 = Kannan, centered
      B0 = CVP,    uncentered
      B1 = CVP,    centered

Experiments:
  E1  correctness of the CVP formulation (recovers d where A0 does)
  E2  the 2026-07-29 falsifier: does the target become lambda_1 / does the
      closest-vector gap ||e||/GH drop below 1?
  E3  the K1 wall grid of 2026-07-29 T4, re-run on all four arms
  E4  Gaussian-heuristic vs information-theoretic wall, per arm

Run: python3 glv_hnp_phase2_cvp.py
Deps: fpylll, cysignals, sympy   (pip install cysignals fpylll sympy)
"""

import math
import random
import time

from fpylll import IntegerMatrix, LLL, BKZ, GSO, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) — verbatim from
# glv_hnp_phase2_lambda_threshold.py so results are directly comparable.
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

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures (verbatim from glv_hnp_phase2_20bit.py:235)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) — verbatim from the 2026-06-15 construction."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# ARM A — Kannan embedding, dim 2m+2 (A0 uncentered = historical baseline)
# ---------------------------------------------------------------------------

def build_kannan(sigs, n, lam, k1_bound, k2_bound, center=False):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1 = (k1_bound // 2) if center else 0
    c2 = (k2_bound // 2) if center else 0
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
        # offsetting the Kannan row by -c1*S_K1 centers the k1-errors, and
        # offsetting the k2-block coordinate by -c2*S_K2 centers the k2-errors.
        M[2 * m + 1][i] = (sigs[i]['A'] - c1) * S_K1
        M[2 * m + 1][m + 1 + i] = -c2 * S_K2
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def kannan_planted(sigs, d_secret, n, k1_bound, k2_bound, center=False):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    c1 = (k1_bound // 2) if center else 0
    c2 = (k2_bound // 2) if center else 0
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - c2) * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def kannan_recover(reduced, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# ARM B — CVP reformulation, dim 2m (no d column, no Kannan column)
# ---------------------------------------------------------------------------

def build_cvp(sigs, n, lam, k1_bound, k2_bound, center=False):
    """Returns (generators, target, S_K1, S_K2, c1, c2).

    Coordinates: col i (i<m)      = congruence block, scale S_K1
                 col m+i (i<m)    = k2 block,         scale S_K2
    Generators (2m+1 of them, rank 2m — n*(B row) is in the span of the
    n*e_i rows, which is exactly the statement "d is only defined mod n"):
        g_i     = n*S_K1 * e_i                       i = 0..m-1
        g_d     = (B_i*S_K1)_i                       (coefficient = d)
        g_{k2,i}= -lam*S_K1*e_i + S_K2*e_{m+i}       (coefficient = k2_i)
    Target:
        t_i     = (-A_i + c1) * S_K1,   t_{m+i} = c2 * S_K2
    so that for the true (d, k2_i) the lattice vector v satisfies
        v - t = ((k1_i - c1)*S_K1 ,  (k2_i - c2)*S_K2).
    """
    m = len(sigs)
    dim = 2 * m
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    c1 = (k1_bound // 2) if center else 0
    c2 = (k2_bound // 2) if center else 0

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
        target[i] = (-sigs[i]['A'] + c1) * S_K1
        target[m + i] = c2 * S_K2
    return gens, target, S_K1, S_K2, c1, c2

def cvp_error(sigs, n, k1_bound, k2_bound, center=False):
    """The error vector v_true - t that CVP must find."""
    m = len(sigs)
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    c1 = (k1_bound // 2) if center else 0
    c2 = (k2_bound // 2) if center else 0
    e = [0] * (2 * m)
    for i in range(m):
        e[i] = (sigs[i]['k1'] - c1) * S_K1
        e[m + i] = (sigs[i]['k2'] - c2) * S_K2
    return e

def lll_basis(gens, dim):
    """LLL-reduce a rank-deficient generator set; return a basis (zero rows
    dropped).  fpylll pushes the dependent rows to the top as zeros."""
    A = IntegerMatrix.from_matrix(gens)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    rows = [r for r in rows if any(r)]
    return rows

def recover_d_from_error(err, sigs, n, lam, S_K1, S_K2, c1, c2):
    """d from the CVP error vector: A_i + d*B_i = k1_i + lam*k2_i (mod n)."""
    m = len(sigs)
    if any(err[i] % S_K1 for i in range(m)):
        return None
    if any(err[m + i] % S_K2 for i in range(m)):
        return None
    cands = {}
    for i in range(m):
        k1 = err[i] // S_K1 + c1
        k2 = err[m + i] // S_K2 + c2
        B = sigs[i]['B'] % n
        if math.gcd(B, n) != 1:
            continue
        d = (k1 + lam * k2 - sigs[i]['A']) * modinv(B, n) % n
        cands[d] = cands.get(d, 0) + 1
    if not cands:
        return None
    # self-consistency: the true d is agreed on by every equation
    d_best, votes = max(cands.items(), key=lambda kv: kv[1])
    return d_best if votes == m else None

def run_cvp(curve, m, d_secret, k1_bound, seed=42, center=True,
            exact=False, use_bkz=False, bkz_beta=20, want_stats=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return (False, None) if want_stats else False
    gens, target, S_K1, S_K2, c1, c2 = build_cvp(
        sigs, n, lam, k1_bound, k2_bound, center)
    dim = 2 * m
    basis = lll_basis(gens, dim)
    A = IntegerMatrix.from_matrix(basis)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    if exact:
        v = list(CVP.closest_vector(A, tuple(target)))
    else:
        Gso = GSO.Mat(A)
        Gso.update_gso()
        coeffs = Gso.babai(list(target))
        v = [sum(coeffs[r] * A[r][j] for r in range(A.nrows)) for j in range(dim)]
    err = [v[j] - target[j] for j in range(dim)]
    d_cand = recover_d_from_error(err, sigs, n, lam, S_K1, S_K2, c1, c2)
    ok = (d_cand == d_secret)
    if not want_stats:
        return ok
    e_true = cvp_error(sigs, n, k1_bound, k2_bound, center)
    lam1 = min((norm([A[i][j] for j in range(dim)])
                for i in range(A.nrows)), default=float('inf'))
    stats = {
        'err_norm': norm(e_true),
        'found_norm': norm(err),
        'lam1': lam1,
        'gh': gaussian_heuristic(basis, dim),
        'dim': dim,
    }
    return ok, stats

def run_kannan(curve, m, d_secret, k1_bound, seed=42, center=False,
               use_bkz=False, bkz_beta=20, want_stats=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return (False, None) if want_stats else False
    M = build_kannan(sigs, n, lam, k1_bound, k2_bound, center)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = kannan_recover(reduced, m, n, S_KAN, d_secret) is not None
    if not want_stats:
        return ok
    pv = kannan_planted(sigs, d_secret, n, k1_bound, k2_bound, center)
    sv = min((r for r in reduced if any(r)), key=norm)
    stats = {
        'pv_norm': norm(pv),
        'sv_norm': norm(sv),
        'svpv': norm(sv) / norm(pv),
        'sv_d_frac': abs(sv[m]) / max(1.0, norm(sv)),
        'gh': gaussian_heuristic(M, dim),
        'dim': dim,
    }
    return ok, stats

# ---------------------------------------------------------------------------
# Gaussian heuristic
# ---------------------------------------------------------------------------

def log_det(rows, dim):
    """log(det) of a full-rank integer basis, via GSO norms (float-safe)."""
    A = IntegerMatrix.from_matrix(rows)
    M = GSO.Mat(A)
    M.update_gso()
    return 0.5 * sum(math.log(M.get_r(i, i)) for i in range(A.nrows))

def gaussian_heuristic(rows, dim):
    rows = [r for r in rows if any(r)]
    N = len(rows)
    ld = log_det(rows, dim)
    return math.exp(0.5 * math.log(N / (2 * math.pi * math.e)) + ld / N)

# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

ARMS = [
    ('A0 kannan/uncentered', 'kannan', False),
    ('A1 kannan/centered',   'kannan', True),
    ('B0 cvp/uncentered',    'cvp',    False),
    ('B1 cvp/centered',      'cvp',    True),
]

def success_rate(curve, m, k1_bound, seeds, kind, center, exact=False,
                 use_bkz=False, bkz_beta=20):
    wins = 0
    n = curve[2]
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        if kind == 'kannan':
            ok = run_kannan(curve, m, d, k1_bound, s, center, use_bkz, bkz_beta)
        else:
            ok = run_cvp(curve, m, d, k1_bound, s, center, exact,
                         use_bkz, bkz_beta)
        wins += 1 if ok else 0
    return wins, len(seeds)


# ===========================================================================

print("=" * 78)
print("Thread 23 — CVP reformulation of the GLV-HNP Phase 2 lattice")
print("=" * 78)

HIST = [
    # label,              p,    b, n,    lam,  K1, m
    ("8-bit/199",         211,  2, 199,  106,  2,  6),
    ("12-bit/2557",       2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL",  2677, 2, 2647, 185,  8,  10),
]

curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E1: correctness — does the CVP formulation recover d at all?")
print("-" * 78)
print("Historical K1/m per curve; 5 seeds; Babai nearest-plane on the")
print("LLL-reduced basis, plus exact CVP (fplll enumeration).\n")

print(f"{'curve':<18} {'K1':>3} {'m':>3} "
      f"{'A0':>6} {'A1':>6} {'B0':>6} {'B1':>6} {'B1exact':>8}")
for label, curve, k1, m in curves:
    r = {}
    for name, kind, center in ARMS:
        w, t = success_rate(curve, m, k1, SEEDS, kind, center)
        r[name[:2]] = f"{w}/{t}"
    w, t = success_rate(curve, m, k1, SEEDS, 'cvp', True, exact=True)
    print(f"{label:<18} {k1:>3} {m:>3} "
          f"{r['A0']:>6} {r['A1']:>6} {r['B0']:>6} {r['B1']:>6} "
          f"{str(w)+'/'+str(t):>8}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("E2: the 2026-07-29 falsifier — is the target now the extremal vector?")
print("-" * 78)
print("Kannan arms report sv/pv (2026-07-29 T5 reported 0.603/0.517/0.422,")
print("with 100% of sv energy in the d-column).  CVP arms report the BDD")
print("gap ||e||/lambda_1 and ||e||/GH: BDD is uniquely decodable when")
print("||e|| < lambda_1/2, and heuristically solvable while ||e|| < GH.\n")

print(f"{'curve':<18} {'arm':<22} {'dim':>4} {'sv/pv':>7} {'d-frac':>7} "
      f"{'|e|/L1':>8} {'|e|/GH':>7}")
for label, curve, k1, m in curves:
    n = curve[2]
    d = random.Random(42 + 7777).randint(1, n - 1)
    for name, kind, center in ARMS:
        if kind == 'kannan':
            ok, st = run_kannan(curve, m, d, k1, 42, center, want_stats=True)
            print(f"{label:<18} {name:<22} {st['dim']:>4} {st['svpv']:>7.3f} "
                  f"{st['sv_d_frac']:>7.3f} {'-':>8} {'-':>7}")
        else:
            ok, st = run_cvp(curve, m, d, k1, 42, center, want_stats=True)
            print(f"{label:<18} {name:<22} {st['dim']:>4} {'-':>7} {'-':>7} "
                  f"{st['err_norm']/st['lam1']:>8.3f} "
                  f"{st['err_norm']/st['gh']:>7.3f}")
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("E3: the K1 wall — 2026-07-29 T4 grid re-run on all four arms")
print("-" * 78)
print("2026-07-29 T4 (= arm A0) found the wall at K1~12-16 for lam*=0.340 and")
print("K1~4-6 for lam*=0.070, at m=12, and concluded 'it is a K1 wall'.")
print("Falsifier for Thread 23: the reformulation is a real improvement iff")
print("the wall moves outward.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
M_SIGS = 12
WALL = [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]
wall_curves = []
for label, p, b, n, lam in WALL:
    G = find_generator(p, b, n)
    wall_curves.append((label, (p, b, n, lam, G)))

for label, curve in wall_curves:
    n = curve[2]
    lam = curve[3]
    k2 = math.isqrt(n) + 1
    print(f"{label}  n={n}  lam*={lam_star(lam, n):.3f}  m={M_SIGS}  "
          f"K2={k2}")
    hdr = f"{'arm':<22}" + "".join(f"{k:>6}" for k in K1_GRID)
    print(hdr)
    print(f"{'eff = K1*K2/n':<22}" +
          "".join(f"{k * k2 / n:>6.2f}" for k in K1_GRID))
    for name, kind, center in ARMS:
        cells = []
        for k1 in K1_GRID:
            w, t = success_rate(curve, M_SIGS, k1, SEEDS, kind, center)
            cells.append(f"{w}/{t}")
        print(f"{name:<22}" + "".join(f"{c:>6}" for c in cells))
    # exact CVP on the best arm
    cells = []
    for k1 in K1_GRID:
        w, t = success_rate(curve, M_SIGS, k1, SEEDS, 'cvp', True, exact=True)
        cells.append(f"{w}/{t}")
    print(f"{'B1 + exact CVP':<22}" + "".join(f"{c:>6}" for c in cells))
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("E4: predicted walls — Gaussian heuristic vs information theory")
print("-" * 78)
print("GH wall:   solve while ||e|| < GH(L) = sqrt(N/2*pi*e) * det^(1/N).")
print("Info wall: m*log2(n) bits of equations must cover log2(n) + ")
print("           m*log2(K1*K2) bits of unknowns, i.e. eff < n^(-1/m).\n")

print(f"{'curve':<14} {'arm':<22} {'K1':>4} {'eff':>6} {'|e|':>10} "
      f"{'GH':>10} {'|e|/GH':>7}")
for label, curve in wall_curves:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    d = random.Random(42 + 7777).randint(1, n - 1)
    for k1 in (4, 8, 16, 32):
        for name, kind, center in ARMS:
            if kind == 'kannan':
                ok, st = run_kannan(curve, M_SIGS, d, k1, 42, center,
                                    want_stats=True)
                e = st['pv_norm']
            else:
                ok, st = run_cvp(curve, M_SIGS, d, k1, 42, center,
                                 want_stats=True)
                e = st['err_norm']
            print(f"{label:<14} {name:<22} {k1:>4} {k1 * k2 / n:>6.2f} "
                  f"{e:>10.3g} {st['gh']:>10.3g} {e / st['gh']:>7.3f}")
        print()
    print(f"  info-theoretic wall for m={M_SIGS}, n={n}: "
          f"eff < n^(-1/m) = {n ** (-1.0 / M_SIGS):.3f}  "
          f"(K1 < {n ** (-1.0 / M_SIGS) * n / k2:.1f})")
    print()

# ---------------------------------------------------------------------------
# E5 — does centering also erase the curve-level (lam*, nu_hat) effect?
# ---------------------------------------------------------------------------
# 2026-07-29 reported nu_hat = lambda_1(L2)/sqrt(det L2), with
# L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >, as a separator with AUC 0.935 on
# the June C1/C2 classes, and T4 reported a factor-~3 lam*-dependent shift of
# the K1 wall.  Both were measured on arm A0.  If those effects are artifacts
# of the uncentered embedding, they must shrink or vanish on arm B1.

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

def gauss_reduce_2d(u, v):
    """Exact shortest vector of a 2D integer lattice (Lagrange/Gauss)."""
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

def nu_hat(n, lam, S_K1, S_K2):
    """lambda_1(L2)/sqrt(det L2), the 2026-07-29 separator."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    l1 = math.sqrt(w[0] * w[0] + w[1] * w[1])
    return l1 / math.sqrt(n * S_K1 * S_K2)

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

def search_curves(lo, hi, per_bin=2, nbins=8):
    """j=0 GLV curves with p in [lo,hi), n prime, bucketed by lam*."""
    import sympy
    bins = {i: [] for i in range(nbins)}
    p = int(sympy.nextprime(lo))
    while p < hi:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_c = p + 1 - t
                    if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                        continue
                    roots = glv_roots(n_c)
                    if roots is None:
                        continue
                    lam = roots[0]
                    idx = min(nbins - 1, int(lam_star(lam, n_c) / (0.5 / nbins)))
                    if len(bins[idx]) >= per_bin:
                        continue
                    cur = build_curve(p, n_c)
                    if cur is None:
                        continue
                    bins[idx].append((p, cur[1], n_c, lam, cur[4]))
        if all(len(v) >= per_bin for v in bins.values()):
            break
        p = int(sympy.nextprime(p))
    out = []
    for i in range(nbins):
        out.extend(bins[i])
    return out

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    nn = len(xs)
    mx, my = sum(rx) / nn, sum(ry) / nn
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(nn))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(nn)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(nn)))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0

print("-" * 78)
print("E5: does centering erase the curve-level (lam*, nu_hat) effect?")
print("-" * 78)

LO, HI = 1 << 13, 1 << 14
print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
e5_curves = search_curves(LO, HI, per_bin=2, nbins=8)
print(f"Found {len(e5_curves)} curves, lam* spread over [0, 0.5].\n")

E5_SEEDS = [42, 1234, 9999, 555, 31337, 6, 77, 808]
for eff_t in (0.15, 0.30, 0.45, 0.60):
    print(f"=== eff target = {eff_t}   (m = {M_SIGS}, {len(E5_SEEDS)} seeds) ===")
    print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} {'nu_hat':>7} "
          f"{'A0':>6} {'B1':>6}")
    rows = []
    for (p, b, n, lam, G) in e5_curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        S_K1, _, S_K2, _ = scales(n, k1, k2)
        nh = nu_hat(n, lam, S_K1, S_K2)
        curve = (p, b, n, lam, G)
        wA, tA = success_rate(curve, M_SIGS, k1, E5_SEEDS, 'kannan', False)
        wB, tB = success_rate(curve, M_SIGS, k1, E5_SEEDS, 'cvp', True)
        rows.append({'lam_star': lam_star(lam, n), 'nu': nh,
                     'pA': wA / tA, 'pB': wB / tB})
        print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1:>4} "
              f"{k1 * k2 / n:>6.3f} {nh:>7.3f} "
              f"{str(wA)+'/'+str(tA):>6} {str(wB)+'/'+str(tB):>6}")
    for arm in ('pA', 'pB'):
        nm = 'A0 (uncentered)' if arm == 'pA' else 'B1 (centered CVP)'
        ps = [r[arm] for r in rows]
        print(f"  {nm:<20} mean p_hat = {sum(ps)/len(ps):.3f}   "
              f"spearman(nu_hat,p) = {spearman([r['nu'] for r in rows], ps):+.3f}   "
              f"spearman(lam*,p) = {spearman([r['lam_star'] for r in rows], ps):+.3f}")
    print()

print("=" * 78)
print("done")
print("=" * 78)
