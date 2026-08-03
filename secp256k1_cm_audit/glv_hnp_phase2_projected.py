"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is lambda_1.

Background (2026-07-29 log, T5): in the Phase-2 lattice of
`glv_hnp_phase2_20bit.py::build_glv_lattice` the shortest vector after LLL is ALWAYS the
trivial vector n*S_D*e_m (100% of its energy in the d-column, Kannan coord 0). It carries
no information (d is only defined mod n) and no choice of S_D removes it, because both it
and the planted vector scale linearly in S_D. Recovery was therefore a BDD/coset condition
("shortest among vectors with last coord +-S_KANNAN"), not an SVP condition.

This script tests two independent fixes, in a 2x2 design:

  PROJECTION  -- delete the d-column. The trivial vector n*S_D*e_m maps to 0, so it is
                 gone from the lattice entirely. d is recovered afterwards in closed form
                 from the k1/k2 blocks:  d = (k1_0 + lam*k2_0 - A_0) * B_0^-1  (mod n).
                 Basis care: {B*S_K1} u {n*S_K1*e_i, i=0..m-1} is m+1 generators of a
                 rank-m lattice. Normalising B' = B*B_0^-1 (mod n), so B'_0 = 1, the set
                 {B'*S_K1} u {n*S_K1*e_i, i=1..m-1} is an honest basis (det S_K1^m n^(m-1)),
                 and n*e_0 = n*B' - sum_{i>=1} n*B'_i*e_i stays in the span.

  CENTERING   -- k1_i in [0,K1) has E[(k1_i*S_K1)^2] ~ n^2/3. Substituting
                 k1_i = k1'_i + K1//2, k2_i = k2'_i + K2//2 (so k1'_i in [-K1/2, K1/2))
                 costs only a shift of the A-row,
                     A'_i = A_i - K1//2 - lam*(K2//2)  (mod n),
                 and drops the per-coordinate second moment to n^2/12. Expected:
                 ||pv||^2 : n^2*(2m/3 + 4/3) -> n^2*(m/6 + 1)  (projected+centered).

Variants: V0 = baseline (verbatim from glv_hnp_phase2_20bit.py), V1 = baseline+centered,
V2 = projected, V3 = projected+centered.

Experiments:
  E1  membership/determinant sanity check on the projected basis
  E2  sv/pv (is the planted vector lambda_1?) for all four variants
  E3  the T4 K1-grid, replicated across all four variants -- does the K1 wall move?
  E4  Gaussian-heuristic eff-ceiling prediction vs. empirical success

Run: python3 glv_hnp_phase2_projected.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- copied verbatim from
# glv_hnp_phase2_20bit.py so the comparison to 2026-07-26/07-29 is exact.
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
        if y is not None and y != 0:
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
        for num in [a + s, a - s]:
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_eigenvalue(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    assert (r1 * r1 + r1 + 1) % n == 0, f"GLV eigenvalue check failed for n={n}"
    lam = min(r1, r2)
    return lam, n - 1 - lam

# ---------------------------------------------------------------------------
# Signature generation (verbatim)
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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

# ---------------------------------------------------------------------------
# V0/V1: baseline lattice, dim 2m+2 (verbatim build_glv_lattice + centering hook)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def shifted_A(sig, n, lam, k1_bound, k2_bound, centered):
    """A'_i = A_i - K1//2 - lam*(K2//2) mod n when centering, else A_i."""
    if not centered:
        return sig['A']
    return (sig['A'] - k1_bound // 2 - lam * (k2_bound // 2)) % n

def build_baseline(sigs, n, lam, k1_bound, k2_bound, centered=False):
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
        M[2 * m + 1][i] = shifted_A(sigs[i], n, lam, k1_bound, k2_bound, centered) * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M

def planted_baseline(sigs, n, lam, k1_bound, k2_bound, d_secret, centered):
    """The vector LLL is supposed to find, in the dim-(2m+2) baseline lattice."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    off1 = k1_bound // 2 if centered else 0
    off2 = k2_bound // 2 if centered else 0
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - off1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - off2) * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KANNAN
    return v

def recover_baseline(row, m, n, lam, sigs, k1_bound, k2_bound, S, centered):
    S_K1, S_D, S_K2, S_KANNAN = S
    if abs(row[2 * m + 1]) != S_KANNAN:
        return None
    r = row if row[2 * m + 1] > 0 else [-x for x in row]
    return (r[m] * modinv(S_D, n)) % n

# ---------------------------------------------------------------------------
# V2/V3: projected lattice, dim 2m+1 (d-column deleted)
# ---------------------------------------------------------------------------

def build_projected(sigs, n, lam, k1_bound, k2_bound, centered=False):
    """
    Columns: 0..m-1 = k1 block, m..2m-1 = k2 block, 2m = Kannan.
    Rows (2m+1, an honest basis -- see module docstring):
        B'*S_K1                          (B'_0 = 1)
        n*S_K1*e_i,           i=1..m-1
        -lam*S_K1*e_i + S_K2*e_{m+i},    i=0..m-1
        A'_i*S_K1 (+ S_KANNAN at col 2m)
    """
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    Bp = [sigs[i]['B'] * B0inv % n for i in range(m)]
    assert Bp[0] == 1

    M = []
    row = [0] * dim
    for i in range(m):
        row[i] = Bp[i] * S_K1
    M.append(row)
    for i in range(1, m):
        row = [0] * dim
        row[i] = n * S_K1
        M.append(row)
    for i in range(m):
        row = [0] * dim
        row[i] = -lam * S_K1
        row[m + i] = S_K2
        M.append(row)
    row = [0] * dim
    for i in range(m):
        row[i] = shifted_A(sigs[i], n, lam, k1_bound, k2_bound, centered) * S_K1
    row[2 * m] = S_KANNAN
    M.append(row)
    assert len(M) == dim
    return M

def planted_projected(sigs, n, lam, k1_bound, k2_bound, d_secret, centered):
    """d_secret is unused: the projection deletes the d-coordinate."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    off1 = k1_bound // 2 if centered else 0
    off2 = k2_bound // 2 if centered else 0
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - off1) * S_K1
        v[m + i] = (sigs[i]['k2'] - off2) * S_K2
    v[2 * m] = S_KANNAN
    return v

def recover_projected(row, m, n, lam, sigs, k1_bound, k2_bound, S, centered):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^-1 (mod n), read off the k1/k2 blocks."""
    S_K1, S_D, S_K2, S_KANNAN = S
    if abs(row[2 * m]) != S_KANNAN:
        return None
    r = row if row[2 * m] > 0 else [-x for x in row]
    if r[0] % S_K1 != 0 or r[m] % S_K2 != 0:
        return None
    off1 = k1_bound // 2 if centered else 0
    off2 = k2_bound // 2 if centered else 0
    k1_0 = r[0] // S_K1 + off1
    k2_0 = r[m] // S_K2 + off2
    return (k1_0 + lam * k2_0 - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n

# ---------------------------------------------------------------------------
# Variant dispatch
# ---------------------------------------------------------------------------

VARIANTS = {
    'V0 base':      (build_baseline,  planted_baseline,  recover_baseline,  False),
    'V1 base+ctr':  (build_baseline,  planted_baseline,  recover_baseline,  True),
    'V2 proj':      (build_projected, planted_projected, recover_projected, False),
    'V3 proj+ctr':  (build_projected, planted_projected, recover_projected, True),
}

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def run_one(variant, curve_params, m, d_secret, k1_bound, seed):
    """Returns dict with recovered flag, sv/pv ratio, dims."""
    p, b, n, lam, G = curve_params
    build, planted, recover, centered = VARIANTS[variant]
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S = scales(n, k1_bound, k2_bound)
    M = build(sigs, n, lam, k1_bound, k2_bound, centered)
    dim = len(M)
    pv = planted(sigs, n, lam, k1_bound, k2_bound, d_secret, centered)

    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]

    ok = any(recover(r, m, n, lam, sigs, k1_bound, k2_bound, S, centered) == d_secret
             for r in reduced)
    nz = [r for r in reduced if any(r)]
    sv = min(nz, key=norm) if nz else [0] * dim
    return {'ok': ok, 'sv': norm(sv), 'pv': norm(pv), 'dim': dim,
            'ratio': norm(sv) / norm(pv) if norm(pv) else float('inf')}

def sweep(variant, curve_params, m, k1_bound, seeds):
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, curve_params[2] - 1)
        res = run_one(variant, curve_params, m, d_trial, k1_bound, seed)
        if res is None:
            continue
        wins += res['ok']
        ratios.append(res['ratio'])
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))

# ---------------------------------------------------------------------------
# Reference curves (identical to 2026-07-26 / 2026-07-29)
# ---------------------------------------------------------------------------

def make_curve(p, b, n):
    lam, _ = glv_eigenvalue(n)
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for p={p}"
    return (p, b, n, lam, G)

print("=" * 78)
print("GLV-HNP Phase 2 -- Thread 23: projected + centered lattice reformulation")
print("=" * 78)

C_8    = make_curve(211, 2, 199)
C_2557 = make_curve(2557, 2, 2659)
C_2677 = make_curve(2677, 2, 2647)
REFS = [("8-bit/199", C_8, 2), ("12-bit/2557", C_2557, 8), ("12-bit/2677", C_2677, 8)]
for label, c, _ in REFS:
    lam = c[3]
    print(f"  {label:14s} p={c[0]:5d} n={c[2]:5d} lam={lam:5d} "
          f"lam*={min(lam, c[2]-lam)/c[2]:.4f}")

# ---- E1: membership + determinant sanity ----------------------------------
print("\n" + "=" * 78)
print("E1: projected-basis sanity (membership of planted vector, determinant)")
print("=" * 78)

def lattice_det(M):
    A = IntegerMatrix.from_matrix([r[:] for r in M])
    LLL.reduction(A)
    dim = len(M)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    nz = [r for r in rows if any(r)]
    return len(nz), sympy.Matrix(nz).det() if len(nz) == dim else None

for label, c, k1b in REFS:
    p, b, n, lam, G = c
    m = 6
    k2b = math.isqrt(n) + 1
    d = random.Random(99).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, b, k1b, k2b, 42)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
    for centered in (False, True):
        M = build_projected(sigs, n, lam, k1b, k2b, centered)
        pv = planted_projected(sigs, n, lam, k1b, k2b, d, centered)
        # membership: solve M^T x = pv over Q, check x integral
        x = sympy.Matrix(M).T.solve(sympy.Matrix(pv))
        integral = all(v.q == 1 for v in x)
        rank, det = lattice_det(M)
        pred = S_K1**m * n**(m - 1) * S_K2**m * S_KAN
        print(f"  {label:14s} ctr={int(centered)}  rank={rank}/{2*m+1}  "
              f"planted in L: {integral}  det==S_K1^m*n^(m-1)*S_K2^m*S_KAN: "
              f"{abs(det) == pred if det is not None else 'n/a'}")

# ---- E2: is the planted vector lambda_1? ----------------------------------
print("\n" + "=" * 78)
print("E2: sv/pv  (>1.0 means the planted vector is at or below lambda_1)")
print("=" * 78)
print(f"  {'curve':14s} {'K1':>3s} {'m':>3s} " + " ".join(f"{v:>13s}" for v in VARIANTS))
E2_ROWS = [("8-bit/199", C_8, 2, 6), ("12-bit/2557", C_2557, 8, 8),
           ("12-bit/2677", C_2677, 8, 8), ("12-bit/2677", C_2677, 4, 8)]
for label, c, k1b, m in E2_ROWS:
    cells = []
    for v in VARIANTS:
        _, _, ratio = sweep(v, c, m, k1b, [42, 1234, 9999])
        cells.append(f"{ratio:13.3f}")
    print(f"  {label:14s} {k1b:3d} {m:3d} " + " ".join(cells))

# ---- E3: the T4 K1 grid, all variants -------------------------------------
print("\n" + "=" * 78)
print("E3: T4 K1-wall replication (m=12, 5 seeds).  2026-07-29 baseline in brackets.")
print("=" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
SEEDS5 = [42, 1234, 9999, 31337, 271828]
T4_PRIOR = {
    "12-bit/2557": {2: 5, 3: 5, 4: 5, 6: 5, 8: 5, 12: 4, 16: 1, 24: 0},
    "12-bit/2677": {2: 5, 3: 5, 4: 5, 6: 2, 8: 0, 12: 0, 16: 0, 24: 0},
}
for label, c in [("12-bit/2557", C_2557), ("12-bit/2677", C_2677)]:
    print(f"\n  {label}  (lam*={min(c[3], c[2]-c[3])/c[2]:.3f})")
    print("    variant        " + " ".join(f"K1={k:<3d}" for k in K1_GRID))
    prior = " ".join(f"{T4_PRIOR[label].get(k, '-'):>3}   " for k in K1_GRID)
    print(f"    {'[2026-07-29]':15s}" + prior)
    walls = {}
    for v in VARIANTS:
        cells, wall = [], None
        for k1b in K1_GRID:
            w, tot, _ = sweep(v, c, 12, k1b, SEEDS5)
            cells.append(f"{w}/{tot}  ")
            if w == tot:
                wall = k1b
        walls[v] = wall
        print(f"    {v:15s}" + " ".join(cells))
    print("    last K1 with 5/5: " +
          ", ".join(f"{v}={walls[v]}" for v in VARIANTS))

# ---- E4: Gaussian-heuristic eff ceiling ------------------------------------
print("\n" + "=" * 78)
print("E4: Gaussian-heuristic prediction for the projected+centered lattice")
print("=" * 78)
print("""  Projected lattice: dim 2m+1, det = S_K1^m * n^(m-1) * S_K2^m * S_KANNAN.
  With S_K1=n/K1, S_K2=n/K2, S_KANNAN=n and eff=K1*K2/n:
      det = n^(3m) / (K1*K2)^m = n^(2m) / eff^m
  The competitor is the last-coord-0 sublattice (rank 2m, det = det/S_KANNAN):
      GH_0 = sqrt(2m/(2*pi*e)) * (n^(2m-1)/eff^m)^(1/(2m))
           = sqrt(2m/(2*pi*e)) * n^(1-1/(2m)) / sqrt(eff)
  Planted norm (centered):  ||pv|| = n*sqrt(2m/12 + 1) = n*sqrt(m/6 + 1)
  Recovery predicted when GH_0 > ||pv||, i.e.
      sqrt(eff) < sqrt(2m/(2*pi*e)) / (sqrt(m/6+1) * n^(1/(2m)))
  As m -> infinity this tends to eff < 3/(pi*e) = 0.351 (centered),
  and to eff < 3/(4*pi*e) = 0.0878 (uncentered, ||pv||=n*sqrt(2m/3)).

  CAUTION: those are asymptotes in m and are NOT attained at m=12 -- the n^(1/m)
  factor is still 1.9-2.5 there. Use the finite-m form
      eff < [2m/(2*pi*e)] / [C(m) * n^(1/m)],  C = m/6+1 centered, 2m/3+1 uncentered,
  which is what glv_hnp_phase2_eff_ceiling.py tabulates against observation.
""")

def gh_ratio(n, m, k1b, centered):
    k2b = math.isqrt(n) + 1
    eff = k1b * k2b / n
    gh0 = math.sqrt(2 * m / (2 * math.pi * math.e)) * n ** (1 - 1 / (2 * m)) / math.sqrt(eff)
    pv = n * math.sqrt(m / 6 + 1) if centered else n * math.sqrt(2 * m / 3 + 1)
    return eff, gh0 / pv

print(f"  {'curve':14s} {'K1':>3s} {'eff':>7s} {'GH0/pv V2':>10s} {'V2 emp':>8s}"
      f" {'GH0/pv V3':>10s} {'V3 emp':>8s}")
for label, c in [("12-bit/2557", C_2557), ("12-bit/2677", C_2677)]:
    for k1b in [2, 4, 8, 16, 32, 48]:
        eff, r_u = gh_ratio(c[2], 12, k1b, False)
        _, r_c = gh_ratio(c[2], 12, k1b, True)
        w_u, t, _ = sweep('V2 proj', c, 12, k1b, SEEDS5)
        w_c, _, _ = sweep('V3 proj+ctr', c, 12, k1b, SEEDS5)
        print(f"  {label:14s} {k1b:3d} {eff:7.3f} {r_u:10.3f} {w_u:>4d}/{t:<3d}"
              f" {r_c:10.3f} {w_c:>4d}/{t:<3d}")

print("\n" + "=" * 78)
print("DONE")
print("=" * 78)
