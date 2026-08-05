"""
Thread 23 — GLV-HNP Phase 2: reformulate the lattice so the planted vector is
the target of the search, then ask whether the K1 wall survives.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5): in the Kannan lattice
built by `build_glv_lattice` (glv_hnp_phase2_20bit.py:262) the shortest vector is
always the trivial `n*S_D*e_m`, never the planted vector.  The 2026-07-29 entry
proposed two fixes: (i) quotient out the trivial direction, (ii) replace the
Kannan embedding by an explicit CVP.  This script does both and adds the
decisive measurement the earlier runs could not make.

Experiments
  T23-A  S_D sweep.  Corrects the 2026-07-29 claim that "no choice of S_D
         removes the trivial vector"; derives and verifies the crossover.
  T23-B  BDD reformulation.  d is eliminated algebraically (not given a
         coordinate), giving a 2m-dimensional lattice of determinant n^(m-1).
         Solvers: Babai nearest-plane, and *exact* CVP (fpylll, method="proved").
  T23-C  Is the planted vector the closest lattice point at all?  This is an
         information-theoretic question about the instance, independent of any
         reduction algorithm.  ||v_planted - t|| vs ||v_CVP - t||.
  T23-D  K1-wall head-to-head on the 2026-07-29 T4 grid (the falsifier).

Run: python3 glv_hnp_phase2_bdd.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic + helpers.  Copied verbatim from
# glv_hnp_phase2_20bit.py so the comparison against 2026-07-26 / 2026-07-29 is
# exact (that file runs its experiments at import time, so it is not importable).
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

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
    """Verbatim from glv_hnp_phase2_20bit.py:236."""
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
# Baseline: the Kannan lattice under study (verbatim, glv_hnp_phase2_20bit.py:262)
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    """Verbatim except that S_D is now a parameter (was hard-coded to 1)."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n

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
    return M, S_K1, S_D, S_K2, S_KANNAN

def recover_d(M_reduced, m, n, S_KANNAN, S_D, d_secret=None):
    """Verbatim from glv_hnp_phase2_20bit.py:295, generalised for S_D != 1."""
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % S_D != 0: continue
        d_cand = (val // S_D) % n
        if d_cand == 0: continue
        if d_secret is not None and d_cand == d_secret:
            return d_cand
    return None

def kannan_attack(curve, m, d_secret, k1_bound, seed, S_D=1, beta=None):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, None
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice(sigs, n, lam, k1_bound,
                                                     k2_bound, S_D=S_D)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = recover_d(reduced, m, n, S_KANNAN, S_D, d_secret) is not None
    return ok, (sigs, reduced, S_K1, S_D, S_K2, S_KANNAN)

# ---------------------------------------------------------------------------
# T23-B: the BDD reformulation.  d is eliminated, not embedded.
#
# From  k1_i = A_i + B_i*d - lam*k2_i + c_i*n  set  u_i = k1_i - A_i, v_i = k2_i:
#     u_i + lam*v_i  ==  B_i * d           (mod n)
# Multiply the i=0 relation by B_0^{-1} and substitute d:
#     u_i + lam*v_i  ==  beta_i * (u_0 + lam*v_0)   (mod n),  beta_i = B_i/B_0
# for i = 1..m-1.  That is m-1 congruences on the 2m free integers (u, v), so
# the solution set is a rank-2m lattice L in Z^{2m} of determinant n^(m-1),
# with NO free direction: the trivial vector n*S_D*e_m of the Kannan lattice
# does not exist here because d is no longer a coordinate.
#
# The secret is the lattice point closest to  t = (-A_i*S_K1, 0),  at distance
# ~ n*sqrt(2m/3)  (uncentered) or  ~ n*sqrt(m/6)  (centered).
# ---------------------------------------------------------------------------

def build_bdd_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    Bc = [s['B'] for s in sigs]
    B0inv = modinv(Bc[0], n)
    beta = [Bc[i] * B0inv % n for i in range(m)]     # beta[0] == 1
    dim = 2 * m
    rows = []
    # (a) c_j generators: u_j -> u_j + n   (j = 1..m-1)
    for j in range(1, m):
        r = [0] * dim
        r[j] = n * S_K1
        rows.append(r)
    # (b) u_0 = 1  =>  u_i = beta_i
    r = [0] * dim
    for i in range(m):
        r[i] = (1 if i == 0 else beta[i]) * S_K1
    rows.append(r)
    # (c) v_0 = 1  =>  u_0 = 0, u_i = beta_i*lam (i >= 1)
    r = [0] * dim
    for i in range(1, m):
        r[i] = (beta[i] * lam % n) * S_K1
    r[m] = S_K2
    rows.append(r)
    # (d) v_j = 1 (j >= 1)  =>  u_j = -lam, all other u = 0
    for j in range(1, m):
        r = [0] * dim
        r[j] = ((-lam) % n) * S_K1
        r[m + j] = S_K2
        rows.append(r)
    assert len(rows) == dim
    return rows, S_K1, S_K2, B0inv

def bdd_target(sigs, S_K1, S_K2, k1_bound, k2_bound, center):
    m = len(sigs)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
        if center:
            t[i] += (k1_bound * S_K1) // 2
    if center:
        for i in range(m):
            t[m + i] = (k2_bound * S_K2) // 2
    return t

def planted_vector(sigs, S_K1, S_K2):
    """The lattice point the attacker is after: u_i = k1_i - A_i, v_i = k2_i."""
    m = len(sigs)
    return ([(s['k1'] - s['A']) * S_K1 for s in sigs] +
            [s['k2'] * S_K2 for s in sigs])

def in_lattice(vec, sigs, n, lam, S_K1, S_K2):
    """Check the m-1 congruences directly (independent of the basis)."""
    m = len(sigs)
    Bc = [s['B'] for s in sigs]
    B0inv = modinv(Bc[0], n)
    u = [vec[i] // S_K1 for i in range(m)]
    v = [vec[m + i] // S_K2 for i in range(m)]
    if any(vec[i] % S_K1 for i in range(m)): return False
    if any(vec[m + i] % S_K2 for i in range(m)): return False
    rhs0 = (u[0] + lam * v[0]) % n
    for i in range(1, m):
        beta_i = Bc[i] * B0inv % n
        if (u[i] + lam * v[i] - beta_i * rhs0) % n != 0:
            return False
    return True

def gram_schmidt(B):
    nr, nc = len(B), len(B[0])
    Bs, mu = [], [[0.0] * nr for _ in range(nr)]
    for i in range(nr):
        v = [float(x) for x in B[i]]
        for j in range(i):
            nj = sum(x * x for x in Bs[j])
            if nj == 0.0:
                continue
            mu[i][j] = sum(float(B[i][k]) * Bs[j][k] for k in range(nc)) / nj
            for k in range(nc):
                v[k] -= mu[i][j] * Bs[j][k]
        Bs.append(v)
    return Bs

def babai_nearest_plane(B, t):
    Bs = gram_schmidt(B)
    nr, nc = len(B), len(B[0])
    b = [float(x) for x in t]
    coeffs = [0] * nr
    for i in range(nr - 1, -1, -1):
        nj = sum(x * x for x in Bs[i])
        if nj == 0.0:
            continue
        c = sum(b[k] * Bs[i][k] for k in range(nc)) / nj
        ci = int(round(c))
        coeffs[i] = ci
        for k in range(nc):
            b[k] -= ci * float(B[i][k])
    return [sum(coeffs[i] * B[i][k] for i in range(nr)) for k in range(nc)]

def decode_candidate(w, sigs, n, lam, S_K1, S_K2, k1_bound, k2_bound, B0inv):
    """Oracle-free: turn a lattice point into d, or None if it is inconsistent."""
    m = len(sigs)
    t0 = [-s['A'] * S_K1 for s in sigs] + [0] * m
    e = [w[i] - t0[i] for i in range(2 * m)]
    if any(e[i] % S_K1 for i in range(m)):
        return None
    if any(e[m + i] % S_K2 for i in range(m)):
        return None
    k1 = [e[i] // S_K1 for i in range(m)]
    k2 = [e[m + i] // S_K2 for i in range(m)]
    if any(not (0 <= x < k1_bound) for x in k1):
        return None
    if any(not (0 <= x < k2_bound) for x in k2):
        return None
    d_cand = B0inv * (k1[0] + lam * k2[0] - sigs[0]['A']) % n
    # full self-check against every signature (no knowledge of the secret)
    for i in range(m):
        if (sigs[i]['A'] + sigs[i]['B'] * d_cand - k1[i] - lam * k2[i]) % n != 0:
            return None
    return d_cand

def norm2(v):
    return sum(int(x) * int(x) for x in v)

def bdd_attack(curve, m, d_secret, k1_bound, seed, center=True, exact=True):
    """Returns dict with per-solver outcomes and the T23-C distance diagnostic."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    rows, S_K1, S_K2, B0inv = build_bdd_lattice(sigs, n, lam, k1_bound, k2_bound)
    t = bdd_target(sigs, S_K1, S_K2, k1_bound, k2_bound, center)
    vp = planted_vector(sigs, S_K1, S_K2)
    assert in_lattice(vp, sigs, n, lam, S_K1, S_K2), "planted vector not in L"

    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    red = [[A[i][j] for j in range(2 * m)] for i in range(2 * m)]

    out = {'n': n, 'm': m, 'K1': k1_bound, 'K2': k2_bound,
           'd_planted': math.sqrt(norm2([vp[i] - t[i] for i in range(2 * m)]))}

    wb = babai_nearest_plane(red, t)
    db = decode_candidate(wb, sigs, n, lam, S_K1, S_K2, k1_bound, k2_bound, B0inv)
    out['babai'] = (db is not None and db == d_secret)
    out['babai_selfok'] = db is not None
    out['d_babai'] = math.sqrt(norm2([wb[i] - t[i] for i in range(2 * m)]))

    if exact:
        wc = list(CVP.closest_vector(A, [int(x) for x in t], method="proved"))
        dc = decode_candidate(wc, sigs, n, lam, S_K1, S_K2, k1_bound, k2_bound, B0inv)
        out['cvp'] = (dc is not None and dc == d_secret)
        out['cvp_selfok'] = dc is not None
        out['d_cvp'] = math.sqrt(norm2([wc[i] - t[i] for i in range(2 * m)]))
        out['planted_is_closest'] = out['d_cvp'] >= out['d_planted'] - 1e-9
    return out

# ===========================================================================
# Setup: the historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 /
# 2026-07-26 / 2026-07-29).  lam* = min(lam, n-lam)/n is the root-convention-
# free eigenvalue statistic fixed by the 2026-07-29 entry (exp T1).
# ===========================================================================

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

def lam_star(lam, n):
    return min(lam, n - lam) / n

curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    curves.append((label, (p, b, n, lam, G), k1, m))

print("=" * 78)
print("Thread 23 — reformulating the Phase-2 lattice so the planted vector")
print("            is the search target (BDD/CVP instead of Kannan SVP)")
print("=" * 78)
for label, cur, k1, m in curves:
    p, b, n, lam, G = cur
    print(f"  {label:18s} p={p:5d} n={n:5d} lam={lam:5d} "
          f"lam*={lam_star(lam, n):.4f} K1={k1} m={m}")

# ---------------------------------------------------------------------------
# T23-A — the S_D claim from 2026-07-29 is quantitatively wrong
# ---------------------------------------------------------------------------

print("\n" + "-" * 78)
print("T23-A: does S_D scaling remove the trivial vector n*S_D*e_m?")
print("-" * 78)
print("2026-07-29 T5 claimed: 'no choice of S_D removes it - both vectors scale")
print("linearly in S_D.'  Only the trivial vector scales fully:")
print("  ||n*S_D*e_m||^2 = n^2 * S_D^2")
print("  ||v_planted||^2 ~ n^2 * (2m/3 + S_D^2/3 + 1)   (only ONE coord has S_D)")
print("  => trivial is longer once S_D^2/1 > 2m/3 + S_D^2/3 + 1, i.e. S_D > sqrt(m+3/2)")
print("The cost is only det^(1/dim): S_D^(1/(2m+2)).\n")

print(f"{'curve':18s} {'m':>3s} {'S_D*':>6s} | "
      + " ".join(f"S_D={s:<4d}" for s in [1, 2, 4, 8, 16, 64]))
print(f"{'':18s} {'':>3s} {'':>6s} | " + " ".join(f"{'w/5':8s}" for _ in range(6)))
sd_grid = [1, 2, 4, 8, 16, 64]
sd_table = []
for label, cur, k1, m in curves:
    p, b, n, lam, G = cur
    sd_star = math.sqrt(m + 1.5)
    cells = []
    for S_D in sd_grid:
        wins = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok, _ = kannan_attack(cur, m, d_trial, k1, seed, S_D=S_D)
            wins += bool(ok)
        cells.append(wins)
    sd_table.append((label, m, sd_star, cells))
    print(f"{label:18s} {m:3d} {sd_star:6.2f} | "
          + " ".join(f"{c}/5{'':5s}" for c in cells))

# Confirm the crossover empirically: is n*S_D*e_m still the shortest vector?
print("\nShortest-vector identity after LLL (is it the trivial n*S_D*e_m?):")
print(f"{'curve':18s} " + " ".join(f"S_D={s:<4d}" for s in sd_grid))
for label, cur, k1, m in curves:
    p, b, n, lam, G = cur
    marks = []
    for S_D in sd_grid:
        d_trial = random.Random(42 + 7777).randint(1, n - 1)
        ok, aux = kannan_attack(cur, m, d_trial, k1, 42, S_D=S_D)
        sigs, reduced, S_K1, S_Dr, S_K2, S_KANNAN = aux
        dim = 2 * m + 2
        sv = min(reduced, key=lambda r: norm2(r))
        triv = (abs(sv[m]) == n * S_D and
                all(sv[j] == 0 for j in range(dim) if j != m))
        marks.append("triv" if triv else "othr")
    print(f"{label:18s} " + " ".join(f"{x:<8s}" for x in marks))

# ---------------------------------------------------------------------------
# T23-B / T23-C — BDD reformulation and the closest-point diagnostic
# ---------------------------------------------------------------------------

print("\n" + "-" * 78)
print("T23-B/C: BDD reformulation (d eliminated; dim 2m, det n^(m-1))")
print("-" * 78)
print("Solvers: Babai nearest-plane on the LLL basis, and EXACT CVP")
print("(fpylll method='proved').  Exact CVP is the ceiling of this formulation:")
print("if it fails, no amount of reduction strength can help.")
print("'closest' = is the planted vector the true closest lattice point to t?\n")

hdr = (f"{'curve':18s} {'ctr':>4s} {'sd':>6s} {'kannan':>7s} {'babai':>7s} "
       f"{'cvp':>7s} {'closest':>8s} {'d_pl/d_cvp':>11s}")
print(hdr)
print("-" * len(hdr))

bdd_rows = []
for label, cur, k1, m in curves:
    p, b, n, lam, G = cur
    for center in (False, True):
        kan = bab = cvp = clo = 0
        ratios = []
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            okk, _ = kannan_attack(cur, m, d_trial, k1, seed)
            kan += bool(okk)
            r = bdd_attack(cur, m, d_trial, k1, seed, center=center, exact=True)
            bab += bool(r['babai'])
            cvp += bool(r['cvp'])
            clo += bool(r['planted_is_closest'])
            ratios.append(r['d_planted'] / r['d_cvp'] if r['d_cvp'] else float('inf'))
        rr = sum(ratios) / len(ratios)
        bdd_rows.append((label, center, kan, bab, cvp, clo, rr))
        print(f"{label:18s} {str(center):>4s} {'':>6s} {kan}/5{'':4s} "
              f"{bab}/5{'':4s} {cvp}/5{'':4s} {clo}/5{'':5s} {rr:11.4f}")

# ---------------------------------------------------------------------------
# T23-D — the falsifier: does the K1 wall move?
# ---------------------------------------------------------------------------

print("\n" + "-" * 78)
print("T23-D: K1 wall, Kannan vs BDD-exact-CVP  (2026-07-29 T4 grid)")
print("-" * 78)
print("2026-07-29 T4: 12-bit/2557 wall at K1~12-16; 12-bit/2677 wall at K1~4-6.")
print("Falsifier: if the wall moves outward under exact CVP the reformulation")
print("is a real improvement; if it does not, the wall is information-theoretic.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
wall_rows = []
for label, cur, k1_hist, m in curves[1:]:      # the two 12-bit T4 curves
    p, b, n, lam, G = cur
    print(f"{label}  (m={m}, lam*={lam_star(lam, n):.3f})")
    print(f"  {'K1':>4s} " + " ".join(f"{k:>6d}" for k in K1_GRID))
    for arm in ("kannan", "bdd-cvp"):
        cells = []
        for K1 in K1_GRID:
            wins = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                if arm == "kannan":
                    ok, _ = kannan_attack(cur, m, d_trial, K1, seed)
                    wins += bool(ok)
                else:
                    r = bdd_attack(cur, m, d_trial, K1, seed,
                                   center=True, exact=True)
                    wins += bool(r['cvp'])
            cells.append(wins)
        wall_rows.append((label, arm, cells))
        print(f"  {arm:>7s} " + " ".join(f"{c:>4d}/5" for c in cells))
    # closest-point diagnostic across the same grid
    cells = []
    for K1 in K1_GRID:
        c = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = bdd_attack(cur, m, d_trial, K1, seed, center=True, exact=True)
            c += bool(r['planted_is_closest'])
        cells.append(c)
    wall_rows.append((label, "closest", cells))
    print(f"  {'closest':>7s} " + " ".join(f"{c:>4d}/5" for c in cells))
    print()

# ---------------------------------------------------------------------------
# T23-E — CONTROL: is the gain from eliminating d, or just from centering?
# The BDD arm changed two things at once.  Centering can be applied to the
# Kannan embedding too: shift the last row by half the box in every column
# (including the d column, whose planted entry is d*S_D with d ~ U(0,n)).
# ---------------------------------------------------------------------------

def build_glv_lattice_centered(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    S_KANNAN = n
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    # NB: use one explicit half and negate it.  Writing -(x) // 2 inline floors
    # the *negative*, which desynchronises the offset added back in recovery.
    half_d = (n * S_D) // 2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1 - (k1_bound * S_K1) // 2
    M[2 * m + 1][m] = -half_d
    for i in range(m):
        M[2 * m + 1][m + 1 + i] = -((k2_bound * S_K2) // 2)
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, S_K1, S_D, S_K2, S_KANNAN

def recover_d_centered(M_reduced, m, n, S_KANNAN, S_D, d_secret=None):
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m] + (n * S_D) // 2
        if val % S_D != 0: continue
        d_cand = (val // S_D) % n
        if d_cand == 0: continue
        if d_secret is not None and d_cand == d_secret:
            return d_cand
    return None

def kannan_attack_centered(curve, m, d_secret, k1_bound, seed, S_D=1):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M, S_K1, S_D, S_K2, S_KANNAN = build_glv_lattice_centered(
        sigs, n, lam, k1_bound, k2_bound, S_D=S_D)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d_centered(reduced, m, n, S_KANNAN, S_D, d_secret) is not None

print("\n" + "-" * 78)
print("T23-E CONTROL: centering alone vs. d-elimination + centering")
print("-" * 78)
print("The BDD arm changed two things at once.  'kan-ctr' is the ORIGINAL Kannan")
print("embedding with only the target centered -- it isolates the centering effect.\n")

ctrl_rows = []
for label, cur, k1_hist, m in curves[1:]:
    p, b, n, lam, G = cur
    print(f"{label}  (m={m}, lam*={lam_star(lam, n):.3f})")
    print(f"  {'K1':>8s} " + " ".join(f"{k:>6d}" for k in K1_GRID))
    for arm in ("kannan", "kan-ctr", "bdd-cvp"):
        cells = []
        for K1 in K1_GRID:
            wins = 0
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                if arm == "kannan":
                    ok, _ = kannan_attack(cur, m, d_trial, K1, seed)
                    wins += bool(ok)
                elif arm == "kan-ctr":
                    wins += bool(kannan_attack_centered(cur, m, d_trial, K1, seed))
                else:
                    r = bdd_attack(cur, m, d_trial, K1, seed, center=True, exact=True)
                    wins += bool(r['cvp'])
            cells.append(wins)
        ctrl_rows.append((label, arm, cells))
        print(f"  {arm:>8s} " + " ".join(f"{c:>4d}/5" for c in cells))
    print()

# ---------------------------------------------------------------------------
# T23-F — list decoding.  On 12-bit/2557 exact CVP under-performs the Kannan
# row-scan, because the Kannan procedure tests ~dim candidate rows while CVP
# returns only the single closest point.  Give BDD the same courtesy: scan the
# neighbourhood w +/- b_i (+/- b_j) of the CVP solution, self-checking each.
# ---------------------------------------------------------------------------

def bdd_attack_list(curve, m, d_secret, k1_bound, seed, depth=2):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return False, 0
    rows, S_K1, S_K2, B0inv = build_bdd_lattice(sigs, n, lam, k1_bound, k2_bound)
    t = bdd_target(sigs, S_K1, S_K2, k1_bound, k2_bound, True)
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    dim = 2 * m
    red = [[A[i][j] for j in range(dim)] for i in range(dim)]
    w0 = list(CVP.closest_vector(A, [int(x) for x in t], method="proved"))

    def try_vec(w):
        return decode_candidate(w, sigs, n, lam, S_K1, S_K2,
                                k1_bound, k2_bound, B0inv)

    tested = 0
    d = try_vec(w0); tested += 1
    if d is not None:
        return d == d_secret, tested
    offs = [[0] * dim]
    for i in range(dim):
        for s in (1, -1):
            offs.append([s * red[i][k] for k in range(dim)])
    if depth >= 2:
        for i in range(dim):
            for j in range(i + 1, dim):
                for si in (1, -1):
                    for sj in (1, -1):
                        offs.append([si * red[i][k] + sj * red[j][k]
                                     for k in range(dim)])
    for o in offs:
        d = try_vec([w0[k] + o[k] for k in range(dim)]); tested += 1
        if d is not None:
            return d == d_secret, tested
    return False, tested

print("-" * 78)
print("T23-F: list decoding -- CVP solution + neighbourhood scan (self-checked,")
print("       no secret oracle).  Matches the candidate budget the Kannan")
print("       row-scan already spends.")
print("-" * 78)
for label, cur, k1_hist, m in curves[1:]:
    p, b, n, lam, G = cur
    print(f"{label}  (m={m})")
    print(f"  {'K1':>8s} " + " ".join(f"{k:>6d}" for k in K1_GRID))
    cells, costs = [], []
    for K1 in K1_GRID:
        wins, tot = 0, 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok, tested = bdd_attack_list(cur, m, d_trial, K1, seed)
            wins += bool(ok); tot += tested
        cells.append(wins); costs.append(tot // len(SEEDS))
    print(f"  {'bdd-list':>8s} " + " ".join(f"{c:>4d}/5" for c in cells))
    print(f"  {'cands':>8s} " + " ".join(f"{c:>6d}" for c in costs))
    print()

# ---------------------------------------------------------------------------
# T23-G — does centering generalise?  Replicates the 2026-07-29 T3 sweep
# (20 fresh 17-bit j=0 GLV curves, m=12, 5 seeds, eff = K1*K2/n) with the
# centered target.  T3 reported 19/20, 3/20, 0/20 at eff = 0.05, 0.15, 0.25.
# ---------------------------------------------------------------------------

import sympy

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
                bb = num // 2
                if bb >= 0 and a * a - a * bb + bb * bb == p:
                    return (a, bb)
    return None

def j0_traces(a, b):
    return [2 * a - b, -2 * a + b, -(a + b), a + b, 2 * b - a, a - 2 * b]

def glv_lambda(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return min(r1, n - 1 - r1)

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        y = tonelli_shanks((pow(x, 3, p) + b) % p, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n)
            if G is not None:
                return (p, b, n, glv_lambda(n), G)
    return None

def sweep_curves(lo, hi, want=20):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_c = p + 1 - t
                    if n_c < 2 or n_c % 3 != 1 or not sympy.isprime(n_c):
                        continue
                    lam = glv_lambda(n_c)
                    if lam is None:
                        continue
                    cur = build_curve(p, n_c)
                    if cur is not None and cur[3] is not None:
                        out.append(cur)
                        break
        p = int(sympy.nextprime(p))
    return out

print("-" * 78)
print("T23-G: does centering generalise?  17-bit sweep, m=12, 5 seeds.")
print("       (replicates 2026-07-29 T3: 19/20, 3/20, 0/20 at eff .05/.15/.25)")
print("-" * 78)

sweep = sweep_curves(2**16, 2**17, want=20)
M_SWEEP = 12
print(f"collected {len(sweep)} fresh 17-bit j=0 GLV curves\n")
print(f"{'eff':>5s} | {'kannan 5/5':>11s} {'kan-ctr 5/5':>12s} | "
      f"{'kannan tot':>11s} {'kan-ctr tot':>12s}")
sweep_rows = []
for eff in (0.05, 0.15, 0.25, 0.35):
    kan_full = ctr_full = 0
    kan_tot = ctr_tot = 0
    for cur in sweep:
        p, b, n, lam, G = cur
        K2 = math.isqrt(n) + 1
        K1 = max(2, round(eff * n / K2))
        kw = cw = 0
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok, _ = kannan_attack(cur, M_SWEEP, d_trial, K1, seed)
            kw += bool(ok)
            cw += bool(kannan_attack_centered(cur, M_SWEEP, d_trial, K1, seed))
        kan_full += (kw == len(SEEDS)); ctr_full += (cw == len(SEEDS))
        kan_tot += kw; ctr_tot += cw
    denom = len(sweep) * len(SEEDS)
    sweep_rows.append((eff, kan_full, ctr_full, kan_tot, ctr_tot, len(sweep), denom))
    print(f"{eff:5.2f} | {kan_full:>4d}/{len(sweep):<6d} {ctr_full:>5d}/{len(sweep):<6d} | "
          f"{kan_tot:>4d}/{denom:<6d} {ctr_tot:>5d}/{denom:<6d}")

print("\n" + "-" * 78)
print("T23-H: seed-robustness of the headline cells (20 seeds, not 5)")
print("-" * 78)
BIG_SEEDS = [42, 1234, 9999, 555, 31337, 7, 13, 101, 2024, 88,
             314, 271, 1618, 4242, 999, 5150, 60, 77, 123, 456]
print(f"{'curve':18s} {'K1':>4s} {'kannan':>9s} {'kan-ctr':>9s} {'bdd-cvp':>9s}")
for label, cur, k1_hist, m in curves[1:]:
    p, b, n, lam, G = cur
    for K1 in (6, 8, 12, 16):
        kw = cw = bw = 0
        for seed in BIG_SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            ok, _ = kannan_attack(cur, m, d_trial, K1, seed)
            kw += bool(ok)
            cw += bool(kannan_attack_centered(cur, m, d_trial, K1, seed))
            r = bdd_attack(cur, m, d_trial, K1, seed, center=True, exact=True)
            bw += bool(r['cvp'])
        nb = len(BIG_SEEDS)
        print(f"{label:18s} {K1:>4d} {kw:>5d}/{nb:<3d} {cw:>5d}/{nb:<3d} {bw:>5d}/{nb:<3d}")

print("=" * 78)
print("done")
print("=" * 78)
