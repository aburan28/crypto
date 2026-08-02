#!/usr/bin/env python3
"""
Thread 23 (autolab 2026-08-02) -- Phase-2 GLV-HNP lattice reformulation.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, finding T5):
  The Phase-2 lattice of `glv_hnp_phase2_20bit.py:262` (`build_glv_lattice`)
  carries an explicit d-coordinate scaled by S_D = 1.  Because
      n*row_d - sum_i B_i*row_i  =  (0,...,0, n*S_D, 0,...,0)
  is a lattice vector of norm exactly n*S_D, and because
      ||v_planted|| ~ n*sqrt(2m/3 + 4/3)  >  n*S_D  for every m >= 1,
  the planted vector is NEVER lambda_1.  Every recovery in Phase 2 has been a
  BDD/coset event, not an SVP event.  T5 also showed the trivial vector cannot
  be scaled away: both it and the planted vector's d-component are linear in S_D.

This script tests the obvious structural fix -- DELETE the d-column -- and two
further reformulations, and measures whether the K1 wall of T4 moves.

Variants
--------
V1   baseline               dim (2m+2) x (2m+2), d-column, Kannan embedding.
V2   d-column dropped       (2m+2) generators in Z^(2m+1); rank 2m+1, so the
                            generating set is rank-deficient by exactly 1 --
                            the trivial vector n*e_d IS that deficiency and
                            LLL returns it as a zero row.
V2c  V2 + centering         Kannan row shifted by -(K1/2)*S_K1 and -(K2/2)*S_K2
                            so residuals lie in [-X/2, X/2); shrinks the planted
                            vector by ~sqrt(2).
V3   V2c minus Kannan       Babai nearest-plane CVP in Z^(2m) against the
                            centered target; no embedding coordinate at all.

d is recovered a posteriori from the k1_i, k2_i:  d = B_i^{-1}(k1_i + lam*k2_i - A_i).

Usage:  python3 glv_hnp_phase2_dcol.py
"""

import math
import random
import sys
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py:87)
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

# ---------------------------------------------------------------------------
# Signatures (verbatim from glv_hnp_phase2_lambda_threshold.py:230)
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- identical to the Phase-2 convention."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# V1 -- baseline (verbatim)
# ---------------------------------------------------------------------------
# columns:  0..m-1 = k1 residuals | m = d | m+1..2m = k2 | 2m+1 = Kannan

def build_v1(sigs, n, lam, K1, K2, centered=False):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    c1 = (K1 // 2) * S_K1 if centered else 0
    c2 = (K2 // 2) * S_K2 if centered else 0
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1 - c1
        M[2 * m + 1][m + 1 + i] = -c2
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def planted_v1(sigs, d, n, K1, K2, centered=False):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    c1 = (K1 // 2) if centered else 0
    c2 = (K2 // 2) if centered else 0
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + 1 + i] = (sigs[i]['k2'] - c2) * S_K2
    v[m] = d * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_v1(rows, sigs, n, K1, K2, d_secret, centered=False):
    """The d-column still carries d verbatim, centered or not."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
    for row in rows:
        if abs(row[dim - 1]) != S_KAN: continue
        sgn = 1 if row[dim - 1] > 0 else -1
        d_cand = (sgn * row[m]) % n
        if d_cand and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# V2 / V2c -- d-column dropped
# ---------------------------------------------------------------------------
# columns:  0..m-1 = k1 residuals | m..2m-1 = k2 | 2m = Kannan   (2m+1 cols)
# rows:     m x (n*S_K1*e_i) | row_d | m x row_k2 | row_Kannan    (2m+2 rows)

def build_v2(sigs, n, lam, K1, K2, centered=False):
    m = len(sigs)
    ncol = 2 * m + 1
    S_K1, _S_D, S_K2, S_KAN = scales(n, K1, K2)
    rows = []
    for i in range(m):
        r = [0] * ncol; r[i] = n * S_K1; rows.append(r)
    r = [0] * ncol
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)                                   # row_d  (no d coordinate)
    for i in range(m):
        r = [0] * ncol; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    r = [0] * ncol
    c1 = (K1 // 2) * S_K1 if centered else 0
    c2 = (K2 // 2) * S_K2 if centered else 0
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1 - c1
        r[m + i] = -c2
    r[2 * m] = S_KAN
    rows.append(r)                                   # Kannan row
    return rows

def planted_v2(sigs, n, K1, K2, centered=False):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, K1, K2)
    c1 = (K1 // 2) if centered else 0
    c2 = (K2 // 2) if centered else 0
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = (sigs[i]['k1'] - c1) * S_K1
        v[m + i] = (sigs[i]['k2'] - c2) * S_K2
    v[2 * m] = S_KAN
    return v

def _d_from_residuals(k1s, k2s, sigs, n, lam):
    """d = B_i^{-1}(k1_i + lam*k2_i - A_i) mod n, majority-consistent over i."""
    cands = {}
    for i, s in enumerate(sigs):
        try:
            d = (k1s[i] + lam * k2s[i] - s['A']) * modinv(s['B'], n) % n
        except ValueError:
            continue
        cands[d] = cands.get(d, 0) + 1
    if not cands:
        return None
    return max(cands.items(), key=lambda kv: kv[1])[0]

def recover_v2(rows, sigs, n, lam, K1, K2, d_secret, centered=False):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, K1, K2)
    c1 = (K1 // 2) if centered else 0
    c2 = (K2 // 2) if centered else 0
    for row in rows:
        if abs(row[2 * m]) != S_KAN: continue
        sgn = 1 if row[2 * m] > 0 else -1
        k1s, k2s, ok = [], [], True
        for i in range(m):
            a, b = sgn * row[i], sgn * row[m + i]
            if a % S_K1 or b % S_K2:
                ok = False; break
            k1s.append(a // S_K1 + c1)
            k2s.append(b // S_K2 + c2)
        if not ok: continue
        d = _d_from_residuals(k1s, k2s, sigs, n, lam)
        if d and d == d_secret:
            return d
    return None

# ---------------------------------------------------------------------------
# V3 -- no Kannan coordinate; Babai nearest-plane CVP in Z^(2m)
# ---------------------------------------------------------------------------

def build_v3(sigs, n, lam, K1, K2):
    """Generators of L3 subset Z^(2m): columns 0..m-1 = k1, m..2m-1 = k2."""
    m = len(sigs)
    ncol = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, K1, K2)
    rows = []
    for i in range(m):
        r = [0] * ncol; r[i] = n * S_K1; rows.append(r)
    r = [0] * ncol
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * ncol; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    return rows

def target_v3(sigs, n, K1, K2):
    """t with  v_planted - t  centered:  v_planted = ((k1_i-A_i)S_K1, k2_i S_K2)."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, K1, K2)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = (K1 // 2) * S_K1 - sigs[i]['A'] * S_K1
        t[m + i] = (K2 // 2) * S_K2
    return t

def gram_schmidt(B):
    """Exact rational GS of a list of integer row vectors."""
    Bs, mu = [], []
    for i, b in enumerate(B):
        v = [Fraction(x) for x in b]
        row = []
        for j, bs in enumerate(Bs):
            den = sum(x * x for x in bs)
            c = Fraction(0) if den == 0 else sum(Fraction(x) * y for x, y in zip(b, bs)) / den
            row.append(c)
            if c:
                v = [vi - c * bj for vi, bj in zip(v, bs)]
        Bs.append(v); mu.append(row)
    return Bs

def babai_nearest_plane(B, t):
    """Babai nearest-plane on an LLL-reduced integer basis B (list of rows).
    Returns the lattice vector closest to t found by the procedure."""
    Bs = gram_schmidt(B)
    w = [Fraction(x) for x in t]
    coeffs = [0] * len(B)
    for i in range(len(B) - 1, -1, -1):
        den = sum(x * x for x in Bs[i])
        if den == 0:
            continue
        c = sum(x * y for x, y in zip(w, Bs[i])) / den
        ci = int(math.floor(c + Fraction(1, 2)))
        coeffs[i] = ci
        if ci:
            w = [wi - ci * bj for wi, bj in zip(w, B[i])]
    v = [0] * len(t)
    for ci, b in zip(coeffs, B):
        if ci:
            v = [vi + ci * bj for vi, bj in zip(v, b)]
    return v

def lll_rows(rows, strip_zero=True, bkz_beta=0):
    A = IntegerMatrix.from_matrix(rows)
    if bkz_beta:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    out = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    if strip_zero:
        out = [r for r in out if any(r)]
    return out

def run_v3(sigs, n, lam, K1, K2, d_secret):
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, K1, K2)
    B = lll_rows(build_v3(sigs, n, lam, K1, K2))
    if len(B) != 2 * m:
        return None, None
    t = target_v3(sigs, n, K1, K2)
    v = babai_nearest_plane(B, t)
    dist = norm([a - b for a, b in zip(v, t)])
    k1s, k2s = [], []
    for i in range(m):
        if v[i] % S_K1 or v[m + i] % S_K2:
            return None, dist
        k1s.append(v[i] // S_K1 + sigs[i]['A'])
        k2s.append(v[m + i] // S_K2)
    d = _d_from_residuals(k1s, k2s, sigs, n, lam)
    return (d if d == d_secret else None), dist

# ---------------------------------------------------------------------------
# Per-variant single experiment
# ---------------------------------------------------------------------------

VARIANTS = ["V1", "V1c", "V2", "V2c", "V3"]

def experiment(curve, m, d_secret, K1, seed, variant, bkz_beta=0):
    """Returns (recovered: bool, sv_over_pv: float|None)."""
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, seed)
    if len(sigs) < m:
        return False, None

    if variant in ("V1", "V1c"):
        cen = (variant == "V1c")
        rows = lll_rows(build_v1(sigs, n, lam, K1, K2, centered=cen), bkz_beta=bkz_beta)
        pv = norm(planted_v1(sigs, d_secret, n, K1, K2, centered=cen))
        sv = min(norm(r) for r in rows)
        ok = recover_v1(rows, sigs, n, K1, K2, d_secret, centered=cen) is not None
        return ok, sv / pv

    if variant in ("V2", "V2c"):
        cen = (variant == "V2c")
        rows = lll_rows(build_v2(sigs, n, lam, K1, K2, centered=cen), bkz_beta=bkz_beta)
        pv = norm(planted_v2(sigs, n, K1, K2, centered=cen))
        sv = min(norm(r) for r in rows)
        ok = recover_v2(rows, sigs, n, lam, K1, K2, d_secret, centered=cen) is not None
        return ok, sv / pv

    if variant == "V3":
        d, dist = run_v3(sigs, n, lam, K1, K2, d_secret)
        # "sv/pv" analogue for a CVP formulation: Babai distance / true error norm
        m_ = len(sigs)
        S_K1, _S, S_K2, _K = scales(n, K1, K2)
        err = [0] * (2 * m_)
        for i in range(m_):
            err[i] = (sigs[i]['k1'] - K1 // 2) * S_K1
            err[m_ + i] = (sigs[i]['k2'] - K2 // 2) * S_K2
        pv = norm(err)
        return d is not None, (dist / pv if dist and pv else None)

    raise ValueError(variant)

def rate(curve, m, K1, seeds, variant, bkz_beta=0):
    n = curve[2]
    wins, ratios = 0, []
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        ok, r = experiment(curve, m, d_trial, K1, seed, variant, bkz_beta=bkz_beta)
        wins += ok
        if r is not None:
            ratios.append(r)
    return wins, (sum(ratios) / len(ratios) if ratios else float('nan'))

# ---------------------------------------------------------------------------
# Curves -- identical to the 2026-07-29 T4/T5 table
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]
HIST = [
    # label,           p,    b, n,    lam,  K1_hist, m
    ("8-bit/199",      211,  2, 199,  106,  2,  6),
    ("12-bit/2557",    2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",    2677, 2, 2647, 185,  8,  10),
]
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]

def load_curves():
    out = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        out.append((label, (p, b, n, lam, G), k1, m))
    return out

def hdr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)

# ---------------------------------------------------------------------------
# E1 -- does the trivial vector actually disappear?
# ---------------------------------------------------------------------------

def E1(curves):
    hdr("E1  Shortest-vector structure: is the planted vector lambda_1 now?")
    print("sv/pv  = ||shortest LLL row|| / ||v_planted||   (V3: Babai dist / ||err||)")
    print("triv   = 1 if the shortest row is exactly the trivial n*S_D*e_d vector\n")
    print(f"{'curve':<14}{'K1':>4}{'m':>4}  " + "".join(f"{v+' sv/pv':>12}" for v in VARIANTS)
          + f"{'V1 triv':>9}")
    for label, curve, K1, m in curves:
        p, b, n, lam, G = curve
        d_secret = random.Random(7777 + 42).randint(1, n - 1)
        K2 = math.isqrt(n) + 1
        sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, 42)
        S_K1, S_D, S_K2, S_KAN = scales(n, K1, K2)
        ratios = []
        for v in VARIANTS:
            _ok, r = experiment(curve, m, d_secret, K1, 42, v)
            ratios.append(r)
        rows1 = lll_rows(build_v1(sigs, n, lam, K1, K2))
        sv1 = min(rows1, key=norm)
        triv = int(abs(sv1[m]) == n * S_D and all(x == 0 for j, x in enumerate(sv1) if j != m))
        print(f"{label:<14}{K1:>4}{m:>4}  "
              + "".join(f"{(f'{r:.3f}' if r is not None and r == r else 'n/a'):>12}" for r in ratios)
              + f"{triv:>9}")

# ---------------------------------------------------------------------------
# E2 -- rank check: the deficiency of the V2 generating set is exactly 1
# ---------------------------------------------------------------------------

def E2(curves):
    hdr("E2  Rank of the V2 generating set (deficiency should be exactly 1)")
    print(f"{'curve':<14}{'m':>4}{'gens':>6}{'cols':>6}{'rank(LLL)':>11}{'defic':>7}")
    for label, curve, K1, m in curves:
        p, b, n, lam, G = curve
        d_secret = random.Random(7777 + 42).randint(1, n - 1)
        K2 = math.isqrt(n) + 1
        sigs = gen_signatures(G, d_secret, m, n, lam, p, K1, K2, 42)
        rows = build_v2(sigs, n, lam, K1, K2)
        red = lll_rows(rows)
        print(f"{label:<14}{m:>4}{len(rows):>6}{len(rows[0]):>6}{len(red):>11}"
              f"{len(rows) - len(red):>7}")

# ---------------------------------------------------------------------------
# E3 -- the T4 grid, re-run for every variant.  Does the K1 wall move?
# ---------------------------------------------------------------------------

def E3(curves):
    hdr("E3  T4 grid re-run: successes out of 5 seeds, per K1, per variant")
    for label, curve, _K1, m in curves:
        p, b, n, lam, G = curve
        print(f"\n{label}  p={p} n={n} lam={lam} lam*={lam_star(lam, n):.4f} m={m}")
        print(f"  {'variant':<9}" + "".join(f"{('K1='+str(k)):>7}" for k in K1_GRID)
              + f"{'wall':>7}")
        for v in VARIANTS:
            cells, wall = [], None
            for K1 in K1_GRID:
                w, _ = rate(curve, m, K1, SEEDS, v)
                cells.append(w)
                if w == len(SEEDS):
                    wall = K1
            print(f"  {v:<9}" + "".join(f"{c:>7}" for c in cells)
                  + f"{(str(wall) if wall else '-'):>7}")
        sys.stdout.flush()

# ---------------------------------------------------------------------------
# E4 -- m-sweep at the historical failure point (K1=8 on 12-bit/2677)
# ---------------------------------------------------------------------------

def E4(curves):
    hdr("E4  m-sweep at K1=8 (the 2026-07-29 T4b point where more data did not help)")
    M_GRID = [8, 12, 16, 24, 32]
    for label, curve, _K1, _m in curves:
        if not label.startswith("12-bit"):
            continue
        print(f"\n{label}")
        print(f"  {'variant':<9}" + "".join(f"{('m='+str(mm)):>7}" for mm in M_GRID))
        for v in VARIANTS:
            cells = []
            for mm in M_GRID:
                w, _ = rate(curve, mm, 8, SEEDS, v)
                cells.append(w)
            print(f"  {v:<9}" + "".join(f"{c:>7}" for c in cells))
        sys.stdout.flush()

# ---------------------------------------------------------------------------
# E5 -- replicate the 2026-07-29 T3 sweep (20 fresh 17-bit curves) per variant
# ---------------------------------------------------------------------------
# Reference (T3, variant V1):  eff=0.05 -> 19/20 curves at 5/5,
#                              eff=0.15 ->  3/20,  eff=0.25 -> 0/20.

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
    """Verbatim from glv_hnp_phase2_lambda_threshold.py:373 -- j=0 GLV curves
    with p in [lo,hi), n prime, n = 1 mod 3, bucketed by lam*."""
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

def E5(m_sigs=12):
    hdr("E5  T3 sweep replication on 20 fresh 17-bit curves, V1 vs V2c vs V3")
    print("2026-07-29 reference (V1): eff=0.05 -> 19/20 at 5/5; 0.15 -> 3/20; 0.25 -> 0/20\n")
    print(f"Searching j=0 GLV curves with p in [2^16, 2^17) ...")
    curves = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
    print(f"Found {len(curves)} curves.  m = {m_sigs}, seeds = {SEEDS}")
    for eff_t in (0.05, 0.15, 0.25, 0.35):
        print(f"\n--- eff target = {eff_t} ---")
        print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6}"
              + "".join(f"{v:>8}" for v in ("V1", "V1c", "V2c", "V3")))
        tally = {v: 0 for v in ("V1", "V1c", "V2c", "V3")}
        for (p, b, n, lam, G) in curves:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            eff = k1 * k2 / n
            curve = (p, b, n, lam, G)
            cells = []
            for v in ("V1", "V1c", "V2c", "V3"):
                w, _ = rate(curve, m_sigs, k1, SEEDS, v)
                cells.append(w)
                if w == len(SEEDS):
                    tally[v] += 1
            print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1:>4} {eff:>6.3f}"
                  + "".join(f"{c:>8}" for c in cells))
            sys.stdout.flush()
        print(f"{'CURVES AT 5/5':>32}"
              + "".join(f"{str(tally[v]) + '/' + str(len(curves)):>8}"
                        for v in ("V1", "V1c", "V2c", "V3")))

# ---------------------------------------------------------------------------

def main():
    curves = load_curves()
    print("Thread 23 -- Phase-2 lattice reformulation (d-column removal)")
    print(f"seeds = {SEEDS}")
    E1(curves)
    E2(curves)
    E3(curves)
    E4(curves)
    E5()
    E6()
    print("\nDone.")

if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------------
# E6 -- where is the new ceiling?  BKZ-20 on top of the centered lattice.
# ---------------------------------------------------------------------------

def E6(m_sigs=12):
    hdr("E6  New ceiling: V1c with LLL vs BKZ-20, on the same 20 curves")
    curves = search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=10)
    print(f"{len(curves)} curves, m = {m_sigs}, seeds = {SEEDS}\n")
    print(f"{'eff':>6}{'V1 LLL':>10}{'V1c LLL':>10}{'V1c BKZ20':>12}")
    for eff_t in (0.25, 0.35, 0.45):
        tal = {'a': 0, 'b': 0, 'c': 0}
        for (p, b, n, lam, G) in curves:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            curve = (p, b, n, lam, G)
            if rate(curve, m_sigs, k1, SEEDS, "V1")[0] == len(SEEDS):   tal['a'] += 1
            if rate(curve, m_sigs, k1, SEEDS, "V1c")[0] == len(SEEDS):  tal['b'] += 1
            if rate(curve, m_sigs, k1, SEEDS, "V1c", bkz_beta=20)[0] == len(SEEDS):
                tal['c'] += 1
        N = len(curves)
        print(f"{eff_t:>6.2f}{f'{tal[chr(97)]}/{N}':>10}{f'{tal[chr(98)]}/{N}':>10}"
              f"{f'{tal[chr(99)]}/{N}':>12}")
        sys.stdout.flush()
