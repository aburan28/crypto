"""
GLV-HNP Thread 23: remove the parasitic short vector from the Phase-2 lattice.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, finding T5)
----------------------------------------------------------------
The Phase-2 lattice `build_glv_lattice()` in
`glv_hnp_phase2_lambda_threshold.py:254` is (2m+2)-dimensional with columns

    [ k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan ]

and scales S_K1 = n//K1, S_D = 1, S_K2 = n//K2, S_KANNAN = n.

T5 showed that on EVERY curve tested the vector found by LLL is exactly the
trivial lattice vector  n*S_D*e_m  ("the d-column vector"): 100% of its energy
sits in the d coordinate, |sv[m]|/n = 1.0000 exactly.  It is in the lattice
because d is only determined mod n (n*row_m - sum_i B_i*row_i = n*S_D*e_m).
Its norm is n*S_D, while

    ||v_planted|| ~ n*sqrt(2m/3 + 4/3)   >   n*S_D   for every m >= 1,

so the planted vector is NEVER lambda_1 and no choice of S_D fixes this (both
vectors scale linearly in S_D).  Recovery is therefore a BDD/coset condition,
not an SVP condition.

The Thread 23 proposal was: reformulate so that the target IS lambda_1.

Key observation this script implements
--------------------------------------
The d-column is REDUNDANT.  Once the small unknowns k1_i, k2_i are known,

    d = (k1_i + lam*k2_i - A_i) * B_i^{-1}   (mod n)

recovers d by back-substitution.  So we may DELETE column m entirely.  Deleting
it is exactly quotienting L by the rank-1 sublattice Z*(n*S_D*e_m), i.e. the
parasitic vector is not merely deprioritised, it is annihilated (it maps to 0),
and the dimension DROPS by one.

Three lattices are compared on identical signature sets and seeds:

  V0  baseline           dim 2m+2, Kannan embedding, d read off column m.
                         (verbatim `build_glv_lattice`, reproduces the
                         2026-07-29 T4 grid as a control.)

  V1  d-projected Kannan dim 2m+1.  Column m deleted.  Basis is written in
                         Hermite form so the row set is a true basis:
                             r_0    = S_K1 * (1, C_1, ..., C_{m-1} | 0 | 0)
                             r_i    = n*S_K1*e_i              (i = 1..m-1)
                             q_i    = -lam*S_K1*e_i + S_K2*e_{m+i}
                             r_kan  = (A_i*S_K1 | 0 | S_KANNAN)
                         with C_i = B_i * B_0^{-1} mod n.
                         det(V1) = det(V0)/n, dim(V1) = dim(V0)-1.

  V2  d-projected CVP    dim 2m.  As V1 but the Kannan row is dropped too and
                         the problem is solved as an explicit CVP by Babai
                         nearest-plane against the target
                             t = (-A_i*S_K1 | 0).
                         The offset  v - t = (k1_i*S_K1 | k2_i*S_K2)  has no
                         d term and no Kannan term at all, so the target norm
                         is the smallest of the three.

Falsifier stated in the 2026-07-29 log
--------------------------------------
  "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a
   real improvement; if the wall stays at K1 ~ 4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

E1 answers the first half analytically (Gaussian-heuristic ratios), E2 answers
the second half empirically on the exact T4 grid, E3 measures sv/pv directly.

Run:  python3 glv_hnp_thread23_projected_lattice.py
Deps: fpylll (+ cysignals), sympy.
"""

import math
import os
import random
from fractions import Fraction

from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Reuse the Thread 20 primitives verbatim so the comparison to the 2026-07-29
# numbers is exact.  That module runs its experiments at import time, so only
# the definition prefix (everything before the "EXPERIMENT T1" banner) is
# exec'd here.
# ---------------------------------------------------------------------------
_T20_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "glv_hnp_phase2_lambda_threshold.py")
with open(_T20_PATH) as _fh:
    _prefix = _fh.read().split("# EXPERIMENT T1")[0]
_T20 = {}
exec(compile(_prefix, _T20_PATH, "exec"), _T20)

modinv = _T20["modinv"]
find_generator = _T20["find_generator"]
lam_star = _T20["lam_star"]
scales = _T20["scales"]
gen_signatures = _T20["gen_signatures"]
build_glv_lattice = _T20["build_glv_lattice"]
planted_vector = _T20["planted_vector"]
norm = _T20["norm"]
recover_d = _T20["recover_d"]

TWO_PI_E = 2.0 * math.pi * math.e


# ---------------------------------------------------------------------------
# V1 / V2 lattice construction
# ---------------------------------------------------------------------------

def build_projected_basis(sigs, n, lam, k1_bound, k2_bound, kannan=True):
    """d-projected lattice.  dim = 2m+1 (kannan=True) or 2m (kannan=False).

    Columns: [k1_0..k1_{m-1} | k2_0..k2_{m-1} ( | kannan )].
    """
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + (1 if kannan else 0)
    B0inv = modinv(sigs[0]['B'] % n, n)
    C = [sigs[i]['B'] * B0inv % n for i in range(m)]     # C[0] == 1

    rows = []
    # r_0: Hermite pivot row of the d-image sublattice in the k1 block.
    r0 = [0] * dim
    for i in range(m):
        r0[i] = C[i] * S_K1
    rows.append(r0)
    # r_i = n*S_K1*e_i for i = 1..m-1  (n*e_0 is generated as n*r_0 - sum C_i r_i)
    for i in range(1, m):
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    # k2 rows
    for i in range(m):
        r = [0] * dim
        r[i] = -(lam % n) * S_K1
        r[m + i] = S_K2
        rows.append(r)
    if kannan:
        r = [0] * dim
        for i in range(m):
            r[i] = sigs[i]['A'] * S_K1
        r[dim - 1] = S_KANNAN
        rows.append(r)
    return rows


def projected_planted(sigs, n, k1_bound, k2_bound, kannan=True):
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    dim = 2 * m + (1 if kannan else 0)
    v = [0] * dim
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    if kannan:
        v[dim - 1] = S_KANNAN
    return v


def cvp_target(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _S_D, _S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return t


def d_from_smalls(k1s, k2s, sigs, n, lam):
    """Back-substitute d from every index; return the set of candidates."""
    out = set()
    for i, sg in enumerate(sigs):
        try:
            Binv = modinv(sg['B'] % n, n)
        except Exception:
            continue
        out.add((k1s[i] + lam * k2s[i] - sg['A']) * Binv % n)
    return out


def recover_from_projected_row(row, sigs, n, k1_bound, k2_bound, lam,
                               d_secret, kannan=True):
    """Read (k1,k2) out of a candidate lattice vector and test d."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    if kannan:
        last = row[2 * m]
        if abs(last) != S_KANNAN:
            return False
        sign = 1 if last > 0 else -1
        row = [sign * x for x in row]
    k1s, k2s = [], []
    for i in range(m):
        if row[i] % S_K1 or row[m + i] % S_K2:
            return False
        k1s.append(row[i] // S_K1)
        k2s.append(row[m + i] // S_K2)
    return d_secret in d_from_smalls(k1s, k2s, sigs, n, lam)


# ---------------------------------------------------------------------------
# Exact Babai nearest-plane (Fraction GSO -- dim <= ~30, entries ~1e6)
# ---------------------------------------------------------------------------

def babai_nearest_plane(basis, target):
    """Exact nearest-plane against `basis` (assumed already LLL-reduced)."""
    k = len(basis)
    B = [[Fraction(x) for x in row] for row in basis]
    Bs, mu, nrm = [], [], []
    for i in range(k):
        v = list(B[i])
        mu_i = []
        for j in range(len(Bs)):
            c = sum(B[i][t] * Bs[j][t] for t in range(len(v))) / nrm[j]
            mu_i.append(c)
            for t in range(len(v)):
                v[t] -= c * Bs[j][t]
        Bs.append(v)
        mu.append(mu_i)
        nrm.append(sum(x * x for x in v))

    w = [Fraction(x) for x in target]
    coeffs = [0] * k
    for i in range(k - 1, -1, -1):
        if nrm[i] == 0:
            continue
        c = sum(w[t] * Bs[i][t] for t in range(len(w))) / nrm[i]
        ci = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        coeffs[i] = ci
        for t in range(len(w)):
            w[t] -= ci * B[i][t]
    v = [0] * len(target)
    for i in range(k):
        if not coeffs[i]:
            continue
        for t in range(len(v)):
            v[t] += coeffs[i] * basis[i][t]
    return v


# ---------------------------------------------------------------------------
# One trial of each variant
# ---------------------------------------------------------------------------

def reduce_rows(rows, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]


def trial(curve, m, d_secret, k1_bound, seed, variant, use_bkz=False, beta=20):
    """variant in {'V0','V1','V2'}.  Returns dict or None."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    _S_K1, _S_D, _S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)

    if variant == 'V0':
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        red = reduce_rows(M, use_bkz, beta)
        ok = recover_d(red, m, n, S_KANNAN, d_secret) is not None
        pv = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
        sv = min(norm(r) for r in red if any(r))
        return {'ok': ok, 'pv': pv, 'sv': sv, 'dim': 2 * m + 2}

    if variant == 'V1':
        M = build_projected_basis(sigs, n, lam, k1_bound, k2_bound, kannan=True)
        red = reduce_rows(M, use_bkz, beta)
        ok = any(recover_from_projected_row(r, sigs, n, k1_bound, k2_bound,
                                            lam, d_secret, kannan=True)
                 for r in red)
        pv = norm(projected_planted(sigs, n, k1_bound, k2_bound, kannan=True))
        sv = min(norm(r) for r in red if any(r))
        return {'ok': ok, 'pv': pv, 'sv': sv, 'dim': 2 * m + 1}

    if variant == 'V2':
        M = build_projected_basis(sigs, n, lam, k1_bound, k2_bound, kannan=False)
        red = reduce_rows(M, use_bkz, beta)
        t = cvp_target(sigs, n, k1_bound, k2_bound)
        v = babai_nearest_plane(red, t)
        off = [v[i] - t[i] for i in range(len(t))]
        ok = recover_from_projected_row(off, sigs, n, k1_bound, k2_bound,
                                        lam, d_secret, kannan=False)
        pv = norm(projected_planted(sigs, n, k1_bound, k2_bound, kannan=False))
        sv = norm(off)          # for V2 this is ||Babai error||, not lambda_1
        return {'ok': ok, 'pv': pv, 'sv': sv, 'dim': 2 * m}

    raise ValueError(variant)


def rate(curve, m, k1_bound, seeds, variant, use_bkz=False, beta=20):
    wins, tot, ratios = 0, 0, []
    for s in seeds:
        d_trial = random.Random(s + 7777).randint(1, curve[2] - 1)
        r = trial(curve, m, d_trial, k1_bound, s, variant, use_bkz, beta)
        if r is None:
            continue
        tot += 1
        wins += bool(r['ok'])
        if r['pv']:
            ratios.append(r['sv'] / r['pv'])
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, tot, mr


# ---------------------------------------------------------------------------
# Analytic Gaussian-heuristic predictions
# ---------------------------------------------------------------------------

def gh_table(n, m, k1_bound):
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    ln = math.log
    out = {}

    # V0
    dim = 2 * m + 2
    logdet = m * ln(n * S_K1) + ln(S_D) + m * ln(S_K2) + ln(S_KAN)
    tgt = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3 + (n * S_D) ** 2 / 3
                    + m * (k2_bound * S_K2) ** 2 / 3 + S_KAN ** 2)
    out['V0'] = (dim, logdet, tgt)

    # V1  (det(V0)/n, one fewer dimension, no d term in the target)
    dim = 2 * m + 1
    logdet = (m - 1) * ln(n) + m * ln(S_K1) + m * ln(S_K2) + ln(S_KAN)
    tgt = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3
                    + m * (k2_bound * S_K2) ** 2 / 3 + S_KAN ** 2)
    out['V1'] = (dim, logdet, tgt)

    # V2  (no Kannan coordinate; target is a BDD distance)
    dim = 2 * m
    logdet = (m - 1) * ln(n) + m * ln(S_K1) + m * ln(S_K2)
    tgt = math.sqrt(m * (k1_bound * S_K1) ** 2 / 3
                    + m * (k2_bound * S_K2) ** 2 / 3)
    out['V2'] = (dim, logdet, tgt)

    res = {}
    for k, (dim, logdet, tgt) in out.items():
        gh = math.sqrt(dim / TWO_PI_E) * math.exp(logdet / dim)
        res[k] = {'dim': dim, 'gh': gh, 'tgt': tgt, 'ratio': tgt / gh}
    return res


# ===========================================================================
if __name__ == "__main__":
    print("=" * 78)
    print("Thread 23 - does deleting the d-column make the planted vector lambda_1?")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337]
    HIST = [
        # label,              p,    b, n,    lam
        ("12-bit/2557 lam*=0.340", 2557, 2, 2659, 1755),
        ("12-bit/2677 lam*=0.070", 2677, 2, 2647, 185),
    ]
    curves = []
    for label, p, b, n, lam in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        assert (lam * lam + lam + 1) % n == 0, label
        curves.append((label, (p, b, n, lam, G)))

    M_SIGS = 12                      # matches the 2026-07-29 T4 grid
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("E1  ANALYTIC: Gaussian-heuristic ratio  ||target|| / GH(L)")
    print("-" * 78)
    print("ratio < 1  =>  the target is (heuristically) the unique shortest /")
    print("closest vector and lattice reduction should find it.\n")
    for label, cur in curves:
        n = cur[2]
        print(f"{label}   n={n}  m={M_SIGS}  K2={math.isqrt(n)+1}")
        print(f"  {'K1':>3} {'eff=K1K2/n':>11} | "
              f"{'V0 dim':>6} {'V0 ratio':>9} | {'V1 dim':>6} {'V1 ratio':>9} | "
              f"{'V2 dim':>6} {'V2 ratio':>9}")
        for k1 in K1_GRID:
            g = gh_table(n, M_SIGS, k1)
            eff = k1 * (math.isqrt(n) + 1) / n
            print(f"  {k1:>3} {eff:>11.3f} | "
                  f"{g['V0']['dim']:>6} {g['V0']['ratio']:>9.3f} | "
                  f"{g['V1']['dim']:>6} {g['V1']['ratio']:>9.3f} | "
                  f"{g['V2']['dim']:>6} {g['V2']['ratio']:>9.3f}")
        print()

    # -------------------------------------------------------------------
    print("-" * 78)
    print("E2  EMPIRICAL: the 2026-07-29 T4 K1 grid, three lattice variants")
    print("-" * 78)
    print("V0 row must reproduce the 2026-07-29 T4 table (control).")
    print(f"m={M_SIGS}, {len(SEEDS)} seeds {SEEDS}\n")
    e2 = {}
    for label, cur in curves:
        print(label)
        header = "  " + f"{'variant':<10}" + "".join(f"{('K1=' + str(k)):>8}"
                                                     for k in K1_GRID)
        print(header)
        for variant in ('V0', 'V1', 'V2'):
            cells = []
            for k1 in K1_GRID:
                w, t, _ = rate(cur, M_SIGS, k1, SEEDS, variant)
                e2[(label, variant, k1)] = (w, t)
                cells.append(f"{w}/{t}")
            print("  " + f"{variant:<10}" + "".join(f"{c:>8}" for c in cells))
        # BKZ-20 retry on the projected CVP variant
        cells = []
        for k1 in K1_GRID:
            w, t, _ = rate(cur, M_SIGS, k1, SEEDS, 'V2', use_bkz=True, beta=20)
            e2[(label, 'V2-BKZ20', k1)] = (w, t)
            cells.append(f"{w}/{t}")
        print("  " + f"{'V2-BKZ20':<10}" + "".join(f"{c:>8}" for c in cells))
        print()

    # -------------------------------------------------------------------
    print("-" * 78)
    print("E3  sv/pv DIAGNOSTIC: is the planted vector lambda_1 now?")
    print("-" * 78)
    print("V0/V1: sv = ||shortest reduced row||, pv = ||planted||.")
    print("       sv/pv >= 1 means nothing in the reduced basis beats the")
    print("       planted vector, i.e. the parasitic vector is gone.")
    print("V2:    sv = ||Babai error||; sv/pv == 1.000 means Babai landed")
    print("       exactly on the planted coset representative.\n")
    print(f"  {'curve':<24}{'K1':>4}{'V0 sv/pv':>10}{'V1 sv/pv':>10}"
          f"{'V2 err/pv':>11}")
    for label, cur in curves:
        for k1 in (2, 4, 8, 16):
            row = [label, k1]
            for variant in ('V0', 'V1', 'V2'):
                _, _, mr = rate(cur, M_SIGS, k1, SEEDS, variant)
                row.append(mr)
            print(f"  {row[0]:<24}{row[1]:>4}{row[2]:>10.3f}{row[3]:>10.3f}"
                  f"{row[4]:>11.3f}")
    print()

    # -------------------------------------------------------------------
    print("-" * 78)
    print("E4  WALL LOCATION (largest K1 with >= 4/5 recoveries)")
    print("-" * 78)
    for label, _cur in curves:
        line = [f"  {label:<24}"]
        for variant in ('V0', 'V1', 'V2', 'V2-BKZ20'):
            wall = 0
            for k1 in K1_GRID:
                w, t = e2[(label, variant, k1)]
                if t and w >= 4:
                    wall = k1
            line.append(f"{variant}={wall}")
        print("  ".join(line))
    print()

    # -------------------------------------------------------------------
    # E5 asks what the wall IS, given that E2/E4 showed it is not a
    # consequence of the parasitic d-vector.  Two candidate explanations:
    #   (a) information-theoretic: the m signatures do not determine d.
    #       Bound: (m-1)*log n >= m*log(K1*K2), i.e. K1 <= n^{(m-1)/m}/K2.
    #   (b) lattice-reduction: the wall sits where the Gaussian-heuristic
    #       ratio ||target||/GH(L) crosses 1.
    # These make very different predictions, so the wall location settles it.
    # -------------------------------------------------------------------
    print("-" * 78)
    print("E5  WHAT IS THE WALL?  info-theoretic bound vs Gaussian heuristic")
    print("-" * 78)
    print(f"  {'curve':<24}{'K1 obs':>8}{'K1 info':>9}{'K1 GH':>8}")
    for label, cur in curves:
        n = cur[2]
        k2 = math.isqrt(n) + 1
        obs = 0
        for k1 in K1_GRID:
            w, t = e2[(label, 'V0', k1)]
            if t and w >= 4:
                obs = k1
        k1_info = n ** ((M_SIGS - 1) / M_SIGS) / k2
        k1_gh = 0
        for k1 in range(2, 200):
            if gh_table(n, M_SIGS, k1)['V0']['ratio'] <= 1.0:
                k1_gh = k1
        print(f"  {label:<24}{obs:>8}{k1_info:>9.1f}{k1_gh:>8}")
    print()

    # -------------------------------------------------------------------
    print("-" * 78)
    print("E6  IS THE lam* WALL SHIFT REAL?  fresh-curve sweep, V0")
    print("-" * 78)
    print("E1 shows GH(L) is essentially identical for the two historical")
    print("curves (ratios agree to 3 decimals), yet E4 puts their walls at")
    print("K1=12 vs K1=4.  If lam* drives that, the wall should track lam*")
    print("across fresh curves; if it is per-curve noise, it should not.\n")
    def pearson(xs, ys):
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        vx = math.sqrt(sum((x - mx) ** 2 for x in xs))
        vy = math.sqrt(sum((y - my) ** 2 for y in ys))
        return cov / (vx * vy) if vx and vy else float('nan')

    search_curves = _T20["search_curves"]
    sweep = search_curves(2 ** 12, 2 ** 13, per_bin=2, nbins=8)
    print(f"  found {len(sweep)} fresh 12-13 bit j=0 GLV curves")
    print("  'wall'  = largest K1 in the grid with >= 4/5 recoveries")
    print("  'score' = total recoveries over the whole grid (0..40); a smoother")
    print("            statistic than the wall, which is a noisy argmax.\n")
    print(f"  {'p':>7}{'n':>7}{'lam*':>8}{'wall':>6}{'score':>7}{'GH K1':>7}   "
          + "".join(f"{('K'+str(k)):>5}" for k in K1_GRID))
    rows = []
    for (p, b, n, lam, G) in sweep:
        cur = (p, b, n, lam, G)
        cells, wall, score = [], 0, 0
        for k1 in K1_GRID:
            w, t, _ = rate(cur, M_SIGS, k1, SEEDS, 'V0')
            cells.append(f"{w}/{t}")
            score += w
            if t and w >= 4:
                wall = k1
        k1_gh = 0
        for k1 in range(2, 200):
            if gh_table(n, M_SIGS, k1)['V0']['ratio'] <= 1.0:
                k1_gh = k1
        ls = lam_star(lam, n)
        rows.append((ls, wall, score, p, n))
        print(f"  {p:>7}{n:>7}{ls:>8.3f}{wall:>6}{score:>7}{k1_gh:>7}   "
              + "".join(f"{c:>5}" for c in cells))
    if len(rows) >= 3:
        xs = [r[0] for r in rows]
        rw = pearson(xs, [math.log2(r[1]) if r[1] else 0.0 for r in rows])
        rs = pearson(xs, [r[2] for r in rows])
        print(f"\n  Pearson r(lam*, log2 wall)  = {rw:+.3f}   over {len(rows)} curves")
        print(f"  Pearson r(lam*, score)      = {rs:+.3f}")
        lo = [r_ for r_ in rows if r_[0] < 0.25]
        hi = [r_ for r_ in rows if r_[0] >= 0.25]
        if lo and hi:
            print(f"  mean wall   lam*<0.25: {sum(r_[1] for r_ in lo)/len(lo):.2f}"
                  f" (n={len(lo)})   lam*>=0.25:"
                  f" {sum(r_[1] for r_ in hi)/len(hi):.2f} (n={len(hi)})")
            print(f"  mean score  lam*<0.25: {sum(r_[2] for r_ in lo)/len(lo):.2f}"
                  f"              lam*>=0.25:"
                  f" {sum(r_[2] for r_ in hi)/len(hi):.2f}")
        print(f"  wall range over the sample: {min(r_[1] for r_ in rows)}"
              f" .. {max(r_[1] for r_ in rows)}  (GH predicts the same K1 for all)")
    print()

    # -------------------------------------------------------------------
    # Control.  E6's wall is an argmax over 5 seeds, so a 6-vs-16 spread
    # could be sampling noise rather than a curve property.  Re-measure the
    # extreme curves on three DISJOINT seed blocks: if the wall is a curve
    # property the three blocks agree, if it is noise they scatter.
    # -------------------------------------------------------------------
    print("-" * 78)
    print("E7  SEED-STABILITY CONTROL on the extreme curves of E6")
    print("-" * 78)
    BLOCKS = [[42, 1234, 9999, 555, 31337],
              [7, 71, 701, 7001, 70001],
              [2024, 2025, 2026, 2027, 2028]]
    extremes = sorted(rows, key=lambda r_: r_[1])
    picks = extremes[:2] + extremes[-2:]
    print(f"  {'p':>7}{'n':>7}{'lam*':>8}{'blk1':>6}{'blk2':>6}{'blk3':>6}"
          f"{'wall(15 seeds)':>16}")
    for (ls, _w, _sc, p, n) in picks:
        cur = next(c for c in sweep if c[0] == p and c[2] == n)
        cur = (cur[0], cur[1], cur[2], cur[3], cur[4])
        walls = []
        for blk in BLOCKS:
            wl = 0
            for k1 in K1_GRID:
                w, t, _ = rate(cur, M_SIGS, k1, blk, 'V0')
                if t and w >= 4:
                    wl = k1
            walls.append(wl)
        allseeds = BLOCKS[0] + BLOCKS[1] + BLOCKS[2]
        wl_all = 0
        for k1 in K1_GRID:
            w, t, _ = rate(cur, M_SIGS, k1, allseeds, 'V0')
            if t and w >= 0.8 * t:
                wl_all = k1
        print(f"  {p:>7}{n:>7}{ls:>8.3f}{walls[0]:>6}{walls[1]:>6}"
              f"{walls[2]:>6}{wl_all:>16}")
    print()
