"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
the target of a BDD/CVP instance rather than a (hopeless) SVP instance.

Motivation (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 lattice always contains the trivial vector n*S_D*e_m (norm
  n*S_D), because  n*row_m - sum_i B_i*row_i  cancels every k1-column.  The
  planted vector has norm ~ n*sqrt(2m/3 + 4/3) > n, so it is NEVER lambda_1.
  Recovery is therefore not an SVP condition but a BDD/coset condition.
  That run proposed Thread 23: quotient out the trivial direction and/or
  replace the Kannan embedding by an explicit CVP call, then ask whether the
  K1 wall measured in T4 moves outward.

  Falsifier stated on 2026-07-29:
    if the wall moves outward on the lam*=0.07 curve (currently K1 ~ 4-6),
    the reformulation is a real improvement; if the wall stays, the wall is
    information-theoretic and Phase 2 is at its ceiling.

Four formulations are compared on identical signature sets:

  full   the original (2m+2)-dim Kannan lattice          [baseline, = T4]
  proj   the same lattice with the d-column DELETED.  Deleting column m is
         exactly the orthogonal projection along e_m, and since
         L cap R*e_m = n*S_D*Z*e_m (proved below in check_A2), the image is
         a rank-(2m+1) integer lattice in which the trivial vector is 0.
         d is recovered afterwards from signature 0.
  bdd    explicit CVP.  Drop the Kannan row AND the Kannan column AND the
         d-column: L'' = rowspan(rows 0..2m) restricted to the k1- and
         k2-columns, rank 2m, det = S_K1^m * n^(m-1) * S_K2^m.  The planted
         relation becomes  (k1_i*S_K1, k2_i*S_K2) = w + a  with w in L'' and
         a = (A_i*S_K1, 0), i.e. w is the lattice point nearest to -a and the
         error IS the planted small vector.
  bdd_c  same lattice, CENTERED target: shift by (K1/2*S_K1, K2/2*S_K2) so the
         error is symmetric about 0.  E[x^2] drops from X^2/3 to X^2/12, i.e.
         ||error|| halves.  This costs nothing and is the only one of the four
         changes that alters the error norm rather than just the dimension.

Analytic prediction being tested (derived in the log entry for 2026-08-03):
  det(L'') = n^(3m-1)/(K1*K2)^m ,  ||error||_uncentered = n*sqrt(2m/3), so
      rho_GH := ||error|| / GH(L'') = sqrt(2*pi*e/3) * sqrt(eff) * n^(1/(2m))
  with eff = K1*K2/n.  The leading factor sqrt(2*pi*e/3) = 2.386, so
  rho_GH = 1 at eff = 0.176 * n^(-1/m).  Centering divides rho_GH by 2, which
  should move the eff wall by a factor of 4 (K1 wall by 4x) if rho_GH is the
  governing quantity.  Note rho_GH does NOT depend on lambda at all, so if the
  lam*=0.34 and lam*=0.07 curves keep their factor-~3 gap under bdd_c, then
  lambda acts through the GS profile of L'' and not through its determinant.

Usage:  python3 glv_hnp_thread23_bdd.py
"""

import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fpylll import IntegerMatrix, LLL, BKZ, GSO, CVP

import glv_hnp_common as C

TWO_PI_E = 2.0 * math.pi * math.e


# ---------------------------------------------------------------------------
# Lattice surgery
# ---------------------------------------------------------------------------

def rows_of(M):
    return [list(r) for r in M]


def delete_column(M, j):
    return [[v for k, v in enumerate(row) if k != j] for row in M]


def lll_rows(M, drop_zero=True):
    """LLL-reduce a list-of-lists (possibly rank-deficient).  fpylll pushes the
    zero rows to the top; drop them so the result is a basis."""
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    out = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    if drop_zero:
        out = [r for r in out if any(r)]
    return out


def gs_profile(basis):
    """||b*_i|| for an integer basis, via fpylll GSO."""
    A = IntegerMatrix.from_matrix(basis)
    M = GSO.Mat(A)
    M.update_gso()
    return [math.sqrt(M.get_r(i, i)) for i in range(A.nrows)]


def log_det(basis):
    """log(det) of a full-rank integer basis, from the GS profile (stable)."""
    return sum(math.log(x) for x in gs_profile(basis))


def gaussian_heuristic(basis):
    d = len(basis)
    return math.exp(math.log(d / TWO_PI_E) / 2.0 + log_det(basis) / d)


def babai_nearest_plane(basis, target):
    """Classical Babai nearest-plane on an LLL-reduced basis.  Returns the
    lattice point.  Pure python w/ floats -- fine at these sizes."""
    A = IntegerMatrix.from_matrix(basis)
    M = GSO.Mat(A)
    M.update_gso()
    d = A.nrows
    b = [list(r) for r in basis]
    bstar = []
    # explicit GS in floats
    for i in range(d):
        v = [float(x) for x in b[i]]
        for j in range(len(bstar)):
            num = sum(v[k] * bstar[j][k] for k in range(len(v)))
            den = sum(x * x for x in bstar[j])
            if den:
                mu = num / den
                v = [v[k] - mu * bstar[j][k] for k in range(len(v))]
        bstar.append(v)
    w = [float(x) for x in target]
    coeffs = [0] * d
    for i in range(d - 1, -1, -1):
        den = sum(x * x for x in bstar[i])
        if den == 0:
            continue
        num = sum(w[k] * bstar[i][k] for k in range(len(w)))
        c = int(math.floor(num / den + 0.5))
        coeffs[i] = c
        w = [w[k] - c * b[i][k] for k in range(len(w))]
    pt = [0] * len(target)
    for i in range(d):
        if coeffs[i]:
            for k in range(len(pt)):
                pt[k] += coeffs[i] * b[i][k]
    return pt


# ---------------------------------------------------------------------------
# The three reformulated attacks.  All take the SAME signature set.
# ---------------------------------------------------------------------------

def attack_full(sigs, n, lam, K1, K2, d_secret, use_bkz=False, beta=20):
    """Baseline: original Kannan lattice, SVP-style read-off (as in T4)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = C.scales(n, K1, K2)
    M = C.build_glv_lattice(sigs, n, lam, K1, K2)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    red = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    return C.recover_d(red, m, n, S_KAN, d_secret) is not None


def _d_from_k(sigs, k1_0, k2_0, n, lam):
    """d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1} mod n."""
    B0 = sigs[0]['B'] % n
    if B0 == 0:
        return None
    kfull = (k1_0 + lam * k2_0) % n
    return (kfull - sigs[0]['A']) * C.modinv(B0, n) % n


def attack_proj(sigs, n, lam, K1, K2, d_secret, use_bkz=False, beta=20):
    """d-column deleted; read the Kannan row off the reduced basis."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = C.scales(n, K1, K2)
    M = delete_column(C.build_glv_lattice(sigs, n, lam, K1, K2), m)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    dim = 2 * m + 1
    for i in range(A.nrows):
        row = [A[i][j] for j in range(dim)]
        last = row[dim - 1]
        if abs(last) != S_KAN:
            continue
        if last < 0:
            row = [-x for x in row]
        k1_0 = row[0] // S_K1 if row[0] % S_K1 == 0 else None
        k2_0 = row[m] // S_K2 if row[m] % S_K2 == 0 else None
        if k1_0 is None or k2_0 is None:
            continue
        d_cand = _d_from_k(sigs, k1_0, k2_0, n, lam)
        if d_cand is not None and d_cand == d_secret:
            return True
    return False


def build_Lpp(sigs, n, lam, K1, K2):
    """L'' : rows 0..2m of the Phase-2 lattice, restricted to the k1- and
    k2-columns (d-column and Kannan column both removed).  Rank 2m."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = C.scales(n, K1, K2)
    dim = 2 * m
    rows = []
    for i in range(m):                       # n*S_K1*e_i
        r = [0] * dim
        r[i] = n * S_K1
        rows.append(r)
    r = [0] * dim                            # the d-row
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                       # the k2-rows
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    return rows


def attack_bdd(sigs, n, lam, K1, K2, d_secret, centered=False, exact=True):
    """Explicit CVP in L''.  error == the planted small vector."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = C.scales(n, K1, K2)
    dim = 2 * m
    basis = lll_rows(build_Lpp(sigs, n, lam, K1, K2))
    if len(basis) != dim:
        return False, None
    a = [0] * dim
    for i in range(m):
        a[i] = sigs[i]['A'] * S_K1
    # (k1_i*S_K1, k2_i*S_K2) = w + a  with w in L''  =>  w near -a
    t = [-x for x in a]
    if centered:
        for i in range(m):
            t[i] += (K1 // 2) * S_K1
            t[m + i] += (K2 // 2) * S_K2
    if exact:
        try:
            w = list(CVP.closest_vector(IntegerMatrix.from_matrix(basis), tuple(t)))
        except Exception:
            w = babai_nearest_plane(basis, t)
    else:
        w = babai_nearest_plane(basis, t)
    # planted small vector = w + a  (recover k1_0, k2_0 -> d)
    v = [w[k] + a[k] for k in range(dim)]
    if v[0] % S_K1 or v[m] % S_K2:
        return False, None
    d_cand = _d_from_k(sigs, v[0] // S_K1, v[m] // S_K2, n, lam)
    return (d_cand is not None and d_cand == d_secret), v


# ---------------------------------------------------------------------------
# Diagnostics for one (curve, K1, m) cell -- lambda-independent and
# lambda-dependent quantities side by side.
# ---------------------------------------------------------------------------

def diagnostics(sigs, n, lam, K1, K2):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = C.scales(n, K1, K2)
    basis = lll_rows(build_Lpp(sigs, n, lam, K1, K2))
    prof = gs_profile(basis)
    gh = math.exp(math.log(len(basis) / TWO_PI_E) / 2.0
                  + sum(math.log(x) for x in prof) / len(basis))
    # actual planted error, uncentered and centered
    e_u = 0.0
    e_c = 0.0
    for s in sigs:
        e_u += (s['k1'] * S_K1) ** 2 + (s['k2'] * S_K2) ** 2
        e_c += ((s['k1'] - K1 // 2) * S_K1) ** 2 + ((s['k2'] - K2 // 2) * S_K2) ** 2
    e_u, e_c = math.sqrt(e_u), math.sqrt(e_c)
    # shortest vector actually present after LLL (upper bound on lambda_1)
    lam1_ub = min(math.sqrt(sum(x * x for x in r)) for r in basis)
    mu, _ = C.lambda_block_mu(n, lam, S_K1, S_K2)
    return {
        'dim': len(basis), 'gh': gh, 'e_u': e_u, 'e_c': e_c,
        'rho_gh_u': e_u / gh, 'rho_gh_c': e_c / gh,
        'lam1_ub': lam1_ub, 'mu': mu,
        'gs_min': min(prof), 'gs_max': max(prof),
        'np_margin': min(prof) / (2.0 * e_c) if e_c else float('inf'),
    }


# ===========================================================================

def hr(t=''):
    print('\n' + '=' * 78)
    if t:
        print(t)
        print('=' * 78)


SEEDS = [42, 1234, 9999, 555, 31337, 7, 20250803, 65521]

def secret_for(seed, n):
    """d in [1, n-1].  NOTE: earlier drafts used d = 1234 + seed unreduced,
    which silently capped every E1/E3 cell at 3/5 for the 12-bit curves
    (n ~ 2650) because two seeds gave d > n and could never be matched."""
    d = (1234 + seed) % n
    return d if d else 1

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,            p,    b, n,    lam,   K1, m
    ("12-bit/2557", 2557, 2, 2659, 1755, 8, 8),    # lam* = 0.340
    ("12-bit/2677", 2677, 2, 2647, 185, 8, 10),    # lam* = 0.070
]

print(__doc__.split('Usage:')[0])
hr('Thread 23 -- BDD reformulation of the GLV-HNP Phase 2 lattice')

curves = []
for label, p, b, n, lam, k1, m in HIST:
    G = C.find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    curves.append((label, (p, b, n, lam, G), k1, m))
    print(f"  {label:<14} p={p} n={n} lam={lam} lam*={C.lam_star(lam, n):.4f} "
          f"K2={math.isqrt(n)+1}")


# ---------------------------------------------------------------------------
hr('CHECK A -- structural claims about the lattice')
# ---------------------------------------------------------------------------
label, curve, _, m0 = curves[0]
p, b, n, lam, G = curve
K2_0 = math.isqrt(n) + 1
sigs0 = C.gen_signatures(G, 1234, m0, n, lam, p, 8, K2_0, seed=42)
S_K1, S_D, S_K2, S_KAN = C.scales(n, 8, K2_0)

# A1: the trivial vector n*S_D*e_m really is in L
M0 = C.build_glv_lattice(sigs0, n, lam, 8, K2_0)
triv = [0] * (2 * m0 + 2)
for k in range(2 * m0 + 2):
    triv[k] = n * M0[m0][k] - sum(sigs0[i]['B'] * M0[i][k] for i in range(m0))
ok_a1 = (triv[m0] == n * S_D and all(v == 0 for k, v in enumerate(triv) if k != m0))
print(f"A1  n*row_m - sum B_i*row_i == n*S_D*e_m                 : {ok_a1}"
      f"   (norm {abs(triv[m0])} = n*S_D = {n*S_D})")

# A2: L cap R*e_m == n*S_D*Z  (so deleting column m gives a lattice)
#     a*B_i = 0 mod n for all i and n prime, B_i != 0  =>  a = 0 mod n
ok_a2 = all(s['B'] % n != 0 for s in sigs0)
print(f"A2  all B_i != 0 mod n (=> L cap R*e_m = n*S_D*Z)        : {ok_a2}")

# A3: planted vector norm vs trivial vector norm
vp = C.planted_vector(sigs0, 1234, n, 8, K2_0)
print(f"A3  ||v_planted|| / ||trivial|| = {C.norm(vp)/(n*S_D):.3f}"
      f"   (>1 => planted is never lambda_1)")

# A4: det(L'') == S_K1^m * n^(m-1) * S_K2^m   (predicted)
bas = lll_rows(build_Lpp(sigs0, n, lam, 8, K2_0))
pred = m0 * math.log(S_K1) + (m0 - 1) * math.log(n) + m0 * math.log(S_K2)
print(f"A4  log det(L'')  measured {log_det(bas):.6f}   predicted {pred:.6f}"
      f"   rank {len(bas)} (expect {2*m0})")

# A5: the CVP error really is the planted vector
ok5, v5 = attack_bdd(sigs0, n, lam, 8, K2_0, 1234, centered=True)
exp5 = [sigs0[i]['k1'] * S_K1 for i in range(m0)] + \
       [sigs0[i]['k2'] * S_K2 for i in range(m0)]
print(f"A5  CVP error == planted (k1_i*S_K1, k2_i*S_K2)          : {v5 == exp5}"
      f"   (d recovered: {ok5})")


# ---------------------------------------------------------------------------
hr('EXP E1 -- does the reformulation move the K1 wall?  (T4 grid)')
# ---------------------------------------------------------------------------
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
FORMS = ['full', 'proj', 'bdd', 'bdd_c']

e1_rows = []
for label, curve, _k1, m in curves:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    print(f"\n{label}   lam*={C.lam_star(lam,n):.4f}  n={n}  K2={K2}  m={m}")
    print(f"  {'K1':>4} {'eff':>7} | " +
          ' '.join(f"{f:>6}" for f in FORMS) + " |  rho_gh_u rho_gh_c")
    for K1 in K1_GRID:
        eff = K1 * K2 / n
        if eff > 1.2:
            continue
        wins = {f: 0 for f in FORMS}
        diag_acc = None
        for seed in SEEDS:
            d_sec = secret_for(seed, n)
            sigs = C.gen_signatures(G, d_sec, m, n, lam, p, K1, K2, seed=seed)
            if len(sigs) < m:
                continue
            if attack_full(sigs, n, lam, K1, K2, d_sec):
                wins['full'] += 1
            if attack_proj(sigs, n, lam, K1, K2, d_sec):
                wins['proj'] += 1
            if attack_bdd(sigs, n, lam, K1, K2, d_sec, centered=False)[0]:
                wins['bdd'] += 1
            if attack_bdd(sigs, n, lam, K1, K2, d_sec, centered=True)[0]:
                wins['bdd_c'] += 1
            if diag_acc is None:
                diag_acc = diagnostics(sigs, n, lam, K1, K2)
        row = {'curve': label, 'lam_star': C.lam_star(lam, n), 'K1': K1, 'm': m,
               'eff': eff, 'n': n, **{f'w_{f}': wins[f] for f in FORMS},
               **diag_acc}
        e1_rows.append(row)
        print(f"  {K1:>4} {eff:>7.4f} | " +
              ' '.join(f"{wins[f]:>3}/{len(SEEDS)}" for f in FORMS) +
              f" |  {diag_acc['rho_gh_u']:>7.3f} {diag_acc['rho_gh_c']:>7.3f}")


def wall_of(rows, form, curve):
    """largest eff with >=3/5 successes, and smallest eff with 0/5."""
    rs = [r for r in rows if r['curve'] == curve]
    ok = [r['eff'] for r in rs if r[f'w_{form}'] >= 0.75 * len(SEEDS)]
    bad = [r['eff'] for r in rs if r[f'w_{form}'] == 0]
    return (max(ok) if ok else None, min(bad) if bad else None)


print("\n  Wall summary (last eff with >=75% -> first eff with 0):")
print(f"  {'curve':<14} {'form':<7} {'last_ok':>9} {'first_0':>9}  {'K1 last_ok':>11}")
for label, curve, _k1, m in curves:
    n = curve[2]
    K2 = math.isqrt(n) + 1
    for f in FORMS:
        lo, hi = wall_of(e1_rows, f, label)
        k1s = f"{lo*n/K2:.0f}" if lo else "-"
        print(f"  {label:<14} {f:<7} {lo if lo is None else f'{lo:9.4f}':>9} "
              f"{hi if hi is None else f'{hi:9.4f}':>9}  {k1s:>11}")


# ---------------------------------------------------------------------------
hr('EXP E2 -- is rho_GH = ||error||/GH(L'') the governing quantity?')
# ---------------------------------------------------------------------------
print("rho_GH is lambda-INDEPENDENT (det(L'') has no lam in it).  If the two")
print("curves keep different walls at equal rho_GH, lambda acts through the GS")
print("profile, not the determinant.\n")
print(f"  {'curve':<14} {'K1':>4} {'eff':>7} {'rho_gh_c':>9} {'bdd_c':>6} "
      f"{'mu/e_c':>8} {'lam1/e_c':>9} {'np_marg':>8}")
for r in e1_rows:
    print(f"  {r['curve']:<14} {r['K1']:>4} {r['eff']:>7.4f} {r['rho_gh_c']:>9.3f} "
          f"{r['w_bdd_c']:>3}/{len(SEEDS)} {r['mu']/r['e_c']:>8.3f} "
          f"{r['lam1_ub']/r['e_c']:>9.3f} {r['np_margin']:>8.3f}")


def separation(rows, key, form='bdd_c', invert=False):
    """best single threshold on `key`; returns (acc, thr, n_pos, n_neg)."""
    pos = [r[key] for r in rows if r[f'w_{form}'] >= 0.75 * len(SEEDS)]
    neg = [r[key] for r in rows if r[f'w_{form}'] <= 0.25 * len(SEEDS)]
    if not pos or not neg:
        return None
    cands = sorted(set(pos + neg))
    best = (0.0, None)
    for t in cands:
        if invert:
            acc = (sum(1 for x in pos if x >= t) + sum(1 for x in neg if x < t))
        else:
            acc = (sum(1 for x in pos if x <= t) + sum(1 for x in neg if x > t))
        acc /= (len(pos) + len(neg))
        if acc > best[0]:
            best = (acc, t)
    base = max(len(pos), len(neg)) / (len(pos) + len(neg))
    return best[0], best[1], base, len(pos), len(neg)


print("\n  Single-threshold separation of bdd_c success (>=75%) vs failure (<=25%):")
for key, inv in [('rho_gh_c', False), ('eff', False), ('lam_star', True)]:
    s = separation(e1_rows, key, 'bdd_c', inv)
    if s:
        print(f"    {key:<10} acc={s[0]:.3f}  thr={s[1]:.4f}  "
              f"(majority baseline {s[2]:.3f}; {s[3]} pos / {s[4]} neg)")
for key in ['mu_over_e', 'lam1_over_e']:
    for r in e1_rows:
        r['mu_over_e'] = r['mu'] / r['e_c']
        r['lam1_over_e'] = r['lam1_ub'] / r['e_c']
    s = separation(e1_rows, key, 'bdd_c', invert=True)
    if s:
        print(f"    {key:<10} acc={s[0]:.3f}  thr={s[1]:.4f}  "
              f"(majority baseline {s[2]:.3f}; {s[3]} pos / {s[4]} neg)")


# ---------------------------------------------------------------------------
hr('EXP E3 -- m-dependence.  rho_GH ~ 2.386*sqrt(eff)*n^(1/(2m)) predicts')
# ---------------------------------------------------------------------------
print("that more signatures help only logarithmically (n^(1/2m) -> 1).")
print("T4b (2026-07-29) found m = 8..32 does NOT rescue K1=8 in the `full`")
print("formulation.  Re-test under bdd_c.\n")
M_GRID = [6, 8, 12, 16, 24, 32]
print(f"  {'curve':<14} {'K1':>4} " + ' '.join(f"m={mm:<5}" for mm in M_GRID))
e3_rows = []
for label, curve, _k1, _m in curves:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    for K1 in (8, 16):
        cells = []
        for mm in M_GRID:
            w = 0
            for seed in SEEDS:
                d_sec = secret_for(seed, n)
                sigs = C.gen_signatures(G, d_sec, mm, n, lam, p, K1, K2, seed=seed)
                if len(sigs) < mm:
                    continue
                if attack_bdd(sigs, n, lam, K1, K2, d_sec, centered=True)[0]:
                    w += 1
            sigs = C.gen_signatures(G, 1234, mm, n, lam, p, K1, K2, seed=42)
            dg = diagnostics(sigs, n, lam, K1, K2)
            cells.append((mm, w, dg['rho_gh_c']))
            e3_rows.append({'curve': label, 'K1': K1, 'm': mm, 'w': w,
                            'rho_gh_c': dg['rho_gh_c']})
        print(f"  {label:<14} {K1:>4} " +
              ' '.join(f"{w}/{len(SEEDS)}({rg:.2f})" for _mm, w, rg in cells))


# ---------------------------------------------------------------------------
hr('EXP E4 -- where is the eff wall on fresh curves?  predicted 0.176*n^(-1/m)')
# ---------------------------------------------------------------------------
t0 = time.time()
fresh = C.search_curves(2 ** 16, 2 ** 17, per_bin=2, nbins=6)
print(f"  {len(fresh)} fresh 17-bit j=0 GLV curves "
      f"(lam* spread over 6 bins, {time.time()-t0:.1f}s)")
EFF_GRID = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.55, 0.70]
M_E4 = 12
print(f"\n  m={M_E4}, {len(fresh)} curves x {len(SEEDS)} seeds per eff cell")
print(f"  {'eff':>6} | {'full':>9} {'bdd_c':>9} | {'rho_gh_c':>9}")
e4_rows = []
for eff_t in EFF_GRID:
    tot = wf = wb = 0
    rg = []
    for (p, b, n, lam, G) in fresh:
        K2 = math.isqrt(n) + 1
        K1 = max(2, int(round(eff_t * n / K2)))
        for seed in SEEDS:
            d_sec = secret_for(seed, n)
            sigs = C.gen_signatures(G, d_sec, M_E4, n, lam, p, K1, K2, seed=seed)
            if len(sigs) < M_E4:
                continue
            tot += 1
            if attack_full(sigs, n, lam, K1, K2, d_sec):
                wf += 1
            if attack_bdd(sigs, n, lam, K1, K2, d_sec, centered=True)[0]:
                wb += 1
        sigs = C.gen_signatures(G, 1234, M_E4, n, lam, p, K1, K2, seed=42)
        rg.append(diagnostics(sigs, n, lam, K1, K2)['rho_gh_c'])
    e4_rows.append({'eff': eff_t, 'tot': tot, 'full': wf, 'bdd_c': wb,
                    'rho': sum(rg) / len(rg)})
    print(f"  {eff_t:>6.3f} | {wf:>4}/{tot:<4} {wb:>4}/{tot:<4} | "
          f"{sum(rg)/len(rg):>9.3f}")

hr('SUMMARY')
for label, curve, _k1, m in curves:
    n = curve[2]
    K2 = math.isqrt(n) + 1
    lo_f, _ = wall_of(e1_rows, 'full', label)
    lo_c, _ = wall_of(e1_rows, 'bdd_c', label)
    if lo_f and lo_c:
        print(f"  {label:<14} lam*={C.lam_star(curve[3],n):.4f}  "
              f"eff wall  full {lo_f:.4f} -> bdd_c {lo_c:.4f}   "
              f"({lo_c/lo_f:.2f}x)   K1 wall {lo_f*n/K2:.0f} -> {lo_c*n/K2:.0f}")
print("\nDone.")
