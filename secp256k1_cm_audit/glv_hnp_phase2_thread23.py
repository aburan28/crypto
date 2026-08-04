"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the target IS the
short vector, and find a predictor that actually separates success from
failure.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, commits d525931 /
e845207):

  T5 established that the shortest vector of the Phase-2 lattice is ALWAYS
  the parasitic vector n*S_D*e_m (100% of its energy in the d column).  The
  planted vector is therefore never lambda_1, and recovery is a BDD/coset
  condition, not an SVP condition: the planted vector must be shortest
  among lattice vectors whose last coordinate is +-S_KANNAN.

  That run proposed Thread 23:
    "project the lattice along e_m (quotient out the trivial n*e_m
     direction) and solve BDD in the projection, or replace the Kannan
     embedding with an explicit CVP call.  Falsifier: if sv/pv rises above
     1 after the reformulation AND the K1 wall in T4 moves outward on the
     lam*=0.07 curve (currently K1 ~ 4-6), the reformulation is a real
     improvement; if the wall stays at K1 ~ 4-6, then the wall is
     information-theoretic and Phase 2 is at its ceiling."

Experiments:

  E1  Geometry of the coset sublattice L_0 (= all lattice vectors with last
      coordinate 0; spanned by rows 0..2m).  Recovery is BDD of the target
      in L_0, so the governing quantity should be

          gamma = lambda_1(L_0) / ||v'||

      where v' is the planted vector minus its Kannan coordinate.  Two
      variants: gamma_gh uses the Gaussian heuristic for lambda_1(L_0),
      gamma_emp uses the measured LLL shortest vector of L_0.  Compared
      head-to-head against the three predictors already falsified
      (lam*, rho, eff) on the historical curves.

  E2  Same predictors on a fresh 17-bit multi-curve sweep (AUC + best
      threshold accuracy), to check E1 is not a 3-curve artefact.

  E3  THE THREAD-23 FALSIFIER.  A genuinely reformulated lattice in which
      d is eliminated algebraically (pivot signature 0), so no d column and
      no parasitic n*S_D*e_m vector exists at all.  dim drops 2m+2 -> 2m+1.
      Head-to-head against the original lattice at identical
      (curve, K1, m, seed), plus sv/pv, plus the K1 wall on the lam*=0.07
      curve.

  E4  Closed form for gamma_gh and the eff* wall it predicts:

          gamma_gh = sqrt(3/(2*pi*e)) * (n * eff^m)^(-1/(2m+1)),  eff=K1*K2/n

      so gamma_gh >= 1 requires m >= (ln n + 0.868)/(ln(1/eff) - 1.736),
      which has no solution at all once eff >= exp(-1.736) = 0.1763.
      Verified against the exact determinant, then tested by sweeping m at
      fixed eff on both sides of eff*.

Run: python3 glv_hnp_phase2_thread23.py
"""

import math
import random

from fpylll import IntegerMatrix, LLL, BKZ

from glv_hnp_common import (
    modinv, ec_mul, find_generator, glv_roots, lam_star,
    lambda_block_mu, scales, gen_signatures, build_glv_lattice,
    planted_vector, norm, recover_d, build_curve, search_curves,
)

TWO_PI_E = 2.0 * math.pi * math.e
GH_CONST = math.sqrt(3.0 / TWO_PI_E)          # 0.419766...
EFF_STAR = math.exp(-2.0 * math.log(1.0 / GH_CONST))


# ---------------------------------------------------------------------------
# Geometry of the coset sublattice L_0
# ---------------------------------------------------------------------------

def coset_sublattice(M, m):
    """Rows 0..2m of the Phase-2 lattice restricted to columns 0..2m.

    Those rows all have last coordinate 0, so they span exactly the set of
    lattice vectors with vanishing Kannan coordinate.  Recovery of d is a
    BDD instance in this lattice.
    """
    return [[M[i][j] for j in range(2 * m + 1)] for i in range(2 * m + 1)]


def log_det_L0(n, k1_bound, k2_bound, m):
    """ln det(L_0) = m*ln(n*S_K1) + ln(S_D) + m*ln(S_K2)  (block triangular)."""
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound)
    return m * math.log(n * S_K1) + math.log(S_D) + m * math.log(S_K2)


def gh(log_det, dim):
    return math.sqrt(dim / TWO_PI_E) * math.exp(log_det / dim)


def lll_lambda1(basis):
    """LLL-reduced shortest row norm (an upper bound on lambda_1)."""
    A = IntegerMatrix.from_matrix(basis)
    LLL.reduction(A)
    dim = len(basis)
    best = None
    for i in range(dim):
        row = [A[i][j] for j in range(len(basis[0]))]
        if all(x == 0 for x in row):
            continue
        nr = norm(row)
        if best is None or nr < best:
            best = nr
    return best


def geometry(curve, m, d_secret, k1_bound, seed):
    """All geometric quantities for one (curve, K1, m, seed) instance."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)

    v = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    v_full = norm(v)
    v_proj = norm(v[:2 * m + 1])          # planted vector minus Kannan coord

    L0 = coset_sublattice(M, m)
    lam1_emp = lll_lambda1(L0)
    gh_L0 = gh(log_det_L0(n, k1_bound, k2_bound, m), 2 * m + 1)

    mu, _mu_vec = lambda_block_mu(n, lam, S_K1, S_K2)
    eff = k1_bound * k2_bound / n

    return {
        'n': n, 'lam': lam, 'K1': k1_bound, 'K2': k2_bound, 'm': m,
        'eff': eff, 'lam_star': lam_star(lam, n),
        'v_full': v_full, 'v_proj': v_proj,
        'lam1_emp': lam1_emp, 'gh_L0': gh_L0, 'mu': mu,
        'gamma_emp': lam1_emp / v_proj,
        'gamma_gh': gh_L0 / v_proj,
        'rho': mu / v_full,
        'sigs': sigs, 'M': M, 'S_KAN': S_KAN,
    }


def lll_recovers(M, m, n, S_KAN, d_secret, beta=None):
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    return recover_d(reduced, m, n, S_KAN, d_secret) is not None


# ---------------------------------------------------------------------------
# E3: the d-eliminated ("reformulated") lattice
# ---------------------------------------------------------------------------
#
# From  A_i + B_i*d = k1_i + lam*k2_i  (mod n),  pivot on signature 0:
#     R_i = B_i * B_0^{-1},   C_i = A_i - R_i*A_0
#     k1_i = C_i + R_i*k1_0 + R_i*lam*k2_0 - lam*k2_i   (mod n),  i = 1..m-1
#
# so d never appears.  Free unknowns: k1_0, k2_0, ..., k2_{m-1} (m+1 of them);
# determined unknowns: k1_1..k1_{m-1} (m-1 of them, one column each).
# dim = (m-1) + 1 + m + 1 = 2m+1, one lower than the original, and there is
# no d column, hence no parasitic n*S_D*e_m vector.

def build_elim_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    S_1 = n // k1_bound
    S_2 = max(1, n // k2_bound)
    S_K = n
    dim = 2 * m + 1

    B0inv = modinv(sigs[0]['B'] % n, n)
    R = [0] * m
    C = [0] * m
    for i in range(1, m):
        R[i] = sigs[i]['B'] * B0inv % n
        C[i] = (sigs[i]['A'] - R[i] * sigs[0]['A']) % n

    # column layout
    col_c = {i: i - 1 for i in range(1, m)}        # 0 .. m-2
    col_k1_0 = m - 1
    col_k2 = {i: m + i for i in range(m)}          # m .. 2m-1
    col_kan = 2 * m

    M = [[0] * dim for _ in range(dim)]
    r = 0
    for i in range(1, m):                          # modulus rows
        M[r][col_c[i]] = n * S_1
        r += 1
    for i in range(1, m):                          # k1_0 row coefficients
        M[r][col_c[i]] = R[i] * S_1 % (n * S_1)
    M[r][col_k1_0] = S_1
    r += 1
    for i in range(1, m):                          # k2_0 row coefficients
        M[r][col_c[i]] = R[i] * lam % n * S_1
    M[r][col_k2[0]] = S_2
    r += 1
    for j in range(1, m):                          # k2_j rows
        M[r][col_c[j]] = (-lam) % n * S_1
        M[r][col_k2[j]] = S_2
        r += 1
    for i in range(1, m):                          # target row
        M[r][col_c[i]] = C[i] * S_1
    M[r][col_kan] = S_K
    assert r == dim - 1
    return M, (S_1, S_2, S_K), (col_k1_0, col_k2, col_kan), (R, C)


def elim_planted(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_1 = n // k1_bound
    S_2 = max(1, n // k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(1, m):
        v[i - 1] = sigs[i]['k1'] * S_1
    v[m - 1] = sigs[0]['k1'] * S_1
    for i in range(m):
        v[m + i] = sigs[i]['k2'] * S_2
    v[2 * m] = n
    return v


def elim_recover(reduced, sigs, n, lam, scl, cols, d_secret):
    S_1, S_2, S_K = scl
    col_k1_0, col_k2, col_kan = cols
    B0inv = modinv(sigs[0]['B'] % n, n)
    for row in reduced:
        if abs(row[col_kan]) != S_K:
            continue
        sign = 1 if row[col_kan] > 0 else -1
        if row[col_k1_0] % S_1 or row[col_k2[0]] % S_2:
            continue
        k1_0 = sign * row[col_k1_0] // S_1
        k2_0 = sign * row[col_k2[0]] // S_2
        k0 = (k1_0 + lam * k2_0) % n
        d_cand = (k0 - sigs[0]['A']) * B0inv % n
        if d_cand == d_secret:
            return d_cand
    return None


def elim_experiment(curve, m, d_secret, k1_bound, seed, beta=None):
    """Returns (recovered, sv_over_pv)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None, None
    M, scl, cols, _ = build_elim_lattice(sigs, n, lam, k1_bound, k2_bound)
    dim = 2 * m + 1
    pv = norm(elim_planted(sigs, n, k1_bound, k2_bound))
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    sv = min((norm(r) for r in reduced if any(r)), default=None)
    ok = elim_recover(reduced, sigs, n, lam, scl, cols, d_secret) is not None
    return ok, sv / pv


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------

def auc(rows, key, positive='ok'):
    """Rank-based AUC of `key` as a predictor of `positive` (higher = better)."""
    pos = [r[key] for r in rows if r[positive]]
    neg = [r[key] for r in rows if not r[positive]]
    if not pos or not neg:
        return float('nan')
    wins = 0.0
    for a in pos:
        for b in neg:
            wins += 1.0 if a > b else (0.5 if a == b else 0.0)
    return wins / (len(pos) * len(neg))


def best_threshold(rows, key, positive='ok'):
    """Best accuracy achievable by a single threshold on `key`."""
    vals = sorted({r[key] for r in rows})
    n_tot = len(rows)
    best = (0.0, None)
    for t in vals:
        acc = sum(1 for r in rows if (r[key] >= t) == bool(r[positive])) / n_tot
        if acc > best[0]:
            best = (acc, t)
    return best


def majority(rows, positive='ok'):
    k = sum(1 for r in rows if r[positive])
    return max(k, len(rows) - k) / len(rows)


# ===========================================================================

print("=" * 78)
print("Thread 23 — reformulate the Phase-2 lattice so the target is findable")
print("=" * 78)
print(f"GH constant sqrt(3/(2*pi*e)) = {GH_CONST:.5f}")
print(f"predicted hard wall eff* = exp(-2*ln(1/GH_CONST)) = {EFF_STAR:.5f}")

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
    assert G is not None and (lam * lam + lam + 1) % n == 0
    hist.append((label, (p, b, n, lam, G), k1, m))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E1: BDD geometry of the coset sublattice L_0, historical curves")
print("-" * 78)
print("L_0 = lattice vectors with Kannan coordinate 0 (rows 0..2m, cols 0..2m).")
print("v' = planted vector minus its Kannan coordinate.  Recovery is BDD of")
print("the target in L_0, so the governing ratio should be lambda_1(L_0)/||v'||.\n")
print(f"{'curve':<18} {'K1':>3} {'eff':>6} {'lam*':>6} {'g_emp':>7} {'g_gh':>6} "
      f"{'rho':>6} {'ok':>5}")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_E1 = 12
e1_rows = []
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    for k1 in K1_GRID:
        succ = 0
        rec = None
        for seed in SEEDS:
            d_secret = random.Random(seed * 7919).randrange(2, n - 1)
            g = geometry(curve, M_E1, d_secret, k1, seed)
            if g is None:
                continue
            rec = g
            if lll_recovers(g['M'], M_E1, n, g['S_KAN'], d_secret):
                succ += 1
        if rec is None:
            continue
        row = {k: rec[k] for k in
               ('eff', 'lam_star', 'gamma_emp', 'gamma_gh', 'rho')}
        row.update(label=label, K1=k1, succ=succ, ok=(succ >= 3))
        e1_rows.append(row)
        print(f"{label:<18} {k1:>3} {rec['eff']:>6.3f} {rec['lam_star']:>6.3f} "
              f"{rec['gamma_emp']:>7.3f} {rec['gamma_gh']:>6.3f} "
              f"{rec['rho']:>6.3f} {succ:>3}/5")

print("\nPredictor comparison on E1 (ok = >=3/5 seeds recovered):")
print(f"  majority baseline           {majority(e1_rows):.3f}")
for key in ('gamma_emp', 'gamma_gh', 'eff', 'lam_star', 'rho'):
    sign = -1 if key == 'eff' else 1
    tmp = [dict(r, _k=sign * r[key]) for r in e1_rows]
    acc, thr = best_threshold(tmp, '_k')
    a = auc(tmp, '_k')
    lbl = f"-{key}" if sign < 0 else key
    print(f"  {lbl:<26} AUC={a:.3f}  best-threshold acc={acc:.3f}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E2: same predictors on a fresh 17-bit multi-curve sweep")
print("-" * 78)
LO, HI = 2 ** 16, 2 ** 17
M_E2 = 12
print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
curves = search_curves(LO, HI, per_bin=2, nbins=8)
print(f"Found {len(curves)} curves.  m={M_E2}, seeds={SEEDS[:3]}\n")

e2_rows = []
for eff_target in (0.05, 0.12, 0.18, 0.25):
    for (p, b, n, lam, G) in curves:
        k2_bound = math.isqrt(n) + 1
        k1 = max(2, int(eff_target * n / k2_bound))
        curve = (p, b, n, lam, G)
        succ = 0
        rec = None
        for seed in SEEDS[:3]:
            d_secret = random.Random(seed * 7919).randrange(2, n - 1)
            g = geometry(curve, M_E2, d_secret, k1, seed)
            if g is None:
                continue
            rec = g
            if lll_recovers(g['M'], M_E2, n, g['S_KAN'], d_secret):
                succ += 1
        if rec is None:
            continue
        row = {k: rec[k] for k in
               ('eff', 'lam_star', 'gamma_emp', 'gamma_gh', 'rho')}
        row.update(n=n, K1=k1, succ=succ, ok=(succ >= 2), eff_target=eff_target)
        e2_rows.append(row)

for eff_target in (0.05, 0.12, 0.18, 0.25):
    sub = [r for r in e2_rows if r['eff_target'] == eff_target]
    nok = sum(1 for r in sub if r['ok'])
    gg = [r['gamma_gh'] for r in sub]
    ge = [r['gamma_emp'] for r in sub]
    print(f"  eff~{eff_target:.2f}: {nok:>2}/{len(sub)} curves recover  "
          f"gamma_gh in [{min(gg):.3f},{max(gg):.3f}]  "
          f"gamma_emp in [{min(ge):.3f},{max(ge):.3f}]")

print("\nPredictor comparison on E2 (pooled, ok = >=2/3 seeds):")
print(f"  majority baseline           {majority(e2_rows):.3f}")
for key in ('gamma_emp', 'gamma_gh', 'eff', 'lam_star', 'rho'):
    sign = -1 if key == 'eff' else 1
    tmp = [dict(r, _k=sign * r[key]) for r in e2_rows]
    acc, thr = best_threshold(tmp, '_k')
    a = auc(tmp, '_k')
    lbl = f"-{key}" if sign < 0 else key
    print(f"  {lbl:<26} AUC={a:.3f}  best-threshold acc={acc:.3f}  thr={thr:.4f}")

pooled = e1_rows + e2_rows
print("\nPooled E1+E2:")
for key in ('gamma_emp', 'gamma_gh'):
    tmp = [dict(r, _k=r[key]) for r in pooled]
    acc, thr = best_threshold(tmp, '_k')
    print(f"  {key:<26} AUC={auc(tmp,'_k'):.3f}  acc={acc:.3f}  thr={thr:.4f}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E3: THE FALSIFIER — d-eliminated reformulation, head to head")
print("-" * 78)
print("Reformulated lattice: d removed algebraically (pivot on signature 0),")
print("dim 2m+2 -> 2m+1, no d column, so the parasitic n*S_D*e_m vector that")
print("T5 (2026-07-29) identified as lambda_1 cannot exist.")
print("Falsifier: does sv/pv rise above 1, and does the K1 wall move?\n")
print(f"{'curve':<18} {'K1':>3} {'eff':>6} | {'orig ok':>8} {'sv/pv':>7} | "
      f"{'elim ok':>8} {'sv/pv':>7}")

M_E3 = 12
e3_rows = []
for label, curve, _k1, _m in hist:
    p, b, n, lam, G = curve
    for k1 in K1_GRID:
        o_succ = e_succ = 0
        o_svpv = e_svpv = []
        o_svpv, e_svpv = [], []
        for seed in SEEDS:
            d_secret = random.Random(seed * 7919).randrange(2, n - 1)
            g = geometry(curve, M_E3, d_secret, k1, seed)
            if g is None:
                continue
            if lll_recovers(g['M'], M_E3, n, g['S_KAN'], d_secret):
                o_succ += 1
            # original sv/pv
            A = IntegerMatrix.from_matrix(g['M'])
            LLL.reduction(A)
            dim = 2 * M_E3 + 2
            sv = min(norm([A[i][j] for j in range(dim)])
                     for i in range(dim)
                     if any(A[i][j] for j in range(dim)))
            o_svpv.append(sv / g['v_full'])
            ok, svpv = elim_experiment(curve, M_E3, d_secret, k1, seed)
            if ok:
                e_succ += 1
            if svpv is not None:
                e_svpv.append(svpv)
        if not o_svpv:
            continue
        eff = k1 * (math.isqrt(n) + 1) / n
        o_m = sum(o_svpv) / len(o_svpv)
        e_m = sum(e_svpv) / len(e_svpv)
        e3_rows.append(dict(label=label, K1=k1, eff=eff, orig=o_succ,
                            elim=e_succ, o_svpv=o_m, e_svpv=e_m))
        print(f"{label:<18} {k1:>3} {eff:>6.3f} | {o_succ:>6}/5 {o_m:>7.3f} | "
              f"{e_succ:>6}/5 {e_m:>7.3f}")

print("\nK1 wall (largest K1 with >=3/5 recovery):")
for label, _c, _k, _m in hist:
    sub = [r for r in e3_rows if r['label'] == label]
    wo = max([r['K1'] for r in sub if r['orig'] >= 3], default=0)
    we = max([r['K1'] for r in sub if r['elim'] >= 3], default=0)
    print(f"  {label:<18} original K1<={wo:<3} reformulated K1<={we:<3} "
          f"{'MOVED' if we > wo else 'unchanged' if we == wo else 'WORSE'}")

n_svpv_gt1 = sum(1 for r in e3_rows if r['e_svpv'] > 1.0)
print(f"\nsv/pv > 1 in the reformulated lattice: {n_svpv_gt1}/{len(e3_rows)} "
      f"configurations  (original: "
      f"{sum(1 for r in e3_rows if r['o_svpv'] > 1.0)}/{len(e3_rows)})")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E4: closed form for gamma_gh and the eff* wall")
print("-" * 78)
print("Claim:  gamma_gh = sqrt(3/(2*pi*e)) * (n * eff^m)^(-1/(2m+1))")
print("Check against the exact determinant of L_0:\n")
print(f"{'n':>7} {'K1':>3} {'m':>3} {'exact':>8} {'closed':>8} {'rel.err':>9}")
for (nn, k1v, mv) in [(2659, 4, 12), (2659, 8, 12), (2647, 6, 12),
                      (2647, 4, 20), (199, 2, 6), (100003, 10, 16)]:
    k2v = math.isqrt(nn) + 1
    eff = k1v * k2v / nn
    exact = gh(log_det_L0(nn, k1v, k2v, mv), 2 * mv + 1) / \
        (nn * math.sqrt((2 * mv + 1) / 3.0))
    closed = GH_CONST * (nn * eff ** mv) ** (-1.0 / (2 * mv + 1))
    print(f"{nn:>7} {k1v:>3} {mv:>3} {exact:>8.4f} {closed:>8.4f} "
          f"{abs(exact-closed)/closed:>8.2%}")

print("\nm required for gamma_gh >= 1:  m >= (ln n + 0.868)/(ln(1/eff) - 1.736)")
print(f"(no solution once eff >= eff* = {EFF_STAR:.4f})\n")
print(f"{'eff':>6} {'ln(1/eff)':>10} {'m_required':>12}   {'naive m_thresh':>14}")
for eff in (0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.1763, 0.20, 0.25):
    nn = 2659
    denom = math.log(1.0 / eff) - 2.0 * math.log(1.0 / GH_CONST)
    naive = math.ceil(math.log(nn) / math.log(1.0 / eff))
    req = "infeasible" if denom <= 0 else \
        f"{math.ceil((math.log(nn) + math.log(1.0/GH_CONST))/denom):d}"
    print(f"{eff:>6.4f} {math.log(1.0/eff):>10.3f} {req:>12}   {naive:>14}")

print("\nE4b: does more data help?  m-sweep at fixed eff, curve 12-bit/2677")
curve2677 = hist[2][1]
n2677 = curve2677[2]
k2b = math.isqrt(n2677) + 1
print(f"{'eff':>6} {'K1':>3} " + " ".join(f"m={mm:<4}" for mm in (8, 12, 16, 24, 32)))
for k1v in (2, 4, 6, 8):
    eff = k1v * k2b / n2677
    cells = []
    for mm in (8, 12, 16, 24, 32):
        s = 0
        for seed in SEEDS[:3]:
            d_secret = random.Random(seed * 7919).randrange(2, n2677 - 1)
            k2_bound = k2b
            sigs = gen_signatures(curve2677[4], d_secret, mm, n2677,
                                  curve2677[3], curve2677[0], k1v, k2_bound, seed)
            if len(sigs) < mm:
                continue
            M = build_glv_lattice(sigs, n2677, curve2677[3], k1v, k2_bound)
            _, _, _, S_KAN = scales(n2677, k1v, k2_bound)
            if lll_recovers(M, mm, n2677, S_KAN, d_secret):
                s += 1
        cells.append(f"{s}/3 ")
    print(f"{eff:>6.3f} {k1v:>3} " + " ".join(f"{c:<6}" for c in cells))


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E5: is eff* a real ceiling, or just an LLL limitation?")
print("-" * 78)
print("As m -> inf the closed form has the finite limit")
print(f"    gamma_inf = {GH_CONST:.5f} / sqrt(eff)")
print("so eff* is where gamma_inf = 1: beyond it the planted vector is longer")
print("than the Gaussian-heuristic lambda_1(L_0) NO MATTER how much data is")
print("used, and no amount of reduction effort can help.  Sharp prediction:")
print("BKZ can push the achievable eff up towards eff* but never past it.\n")
print(f"{'K1':>3} {'eff':>6} {'g_inf':>6} {'m':>3}  " +
      "  ".join(f"{lbl:>7}" for lbl in ("LLL", "BKZ20", "BKZ30", "BKZ40")))

for k1v in (6, 8, 12, 16):
    eff = k1v * k2b / n2677
    g_inf = GH_CONST / math.sqrt(eff)
    for mm in (12, 20):
        cells = []
        for beta in (None, 20, 30, 40):
            s = 0
            for seed in SEEDS[:3]:
                d_secret = random.Random(seed * 7919).randrange(2, n2677 - 1)
                sigs = gen_signatures(curve2677[4], d_secret, mm, n2677,
                                      curve2677[3], curve2677[0], k1v, k2b, seed)
                if len(sigs) < mm:
                    continue
                M = build_glv_lattice(sigs, n2677, curve2677[3], k1v, k2b)
                _, _, _, S_KAN = scales(n2677, k1v, k2b)
                b_use = None if beta is None else min(beta, 2 * mm + 1)
                if lll_recovers(M, mm, n2677, S_KAN, d_secret, beta=b_use):
                    s += 1
            cells.append(f"{s}/3")
        print(f"{k1v:>3} {eff:>6.3f} {g_inf:>6.3f} {mm:>3}  " +
              "  ".join(f"{c:>7}" for c in cells))

print(f"\nPrediction: rows with gamma_inf < 1 (eff > eff* = {EFF_STAR:.4f}) stay 0/3")
print("across every block size and every m.  Rows with gamma_inf > 1 may be")
print("rescued by larger beta.")

print("\nDone.")
