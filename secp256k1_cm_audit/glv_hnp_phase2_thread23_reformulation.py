"""
GLV-HNP Thread 23 — reformulate the Phase-2 lattice so the planted vector is
the shortest vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5)
--------------------------------------------------------
The Phase-2 lattice built by `glv_hnp_phase2_lambda_threshold.py:254`
(`build_glv_lattice`) is (2m+2)-dimensional with column scales

    S_K1 = n // K1,   S_D = 1,   S_K2 = n // K2,   S_KANNAN = n

and columns [ k1_0..k1_{m-1} | d | k2_0..k2_{m-1} | kannan ].  T5 showed that on
*every* curve tested the shortest vector after LLL is the trivial vector

    v_triv = n * S_D * e_m         (norm n*S_D)

which carries no information (d is only defined mod n).  Meanwhile

    ||v_planted||^2 ~ n^2 (2m/3 + 4/3)

so v_planted is never lambda_1 and recovery is a *coset/BDD* condition, not an
SVP condition.  T5's proposed remedy (Thread 23) is implemented and tested here.

The three formulations compared
-------------------------------
V0  baseline   (2m+2)-dim Kannan lattice, verbatim from the T4 run.
V1  quotient   drop the d column entirely -> (2m+1) columns, rank 2m+1.
               This is exactly L / <n*S_D*e_m>, so v_triv is quotiented away.
               d is no longer a coordinate; it is recovered arithmetically from
               the (k1_0, k2_0) block of the surviving Kannan row via
               d = (k1_0 - A_0 + lam*k2_0) * B_0^{-1} mod n.
V2  CVP        drop the d column AND the Kannan row -> 2m columns, rank 2m
               homogeneous lattice; solve the CVP directly with Babai
               nearest-plane against target -A_vec = (-A_i*S_K1, 0).
               Same d-recovery step as V1.

All three use identical curves, identical signatures (same seeds), identical
scales, so the comparison is exact.

Falsifier (stated in the 2026-07-29 log entry)
----------------------------------------------
  "if sv/pv rises above 1 after the reformulation and the K1 wall in T4 moves
   outward on the lam*=0.07 curve (currently K1~4-6), the reformulation is a
   real improvement; if the wall stays at K1~4-6, then the wall is
   information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_thread23_reformulation.py
"""

import math
import os
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Load the Phase-2 primitives VERBATIM from the Thread-20 script.
# That file has no __main__ guard (it runs its own experiments at import), so
# we exec only the prefix that precedes its first top-level print statement.
# ---------------------------------------------------------------------------

_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "glv_hnp_phase2_lambda_threshold.py")
with open(_SRC) as _fh:
    _lines = _fh.readlines()
_cut = next(i for i, l in enumerate(_lines) if l.startswith('print('))
_ns = {'__name__': '_t20_prims', '__file__': _SRC}
exec("".join(_lines[:_cut]), _ns)

modinv = _ns['modinv']
scales = _ns['scales']
gen_signatures = _ns['gen_signatures']
build_glv_lattice = _ns['build_glv_lattice']
planted_vector = _ns['planted_vector']
recover_d = _ns['recover_d']
norm = _ns['norm']
find_generator = _ns['find_generator']
lam_star = _ns['lam_star']
search_curves = _ns['search_curves']
lambda_block_mu = _ns['lambda_block_mu']

# ---------------------------------------------------------------------------
# V1 — quotient lattice (d column dropped)
# ---------------------------------------------------------------------------

def build_quotient_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Same rows as build_glv_lattice with the d column deleted.
    Columns: [k1_0..k1_{m-1} | k2_0..k2_{m-1} | kannan],  dim = 2m+1.
    Rows: 2m+2 generators of rank 2m+1 (the dependency is exactly the
    quotiented-out direction n*e_d)."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(2 * m + 2)]
    for i in range(m):                       # n*S_K1*e_i
        M[i][i] = n * S_K1
    for i in range(m):                       # d-row (no d marker any more)
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):                       # k2 rows
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    for i in range(m):                       # Kannan row
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M


def quotient_planted_vector(sigs, n, k1_bound, k2_bound):
    """Image of v_planted in the quotient: (k1_i*S_K1, k2_i*S_K2, S_KANNAN)."""
    m = len(sigs)
    S_K1, _S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v


def d_from_residual(u, sigs, n, lam, k1_bound, k2_bound):
    """u = (k1_i*S_K1, k2_i*S_K2, ...).  Recover d from equation 0:
         k1 = A + B*d - lam*k2  (mod n)  =>  d = (k1 - A + lam*k2) * B^{-1}."""
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    if u[0] % S_K1 or u[m] % S_K2:
        return None
    k1_0, k2_0 = u[0] // S_K1, u[m] // S_K2
    B = sigs[0]['B'] % n
    if B == 0:
        return None
    return (k1_0 - sigs[0]['A'] + lam * k2_0) * modinv(B, n) % n


def recover_d_quotient(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Scan reduced rows for a Kannan row (|last| == S_KANNAN) and decode."""
    m = len(sigs)
    dim = 2 * m + 1
    _S_K1, _S_D, _S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        if abs(row[dim - 1]) != S_KAN:
            continue
        u = [x if row[dim - 1] > 0 else -x for x in row]
        cand = d_from_residual(u, sigs, n, lam, k1_bound, k2_bound)
        if cand is not None and cand == d_secret:
            return cand
    return None


# ---------------------------------------------------------------------------
# V2 — homogeneous lattice + Babai nearest plane (true CVP, no Kannan row)
# ---------------------------------------------------------------------------

def build_homogeneous_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Columns [k1_0..k1_{m-1} | k2_0..k2_{m-1}], dim 2m, 2m+1 generators."""
    m = len(sigs)
    dim = 2 * m
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M = [[0] * dim for _ in range(2 * m + 1)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + i] = S_K2
    return M


def gram_schmidt(basis):
    """Float GS of a list of integer rows.  Returns (bstar, mu)."""
    bstar, mu = [], []
    for v in basis:
        w = [float(x) for x in v]
        coeffs = []
        for bs in bstar:
            nb = sum(x * x for x in bs)
            c = (sum(a * b for a, b in zip(v, bs)) / nb) if nb > 0 else 0.0
            coeffs.append(c)
            for j in range(len(w)):
                w[j] -= c * bs[j]
        bstar.append(w)
        mu.append(coeffs)
    return bstar, mu


def babai_nearest_plane(basis, target):
    """Babai nearest plane.  basis: LLL-reduced integer rows.  Returns the
    lattice vector (exact integers) closest to `target`."""
    bstar, _ = gram_schmidt(basis)
    b = [float(x) for x in target]
    coeffs = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        nb = sum(x * x for x in bstar[i])
        if nb == 0:
            continue
        c = sum(x * y for x, y in zip(b, bstar[i])) / nb
        q = int(math.floor(c + 0.5))
        coeffs[i] = q
        for j in range(len(b)):
            b[j] -= q * basis[i][j]
    out = [0] * len(target)
    for i, q in enumerate(coeffs):
        if q:
            for j in range(len(out)):
                out[j] += q * basis[i][j]
    return out


# ---------------------------------------------------------------------------
# Trial drivers.  All three variants share signatures for a given (curve,seed).
# ---------------------------------------------------------------------------

def reduce_rows(M, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]


def one_trial(curve, m, d_secret, k1_bound, seed, use_bkz=False, beta=20):
    """Returns dict with per-variant success flags and geometry, or None."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    out = {}

    # --- V0 baseline -------------------------------------------------------
    rows0 = reduce_rows(build_glv_lattice(sigs, n, lam, k1_bound, k2_bound),
                        use_bkz, beta)
    out['v0'] = recover_d(rows0, m, n, S_KAN, d_secret) is not None
    pv0 = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
    sv0 = min(norm(r) for r in rows0 if any(r))
    out['svpv0'] = sv0 / pv0

    # --- V1 quotient -------------------------------------------------------
    rows1 = reduce_rows(build_quotient_lattice(sigs, n, lam, k1_bound, k2_bound),
                        use_bkz, beta)
    out['v1'] = recover_d_quotient(rows1, sigs, n, lam, k1_bound, k2_bound,
                                   d_secret) is not None
    pv1 = norm(quotient_planted_vector(sigs, n, k1_bound, k2_bound))
    sv1 = min(norm(r) for r in rows1 if any(r))
    out['svpv1'] = sv1 / pv1
    out['pv1'] = pv1

    # --- V2 CVP ------------------------------------------------------------
    rowsH = reduce_rows(build_homogeneous_lattice(sigs, n, lam, k1_bound,
                                                  k2_bound), use_bkz, beta)
    basis = [r for r in rowsH if any(r)]
    target = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * m
    w = babai_nearest_plane(basis, target)
    u = [w[j] - target[j] for j in range(2 * m)]   # = (k1_i*S_K1, k2_i*S_K2)
    ok2 = False
    for cand_u in (u, [-x for x in u]):
        cand = d_from_residual(cand_u, sigs, n, lam, k1_bound, k2_bound)
        if cand is not None and cand == d_secret:
            ok2 = True
            break
    out['v2'] = ok2
    out['babai_res'] = norm(u)

    # --- geometry: Gaussian heuristic for the quotient lattice -------------
    # det L' = (n*S_K1)^m * S_K2^m * S_KANNAN / n     (rank 2m+1)
    dim1 = 2 * m + 1
    logdet = (m * math.log(n * S_K1) + m * math.log(S_K2)
              + math.log(S_KAN) - math.log(n))
    gh = math.sqrt(dim1 / (2 * math.pi * math.e)) * math.exp(logdet / dim1)
    out['gh1'] = gh
    out['gap1'] = gh / pv1            # >1 predicts planted is below lambda_1^GH
    out['sv1'] = sv1                  # absolute shortest length (for EXP F)
    return out


def sweep(curve, m, k1_bound, seeds, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    agg = {'v0': 0, 'v1': 0, 'v2': 0, 'tot': 0,
           'svpv0': [], 'svpv1': [], 'gap1': [], 'sv1': [], 'pv1': []}
    for seed in seeds:
        d_trial = random.Random(seed + 7777).randint(1, n - 1)
        r = one_trial(curve, m, d_trial, k1_bound, seed, use_bkz, beta)
        if r is None:
            continue
        agg['tot'] += 1
        for k in ('v0', 'v1', 'v2'):
            agg[k] += 1 if r[k] else 0
        for k in ('svpv0', 'svpv1', 'gap1', 'sv1', 'pv1'):
            agg[k].append(r[k])
    return agg


def mean(xs):
    return sum(xs) / len(xs) if xs else float('nan')


# ===========================================================================
print("=" * 78)
print("Thread 23 — does quotienting out the trivial vector fix the Phase-2")
print("            lattice?  (V0 baseline / V1 quotient / V2 CVP-Babai)")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]          # verbatim from the T4 run
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
M_SIGS = 12

HIST = [
    # label,             p,    b, n,    lam
    ("12-bit/2557",      2557, 2, 2659, 1755),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185),
]
hist_curves = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None, f"no generator for {label}"
    assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
    hist_curves.append((label, (p, b, n, lam, G)))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP A: K1 wall, three formulations, m = 12, 5 seeds")
print("-" * 78)
print("V0 row must reproduce the 2026-07-29 T4 table exactly (sanity check).")
print("K2 = 52 for both curves; eff = K1*K2/n.\n")

expA = {}
for label, curve in hist_curves:
    p, b, n, lam, G = curve
    print(f"{label}   lam* = {lam_star(lam, n):.4f}   n = {n}")
    print(f"  {'var':<4} " + " ".join(f"K1={k:<4}" for k in K1_GRID))
    per_var = {'v0': [], 'v1': [], 'v2': []}
    stats = {}
    for k1 in K1_GRID:
        agg = sweep(curve, M_SIGS, k1, SEEDS)
        stats[k1] = agg
        for v in ('v0', 'v1', 'v2'):
            per_var[v].append(f"{agg[v]}/{agg['tot']}  ")
    for v in ('v0', 'v1', 'v2'):
        print(f"  {v:<4} " + " ".join(per_var[v]))
    expA[label] = stats
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP B: is the planted vector now the shortest vector?")
print("-" * 78)
print("sv/pv = ||shortest reduced vector|| / ||planted vector||.")
print("sv/pv >= 1  <=>  nothing in the lattice is shorter than the planted")
print("vector, i.e. recovery has become a genuine SVP problem.\n")
print(f"{'curve':<18} {'K1':>4} {'eff':>6} {'sv/pv V0':>9} {'sv/pv V1':>9} "
      f"{'GH/pv V1':>9} {'V0':>5} {'V1':>5} {'V2':>5}")
for label, curve in hist_curves:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    for k1 in K1_GRID:
        a = expA[label][k1]
        print(f"{label:<18} {k1:>4} {k1*k2/n:>6.3f} "
              f"{mean(a['svpv0']):>9.3f} {mean(a['svpv1']):>9.3f} "
              f"{mean(a['gap1']):>9.3f} "
              f"{str(a['v0'])+'/'+str(a['tot']):>5} "
              f"{str(a['v1'])+'/'+str(a['tot']):>5} "
              f"{str(a['v2'])+'/'+str(a['tot']):>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP C: does the Gaussian-heuristic gap predict success?")
print("-" * 78)
print("gap = lambda_1^GH(L') / ||v'_planted||, L' = quotient lattice (rank 2m+1).")
print("Prediction rule: success iff gap >= 1 (planted is below the generic")
print("shortest length, so it is findable).  Scored against V1 outcomes.\n")

rows_all = []
for label, curve in hist_curves:
    for k1 in K1_GRID:
        a = expA[label][k1]
        if a['tot'] == 0:
            continue
        rows_all.append({'label': label, 'k1': k1, 'gap': mean(a['gap1']),
                         'ok_v0': a['v0'] == a['tot'],
                         'ok_v1': a['v1'] == a['tot'],
                         'ok_v2': a['v2'] == a['tot']})

for key in ('ok_v0', 'ok_v1', 'ok_v2'):
    tp = sum(1 for r in rows_all if r['gap'] >= 1.0 and r[key])
    tn = sum(1 for r in rows_all if r['gap'] < 1.0 and not r[key])
    fp = sum(1 for r in rows_all if r['gap'] >= 1.0 and not r[key])
    fn = sum(1 for r in rows_all if r['gap'] < 1.0 and r[key])
    acc = (tp + tn) / len(rows_all) if rows_all else float('nan')
    print(f"  gap>=1 vs {key:<6}  TP={tp} TN={tn} FP={fp} FN={fn}  acc={acc:.3f}")

print("\n  per-cell detail:")
print(f"  {'curve':<18} {'K1':>4} {'gap':>7} {'pred':>6} {'V0':>5} {'V1':>5} {'V2':>5}")
for r in rows_all:
    print(f"  {r['label']:<18} {r['k1']:>4} {r['gap']:>7.3f} "
          f"{('win' if r['gap'] >= 1 else 'lose'):>6} "
          f"{str(r['ok_v0']):>5} {str(r['ok_v1']):>5} {str(r['ok_v2']):>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP D: are V0 and V1 the same experiment?  (per-seed agreement)")
print("-" * 78)
print("If the quotient map is injective on the Kannan coset then V0 and V1 must")
print("agree trial by trial, not merely in aggregate counts.\n")

agree = disagree = 0
mismatches = []
for label, curve in hist_curves:
    p, b, n, lam, G = curve
    for k1 in K1_GRID:
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            r = one_trial(curve, M_SIGS, d_trial, k1, seed)
            if r is None:
                continue
            if r['v0'] == r['v1']:
                agree += 1
            else:
                disagree += 1
                mismatches.append((label, k1, seed, r['v0'], r['v1']))
print(f"  per-trial agreement V0 == V1: {agree}/{agree+disagree}")
if mismatches:
    for mm in mismatches:
        print(f"    mismatch {mm}")
else:
    print("  no mismatches: the quotient changes the lattice but not the outcome.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP E: the gap law as a falsifiable prediction — m-scaling at fixed eff")
print("-" * 78)
print("Closed form for the quotient lattice (S_K1=n/K1, S_K2=n/K2, S_KAN=n):")
print("    det L' = n^(2m) / eff^m,   dim = 2m+1,   ||v'|| = n*sqrt(2m/3+1)")
print("    gap(m, eff, n) = sqrt((2m+1)/(2*pi*e*(2m/3+1)))")
print("                     * n^(-1/(2m+1)) * eff^(-m/(2m+1))")
print("As m -> infinity this tends to eff^(-1/2)*sqrt(3/(2*pi*e)), so gap >= 1")
print("iff eff <= 3/(2*pi*e) = %.4f -- a curve-INDEPENDENT ceiling strictly"
      % (3.0 / (2 * math.pi * math.e)))
print("below the information-theoretic bound eff <= n^(-1/m) -> 1.\n")
print("2026-07-29 T4b tested m = 8,12,16,24,32 at K1=8 (eff=0.157) on the")
print("lam*=0.07 curve, all failed, and concluded 'the K1 wall is genuine, not")
print("sample size'.  The gap law says gap(32) = %.3f < 1 so those m were ALL"
      % (math.sqrt(65 / (2 * math.pi * math.e * (64 / 3 + 1)))
         * 2647 ** (-1 / 65) * 0.157 ** (-32 / 65)))
print("too small, and predicts the crossing near m ~ 100.  Testing that now.\n")

def gap_closed_form(m, n, k1, k2):
    eff = k1 * k2 / n
    return (math.sqrt((2 * m + 1) / (2 * math.pi * math.e * (2 * m / 3.0 + 1)))
            * n ** (-1.0 / (2 * m + 1)) * eff ** (-m / (2.0 * m + 1)))

fail_curve = hist_curves[1][1]           # 12-bit/2677, lam* = 0.07
p_f, b_f, n_f, lam_f, G_f = fail_curve
K1_FIXED = 8
K2_f = math.isqrt(n_f) + 1
M_GRID = [12, 32, 48, 64, 80, 100, 120]
E_SEEDS = [42, 1234, 9999]

print(f"  curve n = {n_f}, K1 = {K1_FIXED}, K2 = {K2_f}, "
      f"eff = {K1_FIXED*K2_f/n_f:.3f}, lam* = {lam_star(lam_f, n_f):.4f}")
print(f"  {'m':>5} {'dim':>5} {'gap_cf':>7} {'gap_obs':>8} {'pred':>6} "
      f"{'V1':>5} {'sv/pv':>7}")
expE = {}
for m_try in M_GRID:
    agg = sweep(fail_curve, m_try, K1_FIXED, E_SEEDS)
    expE[m_try] = agg
    g_cf = gap_closed_form(m_try, n_f, K1_FIXED, K2_f)
    print(f"  {m_try:>5} {2*m_try+1:>5} {g_cf:>7.3f} {mean(agg['gap1']):>8.3f} "
          f"{('win' if g_cf >= 1 else 'lose'):>6} "
          f"{str(agg['v1'])+'/'+str(agg['tot']):>5} {mean(agg['svpv1']):>7.3f}")
    sys.stdout.flush()

print("\nIf V1 turns from 0/3 to 3/3 somewhere near the gap_cf = 1 crossing, the")
print("gap law is confirmed and the 2026-07-29 'not sample size' reading of T4b")
print("is corrected: it was sample size, just far more of it than was tested.")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP F: why the Gaussian heuristic does not apply — the m-independent")
print("       lambda-block parasite")
print("-" * 78)
print("The quotient lattice contains m disjoint copies of the 2-D sublattice")
print("    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >   (columns k1_i, k2_i)")
print("whose exact shortest vector mu = lambda_1(L2) is computed in O(log n) by")
print("Lagrange-Gauss (`lambda_block_mu`, glv_hnp_phase2_lambda_threshold.py:209)")
print("and does NOT depend on m.  If the observed shortest vector equals mu then")
print("lambda_1(L') is pinned at an m-independent value while ||v'|| grows like")
print("sqrt(m), so the Gaussian heuristic over-estimates lambda_1 and the gap law")
print("of EXP E must fail.  mu is the same quantity that gives the ratio")
print("nu_hat = mu/sqrt(det L2), the 2026-07-29 separator (AUC 0.935).\n")

S_K1_f, _sd, S_K2_f, _sk = scales(n_f, K1_FIXED, K2_f)
mu_f, _w = lambda_block_mu(n_f, lam_f, S_K1_f, S_K2_f)
nu_hat_f = mu_f / math.sqrt(n_f * S_K1_f * S_K2_f)
print(f"  fail curve, K1 = {K1_FIXED}:  mu = {mu_f:.1f}  "
      f"nu_hat = {nu_hat_f:.4f}   (mu is independent of m)")
print(f"  {'m':>5} {'sv_obs':>10} {'mu':>10} {'sv/mu':>7} {'||v'+chr(39)+'||':>10} "
      f"{'sv*sqrt(m)':>11} {'mu/||v'+chr(39)+'||':>11} {'V1':>5}")
for m_try in M_GRID:
    a = expE[m_try]
    sv = mean(a['sv1'])
    pv = mean(a['pv1'])
    print(f"  {m_try:>5} {sv:>10.1f} {mu_f:>10.1f} {sv/mu_f:>7.3f} "
          f"{pv:>10.1f} {sv*math.sqrt(m_try):>11.1f} {mu_f/pv:>11.3f} "
          f"{str(a['v1'])+'/'+str(a['tot']):>5}")

print("\n  Same check across the EXP A K1 grid at m = 12 (both curves):")
print(f"  {'curve':<18} {'K1':>4} {'sv_obs':>10} {'mu':>10} {'sv/mu':>7} "
      f"{'nu_hat':>7} {'mu/||v'+chr(39)+'||':>11} {'V1':>5}")
for label, curve in hist_curves:
    p, b, n, lam, G = curve
    k2 = math.isqrt(n) + 1
    for k1 in K1_GRID:
        a = expA[label][k1]
        S1, _s, S2, _k = scales(n, k1, k2)
        mu, _ = lambda_block_mu(n, lam, S1, S2)
        nh = mu / math.sqrt(n * S1 * S2)
        sv, pv = mean(a['sv1']), mean(a['pv1'])
        print(f"  {label:<18} {k1:>4} {sv:>10.1f} {mu:>10.1f} {sv/mu:>7.3f} "
              f"{nh:>7.4f} {mu/pv:>11.3f} "
              f"{str(a['v1'])+'/'+str(a['tot']):>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP G: out-of-sample test of  lambda_1(L') = min(mu, ||v'||)")
print("-" * 78)
print("EXP F saw sv = mu exactly wherever mu < ||v'||, and sv = ||v'|| exactly")
print("in the three cells where mu > ||v'||.  That is the law")
print("      lambda_1(L') = min( mu, ||v'_planted|| ) ,   mu independent of m.")
print("Tested here on fresh curves not used in EXP A-F.\n")

LO, HI = 2 ** 14, 2 ** 15
print(f"Searching fresh j=0 GLV curves with p in [{LO}, {HI}) ...")
fresh = search_curves(LO, HI, per_bin=1, nbins=8)
print(f"Found {len(fresh)} curves.\n")

G_SEEDS = [42, 1234, 9999]
G_K1 = [4, 12, 32]
hits = total = 0
viol = []
print(f"  {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} {'nu_hat':>7} {'sv_obs':>12} "
      f"{'min(mu,pv)':>12} {'law':>5} {'mu/||v'+chr(39)+'||':>11} {'V1':>5}")
by_nu = {'lo': [0, 0], 'hi': [0, 0]}
for (p, b, n, lam, G) in fresh:
    k2 = math.isqrt(n) + 1
    for k1 in G_K1:
        agg = sweep((p, b, n, lam, G), M_SIGS, k1, G_SEEDS)
        if agg['tot'] == 0:
            continue
        S1, _s, S2, _k = scales(n, k1, k2)
        mu, _ = lambda_block_mu(n, lam, S1, S2)
        nh = mu / math.sqrt(n * S1 * S2)
        sv, pv = mean(agg['sv1']), mean(agg['pv1'])
        pred = min(mu, pv)
        ok = abs(sv - pred) / pred < 1e-9
        hits += 1 if ok else 0
        total += 1
        if not ok:
            # 'shorter' = a cross-block combination beats both mu and ||v'||;
            # 'LLLmiss' = LLL simply did not return the known short vector.
            viol.append((n, k1, mu / pv,
                         'shorter' if sv < pred else 'LLLmiss'))
        bucket = 'lo' if nh < 0.55 else 'hi'
        by_nu[bucket][0] += agg['v1']
        by_nu[bucket][1] += agg['tot']
        print(f"  {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {k1*k2/n:>6.3f} "
              f"{nh:>7.4f} {sv:>12.1f} {pred:>12.1f} {str(ok):>5} "
              f"{mu/pv:>11.3f} {str(agg['v1'])+'/'+str(agg['tot']):>5}")
print(f"\n  law lambda_1(L') = min(mu, ||v'||) holds in {hits}/{total} cells")
nsh = sum(1 for v in viol if v[3] == 'shorter')
print(f"  violations: {nsh} 'shorter vector exists', {len(viol)-nsh} 'LLL miss'")
if viol:
    rr = [v[2] for v in viol]
    print(f"  all violations sit in the crossover band mu/||v'|| in "
          f"[{min(rr):.2f}, {max(rr):.2f}]")
print("  NOTE: mu is a lattice vector by construction, so lambda_1(L') <= mu")
print("  holds unconditionally; only the equality is approximate.  The")
print("  m-independence conclusion rests on the inequality, not the equality.")
for k in ('lo', 'hi'):
    w, t = by_nu[k]
    lbl = 'nu_hat < 0.55' if k == 'lo' else 'nu_hat >= 0.55'
    print(f"  {lbl:<16} recovery {w}/{t}"
          + (f" = {w/t:.3f}" if t else ""))

print("\nDone.")
