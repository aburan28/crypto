"""
GLV-HNP Phase 2 — Thread 23: reformulate the lattice so the target is recoverable
by CVP instead of SVP.

Background (2026-07-29, RESEARCH_AUTOLAB_LOG.md T5)
--------------------------------------------------
The Phase-2 lattice built by `glv_hnp_phase2_lambda_threshold.py:254`
(`build_glv_lattice`) is a (2m+2)-dimensional Kannan embedding with a d-column
scaled by S_D = 1.  It contains the vector `n * S_D * e_m` (multiply the d-row by
n and subtract n*B_i times each modular row), which is 2-3x SHORTER than the
planted vector on every curve tested and carries no information at all (d is only
defined mod n).  So the planted vector is never lambda_1, and recovery is a
BDD/coset condition, not SVP.  T5 measured |sv[m]|/n = 1.0000 exactly on 3/3
historical curves and sv/pv in [0.337, 0.368] across a 17-bit sweep.

This script implements the proposed fix and runs the falsifier stated in the
2026-07-29 "Next step proposal".

The reformulation
-----------------
Drop BOTH the Kannan row and the d-column.  Work in dimension 2m with

    rows   r_i = n*S1*e_i                          (i = 0..m-1)  [the c_i freedom]
           g_i = -lam*S1*e_i + S2*e_{m+i}          (i = 0..m-1)  [the k2_i freedom]
           b   = sum_i B_i*S1*e_i                                [the d freedom]

    target t   = (-A_i*S1)_i  ||  (0)_i

A lattice vector v = d*b + sum k2_i*g_i - sum c_i*r_i has

    v_i     - t_i     = S1*(A_i + B_i*d - lam*k2_i - c_i*n) = S1*k1_i
    v_{m+i} - t_{m+i} = S2*k2_i

so the planted point is the lattice point closest to t, with error norm
~ n*sqrt(2m/3) under the standard scaling S1 = n/K1, S2 = n/K2.  There is no
trivial short vector: `n*e_m` does not exist because there is no d-column, and d
is read back out of the recovered (k1_i, k2_i) via d = (k1_0 + lam*k2_0 - A_0)/B_0
mod n.  The dimension also drops from 2m+2 to 2m.

Three solvers are compared on the SAME signatures:
  KAN   = the 2026-07-29 baseline (Kannan embedding + LLL, read d off the row
          whose last coordinate is +-S_KANNAN)
  BAB   = this lattice + LLL + Babai nearest plane
  CVPx  = this lattice + LLL + exact CVP (fpylll enumeration)

CVPx is the decisive one: if exact CVP does not return the planted point, the
planted point is genuinely not the closest lattice vector and the wall is
information-theoretic, not an artifact of the reduction.

Falsifier (from the 2026-07-29 log entry)
-----------------------------------------
"if the K1 wall in T4 moves outward on the lam*=0.07 curve (currently K1~4-6),
the reformulation is a real improvement; if the wall stays at K1~4-6, then the
wall is information-theoretic and Phase 2 is at its ceiling."

Run: python3 glv_hnp_phase2_cvp.py
"""

import importlib.util
import math
import os
import random
import sys
import time

import sympy
from fpylll import CVP, GSO, IntegerMatrix, LLL

_HERE = os.path.dirname(os.path.abspath(__file__))

# glv_hnp_phase2_lambda_threshold.py runs its own experiments at module level
# (from the "EXPERIMENT T1" banner onwards), so import only the helper prefix.
# Everything below is therefore the 2026-07-29 code verbatim.
_SRC = open(os.path.join(_HERE, "glv_hnp_phase2_lambda_threshold.py")).read()
_CUT = _SRC.index("# EXPERIMENT T1")
_t20a = importlib.util.module_from_spec(
    importlib.util.spec_from_loader("_t20a", loader=None))
exec(compile(_SRC[:_CUT], "glv_hnp_phase2_lambda_threshold.py(prefix)", "exec"),
     _t20a.__dict__)

modinv = _t20a.modinv
scales = _t20a.scales
gen_signatures = _t20a.gen_signatures
build_glv_lattice = _t20a.build_glv_lattice
planted_vector = _t20a.planted_vector
recover_d = _t20a.recover_d
norm = _t20a.norm
lam_star = _t20a.lam_star
build_curve = _t20a.build_curve
search_curves = _t20a.search_curves
find_generator = _t20a.find_generator

SEEDS = [42, 1234, 9999, 555, 31337]


# ---------------------------------------------------------------------------
# The reformulated (2m-dimensional) lattice
# ---------------------------------------------------------------------------

def build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound):
    """Generators (2m+1 rows, rank 2m) and the CVP target."""
    m = len(sigs)
    S1, _S_D, S2, _S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                       # c_i freedom
        r = [0] * (2 * m)
        r[i] = n * S1
        rows.append(r)
    for i in range(m):                       # k2_i freedom
        r = [0] * (2 * m)
        r[i] = -lam * S1
        r[m + i] = S2
        rows.append(r)
    b = [0] * (2 * m)                        # d freedom
    for i in range(m):
        b[i] = sigs[i]['B'] * S1
    rows.append(b)
    target = [(-sigs[i]['A'] * S1) % (n * S1) for i in range(m)] + [0] * m
    return rows, target, S1, S2


def planted_error(sigs, n, k1_bound, k2_bound):
    """The error vector v - t of the planted point: (S1*k1_i, S2*k2_i)."""
    m = len(sigs)
    S1, _, S2, _ = scales(n, k1_bound, k2_bound)
    e = [0] * (2 * m)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S1
        e[m + i] = sigs[i]['k2'] * S2
    return e


def reduced_basis(rows):
    """LLL-reduce a generating set; return a full-rank basis (zero rows stripped)."""
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    keep = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)
            if any(A[i][j] for j in range(A.ncols))]
    return IntegerMatrix.from_matrix(keep)


def d_from_error(err, sigs, n, lam, k1_bound, k2_bound):
    """Recover d from a candidate error vector (S1*k1_i, S2*k2_i)."""
    m = len(sigs)
    S1, _, S2, _ = scales(n, k1_bound, k2_bound)
    if any(err[i] % S1 for i in range(m)) or any(err[m + i] % S2 for i in range(m)):
        return None
    k1 = [err[i] // S1 for i in range(m)]
    k2 = [err[m + i] // S2 for i in range(m)]
    for i in range(m):
        B = sigs[i]['B'] % n
        if math.gcd(B, n) != 1:
            continue
        k = (k1[i] + lam * k2[i]) % n
        d = (k - sigs[i]['A']) * modinv(B, n) % n
        # cross-check against every signature
        ok = True
        for j in range(m):
            kj = (sigs[j]['A'] + sigs[j]['B'] * d) % n
            if kj != (k1[j] + lam * k2[j]) % n:
                ok = False
                break
        if ok:
            return d
    return None


def box_center(m, k1_bound, k2_bound, n):
    """Centre of the error box, as an exact lattice-coordinate-compatible shift.

    The unknowns satisfy 0 <= k1_i < K1 and 0 <= k2_i < K2, so the error vector
    v - t has non-negative entries with mean (S1*K1/2, S2*K2/2).  Solving CVP
    against t therefore aims at a corner of the box rather than its centre, and
    the planted error is ~2x longer than it needs to be:
        uncentred  E|e|^2 = 2m*n^2/3      (E[u^2] = 1/3 for u ~ U[0,1])
        centred    E|e|^2 = 2m*n^2/12     (E[u^2] = 1/12 for u ~ U[-1/2,1/2])
    The shift is taken as an exact multiple of S1 (resp. S2) so that the
    divisibility test in d_from_error still holds.
    """
    S1, _, S2, _ = scales(n, k1_bound, k2_bound)
    return [S1 * ((k1_bound - 1) // 2)] * m + [S2 * ((k2_bound - 1) // 2)] * m


def centered_sigs(sigs, n, lam, k1_bound, k2_bound):
    """Re-centre the unknowns by rewriting the PUBLIC data.

    k1_i = A_i + B_i*d - lam*k2_i - c_i*n with k1_i in [0,K1), k2_i in [0,K2).
    Put h1 = (K1-1)//2, h2 = (K2-1)//2 and A'_i = A_i - h1 - lam*h2 mod n.  Then
    k1'_i = k1_i - h1 and k2'_i = k2_i - h2 satisfy k1' + lam*k2' = A' + B*d with
    the SAME d, and now range over [-h1, h1] / [-h2, h2].  So centring is a
    rewrite of the instance, not a change of lattice, and it can be applied to
    any solver -- including the original Kannan embedding.  This is the control
    that separates "CVP reformulation" from "centring".
    """
    h1 = (k1_bound - 1) // 2
    h2 = (k2_bound - 1) // 2
    return [{'A': (s['A'] - h1 - lam * h2) % n, 'B': s['B'],
             'k1': s['k1'] - h1, 'k2': s['k2'] - h2} for s in sigs]


def solve_cvp(sigs, n, lam, k1_bound, k2_bound, exact=False, center=False):
    """Returns (d_candidate, err_norm_vs_true_target, lambda1_of_lattice)."""
    m = len(sigs)
    rows, target, S1, S2 = build_cvp_lattice(sigs, n, lam, k1_bound, k2_bound)
    B = reduced_basis(rows)
    dim = B.ncols
    lam1 = min(norm([B[i][j] for j in range(dim)]) for i in range(B.nrows))
    if center:
        c = box_center(m, k1_bound, k2_bound, n)
        aim = [target[j] + c[j] for j in range(dim)]
    else:
        aim = list(target)
    if exact:
        v = list(CVP.closest_vector(B, tuple(aim)))
    else:
        M = GSO.Mat(B)
        M.update_gso()
        coords = M.babai(list(aim))
        v = [0] * dim
        for i, co in enumerate(coords):
            if co:
                for j in range(dim):
                    v[j] += int(co) * B[i][j]
    # d is always read off the error against the TRUE target
    err = [v[j] - target[j] for j in range(dim)]
    d = d_from_error(err, sigs, n, lam, k1_bound, k2_bound)
    return d, norm([v[j] - aim[j] for j in range(dim)]), lam1


# ---------------------------------------------------------------------------
# One trial: all three solvers on identical signatures
# ---------------------------------------------------------------------------

def trial(curve, m, d_secret, k1_bound, seed, do_exact=True):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S1, S_D, S2, S_KAN = scales(n, k1_bound, k2_bound)

    # --- KAN baseline (2026-07-29 construction, verbatim) ---
    M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dimK = 2 * m + 2
    red = [[A[i][j] for j in range(dimK)] for i in range(dimK)]
    kan_ok = recover_d(red, m, n, S_KAN, d_secret) is not None
    pv = norm(planted_vector(sigs, d_secret, n, k1_bound, k2_bound))
    sv = min(norm(r) for r in red)

    # --- KANc control: same Kannan embedding, centred instance ---
    csigs = centered_sigs(sigs, n, lam, k1_bound, k2_bound)
    Mc = build_glv_lattice(csigs, n, lam, k1_bound, k2_bound)
    Ac = IntegerMatrix.from_matrix(Mc)
    LLL.reduction(Ac)
    redc = [[Ac[i][j] for j in range(dimK)] for i in range(dimK)]
    kanc_ok = recover_d(redc, m, n, S_KAN, d_secret) is not None

    # --- reformulated lattice ---
    pe_vec = planted_error(sigs, n, k1_bound, k2_bound)
    pe = norm(pe_vec)
    ctr = box_center(m, k1_bound, k2_bound, n)
    pe_c = norm([pe_vec[j] - ctr[j] for j in range(2 * m)])
    d_bab, err_bab, lam1 = solve_cvp(sigs, n, lam, k1_bound, k2_bound, exact=False)
    bab_ok = (d_bab == d_secret)
    d_babc, err_babc, _ = solve_cvp(sigs, n, lam, k1_bound, k2_bound,
                                    exact=False, center=True)
    babc_ok = (d_babc == d_secret)
    if do_exact:
        d_cvp, err_cvp, _ = solve_cvp(sigs, n, lam, k1_bound, k2_bound, exact=True)
        cvp_ok = (d_cvp == d_secret)
        d_cvpc, err_cvpc, _ = solve_cvp(sigs, n, lam, k1_bound, k2_bound,
                                        exact=True, center=True)
        cvpc_ok = (d_cvpc == d_secret)
    else:
        cvp_ok, err_cvp = None, float('nan')
        cvpc_ok, err_cvpc = None, float('nan')

    # Classification of the exact-CVP outcome.
    #   'planted'  : the CVP solution IS the planted point (err == pe) -> d found
    #   'alt'      : a strictly closer point that still yields d (the GLV
    #                decomposition k = k1 + lam*k2 is not unique, so a shorter
    #                (k1,k2) for the same k is an equally valid solution)
    #   'decoy'    : a strictly closer point that does NOT yield d.  This is
    #                information-theoretic failure: no CVP oracle can help,
    #                because the true closest vector is not the secret.
    def classify(ok, err, planted_err):
        if ok is None:
            return None
        if ok:
            return 'planted' if err >= planted_err - 1e-9 else 'alt'
        return 'decoy' if err < planted_err - 1e-9 else 'miss'

    cvp_kind = classify(cvp_ok, err_cvp, pe)
    cvpc_kind = classify(cvpc_ok, err_cvpc, pe_c)

    return {
        'kan': kan_ok, 'kanc': kanc_ok, 'bab': bab_ok, 'cvp': cvp_ok,
        'babc': babc_ok, 'cvpc': cvpc_ok,
        'cvp_kind': cvp_kind, 'cvpc_kind': cvpc_kind,
        'sv_pv': sv / pv if pv else float('nan'),
        'pe': pe, 'pe_c': pe_c,
        'lam1': lam1, 'pe_lam1': pe / lam1 if lam1 else float('nan'),
        'pec_lam1': pe_c / lam1 if lam1 else float('nan'),
        'err_bab': err_bab, 'err_cvp': err_cvp, 'err_cvpc': err_cvpc,
        'gap': err_cvp / pe if (do_exact and pe) else float('nan'),
        'gap_c': err_cvpc / pe_c if (do_exact and pe_c) else float('nan'),
    }


def rates(curve, m, k1_bound, seeds, do_exact=True):
    agg = {'kan': 0, 'kanc': 0, 'bab': 0, 'cvp': 0, 'babc': 0, 'cvpc': 0, 'n': 0,
           'planted': 0, 'alt': 0, 'decoy': 0, 'miss': 0,
           'c_planted': 0, 'c_alt': 0, 'c_decoy': 0, 'c_miss': 0}
    pe_lam1, pec_lam1, svpv, gaps, gaps_c = [], [], [], [], []
    for s in seeds:
        d_trial = random.Random(s + 7777).randint(1, curve[2] - 1)
        r = trial(curve, m, d_trial, k1_bound, s, do_exact=do_exact)
        if r is None:
            continue
        agg['n'] += 1
        agg['kan'] += bool(r['kan'])
        agg['kanc'] += bool(r['kanc'])
        agg['bab'] += bool(r['bab'])
        agg['babc'] += bool(r['babc'])
        if do_exact:
            agg['cvp'] += bool(r['cvp'])
            agg['cvpc'] += bool(r['cvpc'])
            agg[r['cvp_kind']] += 1
            agg['c_' + r['cvpc_kind']] += 1
            gaps.append(r['gap'])
            gaps_c.append(r['gap_c'])
        pe_lam1.append(r['pe_lam1'])
        pec_lam1.append(r['pec_lam1'])
        svpv.append(r['sv_pv'])
    mean = lambda xs: sum(xs) / len(xs) if xs else float('nan')
    agg['pe_lam1'] = mean(pe_lam1)
    agg['pec_lam1'] = mean(pec_lam1)
    agg['sv_pv'] = mean(svpv)
    agg['gap'] = mean(gaps)
    agg['gap_c'] = mean(gaps_c)
    return agg


# ===========================================================================
print("=" * 78)
print("Thread 23 — CVP reformulation of the GLV-HNP Phase 2 lattice")
print("=" * 78)
t_start = time.time()

# Historical curves, identical to glv_hnp_phase2_lambda_threshold.py:423
HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
    ("8-bit/199",   211,  2, 199,  106),
]
hist = []
for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    hist.append((label, (p, b, n, lam, G)))
    print(f"  {label}: p={p} n={n} lam={lam} lam*={lam_star(lam, n):.4f}")


# ===========================================================================
# E0 — sanity: the reformulated lattice has no trivial short vector
# ===========================================================================
print("\n" + "=" * 78)
print("E0 — geometry check: planted error vs lambda_1 (m=12, K1=8)")
print("=" * 78)
print(f"{'curve':>12} {'|pe|':>12} {'lam1(L)':>12} {'|pe|/lam1':>10} "
      f"{'sv/pv (KAN)':>12}")
for label, curve in hist:
    n = curve[2]
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    r = trial(curve, 12, d_trial, 8, 42)
    print(f"{label:>12} {r['pe']:>12.4g} {r['lam1']:>12.4g} "
          f"{r['pe_lam1']:>10.3f} {r['sv_pv']:>12.3f}")
print("\n(KAN sv/pv < 1 means the planted vector was never lambda_1 — the T5")
print(" finding.  |pe|/lam1 > 1 in the new lattice is expected and harmless:")
print(" BDD does not require the error to be shorter than lambda_1, only that")
print(" it lie inside the Voronoi cell / within lambda_1/2 for uniqueness.)")


# ===========================================================================
# E0b — the centring rewrite is exact (assertion, not a measurement)
# ===========================================================================
print("\n" + "=" * 78)
print("E0b — centring rewrite faithfulness check")
print("=" * 78)
_checked = 0
for _label, _curve in hist:
    _p, _b, _n, _lam, _G = _curve
    _k2 = math.isqrt(_n) + 1
    for _k1 in (4, 8, 24):
        _d = random.Random(999).randint(1, _n - 1)
        _sg = gen_signatures(_G, _d, 8, _n, _lam, _p, _k1, _k2, 42)
        _cs = centered_sigs(_sg, _n, _lam, _k1, _k2)
        _h1, _h2 = (_k1 - 1) // 2, (_k2 - 1) // 2
        for _o, _c in zip(_sg, _cs):
            assert (_c['k1'] + _lam * _c['k2']) % _n == (_c['A'] + _c['B'] * _d) % _n
            assert _c['k1'] == _o['k1'] - _h1 and _c['k2'] == _o['k2'] - _h2
            assert -_h1 <= _c['k1'] < _k1 - _h1
            _checked += 1
print(f"  OK — {_checked} signatures: k1' + lam*k2' == A' + B*d (mod n) with the")
print("  SAME d, and k1' in [-h1, K1-1-h1].  Centring is a rewrite of the public")
print("  data only; it uses no secret information.")


# ===========================================================================
# E1 — the T4 K1 grid: does the wall move?
# ===========================================================================
print("\n" + "=" * 78)
print("E1 — K1 wall (T4 replication + reformulation), m=12, 5 seeds")
print("=" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
M_E1 = 12
e1 = {}
for label, curve in hist[:2]:
    n = curve[2]
    k2 = math.isqrt(n) + 1
    print(f"\n--- {label}  (lam*={lam_star(curve[3], n):.3f}, K2={k2}) ---")
    print(f"{'K1':>4} {'eff':>7} {'KAN':>6} {'BAB':>6} {'CVPx':>6} "
          f"{'KANc':>6} {'BABc':>6} {'CVPc':>6} {'gap':>6} {'gap_c':>6}  "
          f"{'pl/alt/dec/mis':>14}  {'centred kinds':>14}")
    e1[label] = []
    for k1 in K1_GRID:
        a = rates(curve, M_E1, k1, SEEDS)
        eff = k1 * k2 / n
        kinds = f"{a['planted']}/{a['alt']}/{a['decoy']}/{a['miss']}"
        kinds_c = (f"{a['c_planted']}/{a['c_alt']}/{a['c_decoy']}/"
                   f"{a['c_miss']}")
        print(f"{k1:>4} {eff:>7.3f} {str(a['kan'])+'/'+str(a['n']):>6} "
              f"{str(a['bab'])+'/'+str(a['n']):>6} "
              f"{str(a['cvp'])+'/'+str(a['n']):>6} "
              f"{str(a['kanc'])+'/'+str(a['n']):>6} "
              f"{str(a['babc'])+'/'+str(a['n']):>6} "
              f"{str(a['cvpc'])+'/'+str(a['n']):>6} "
              f"{a['gap']:>6.3f} {a['gap_c']:>6.3f}  {kinds:>14}  {kinds_c:>14}")
        e1[label].append((k1, eff, a))

print("\nWall = largest K1 with a full 5/5 sweep:")
for label, rows in e1.items():
    walls = {}
    for key in ('kan', 'bab', 'cvp', 'kanc', 'babc', 'cvpc'):
        w = None
        for k1, eff, a in rows:
            if a[key] == a['n']:
                w = k1
        walls[key] = w
    print(f"  {label}: KAN K1<={walls['kan']}  BAB K1<={walls['bab']}  "
          f"CVPx K1<={walls['cvp']}  |  KANc K1<={walls['kanc']}  "
          f"BABc K1<={walls['babc']}  CVPc K1<={walls['cvpc']}")


# ===========================================================================
# E2 — m-scaling at a fixed hard K1
# ===========================================================================
print("\n" + "=" * 78)
print("E2 — more signatures at the KAN wall (K1=8, the 2026-07-26 setting)")
print("=" * 78)
print(f"{'curve':>12} {'m':>4} {'KAN':>6} {'BAB':>6} {'CVPx':>6} "
      f"{'KANc':>6} {'BABc':>6} {'CVPc':>6}")
for label, curve in hist[:2]:
    for m in (8, 12, 16, 24, 32):
        a = rates(curve, m, 8, SEEDS, do_exact=(m <= 16))
        cvp_s = f"{a['cvp']}/{a['n']}" if m <= 16 else "-"
        cvpc_s = f"{a['cvpc']}/{a['n']}" if m <= 16 else "-"
        print(f"{label:>12} {m:>4} {str(a['kan'])+'/'+str(a['n']):>6} "
              f"{str(a['bab'])+'/'+str(a['n']):>6} {cvp_s:>6} "
              f"{str(a['kanc'])+'/'+str(a['n']):>6} "
              f"{str(a['babc'])+'/'+str(a['n']):>6} {cvpc_s:>6}")


# ===========================================================================
# E3 — fresh 17-bit curves: does the reformulation extend the usable eff?
# ===========================================================================
print("\n" + "=" * 78)
print("E3 — fresh 17-bit j=0 GLV curves, eff sweep (m=12, 5 seeds)")
print("=" * 78)
LO, HI = 2 ** 16, 2 ** 17
print(f"Searching j=0 GLV curves with p in [{LO}, {HI}) ...")
curves = search_curves(LO, HI, per_bin=1, nbins=10)
print(f"Found {len(curves)} curves.")

for eff_t in (0.05, 0.15, 0.25, 0.40):
    print(f"\n=== eff target = {eff_t} ===")
    print(f"{'p':>7} {'n':>7} {'lam*':>6} {'K1':>5} {'eff':>6} "
          f"{'KAN':>6} {'BAB':>6} {'CVPx':>6} {'KANc':>6} {'BABc':>6} "
          f"{'CVPc':>6} {'gap_c':>6}  {'centred kinds':>14}")
    tot = {'kan': 0, 'kanc': 0, 'bab': 0, 'cvp': 0, 'babc': 0, 'cvpc': 0, 'n': 0,
           'planted': 0, 'alt': 0, 'decoy': 0, 'miss': 0,
           'c_planted': 0, 'c_alt': 0, 'c_decoy': 0, 'c_miss': 0}
    for (p, b, n, lam, G) in curves:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(eff_t * n / k2))
        a = rates((p, b, n, lam, G), 12, k1, SEEDS)
        for key in tot:
            tot[key] += a[key]
        kinds_c = f"{a['c_planted']}/{a['c_alt']}/{a['c_decoy']}/{a['c_miss']}"
        print(f"{p:>7} {n:>7} {lam_star(lam, n):>6.3f} {k1:>5} "
              f"{k1 * k2 / n:>6.3f} {str(a['kan'])+'/'+str(a['n']):>6} "
              f"{str(a['bab'])+'/'+str(a['n']):>6} "
              f"{str(a['cvp'])+'/'+str(a['n']):>6} "
              f"{str(a['kanc'])+'/'+str(a['n']):>6} "
              f"{str(a['babc'])+'/'+str(a['n']):>6} "
              f"{str(a['cvpc'])+'/'+str(a['n']):>6} "
              f"{a['gap_c']:>6.3f}  {kinds_c:>14}")
    print(f"  TOTAL eff={eff_t}: KAN {tot['kan']}/{tot['n']}  "
          f"BAB {tot['bab']}/{tot['n']}  CVPx {tot['cvp']}/{tot['n']}  |  "
          f"KANc {tot['kanc']}/{tot['n']}  BABc {tot['babc']}/{tot['n']}  "
          f"CVPc {tot['cvpc']}/{tot['n']}")
    print(f"    uncentred kinds: planted {tot['planted']}, alt {tot['alt']}, "
          f"decoy {tot['decoy']}, miss {tot['miss']}")
    print(f"    centred   kinds: planted {tot['c_planted']}, alt {tot['c_alt']}, "
          f"decoy {tot['c_decoy']}, miss {tot['c_miss']}")


# ===========================================================================
# E4 — does the centred wall move outward with more signatures?
# ===========================================================================
print("\n" + "=" * 78)
print("E4 — eff wall vs number of signatures m (centred solvers, 5 seeds)")
print("=" * 78)
EFF_GRID = [0.05, 0.15, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60]
E4_CURVES = curves[:5]
print(f"{'m':>4} " + " ".join(f"{e:>11.2f}" for e in EFF_GRID))
print("     " + " ".join(f"{'KANc/CVPc':>11}" for _ in EFF_GRID))
for m in (6, 8, 12, 16, 24, 32):
    cells = []
    for eff_t in EFF_GRID:
        kc = cc = tot = 0
        for (p, b, n, lam, G) in E4_CURVES:
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            a = rates((p, b, n, lam, G), m, k1, SEEDS, do_exact=(m <= 16))
            kc += a['kanc']
            cc += a['cvpc']
            tot += a['n']
        cc_s = str(cc) if m <= 16 else "-"
        cells.append(f"{kc}/{cc_s} of {tot}")
    print(f"{m:>4} " + " ".join(f"{c:>11}" for c in cells))
print("\n(cell = KANc / CVPc successes out of 25 trials; '-' = exact CVP skipped)")

# ===========================================================================
# E5 — does the same defect sit in the plain (non-GLV) BV lattice?
# ===========================================================================
# src/cryptanalysis/hnp_ecdsa.rs:207-229 builds the classic Boneh-Venkatesan
# basis
#       n^2 * e_i                       (i < m)
#       (n*a_i , 2^l)                   row m
#       (n*t_i , 0.., n*2^l)            row m+1
# and reads d off row[m].  Its planted vector has entries n*k_i with
# k_i in [0, 2^l) -- uncentred in exactly the same way as the GLV lattice.
# The fix is the same one-line rewrite: t_i <- t_i - 2^(l-1) mod n.
print("\n" + "=" * 78)
print("E5 — plain BV lattice (the Rust code path), uncentred vs centred")
print("=" * 78)

def bv_trial(n, m, l, seed, center):
    """Classic BV HNP: k_i = t_i + a_i*d mod n with 0 <= k_i < 2^l."""
    rng = random.Random(seed)
    d = rng.randrange(1, n)
    two_l = 1 << l
    a, t, ks = [], [], []
    for _ in range(m):
        ai = rng.randrange(1, n)
        ki = rng.randrange(0, two_l)
        ti = (ki - ai * d) % n
        a.append(ai); t.append(ti); ks.append(ki)
    if center:
        t = [(ti - (two_l >> 1)) % n for ti in t]
    dim = m + 2
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * n
    for i in range(m):
        M[m][i] = n * a[i]
    M[m][m] = two_l
    for i in range(m):
        M[m + 1][i] = n * t[i]
    M[m + 1][m + 1] = n * two_l
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    inv = modinv(two_l % n, n)
    for i in range(dim):
        for sign in (1, -1):
            val = (sign * A[i][m]) % n
            if val == 0:
                continue
            if (val * inv) % n == d:
                return True
    return False

M_BV = 40
N_BV = int(sympy.nextprime(2 ** 40)) if False else 1099511627791  # 41-bit prime
print(f"  n = {N_BV} ({N_BV.bit_length()} bits), m = {M_BV} signatures, 10 seeds")
print(f"{'l (nonce bits)':>15} {'bias bits':>10} {'uncentred':>10} {'centred':>10}")
BV_SEEDS = list(range(10))
for l in range(N_BV.bit_length() - 1, N_BV.bit_length() - 9, -1):
    u = sum(bv_trial(N_BV, M_BV, l, s, False) for s in BV_SEEDS)
    c = sum(bv_trial(N_BV, M_BV, l, s, True) for s in BV_SEEDS)
    print(f"{l:>15} {N_BV.bit_length() - l:>10} "
          f"{str(u) + '/' + str(len(BV_SEEDS)):>10} "
          f"{str(c) + '/' + str(len(BV_SEEDS)):>10}")

print(f"\nElapsed: {time.time() - t_start:.1f}s")
print("Done.")
