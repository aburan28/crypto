"""
GLV-HNP Thread 23 — is the Phase-2 wall a uSVP-geometry wall?

Background
----------
Thread 20 (2026-07-29, commits d525931 / e845207) established:
  T4   the K1 wall is real and is not a lam wall (lam* only shifts it ~3x);
  T4b  at K1=8 the lam*=0.070 curve is not rescued by more signatures
       (m = 8/12/16/24/32 -> 0,0,1,0,1 of 5);
  T5   the planted vector is NEVER lambda_1: the trivial vector n*S_D*e_m is
       always shorter, so recovery is a BDD/coset condition, not SVP;
  20c  nu_hat = lambda_1(L2)/sqrt(det L2) of the 2D rival block correlates with
       success, with the *counter-intuitive* sign (small nu_hat = easier).
and proposed Thread 23: "project the lattice along e_m (quotient out the
trivial n*e_m direction) ... Falsifier: if sv/pv rises above 1 after the
reformulation and the K1 wall moves outward, the reformulation is a real
improvement; if the wall stays at K1 ~ 4-6, the wall is information-theoretic
and Phase 2 is at its ceiling."

This script executes that falsifier and adds the missing quantitative frame.

Hypotheses under test
---------------------
H23a  Removing the d column (kill the trivial vector, dim 2m+2 -> 2m+1) moves
      the K1 wall outward.
      Falsifier: identical wall position on the T4 grid.

H23b  The wall is a *Gaussian-heuristic* wall, i.e. the planted vector is
      simply longer than GH(L), so no amount of reduction can help.  The
      closed form derived here is

          nu_gh = ||v||/GH(L) ~= sqrt(2*pi*e/3) * (n * eff^m)^(1/(2m+2))

      with eff = K1*K2/n.  Note lim_{m->inf} nu_gh = sqrt(2*pi*e/3)*sqrt(eff),
      so nu_gh < 1 is achievable at large m only when eff < 3/(2*pi*e) = 0.1757.
      Falsifier: nu_gh does not separate success from failure, or the observed
      successes sit at nu_gh > 1.

H23c  Residual (lam-dependent) variation not explained by nu_gh is explained by
      the local Gram-Schmidt profile via the ADPS visibility margin
          margin = max_i ||b*_i|| / ||pi_i(v)||,
      which should beat nu_gh, nu_hat and lam* on AUC.
      Falsifier: margin's AUC is not above nu_gh's.

Run: python3 glv_hnp_usvp_gap.py
"""

import math
import random
import sys
import time

import glv_hnp_lib as L
from fpylll import IntegerMatrix, LLL

t_start = time.time()

# (label, p, b, n, lam, K1, m) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py:425-427
ANCHORS = [
    ("8-bit/199",   211,  2, 199,  106,  2,  6),
    ("12-bit/2557", 2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677", 2677, 2, 2647, 185,  8,  10),
]

SEEDS = [42, 1234, 9999, 555, 31337]


def anchor_curves():
    out = []
    for label, p, b, n, lam, k1, m in ANCHORS:
        G = L.find_generator(p, b, n)
        assert G is not None, label
        out.append((label, (p, b, n, lam, G), k1, m))
    return out


def instance_geometry(curve, m, d_secret, k1_bound, seed):
    """Everything measurable about one (curve, m, K1, seed) instance."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = L.gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = L.scales(n, k1_bound, k2_bound)
    dim = 2 * m + 2

    M = L.build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
    v = L.planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    pn = L.norm(v)

    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    reduced = [[A[i][j] for j in range(dim)] for i in range(dim)]
    ok = L.recover_d(reduced, m, n, S_KAN, d_secret) is not None
    sv = min(L.norm(r) for r in reduced)

    margin, mi = L.adps_margin(reduced, v)
    _w, mu, nu_hat = L.rival_sublattice_nu(n, lam, k1_bound, k2_bound)
    nu_proj, proj_nrm, _mu = L.usvp_gap_projected(
        n, lam, m, k1_bound, k2_bound, sigs, d_secret)

    return {
        'n': n, 'lam': lam, 'm': m, 'K1': k1_bound, 'K2': k2_bound,
        'eff': k1_bound * k2_bound / n,
        'ok': ok,
        'pn': pn, 'sv': sv, 'sv_over_pn': sv / pn,
        'nu_gh': L.usvp_gap(n, m, k1_bound, k2_bound, pn),
        'nu_cf': L.usvp_gap_closed_form(n, m, k1_bound, k2_bound),
        'nu_hat': nu_hat,
        'nu_proj': nu_proj, 'proj_nrm': proj_nrm,
        'lam_star': L.lam_star(lam, n),
        'margin': margin, 'margin_i': mi,
    }


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ===========================================================================
# E1 -- exact geometry of the anchors: is the planted vector above GH?
# ===========================================================================
banner("E1  Exact uSVP geometry of the three anchor curves")

print("nu_gh = ||v_planted|| / GH(L).  nu_gh > 1 means the Gaussian heuristic "
      "predicts\n~nu_gh^dim lattice vectors SHORTER than the target, i.e. uSVP "
      "is ill-posed.\n")
print(f"{'curve':<14}{'m':>3}{'K1':>4}{'eff':>7}{'lam*':>7}{'nu_hat':>8}"
      f"{'nu_gh':>8}{'nu_cf':>8}{'nu_proj':>9}{'sv/pv':>8}{'margin':>8}"
      f"{'i*':>4}{'LLL':>6}")

anchors = anchor_curves()
e1_rows = []
for label, curve, k1, m in anchors:
    wins = 0
    agg = None
    for seed in SEEDS:
        d = random.Random(seed + 7777).randint(1, curve[2] - 1)
        g = instance_geometry(curve, m, d, k1, seed)
        if g is None:
            continue
        wins += g['ok']
        agg = g if agg is None else agg
        e1_rows.append(g)
    g = agg
    print(f"{label:<14}{m:>3}{k1:>4}{g['eff']:>7.3f}{g['lam_star']:>7.3f}"
          f"{g['nu_hat']:>8.3f}{g['nu_gh']:>8.3f}{g['nu_cf']:>8.3f}"
          f"{g['nu_proj']:>9.3f}{g['sv_over_pn']:>8.3f}{g['margin']:>8.3f}"
          f"{g['margin_i']:>4}{str(wins) + '/' + str(len(SEEDS)):>6}")

print("\nClosed-form check: nu_cf is computed from (n, m, K1, K2) alone; nu_gh "
      "uses the\nexact integer det and the exact planted norm of seed 42.")
print("nu_proj is the gap in the projection orthogonal to the rank-(m+1)")
print("nuisance sublattice T = <n*S_D*e_m> + <rival block shortest vectors>.")
print("T is primitive in L (index [sat(T):T] = 1, verified exactly by HNF on")
print("the 8-bit and 12-bit/2557 anchors), so det(pi(L)) = det(L)/det(T) holds")
print("with equality, not as a heuristic.")

# ===========================================================================
# E2 -- H23a: does deleting the d column move the K1 wall?
# ===========================================================================
banner("E2  H23a -- Thread 23's proposed reformulation (delete the d column)")

print("Left block: the (2m+2)-dim Phase-2 lattice (Thread 20 baseline).")
print("Right block: the (2m+1)-dim lattice with column m deleted; the trivial")
print("vector n*S_D*e_m maps to 0 and d is recovered post hoc from")
print("d = B_0^{-1}(k1_0 - A_0 + lam*k2_0) mod n.\n")

K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
print(f"{'curve':<14}{'lat':>6}" + "".join(f"{('K1=' + str(k)):>7}" for k in K1_GRID))

wall_rows = []
for label, curve, _k1, m in anchors:
    n = curve[2]
    for kind in ('2m+2', '2m+1'):
        cells = []
        for k1 in K1_GRID:
            wins = 0
            for seed in SEEDS:
                d = random.Random(seed + 7777).randint(1, n - 1)
                if kind == '2m+2':
                    res = L.run_experiment(curve, m, d, k1, seed)
                else:
                    res = L.run_experiment_nod(curve, m, d, k1, seed)
                if res is not None:
                    wins += res[0]
            cells.append(wins)
        wall_rows.append((label, kind, cells))
        print(f"{label:<14}{kind:>6}" + "".join(f"{str(c) + '/5':>7}" for c in cells))
    print()

print("Wall position = largest K1 with >= 3/5 recoveries:")
for label, kind, cells in wall_rows:
    wall = max([k for k, c in zip(K1_GRID, cells) if c >= 3] + [0])
    print(f"  {label:<14} {kind:>6}  K1_wall = {wall}")

# ===========================================================================
# E3 -- predictor comparison on a fresh out-of-sample sweep
# ===========================================================================
banner("E3  H23b/H23c -- predictor comparison (fresh 17-bit curves)")

print("Searching fresh j=0 GLV curves, 17-bit, lam* spread over 6 bins ...")
sys.stdout.flush()
curves = L.search_curves(2 ** 16, 2 ** 17, per_bin=1, nbins=6)
print(f"  found {len(curves)} curves in {time.time() - t_start:.2f}s")

M_SIGS = 10
EFFS = [0.05, 0.12, 0.20, 0.30]
SWEEP_SEEDS = [42, 1234, 9999]

rows = []
print(f"\n{'n':>8}{'lam*':>7}{'nu_hat':>8}{'eff':>7}{'K1':>5}"
      f"{'nu_gh':>8}{'nu_cf':>8}{'nu_proj':>9}{'margin':>9}{'LLL':>6}")
for (p, b, n, lam, G) in curves:
    curve = (p, b, n, lam, G)
    k2_bound = math.isqrt(n) + 1
    for eff in EFFS:
        k1 = max(2, int(round(eff * n / k2_bound)))
        wins, gs = 0, []
        for seed in SWEEP_SEEDS:
            d = random.Random(seed + 7777).randint(1, n - 1)
            g = instance_geometry(curve, M_SIGS, d, k1, seed)
            if g is None:
                continue
            gs.append(g)
            rows.append(g)
            wins += g['ok']
        if not gs:
            continue
        g = gs[0]
        print(f"{n:>8}{g['lam_star']:>7.3f}{g['nu_hat']:>8.3f}{g['eff']:>7.3f}"
              f"{k1:>5}{g['nu_gh']:>8.3f}{g['nu_cf']:>8.3f}{g['nu_proj']:>9.3f}"
              f"{g['margin']:>9.3f}{str(wins) + '/' + str(len(gs)):>6}")
        sys.stdout.flush()

# pool anchors' grid + sweep for the AUC table
pool = rows + e1_rows
labels = [r['ok'] for r in pool]
print(f"\nPooled instances: {len(pool)}  (successes {sum(labels)}, "
      f"failures {len(labels) - sum(labels)})")

print(f"\n{'predictor':<34}{'AUC':>8}{'Spearman':>10}   direction")
for name, key, flip in [
        ("nu_proj (Thread 23, unified)",   'nu_proj', True),
        ("nu_gh   = ||v||/GH(L)",          'nu_gh',   True),
        ("nu_cf   (closed form)",          'nu_cf',   True),
        ("margin  = max||b*_i||/||pi_i(v)||", 'margin', False),
        ("nu_hat  (Thread 20c)",           'nu_hat',  True),
        ("lam*    (Thread 20a)",           'lam_star', False),
        ("eff     = K1*K2/n",              'eff',     True),
        ("sv/pv   (Thread 20 T5)",         'sv_over_pn', False),
]:
    xs = [(-r[key] if flip else r[key]) for r in pool]
    a = L.auc(xs, labels)
    s = L.spearman([r[key] for r in pool], [1.0 if y else 0.0 for y in labels])
    print(f"{name:<34}{a:>8.3f}{s:>10.3f}   "
          f"{'smaller is better' if flip else 'larger is better'}")

# best achievable threshold for each of the two lattice-geometric gaps
for key in ('nu_proj', 'nu_gh'):
    vals = sorted(set(r[key] for r in pool))
    best = (0.0, None)
    for t in vals:
        acc = sum(1 for r in pool if (r[key] < t) == r['ok']) / len(pool)
        if acc > best[0]:
            best = (acc, t)
    print(f"\nBest single threshold on {key}: '{key} < {best[1]:.3f}' "
          f"=> accuracy {100 * best[0]:.1f}%")

# the decision rule nu_gh < 1
tp = sum(1 for r in pool if r['nu_gh'] < 1 and r['ok'])
fp = sum(1 for r in pool if r['nu_gh'] < 1 and not r['ok'])
fn = sum(1 for r in pool if r['nu_gh'] >= 1 and r['ok'])
tn = sum(1 for r in pool if r['nu_gh'] >= 1 and not r['ok'])
print(f"\nRule 'nu_gh < 1 => recoverable':  TP={tp} FP={fp} FN={fn} TN={tn}"
      f"   accuracy={100.0 * (tp + tn) / len(pool):.1f}%")
succ = [r['nu_gh'] for r in pool if r['ok']]
fail = [r['nu_gh'] for r in pool if not r['ok']]
if succ:
    print(f"  nu_gh over successes: min {min(succ):.3f}  max {max(succ):.3f}")
if fail:
    print(f"  nu_gh over failures : min {min(fail):.3f}  max {max(fail):.3f}")

tp = sum(1 for r in pool if r['margin'] >= 1 and r['ok'])
fp = sum(1 for r in pool if r['margin'] >= 1 and not r['ok'])
fn = sum(1 for r in pool if r['margin'] < 1 and r['ok'])
tn = sum(1 for r in pool if r['margin'] < 1 and not r['ok'])
print(f"Rule 'margin >= 1 => recoverable': TP={tp} FP={fp} FN={fn} TN={tn}"
      f"   accuracy={100.0 * (tp + tn) / len(pool):.1f}%")

# ===========================================================================
# E5 -- held-out validation of the nu_proj threshold
# ===========================================================================
banner("E5  Held-out validation: fit the threshold on 17-bit, test on 20-bit")

fit_pool = rows                       # 17-bit sweep only
vals = sorted(set(r['nu_proj'] for r in fit_pool))
best = (0.0, None)
for t in vals:
    acc = sum(1 for r in fit_pool if (r['nu_proj'] < t) == r['ok']) / len(fit_pool)
    if acc > best[0]:
        best = (acc, t)
THRESH = best[1]
print(f"Threshold fitted on the 17-bit sweep only ({len(fit_pool)} instances): "
      f"nu_proj < {THRESH:.3f}  (in-sample accuracy {100 * best[0]:.1f}%)\n")

print("Searching fresh 20-bit curves (never used to fit) ...")
sys.stdout.flush()
ho_curves = L.search_curves(2 ** 19, 2 ** 20, per_bin=1, nbins=5)
print(f"  found {len(ho_curves)} curves\n")

print(f"{'n':>9}{'lam*':>7}{'nu_hat':>8}{'eff':>7}{'K1':>5}{'nu_gh':>8}"
      f"{'nu_proj':>9}{'pred':>6}{'LLL':>6}{'':>4}")
ho = []
for (p, b, n, lam, G) in ho_curves:
    curve = (p, b, n, lam, G)
    k2_bound = math.isqrt(n) + 1
    for eff in (0.06, 0.14, 0.24):
        k1 = max(2, int(round(eff * n / k2_bound)))
        for seed in SWEEP_SEEDS:
            d = random.Random(seed + 7777).randint(1, n - 1)
            g = instance_geometry(curve, M_SIGS, d, k1, seed)
            if g is not None:
                ho.append(g)
        gs = ho[-len(SWEEP_SEEDS):]
        wins = sum(x['ok'] for x in gs)
        g = gs[0]
        pred = g['nu_proj'] < THRESH
        print(f"{n:>9}{g['lam_star']:>7.3f}{g['nu_hat']:>8.3f}{g['eff']:>7.3f}"
              f"{k1:>5}{g['nu_gh']:>8.3f}{g['nu_proj']:>9.3f}"
              f"{('yes' if pred else 'no'):>6}"
              f"{str(wins) + '/' + str(len(gs)):>6}"
              f"{('  ok' if pred == (wins > len(gs) // 2) else '  MISS'):>4}")
        sys.stdout.flush()

ho_labels = [r['ok'] for r in ho]
correct = sum(1 for r in ho if (r['nu_proj'] < THRESH) == r['ok'])
print(f"\nHeld-out: {len(ho)} instances "
      f"({sum(ho_labels)} successes, {len(ho_labels) - sum(ho_labels)} failures)")
print(f"  AUC(nu_proj)               = {L.auc([-r['nu_proj'] for r in ho], ho_labels):.3f}")
print(f"  AUC(nu_gh)                 = {L.auc([-r['nu_gh'] for r in ho], ho_labels):.3f}")
print(f"  AUC(nu_hat)                = {L.auc([-r['nu_hat'] for r in ho], ho_labels):.3f}")
print(f"  accuracy of the fitted rule = {100.0 * correct / len(ho):.1f}%")

# ===========================================================================
# E6 -- the controlled contrast: fix (n-size, m, eff), vary only lam
# ===========================================================================
banner("E6  Controlled contrast -- within a stratum of near-constant nu_gh")

print("nu_gh depends only on (n, m, K1, K2), so inside a stratum of fixed")
print("bit-size and fixed eff it is constant up to the integer rounding of K1")
print("and cannot express any lam dependence.  nu_proj can.  AUC computed")
print("*within* each stratum therefore isolates exactly the lam-geometry term.\n")

strata = {}
for r in rows + ho:
    key = (r['n'].bit_length(), round(r['eff'], 2))
    strata.setdefault(key, []).append(r)

print(f"{'stratum (bits, eff)':<22}{'n':>4}{'wins':>6}{'nu_gh range':>18}"
      f"{'AUC nu_proj':>13}{'AUC nu_hat':>12}")
tot_w, tot_n = 0.0, 0
for key in sorted(strata):
    grp = strata[key]
    labs = [r['ok'] for r in grp]
    if not (0 < sum(labs) < len(labs)):
        continue
    a_p = L.auc([-r['nu_proj'] for r in grp], labs)
    a_h = L.auc([-r['nu_hat'] for r in grp], labs)
    lo = min(r['nu_gh'] for r in grp)
    hi = max(r['nu_gh'] for r in grp)
    npair = sum(labs) * (len(labs) - sum(labs))
    tot_w += a_p * npair
    tot_n += npair
    print(f"{str(key):<22}{len(grp):>4}{sum(labs):>6}"
          f"{f'[{lo:.3f}, {hi:.3f}]':>18}{a_p:>13.3f}{a_h:>12.3f}")
if tot_n:
    print(f"\nPair-weighted within-stratum AUC(nu_proj) = {tot_w / tot_n:.3f}"
          f"   (nu_gh is uninformative here by construction)")

# ===========================================================================
# E4 -- the asymptotic consequence
# ===========================================================================
banner("E4  Asymptotic ceiling implied by H23b")

print("lim_{m->inf} nu_gh = sqrt(2*pi*e/3) * sqrt(eff),  so nu_gh < 1 is")
print("reachable at large m only if eff < 3/(2*pi*e) = %.4f.\n"
      % (3.0 / (2 * math.pi * math.e)))
print(f"{'eff':>7}{'nu_gh(m=inf)':>14}   m needed for nu_gh<1 (n=2^12 / 2^17 / 2^256)")
for eff in [0.03, 0.05, 0.08, 0.1757, 0.20, 0.25]:
    lim = math.sqrt(2 * math.pi * math.e / 3.0) * math.sqrt(eff)
    ms = []
    for nb in (12, 17, 256):
        nn = 2 ** nb
        found = None
        for m in range(2, 4001):
            if L.usvp_gap_closed_form(nn, m, max(2, int(eff * math.isqrt(nn))),
                                      math.isqrt(nn) + 1) < 1.0:
                found = m
                break
        ms.append(str(found) if found else ">4000")
    print(f"{eff:>7.4f}{lim:>14.3f}   {ms[0]:>6} {ms[1]:>6} {ms[2]:>6}")

print(f"\nTotal runtime {time.time() - t_start:.2f}s")
print("Done.")
