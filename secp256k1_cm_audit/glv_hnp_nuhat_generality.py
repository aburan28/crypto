"""
GLV-HNP -- Thread 23, parts F/G.

F  Does the C2 ceiling nu_hat ~ 0.645 survive a change of (K1, m)?
   20d used ONE cell (K1=72, m=12).  Sweep K1 in {36,72,144} x m in {10,12,14}.
   If the ceiling is a property of the lattice family the cut stays put; if it
   is a property of the sample it moves.

G  How local is nu_hat in the CF expansion?  Part D showed the scale-local
   partial quotient a_{j*+1} has AUC 0.616 -- WORSE than scale-free max_a
   (0.748) -- and a non-monotone response.  Test whether nu_hat is recovered
   by a WINDOW of convergents around the critical index, and decompose the two
   mechanisms that make nu_hat small.

Run: python3 glv_hnp_nuhat_generality.py
"""

import importlib.util
import math
import random
import sys
import time

import sympy

_d = __file__.rsplit("/", 1)[0]
_spec = importlib.util.spec_from_file_location("_core", _d + "/glv_hnp_nuhat_core.py")
_core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_core)
_spec2 = importlib.util.spec_from_file_location(
    "_cf", _d + "/glv_hnp_nuhat_closed_form.py")
_cf = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(_cf)

scales = _core.scales
rival_sublattice_nu = _core.rival_sublattice_nu
run_experiment = _core.run_experiment
mu_of = _core.mu_of

cf_convergents = _cf.cf_convergents
centered_residue = _cf.centered_residue
objective = _cf.objective
critical_index = _cf.critical_index
max_partial_quotient = _cf.max_partial_quotient
harvest = _cf.harvest
auc = _cf.auc
best_cut_accuracy = _cf.best_cut_accuracy
spearman = _cf.spearman
pearson = _cf.pearson

N_CURVES = 100
SEEDS_N = 6


def nu_hat_window(n, lam, k1_bound, k2_bound, half_width):
    """nu_hat restricted to convergents within +/- half_width of the critical
    index j*.  half_width = None means all convergents (the exact value)."""
    s_k1, _, s_k2, _ = scales(n, k1_bound, k2_bound)
    j_star, _, _, _ = critical_index(n, lam, k1_bound, k2_bound)
    _, convs = cf_convergents(lam, n)
    best = (n * s_k1) ** 2
    for j, (_, q_j) in enumerate(convs):
        if q_j <= 0:
            continue
        if half_width is not None and abs(j - j_star) > half_width:
            continue
        if s_k2 * q_j >= n * s_k1:
            break
        best = min(best, objective(q_j, lam, n, s_k1, s_k2))
    return math.sqrt(best / (n * s_k1 * s_k2))


def part_f():
    print("=" * 78)
    print("PART F -- does the nu_hat C2 ceiling survive a change of (K1, m)?")
    print("=" * 78)
    curves = harvest(N_CURVES)
    print(f"  {len(curves)} curves (same deterministic 20-bit sample as 20d)\n")

    print(f"{'K1':>5} {'m':>4} {'C2/n':>7} {'base':>7} {'AUC nu':>8} "
          f"{'acc nu':>8} {'C2 ceil':>8} {'C1 floor':>9} {'AUC max_a':>10}")
    grid = []
    for k1 in (36, 72, 144):
        for m in (10, 12, 14):
            rows = []
            for c in curves:
                n = c['n']
                k2 = math.isqrt(n) + 1
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, n - 1)
                    wins += run_experiment(c['p'], n, c['lam'], c['G'], m, d,
                                           k1, k2, s)
                _, _, nh = rival_sublattice_nu(n, c['lam'], k1, k2)
                rows.append({'nu': nh, 'C2': wins == SEEDS_N, 'wins': wins,
                             'max_a': max_partial_quotient(c['lam'], n)})
            lab = [r['C2'] for r in rows]
            n_c2 = sum(lab)
            base = max(n_c2, len(rows) - n_c2) / len(rows)
            a_nu = auc([r['nu'] for r in rows], lab)
            a_ma = auc([float(r['max_a']) for r in rows], lab)
            ceil_ = max((r['nu'] for r in rows if r['C2']), default=float('nan'))
            floor_ = min((r['nu'] for r in rows if not r['C2']),
                         default=float('nan'))
            acc = best_cut_accuracy([r['nu'] for r in rows], lab)
            grid.append((k1, m, n_c2, base, max(a_nu, 1 - a_nu), acc, ceil_,
                         floor_, max(a_ma, 1 - a_ma)))
            print(f"{k1:>5} {m:>4} {n_c2:>6}% {base*100:>6.0f}% "
                  f"{max(a_nu,1-a_nu):>8.3f} {acc*100:>7.1f}% "
                  f"{ceil_:>8.3f} {floor_:>9.3f} {max(a_ma,1-a_ma):>10.3f}")

    live = [g for g in grid if 0 < g[2] < N_CURVES]
    if live:
        ceils = [g[6] for g in live]
        aucs = [g[4] for g in live]
        print(f"\n  cells with both classes present : {len(live)}/9")
        print(f"  C2 ceiling  range               : "
              f"[{min(ceils):.3f}, {max(ceils):.3f}]  "
              f"(20d single-cell value 0.645)")
        print(f"  AUC(nu_hat) range               : "
              f"[{min(aucs):.3f}, {max(aucs):.3f}]")
        print(f"  AUC(nu_hat) beats AUC(max_a) in : "
              f"{sum(1 for g in live if g[4] > g[8])}/{len(live)} cells")
    return curves


def part_g(curves):
    print()
    print("=" * 78)
    print("PART G -- how local is nu_hat in the CF expansion?")
    print("=" * 78)
    k1 = 72
    for c in curves:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        s_k1, _, s_k2, _ = scales(n, k1, k2)
        _, _, c['nu'] = rival_sublattice_nu(n, lam, k1, k2)
        for w in (0, 1, 2, 3):
            c[f'w{w}'] = nu_hat_window(n, lam, k1, k2, w)
        j_star, q_star, a_star, target = critical_index(n, lam, k1, k2)
        c['a_star'], c['q_star'], c['target'] = a_star, q_star, target
        # the two components of the winning objective
        _, convs = cf_convergents(lam, n)
        best, comp = float('inf'), None
        for j, (_, q_j) in enumerate(convs):
            if q_j <= 0 or s_k2 * q_j >= n * s_k1:
                continue
            v = objective(q_j, lam, n, s_k1, s_k2)
            if v < best:
                best = v
                comp = (s_k1 * centered_residue(q_j, lam, n), s_k2 * q_j)
        c['u'], c['v'] = comp if comp else (0, 1)

    print("  nu_hat computed from convergents within +/- w of j*:\n")
    print(f"{'window w':>10} {'exact match':>13} {'max rel err':>13} "
          f"{'pearson r':>11}")
    for w in (0, 1, 2, 3):
        vals = [c[f'w{w}'] for c in curves]
        ex = sum(1 for c in curves
                 if abs(c[f'w{w}'] - c['nu']) <= 1e-12 * max(1.0, c['nu']))
        rel = max(abs(c[f'w{w}'] - c['nu']) / c['nu'] for c in curves)
        print(f"{'+/-' + str(w):>10} {ex:>10}/100 {rel:>13.4f} "
              f"{pearson(vals, [c['nu'] for c in curves]):>11.4f}")

    print("\n  Two mechanisms for small nu_hat (winning vector components):")
    print("    u = S_K1 * delta_j  (approximation quality)")
    print("    v = S_K2 * q_j      (denominator cost)")
    print(f"\n{'nu_hat band':>14} {'curves':>7} {'mean u/v':>9} "
          f"{'mean a_{j*+1}':>14} {'mean q*/scale':>14}")
    bands = [(0.0, 0.45), (0.45, 0.645), (0.645, 0.80), (0.80, 2.0)]
    for lo, hi in bands:
        grp = [c for c in curves if lo <= c['nu'] < hi]
        if not grp:
            continue
        print(f"{'[' + f'{lo:.2f},{hi:.2f}' + ')':>14} {len(grp):>7} "
              f"{sum(c['u']/c['v'] for c in grp)/len(grp):>9.3f} "
              f"{sum(c['a_star'] for c in grp)/len(grp):>14.2f} "
              f"{sum(c['q_star']/c['target'] for c in grp)/len(grp):>14.3f}")

    r_a = spearman([float(c['a_star']) for c in curves],
                   [c['nu'] for c in curves])
    r_q = spearman([c['q_star'] / c['target'] for c in curves],
                   [c['nu'] for c in curves])
    print(f"\n  spearman(a_{{j*+1}}, nu_hat) = {r_a:+.3f}")
    print(f"  spearman(q*/scale,  nu_hat) = {r_q:+.3f}")

    # integer-exact check: nu_hat^2 is rational, so the comparison
    # nu_hat <= cut can be done in exact integer arithmetic
    print("\n  Integer-exactness: nu_hat^2 = N/D with "
          "N = min_j(S_K1^2 delta_j^2 + S_K2^2 q_j^2), D = n S_K1 S_K2.")
    t0 = time.time()
    for c in curves:
        k2 = math.isqrt(c['n']) + 1
        nu_hat_window(c['n'], c['lam'], k1, k2, None)
    t_cf = time.time() - t0
    t0 = time.time()
    for c in curves:
        k2 = math.isqrt(c['n']) + 1
        rival_sublattice_nu(c['n'], c['lam'], k1, k2)
    t_lg = time.time() - t0
    print(f"  100 evaluations: CF {t_cf*1e3:.1f} ms   "
          f"Lagrange-Gauss(float) {t_lg*1e3:.1f} ms")


def main():
    t0 = time.time()
    curves = part_f()
    part_g(curves)
    print(f"\nDone in {time.time()-t0:.1f}s.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
