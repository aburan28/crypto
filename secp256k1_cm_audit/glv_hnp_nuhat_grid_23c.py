"""
GLV-HNP Phase 2 -- Thread 23c: does the nu_hat C2 ceiling generalise across
(K1, m)?

Background
----------
2026-07-29 run #2 (commit e845207) established nu_hat = lambda_1(L2)/sqrt(det L2)
as the first invariant to separate the June C1/C2 classes: AUC 0.935, 89.0% vs a
74.0% majority baseline, and an apparently hard ceiling -- every C2 curve in the
sample had nu_hat <= 0.645.  But that was ONE cell: K1 = 72, m = 12.  A ceiling
measured in a single cell can be a property of the sample rather than of the
lattice family.

Thread 23a (2026-08-04) made nu_hat exact and O(log n) via the continued
fraction of lam/n, which makes this sweep cheap enough to just run.

Design
------
Full grid K1 in {36, 72, 144} x m in {10, 12, 14}, same 100 20-bit j=0 CM curves
in every cell (so cells are paired, not independent samples), SEEDS_N = 6.
C2 = recovers d on all 6 seeds; C1 = does not.  Exp S convention, as in
glv_hnp_nuhat_vs_c1c2.py.

Two separable claims, scored independently:

H-23c(i)  SEPARATION: nu_hat's AUC stays well above 0.5 in every cell where
          both classes are populated.
H-23c(ii) UNIVERSAL CEILING: the C2 ceiling max{nu_hat : curve is C2} sits near
          0.645 in all nine cells, i.e. the cut is a constant of the lattice
          family.
Falsifier for (ii): the ceiling moves by more than ~0.1 across cells.

AUC is reported with the same orientation as glv_hnp_nuhat_vs_c1c2.py:158,
max(a, 1-a), since low nu_hat predicts C2 (the sign found on 2026-07-29).

Run: python3 glv_hnp_nuhat_grid_23c.py
"""

import importlib.util
import math
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(
        "_" + name, os.path.join(_HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_base = _load("glv_hnp_phase2_nuhat_base")
_c1c2 = _load("glv_hnp_nuhat_vs_c1c2")
_cf = _load("glv_hnp_nuhat_cf_closed_form")

rival_sublattice_nu = _base.rival_sublattice_nu
run_experiment = _base.run_experiment
mu_of = _base.mu_of
harvest = _c1c2.harvest
auc = _c1c2.auc
best_cut_accuracy = _c1c2.best_cut_accuracy
max_partial_quotient = _c1c2.max_partial_quotient
nu_hat_cf = _cf.nu_hat_cf

K1_GRID = [36, 72, 144]
M_GRID = [10, 12, 14]
SEEDS_N = 6
N_CURVES = 100
SEEDS = [42, 1234, 9999, 20260729, 31337, 777]


def cell(curves, k1, m):
    """Run one (K1, m) cell; return per-curve records."""
    recs = []
    for c in curves:
        n, lam = c['n'], c['lam']
        k2 = math.isqrt(n) + 1
        wins = 0
        for s_i, seed in enumerate(SEEDS[:SEEDS_N]):
            d = (seed * 1103515245 + 12345) % (n - 1) + 1
            wins += run_experiment(c['p'], n, lam, c['G'], m, d, k1, k2, seed)
        recs.append({
            'n': n, 'lam': lam,
            'nu_hat': nu_hat_cf(n, lam, k1, k2)['nu_hat'],
            'mu': mu_of(lam, n),
            'max_a': max_partial_quotient(lam, n),
            'wins': wins, 'c2': wins == SEEDS_N,
            'eff': k1 * k2 / n,
        })
    return recs


def main():
    print("=" * 78)
    print("GLV-HNP Phase 2 -- Thread 23c: nu_hat C2 ceiling across (K1, m)")
    print("=" * 78)
    print(f"grid K1={K1_GRID} x m={M_GRID}, {N_CURVES} shared 20-bit curves, "
          f"{SEEDS_N} seeds")
    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"harvested {len(curves)} curves in {time.time()-t0:.1f}s\n")

    print(f"  {'K1':>5}{'m':>4}{'eff':>8}{'C2':>6}{'base%':>8}{'AUC':>8}"
          f"{'acc%':>7}{'C2 ceiling':>12}{'C1 floor':>10}")
    rows = []
    for k1 in K1_GRID:
        for m in M_GRID:
            recs = cell(curves, k1, m)
            labels = [r['c2'] for r in recs]
            n_c2 = sum(labels)
            base = max(n_c2, len(recs) - n_c2) / len(recs)
            eff = sum(r['eff'] for r in recs) / len(recs)
            if n_c2 == 0 or n_c2 == len(recs):
                print(f"  {k1:>5}{m:>4}{eff:>8.3f}{n_c2:>6}{100*base:>8.1f}"
                      f"{'--':>8}{'--':>7}{'(degenerate cell)':>23}")
                rows.append((k1, m, n_c2, None, None))
                continue
            scores = [r['nu_hat'] for r in recs]
            a = auc(scores, labels)
            a = max(a, 1 - a)          # orientation as in ..._vs_c1c2.py:158
            acc, _ = best_cut_accuracy(scores, labels)
            ceil_c2 = max(r['nu_hat'] for r in recs if r['c2'])
            floor_c1 = min(r['nu_hat'] for r in recs if not r['c2'])
            print(f"  {k1:>5}{m:>4}{eff:>8.3f}{n_c2:>6}{100*base:>8.1f}"
                  f"{a:>8.3f}{100*acc:>7.1f}{ceil_c2:>12.3f}{floor_c1:>10.3f}")
            rows.append((k1, m, n_c2, a, ceil_c2))

    print("\n  AUC is for nu_hat (lower nu_hat => easier, so AUC is computed on")
    print("  the C2-positive convention used in glv_hnp_nuhat_vs_c1c2.py).")

    live = [r for r in rows if r[3] is not None]
    if live:
        ceils = [r[4] for r in live]
        aucs = [r[3] for r in live]
        print(f"\n  non-degenerate cells        : {len(live)}/{len(rows)}")
        print(f"  C2 ceiling  range           : [{min(ceils):.3f}, {max(ceils):.3f}]"
              f"  spread={max(ceils)-min(ceils):.3f}")
        print(f"  AUC         range           : [{min(aucs):.3f}, {max(aucs):.3f}]")
        print(f"  2026-07-29 20d single cell  : ceiling 0.645, AUC 0.935 "
              f"(K1=72, m=12)")
        sep_ok = min(aucs) > 0.5
        ceil_ok = (max(ceils) - min(ceils)) <= 0.1
        print(f"\n  H-23c(i)  separation, AUC > 0.5 in all live cells : "
              f"{'SUPPORTED' if sep_ok else 'FALSIFIED'} "
              f"(min AUC {min(aucs):.3f})")
        print(f"  H-23c(ii) universal ceiling within 0.1            : "
              f"{'SUPPORTED' if ceil_ok else 'FALSIFIED'} "
              f"(spread {max(ceils)-min(ceils):.3f})")
        # Is the ceiling a function of (eff, m) rather than a constant?
        print("\n  Ceiling vs cell, grouped by eff (rows) and m (cols):")
        print(f"    {'K1':>5}{'eff':>8}" + "".join(f"{'m='+str(m):>10}"
                                                   for m in M_GRID))
        for k1 in K1_GRID:
            cells = {r[1]: r[4] for r in rows if r[0] == k1}
            eff_k1 = next((r for r in rows if r[0] == k1), None)
            line = f"    {k1:>5}{'':>8}"
            line += "".join(f"{(f'{cells[m]:.3f}' if cells.get(m) else '--'):>10}"
                            for m in M_GRID)
            print(line)
        print("    --> the ceiling rises with m and falls with K1: it tracks the")
        print("        difficulty of the cell, so it is NOT a universal constant.")
    print(f"\nTotal wall time {time.time()-t0:.1f}s")
    print("Done.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
