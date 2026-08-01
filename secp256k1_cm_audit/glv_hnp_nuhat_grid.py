"""
GLV-HNP Thread 23, second sub-task: does the nu_hat cut generalise off the
single (K1, m) cell of Thread 20d?

Thread 20d (2026-07-29 run #2) measured the C1/C2 separation at exactly one
operating point -- K1 = 72, m = 12, 6 seeds -- and found

    AUC 0.935, every C2 curve has nu_hat <= 0.645, no curve with
    nu_hat > 0.645 is ever C2.

The proposed check was: re-run at K1 in {36, 72, 144} and m in {10, 12, 14}.
If the C2 ceiling stays near nu_hat ~ 0.65 across all nine cells the cut is a
property of the lattice family; if it drifts with (K1, m) it is a property of
that one sample and the "ceiling" is not a constant.

Note nu_hat itself depends on K1 (via S_K1 = n // K1), so the three K1 columns
are genuinely different predictors, not a relabelling.

Run: python3 glv_hnp_nuhat_grid.py
"""

import math
import random
import sys
import time

import importlib.util
_d = __file__.rsplit("/", 1)[0]

def _load(name, fname):
    s = importlib.util.spec_from_file_location(name, _d + "/" + fname)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

_t20a = _load("_t20a", "glv_hnp_phase2_nuhat_base.py")
_t20d = _load("_t20d", "glv_hnp_nuhat_vs_c1c2.py")

run_experiment = _t20a.run_experiment
rival_sublattice_nu = _t20a.rival_sublattice_nu
mu_of = _t20a.mu_of
harvest = _t20d.harvest
auc = _t20d.auc

K1_GRID = [36, 72, 144]
M_GRID = [10, 12, 14]
SEEDS_N = 6
N_CURVES = 60

def main():
    print("=" * 78)
    print("GLV-HNP Thread 23b: does the nu_hat cut hold across (K1, m)?")
    print("=" * 78)
    print(f"{N_CURVES} fresh 20-bit j=0 CM curves, {SEEDS_N} seeds, "
          f"K1 in {K1_GRID}, m in {M_GRID}")
    print("C2 = recovers d on all seeds (Exp S convention)")

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"\nharvested {len(curves)} curves in {time.time()-t0:.1f}s")
    for c in curves:
        c['k2'] = math.isqrt(c['n']) + 1

    print(f"\n{'K1':>5} {'m':>4} {'C2/N':>8} {'AUC':>7} {'C2 ceiling':>11}"
          f" {'C1 floor':>9} {'overlap':>15} {'sec/cell':>9}")
    rows = []
    for k1 in K1_GRID:
        for c in curves:
            _, _, c[f'nu{k1}'] = rival_sublattice_nu(c['n'], c['lam'],
                                                     k1, c['k2'])
        for m in M_GRID:
            t1 = time.time()
            for c in curves:
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, c['n'] - 1)
                    wins += run_experiment(c['p'], c['n'], c['lam'], c['G'],
                                           m, d, k1, c['k2'], s)
                c['C2'] = (wins == SEEDS_N)
            nus = [c[f'nu{k1}'] for c in curves]
            labs = [c['C2'] for c in curves]
            a = auc(nus, labs)
            a_or = max(a, 1 - a)
            c2 = [v for v, l in zip(nus, labs) if l]
            c1 = [v for v, l in zip(nus, labs) if not l]
            ceil_ = max(c2) if c2 else float('nan')
            floor_ = min(c1) if c1 else float('nan')
            ov = (f"[{floor_:.3f},{ceil_:.3f}]"
                  if c1 and c2 and floor_ <= ceil_ else "none")
            print(f"{k1:>5} {m:>4} {sum(labs):>4}/{len(curves):<3}"
                  f" {a_or:>7.3f} {ceil_:>11.3f} {floor_:>9.3f} {ov:>15}"
                  f" {time.time()-t1:>9.1f}")
            rows.append((k1, m, sum(labs), a_or, ceil_, floor_))

    print("\n--- summary ---")
    valid = [r for r in rows if 0 < r[2] < len(curves)]
    if valid:
        ceilings = [r[4] for r in valid]
        aucs = [r[3] for r in valid]
        print(f"  {len(valid)}/9 cells non-degenerate (both classes present)")
        print(f"  C2 ceiling: min={min(ceilings):.3f} max={max(ceilings):.3f} "
              f"mean={sum(ceilings)/len(ceilings):.3f}")
        print(f"  AUC:        min={min(aucs):.3f} max={max(aucs):.3f} "
              f"mean={sum(aucs)/len(aucs):.3f}")
        print("  Thread 20d reference (K1=72, m=12, 100 curves): "
              "ceiling 0.645, AUC 0.935")
    else:
        print("  all cells degenerate -- (K1, m) grid too easy or too hard")

    print("\nDone.")

if __name__ == '__main__':
    sys.exit(main())
