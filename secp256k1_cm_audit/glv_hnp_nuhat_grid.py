"""
GLV-HNP -- Thread 23-CF Part H: does the nu_hat cut generalise across (K1, m)?

The 2026-07-29 entry (recovered from commit e845207) reported a C2 ceiling at
nu_hat ~ 0.645 -- every curve that recovered d on all seeds had nu_hat below it
-- but measured it at ONE cell only (K1=72, m=12).  Its own next-step proposal
was to check K1 in {36,72,144} x m in {10,12,14}: if the ceiling is a property
of the lattice family rather than of that sample, it should persist.

This also tests the theory-derived a_crit (see glv_hnp_nuhat_cf.py) in each
cell, and reports whether the cut value drifts with eff = K1*K2/n.

FALSIFIER: if the C2 ceiling moves by more than ~0.15 across the nine cells, the
"nu_hat < ceiling" rule is sample-specific and only the ordering (AUC) is real.

Run: python3 glv_hnp_nuhat_grid.py
"""

import math
import random
import sys
import time

import importlib.util
_here = __file__.rsplit("/", 1)[0]

def _load(name, path):
    s = importlib.util.spec_from_file_location(name, _here + "/" + path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

_t20a = _load("_t20a", "glv_hnp_phase2_lambda_threshold_20a.py")
_t20d = _load("_t20d", "glv_hnp_nuhat_vs_c1c2.py")
_cf = _load("_cf", "glv_hnp_nuhat_cf.py")

run_experiment = _t20a.run_experiment
rival_sublattice_nu = _t20a.rival_sublattice_nu
harvest = _t20d.harvest
auc = _t20d.auc
best_cut_accuracy = _t20d.best_cut_accuracy
max_partial_quotient = _t20d.max_partial_quotient
nu_hat_cf = _cf.nu_hat_cf

K1_GRID = (36, 72, 144)
M_GRID = (10, 12, 14)
SEEDS_N = 6
N_CURVES = 100


def main():
    print("=" * 78)
    print("Thread 23-CF Part H: nu_hat / a_crit across a (K1, m) grid")
    print("=" * 78)
    print(f"{N_CURVES} fresh 20-bit j=0 CM curves, {SEEDS_N} seeds per cell, "
          f"K2 = isqrt(n)+1")
    print("C2 = recovers d on all seeds.  Ceiling = max nu_hat over C2 curves.")

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"\n  harvested {len(curves)} curves in {time.time()-t0:.1f}s\n")

    print(f"{'K1':>5} {'m':>4} {'eff':>7} {'C2/N':>8} {'base':>7} "
          f"{'nu AUC':>7} {'nu acc':>7} {'ceiling':>8} {'ac AUC':>7} {'ma AUC':>7}")
    rows = []
    for K1 in K1_GRID:
        for M in M_GRID:
            recs = []
            for c in curves:
                n, lam = c['n'], c['lam']
                k2 = math.isqrt(n) + 1
                wins = 0
                for s in range(SEEDS_N):
                    d = random.Random(s + 7777).randint(1, n - 1)
                    wins += run_experiment(c['p'], n, lam, c['G'], M, d,
                                           K1, k2, s)
                _, _, nh = rival_sublattice_nu(n, lam, K1, k2)
                _, det = nu_hat_cf(n, lam, K1, k2, want_detail=True)
                recs.append({
                    'nu_hat': nh,
                    'a_crit': det['a_next'] if det['a_next'] is not None else 0,
                    'max_a': max_partial_quotient(lam, n),
                    'C2': wins == SEEDS_N,
                    'wins': wins,
                })
            n0 = curves[0]['n']
            eff = K1 * (math.isqrt(n0) + 1) / n0
            c2 = [r for r in recs if r['C2']]
            c1 = [r for r in recs if not r['C2']]
            base = max(len(c1), len(c2)) / len(recs)
            labels = [r['C2'] for r in recs]
            a_nu = auc([r['nu_hat'] for r in recs], labels)
            a_nu = max(a_nu, 1 - a_nu)
            acc_nu, _ = best_cut_accuracy([r['nu_hat'] for r in recs], labels)
            a_ac = auc([r['a_crit'] for r in recs], labels)
            a_ac = max(a_ac, 1 - a_ac)
            a_ma = auc([r['max_a'] for r in recs], labels)
            a_ma = max(a_ma, 1 - a_ma)
            ceil = max((r['nu_hat'] for r in c2), default=float('nan'))
            print(f"{K1:>5} {M:>4} {eff:>7.4f} {len(c2):>4}/{len(recs):<3} "
                  f"{base*100:>6.1f}% {a_nu:>7.3f} {acc_nu*100:>6.1f}% "
                  f"{ceil:>8.3f} {a_ac:>7.3f} {a_ma:>7.3f}")
            rows.append((K1, M, eff, len(c2), base, a_nu, acc_nu, ceil, a_ac, a_ma))

    print()
    valid = [r for r in rows if r[3] >= 3 and not math.isnan(r[7])]
    if valid:
        ceils = [r[7] for r in valid]
        print(f"  cells with >=3 C2 curves: {len(valid)}/{len(rows)}")
        print(f"  C2 ceiling across those cells: min={min(ceils):.3f} "
              f"max={max(ceils):.3f} spread={max(ceils)-min(ceils):.3f}")
        print(f"  nu_hat AUC: min={min(r[5] for r in valid):.3f} "
              f"max={max(r[5] for r in valid):.3f}")
        print(f"  a_crit AUC: min={min(r[8] for r in valid):.3f} "
              f"max={max(r[8] for r in valid):.3f}")
        print(f"  max_a  AUC: min={min(r[9] for r in valid):.3f} "
              f"max={max(r[9] for r in valid):.3f}")
        print()
        verdict = ("CONFIRMED -- ceiling stable"
                   if max(ceils) - min(ceils) <= 0.15
                   else "FALSIFIED -- ceiling is cell-specific; only the ordering is real")
        print(f"  VERDICT (ceiling stability, tolerance 0.15): {verdict}")
    print("\nDone.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
