"""
GLV-HNP -- Thread 23-CF Part G: does the THEORY-DERIVED a_crit separate C1/C2?

`glv_hnp_nuhat_cf.py` proves lambda_1(L2) is attained at a continued-fraction
convergent of lam/n, giving

    nu_hat^2 = min_j (x_j^2 + y_j^2) >= 2*c_{j*},   c_j = q_j*eta_j/n,
    c_j = 1 / (q_{j+1}/q_j + eta_{j+1}/eta_j)

so nu_hat is small exactly when the partial quotient a_{j*+1} at the critical
scale sqrt(n*K2/K1) is large.  That predicts a purely arithmetic separator:

    a_crit = the partial quotient of lam/n immediately after the convergent
             that minimises (S1*eta)^2 + (S2*q)^2.

This script replicates the Exp S protocol (K1=72, m=12, 6 seeds, 20-bit j=0 CM
curves) exactly as `glv_hnp_nuhat_vs_c1c2.py` does, and puts five invariants
head to head on the SAME sample:

    nu_hat  -- the 2026-07-29 empirical separator (AUC 0.935 reported)
    a_crit  -- NEW, derived from the closed form (arithmetic, no lattice)
    c       -- q*eta/n at the winning convergent (the bound quantity)
    max_a   -- June invariant, falsified by Exp S (largest PQ anywhere)
    mu      -- min(lam,n-lam)/n, falsified 2026-07-29

FALSIFIER for the closed-form story: if a_crit's AUC sits at the trivial
baseline while nu_hat's is high, then the "large partial quotient at the
critical scale" reading of nu_hat is wrong and only the full 2-D minimum
matters.  If a_crit tracks nu_hat, nu_hat has been reduced to arithmetic.

Run: python3 glv_hnp_acrit_vs_c1c2.py
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
mu_of = _t20a.mu_of
harvest = _t20d.harvest
auc = _t20d.auc
best_cut_accuracy = _t20d.best_cut_accuracy
max_partial_quotient = _t20d.max_partial_quotient
nu_hat_cf = _cf.nu_hat_cf
spearman = _cf.spearman

# Exp S protocol (2026-06-29), identical to glv_hnp_nuhat_vs_c1c2.py
K1_FIXED = 72
M_FIXED = 12
SEEDS_N = 6
N_CURVES = 100


def main():
    print("=" * 78)
    print("Thread 23-CF Part G: a_crit (theory-derived) vs the June C1/C2 classes")
    print("=" * 78)
    print(f"Exp S protocol: K1={K1_FIXED}, m={M_FIXED}, {SEEDS_N} seeds, "
          f"{N_CURVES} fresh 20-bit j=0 CM curves")
    print("C2 = recovers d on all seeds; C1 = does not (Exp S convention)")

    t0 = time.time()
    curves = harvest(N_CURVES)
    print(f"\n  harvested {len(curves)} curves in {time.time()-t0:.1f}s")

    t1 = time.time()
    for c in curves:
        n = c['n']
        k2 = math.isqrt(n) + 1
        wins = 0
        for s in range(SEEDS_N):
            d = random.Random(s + 7777).randint(1, n - 1)
            wins += run_experiment(c['p'], n, c['lam'], c['G'], M_FIXED, d,
                                   K1_FIXED, k2, s)
        c['wins'] = wins
        c['C2'] = (wins == SEEDS_N)
        _, _, c['nu_hat'] = rival_sublattice_nu(n, c['lam'], K1_FIXED, k2)
        nh_cf, det = nu_hat_cf(n, c['lam'], K1_FIXED, k2, want_detail=True)
        c['nu_hat_cf'] = nh_cf
        c['a_crit'] = det['a_next'] if det['a_next'] is not None else 0
        c['c'] = det['c']
        c['q_star'] = det['q']
        c['mu'] = mu_of(c['lam'], n)
        c['max_a'] = max_partial_quotient(c['lam'], n)
    print(f"  {len(curves)*SEEDS_N} LLL trials in {time.time()-t1:.1f}s")

    # sanity: the closed form must agree with Lagrange-Gauss on this sample too
    worst = max(abs(c['nu_hat_cf'] - c['nu_hat']) for c in curves)
    print(f"  max |nu_hat(CF) - nu_hat(Gauss)| on this sample = {worst:.3e}")

    c1 = [c for c in curves if not c['C2']]
    c2 = [c for c in curves if c['C2']]
    base = max(len(c1), len(c2)) / len(curves)
    print(f"\n  C1 (fails): {len(c1)}/{len(curves)}   C2 (6/6): {len(c2)}/{len(curves)}")
    print(f"  trivial baseline (always predict majority) = {base*100:.1f}%")

    print("\n--- Head to head on the same sample ---")
    print(f"{'invariant':>10} {'C1 range':>20} {'C2 range':>20} {'AUC':>7} {'best acc':>9}")
    for key in ('nu_hat', 'a_crit', 'c', 'max_a', 'mu'):
        s1 = [c[key] for c in c1]
        s2 = [c[key] for c in c2]
        a = auc([c[key] for c in curves], [c['C2'] for c in curves])
        a_or = max(a, 1 - a)
        acc, cut = best_cut_accuracy([c[key] for c in curves],
                                     [c['C2'] for c in curves])
        fmt = (lambda v: f"{v:.3f}") if key in ('nu_hat', 'c', 'mu') else (lambda v: f"{v:g}")
        r1 = f"[{fmt(min(s1))},{fmt(max(s1))}]" if s1 else "-"
        r2 = f"[{fmt(min(s2))},{fmt(max(s2))}]" if s2 else "-"
        print(f"{key:>10} {r1:>20} {r2:>20} {a_or:>7.3f} {acc*100:>8.1f}%")

    print("\n--- a_crit response (C2 rate and mean wins/6) ---")
    print(f"{'a_crit':>8} {'count':>6} {'C2 rate':>9} {'mean wins/6':>12} "
          f"{'mean nu_hat':>12}")
    buckets = {}
    for c in curves:
        a = c['a_crit']
        key = 1 if a <= 1 else 2 if a == 2 else 3 if a <= 4 else 4 if a <= 8 else 5
        buckets.setdefault(key, []).append(c)
    labels = {1: '1', 2: '2', 3: '3-4', 4: '5-8', 5: '>8'}
    for k in sorted(buckets):
        v = buckets[k]
        print(f"{labels[k]:>8} {len(v):>6} {sum(x['C2'] for x in v)/len(v):>9.2f} "
              f"{sum(x['wins'] for x in v)/len(v):>12.2f} "
              f"{sum(x['nu_hat'] for x in v)/len(v):>12.3f}")

    print("\n--- rank correlations with wins/6 ---")
    w = [c['wins'] for c in curves]
    for key in ('nu_hat', 'a_crit', 'c', 'max_a', 'mu'):
        print(f"  spearman({key:>7}, wins) = "
              f"{spearman([c[key] for c in curves], w):+.3f}")

    print("\n--- does a_crit add anything beyond nu_hat (and vice versa)? ---")
    print(f"  spearman(a_crit, nu_hat) = "
          f"{spearman([c['a_crit'] for c in curves], [c['nu_hat'] for c in curves]):+.3f}")
    print(f"  spearman(max_a,  nu_hat) = "
          f"{spearman([c['max_a'] for c in curves], [c['nu_hat'] for c in curves]):+.3f}")

    print("\nDone.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
