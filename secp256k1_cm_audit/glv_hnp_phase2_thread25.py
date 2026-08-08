"""
GLV-HNP Phase 2, Thread 25: is mu a second coordinate, independent of NU?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry.  W5/W6 there
established that NU (exact BDD nearest-plane certificate) and mu =
lambda_1(L2) (equivalently nu_hat = mu/sqrt(det L2)) are UNCORRELATED at
fixed eff, yet both separate recovery from failure.  Two independent
predictors of the same event means recovery = f(NU, X) with X ~ mu-driven.

  H25: within the ambiguous NU band [1.04, 2.20] (the W4 bracket at 17 bits,
       where nearest-plane alone gives no answer: sufficient NU < 1.040,
       necessary NU > 2.199), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  Falsifier: if AUC collapses toward 0.5 inside the band, mu's apparent
  power in W5 is entirely mediated by NU (a stratification artifact of
  eff, since low-eff strata also tend to have small NU), and the closed
  form should be retired as a redundant reparametrization of NU.

Secondary (W1b): the GS profile head is m exact copies of lambda_1(L2) and
the step to the second block vanishes right at the K1 wall.  Define

  step = log2(||b*_{m+1}||) - log2(||b*_1||)     (0-indexed: prof[m]-prof[0])

and test whether step -> 0 predicts the wall better than NU or mu, on the
same U2 K1-grid Thread 24 used (2 curves x 11 K1 values x 5 seeds, exact GS).

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import find_generator, lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, HIST, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

BAND_LO, BAND_HI = 1.040, 2.199   # W4 bracket at 17 bits


def log2_step(row):
    m = row['m']
    prof = row['prof']
    if prof[0] <= 0 or prof[m] <= 0:
        return float('nan')
    return math.log2(prof[m]) - math.log2(prof[0])


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a second coordinate independent of NU?")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # Part A. H25 main test: reproduce the 500-instance 17-bit table
    # (glv_hnp_phase2_gsprofile_strat.py W5/W6) and stratify by NU band.
    # -----------------------------------------------------------------------
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    M17 = 12
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M17, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'm': M17})
                rows.append(r)
    print(f"\n{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s  [reproduces the Thread 24 W5 table]")

    print("\n" + "-" * 78)
    print(f"H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{BAND_LO}, {BAND_HI}]")
    print("-" * 78)
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band N = {len(band)}  ({len(pos)} recovered / {len(neg)} failed)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_eff = auc([r['eff'] for r in pos], [r['eff'] for r in neg])
        print(f"  AUC(-mu)      = {a_mu:.4f}   {'H25 HOLDS' if a_mu >= 0.8 else 'H25 FALSIFIED'}")
        print(f"  AUC(-nu_hat)  = {a_nh:.4f}")
        print(f"  AUC(-NU)      = {a_nu:.4f}   (expect ~0.5: no signal left inside its own band)")
        print(f"  AUC(-eff)     = {a_eff:.4f}   (control: is this just re-reading eff?)")
        print(f"  Spearman(NU, mu) inside band = "
              f"{spearman([r['NU'] for r in band], [r['mu'] for r in band]):.4f}")
    else:
        print("  degenerate band (all one class) -- cannot compute AUC")

    print("\nper-eff-stratum breakdown inside the band:")
    print(f"{'eff':>5} {'N':>4} {'rec':>7} {'AUC mu':>8} {'AUC NU':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p_ = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not sub:
            continue
        if p_ and ng:
            print(f"{eff:>5.2f} {len(sub):>4} {str(len(p_))+'/'+str(len(sub)):>7} "
                  f"{auc([r['mu'] for r in p_],[r['mu'] for r in ng]):>8.4f} "
                  f"{auc([r['NU'] for r in p_],[r['NU'] for r in ng]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>4} {str(len(p_))+'/'+str(len(sub)):>7} "
                  f"{'(degen)':>8} {'(degen)':>8}")

    # -----------------------------------------------------------------------
    # Part B. Secondary: the "step" statistic on the U2 K1-grid (Thread 24).
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("W1b follow-up: step = log2(||b*_{m+1}||) - log2(||b*_1||) vs the wall")
    print("-" * 78)

    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, label
        hist.append((label, (p, b, n, lam, G), k1, m))

    U2 = [("12-bit/2557", hist[1][1], 8, 0.340),
          ("12-bit/2677", hist[2][1], 10, 0.070)]
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

    u2rows = []
    for label, curve, m, ls in U2:
        n = curve[2]
        for k1 in K1_GRID:
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance(curve, m, d_trial, k1, seed, exact=True)
                if r is None:
                    continue
                rk = run_new(curve, m, d_trial, k1, seed)
                r.update({'label': label, 'm': m, 'n': n, 'K1': k1,
                          'seed': seed, 'ok': bool(rk['ok']),
                          'eff': k1 * (math.isqrt(n) + 1) / n})
                r['step'] = log2_step(r)
                u2rows.append(r)

    pos = [r for r in u2rows if r['ok']]
    neg = [r for r in u2rows if not r['ok']]
    print(f"{len(u2rows)} instances (exact GS)  "
          f"({len(pos)} recovered / {len(neg)} failed)\n")
    if pos and neg:
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"pooled (N={len(u2rows)}):")
        print(f"  AUC(-step) = {a_step:.4f}   (step -> 0 predicted to favor recovery)")
        print(f"  AUC(-NU)   = {a_nu:.4f}")
        print(f"  AUC(-mu)   = {a_mu:.4f}")
        print("\nper-curve (does pooling two different step baselines invert the sign?):")
        for label, curve, m, ls in U2:
            sub = [r for r in u2rows if r['label'] == label]
            p_ = [r for r in sub if r['ok']]
            ng = [r for r in sub if not r['ok']]
            if p_ and ng:
                print(f"  {label:>14}: AUC(-step) = "
                      f"{auc([r['step'] for r in p_],[r['step'] for r in ng]):.4f}   "
                      f"AUC(-NU) = {auc([r['NU'] for r in p_],[r['NU'] for r in ng]):.4f}")

    print(f"\n{'label':>14} {'K1':>4} {'eff':>6} {'step':>8} {'mean NU':>8} {'rec':>6}")
    bykey = {}
    for r in u2rows:
        bykey.setdefault((r['label'], r['K1']), []).append(r)
    for (label, k1), g in sorted(bykey.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        print(f"{label:>14} {k1:>4} {g[0]['eff']:>6.3f} "
              f"{sum(x['step'] for x in g)/len(g):>8.3f} "
              f"{sum(x['NU'] for x in g)/len(g):>8.4f} "
              f"{str(sum(1 for x in g if x['ok']))+'/'+str(len(g)):>6}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
