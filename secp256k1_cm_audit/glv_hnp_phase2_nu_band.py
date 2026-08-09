"""
GLV-HNP Phase 2, Thread 25: does mu survive as a SECOND coordinate once NU is
held fixed, and does the GS-profile "step" statistic explain the wall better
than either?

Pre-registered by the 2026-08-07 (autolab run #2) log entry (Thread 24, W5/W6):
  W5/W6 showed recovery = f(NU, X) with X ~ mu-driven and X uncorrelated with
  NU at fixed eff.  H25 asks whether that holds the OTHER way: condition on
  NU instead of eff, and see if mu still separates inside the ambiguous band.

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (17-bit bracket from
       Thread 23b/24 W4, where nearest-plane's own certificate gives no
       answer either way), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  Falsifier: if AUC inside the band drops to ~0.5, mu's apparent power in W5
  was entirely mediated by eff/NU and the closed form should be retired as a
  genuine second coordinate (it may still be useful as a cheap NU proxy).

Secondary (also pre-registered): step = log2(||b*_{m+1}||) - log2(||b*_1||),
the point where the GS profile of L0 transitions from the m exact copies of
lambda_1(L2) (Thread 24 W1b) to the second block.  Test whether step
predicts the wall better than NU or mu, in- and out-of-band.

Data: reuses the EXACT data-collection loop of glv_hnp_phase2_gsprofile_strat.py
(5 eff strata x 20 17-bit curves x 5 seeds, dim 24, float GS -- justified by
W0/W4: max relative NU error ~1e-15 at dim 24) so this is a re-analysis of the
same population, not a new experiment population.

Run: python3 glv_hnp_phase2_nu_band.py [--dump-json OUT.json] [--load-json IN.json]
"""

import argparse
import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

BAND_LO, BAND_HI = 1.040, 2.199  # 17-bit NU bracket, Thread 24 W4
M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)


def collect():
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
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
                m = M17
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                rows.append({
                    'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                    'eff': k1b * k2b / n, 'effq': eff,
                    'lamstar': lam_star(lam, n),
                    'NU': r['NU'], 'mu': r['mu'], 'nuhat': r['nuhat'],
                    'step': step,
                })
    return rows


def load(path):
    with open(path) as f:
        return json.load(f)


def dump(rows, path):
    with open(path, 'w') as f:
        json.dump(rows, f)


def band_report(rows, lo, hi, label):
    band = [r for r in rows if lo <= r['NU'] <= hi]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"\n{label}: N={len(band)}  ({len(pos)} recover / {len(neg)} fail)")
    if not pos or not neg:
        print("  degenerate (one class empty) -- no AUC")
        return None
    a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
    a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
    a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
    a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
    print(f"  AUC(-mu)     = {a_mu:.4f}")
    print(f"  AUC(-nu_hat) = {a_nh:.4f}")
    print(f"  AUC(-NU)     = {a_nu:.4f}   (should be ~0.5: NU is constant-ish in-band)")
    print(f"  AUC(-step)   = {a_st:.4f}")
    return {'N': len(band), 'auc_mu': a_mu, 'auc_nuhat': a_nh,
            'auc_NU': a_nu, 'auc_step': a_st}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--dump-json')
    ap.add_argument('--load-json')
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 -- does mu survive as a 2nd coordinate once NU is fixed?")
    print("=" * 78)

    if args.load_json:
        rows = load(args.load_json)
        print(f"\nloaded {len(rows)} instances from {args.load_json}")
    else:
        t0 = time.time()
        rows = collect()
        print(f"\ncollected {len(rows)} instances (float GS, 17-bit, dim 24) "
              f"in {time.time()-t0:.1f}s")
        if args.dump_json:
            dump(rows, args.dump_json)
            print(f"dumped to {args.dump_json}")

    print("\n" + "-" * 78)
    print(f"H25: within the ambiguous NU band [{BAND_LO}, {BAND_HI}], "
          "does AUC(-mu) stay >= 0.8?")
    print("-" * 78)
    r_band = band_report(rows, BAND_LO, BAND_HI, "IN-BAND")
    r_lo = band_report(rows, 0.0, BAND_LO, "BELOW BAND (NU < 1.04, mostly recovers)")
    r_hi = band_report(rows, BAND_HI, 1e9, "ABOVE BAND (NU > 2.20, mostly fails)")

    print("\n" + "-" * 78)
    print("Pooled (no NU conditioning) for reference:")
    print("-" * 78)
    pos = [r for r in rows if r['ok']]
    neg = [r for r in rows if not r['ok']]
    print(f"  AUC(-mu)     = {auc([r['mu'] for r in pos], [r['mu'] for r in neg]):.4f}")
    print(f"  AUC(-nu_hat) = {auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg]):.4f}")
    print(f"  AUC(-NU)     = {auc([r['NU'] for r in pos], [r['NU'] for r in neg]):.4f}")
    print(f"  AUC(-step)   = {auc([r['step'] for r in pos], [r['step'] for r in neg]):.4f}")

    print("\n" + "-" * 78)
    print("Verdict")
    print("-" * 78)
    if r_band and r_band['auc_mu'] >= 0.8:
        print(f"H25 HOLDS: AUC(-mu) = {r_band['auc_mu']:.4f} >= 0.8 in-band. "
              "mu is a genuine second coordinate.")
    elif r_band:
        print(f"H25 FALSIFIED: AUC(-mu) = {r_band['auc_mu']:.4f} < 0.8 in-band. "
              "mu's apparent power was mediated by NU/eff after all.")
    else:
        print("H25 UNTESTABLE: in-band sample is degenerate (one class empty).")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
