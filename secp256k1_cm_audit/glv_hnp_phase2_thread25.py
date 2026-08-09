"""
GLV-HNP Phase 2, Thread 25: does mu carry information NU does not, and does
the GS-profile step statistic beat both?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry.  W5/W6 of that run
found recovery = f(NU, X) with X ~ mu-driven and Spearman(mu-predictor, NU)
~ 0 inside every eff stratum -- i.e. NU (exact BDD certificate) and mu
(lambda_1(L2), a closed-form, lattice-reduction-free quantity) look like two
independent coordinates, not the same statistic wearing two hats.

  H25: within the ambiguous NU band [1.04, 2.20] (from the parent script's
       W4: sufficient NU < 1.040, necessary NU > 2.199, 17-bit/dim-24), where
       the exact nearest-plane certificate gives no verdict, AUC(-mu ->
       recovery) stays >= 0.8.
  Falsifier: if AUC inside the band drops to ~0.5, mu's apparent power in W5
       is entirely mediated by NU (a stratification artifact of eff), and
       the closed-form nu_hat/mu line of attack should be retired.

Secondary (also pre-registered): step = log2(prof[m]) - log2(prof[0]), the
gap between the m-fold-repeated head of the GS profile (all sitting at
lambda_1(L2), per Thread 24's W1b) and the first vector of the second block.
W1b showed this step visually vanishes right as recovery turns on; here it
is tested as a numeric separator, both pooled and inside the NU band.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json FILE]
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


def collect(m=12, effs=(0.05, 0.10, 0.15, 0.20, 0.25)):
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    rows = []
    for eff in effs:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m, d_trial, k1b, seed)
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                rows.append({
                    'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                    'eff': k1b * k2b / n, 'effq': eff,
                    'lamstar': lam_star(lam, n),
                    'mu': r['mu'], 'NU': r['NU'], 'nuhat': r['nuhat'],
                    'step': step, 'seed': seed,
                })
    return rows, len(curves17)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the collected row table to this path")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 - conditioning on NU: does mu (and the GS step) carry")
    print("independent separating power inside the ambiguous NU band?")
    print("=" * 78)

    M17 = 12
    t0 = time.time()
    rows, ncurves = collect(m=M17)
    print(f"\n{ncurves} 17-bit j=0 GLV curves; {len(rows)} instances "
          f"(float GS, dim {2*M17}) in {time.time()-t0:.1f}s")

    if args.dump_json:
        with open(args.dump_json, "w") as f:
            json.dump(rows, f)
        print(f"wrote {len(rows)} rows to {args.dump_json}")

    # W4 band from the parent (gsprofile_strat) 17-bit run: sufficient
    # NU < 1.040, necessary NU > 2.199.
    BAND_LO, BAND_HI = 1.040, 2.199

    print("\n" + "-" * 78)
    print(f"EXP H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{BAND_LO}, {BAND_HI}]")
    print("-" * 78)
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band population: {len(band)} / {len(rows)} instances "
          f"({len(pos)} recovered, {len(neg)} failed)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu)      = {a_mu:.4f}   (H25 threshold: >= 0.80)")
        print(f"  AUC(-nu_hat)  = {a_nh:.4f}")
        print(f"  AUC(-NU)      = {a_nu:.4f}   (sanity: should be ~0.5, band"
              f" is where NU is uninformative by construction)")
        print(f"  AUC(step)     = {a_st:.4f}   (step -> smaller predicts"
              f" recovery if this direction is right)")
        verdict = "HOLDS" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict} (AUC(-mu) = {a_mu:.4f})")
    else:
        print("degenerate band (one class empty) -- cannot compute AUC; "
              "H25 untestable at this sample size")

    print("\n" + "-" * 78)
    print("EXP secondary: step = log2(prof[m]) - log2(prof[0]) as a")
    print("standalone separator, pooled and per eff-stratum")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} {'AUC mu':>8} "
          f"{'AUC NU':>8}")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in rows if r['effq'] == eff and not math.isnan(r['step'])]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        a_st = auc([r['step'] for r in p], [r['step'] for r in ng])
        a_mu = auc([r['mu'] for r in p], [r['mu'] for r in ng])
        a_nu = auc([r['NU'] for r in p], [r['NU'] for r in ng])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(p))+'/'+str(len(sub)):>7} | {a_st:>9.4f} "
              f"{a_mu:>8.4f} {a_nu:>8.4f}")

    valid = [r for r in rows if not math.isnan(r['step'])]
    pooled_pos = [r for r in valid if r['ok']]
    pooled_neg = [r for r in valid if not r['ok']]
    if pooled_pos and pooled_neg:
        print(f"\npooled (N={len(valid)}): AUC(step -> recovery) = "
              f"{auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg]):.4f}")
        print(f"Spearman(step, NU) = "
              f"{spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
        print(f"Spearman(step, mu) = "
              f"{spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
