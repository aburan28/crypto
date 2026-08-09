"""
GLV-HNP Phase 2, Thread 25: is mu a genuine second coordinate, or is its
power in W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous NU band [1.04, 2.20] (17 bits, where nearest-
       plane gives no answer either way), AUC(-mu -> Kannan-LLL recovery)
       stays >= 0.8.

  If yes: mu is a genuine second coordinate independent of NU, and (NU, mu)
  is a 2-parameter viability test.
  If no: mu's apparent power in W5 is entirely mediated by NU (curves with
  small mu also tend to have small NU), and the W5 result is a
  stratification artifact — the closed form should be retired as a
  predictor and NU treated as sufficient statistic.

Secondary (W1b follow-up): does the profile head/tail step,

    step = log2(||b*_{m+1}||) - log2(||b*_1||)

(vanishing exactly as the K1 wall is crossed, per Thread 24 W1b) predict
recovery better than NU or mu?  One-line addition to the existing sweep.

Data: same 500-instance 17-bit generation as glv_hnp_phase2_gsprofile_strat.py
(20 curves x 5 eff strata x 5 seeds, dim 24, float GS -- justified by W0/W4
of glv_hnp_phase2_gsprofile.py, max relative NU error ~1e-15).

Run: python3 glv_hnp_phase2_thread25.py [--dump-json path]
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


def collect(curves17, m):
    effs = (0.05, 0.10, 0.15, 0.20, 0.25)
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
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                # prof is a list of floats; keep it out of the json dump by
                # default (it is large and reconstructible from seed+curve).
                rows.append(r)
    return rows


def json_safe(row):
    return {k: v for k, v in row.items() if k not in ('prof', 'nus')}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the raw row table to this path")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — does mu survive conditioning on NU?  (H25)")
    print("=" * 78)

    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
    M17 = 12

    t0 = time.time()
    rows = collect(curves17, M17)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        with open(args.dump_json, "w") as f:
            json.dump([json_safe(r) for r in rows], f)
        print(f"dumped {len(rows)} rows to {args.dump_json}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25a: recompute the pooled NU ambiguous band on THIS run")
    print("-" * 78)
    pos_all = [r['NU'] for r in rows if r['ok']]
    neg_all = [r['NU'] for r in rows if not r['ok']]
    band_lo = min(neg_all)   # sufficient: NU < band_lo -> always recovers
    band_hi = max(pos_all)   # necessary:  NU > band_hi -> never recovers
    print(f"pooled (N={len(rows)}): sufficient NU < {band_lo:.4f} , "
          f"necessary NU > {band_hi:.4f}")
    print("(compare 2026-08-07 log: sufficient 1.040, necessary 2.199 -- "
          "this run redraws 500 fresh instances so the band will differ "
          "slightly by sampling noise)")

    band = [r for r in rows if band_lo <= r['NU'] <= band_hi]
    band_pos = [r for r in band if r['ok']]
    band_neg = [r for r in band if not r['ok']]
    print(f"\nband population: {len(band)}/{len(rows)}  "
          f"({len(band_pos)} recover, {len(band_neg)} fail)")

    if band_pos and band_neg:
        auc_mu_band = auc([r['mu'] for r in band_pos], [r['mu'] for r in band_neg])
        auc_nuhat_band = auc([r['nuhat'] for r in band_pos], [r['nuhat'] for r in band_neg])
        auc_NU_band = auc([r['NU'] for r in band_pos], [r['NU'] for r in band_neg])
        print(f"\nWITHIN the ambiguous band, does NU itself still separate?")
        print(f"  AUC(-NU     -> recovery | in band) = {auc_NU_band:.4f}  "
              f"(should be ~0.5 by construction -- band is where NU fails)")
        print(f"  AUC(-mu     -> recovery | in band) = {auc_mu_band:.4f}")
        print(f"  AUC(-nu_hat -> recovery | in band) = {auc_nuhat_band:.4f}")
        print(f"\nH25 verdict: AUC(-mu | in band) = {auc_mu_band:.4f}  "
              f"{'>= 0.8 -> H25 HOLDS, mu is a genuine 2nd coordinate' if auc_mu_band >= 0.8 else '< 0.8 -> H25 FAILS'}")
    else:
        print("band is degenerate (all one class) -- cannot test H25 here")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25b: per-eff-stratum band test (does band population even")
    print("          exist inside single strata, or only pooled?)")
    print("-" * 78)
    print(f"{'eff':>5} {'band N':>7} {'rec':>7} | {'AUC mu|band':>12} "
          f"{'AUC NU|band':>12}")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in rows if r['effq'] == eff and band_lo <= r['NU'] <= band_hi]
        sp = [r for r in sub if r['ok']]
        sn = [r for r in sub if not r['ok']]
        if not sp or not sn:
            print(f"{eff:>5.2f} {len(sub):>7} {'(degenerate)':>7}")
            continue
        a_mu = auc([r['mu'] for r in sp], [r['mu'] for r in sn])
        a_nu = auc([r['NU'] for r in sp], [r['NU'] for r in sn])
        print(f"{eff:>5.2f} {len(sub):>7} "
              f"{str(len(sp))+'/'+str(len(sub)):>7} | {a_mu:>12.4f} "
              f"{a_nu:>12.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25c (secondary): does the head/tail GS step predict better")
    print("          than NU or mu?  step = log2||b*_{m+1}|| - log2||b*_1||")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    print(f"valid step values: {len(valid)}/{len(rows)}")
    pos = [r['step'] for r in valid if r['ok']]
    neg = [r['step'] for r in valid if not r['ok']]
    if pos and neg:
        a_step_pos_smaller = auc(pos, neg)
        a_step_neg_smaller = auc(neg, pos)
        print(f"step | success : mean {sum(pos)/len(pos):.3f}  "
              f"min {min(pos):.3f}  max {max(pos):.3f}")
        print(f"step | failure : mean {sum(neg)/len(neg):.3f}  "
              f"min {min(neg):.3f}  max {max(neg):.3f}")
        print(f"AUC(-step -> recovery) = {a_step_pos_smaller:.4f}   "
              f"AUC(+step -> recovery) = {a_step_neg_smaller:.4f}")
        best = max(a_step_pos_smaller, a_step_neg_smaller)
        print(f"best-direction AUC(step) = {best:.4f}   "
              f"vs pooled AUC(NU) = {auc(pos_all, neg_all):.4f}   "
              f"vs pooled AUC(mu) = "
              f"{auc([r['mu'] for r in rows if r['ok']], [r['mu'] for r in rows if not r['ok']]):.4f}")

    print("\nper-eff-stratum step AUC (best direction):")
    print(f"{'eff':>5} {'N':>5} {'AUC step':>9}")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in valid if r['effq'] == eff]
        sp = [r['step'] for r in sub if r['ok']]
        sn = [r['step'] for r in sub if not r['ok']]
        if not sp or not sn:
            print(f"{eff:>5.2f} {len(sub):>5} {'(degenerate)':>9}")
            continue
        best = max(auc(sp, sn), auc(sn, sp))
        print(f"{eff:>5.2f} {len(sub):>5} {best:>9.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
