"""
GLV-HNP Phase 2, Thread 24b: is the closed-form separator actually doing
cross-curve work, or is it just re-reading the bias strength?

W3 (glv_hnp_phase2_gsprofile.py) found AUC(-nu_hat*sqrt(eff) -> recovery)
= 0.992 on the 22-cell U2 grid, beating the exact BDD certificate NU (0.978).
That grid varies K1 over 11 values on only 2 curves, and recovery is monotone
in K1 within a curve, so a predictor monotone in K1 scores high almost for
free.  The honest test holds eff FIXED and asks whether nu_hat still ranks
curves correctly.

W5  per-eff-stratum AUC at 17 bits (20 curves x 5 seeds per stratum), for
    three scores: nu_hat (no lattice work), NU (exact BDD certificate),
    and lam* (the quantity Thread 20 falsified).
W6  is C = NU/(nu_hat*sqrt(eff)) stable across strata, i.e. does the closed
    form carry an absolute scale or only a ranking?

Gram-Schmidt is float here, justified by W0/W4 of the parent script
(max relative NU error vs exact Fractions ~1e-15 at dim 20 and dim 24).

Run: python3 glv_hnp_phase2_gsprofile_strat.py
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

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the raw per-instance rows (incl. the "
                          "per-index GS profile 'prof' and nu-list 'nus') "
                          "to this path so later runs can re-analyze "
                          "without regenerating the 500 instances")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 24b — cross-curve test of the closed-form separator (eff fixed)")
    print("=" * 78)

    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
    M17 = 12
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

    t0 = time.time()
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
                          'lamstar': lam_star(lam, n)})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        with open(args.dump_json, "w") as f:
            json.dump(rows, f)
        print(f"wrote {len(rows)} rows to {args.dump_json}")

    print("\n" + "-" * 78)
    print("EXP W5: AUC within each eff stratum — eff is CONSTANT, so the only")
    print("        signal left is the cross-curve geometry.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC nu_hat':>11} {'AUC NU':>8} "
          f"{'AUC lam*':>9} | {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | "
                  f"{'(degenerate)':>11}")
            continue
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_nh:>11.4f} "
              f"{a_nu:>8.4f} {a_ls:>9.4f} | {a_mu:>8.4f}")

    print("\nAUC > 0.5 means SMALLER score -> more likely to recover.")
    print("Thread 20 falsified lam* as a predictor; it is the control column.")

    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    print(f"\npooled over all strata (N={len(rows)}):")
    print(f"  AUC(-nu_hat*sqrt(eff)) = "
          f"{auc([r['nuhat']*math.sqrt(r['eff']) for r in pooled_pos], [r['nuhat']*math.sqrt(r['eff']) for r in pooled_neg]):.4f}")
    print(f"  AUC(-NU)               = "
          f"{auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-nu_hat alone)     = "
          f"{auc([r['nuhat'] for r in pooled_pos], [r['nuhat'] for r in pooled_neg]):.4f}")
    print(f"  AUC(-eff alone)        = "
          f"{auc([r['eff'] for r in pooled_pos], [r['eff'] for r in pooled_neg]):.4f}")

    print("\n" + "-" * 78)
    print("EXP W6: is C = NU / (nu_hat*sqrt(eff)) an absolute constant?")
    print("-" * 78)
    print(f"{'eff':>5} {'mean C':>9} {'min':>8} {'max':>8} {'spread':>8} "
          f"{'Spearman(pred,NU)':>19}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        cs = [r['NU'] / (r['nuhat'] * math.sqrt(r['eff'])) for r in sub
              if r['nuhat'] > 0]
        pr = [r['nuhat'] * math.sqrt(r['eff']) for r in sub]
        nu = [r['NU'] for r in sub]
        print(f"{eff:>5.2f} {sum(cs)/len(cs):>9.3f} {min(cs):>8.3f} "
              f"{max(cs):>8.3f} {max(cs)/min(cs):>8.2f}x "
              f"{spearman(pr, nu):>19.4f}")

    print("\n" + "-" * 78)
    print("EXP W7: per-curve detail at the discriminating stratum")
    print("-" * 78)
    disc = None
    best = -1
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        w = sum(1 for r in sub if r['ok'])
        bal = min(w, len(sub) - w)
        if bal > best:
            best, disc = bal, eff
    sub = [r for r in rows if r['effq'] == disc]
    print(f"stratum eff = {disc:.2f}  ({best} of the minority class)\n")
    print(f"{'n':>8} {'lam*':>7} {'nu_hat':>8} {'mean NU':>9} {'rec':>6}")
    bycurve = {}
    for r in sub:
        bycurve.setdefault(r['n'], []).append(r)
    for n in sorted(bycurve, key=lambda x: bycurve[x][0]['nuhat']):
        g = bycurve[n]
        print(f"{n:>8} {g[0]['lamstar']:>7.4f} {g[0]['nuhat']:>8.4f} "
              f"{sum(x['NU'] for x in g)/len(g):>9.4f} "
              f"{str(sum(1 for x in g if x['ok']))+'/'+str(len(g)):>6}")

    print("\n" + "-" * 78)
    print("Thread 25 — H25: does mu separate INSIDE the NU-ambiguous band?")
    print("-" * 78)
    print("Band from the 2026-08-07 log (W4, 17 bits): sufficient NU < 1.040,")
    print("necessary NU > 2.199 -> nearest-plane gives no answer in between.")
    NU_LO, NU_HI = 1.040, 2.199
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    print(f"band size: {len(band)} / {len(rows)}  "
          f"({len(pos_b)} recovered, {len(neg_b)} not)")
    if pos_b and neg_b:
        a_mu_band = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        a_nh_band = auc([r['nuhat'] for r in pos_b], [r['nuhat'] for r in neg_b])
        print(f"  AUC(-mu     -> recovery | band) = {a_mu_band:.4f}")
        print(f"  AUC(-nu_hat -> recovery | band) = {a_nh_band:.4f}")
        verdict = "H25 HOLDS" if a_mu_band >= 0.8 else "H25 FAILS"
        print(f"  {verdict}: threshold is AUC >= 0.8000")
    else:
        print("  degenerate: band is all-positive or all-negative, no AUC defined")

    print("\n" + "-" * 78)
    print("Thread 25 secondary — W1b step statistic: does step -> 0 predict")
    print("the wall better than NU or mu?")
    print("-" * 78)
    print("step := log2(||b*_{m+1}||) - log2(||b*_1||) = log2(prof[M]) - log2(prof[0])")
    step_rows = [r for r in rows if r['prof'][0] > 0 and r['prof'][M17] > 0]
    for r in step_rows:
        r['step'] = math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
    pos_s = [r for r in step_rows if r['ok']]
    neg_s = [r for r in step_rows if not r['ok']]
    if pos_s and neg_s:
        a_step = auc([r['step'] for r in pos_s], [r['step'] for r in neg_s])
        a_nu_s = auc([r['NU'] for r in pos_s], [r['NU'] for r in neg_s])
        a_mu_s = auc([r['mu'] for r in pos_s], [r['mu'] for r in neg_s])
        print(f"  N = {len(step_rows)} ({len(pos_s)} recovered)")
        print(f"  AUC(-step -> recovery) = {a_step:.4f}   "
              f"(step -> 0 predicts recovery if this is >> 0.5)")
        print(f"  AUC(-NU   -> recovery) = {a_nu_s:.4f}  (same pool, for comparison)")
        print(f"  AUC(-mu   -> recovery) = {a_mu_s:.4f}  (same pool, for comparison)")
        print(f"  Spearman(step, NU)     = {spearman([r['step'] for r in step_rows], [r['NU'] for r in step_rows]):.4f}")
        print(f"  Spearman(step, mu)     = {spearman([r['step'] for r in step_rows], [r['mu'] for r in step_rows]):.4f}")
    else:
        print("  degenerate pool, skipped")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
