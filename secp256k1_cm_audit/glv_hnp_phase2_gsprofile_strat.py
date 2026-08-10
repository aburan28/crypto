"""
GLV-HNP Phase 2, Thread 24b/25: is the closed-form separator actually doing
cross-curve work, or is it just re-reading the bias strength?  And once NU
(the exact BDD certificate) is ambiguous, does mu still carry information?

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

Thread 25 (2026-08-10) — pre-registered by the 2026-08-07 #2 log entry:

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (where the exact
       nearest-plane certificate gives no answer either way), does
       AUC(-mu -> Kannan-LLL recovery) stay >= 0.8?

  If yes, mu is a genuine second coordinate independent of NU and (NU, mu)
  is a 2-parameter viability test.  If no, mu's apparent power (W5) is
  entirely mediated by NU and the closed form should be retired.

  Secondary: step = log2(||b*_{m+1}||) - log2(||b*_1||), the jump from the
  flat head of the GS profile (m copies of lambda_1(L2), per Thread 24 W1b)
  into the second block.  Test whether step -> 0 predicts the wall better
  than NU or mu.

Gram-Schmidt is float here, justified by W0/W4 of the parent script
(max relative NU error vs exact Fractions ~1e-15 at dim 20 and dim 24).

Run: python3 glv_hnp_phase2_gsprofile_strat.py [--dump-json FILE]
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
                     help="write the collected instance table to this path")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 24b/25 — cross-curve test of the closed-form separator, "
          "and mu inside the NU-ambiguous band")
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
                prof = r['prof']
                step = (math.log2(prof[M17]) - math.log2(prof[0])
                        if prof[0] > 0 and prof[M17] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'm': M17, 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        dump = [{k: v for k, v in r.items() if k not in ('nus',)}
                for r in rows]
        with open(args.dump_json, "w") as f:
            json.dump(dump, f)
        print(f"[dumped {len(dump)} rows to {args.dump_json}]")

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
    print("EXP W8 (Thread 25, H25): inside the NU-ambiguous band, does mu")
    print("        still separate recovery from failure?")
    print("-" * 78)
    lo, hi = 1.040, 2.199  # bracket from Thread 24 EXP W4 (17-bit, this table's size)
    print(f"ambiguous band taken from the 17-bit bracket (Thread 24 EXP W4): "
          f"NU in [{lo:.3f}, {hi:.3f}]")
    band = [r for r in rows if lo <= r['NU'] <= hi]
    band_pos = [r for r in band if r['ok']]
    band_neg = [r for r in band if not r['ok']]
    print(f"band size N={len(band)}  ({len(band_pos)} recover / "
          f"{len(band_neg)} fail)")
    if band_pos and band_neg:
        a_mu_band = auc([r['mu'] for r in band_pos], [r['mu'] for r in band_neg])
        a_nh_band = auc([r['nuhat'] for r in band_pos],
                         [r['nuhat'] for r in band_neg])
        a_nu_band = auc([r['NU'] for r in band_pos], [r['NU'] for r in band_neg])
        print(f"  AUC(-mu     -> recovery | band) = {a_mu_band:.4f}")
        print(f"  AUC(-nu_hat -> recovery | band) = {a_nh_band:.4f}")
        print(f"  AUC(-NU     -> recovery | band) = {a_nu_band:.4f}  "
              f"(sanity: NU is nearly constant inside its own band)")
        verdict = "HOLDS" if a_mu_band >= 0.8 else "FAILS"
        print(f"\nH25 verdict (pooled): AUC(-mu | band) = {a_mu_band:.4f}  "
              f"({'>=' if a_mu_band >= 0.8 else '<'} 0.8)  -> H25 {verdict}")
        print("\nEXP W8b: is the pooled band confounded by eff?  Spearman(eff, NU | band) "
              f"= {spearman([r['eff'] for r in band], [r['NU'] for r in band]):.4f}")
        print("(NU still correlates with eff inside its own 'ambiguous' band -> "
              "pooling across eff strata is not a clean test of mu.)")
        print(f"\n{'eff':>5} {'N':>4} {'rec':>7} | {'AUC -mu':>8} {'AUC -nu_hat':>11}")
        band_strata_aucs = []
        for eff in EFFS:
            sub = [r for r in band if r['effq'] == eff]
            pos = [r for r in sub if r['ok']]
            neg = [r for r in sub if not r['ok']]
            if not pos or not neg:
                print(f"{eff:>5.2f} {len(sub):>4} "
                      f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
                continue
            a_mu_s = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
            a_nh_s = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
            band_strata_aucs.append(a_mu_s)
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu_s:>8.4f} "
                  f"{a_nh_s:>11.4f}")
        if band_strata_aucs:
            mean_strat = sum(band_strata_aucs) / len(band_strata_aucs)
            v2 = "HOLDS" if mean_strat >= 0.8 else "FAILS"
            print(f"\nH25 verdict (eff-controlled): mean per-stratum AUC(-mu | "
                  f"band, eff fixed) = {mean_strat:.4f}  -> H25 {v2}")
    else:
        print("band is degenerate (all-pos or all-neg) at this table's "
              "resolution — H25 untestable here.")

    print("\n" + "-" * 78)
    print("EXP W9 (Thread 25, secondary): step = log2||b*_{m+1}|| - "
          "log2||b*_1|| vs the wall")
    print("-" * 78)
    print("W1b (Thread 24) showed the GS-profile head is m exact copies of "
          "lambda_1(L2)\nand the step to the second block vanishes as the "
          "K1 wall is crossed.\nHypothesis: step -> 0 predicts failure at "
          "least as well as NU or mu.\n")
    step_pos = [r['step'] for r in rows if r['ok']]
    step_neg = [r['step'] for r in rows if not r['ok']]
    a_step_hi = auc(step_neg, step_pos)   # larger step -> recovery
    print(f"step | success : mean {sum(step_pos)/len(step_pos):.3f}  "
          f"min {min(step_pos):.3f}  max {max(step_pos):.3f}")
    print(f"step | failure : mean {sum(step_neg)/len(step_neg):.3f}  "
          f"min {min(step_neg):.3f}  max {max(step_neg):.3f}")
    print(f"\nAUC(+step -> recovery)  = {a_step_hi:.4f}   (larger step -> "
          f"more likely to recover)")
    a_nu_pool = auc([r['NU'] for r in rows if r['ok']],
                     [r['NU'] for r in rows if not r['ok']])
    a_mu_pool = auc([r['mu'] for r in rows if r['ok']],
                     [r['mu'] for r in rows if not r['ok']])
    print(f"compare pooled AUC(-NU -> recovery) = {a_nu_pool:.4f}, "
          f"AUC(-mu -> recovery) = {a_mu_pool:.4f}")
    if band_pos and band_neg:
        step_band_pos = [r['step'] for r in band_pos]
        step_band_neg = [r['step'] for r in band_neg]
        a_step_band = auc(step_band_neg, step_band_pos)
        print(f"AUC(+step -> recovery | NU-ambiguous band) = {a_step_band:.4f}")
    print(f"\nSpearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}   "
          f"Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
