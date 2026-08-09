"""
GLV-HNP Phase 2, Thread 25: does mu separate INSIDE the NU ambiguous band?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  W5/W6 established recovery = f(NU, X), X ~ mu-driven and independent of NU
  (Spearman(pred, NU) ~ 0 or negative in every eff stratum; W6).  NU is a
  sound but loose certificate (TP 71 / FP 0 at NU<=1, but ambiguous band
  [1.040, 2.199] holds most of the mass at 17 bits -- Thread 24 W4).

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate and (NU, mu) is a 2-parameter
  viability test.  If no: mu's apparent power (W5) is entirely mediated by
  NU and the closed form should be retired.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||) on the
Kannan-free L0.  Test whether step -> 0 predicts the K1 wall better than NU
or mu (one-line addition to the existing GS-profile machinery).

Data: reuses glv_hnp_phase2_gsprofile_strat.py's exact generation path
(same search_curves/instance/run_new calls, same SEEDS, same M17=12,
same EFFS) so this is a re-analysis, not new data -- search_curves is a
deterministic prime sweep (glv_hnp_common.py:306) and SEEDS/d_trial draws
are fixed, so the row set reproduces byte-for-byte.  --dump-json writes the
row table so future threads can re-slice without re-running the sweep.

Run: python3 glv_hnp_phase2_h25_nuband.py [--dump-json path]
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

# NU ambiguous band from Thread 24 W4 (17 bits, dim 24): sufficient < 1.040,
# necessary > 2.199.
NU_LO, NU_HI = 1.040, 2.199


def build_rows():
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
                # secondary: GS-profile step at the K1/K2 block boundary.
                # prof[i] = ||b*_i|| (not logged; instance() in
                # glv_hnp_phase2_gsprofile.py:153/156 returns raw sqrt-norms).
                # dim is 2m (k1-block = indices [0,m), k2-block = [m,2m)).
                m = M17
                if len(r['prof']) > m and r['prof'][0] > 0 and r['prof'][m] > 0:
                    step = math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                else:
                    step = None
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None)
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — does mu separate INSIDE the NU ambiguous band?")
    print("=" * 78)

    t0 = time.time()
    rows = build_rows()
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        # 'prof'/'nus' are per-index lists; keep them, everything else is
        # already JSON-safe (float/int/bool).
        with open(args.dump_json, "w") as f:
            json.dump(rows, f)
        print(f"wrote {len(rows)} rows to {args.dump_json}")

    print("\n" + "-" * 78)
    print("EXP H25: stratify by NU band, AUC(-mu -> recovery) inside each band")
    print("-" * 78)
    bands = [
        ("NU < lo (sufficient)", lambda r: r['NU'] < NU_LO),
        ("lo <= NU <= hi (AMBIGUOUS)", lambda r: NU_LO <= r['NU'] <= NU_HI),
        ("NU > hi (necessary-fail)", lambda r: r['NU'] > NU_HI),
    ]
    print(f"band = [{NU_LO}, {NU_HI}] (Thread 24 W4, 17-bit dim-24 bracket)\n")
    print(f"{'band':<28} {'N':>5} {'rec':>9} | {'AUC mu':>8} {'AUC nuhat':>10} "
          f"{'AUC step':>9}")
    band_result = {}
    for name, pred in bands:
        sub = [r for r in rows if pred(r)]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not sub:
            print(f"{name:<28} {'0':>5}   (empty)")
            continue
        rec = f"{len(pos)}/{len(sub)}"
        if not pos or not neg:
            print(f"{name:<28} {len(sub):>5} {rec:>9} | (degenerate: "
                  f"single class)")
            band_result[name] = None
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        steps_pos = [r['step'] for r in pos if r['step'] is not None]
        steps_neg = [r['step'] for r in neg if r['step'] is not None]
        a_st = (auc(steps_pos, steps_neg)
                if steps_pos and steps_neg else float('nan'))
        print(f"{name:<28} {len(sub):>5} {rec:>9} | {a_mu:>8.4f} "
              f"{a_nh:>10.4f} {a_st:>9.4f}")
        band_result[name] = a_mu

    amb = band_result.get("lo <= NU <= hi (AMBIGUOUS)")
    print()
    if amb is None:
        print("H25 UNTESTABLE: ambiguous band is single-class in this sample "
              "(need more instances or a wider band).")
    elif amb >= 0.8 or amb <= 0.2:
        # AUC<=0.2 means mu separates with the opposite-than-hypothesized
        # sign but still separates -- still a genuine second coordinate.
        print(f"H25 SURVIVES: AUC(-mu) = {amb:.4f} inside the ambiguous band "
              f"(threshold 0.8/0.2). mu carries information NU does not.")
    else:
        print(f"H25 FALSIFIED: AUC(-mu) = {amb:.4f} inside the ambiguous "
              f"band, below threshold. mu's W5 power looks mediated by NU.")

    print("\n" + "-" * 78)
    print("Control: same three-band split, but with lam* (Thread 20's "
          "falsified predictor)")
    print("-" * 78)
    for name, pred in bands:
        sub = [r for r in rows if pred(r)]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            continue
        a_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        print(f"{name:<28} AUC(-lam*) = {a_ls:.4f}")

    print("\n" + "-" * 78)
    print("EXP H25-secondary: does step -> 0 predict the wall pooled "
          "(all instances, not just the ambiguous band)?")
    print("-" * 78)
    pos_all = [r for r in rows if r['ok'] and r['step'] is not None]
    neg_all = [r for r in rows if not r['ok'] and r['step'] is not None]
    a_step_pool = auc([r['step'] for r in pos_all], [r['step'] for r in neg_all])
    a_mu_pool = auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all])
    a_nu_pool = auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all])
    print(f"pooled AUC(-step) = {a_step_pool:.4f}   "
          f"AUC(-mu) = {a_mu_pool:.4f}   AUC(-NU) = {a_nu_pool:.4f}")
    print(f"Spearman(step, mu)   = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")
    print(f"Spearman(step, NU)   = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
