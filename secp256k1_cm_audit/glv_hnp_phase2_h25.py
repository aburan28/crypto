"""
GLV-HNP Phase 2, Thread 25: is mu a genuine SECOND coordinate once NU is
held fixed, or is its apparent power in W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (W4's bracket, where
       nearest-plane gives no answer either way), AUC(-mu -> Kannan-LLL
       recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate and (NU, mu) is a 2-parameter
    viability test.
  If no: mu's apparent power in W5 is entirely mediated by NU, and the W5
    stratification (eff held fixed) was not enough to isolate it.

Secondary (W1b): the GS profile of L0 (dim 2m) is two blocks -- the first m
entries pinned at lambda_1(L2), the second m moving with K1. Define
    step = log2(||b*_{m}||) - log2(||b*_1||)     (0-indexed: prof[m]-prof[0])
and test whether step alone separates recovery as well as NU or mu.

Same 17-bit, 5-strata x 20-curve x 5-seed dataset as
glv_hnp_phase2_gsprofile_strat.py (W5/W6), regenerated here (deterministic:
search_curves has no unseeded randomness) so 'prof' is retained per-row --
the strat script discards it. Float GS, justified by W0 of the parent script
(rel. error ~1e-15 at this dimension).

Run: python3 glv_hnp_phase2_h25.py [--dump-json out.json]
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

NU_LO, NU_HI = 1.040, 2.199  # W4's empirical bracket at 17 bits (log line ~6275)


def build_dataset():
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
                m = r['k'] // 2
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else 0.0)
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step,
                          'seed': seed})
                rows.append(r)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the row table (minus 'prof'/'nus' arrays) to this path")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — does mu survive as a 2nd coordinate once NU is fixed?")
    print("=" * 78)

    t0 = time.time()
    rows = build_dataset()
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        slim = [{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                for r in rows]
        with open(args.dump_json, "w") as f:
            json.dump(slim, f)
        print(f"wrote {len(slim)} rows to {args.dump_json}")

    print("\n" + "-" * 78)
    print(f"EXP H25: stratify by NU band. Ambiguous band = [{NU_LO}, {NU_HI}]")
    print("-" * 78)

    below = [r for r in rows if r['NU'] < NU_LO]
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    above = [r for r in rows if r['NU'] > NU_HI]
    for name, sub in (("NU < lo (sufficient)", below),
                       ("lo <= NU <= hi (ambiguous)", band),
                       ("NU > hi (necessary-fail)", above)):
        if not sub:
            print(f"{name:<28} N=0")
            continue
        rec = sum(1 for r in sub if r['ok'])
        print(f"{name:<28} N={len(sub):>4}  recovery {rec}/{len(sub)} "
              f"({100.0*rec/len(sub):.1f}%)")

    print(f"\nSanity: bracket transported from the 12-bit/17-bit run "
          f"(different lattice/dataset instantiation than W4 -- NU should\n"
          f"still be near-sound at the edges, but is not guaranteed to be "
          f"identical since W4 used the 22-cell K1 grid, not this eff grid).")

    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    if pos_b and neg_b:
        a_mu = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        a_step = auc([r['step'] for r in pos_b], [r['step'] for r in neg_b])
        a_nuhat = auc([r['nuhat'] for r in pos_b], [r['nuhat'] for r in neg_b])
        print(f"\nWithin ambiguous band (N={len(band)}, "
              f"{len(pos_b)} recover / {len(neg_b)} fail):")
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}   "
              f"(H25 predicts >= 0.80)")
        print(f"  AUC(-step   -> recovery) = {a_step:.4f}")
        print(f"  AUC(-nu_hat -> recovery) = {a_nuhat:.4f}")
        verdict = "SURVIVES" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict} (AUC(-mu) = {a_mu:.4f} "
              f"{'>=' if a_mu >= 0.80 else '<'} 0.80)")
    else:
        print(f"\nband is degenerate (pos={len(pos_b)}, neg={len(neg_b)}); "
              f"H25 cannot be tested against this bracket.")

    print("\n" + "-" * 78)
    print("EXP H25b: same test, but per-eff-stratum WITHIN the band -- guards")
    print("against band membership itself being confounded with eff")
    print("-" * 78)
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"  eff={eff:.2f}  N={len(sub):>4}  degenerate "
                  f"(pos={len(pos)}, neg={len(neg)})")
            continue
        a_mu_s = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh_s = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        print(f"  eff={eff:.2f}  N={len(sub):>4}  rec {len(pos)}/{len(sub)}  "
              f"AUC(-mu)={a_mu_s:.4f}  AUC(-nu_hat)={a_nh_s:.4f}")

    print("\n" + "-" * 78)
    print("EXP H25-secondary: step = log2||b*_{m}|| - log2||b*_1|| as a")
    print("standalone predictor, pooled over all strata")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    a_step_all = auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg])
    a_NU_all = auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg])
    a_mu_all = auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg])
    print(f"pooled (N={len(rows)}):")
    print(f"  AUC(-step -> recovery) = {a_step_all:.4f}")
    print(f"  AUC(-NU   -> recovery) = {a_NU_all:.4f}")
    print(f"  AUC(-mu   -> recovery) = {a_mu_all:.4f}")

    print("\nper-eff-stratum AUC(-step):")
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"  eff={eff:.2f}  degenerate")
            continue
        a = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  eff={eff:.2f}  N={len(sub):>4}  AUC(-step) = {a:.4f}")

    print("\nSpearman(step, NU) pooled: "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print("Spearman(step, mu) pooled: "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
