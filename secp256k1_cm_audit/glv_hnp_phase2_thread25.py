"""
GLV-HNP Phase 2, Thread 25: is there a second mechanism beyond NU?

Pre-registered by the 2026-08-07 (Thread 24) log entry.  W5/W6 there showed
recovery = f(NU, X) with X ~ mu-driven (mu = lambda_1(L2)) and X uncorrelated
with NU at fixed eff.  This script runs the concrete falsifier:

  H25: within the ambiguous NU band [1.04, 2.20] (17 bits, ambiguous because
       nearest-plane gives no answer there per Thread 24's W4 bracket),
       AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate; (NU, mu) is a 2-parameter
  viability test.  If no: mu's apparent power in W5 was mediated by NU and
  the closed form should be retired as a stratification artifact.

Secondary (W1b follow-up): does the GS-profile "step" at the block boundary,

  step = log2(||b*_{m+1}||) - log2(||b*_1||)   (prof[m] - prof[0], log2)

predict recovery better than NU or mu?  One-line addition to the same sweep.

Uses the same 17-bit generation as glv_hnp_phase2_gsprofile_strat.py (float
GS, justified there by W0/W4: max relative NU error vs exact ~1e-15 at dim
24).  Adds --dump-json to persist the row table so re-analysis doesn't need
to regenerate data.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json FILE] [--load-json FILE]
"""

import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc

NU_BAND = (1.040, 2.199)  # Thread 24 W4: sufficient/necessary bracket at 17 bits


def collect_rows():
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"{len(curves17)} 17-bit j=0 GLV curves")
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
                          'eff': k1b * k2b / n, 'effq': eff})
                r['step'] = (math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                             if r['prof'][0] > 0 and r['prof'][M17] > 0 else float('nan'))
                del r['prof']  # not JSON-critical, keeps dump small; nus kept
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim 24) in {time.time()-t0:.1f}s")
    return rows


def main():
    dump_path = None
    load_path = None
    args = sys.argv[1:]
    while args:
        a = args.pop(0)
        if a == "--dump-json":
            dump_path = args.pop(0)
        elif a == "--load-json":
            load_path = args.pop(0)

    print("=" * 78)
    print("Thread 25 — second-mechanism test: does mu separate inside the")
    print("NU-ambiguous band?  (H25)")
    print("=" * 78)

    if load_path:
        with open(load_path) as f:
            rows = json.load(f)
        print(f"loaded {len(rows)} rows from {load_path}")
    else:
        rows = collect_rows()
        if dump_path:
            with open(dump_path, "w") as f:
                json.dump(rows, f)
            print(f"dumped {len(rows)} rows to {dump_path}")

    print("\n" + "-" * 78)
    print(f"H25: AUC(-mu -> recovery) inside NU in [{NU_BAND[0]}, {NU_BAND[1]}]")
    print("-" * 78)
    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    print(f"band N = {len(band)} / {len(rows)} total "
          f"({100.0*len(band)/len(rows):.1f}%)")
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band recovery: {len(pos)}/{len(band)}")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        print(f"AUC(-mu    -> recovery) inside band = {a_mu:.4f}   "
              f"(H25 threshold: >= 0.80)")
        print(f"AUC(-NU    -> recovery) inside band = {a_nu:.4f}   "
              f"(expect ~0.5: band is defined to be ambiguous in NU)")
        print(f"AUC(-step  -> recovery) inside band = {a_step:.4f}")
        print(f"AUC(-nuhat -> recovery) inside band = {a_nh:.4f}")
        verdict = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict}  (AUC={a_mu:.4f})")
    else:
        print("band has only one class -- AUC undefined, cannot test H25 "
              "(degenerate band)")

    print("\n" + "-" * 78)
    print("Pooled (unconditional) comparison: is step a better global "
          "predictor than NU or mu?")
    print("-" * 78)
    allpos = [r for r in rows if r['ok']]
    allneg = [r for r in rows if not r['ok']]
    print(f"AUC(-step -> recovery), all rows   = "
          f"{auc([r['step'] for r in allpos], [r['step'] for r in allneg]):.4f}")
    print(f"AUC(-mu   -> recovery), all rows   = "
          f"{auc([r['mu'] for r in allpos], [r['mu'] for r in allneg]):.4f}")
    print(f"AUC(-NU   -> recovery), all rows   = "
          f"{auc([r['NU'] for r in allpos], [r['NU'] for r in allneg]):.4f}")

    print("\n" + "-" * 78)
    print("Per-eff-stratum breakdown of the band test (does H25 hold at "
          "every bias strength, or only pooled?)")
    print("-" * 78)
    print(f"{'eff':>5} {'bandN':>6} {'rec':>7} {'AUC mu':>8} {'AUC step':>9}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>6} "
                  f"{str(len(p))+'/'+str(len(sub)):>7}   (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>6} {str(len(p))+'/'+str(len(sub)):>7} "
              f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f} "
              f"{auc([r['step'] for r in p], [r['step'] for r in ng]):>9.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
