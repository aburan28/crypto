"""
GLV-HNP Phase 2, Thread 25: does mu = lambda_1(L2) separate recovery
INSIDE the NU-ambiguous band, or is its apparent power in Thread 24's W5
entirely mediated by NU?

Pre-registered (RESEARCH_AUTOLAB_LOG.md, 2026-08-07 run #2, Thread 24 "Next
step proposal"):

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  Falsifier: if mu's apparent power is entirely mediated by NU, AUC inside
  the band collapses toward 0.5 and the W5 result is a stratification
  artifact.

Secondary (also pre-registered): W1b showed the GS profile head is m exact
copies of lambda_1(L2) and the step to the second block vanishes as the K1
wall is crossed. Define
      step = log2(||b*_{m+1}||) - log2(||b*_1||)     (0-indexed: prof[m]-prof[0])
and test whether step -> 0 predicts recovery better than NU or mu.

This regenerates the 24b dataset (5 strata x 20 curves x 5 seeds, 17 bits,
dim 24) because that run did not persist rows to disk (see 2026-08-07 run #2
log: "Add a --dump-json flag so the table survives the run"). Float GS is
used throughout per W0 (max relative NU error ~1e-15 at this dimension).

Run: python3 glv_hnp_phase2_nu_stratify.py [--dump-json out.json]
"""

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

M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
# Bracket from Thread 24 W4 (17-bit): sufficient NU < 1.040, necessary NU > 2.199.
NU_LO, NU_HI = 1.040, 2.199


def build_rows():
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

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
                step = math.log2(r['prof'][M17]) - math.log2(r['prof'][0]) \
                    if r['prof'][0] > 0 and r['prof'][M17] > 0 else float('nan')
                rows.append({
                    'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                    'eff': k1b * k2b / n, 'effq': eff,
                    'lamstar': lam_star(lam, n),
                    'NU': r['NU'], 'mu': r['mu'], 'nuhat': r['nuhat'],
                    'step': step, 'seed': seed,
                })
    print(f"{len(rows)} instances in {time.time()-t0:.1f}s")
    return rows


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        dump_path = sys.argv[sys.argv.index("--dump-json") + 1]

    print("=" * 78)
    print("Thread 25 — does mu separate recovery INSIDE the NU-ambiguous band?")
    print("=" * 78)

    rows = build_rows()
    if dump_path:
        with open(dump_path, "w") as f:
            json.dump(rows, f)
        print(f"\ndumped {len(rows)} rows -> {dump_path}")

    print("\n" + "-" * 78)
    print(f"H25: stratify by NU. Ambiguous band = [{NU_LO}, {NU_HI}] (Thread 24 W4)")
    print("-" * 78)

    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    below = [r for r in rows if r['NU'] < NU_LO]
    above = [r for r in rows if r['NU'] > NU_HI]
    print(f"below band (NU < {NU_LO}): {len(below)}  "
          f"[{sum(1 for r in below if r['ok'])} recovered]")
    print(f"in band                : {len(band)}  "
          f"[{sum(1 for r in band if r['ok'])} recovered]")
    print(f"above band (NU > {NU_HI}): {len(above)}  "
          f"[{sum(1 for r in above if r['ok'])} recovered]")

    for name, sub in (("BELOW", below), ("IN-BAND", band), ("ABOVE", above)):
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{name:>8}: N={len(sub):<4} (degenerate: "
                  f"{len(pos)} pos / {len(neg)} neg)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{name:>8}: N={len(sub):<4} rec={len(pos)}/{len(sub)}  "
              f"AUC(-mu)={a_mu:.4f}  AUC(-nu_hat)={a_nh:.4f}  "
              f"AUC(-step)={a_step:.4f}")

    verdict = "UNKNOWN"
    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    if pos_b and neg_b:
        a_mu_band = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        verdict = "H25 HOLDS" if a_mu_band >= 0.8 else "H25 FALSIFIED"
        print(f"\nAUC(-mu) inside band = {a_mu_band:.4f}  -> {verdict}")
    else:
        print("\nband is degenerate (all-pos or all-neg): H25 untestable "
              "at this bracket/sample size")

    print("\n" + "-" * 78)
    print("Secondary: does the block-step statistic beat NU / mu?")
    print("-" * 78)
    pos_all = [r for r in rows if r['ok']]
    neg_all = [r for r in rows if not r['ok']]
    print(f"pooled AUC(-step) = "
          f"{auc([r['step'] for r in pos_all], [r['step'] for r in neg_all]):.4f}")
    print(f"pooled AUC(-mu)   = "
          f"{auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all]):.4f}")
    print(f"pooled AUC(-NU)   = "
          f"{auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all]):.4f}")
    print(f"Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
