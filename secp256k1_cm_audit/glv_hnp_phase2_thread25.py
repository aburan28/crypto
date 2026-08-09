"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry.  W5/W6 of the
parent script established that recovery = f(NU, X), with NU the exact BDD
(nearest-plane) certificate and X ~ mu = lambda_1(L2)-driven, and that NU and
nu_hat*sqrt(eff) are mutually UNCORRELATED at fixed bias strength (Spearman
in [-0.28, +0.16] across every eff stratum).  NU is a sound sufficient
certificate (AUC 0.978, zero false positives) but a size-degrading separator
(bracket widens from [1.19, 1.87] at 12 bits to [1.04, 2.20] at 17 bits), so
~24% of instances at 17 bits are NU-ambiguous.  This script asks whether mu
resolves that ambiguity.

H25 (primary): within the NU-ambiguous band 1.04 <= NU <= 2.20 (17-bit
     bracket from Thread 23's exp V3/V4), AUC(-mu -> Kannan-LLL recovery)
     stays >= 0.8.
     If YES: mu is a genuine second coordinate and (NU, mu) is a real
     2-parameter viability test.
     If NO: mu's apparent power (W5: AUC 0.42-0.92 across eff strata) is
     entirely mediated by NU, and the W5 result is a stratification
     artifact of eff, not evidence for a second mechanism.

W9 (secondary, one-line add suggested by Thread 24's W1b): W1b showed the L0
   Gram-Schmidt profile is m exact copies of lambda_1(L2) in the head block,
   and the step to the second block vanishes right as the K1 wall is
   crossed.  Define
       step = log2(||b*_{m+1}||) - log2(||b*_1||)          (prof[m] - prof[0]
                                                              in log2, 0-indexed)
   and test whether step -> 0 predicts recovery better than NU or mu alone,
   and whether it is a distinct signal (Spearman vs NU, vs mu).

Same 17-bit construction as glv_hnp_phase2_gsprofile_strat.py (5 eff strata x
20 curves x 5 seeds, M17=12, dim 24), so results are directly comparable to
the W5/W6 table.  Adds --dump-json so the raw table survives the run for
future re-analysis without rebuilding it (requested by the 2026-08-07 #2 log
entry).

Run: python3 glv_hnp_phase2_thread25.py [--dump-json PATH]
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

# 17-bit bracket from the 2026-08-07 (Thread 23) log entry, exp V3/V4:
# NU < 1.040 -> 33/33 recover; NU > 2.199 -> 50/50 fail.
NU_LO, NU_HI = 1.040, 2.199


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-json", default=None,
                     help="write the raw instance table to this path")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — does mu resolve the NU-ambiguous band?")
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
                step = math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if args.dump_json:
        dump = [{k: v for k, v in r.items() if k != 'prof' and k != 'nus'}
                for r in rows]
        with open(args.dump_json, "w") as f:
            json.dump(dump, f)
        print(f"wrote {len(dump)} rows to {args.dump_json}")

    print("\n" + "-" * 78)
    print(f"EXP H25: AUC(-mu -> recovery) inside the NU-ambiguous band "
          f"[{NU_LO}, {NU_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos_b = [r for r in band if r['ok']]
    neg_b = [r for r in band if not r['ok']]
    print(f"band N = {len(band)} of {len(rows)} ({100*len(band)/len(rows):.1f}%), "
          f"recover {len(pos_b)}/{len(band)}")
    if pos_b and neg_b:
        a_mu_band = auc([r['mu'] for r in pos_b], [r['mu'] for r in neg_b])
        a_nh_band = auc([r['nuhat'] for r in pos_b], [r['nuhat'] for r in neg_b])
        a_nu_band = auc([r['NU'] for r in pos_b], [r['NU'] for r in neg_b])
        print(f"  AUC(-mu)      inside band = {a_mu_band:.4f}   "
              f"(H25 target: >= 0.80)")
        print(f"  AUC(-nu_hat)  inside band = {a_nh_band:.4f}")
        print(f"  AUC(-NU)      inside band = {a_nu_band:.4f}  "
              f"(expect ~0.5 — NU is constant-ish inside its own ambiguous band)")
        verdict = "CONFIRMED" if a_mu_band >= 0.80 else "FALSIFIED"
        print(f"  H25: {verdict}")
    else:
        print("  degenerate band (all one class) — cannot compute AUC")

    print("\nPer-stratum band breakdown (does the band result hold within eff too):")
    print(f"{'eff':>5} {'band N':>7} {'rec':>7} | {'AUC -mu':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if len(sub) < 4 or not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | {'(n/a)':>8}")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>7} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu:>8.4f}")

    print("\n" + "-" * 78)
    print("EXP W9: step = log2||b*_{m+1}|| - log2||b*_1||  vs. recovery")
    print("-" * 78)
    print("NOTE ON SIGN: auc(pos,neg) = P(pos < neg). H24/W1b's hypothesis is")
    print("'step -> 0 predicts the wall', i.e. recovery goes with LARGER step")
    print("-- the OPPOSITE convention from NU/mu/nu_hat (where smaller predicts")
    print("recovery). So auc(pos,neg) on raw step undershoots 0.5 exactly when")
    print("step is doing the predicted job; report 1-auc as 'AUC(+step)'.")
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    a_step_raw = auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg])
    a_nu_all = auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg])
    a_mu_all = auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg])
    print(f"pooled (N={len(rows)}):")
    print(f"  AUC(+step) = {1-a_step_raw:.4f}   AUC(-NU) = {a_nu_all:.4f}   "
          f"AUC(-mu) = {a_mu_all:.4f}")
    print(f"  Spearman(step, NU) = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu) = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\nstep inside the NU-ambiguous band:")
    if pos_b and neg_b:
        a_step_band_raw = auc([r['step'] for r in pos_b], [r['step'] for r in neg_b])
        print(f"  AUC(+step) inside band = {1-a_step_band_raw:.4f}")

    print(f"\n{'eff':>5} {'mean step|ok':>13} {'mean step|fail':>15}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            continue
        print(f"{eff:>5.2f} {sum(r['step'] for r in pos)/len(pos):>13.4f} "
              f"{sum(r['step'] for r in neg)/len(neg):>15.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
