"""
GLV-HNP Phase 2, Thread 25: conditioning on NU -- is mu a second coordinate?

Pre-registered by the 2026-08-07 #2 log entry (end of the Thread 24 session,
"Next step proposal"). W5/W6 of Thread 24 established:

  - recovery = f(NU, X) for some X that behaves like mu = lambda_1(L2): mu
    (equiv. nu_hat) separates curves cross-curve with eff HELD FIXED
    (AUC 0.75-0.93 per stratum), which NU does not (AUC 0.35-0.73, one
    stratum inverted).
  - NU and nu_hat*sqrt(eff) are mutually uncorrelated at fixed eff
    (Spearman in [-0.28, +0.16] across strata) -- so if mu really is a
    second, NU-independent coordinate, conditioning on NU (rather than on
    eff) should ALSO fail to erase mu's power.

  H25: within the ambiguous NU band [1.04, 2.20] (the 17-bit bracket from
       Thread 24 W4, where NU's nearest-plane certificate gives no answer
       either way), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if AUC(-mu) inside the band drops toward 0.5, mu's apparent
       cross-curve power is entirely mediated by NU (i.e. a monotone
       function of NU that just happens to also separate at fixed eff) and
       the closed form should be retired as a stratification artifact.

Secondary (also pre-registered, from W1b): the GS profile of L0 has a flat
head of m copies of lambda_1(L2) followed by a tail block that starts to
move once the K1 wall is crossed. Quantify the jump with

    step = log2(||b*_{m+1}||) - log2(||b*_1||)     (1-based, so prof[m] - prof[0]
                                                      in the 0-based array)

and test whether step -> 0 predicts recovery better than NU or mu, both
pooled and inside the H25 band.

Data: identical generator to glv_hnp_phase2_gsprofile_strat.py (5 eff strata
x 20 curves x 5 seeds at 17 bits, dim 24). Float GS, justified safe by
Thread 24 W0/W4 (max relative NU error ~6.6e-16 at dim 24 vs exact Fractions).

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# 17-bit bracket from Thread 24 W4 (300 instances, dim 24):
#   sufficient NU < 1.040 , necessary NU > 2.199
BAND_LO, BAND_HI = 1.040, 2.199


def summarize(label, sub):
    pos = [r for r in sub if r['ok']]
    neg = [r for r in sub if not r['ok']]
    if not pos or not neg:
        print(f"{label:<28} N={len(sub):<4} rec={len(pos)}/{len(sub):<5} "
              f"(degenerate -- one class empty, AUC undefined)")
        return None
    a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
    a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
    a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
    print(f"{label:<28} N={len(sub):<4} rec={len(pos)}/{len(sub):<5} "
          f"AUC(-mu)={a_mu:.4f}  AUC(-NU)={a_nu:.4f}  AUC(-step)={a_st:.4f}")
    return {'auc_mu': a_mu, 'auc_nu': a_nu, 'auc_step': a_st, 'n': len(sub)}


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 -- conditioning on NU: is mu a second coordinate? (H25)")
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
                m = r['k'] // 2
                b1, bm1 = r['prof'][0], r['prof'][m]
                step = (math.log2(bm1) - math.log2(b1)
                        if b1 > 0 and bm1 > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    nan_steps = sum(1 for r in rows if math.isnan(r['step']))
    if nan_steps:
        print(f"WARNING: {nan_steps}/{len(rows)} rows have degenerate step "
              "(prof[0] or prof[m] == 0) -- dropped from step AUCs below.")
        rows_step_ok = [r for r in rows if not math.isnan(r['step'])]
    else:
        rows_step_ok = rows

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SANITY: is the head really flat?  prof[0] vs prof[m-1] (last head)")
    print("-" * 78)
    head_gap = [abs(math.log2(r['prof'][r['k']//2 - 1]) - math.log2(r['prof'][0]))
                for r in rows if r['prof'][0] > 0]
    print(f"|log2 prof[m-1] - log2 prof[0]| over {len(head_gap)} rows: "
          f"mean {sum(head_gap)/len(head_gap):.4f} bits, "
          f"max {max(head_gap):.4f} bits")
    print("(near 0 confirms the head is a flat block, so step as specified")
    print(" -- b*_1 vs b*_{m+1} -- is interchangeable with the adjacent-pair")
    print(" step b*_m vs b*_{m+1} used informally in the Thread 24 log.)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("H25 PRIMARY: AUC(-mu -> recovery) inside vs outside the NU band")
    print(f"band = [{BAND_LO}, {BAND_HI}]  (Thread 24 W4, 17-bit bracket)")
    print("-" * 78)
    below = [r for r in rows if r['NU'] < BAND_LO]
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    above = [r for r in rows if r['NU'] > BAND_HI]
    summarize("NU < 1.040 (certificate zone)", below)
    band_stats = summarize("1.040 <= NU <= 2.199 (BAND)", band)
    summarize("NU > 2.199 (should mostly fail)", above)
    summarize("pooled (all 500)", rows)

    print("\nband composition by eff stratum (does the band just re-select")
    print("one eff value, making this equivalent to the W5 test?):")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        print(f"  eff={eff:.2f}: {len(sub)}/{len([r for r in rows if r['effq']==eff])} "
              f"of that stratum falls in the band")

    print("\nband composition by curve (n) -- top 5 by count:")
    from collections import Counter
    cnt = Counter(r['n'] for r in band)
    for n_, c in cnt.most_common(5):
        print(f"  n={n_}: {c} instances in band")

    if band_stats is not None:
        verdict = "HOLDS" if band_stats['auc_mu'] >= 0.8 else "FALSIFIED"
        print(f"\nH25 verdict: AUC(-mu) in band = {band_stats['auc_mu']:.4f} "
              f"(threshold 0.80) -> H25 {verdict}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECONDARY: does step = log2||b*_{m+1}|| - log2||b*_1|| predict the")
    print("wall better than NU or mu?  (sign unknown a priori -- report both)")
    print("-" * 78)
    pos = [r['step'] for r in rows_step_ok if r['ok']]
    neg = [r['step'] for r in rows_step_ok if not r['ok']]
    print(f"AUC(-step -> recovery), pooled N={len(rows_step_ok)}: "
          f"{auc(pos, neg):.4f}   (>0.5 => smaller/more-negative step "
          "=> more likely to recover)")
    print(f"AUC(+step -> recovery), pooled: {auc(neg, pos):.4f}")
    print(f"step | success: mean {sum(pos)/len(pos):.3f}  "
          f"min {min(pos):.3f}  max {max(pos):.3f}")
    print(f"step | failure: mean {sum(neg)/len(neg):.3f}  "
          f"min {min(neg):.3f}  max {max(neg):.3f}")
    print(f"Spearman(step, NU)  = {spearman([r['step'] for r in rows_step_ok], [r['NU'] for r in rows_step_ok]):.4f}")
    print(f"Spearman(step, mu)  = {spearman([r['step'] for r in rows_step_ok], [r['mu'] for r in rows_step_ok]):.4f}")

    band_step_ok = [r for r in band if not math.isnan(r['step'])]
    posb = [r['step'] for r in band_step_ok if r['ok']]
    negb = [r['step'] for r in band_step_ok if not r['ok']]
    if posb and negb:
        print(f"\ninside the H25 band: AUC(-step -> recovery) = "
              f"{auc(posb, negb):.4f}   AUC(+step -> recovery) = "
              f"{auc(negb, posb):.4f}   N={len(band_step_ok)}")
    else:
        print("\ninside the H25 band: one class empty, step AUC undefined")

    print("\nrobustness: is AUC(+step) inside the band an eff-stratum mix")
    print("artifact, or does it hold within each eff stratum separately?")
    for eff in EFFS:
        sub = [r for r in band_step_ok if r['effq'] == eff]
        p = [r['step'] for r in sub if r['ok']]
        ng = [r['step'] for r in sub if not r['ok']]
        if p and ng:
            print(f"  eff={eff:.2f}  N={len(sub):<4} rec={len(p)}/{len(sub):<4} "
                  f"AUC(+step)={auc(ng, p):.4f}  AUC(-mu)={auc([r['mu'] for r in sub if r['ok']], [r['mu'] for r in sub if not r['ok']]):.4f}")
        else:
            print(f"  eff={eff:.2f}  N={len(sub):<4} rec={len(p)}/{len(sub):<4} "
                  "(degenerate)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
