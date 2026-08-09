"""
GLV-HNP Phase 2, Thread 25: is mu a genuine SECOND coordinate, orthogonal to
NU, or is its apparent power in W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  H25: within the ambiguous NU band (NU theory gives no answer -- neither
       NU<=1 sufficiency nor NU>bracket necessity fires), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

  Falsifier: if AUC(-mu) inside the band collapses to ~0.5, mu's apparent
       power in W5 was a stratification artifact riding on NU after all,
       and the closed form nu_hat should be retired as "NU in disguise"
       rather than treated as a second mechanism.

Band definition, taken directly from the W4 bracket (17 bits, 300 instances,
RESEARCH_AUTOLAB_LOG.md 2026-08-07 #2): sufficient NU < 1.040 (min over
failures), necessary NU > 2.199 (max over successes).  So the ambiguous band
is 1.040 <= NU <= 2.199 -- nearest-plane theory abstains there.

Secondary (also pre-registered): the two-block GS profile of Thread 24 (W1b)
showed the step from the head block (m copies of lambda_1(L2)) to the second
block vanishes exactly as the K1 wall is crossed.  Define
    step = log2(||b*_{m+1}||) - log2(||b*_1||)   (prof[m] vs prof[0])
and test whether -step (small step -> more likely to recover) beats NU and mu
as a separator, using the SAME 500-instance table so the comparison is exact.

Data: regenerates the identical 17-bit, 5-eff-stratum x 20-curve x 5-seed
table as glv_hnp_phase2_gsprofile_strat.py (same search_curves call, same
d_trial derivation, same seeds) -- no new lattice experiments, pure
re-analysis as proposed.  exact=False (float GS) per W0/W4: max relative NU
error vs exact Fractions was 6.6e-16 at this dimension, so float is safe and
much faster to regenerate.

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 -- does mu separate INSIDE the NU-ambiguous band?")
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
                          'step': math.log2(r['prof'][r['k'] // 2])
                                  - math.log2(r['prof'][0])
                                  if r['prof'][0] > 0 and r['prof'][r['k'] // 2] > 0
                                  else float('nan')})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("W4 bracket (from the 2026-08-07 #2 log entry, re-derived here "
          "as a sanity check)")
    print("-" * 78)
    pos_all = [r['NU'] for r in rows if r['ok']]
    neg_all = [r['NU'] for r in rows if not r['ok']]
    lo = min(neg_all)   # sufficient: NU below this, only successes seen
    hi = max(pos_all)   # necessary: NU above this, only failures seen
    print(f"N={len(rows)}  rec={len(pos_all)}/{len(rows)}")
    print(f"sufficient NU < {lo:.3f}  (min over {len(neg_all)} failures)")
    print(f"necessary  NU > {hi:.3f}  (max over {len(pos_all)} successes)")
    print(f"ambiguous band: {lo:.3f} <= NU <= {hi:.3f}")

    band = [r for r in rows if lo <= r['NU'] <= hi]
    below = [r for r in rows if r['NU'] < lo]
    above = [r for r in rows if r['NU'] > hi]
    print(f"\npartition: below-band {len(below)} (rec "
          f"{sum(1 for r in below if r['ok'])}/{len(below)}), "
          f"in-band {len(band)} (rec {sum(1 for r in band if r['ok'])}/{len(band)}), "
          f"above-band {len(above)} (rec "
          f"{sum(1 for r in above if r['ok'])}/{len(above)})")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP T25a: H25 -- AUC(-mu -> recovery) INSIDE the ambiguous band")
    print("-" * 78)
    bpos = [r for r in band if r['ok']]
    bneg = [r for r in band if not r['ok']]
    if bpos and bneg:
        a_mu = auc([r['mu'] for r in bpos], [r['mu'] for r in bneg])
        a_nh = auc([r['nuhat'] for r in bpos], [r['nuhat'] for r in bneg])
        a_nu_band = auc([r['NU'] for r in bpos], [r['NU'] for r in bneg])
        print(f"in-band N={len(band)}  rec={len(bpos)}/{len(band)}")
        print(f"AUC(-mu     -> recovery) = {a_mu:.4f}   <-- H25 test statistic")
        print(f"AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"AUC(-NU     -> recovery) = {a_nu_band:.4f}   "
              f"(expected ~0.5: theory abstains here by construction of the band)")
        verdict = "SURVIVES" if a_mu >= 0.8 else (
            "PARTIAL" if a_mu >= 0.65 else "FALSIFIED")
        print(f"\nH25 threshold 0.8: {verdict}  (observed {a_mu:.4f})")
    else:
        print(f"degenerate band: {len(bpos)} pos, {len(bneg)} neg -- cannot compute AUC")

    print("\nper-stratum breakdown inside the band:")
    print(f"{'eff':>5} {'N':>5} {'rec':>7} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degenerate)':>8}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP T25b: pooled AUC comparison, mu vs NU vs nu_hat*sqrt(eff), "
          "ALL instances")
    print("-" * 78)
    ppos = [r for r in rows if r['ok']]
    pneg = [r for r in rows if not r['ok']]
    print(f"AUC(-mu)                = {auc([r['mu'] for r in ppos], [r['mu'] for r in pneg]):.4f}")
    print(f"AUC(-NU)                = {auc([r['NU'] for r in ppos], [r['NU'] for r in pneg]):.4f}")
    print(f"AUC(-nu_hat*sqrt(eff))  = "
          f"{auc([r['nuhat']*math.sqrt(r['eff']) for r in ppos], [r['nuhat']*math.sqrt(r['eff']) for r in pneg]):.4f}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP T25c (secondary): does the GS-profile step predict the wall?")
    print("step = log2||b*_{m+1}|| - log2||b*_1||   (0 == flat profile, no step)")
    print("-" * 78)
    srows = [r for r in rows if not math.isnan(r['step'])]
    spos = [r['step'] for r in srows if r['ok']]
    sneg = [r['step'] for r in srows if not r['ok']]
    print(f"N={len(srows)} (dropped {len(rows)-len(srows)} degenerate profiles)")
    print(f"step | success : mean {sum(spos)/len(spos):.3f}  "
          f"min {min(spos):.3f}  max {max(spos):.3f}")
    print(f"step | failure : mean {sum(sneg)/len(sneg):.3f}  "
          f"min {min(sneg):.3f}  max {max(sneg):.3f}")
    print(f"AUC(-step -> recovery) = {auc(spos, sneg):.4f}   "
          f"(compare AUC(-NU)={auc([r['NU'] for r in srows if r['ok']], [r['NU'] for r in srows if not r['ok']]):.4f}, "
          f"AUC(-mu)={auc([r['mu'] for r in srows if r['ok']], [r['mu'] for r in srows if not r['ok']]):.4f})")
    print(f"Spearman(step, NU) = {spearman([r['step'] for r in srows], [r['NU'] for r in srows]):.4f}")
    print(f"Spearman(step, mu) = {spearman([r['step'] for r in srows], [r['mu'] for r in srows]):.4f}")

    # in-band step check too, since that's the regime that matters
    bsrows = [r for r in band if not math.isnan(r['step'])]
    bspos = [r['step'] for r in bsrows if r['ok']]
    bsneg = [r['step'] for r in bsrows if not r['ok']]
    if bspos and bsneg:
        print(f"\nin-band AUC(-step -> recovery) = {auc(bspos, bsneg):.4f}  "
              f"(N={len(bsrows)})")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
