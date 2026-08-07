"""
GLV-HNP Phase 2, Thread 25: does mu = lambda_1(L2) carry information beyond
NU, or is its apparent power in Thread 24's W5 entirely mediated by NU?

Context (RESEARCH_AUTOLAB_LOG.md, 2026-08-07 #2, Thread 24):
  W5/W6 showed NU (exact BDD certificate) and nu_hat*sqrt(eff) (closed form,
  nu_hat ~ mu/sqrt(det L2)) are mutually uncorrelated at fixed eff, yet both
  predict Kannan-LLL recovery with AUC > 0.75 in every non-degenerate
  17-bit stratum. W4 gave NU a sound but size-degrading bracket:
      sufficient NU < 1.040 (71 TP / 0 FP, n=71)
      necessary  NU > 2.199
  i.e. NU alone gives no answer for 1.040 <= NU <= 2.199.

  H25: within that ambiguous band, AUC(-mu -> recovery) stays >= 0.8.
  If yes: (NU, mu) is a genuine 2-parameter viability test. If no (mu's
  power was entirely mediated by NU): W5's per-stratum AUC was a
  stratification artifact and the closed form should be retired.

Secondary (also proposed 2026-08-07 #2): W1b showed the GS profile head is
m exact copies of lambda_1(L2) and the step to the second block vanishes
right as the K1 wall is crossed. Define
    step = log2(||b*_{m+1}||) - log2(||b*_1||) = log2(prof[m]) - log2(prof[0])
and test whether step predicts recovery, and whether it is redundant with mu.

Reuses glv_hnp_phase2_gsprofile_strat.collect_rows() verbatim (same 17-bit,
5-stratum, 20-curve, 5-seed table as the parent Thread 24 run) so this is
re-analysis of the same instances, not new data.

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman
from glv_hnp_phase2_gsprofile_strat import collect_rows, EFFS, M17

# W4 (17-bit, dim 24) bracket: sufficient NU < 1.040, necessary NU > 2.199.
NU_LO, NU_HI = 1.040, 2.199

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — does mu carry signal beyond NU inside the ambiguous band?")
    print("=" * 78)

    rows, curves17 = collect_rows()
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves, {len(rows)} instances "
          f"(same table as Thread 24 W5/W6)")

    for r in rows:
        r['step'] = math.log2(r['prof'][M17]) - math.log2(r['prof'][0])

    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) restricted to the NU ambiguous band")
    print(f"     band: {NU_LO} <= NU <= {NU_HI}  (W4 sufficient/necessary bracket)")
    print("-" * 78)

    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band N = {len(band)} of {len(rows)}   "
          f"recovered {len(pos)}   failed {len(neg)}")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"AUC(-mu     -> recovery) = {a_mu:.4f}   (H25 threshold: 0.80)")
        print(f"AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"AUC(-NU     -> recovery) = {a_nu:.4f}   (near 0.5 expected: "
              f"NU range is compressed by construction)")
        verdict = "SURVIVES" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 {verdict}: AUC(-mu) = {a_mu:.4f} "
              f"{'>=' if a_mu >= 0.80 else '<'} 0.80")
    else:
        print("degenerate: band has no both-classes data")

    print("\nper-eff-stratum breakdown inside the band:")
    print(f"{'eff':>5} {'N':>4} {'rec':>6} | {'AUC mu':>8} {'AUC nu_hat':>11} "
          f"{'AUC NU':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>6} | (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>4} {str(len(p))+'/'+str(len(sub)):>6} | "
              f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f} "
              f"{auc([r['nuhat'] for r in p], [r['nuhat'] for r in ng]):>11.4f} "
              f"{auc([r['NU'] for r in p], [r['NU'] for r in ng]):>8.4f}")

    print("\n" + "-" * 78)
    print("Secondary: step = log2||b*_{m+1}|| - log2||b*_1|| as a predictor")
    print("-" * 78)
    pos_all = [r for r in rows if r['ok']]
    neg_all = [r for r in rows if not r['ok']]
    a_step = auc([r['step'] for r in pos_all], [r['step'] for r in neg_all])
    print(f"pooled AUC(-step -> recovery)  = {a_step:.4f}   "
          f"(0.5 = uninformative; report as-is, sign not assumed)")
    print(f"pooled AUC(+step -> recovery)  = {1 - a_step:.4f}")
    print(f"Spearman(step, mu)  = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")
    print(f"Spearman(step, NU)  = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")

    print("\nper-eff-stratum AUC(step) (sign chosen to be >=0.5 if any signal):")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            continue
        a = auc([r['step'] for r in p], [r['step'] for r in ng])
        print(f"  eff={eff:.2f}  N={len(sub):<4} AUC(-step)={a:.4f}  "
              f"AUC(+step)={1-a:.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
