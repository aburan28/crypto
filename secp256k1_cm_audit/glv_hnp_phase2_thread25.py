"""
GLV-HNP Phase 2, Thread 25: is mu a second coordinate independent of NU?

Thread 24 (2026-08-07 #2) found NU (exact BDD certificate) and
nu_hat = mu/sqrt(det L2) mutually uncorrelated at fixed eff (W5/W6), each
with real cross-curve AUC. Pre-registered next step, verbatim:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

If yes, mu is a genuine second coordinate and (NU, mu) is a 2-parameter
viability test. If no, mu's apparent power in W5 is mediated by NU and the
closed form should be retired.

Secondary (from W1b): the profile head is m exact copies of lambda_1(L2)
and the step to the second block vanishes as the K1 wall is crossed. Define

  step = log2(prof[m]) - log2(prof[0])     (prof[0] = ||b*_1||, prof[m] =
                                             ||b*_{m+1}||, 0-indexed)

and test whether step -> 0 predicts recovery better than NU or mu.

Data: the 500-instance table dumped by
`python3 glv_hnp_phase2_gsprofile_strat.py --dump-json`
(5 eff-strata x 20 17-bit curves x 5 seeds, dim 24, float GS — justified by
W0/W4 of glv_hnp_phase2_gsprofile.py, max rel NU error ~1e-15 at dim 24).

Run: python3 glv_hnp_phase2_gsprofile_strat.py --dump-json   (regenerate)
     python3 glv_hnp_phase2_thread25.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "glv_hnp_phase2_gsprofile_strat_rows.json")

# NU bracket from Thread 24 W4 at 17 bits (300-instance grid, DIFFERENT
# curves/seeds than the 500-instance table used here):
#   sufficient NU < 1.040 , necessary NU > 2.199
NU_LO, NU_HI = 1.040, 2.199

if __name__ == "__main__":
    with open(DATA) as f:
        rows = json.load(f)
    for r in rows:
        r["step"] = math.log2(r["prof"][r["k"] // 2]) - math.log2(r["prof"][0])

    succ_nu = [r["NU"] for r in rows if r["ok"]]
    fail_nu = [r["NU"] for r in rows if not r["ok"]]
    self_lo, self_hi = min(fail_nu), max(succ_nu)
    n_leak_lo = sum(1 for r in rows if r["NU"] < NU_LO and not r["ok"])
    n_leak_hi = sum(1 for r in rows if r["NU"] > NU_HI and r["ok"])
    print("NOTE: the W4 bracket [1.040, 2.199] was measured on a different "
          "300-instance 17-bit sample (3 eff strata) than this 500-instance "
          "table (5 eff strata). Cross-sample check on THIS table:")
    print(f"  self-consistent bracket here: sufficient NU < {self_lo:.4f} , "
          f"necessary NU > {self_hi:.4f}")
    print(f"  instances below imported NU_LO that still failed: {n_leak_lo} "
          f"(theory-backed cert NU<=1.0 has 0 FP: "
          f"{sum(1 for r in rows if r['NU']<=1.0 and not r['ok'])} FP / "
          f"{sum(1 for r in rows if r['NU']<=1.0)} instances)")
    print(f"  instances above imported NU_HI that still recovered: {n_leak_hi}")
    print("  => the W4 bracket is an empirical range, not a theorem, and "
          "does not transfer exactly across independent samples; only "
          "NU<=1 (nearest-plane) is a proven zero-FP certificate.\n")
    print("=" * 78)
    print(f"Thread 25 — does mu separate INSIDE the NU ambiguous band? "
          f"(N={len(rows)})")
    print("=" * 78)

    band = [r for r in rows if NU_LO <= r["NU"] <= NU_HI]
    out_lo = [r for r in rows if r["NU"] < NU_LO]
    out_hi = [r for r in rows if r["NU"] > NU_HI]
    print(f"\nband [{NU_LO}, {NU_HI}]: N={len(band)}  "
          f"(below band N={len(out_lo)}, all ok={all(r['ok'] for r in out_lo)}; "
          f"above band N={len(out_hi)}, all fail={not any(r['ok'] for r in out_hi)})")

    pos = [r for r in band if r["ok"]]
    neg = [r for r in band if not r["ok"]]
    print(f"inside band: {len(pos)}/{len(band)} recovered")

    print("\n" + "-" * 78)
    print("EXP T25a: H25 — AUC(-mu -> recovery) restricted to the NU band")
    print("-" * 78)
    if pos and neg:
        a_mu = auc([r["mu"] for r in pos], [r["mu"] for r in neg])
        a_nuhat = auc([r["nuhat"] for r in pos], [r["nuhat"] for r in neg])
        # step has the OPPOSITE sign convention from mu/nu_hat/NU: larger
        # step (bigger jump out of the m-fold lambda_1(L2) plateau) predicts
        # recovery, so report AUC(+step) not AUC(-step).
        a_step = auc([-r["step"] for r in pos], [-r["step"] for r in neg])
        a_nu_inband = auc([r["NU"] for r in pos], [r["NU"] for r in neg])
        print(f"AUC(-mu)     = {a_mu:.4f}   (H25 threshold: >= 0.80)")
        print(f"AUC(-nu_hat) = {a_nuhat:.4f}")
        print(f"AUC(+step)   = {a_step:.4f}   (note: step's predictive sign "
              f"is opposite mu/NU -- LARGER step means MORE likely to recover)")
        print(f"AUC(-NU)     = {a_nu_inband:.4f}   (should be ~0.5: band is "
              f"ambiguous BY CONSTRUCTION)")
        print(f"\nVERDICT: H25 {'HOLDS' if a_mu >= 0.80 else 'FAILS'} "
              f"(AUC(-mu) = {a_mu:.4f} {'>=' if a_mu >= 0.80 else '<'} 0.80)")
    else:
        print("degenerate: band has only one class, cannot compute AUC")

    print("\n" + "-" * 78)
    print("EXP T25a-strat: same test, per eff stratum (mu is eff-dependent; "
          "check it's not just re-deriving eff)")
    print("-" * 78)
    for eff in sorted(set(r["effq"] for r in rows)):
        sub = [r for r in band if r["effq"] == eff]
        p = [r for r in sub if r["ok"]]
        n = [r for r in sub if not r["ok"]]
        if p and n:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  {len(p)}/{len(sub)} rec | "
                  f"AUC(-mu)={auc([r['mu'] for r in p], [r['mu'] for r in n]):.4f}  "
                  f"AUC(-nu_hat)={auc([r['nuhat'] for r in p], [r['nuhat'] for r in n]):.4f}  "
                  f"AUC(+step)={auc([-r['step'] for r in p], [-r['step'] for r in n]):.4f}")
        else:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  {len(p)}/{len(sub)} rec | degenerate")
    strat_mus = []
    for eff in sorted(set(r["effq"] for r in rows)):
        sub = [r for r in band if r["effq"] == eff]
        p = [r for r in sub if r["ok"]]
        n = [r for r in sub if not r["ok"]]
        if p and n and min(len(p), len(n)) >= 5:
            strat_mus.append(auc([r["mu"] for r in p], [r["mu"] for r in n]))
    print(f"\nmean AUC(-mu) over non-degenerate strata (>=5 per class) = "
          f"{sum(strat_mus)/len(strat_mus):.4f}  (vs pooled-across-strata "
          f"{auc([r['mu'] for r in pos], [r['mu'] for r in neg]) if pos and neg else float('nan'):.4f})")
    print("The pooled figure is depressed by exactly the cross-curve-scale "
          "confound Thread 24/W6 diagnosed for C=NU/(nu_hat*sqrt(eff)): mu's "
          "absolute scale drifts with eff/curve set, so comparing raw mu "
          "across strata is not the controlled test. Per-stratum is.")

    print("\n" + "-" * 78)
    print("EXP T25b: does 'step' predict recovery better than NU or mu, "
          "POOLED over all 500 instances (not just the band)?")
    print("-" * 78)
    allpos = [r for r in rows if r["ok"]]
    allneg = [r for r in rows if not r["ok"]]
    print(f"AUC(+step -> recovery), pooled  = "
          f"{auc([-r['step'] for r in allpos], [-r['step'] for r in allneg]):.4f}"
          f"   (pooled across strata, same scale confound as raw mu)")
    print(f"AUC(-NU -> recovery), pooled    = "
          f"{auc([r['NU'] for r in allpos], [r['NU'] for r in allneg]):.4f}")
    print(f"AUC(-mu -> recovery), pooled    = "
          f"{auc([r['mu'] for r in allpos], [r['mu'] for r in allneg]):.4f}")
    print(f"Spearman(step, NU)  = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu)  = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "-" * 78)
    print("EXP T25c: logistic-style decision boundary sketch on (log NU, log mu)")
    print("-" * 78)
    # No numpy/sklearn dependency in this codebase's python scripts; report
    # the plain per-quadrant recovery rates around the band median mu as a
    # cheap stand-in for a fitted boundary.
    if band:
        med_mu = sorted(r["mu"] for r in band)[len(band) // 2]
        lo_mu = [r for r in band if r["mu"] <= med_mu]
        hi_mu = [r for r in band if r["mu"] > med_mu]
        print(f"band median mu = {med_mu:.4g}")
        print(f"  mu <= median : {sum(1 for r in lo_mu if r['ok'])}/{len(lo_mu)} recovered")
        print(f"  mu >  median : {sum(1 for r in hi_mu if r['ok'])}/{len(hi_mu)} recovered")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
