"""
GLV-HNP Phase 2, Thread 25: is mu = lambda_1(L2) a genuine SECOND coordinate,
or is its apparent power in Thread 24's W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry:

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (17-bit bracket from
       W4, where the nearest-plane certificate gives no answer either way),
       AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate and (NU, mu) is a 2-parameter
  viability test.  If no: mu's apparent power in W5 is entirely mediated by
  NU (both are large/small together across curves) and the closed-form
  separator should be retired as redundant with NU.

Secondary (also pre-registered): W1b of glv_hnp_phase2_gsprofile.py found
that the GS profile head is m exact copies of lambda_1(L2), and the step to
the second block vanishes as the K1 wall is crossed.  Test whether
  step = log2(||b*_{m+1}||) - log2(||b*_1||)
predicts recovery better than NU or mu, globally and inside the NU band.

Data: reads the JSON dump produced by
  python3 glv_hnp_phase2_gsprofile_strat.py --dump-json <path>
(500 rows: 5 eff strata x 20 17-bit curves x 5 seeds, dim 24, float GS,
 justified safe by W0/W4: relative NU error vs exact Fractions ~1e-15).
No new lattice work is done here -- this is pure re-analysis.

Run: python3 glv_hnp_phase2_gsprofile_strat.py --dump-json rows.json
     python3 glv_hnp_phase2_nu_condition.py rows.json
"""

import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from glv_hnp_phase2_gsprofile import auc, spearman

# 17-bit NU bracket from Thread 24 W4 (glv_hnp_phase2_gsprofile.py output):
#   sufficient NU < 1.040 (all failures above this), necessary NU > 2.199
#   (all successes below this).  Ambiguous band is the open interval between.
NU_LO, NU_HI = 1.040, 2.199


def report_auc(name, rows, key, sign=-1):
    pos = [r[key] for r in rows if r['ok']]
    neg = [r[key] for r in rows if not r['ok']]
    if not pos or not neg:
        print(f"  AUC(-{name:<8}) : degenerate ({len(pos)} pos / {len(neg)} neg)")
        return float('nan')
    a = auc(pos, neg) if sign < 0 else auc(neg, pos)
    print(f"  AUC(-{name:<8}) = {a:.4f}   (N={len(rows)}, "
          f"{len(pos)} rec / {len(neg)} fail)")
    return a


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: glv_hnp_phase2_nu_condition.py <rows.json>")
        sys.exit(1)

    rows = json.load(open(sys.argv[1]))
    print("=" * 78)
    print("Thread 25 — does mu separate INSIDE the NU-ambiguous band? (H25)")
    print("=" * 78)
    print(f"\n{len(rows)} rows loaded, NU band = [{NU_LO}, {NU_HI}] "
          "(17-bit bracket, Thread 24 W4)")

    below = [r for r in rows if r['NU'] < NU_LO]
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    above = [r for r in rows if r['NU'] > NU_HI]
    print(f"\nNU < {NU_LO}:  N={len(below):4d}  rec="
          f"{sum(1 for r in below if r['ok'])}/{len(below)}")
    print(f"band          :  N={len(band):4d}  rec="
          f"{sum(1 for r in band if r['ok'])}/{len(band)}")
    print(f"NU > {NU_HI}:  N={len(above):4d}  rec="
          f"{sum(1 for r in above if r['ok'])}/{len(above)}")

    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) inside the band")
    print("-" * 78)
    if band:
        a_mu_band = report_auc('mu', band, 'mu')
        a_nu_band = report_auc('NU', band, 'NU')
        a_nh_band = report_auc('nu_hat', band, 'nuhat')
        a_step_band = report_auc('step', band, 'step')
        n_pos = sum(1 for r in band if r['ok'])
        n_neg = len(band) - n_pos
        print(f"\nband class balance: {n_pos} recover / {n_neg} fail "
              f"(N={len(band)})")
        verdict = "HOLDS" if (not math.isnan(a_mu_band) and a_mu_band >= 0.8) \
            else ("DEGENERATE (band too pure to test)" if n_pos == 0 or n_neg == 0
                  else "FALSIFIED")
        print(f"\nH25 verdict: {verdict}  (AUC(-mu) in band = {a_mu_band:.4f}, "
              f"threshold 0.8)")
    else:
        print("band is EMPTY -- cannot test H25 on this dataset.")

    print("\n" + "-" * 78)
    print("Global AUCs (for reference / comparison to Thread 24 W5 pooled)")
    print("-" * 78)
    report_auc('mu', rows, 'mu')
    report_auc('NU', rows, 'NU')
    report_auc('nu_hat', rows, 'nuhat')
    report_auc('step', rows, 'step')

    print("\n" + "-" * 78)
    print("Secondary: does step = log2||b*_(m+1)|| - log2||b*_1|| predict the "
          "wall?")
    print("-" * 78)
    steps_pos = [r['step'] for r in rows if r['ok']]
    steps_neg = [r['step'] for r in rows if not r['ok']]
    print(f"step | success : mean {sum(steps_pos)/len(steps_pos):.3f}  "
          f"min {min(steps_pos):.3f}  max {max(steps_pos):.3f}")
    print(f"step | failure : mean {sum(steps_neg)/len(steps_neg):.3f}  "
          f"min {min(steps_neg):.3f}  max {max(steps_neg):.3f}")
    print(f"Spearman(step, NU)     = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu)     = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")
    if band:
        print(f"Spearman(step, NU) in band = "
              f"{spearman([r['step'] for r in band], [r['NU'] for r in band]):.4f}")

    print("\n" + "-" * 78)
    print("Per-stratum breakdown inside the band (does H25 hold uniformly "
          "or only in aggregate?)")
    print("-" * 78)
    effs = sorted(set(r['effq'] for r in rows))
    print(f"{'eff':>5} {'N band':>7} {'rec':>7} {'AUC mu':>8} {'AUC NU':>8}")
    for eff in effs:
        sub = [r for r in band if r['effq'] == eff]
        if not sub:
            print(f"{eff:>5.2f} {0:>7} {'--':>7} {'--':>8} {'--':>8}")
            continue
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        rec = f"{len(pos)}/{len(sub)}"
        if pos and neg:
            am = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
            an = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
            print(f"{eff:>5.2f} {len(sub):>7} {rec:>7} {am:>8.4f} {an:>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} {rec:>7} {'(deg)':>8} {'(deg)':>8}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
