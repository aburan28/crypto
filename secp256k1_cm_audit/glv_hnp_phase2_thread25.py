"""
GLV-HNP Phase 2, Thread 25: is mu a second coordinate once NU is fixed, and
does the head-to-second-block GS "step" beat NU/mu as a wall predictor?

Pre-registered by the 2026-08-07 (Thread 24) log entry, off the back of W5/W6
(glv_hnp_phase2_gsprofile_strat.py): NU and nu_hat (equivalently mu, since W5
showed mu tracks nu_hat to 3 decimal places within a stratum) are MUTUALLY
UNCORRELATED (Spearman ~ 0 or negative, every eff stratum) yet both predict
recovery. Two mutually uncorrelated predictors of the same outcome means they
carry different information; H25 asks whether that information is genuinely
2-dimensional or whether mu's apparent power is a stratification artifact of
NU.

  H25: within the ambiguous NU band (nearest-plane gives no verdict either
       way), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if AUC collapses to ~0.5 inside the band, mu's W5 power was
       entirely mediated by NU's correlation with eff/size, not a genuine
       second coordinate.

Secondary (W1b): the GS profile of L0 is m exact copies of lambda_1(L2) in
its head, and the step to the second block vanishes as the K1 wall is
crossed. Define
    step = log2(||b*_{m+1}||) - log2(||b*_1||)          (0-indexed: prof[m]-prof[0], in log2)
and test whether -step beats NU and mu as a recovery predictor.

Data: reuses the 500-instance table dumped by
`glv_hnp_phase2_gsprofile_strat.py --dump-json` (17 bits, dim 24, 20 curves,
5 eff strata x 5 seeds x 20 curves = 500 rows). No new lattice work.

Run: python3 glv_hnp_phase2_gsprofile_strat.py --dump-json   (once)
     python3 glv_hnp_phase2_thread25.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

DUMP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "glv_hnp_phase2_gsprofile_strat_rows.json")

# 17-bit two-sided NU bracket from Thread 23 (log 2026-08-07 #2, exp V3/V4,
# W4): sufficient NU < 1.040 -> 100% recovery, necessary NU > 2.199 -> 0%
# recovery. The ambiguous band is the interval between them.
NU_LO, NU_HI = 1.040, 2.199


def bandsplit(rows, lo, hi):
    return [r for r in rows if lo <= r['NU'] <= hi]


if __name__ == "__main__":
    if not os.path.exists(DUMP):
        sys.exit(f"missing {DUMP} — run "
                  "'python3 glv_hnp_phase2_gsprofile_strat.py --dump-json' first")
    with open(DUMP) as f:
        rows = json.load(f)
    print("=" * 78)
    print("Thread 25 — is mu a 2nd coordinate given NU, and does the GS-profile")
    print("            step to the second block beat NU/mu?")
    print("=" * 78)
    print(f"\n{len(rows)} instances loaded from {os.path.basename(DUMP)}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1 (H25): AUC(-mu -> recovery) INSIDE the ambiguous NU band")
    print(f"              [{NU_LO}, {NU_HI}], pooled and per eff stratum")
    print("-" * 78)
    band = bandsplit(rows, NU_LO, NU_HI)
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"\nband population: {len(band)}/{len(rows)} instances "
          f"({len(pos)} recover, {len(neg)} fail)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}  "
              f"(NU is constrained to a narrow range by construction here, "
              f"so this is expected to be near chance)")
        print(f"\nH25 verdict: {'HOLDS' if a_mu >= 0.8 else 'FALSIFIED'} "
              f"(threshold 0.8, observed {a_mu:.4f})")
    else:
        print("degenerate: band is pure pos or pure neg, no AUC defined")

    EFFS = sorted(set(r['effq'] for r in rows))
    print(f"\n{'eff':>5} {'N band':>7} {'rec':>7} | {'AUC mu':>8} "
          f"{'AUC nu_hat':>10}")
    for eff in EFFS:
        sub = bandsplit([r for r in rows if r['effq'] == eff], NU_LO, NU_HI)
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not sub:
            print(f"{eff:>5.2f} {'0':>7} {'-':>7} | {'(empty)':>8}")
            continue
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | {'(degenerate)':>8}")
            continue
        am = auc([r['mu'] for r in p], [r['mu'] for r in ng])
        an = auc([r['nuhat'] for r in p], [r['nuhat'] for r in ng])
        print(f"{eff:>5.2f} {len(sub):>7} "
              f"{str(len(p))+'/'+str(len(sub)):>7} | {am:>8.4f} {an:>10.4f}")

    print("\nThe pooled X1 number above mixes 5 eff strata whose mu SCALE\n"
          "differs (mu depends on n, K1, K2 -- W5/W6 already showed pooling\n"
          "across eff confounds mu). Two corrected pooled tests:")
    band_ex = [r for r in band if r['effq'] != 0.05]  # drop near-zero-failure stratum
    p2 = [r for r in band_ex if r['ok']]
    n2 = [r for r in band_ex if not r['ok']]
    a_mu_ex = auc([r['mu'] for r in p2], [r['mu'] for r in n2])
    print(f"  (a) drop degenerate eff=0.05 stratum (N={len(band_ex)}): "
          f"AUC(-mu) = {a_mu_ex:.4f}")
    zpos, zneg = [], []
    for eff in EFFS:
        sub = [r for r in band_ex if r['effq'] == eff]
        if not sub:
            continue
        mus = [r['mu'] for r in sub]
        mn, mx = min(mus), max(mus)
        for r in sub:
            z = (r['mu'] - mn) / (mx - mn) if mx > mn else 0.5
            (zpos if r['ok'] else zneg).append(z)
    a_mu_z = auc(zpos, zneg)
    print(f"  (b) within-stratum min-max normalised mu, pooled: "
          f"AUC = {a_mu_z:.4f}")
    verdict2 = a_mu_ex >= 0.8 or a_mu_z >= 0.8
    print(f"\nH25 verdict (scale-corrected): "
          f"{'HOLDS' if verdict2 else 'FALSIFIED'}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2 (secondary, W1b): does the head->2nd-block GS step predict")
    print("        recovery better than NU or mu?")
    print("-" * 78)
    print("step = log2(||b*_{m+1}||) - log2(||b*_1||) = log2(prof[m]/prof[0])")
    print("(prof is 0-indexed b*_1..b*_2m; dim k=2m so m = k//2)")

    for r in rows:
        m = r['k'] // 2
        p0, pm = r['prof'][0], r['prof'][m]
        r['step'] = (math.log2(pm) - math.log2(p0)) if p0 > 0 and pm > 0 else float('nan')

    bad = sum(1 for r in rows if r['step'] != r['step'])  # NaN check
    print(f"\n{bad}/{len(rows)} rows have a degenerate step (zero GS norm)")
    clean = [r for r in rows if r['step'] == r['step']]
    pos = [r for r in clean if r['ok']]
    neg = [r for r in clean if not r['ok']]
    a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
    a_nu_p = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
    a_mu_p = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
    print(f"\npooled (N={len(clean)}):")
    print(f"  AUC(-step -> recovery) = {a_step:.4f}")
    print(f"  AUC(-NU   -> recovery) = {a_nu_p:.4f}")
    print(f"  AUC(-mu   -> recovery) = {a_mu_p:.4f}")
    print(f"  Spearman(step, NU)  = {spearman([r['step'] for r in clean], [r['NU'] for r in clean]):.4f}")
    print(f"  Spearman(step, mu)  = {spearman([r['step'] for r in clean], [r['mu'] for r in clean]):.4f}")

    print(f"\n{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} "
          f"{'AUC NU':>8} {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in clean if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | {'(degenerate)':>9}")
            continue
        a_s = auc([r['step'] for r in p], [r['step'] for r in ng])
        a_n = auc([r['NU'] for r in p], [r['NU'] for r in ng])
        a_m = auc([r['mu'] for r in p], [r['mu'] for r in ng])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(p))+'/'+str(len(sub)):>7} | {a_s:>9.4f} "
              f"{a_n:>8.4f} {a_m:>8.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
