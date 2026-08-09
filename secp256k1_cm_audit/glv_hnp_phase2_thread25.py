"""
GLV-HNP Phase 2, Thread 25: is mu a second coordinate, or is its apparent
power in W5 entirely mediated by NU?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  H25: within the ambiguous NU band [1.04, 2.20] (17-bit bracket from W4,
       where the nearest-plane certificate gives no answer either way),
       AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate, and (NU, mu) is the missing
  2-parameter viability test.
  If no: mu's apparent power in W5 (AUC 0.75-0.93 per stratum, eff fixed)
  is entirely mediated by NU, and the closed form nu_hat should be retired
  as "NU in disguise, just restricted to the band where it doesn't work."

Secondary (also proposed 2026-08-07 #2): does the GS-profile "step"
  step_i = log2(||b*_{m+1}||) - log2(||b*_1||)   (0-based: prof[m] vs prof[0])
predict the wall better than NU or mu?  W1b observed the step vanishes near
the K1 wall; this tests whether that vanishing is itself a viable predictor.

Data source: glv_hnp_phase2_gsprofile_strat.py --dump-json, run today
(500 rows, 17-bit curves, dim 24, 5 eff strata x 20 curves x 5 seeds).
Re-run: python3 glv_hnp_phase2_gsprofile_strat.py --dump-json
        python3 glv_hnp_phase2_thread25.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

ROWS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "glv_hnp_phase2_gsprofile_strat_rows.json")

# 17-bit NU bracket from Thread 24 W4: sufficient NU < 1.040, necessary
# NU > 2.199.  Ambiguous band is the closed interval between them.
NU_LO, NU_HI = 1.040, 2.199


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (float('nan'), float('nan'))
    ph = k / n
    denom = 1 + z * z / n
    center = (ph + z * z / (2 * n)) / denom
    half = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / denom
    return (center - half, center + half)


if __name__ == "__main__":
    with open(ROWS_PATH) as f:
        d = json.load(f)
    M = d["M"]
    rows = d["rows"]
    for r in rows:
        r["step"] = math.log2(r["prof"][M]) - math.log2(r["prof"][0])
    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a second coordinate?")
    print("=" * 78)
    print(f"loaded {len(rows)} rows (dim {rows[0]['k']}) from {ROWS_PATH}")

    print("\n" + "-" * 78)
    print(f"H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{NU_LO}, {NU_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r["NU"] <= NU_HI]
    print(f"band N = {len(band)} / {len(rows)} "
          f"({100*len(band)/len(rows):.1f}% of all instances)")
    pos = [r for r in band if r["ok"]]
    neg = [r for r in band if not r["ok"]]
    print(f"  recovered {len(pos)} / not {len(neg)}")
    if pos and neg:
        a_mu = auc([r["mu"] for r in pos], [r["mu"] for r in neg])
        a_nh = auc([r["nuhat"] for r in pos], [r["nuhat"] for r in neg])
        a_nu = auc([r["NU"] for r in pos], [r["NU"] for r in neg])
        a_st = auc([r["step"] for r in pos], [r["step"] for r in neg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}   "
              f"[H25 threshold: 0.80]  -> {'HOLDS' if a_mu >= 0.80 else 'FALSIFIED'}")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}   "
              f"(should be ~0.5: NU is constant-ish inside its own band)")
        print(f"  AUC(-step   -> recovery) = {a_st:.4f}")
    else:
        print("  degenerate band (all-pos or all-neg) -- cannot compute AUC")

    print("\nper-stratum breakdown inside the band (eff should still vary "
          "the mix, but mu is compared within a fixed eff too):")
    print(f"{'eff':>5} {'N':>4} {'rec':>6} | {'AUC mu':>8} {'AUC step':>9}")
    for eff in sorted(set(r["effq"] for r in band)):
        sub = [r for r in band if r["effq"] == eff]
        p = [r for r in sub if r["ok"]]
        n = [r for r in sub if not r["ok"]]
        if p and n:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>6} | "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in n]):>8.4f} "
                  f"{auc([r['step'] for r in p], [r['step'] for r in n]):>9.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>6} | (degenerate)")

    print("\n" + "-" * 78)
    print("Control: is NU itself still informative *inside its own band*?")
    print("(it must be close to 0.5 -- if NU still separates here, the band")
    print(" was not actually ambiguous and H25 is not a clean test of mu)")
    print("-" * 78)
    if pos and neg:
        lo, hi = wilson_ci(sum(1 for r in band if r["ok"]), len(band))
        print(f"P(recover | band) = {sum(1 for r in band if r['ok'])}/{len(band)} "
              f"= {sum(1 for r in band if r['ok'])/len(band):.3f}  "
              f"(95% CI [{lo:.3f}, {hi:.3f}])")
        print(f"Spearman(NU, mu) inside band = "
              f"{spearman([r['NU'] for r in band], [r['mu'] for r in band]):.4f}")

    print("\n" + "-" * 78)
    print("SECONDARY: does the GS-profile step predict the wall better than "
          "NU or mu?")
    print(f"step = log2(||b*_{{{M}}}||) - log2(||b*_0||)  (0-based; block "
          "boundary crossing)")
    print("-" * 78)
    allpos = [r for r in rows if r["ok"]]
    allneg = [r for r in rows if not r["ok"]]
    print("pooled over all 500 rows (5 strata mixed -- same caveat as W3, "
          "reported for reference only):")
    print(f"  AUC(-step) = "
          f"{auc([r['step'] for r in allpos], [r['step'] for r in allneg]):.4f}   "
          f"AUC(-NU) = {auc([r['NU'] for r in allpos], [r['NU'] for r in allneg]):.4f}   "
          f"AUC(-mu) = {auc([r['mu'] for r in allpos], [r['mu'] for r in allneg]):.4f}")

    print("\nper-eff-stratum (honest test, eff fixed):")
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} {'AUC NU':>8} "
          f"{'AUC mu':>8}")
    for eff in sorted(set(r["effq"] for r in rows)):
        sub = [r for r in rows if r["effq"] == eff]
        p = [r for r in sub if r["ok"]]
        n = [r for r in sub if not r["ok"]]
        if p and n:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | "
                  f"{auc([r['step'] for r in p], [r['step'] for r in n]):>9.4f} "
                  f"{auc([r['NU'] for r in p], [r['NU'] for r in n]):>8.4f} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in n]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | (degenerate)")

    print("\nstep value ranges (sanity: does step actually approach 0 near "
          "the wall as W1b claimed?):")
    for eff in sorted(set(r["effq"] for r in rows)):
        sub = [r for r in rows if r["effq"] == eff]
        steps = [r["step"] for r in sub]
        print(f"  eff={eff:.2f}: step mean {sum(steps)/len(steps):.3f}  "
              f"min {min(steps):.3f}  max {max(steps):.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
