"""
GLV-HNP Phase 2, Thread 25: does mu = lambda_1(L2) separate INSIDE the NU
ambiguous band, and does the W1b GS-profile "step" predict the wall better
than NU or mu alone?

Background (2026-08-07 log, Thread 24 W5/W6): NU (exact Babai/BDD
certificate) and nu_hat*sqrt(eff) (closed-form, no lattice work) are
uncorrelated at fixed eff, yet both separate recovery. W4 gives a 17-bit
NU bracket: sufficient NU < 1.040, necessary NU > 2.199 — i.e. Babai
nearest-plane gives no verdict for NU in [1.040, 2.199], but Kannan-LLL
still recovers d for many instances in that band (recovery is BDD, not
just nearest-plane -- Thread 23b already showed Kannan beats Babai by
~1.9x in NU).

H25: within the ambiguous band 1.04 <= NU <= 2.20, AUC(-mu -> recovery)
     stays >= 0.8, i.e. mu is a genuine SECOND coordinate independent of
     the NU verdict.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), the
jump from the m-fold-repeated head block (all == lambda_1(L2)) to the
second block.  Test whether step -> 0 predicts the wall better than NU
or mu.

Run: python3 glv_hnp_phase2_gsprofile_strat.py --dump-json rows.json
     python3 glv_hnp_phase2_thread25.py rows.json
"""

import json
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from glv_hnp_phase2_gsprofile import auc, spearman

NU_LO, NU_HI = 1.040, 2.199  # W4 17-bit ambiguous band

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "thread25_rows.json"
    with open(path) as f:
        rows = json.load(f)
    print("=" * 78)
    print("Thread 25 — mu inside the NU band; W1b step vs. the wall")
    print("=" * 78)
    print(f"\n{len(rows)} rows loaded from {path}")

    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) WITHIN the NU ambiguous band "
          f"[{NU_LO}, {NU_HI}]")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    print(f"band population: {len(band)} / {len(rows)} instances "
          f"({sum(1 for r in band if r['ok'])}/{len(band)} recover)")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"  eff={eff:.2f}  N={len(sub):>3}  "
                  f"{len(pos)}/{len(sub)} rec  (degenerate)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        print(f"  eff={eff:.2f}  N={len(sub):>3}  {len(pos)}/{len(sub)} rec  "
              f"AUC(-mu)={a_mu:.4f}  AUC(-nu_hat)={a_nh:.4f}")

    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    if pos and neg:
        a_mu_pool = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh_pool = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        print(f"\npooled over band (N={len(band)}): "
              f"AUC(-mu)={a_mu_pool:.4f}  AUC(-nu_hat)={a_nh_pool:.4f}")
        verdict = "HOLDS" if a_mu_pool >= 0.8 else "FALSIFIED"
        print(f"H25 (AUC(-mu) >= 0.8 inside band): {verdict} "
              f"({a_mu_pool:.4f} {'>=' if a_mu_pool >= 0.8 else '<'} 0.80)")
    else:
        print("\nband is degenerate (all-pos or all-neg) — H25 untestable "
              "pooled; per-stratum numbers above are the only signal.")

    print("\n" + "-" * 78)
    print("Secondary: W1b step = log2||b*_{m+1}|| - log2||b*_1|| vs. the wall")
    print("-" * 78)
    have_step = [r for r in rows if not math.isnan(r.get('step', float('nan')))]
    print(f"{len(have_step)}/{len(rows)} rows have a finite step value")
    print("NOTE: unlike NU/mu/nu_hat, the wall hypothesis predicts step is")
    print("small NEAR failure and large near recovery -- opposite sign. The")
    print("'AUC(+step)' column flips convention so >0.5 also means 'more")
    print("signal', matching the other columns.")
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC(-step)':>10} {'AUC(+step)':>10} "
          f"{'AUC NU':>8} {'AUC mu':>8}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in have_step if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_st:>10.4f} "
              f"{1 - a_st:>10.4f} {a_nu:>8.4f} {a_mu:>8.4f}")

    pooled_pos = [r for r in have_step if r['ok']]
    pooled_neg = [r for r in have_step if not r['ok']]
    a_st_pool = auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg])
    print(f"\npooled AUC(-step) = {a_st_pool:.4f}  ->  AUC(+step) = {1 - a_st_pool:.4f}  "
          f"(recall pooled AUC(-NU) = 0.7996)")

    sp_step_nu = spearman([r['step'] for r in have_step],
                           [r['NU'] for r in have_step])
    sp_step_mu = spearman([r['step'] for r in have_step],
                           [r['mu'] for r in have_step])
    print(f"Spearman(step, NU) = {sp_step_nu:.4f}   "
          f"Spearman(step, mu) = {sp_step_mu:.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
