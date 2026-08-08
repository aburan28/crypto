"""
GLV-HNP Phase 2, Thread 25: does mu carry a second, NU-independent signal, and
does the GS "step" at the block boundary predict the wall directly?

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry.  W5/W6 found
that recovery = f(NU, X) with X ~ mu-driven and X uncorrelated with NU at
fixed eff (Spearman -0.28..+0.16 across strata).  Two questions follow:

  H25: within the ambiguous NU band [1.04, 2.20] (17-bit bracket, ties both
       hands of the nearest-plane certificate), AUC(-mu -> Kannan-LLL
       recovery) stays >= 0.8.  If yes, (NU, mu) is a genuine 2-parameter
       viability test.  If no, W5's mu signal is entirely mediated by NU and
       is a stratification artifact.

  Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||),
       i.e. the jump from the first-block plateau (m copies of lambda_1(L2))
       to the second block.  W1b observed this step visually vanishing right
       as the K1 wall is crossed.  Test AUC(-step -> recovery) against NU
       and mu, pooled and within the ambiguous band.

Reuses the exact same instance()/run_new() machinery and float-GS
justification (W0/W4: relative NU error ~1e-15 at dim 24) as
glv_hnp_phase2_gsprofile_strat.py, so results are directly comparable to W5.
Adds --dump-json (requested in the 2026-08-07 log) so the raw per-instance
table survives the run for later re-analysis.

Run: python3 glv_hnp_thread25_step.py [--dump-json out.json]
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


def collect(m=12, effs=(0.05, 0.10, 0.15, 0.20, 0.25), lo=1 << 16, hi=1 << 17):
    curves = search_curves(lo, hi, per_bin=2, nbins=10)
    rows = []
    for eff in effs:
        for (p, b, n, lam, G) in curves:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m, d_trial, k1b, seed)
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    return rows, curves


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else \
            "glv_hnp_thread25_rows.json"

    print("=" * 78)
    print("Thread 25 — does mu separate inside the NU ambiguous band, and "
          "does the block-step predict the wall?")
    print("=" * 78)

    t0 = time.time()
    M = 12
    rows, curves = collect(m=M)
    print(f"\n{len(curves)} 17-bit j=0 GLV curves; {len(rows)} instances "
          f"(dim {2*M}, float GS) in {time.time()-t0:.1f}s")

    if dump_path:
        with open(dump_path, "w") as f:
            json.dump([{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                       for r in rows], f)
        print(f"raw table dumped to {dump_path} ({len(rows)} rows, "
              f"prof/nus arrays excluded)")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25: AUC(-mu -> recovery) inside vs outside the ambiguous "
          "NU band [1.04, 2.20]  (17-bit bracket, Thread 24 W4)")
    print("-" * 78)
    LO, HI = 1.04, 2.20
    band = [r for r in rows if LO <= r['NU'] <= HI]
    out_lo = [r for r in rows if r['NU'] < LO]
    out_hi = [r for r in rows if r['NU'] > HI]
    print(f"band [{LO},{HI}]: N={len(band)}  (below: {len(out_lo)}, "
          f"above: {len(out_hi)})")
    print(f"  below-band recovery rate: "
          f"{sum(1 for r in out_lo if r['ok'])}/{len(out_lo) or 1}")
    print(f"  above-band recovery rate: "
          f"{sum(1 for r in out_hi if r['ok'])}/{len(out_hi) or 1}")

    def report(sub, name):
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{name}: N={len(sub)}  {len(pos)}/{len(sub)} recovered "
                  f"(degenerate, no AUC)")
            return None
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{name}: N={len(sub)}  rec {len(pos)}/{len(sub)}  "
              f"AUC(-mu)={a_mu:.4f}  AUC(-NU)={a_nu:.4f}  "
              f"AUC(-step)={a_st:.4f}")
        return a_mu, a_nu, a_st

    report(rows, "ALL")
    band_res = report(band, "IN-BAND")
    report(out_lo, "BELOW-BAND")
    report(out_hi, "ABOVE-BAND")

    if band_res:
        a_mu_band = band_res[0]
        verdict = "HOLDS" if a_mu_band >= 0.8 else "FALSIFIED"
        print(f"\nH25 verdict: AUC(-mu) in-band = {a_mu_band:.4f}  "
              f"({verdict}, threshold was >= 0.8000)")
    else:
        print("\nH25 verdict: band degenerate (all-pos or all-neg), "
              "cannot test")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP: per-eff-stratum AUC for mu, NU, step (matches W5 layout)")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC mu':>8} {'AUC NU':>8} "
          f"{'AUC step':>9}")
    EFFS = sorted(set(r['effq'] for r in rows))
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>7} | {a_mu:>8.4f} "
              f"{a_nu:>8.4f} {a_st:>9.4f}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP: step vs NU correlation (is step just re-reading NU?)")
    print("-" * 78)
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        sp_step_nu = spearman([r['step'] for r in sub], [r['NU'] for r in sub])
        sp_step_mu = spearman([r['step'] for r in sub], [r['mu'] for r in sub])
        print(f"eff={eff:.2f}  Spearman(step,NU)={sp_step_nu:+.4f}  "
              f"Spearman(step,mu)={sp_step_mu:+.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
