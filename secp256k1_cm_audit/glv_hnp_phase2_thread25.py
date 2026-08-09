"""
GLV-HNP Phase 2, Thread 25: is mu a genuine second coordinate, conditional on NU?

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry.

W5/W6 (glv_hnp_phase2_gsprofile_strat.py) established recovery = f(NU, X)
with X ~ mu-driven (AUC 0.75-0.93 per eff-stratum) and X uncorrelated with NU
(Spearman ~ -0.2..+0.2 inside every stratum). NU alone is a sound BDD
certificate (0 FP at NU<=1) but a size-degrading separator (AUC 0.978 at 12
bits -> 0.860 at 17 bits), with an ambiguous band [1.040, 2.199] at 17 bits
where nearest-plane gives no answer either way.

  H25: within the ambiguous band 1.04 <= NU <= 2.20, AUC(-mu -> Kannan-LLL
       recovery) stays >= 0.8.
  Falsifier: AUC(-mu) inside the band collapses to ~0.5 (or lam*'s territory,
       <0.5), meaning mu's apparent power in W5 was mediated by NU after all
       and the per-stratum AUC was a stratification artifact.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), the
jump from the head block (m copies of lambda_1(L2)) to the tail block. Test
whether step predicts recovery better than NU or mu.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json path.json]
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

BAND_LO, BAND_HI = 1.040, 2.199  # 17-bit ambiguous band, Thread 24 W4


def collect_rows():
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    M17 = 12
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)

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
                          'lamstar': lam_star(lam, n),
                          'step': math.log2(r['prof'][M17]) - math.log2(r['prof'][0])
                          if r['prof'][0] > 0 and r['prof'][M17] > 0 else float('nan')})
                rows.append(r)
    return rows, M17


def dump_json(rows, path):
    keys = ('n', 'K1', 'eff', 'effq', 'ok', 'NU', 'mu', 'nuhat', 'lamstar',
            'step', 'argmax', 'k')
    with open(path, 'w') as f:
        json.dump([{k: r[k] for k in keys} for r in rows], f)


if __name__ == "__main__":
    dump_path = None
    if '--dump-json' in sys.argv:
        dump_path = sys.argv[sys.argv.index('--dump-json') + 1]

    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, conditional on NU? (H25)")
    print("=" * 78)

    t0 = time.time()
    rows, M17 = collect_rows()
    print(f"\n{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    if dump_path:
        dump_json(rows, dump_path)
        print(f"dumped {len(rows)} rows to {dump_path}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print(f"EXP H25: AUC(-mu -> recovery) inside the ambiguous NU band "
          f"[{BAND_LO}, {BAND_HI}]")
    print("-" * 78)
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band population: N={len(band)}  rec={len(pos)}/{len(band)}  "
          f"(pooled over all strata; NU is fixed-ish, mu is free)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}   <- H25 target (>=0.8?)")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-lam*   -> recovery) = {a_ls:.4f}   (control)")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}   (should be ~0.5, NU is fixed by construction)")
        verdict = "CONFIRMED" if a_mu >= 0.8 else "FALSIFIED"
        print(f"\nH25: {verdict}  (AUC(-mu)={a_mu:.4f} vs threshold 0.8)")
    else:
        print("degenerate band (all-pos or all-neg); cannot compute AUC")

    print("\nper-stratum breakdown inside the band (does mu still separate "
          "when eff is ALSO controlled?):")
    print(f"{'eff':>5} {'N':>4} {'rec':>7} {'AUC mu':>8}")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in band if r['effq'] == eff]
        p2 = [r for r in sub if r['ok']]
        n2 = [r for r in sub if not r['ok']]
        if p2 and n2:
            print(f"{eff:>5.2f} {len(sub):>4} {str(len(p2))+'/'+str(len(sub)):>7} "
                  f"{auc([r['mu'] for r in p2], [r['mu'] for r in n2]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>4} {str(len(p2))+'/'+str(len(sub)):>7} "
                  f"{'(degen)':>8}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP secondary: step = log2||b*_{m+1}|| - log2||b*_1|| vs NU, mu")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    posS = [r for r in valid if r['ok']]
    negS = [r for r in valid if not r['ok']]
    a_step = auc([r['step'] for r in posS], [r['step'] for r in negS])
    a_nu_all = auc([r['NU'] for r in posS], [r['NU'] for r in negS])
    a_mu_all = auc([r['mu'] for r in posS], [r['mu'] for r in negS])
    print(f"pooled over all 5 strata (N={len(valid)}):")
    print(f"  AUC(-step -> recovery) = {a_step:.4f}")
    print(f"  AUC(-NU   -> recovery) = {a_nu_all:.4f}")
    print(f"  AUC(-mu   -> recovery) = {a_mu_all:.4f}")
    print(f"  Spearman(step, NU) = {spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
    print(f"  Spearman(step, mu) = {spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")

    print("\nper-stratum AUC(-step):")
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25):
        sub = [r for r in valid if r['effq'] == eff]
        p2 = [r for r in sub if r['ok']]
        n2 = [r for r in sub if not r['ok']]
        if p2 and n2:
            print(f"  eff={eff:.2f}  AUC(-step) = "
                  f"{auc([r['step'] for r in p2], [r['step'] for r in n2]):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
