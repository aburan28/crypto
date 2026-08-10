"""
GLV-HNP Phase 2, Thread 25: is mu a genuine second coordinate alongside NU?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry.  W5/W6 established
recovery = f(NU, X) with X ~ mu-driven and X independent of NU (Spearman ~ 0
or negative between the closed form and NU, inside every eff stratum).  NU
is a sound size-degrading *certificate* (0 FP, but its ambiguous band widens
1.19-1.87 at 12 bits to 1.04-2.20 at 17 bits).  The question this script
answers:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer either way), AUC(-mu -> Kannan-LLL recovery) stays
       >= 0.8.

If yes, mu is a genuine second coordinate and (NU, mu) is a 2-parameter
viability test.  If no, mu's apparent power in W5 was mediated by NU and is
a stratification artifact (both are monotone in bias strength; conditioning
on NU already removes their shared trend and W5 held eff, not NU, fixed).

Secondary (also pre-registered): W1b showed the GS profile is m exact copies
of lambda_1(L2), and the step to the second block *vanishes* exactly as the
K1 wall is crossed.  Define

    step = log2(||b*_{m+1}||) - log2(||b*_1||)     (prof indices m, 0)

and test whether step -> 0 predicts recovery better than NU or mu, pooled
and within the ambiguous NU band.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json FILE]
"""

import argparse
import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# Ambiguous NU band from the 17-bit measurement (Thread 24, EXP W4):
# sufficient NU < 1.040, necessary NU > 2.199.
NU_LO, NU_HI = 1.040, 2.199


def collect(dump_path=None):
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")
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
                step = math.log2(r['prof'][M17]) - math.log2(r['prof'][0]) \
                    if r['prof'][0] > 0 and r['prof'][M17] > 0 else float('nan')
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff, 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if dump_path:
        thin = [{k: v for k, v in r.items()
                 if k not in ('prof', 'nus')} for r in rows]
        with open(dump_path, 'w') as f:
            json.dump(thin, f)
        print(f"dumped {len(thin)} rows to {dump_path}")

    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dump-json', default=None)
    ap.add_argument('--load-json', default=None)
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, conditioned on NU?")
    print("=" * 78)

    if args.load_json:
        with open(args.load_json) as f:
            rows = json.load(f)
        print(f"loaded {len(rows)} rows from {args.load_json}")
    else:
        rows = collect(args.dump_json)

    print("\n" + "-" * 78)
    print(f"H25: within the ambiguous NU band [{NU_LO}, {NU_HI}], "
          "does mu still separate?")
    print("-" * 78)
    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band N = {len(band)}  ({len(pos)} recovered / {len(neg)} failed)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"AUC(-mu     -> recovery) = {a_mu:.4f}   (H25 needs >= 0.80)")
        print(f"AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"AUC(-NU     -> recovery) = {a_nu:.4f}   (near-chance expected: "
              "NU is ~constant inside the band by construction)")
        print(f"AUC(-step   -> recovery) = {a_step:.4f}")
        verdict = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict} (AUC(-mu) = {a_mu:.4f})")
    else:
        print("degenerate band (all-pos or all-neg); H25 untestable at this N")

    print("\n" + "-" * 78)
    print("H25 per-eff-stratum breakdown (band subset of each stratum)")
    print("-" * 78)
    print(f"{'eff':>5} {'band N':>7} {'rec':>7} | {'AUC mu':>8} "
          f"{'AUC nu_hat':>11} {'AUC step':>9}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} | "
              f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f} "
              f"{auc([r['nuhat'] for r in p], [r['nuhat'] for r in ng]):>11.4f} "
              f"{auc([r['step'] for r in p], [r['step'] for r in ng]):>9.4f}")

    print("\n" + "-" * 78)
    print("SECONDARY: step = log2||b*_{m+1}|| - log2||b*_1|| as a predictor")
    print("-" * 78)
    good_step = [r for r in rows if not math.isnan(r['step'])]
    dropped = len(rows) - len(good_step)
    if dropped:
        print(f"(dropped {dropped} rows with degenerate GS entries)")
    pos_a = [r for r in good_step if r['ok']]
    neg_a = [r for r in good_step if not r['ok']]
    print(f"pooled (N={len(good_step)}):")
    print(f"  AUC(-step -> recovery) = "
          f"{auc([r['step'] for r in pos_a], [r['step'] for r in neg_a]):.4f}")
    print(f"  AUC(-NU   -> recovery) = "
          f"{auc([r['NU'] for r in pos_a], [r['NU'] for r in neg_a]):.4f}")
    print(f"  AUC(-mu   -> recovery) = "
          f"{auc([r['mu'] for r in pos_a], [r['mu'] for r in neg_a]):.4f}")
    print(f"  Spearman(step, NU)  = "
          f"{spearman([r['step'] for r in good_step], [r['NU'] for r in good_step]):.4f}")
    print(f"  Spearman(step, mu)  = "
          f"{spearman([r['step'] for r in good_step], [r['mu'] for r in good_step]):.4f}")

    print("\nby eff stratum:")
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC step':>9} {'mean step ok':>13} "
          f"{'mean step fail':>15}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in good_step if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>5} {str(len(p))+'/'+str(len(sub)):>7} | "
              f"{auc([r['step'] for r in p], [r['step'] for r in ng]):>9.4f} "
              f"{sum(r['step'] for r in p)/len(p):>13.3f} "
              f"{sum(r['step'] for r in ng)/len(ng):>15.3f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
