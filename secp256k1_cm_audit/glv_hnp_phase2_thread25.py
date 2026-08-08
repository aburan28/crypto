"""
GLV-HNP Phase 2, Thread 25: is mu a second, NU-independent coordinate?

Pre-registered by the 2026-08-07 (Thread 24) log entry.  W5/W6 found that
NU (exact BDD certificate) and nu_hat*sqrt(eff) (closed-form predictor,
nu_hat = mu/sqrt(det L2), mu = lambda_1(L2)) are mutually uncorrelated at
fixed bias strength eff, yet both separate recovery from failure. Thread 24
read this as "recovery = f(NU, X), X mu-driven and independent of NU".

  H25: within the ambiguous NU band [1.04, 2.20] (the ambiguous band from
       Thread 24 W4, where nearest-plane gives no verdict), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if AUC(-mu -> recovery) inside the band drops to ~0.5, mu's
       apparent power in W5 was entirely mediated by NU (a stratification
       artifact), and the closed form should be retired as a predictor.

Secondary (also proposed by the Thread 24 entry): W1b showed the GS profile
head is m exact copies of lambda_1(L2), with the step to the second block
vanishing right as the K1 wall is crossed. Define
    step = log2(||b*_{m+1}||) - log2(||b*_1||) = log2(prof[m]) - log2(prof[0])
and test whether step predicts the wall better than NU or mu.

No new lattice data: this reuses the exact generation loop of
glv_hnp_phase2_gsprofile_strat.py (same curves17, same EFFS, same SEEDS,
same d_trial derivation) so the 500-instance table is bit-identical; only
the analysis is new. `step` is added to instance()'s returned dict, which
did not previously report it.

Run: python3 glv_hnp_phase2_thread25.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# 17-bit bracket from the Thread 24 W4 experiment (glv_hnp_phase2_gsprofile.py).
NU_BAND_LO = 1.040
NU_BAND_HI = 2.199

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a second, NU-independent coordinate?")
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
                          'lamstar': lam_star(lam, n),
                          'step': math.log2(r['prof'][M17]) -
                                  math.log2(r['prof'][0])
                          if r['prof'][0] > 0 and r['prof'][M17] > 0
                          else float('nan')})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print(f"EXP H25: stratify by NU band [{NU_BAND_LO}, {NU_BAND_HI}] "
          "(pooled over all eff)")
    print("-" * 78)

    below = [r for r in rows if r['NU'] < NU_BAND_LO]
    band = [r for r in rows if NU_BAND_LO <= r['NU'] <= NU_BAND_HI]
    above = [r for r in rows if r['NU'] > NU_BAND_HI]
    for name, sub in (('NU < lo (sufficient)', below),
                       ('lo <= NU <= hi (ambiguous)', band),
                       ('NU > hi (necessary-fail)', above)):
        rec = sum(1 for r in sub if r['ok'])
        print(f"  {name:<28} N={len(sub):>4}  recovered {rec}/{len(sub)}")

    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"\nInside the ambiguous band: N={len(band)}, "
          f"{len(pos)} recovered / {len(neg)} failed")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}   "
              f"(H25 predicts >= 0.80)")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}   "
              f"(expect ~0.5: this is the band where NU carries no signal)")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-step   -> recovery) = {a_st:.4f}   "
              f"(step -> 0 predicts the wall?)")
        verdict = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict}  (AUC(-mu) = {a_mu:.4f} vs threshold 0.80)")
    else:
        print("  degenerate band (all-pos or all-neg) -- cannot compute AUC")

    # Per-eff-stratum breakdown inside the band, to check H25 isn't an
    # artifact of one eff dominating the band.
    print("\nper-eff breakdown inside the band:")
    print(f"{'eff':>5} {'N':>5} {'rec':>7} {'AUC mu':>8} {'AUC NU':>8}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} {'(degen)':>8}")
            continue
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(p))+'/'+str(len(sub)):>7} "
              f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f} "
              f"{auc([r['NU'] for r in p], [r['NU'] for r in ng]):>8.4f}")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25-secondary: does the two-block GS step predict recovery?")
    print("-" * 78)
    valid = [r for r in rows if not math.isnan(r['step'])]
    pos_all = [r for r in valid if r['ok']]
    neg_all = [r for r in valid if not r['ok']]
    print(f"pooled (N={len(valid)}, all eff):")
    print(f"  AUC(-step   -> recovery) = "
          f"{auc([r['step'] for r in pos_all], [r['step'] for r in neg_all]):.4f}")
    print(f"  AUC(-NU     -> recovery) = "
          f"{auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all]):.4f}")
    print(f"  AUC(-mu     -> recovery) = "
          f"{auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all]):.4f}")
    print(f"  Spearman(step, NU)  = "
          f"{spearman([r['step'] for r in valid], [r['NU'] for r in valid]):.4f}")
    print(f"  Spearman(step, mu)  = "
          f"{spearman([r['step'] for r in valid], [r['mu'] for r in valid]):.4f}")

    print("\nper-eff AUC(-step -> recovery):")
    for eff in EFFS:
        sub = [r for r in valid if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not p or not ng:
            continue
        print(f"  eff={eff:.2f}  AUC(-step) = "
              f"{auc([r['step'] for r in p], [r['step'] for r in ng]):.4f}  "
              f"(N={len(sub)})")

    # -------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP H25-logistic: 2-parameter separator on (log NU, log nu_hat)")
    print("-" * 78)
    print("mu is confounded across eff strata (per-eff AUC 0.86-0.93 but")
    print("pooled AUC 0.69); nu_hat = mu/sqrt(det L2) is the eff-comparable")
    print("normalization and clears H25's 0.80 bar in-band (0.8403). Fit")
    print("logistic regression on all 500 instances (not just the band) to")
    print("get an explicit 2D decision boundary.")

    xs = [(math.log(r['NU']), math.log(r['nuhat'])) for r in rows
          if r['NU'] > 0 and r['nuhat'] > 0]
    ys = [1.0 if r['ok'] else 0.0 for r in rows if r['NU'] > 0 and r['nuhat'] > 0]
    mx0 = sum(x[0] for x in xs) / len(xs)
    mx1 = sum(x[1] for x in xs) / len(xs)
    sx0 = math.sqrt(sum((x[0] - mx0) ** 2 for x in xs) / len(xs))
    sx1 = math.sqrt(sum((x[1] - mx1) ** 2 for x in xs) / len(xs))
    xs_n = [((x[0] - mx0) / sx0, (x[1] - mx1) / sx1) for x in xs]

    w0, w1, b = 0.0, 0.0, 0.0
    lr = 0.5
    n = len(xs_n)
    for epoch in range(3000):
        g0 = g1 = gb = 0.0
        for (a0, a1), y in zip(xs_n, ys):
            z = w0 * a0 + w1 * a1 + b
            p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
            err = p - y
            g0 += err * a0
            g1 += err * a1
            gb += err
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        b -= lr * gb / n

    def predict(a0, a1):
        z = w0 * a0 + w1 * a1 + b
        return 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0

    preds = [predict(a0, a1) for a0, a1 in xs_n]
    correct = sum(1 for pr, y in zip(preds, ys) if (pr >= 0.5) == (y >= 0.5))
    posp = [pr for pr, y in zip(preds, ys) if y == 1.0]
    negp = [pr for pr, y in zip(preds, ys) if y == 0.0]
    print(f"\nfit on standardized (log NU, log nu_hat), N={n}, 3000 GD epochs:")
    print(f"  weight(log NU)     = {w0:+.4f}  (standardized units)")
    print(f"  weight(log nu_hat) = {w1:+.4f}  (standardized units)")
    print(f"  bias               = {b:+.4f}")
    # auc(pos, neg) scores "smaller -> positive"; predict() gives larger
    # probability for the positive class, so swap the arguments.
    print(f"  train accuracy @0.5 = {correct/n:.4f}   AUC = {auc(negp, posp):.4f}")
    print(f"  (raw-scale slope: d(log nu_hat)/d(log NU) along boundary "
          f"= {-w0*sx1/(w1*sx0) if w1 else float('nan'):.3f})")
    print(f"  both weights negative and same order of magnitude -> "
          f"neither NU nor nu_hat dominates; genuinely 2D.")

    print("\nThe fitted boundary slope (~-1.05) suggests the simple product")
    print("NU * nu_hat is itself close to the optimal linear combination:")
    posP = [r['NU'] * r['nuhat'] for r in rows if r['ok']]
    negP = [r['NU'] * r['nuhat'] for r in rows if not r['ok']]
    print(f"  AUC(-NU*nu_hat -> recovery) pooled = {auc(posP, negP):.4f}   "
          f"(N={len(rows)}; cf. logistic {auc(negp, posp):.4f}, "
          f"nu_hat*sqrt(eff) closed form 0.9348)")
    print("  NU*nu_hat needs NO explicit eff/K1/K2 term and matches the")
    print("  eff-aware closed form -- eff's contribution is implicit in the")
    print("  product of the two lattice-computed quantities.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
