"""
GLV-HNP Phase 2, Thread 25: condition on NU and ask whether mu is a genuine
second coordinate, or whether its apparent power (W5/W6 of
glv_hnp_phase2_gsprofile_strat.py) is entirely mediated by NU.

Background (2026-08-07 #2 log entry): NU (exact BDD certificate, sound
nearest-plane bracket) and nu_hat*sqrt(eff) (closed-form, no lattice
reduction) are mutually uncorrelated at fixed eff (Spearman ~ -0.2..+0.2)
yet both predict Kannan-LLL recovery. W4 gives a sound-but-loose NU bracket:
  17 bits: NU <= 1.040 => always recovers (0 FP); NU > 2.199 => never recovers.
  1.040 < NU <= 2.199 is the ambiguous band where nearest-plane gives no answer.

H25 (pre-registered): within the ambiguous band, AUC(-mu -> recovery) >= 0.8.
  If yes: mu is a genuine second coordinate; (NU, mu) is a 2-parameter
    viability test. Deliverable: logistic fit on (log NU, log mu).
  If no: mu's apparent power in W5 is entirely mediated by NU (i.e. NU and mu
    are themselves correlated enough that conditioning removes mu's signal),
    and the closed form should be retired as a *within-band* predictor.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), i.e.
the log-gap between the first "second block" GS norm and the profile head
(which W1b showed sits exactly at lambda_1(L2)). Test whether step -> 0
predicts the wall better than NU or mu.

Run: python3 glv_hnp_phase2_nuband.py [--dump-json out.json]
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

NU_LO, NU_HI = 1.040, 2.199  # W4's 17-bit ambiguous band (log line ~6300)


def logistic_fit_2d(xs, ys, labels, iters=4000, lr=0.3):
    """Plain gradient-descent logistic regression on 2 standardized features.
    No external deps (sklearn/numpy not assumed available)."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n) or 1.0
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n) or 1.0
    X = [((x - mx) / sx, (y - my) / sy) for x, y in zip(xs, ys)]
    w0, w1, w2 = 0.0, 0.0, 0.0
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for (x1, x2), lab in zip(X, labels):
            z = w0 + w1 * x1 + w2 * x2
            p = 1.0 / (1.0 + math.exp(-z)) if z > -50 else 0.0
            err = p - lab
            g0 += err
            g1 += err * x1
            g2 += err * x2
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    correct = 0
    for (x1, x2), lab in zip(X, labels):
        z = w0 + w1 * x1 + w2 * x2
        p = 1.0 / (1.0 + math.exp(-z)) if z > -50 else (0.0 if z < 0 else 1.0)
        correct += (p >= 0.5) == (lab >= 0.5)
    return (w0, w1, w2), (mx, my, sx, sy), correct / n


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        dump_path = sys.argv[sys.argv.index("--dump-json") + 1]

    print("=" * 78)
    print("Thread 25 — condition on NU: is mu a genuine second coordinate?")
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
                m = r['k'] // 2
                step = math.log2(r['prof'][m]) - math.log2(r['prof'][0]) \
                    if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan')
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if dump_path:
        with open(dump_path, "w") as f:
            json.dump([{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                       for r in rows], f)
        print(f"dumped {len(rows)} rows to {dump_path}")

    print("\n" + "-" * 78)
    print(f"H25: within NU band [{NU_LO}, {NU_HI}] (ambiguous, no nearest-")
    print("     plane answer), does mu still separate recovery?")
    print("-" * 78)

    band = [r for r in rows if NU_LO < r['NU'] <= NU_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band size: {len(band)} / {len(rows)} total "
          f"({len(pos)} recover, {len(neg)} fail)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu_within = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu     -> recovery | NU in band) = {a_mu:.4f}")
        print(f"  AUC(-nu_hat -> recovery | NU in band) = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery | NU in band) = {a_nu_within:.4f}"
              "   (sanity: NU restricted to its own ambiguous range)")
        print(f"  AUC(-step   -> recovery | NU in band) = {a_step:.4f}")
        verdict = "HOLDS" if a_mu >= 0.8 else "FALSIFIED"
        print(f"\n  H25 verdict: {verdict}  (threshold 0.8, got {a_mu:.4f})")
    else:
        print("  degenerate band (all-recover or all-fail); H25 untestable here")

    print("\n" + "-" * 78)
    print("Spearman(NU, mu) and Spearman(NU, log mu) — is mu just NU in")
    print("disguise, globally and within-band?")
    print("-" * 78)
    print(f"  global:      Spearman(NU, mu)     = "
          f"{spearman([r['NU'] for r in rows], [r['mu'] for r in rows]):.4f}")
    if band:
        print(f"  within band: Spearman(NU, mu)     = "
              f"{spearman([r['NU'] for r in band], [r['mu'] for r in band]):.4f}")

    print("\n" + "-" * 78)
    print("Per-eff-stratum AUC(-mu) restricted to the NU ambiguous band")
    print("-" * 78)
    print(f"{'eff':>5} {'band N':>7} {'rec':>7} {'AUC mu':>8} {'AUC step':>9}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p_ = [r for r in sub if r['ok']]
        n_ = [r for r in sub if not r['ok']]
        if not p_ or not n_:
            print(f"{eff:>5.2f} {len(sub):>7} "
                  f"{str(len(p_))+'/'+str(len(sub)):>7} {'(degenerate)':>8}")
            continue
        a_mu = auc([r['mu'] for r in p_], [r['mu'] for r in n_])
        a_st = auc([r['step'] for r in p_], [r['step'] for r in n_])
        print(f"{eff:>5.2f} {len(sub):>7} "
              f"{str(len(p_))+'/'+str(len(sub)):>7} {a_mu:>8.4f} {a_st:>9.4f}")

    print("\n" + "-" * 78)
    print("Secondary: step = log2||b*_{m+1}|| - log2||b*_1|| over ALL rows")
    print("(not just the band) -- does step -> 0 predict the K1 wall?")
    print("-" * 78)
    pooled_pos = [r for r in rows if r['ok']]
    pooled_neg = [r for r in rows if not r['ok']]
    a_step_all = auc([r['step'] for r in pooled_pos], [r['step'] for r in pooled_neg])
    a_NU_all = auc([r['NU'] for r in pooled_pos], [r['NU'] for r in pooled_neg])
    a_mu_all = auc([r['mu'] for r in pooled_pos], [r['mu'] for r in pooled_neg])
    print(f"  AUC(-step -> recovery), pooled N={len(rows)}: {a_step_all:.4f}")
    print(f"  AUC(-NU   -> recovery), pooled:               {a_NU_all:.4f}")
    print(f"  AUC(-mu   -> recovery), pooled:                {a_mu_all:.4f}")
    print(f"  Spearman(step, NU) = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu) = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    if pos and neg and band:
        print("\n" + "-" * 78)
        print("Logistic fit on (log NU, log mu) over the FULL 17-bit table")
        print("(deliverable requested by the 2026-08-07 #2 log entry)")
        print("-" * 78)
        xs = [math.log(r['NU']) for r in rows]
        ys = [math.log(r['mu']) for r in rows if r['mu'] > 0]
        # guard: mu should always be > 0 (lambda_1 of a nondegenerate lattice)
        assert len(ys) == len(rows)
        labels = [1.0 if r['ok'] else 0.0 for r in rows]
        (w0, w1, w2), (mx, my, sx, sy), acc = logistic_fit_2d(xs, ys, labels)
        print(f"  standardized weights: w0={w0:.4f} w_logNU={w1:.4f} "
              f"w_logmu={w2:.4f}")
        print(f"  feature means/stds: logNU ~ N({mx:.4f},{sx:.4f})  "
              f"logmu ~ N({my:.4f},{sy:.4f})")
        print(f"  training accuracy @ 0.5 threshold: {acc:.4f}  (N={len(rows)})")
        # AUC of the fitted score itself, as a check that 2D beats either 1D
        score = []
        for x, y in zip(xs, ys):
            z = w0 + w1 * (x - mx) / sx + w2 * (y - my) / sy
            score.append(z)
        pos_s = [s for s, r in zip(score, rows) if r['ok']]
        neg_s = [s for s, r in zip(score, rows) if not r['ok']]
        # score is "recovery propensity" so higher = more likely to recover;
        # auc() expects "smaller = more likely", so negate.
        a_fit = auc([-s for s in pos_s], [-s for s in neg_s])
        print(f"  AUC of fitted (log NU, log mu) score: {a_fit:.4f}  "
              f"(cf. AUC(-NU) pooled = {a_NU_all:.4f}, "
              f"AUC(-mu) pooled = {a_mu_all:.4f})")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
