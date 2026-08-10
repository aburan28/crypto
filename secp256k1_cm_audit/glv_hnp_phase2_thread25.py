"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry.  W5/W6 there
established recovery = f(NU, X) with X ~ mu-driven (mu = lambda_1(L2)) and X
uncorrelated with NU at fixed bias strength eff.  NU alone is a sound
certificate (0 FP at NU<=1) but a size-degrading separator; mu is a better
cross-curve ranker (AUC 0.75-0.93) that NU cannot explain.

  H25: within the ambiguous NU band where nearest-plane gives no answer
       (ballpark [1.04, 2.20] from Thread 24 W4 at 17 bits), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if AUC(-mu) inside the band collapses toward 0.5, mu's power in
       W5 was mediated by NU (a stratification artifact) and should be
       retired as a predictor.

Secondary (Thread 24 W1b): the GS profile of L0 is m exact copies of
lambda_1(L2) followed by a K1-dependent second block; the step between them
vanishes right at the LLL wall.  Test whether
    step = log2(||b*_{m}||) - log2(||b*_{0}||)      (0-indexed, second block
                                                       starts at index m)
predicts recovery better than NU or mu.

Uses the same float-GS instance() generator as Thread 24b (justified by W0:
relative error ~1e-15 at dim 24), same 5-strata x 20-curve x 5-seed 17-bit
grid, but dumps every row to JSON so re-analysis doesn't require a re-run.

Run: python3 glv_hnp_phase2_thread25.py [--dump-json path]
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


def logistic_fit_2d(xs, ys, labels, iters=4000, lr=0.5):
    """Plain gradient-descent logistic regression on 2 features + bias.
    xs, ys: feature vectors (already logged).  labels: 0/1 (1 = recovered).
    Returns (w0, w1, w2) for  sigmoid(w0 + w1*x + w2*y)."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n) or 1.0
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n) or 1.0
    zx = [(x - mx) / sx for x in xs]
    zy = [(y - my) / sy for y in ys]
    w0 = w1 = w2 = 0.0
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for x, y, lab in zip(zx, zy, labels):
            z = w0 + w1 * x + w2 * y
            p = 1.0 / (1.0 + math.exp(-z)) if z > -60 else 0.0
            err = p - lab
            g0 += err
            g1 += err * x
            g2 += err * y
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    # Un-standardize: w0 + w1*(x-mx)/sx + w2*(y-my)/sy
    a1 = w1 / sx
    a2 = w2 / sy
    a0 = w0 - w1 * mx / sx - w2 * my / sy
    return a0, a1, a2


def logistic_auc(a0, a1, a2, xs, ys, labels):
    scores = [a0 + a1 * x + a2 * y for x, y in zip(xs, ys)]
    pos = [s for s, lab in zip(scores, labels) if lab == 1]
    neg = [s for s, lab in zip(scores, labels) if lab == 0]
    # higher score -> more likely recovery, so auc() (which assumes smaller
    # is better) needs negation
    return auc([-s for s in pos], [-s for s in neg])


if __name__ == "__main__":
    dump_path = None
    if "--dump-json" in sys.argv:
        dump_path = sys.argv[sys.argv.index("--dump-json") + 1]

    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a genuine second coordinate?")
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
                m = M17
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    if dump_path:
        keep = ['n', 'K1', 'eff', 'effq', 'ok', 'NU', 'mu', 'nuhat',
                'lamstar', 'step', 'seed' if 'seed' in rows[0] else 'K1']
        with open(dump_path, 'w') as f:
            json.dump([{k: r.get(k) for k in
                        ('n', 'K1', 'eff', 'effq', 'ok', 'NU', 'mu', 'nuhat',
                         'lamstar', 'step')} for r in rows], f)
        print(f"dumped {len(rows)} rows to {dump_path}")

    print("\n" + "-" * 78)
    print("H25: pick the ambiguous NU band from Thread 24 W4 (17-bit bracket)")
    print("     sufficient NU < 1.040 , necessary NU > 2.199")
    print("-" * 78)
    band = [r for r in rows if 1.040 <= r['NU'] <= 2.199]
    outside_easy = [r for r in rows if r['NU'] < 1.040]
    outside_hard = [r for r in rows if r['NU'] > 2.199]
    print(f"in-band: {len(band)}   below (NU<1.040, should all recover): "
          f"{sum(1 for r in outside_easy if r['ok'])}/{len(outside_easy)}   "
          f"above (NU>2.199, should all fail): "
          f"{sum(1 for r in outside_hard if r['ok'])}/{len(outside_hard)}")

    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"in-band recovery: {len(pos)}/{len(band)}")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu     -> recovery) inside band = {a_mu:.4f}   "
              f"(H25 threshold: >= 0.80)")
        print(f"  AUC(-nu_hat -> recovery) inside band = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) inside band = {a_nu:.4f}  "
              f"(expect ~0.5: NU is constant-ish inside its own band)")
        print(f"  AUC(-step   -> recovery) inside band = {a_step:.4f}")
        verdict = "CONFIRMED" if a_mu >= 0.80 else "FALSIFIED"
        print(f"\nH25 verdict: {verdict} (AUC={a_mu:.4f})")
    else:
        print("  degenerate band (all-same label) -- cannot compute AUC")

    print("\n" + "-" * 78)
    print("Per-eff-stratum breakdown of the same in-band AUC(-mu), to check")
    print("the H25 result isn't itself just re-deriving the eff-stratified")
    print("mu column from Thread 24 W5.")
    print("-" * 78)
    print(f"{'eff':>5} {'N band':>7} {'rec':>7} {'AUC -mu (band)':>15}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>15.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} {'degenerate':>7}")

    print("\n" + "-" * 78)
    print("Secondary (W1b step metric): AUC(-step -> recovery), pooled and")
    print("compared to NU and mu, over the FULL 500-instance table.")
    print("-" * 78)
    allpos = [r for r in rows if r['ok'] and not math.isnan(r['step'])]
    allneg = [r for r in rows if not r['ok'] and not math.isnan(r['step'])]
    print(f"AUC(-step -> recovery), full table = "
          f"{auc([r['step'] for r in allpos], [r['step'] for r in allneg]):.4f}")
    print(f"AUC(-NU   -> recovery), full table = "
          f"{auc([r['NU'] for r in allpos], [r['NU'] for r in allneg]):.4f}")
    print(f"AUC(-mu   -> recovery), full table = "
          f"{auc([r['mu'] for r in allpos], [r['mu'] for r in allneg]):.4f}")
    print(f"Spearman(step, NU)  = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"Spearman(step, mu)  = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    if pos and neg and len(band) >= 20:
        print("\n" + "-" * 78)
        print("Logistic fit on (log NU, log mu), full 500-instance table")
        print("-" * 78)
        xs = [math.log(r['NU']) for r in rows]
        ys = [math.log(r['mu']) for r in rows]
        labels = [1 if r['ok'] else 0 for r in rows]
        a0, a1, a2 = logistic_fit_2d(xs, ys, labels)
        auc2d = logistic_auc(a0, a1, a2, xs, ys, labels)
        auc_nu_only = auc([r['NU'] for r in rows if r['ok']],
                           [r['NU'] for r in rows if not r['ok']])
        print(f"logit(recover) = {a0:.4f} + ({a1:.4f})*log(NU) + "
              f"({a2:.4f})*log(mu)")
        print(f"AUC(2D logistic score) = {auc2d:.4f}   "
              f"vs AUC(NU alone) = {auc_nu_only:.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
