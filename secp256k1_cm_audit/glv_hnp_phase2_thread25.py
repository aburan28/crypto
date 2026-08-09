"""
GLV-HNP Phase 2, Thread 25: does mu = lambda_1(L2) separate recovery inside
the ambiguous NU band, i.e. is (NU, mu) a genuine 2-parameter viability test?

Pre-registered by the 2026-08-07 (Thread 24 #2) log entry, following W5/W6:
W5/W6 showed NU (exact BDD certificate) and nu_hat*sqrt(eff) (closed form,
~ mu) are mutually uncorrelated at fixed eff and both predict recovery, so
recovery = f(NU, X) with X mu-driven and independent of NU.

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the NU <= 1
       sufficient certificate and the NU > 2.20 necessary bound both give no
       answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

Falsifier: if AUC(-mu) inside the band drops toward 0.5, mu's apparent power
in W5 is entirely mediated by NU (a stratification artifact along eff), and
the closed form should be retired as a redundant restatement of NU.

Secondary (W1b follow-up): step = log2(||b*_{m+1}||) - log2(||b*_1||), the
GS-profile jump from the first (lambda_1(L2)-flat) block to the second. Test
whether step separates inside the band better than NU or mu alone.

This regenerates the exact same 17-bit / 500-instance grid as
glv_hnp_phase2_gsprofile_strat.py (same curves17 search, same EFFS, same
SEEDS, same M17=12) rather than depending on that script's stdout, since the
prior run did not persist per-instance rows to disk.

Run: python3 glv_hnp_phase2_thread25.py
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

BAND_LO, BAND_HI = 1.040, 2.199  # 17-bit sufficient/necessary bracket, W4


def logistic_fit_2d(xs, ys, labels, iters=20000, lr=0.05):
    """Plain gradient-descent logistic regression on 2 features (+ bias).
    xs, ys: feature vectors (already in the scale to fit, e.g. log NU/log mu).
    labels: 1.0 if recovered, 0.0 if not.
    No numpy/sklearn in this environment; this is small (N ~ few hundred).
    """
    w0, w1, w2 = 0.0, 0.0, 0.0
    n = len(labels)
    mx = sum(xs) / n
    my = sum(ys) / n
    sx = (sum((x - mx) ** 2 for x in xs) / n) ** 0.5 or 1.0
    sy = (sum((y - my) ** 2 for y in ys) / n) ** 0.5 or 1.0
    xs_n = [(x - mx) / sx for x in xs]
    ys_n = [(y - my) / sy for y in ys]
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for x, y, t in zip(xs_n, ys_n, labels):
            z = w0 + w1 * x + w2 * y
            p = 1.0 / (1.0 + math.exp(-z)) if z > -50 else 0.0
            err = p - t
            g0 += err
            g1 += err * x
            g2 += err * y
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    # undo standardization: w0 + w1*(x-mx)/sx + w2*(y-my)/sy
    b1 = w1 / sx
    b2 = w2 / sy
    b0 = w0 - b1 * mx - b2 * my
    return b0, b1, b2


def predict_acc(b0, b1, b2, xs, ys, labels):
    correct = 0
    for x, y, t in zip(xs, ys, labels):
        z = b0 + b1 * x + b2 * y
        pred = 1.0 if z > 0 else 0.0
        correct += (pred == t)
    return correct / len(labels)


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — does mu separate recovery inside the ambiguous NU band?")
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
                step = (math.log2(r['prof'][m]) - math.log2(r['prof'][0])
                        if r['prof'][0] > 0 and r['prof'][m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff,
                          'lamstar': lam_star(lam, n), 'step': step})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'glv_hnp_phase2_thread25_rows.json'), 'w') as f:
        json.dump([{k: v for k, v in r.items() if k not in ('prof', 'nus')}
                   for r in rows], f)
    print("dumped rows -> glv_hnp_phase2_thread25_rows.json")

    print("\n" + "-" * 78)
    print(f"H25: stratify by NU band [{BAND_LO}, {BAND_HI}] "
          "(sufficient/necessary bracket from W4)")
    print("-" * 78)
    below = [r for r in rows if r['NU'] < BAND_LO]
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    above = [r for r in rows if r['NU'] > BAND_HI]
    print(f"below band (NU<{BAND_LO}): N={len(below)}  "
          f"rec={sum(1 for r in below if r['ok'])}/{len(below)}")
    print(f"IN band:                 N={len(band)}  "
          f"rec={sum(1 for r in band if r['ok'])}/{len(band)}")
    print(f"above band (NU>{BAND_HI}): N={len(above)}  "
          f"rec={sum(1 for r in above if r['ok'])}/{len(above)}")

    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    if pos and neg:
        auc_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        auc_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        auc_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        auc_step = auc([r['step'] for r in pos], [r['step'] for r in neg])
        auc_ls = auc([r['lamstar'] for r in pos], [r['lamstar'] for r in neg])
        print(f"\nWITHIN the band (N={len(band)}, {len(pos)} recovered / "
              f"{len(neg)} failed):")
        print(f"  AUC(-mu    -> recovery) = {auc_mu:.4f}   <- H25 target (need >= 0.8)")
        print(f"  AUC(-NU    -> recovery) = {auc_nu:.4f}   (should be ~0.5: NU is "
              "constant-ish inside the band by construction)")
        print(f"  AUC(-nu_hat-> recovery) = {auc_nh:.4f}")
        print(f"  AUC(-step  -> recovery) = {auc_step:.4f}   (W1b secondary)")
        print(f"  AUC(-lam*  -> recovery) = {auc_ls:.4f}   (control, "
              "Thread 20 falsified)")

        verdict = "CONFIRMED" if auc_mu >= 0.8 else "FALSIFIED"
        print(f"\nH25 verdict (raw mu): {verdict}  "
              f"(AUC(-mu) = {auc_mu:.4f} vs threshold 0.8)")
        verdict_nh = "CONFIRMED" if auc_nh >= 0.8 else "FALSIFIED"
        print(f"H25' verdict (nu_hat = mu/sqrt(det L2)): {verdict_nh}  "
              f"(AUC(-nu_hat) = {auc_nh:.4f} vs threshold 0.8)")
        print(f"note: AUC(-step) = {auc_step:.4f} means AUC(+step) = "
              f"{1-auc_step:.4f} -- step is a STRONG signal in the reversed "
              "direction (large step -> recovery), unexamined by H25/H25'.")
    else:
        print("\nband is degenerate (all-recover or all-fail); H25 untestable "
              "on this grid.")
        verdict = verdict_nh = "UNTESTABLE"

    print("\n" + "-" * 78)
    print("Cross-check: does mu separate OUTSIDE the band too (sanity)?")
    print("-" * 78)
    for name, sub in (("below", below), ("above", above)):
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            print(f"  {name}: AUC(-mu -> recovery) = "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):.4f}  "
                  f"(N={len(sub)})")
        else:
            print(f"  {name}: degenerate (N={len(sub)}, "
                  f"rec={sum(1 for r in sub if r['ok'])})")

    if verdict_nh == "CONFIRMED":
        print("\n" + "-" * 78)
        print("Logistic fit on (log NU, log nu_hat) -> recovery, band-restricted "
              "and pooled")
        print("-" * 78)
        for name, sub in (("in-band", band), ("pooled (all rows)", rows)):
            xs = [math.log(r['NU']) for r in sub]
            ys = [math.log(r['nuhat']) for r in sub]
            labels = [1.0 if r['ok'] else 0.0 for r in sub]
            if len(set(labels)) < 2:
                print(f"  {name}: degenerate labels, skip")
                continue
            b0, b1, b2 = logistic_fit_2d(xs, ys, labels)
            acc = predict_acc(b0, b1, b2, xs, ys, labels)
            print(f"  {name} (N={len(sub)}): "
                  f"logit(rec) = {b0:.3f} + {b1:.3f}*log(NU) + {b2:.3f}*log(nu_hat)"
                  f"   train-acc={acc:.3f}")

        print("\n" + "-" * 78)
        print("W1b follow-up: does step alone beat nu_hat in-band? "
              "(reversed sign: large step -> recovery)")
        print("-" * 78)
        xs = [-math.log(r['NU']) for r in band]
        ys = [r['step'] for r in band]  # NOT log; step already a log2-difference
        labels = [1.0 if r['ok'] else 0.0 for r in band]
        b0, b1, b2 = logistic_fit_2d(xs, ys, labels)
        acc = predict_acc(b0, b1, b2, xs, ys, labels)
        print(f"  in-band (N={len(band)}): "
              f"logit(rec) = {b0:.3f} + {b1:.3f}*(-log NU) + {b2:.3f}*step"
              f"   train-acc={acc:.3f}")
        print(f"  (compare nu_hat-based fit train-acc above)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
