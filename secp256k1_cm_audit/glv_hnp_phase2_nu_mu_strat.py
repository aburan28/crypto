"""
GLV-HNP Phase 2, Thread 25: does mu carry a second, NU-independent signal
inside the ambiguous NU band?

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

Background.  Thread 24's W5/W6 (glv_hnp_phase2_gsprofile_strat.py) held eff
FIXED and found nu_hat (equivalently mu = lambda_1(L2)) ranks curves with
AUC 0.75-0.93 per stratum, while the exact BDD certificate NU is markedly
worse (0.35-0.73) and its rank correlation with the nu_hat*sqrt(eff) closed
form is near zero or NEGATIVE inside every stratum (W6). So mu and NU look
like two different mechanisms. The direct test: pick out exactly the
instances where NU is uninformative (the ambiguous band) and ask whether mu
still separates there. If yes, (NU, mu) is a genuine 2-parameter viability
test. If mu's power turns out to be entirely mediated by NU / eff after all,
this is a stratification artifact and the closed form should be retired.

Secondary (from Thread 24's W1b): does the GS-profile step at the K1 wall,
    step = log2(||b*_{m+1}||) - log2(||b*_1||)      (0 -> tail-block onset)
predict the wall better than NU or mu, inside the same ambiguous band?

Run: python3 glv_hnp_phase2_nu_mu_strat.py
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman
from glv_hnp_phase2_gsprofile_strat import build_rows, EFFS, M17

NU_LO, NU_HI = 1.040, 2.199  # Thread 24 W4, 17-bit bracket (300 inst, dim 24)


def logistic_fit_2d(xs, ys, labels, iters=4000, lr=0.5):
    """Plain-gradient-descent logistic regression on 2 standardised features.
    Returns (w0, w1, w2, mx, sx, my, sy) so the boundary can be evaluated on
    raw (x, y) via  w0 + w1*(x-mx)/sx + w2*(y-my)/sy = 0.
    """
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs) / n) or 1.0
    sy = math.sqrt(sum((y - my) ** 2 for y in ys) / n) or 1.0
    zx = [(x - mx) / sx for x in xs]
    zy = [(y - my) / sy for y in ys]
    w0 = w1 = w2 = 0.0
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for xi, yi, li in zip(zx, zy, labels):
            z = w0 + w1 * xi + w2 * yi
            p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
            err = p - li
            g0 += err
            g1 += err * xi
            g2 += err * yi
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    return w0, w1, w2, mx, sx, my, sy


def logistic_predict(w, x, y):
    w0, w1, w2, mx, sx, my, sy = w
    z = w0 + w1 * (x - mx) / sx + w2 * (y - my) / sy
    return 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — mu inside the ambiguous NU band (H25)")
    print("=" * 78)

    import time
    t0 = time.time()
    rows, curves17 = build_rows()
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves, {len(rows)} instances "
          f"(float GS, dim {rows[0]['k']}) in {time.time()-t0:.1f}s")
    print(f"ambiguous band (Thread 24 W4 bracket): {NU_LO} <= NU <= {NU_HI}")

    band = [r for r in rows if NU_LO <= r['NU'] <= NU_HI]
    below = [r for r in rows if r['NU'] < NU_LO]
    above = [r for r in rows if r['NU'] > NU_HI]
    print(f"\nsplit: below {len(below)} (rec "
          f"{sum(1 for r in below if r['ok'])}/{len(below)}), "
          f"band {len(band)} (rec {sum(1 for r in band if r['ok'])}/{len(band)}), "
          f"above {len(above)} (rec {sum(1 for r in above if r['ok'])}/{len(above)})")

    print("\n" + "-" * 78)
    print("EXP X1: AUC(-mu) and AUC(-nu_hat) INSIDE the ambiguous band")
    print("-" * 78)
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_ef = auc([r['eff'] for r in pos], [r['eff'] for r in neg])
        print(f"N={len(band)}  pos={len(pos)}  neg={len(neg)}")
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}   (should be ~0.5: NU is")
        print(f"                                           uninformative BY")
        print(f"                                           CONSTRUCTION in-band)")
        print(f"  AUC(-eff    -> recovery) = {a_ef:.4f}")
        print(f"\nH25 (AUC(-mu) >= 0.8 in-band): "
              f"{'HOLDS' if a_mu >= 0.8 else 'FALSIFIED'} ({a_mu:.4f})")
    else:
        print("degenerate band (all-pos or all-neg) -- cannot compute AUC")

    print("\nper-eff-stratum, restricted to the ambiguous band:")
    print(f"{'eff':>5} {'N':>4} {'rec':>7} | {'AUC mu':>8} {'AUC nu_hat':>11}")
    for eff in EFFS:
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if not sub:
            continue
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f} "
                  f"{auc([r['nuhat'] for r in p], [r['nuhat'] for r in ng]):>11.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | (degenerate)")

    print("\n" + "-" * 78)
    print("EXP X2: is mu's in-band power mediated by eff, or independent of it?")
    print("Spearman(mu, eff) and partial check: AUC(-mu) computed only WITHIN")
    print("each eff stratum (X1 table above) vs pooled across strata (H25).")
    print("-" * 78)
    sp = spearman([r['mu'] for r in band], [r['eff'] for r in band])
    print(f"Spearman(mu, eff) over the band, N={len(band)} = {sp:.4f}")
    print("If the per-stratum AUCs above are all high AND close to the pooled")
    print("H25 figure, mu separates WITHIN fixed eff too -- not just via eff.")

    print("\n" + "-" * 78)
    print("EXP X3: 2-feature logistic fit on (log NU, log mu), in-band only")
    print("-" * 78)
    if pos and neg and len(band) >= 10:
        xs = [math.log(r['NU']) for r in band]
        ys = [math.log(r['mu']) for r in band]
        labels = [1.0 if r['ok'] else 0.0 for r in band]
        w = logistic_fit_2d(xs, ys, labels)
        preds = [logistic_predict(w, x, y) for x, y in zip(xs, ys)]
        auc_fit = auc([p for p, l in zip(preds, labels) if l == 1.0],
                       [p for p, l in zip(preds, labels) if l == 0.0])
        # auc() here scores "-x -> recovery"; predicted prob is already
        # "recovery likely", so flip: AUC(+preds -> recovery) = 1 - auc(-preds)
        print(f"fit: logit(p) = {w[0]:.3f} + {w[1]:.3f}*z(log NU) + "
              f"{w[2]:.3f}*z(log mu)")
        print(f"     [z(x) = (x - {w[3]:.3f}) / {w[4]:.3f} for log NU, "
              f"(x - {w[5]:.3f}) / {w[6]:.3f} for log mu]")
        print(f"in-sample AUC of the fitted 2-feature model (in-band) = "
              f"{1.0 - auc_fit:.4f}")
        print("(in-sample only -- this is a descriptive fit of the band, not")
        print(" a held-out test; compare to the single-feature AUC(-mu) above)")
    else:
        print("band too small / degenerate for a fit")

    print("\n" + "-" * 78)
    print("EXP X4 (secondary): does the GS-profile step at the K1 wall predict")
    print("        recovery inside the band better than NU or mu?")
    print("        step = log2(||b*_{m+1}||) - log2(||b*_1||)")
    print("-" * 78)
    m = M17
    for r in rows:
        prof = r.get('prof')
        r['step'] = (math.log2(prof[m]) - math.log2(prof[0])
                     if prof and prof[0] > 0 and len(prof) > m and prof[m] > 0
                     else None)
    have_step = [r for r in rows if r['step'] is not None]
    print(f"rows with a usable step value: {len(have_step)}/{len(rows)}")
    band_step = [r for r in band if r['step'] is not None]
    posS = [r for r in band_step if r['ok']]
    negS = [r for r in band_step if not r['ok']]
    if posS and negS:
        a_step = auc([r['step'] for r in posS], [r['step'] for r in negS])
        a_mu_s = auc([r['mu'] for r in posS], [r['mu'] for r in negS])
        a_nu_s = auc([r['NU'] for r in posS], [r['NU'] for r in negS])
        print(f"in-band (N={len(band_step)}): AUC(-step) = {a_step:.4f}   "
              f"AUC(-mu) = {a_mu_s:.4f}   AUC(-NU) = {a_nu_s:.4f}")
        print(f"step | success: mean {sum(r['step'] for r in posS)/len(posS):.3f}   "
              f"step | failure: mean {sum(r['step'] for r in negS)/len(negS):.3f}")
    else:
        print("degenerate / no usable step values in-band")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
