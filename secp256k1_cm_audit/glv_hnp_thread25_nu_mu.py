"""
GLV-HNP Phase 2, Thread 25: does mu = lambda_1(L2) carry information beyond
NU, or is its W5 power entirely mediated by NU?

Pre-registered by the 2026-08-07 (Thread 24) log entry, W5/W6:
  - Holding eff fixed, NU (AUC 0.35-0.73) and nu_hat/mu (AUC 0.75-0.93) are
    both cross-curve predictors but are RANK-UNCORRELATED (Spearman ~ -0.2
    to +0.16 against the closed form NU ~ C*nu_hat*sqrt(eff)).
  - W4 gives an NU bracket at 17 bits: sufficient NU < 1.040 (zero FP),
    necessary NU > 2.199 (zero recoveries above). Inside [1.040, 2.199],
    nearest-plane gives no answer either way.

  H25: within the ambiguous NU band, AUC(-mu -> recovery) stays >= 0.8.
  Falsifier: AUC(-mu -> recovery) inside the band collapses to ~0.5, i.e.
    mu's apparent power in W5 was entirely an artifact of also being
    correlated with NU (stratification leakage), not a second coordinate.

Secondary (also proposed by Thread 24's W1b): step = log2(prof[m]) -
log2(prof[0]), i.e. how far the GS profile jumps from the (flat) head block
of m copies of lambda_1(L2) to the first index of the second block. Tests
whether step -> 0 predicts the wall independently of NU/mu.

Data: reads the 500-row dump from glv_hnp_phase2_gsprofile_strat.py
(17 bits, 20 curves x 5 eff strata x 5 seeds, dim 24, float GS -- W0 of the
parent script bounds the float/exact NU gap at ~1e-15 relative at this
dimension, so float GS is not a confound here).

Run:
  python3 glv_hnp_phase2_gsprofile_strat.py --dump-json
  python3 glv_hnp_thread25_nu_mu.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

DUMP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "glv_hnp_phase2_gsprofile_strat_dump.json")

# W4 (17-bit) bracket: sufficient NU < 1.040, necessary NU > 2.199.
BAND_LO, BAND_HI = 1.040, 2.199


def logistic_fit_2d(xs, ys, labels, iters=4000, lr=0.05):
    """Minimal logistic regression on (log NU, log mu) -> recovery, no deps."""
    w0, w1, w2 = 0.0, 0.0, 0.0
    n = len(xs)
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for x, y, t in zip(xs, ys, labels):
            z = w0 + w1 * x + w2 * y
            p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
            err = p - t
            g0 += err
            g1 += err * x
            g2 += err * y
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    return w0, w1, w2


if __name__ == "__main__":
    if not os.path.exists(DUMP):
        print(f"no dump at {DUMP}; run:\n"
              f"  python3 glv_hnp_phase2_gsprofile_strat.py --dump-json")
        sys.exit(1)

    with open(DUMP) as f:
        rows = json.load(f)
    print("=" * 78)
    print(f"Thread 25 — mu vs NU inside the ambiguous band  (N={len(rows)} rows)")
    print("=" * 78)

    for r in rows:
        m = r['k'] // 2
        r['step'] = math.log2(r['prof'][m]) - math.log2(r['prof'][0])

    print("\n" + "-" * 78)
    print(f"H25: AUC(-mu -> recovery) inside NU in [{BAND_LO}, {BAND_HI}]")
    print("-" * 78)
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band N = {len(band)} ({len(pos)} recover, {len(neg)} fail); "
          f"outside band N = {len(rows) - len(band)}")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([r['step'] for r in pos], [r['step'] for r in neg])
        print(f"  AUC(-mu     -> recovery) inside band = {a_mu:.4f}"
              f"  {'CONFIRMS H25 (>=0.8)' if a_mu >= 0.8 else 'FALSIFIES H25 (<0.8)'}")
        print(f"  AUC(-nu_hat -> recovery) inside band = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) inside band = {a_nu:.4f}  "
              f"(expected ~0.5: NU is uninformative by construction of the band)")
        print(f"  AUC(-step   -> recovery) inside band = {a_st:.4f}")
    else:
        print("  degenerate: one class empty inside the band")

    print("\nper-eff breakdown inside the band:")
    print(f"{'eff':>5} {'N':>4} {'rec':>6} {'AUC mu':>8} {'AUC NU':>8}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        ng = [r for r in sub if not r['ok']]
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>4} {str(len(p))+'/'+str(len(sub)):>6} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f} "
                  f"{auc([r['NU'] for r in p], [r['NU'] for r in ng]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>6} (degenerate)")

    print("\n" + "-" * 78)
    print("Secondary: does step = log2(||b*_{m+1}||) - log2(||b*_1||) predict "
          "the wall?")
    print("-" * 78)
    allpos = [r for r in rows if r['ok']]
    allneg = [r for r in rows if not r['ok']]
    print(f"pooled (N={len(rows)}):")
    print(f"  AUC(step -> recovery), larger step = more likely: "
          f"{1 - auc([r['step'] for r in allpos], [r['step'] for r in allneg]):.4f}")
    print(f"  AUC(-NU  -> recovery)                             : "
          f"{auc([r['NU'] for r in allpos], [r['NU'] for r in allneg]):.4f}")
    print(f"  AUC(-mu  -> recovery)                             : "
          f"{auc([r['mu'] for r in allpos], [r['mu'] for r in allneg]):.4f}")
    print(f"  Spearman(step, NU)  = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu)  = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "-" * 78)
    print("Decision boundary: logistic fit on (log NU, log mu) -> recovery")
    print("-" * 78)
    xs = [math.log(r['NU']) for r in rows]
    ys = [math.log(r['mu']) for r in rows]
    ts = [1.0 if r['ok'] else 0.0 for r in rows]
    w0, w1, w2 = logistic_fit_2d(xs, ys, ts)
    preds = []
    for x, y in zip(xs, ys):
        z = w0 + w1 * x + w2 * y
        preds.append(1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0)
    correct = sum(1 for p, t in zip(preds, ts) if (p >= 0.5) == (t >= 0.5))
    print(f"logit(recover) = {w0:.4f} + {w1:.4f}*log(NU) + {w2:.4f}*log(mu)")
    print(f"training accuracy at 0.5 threshold: {correct}/{len(ts)} "
          f"({100*correct/len(ts):.1f}%)")
    a_nu_only = auc([r['NU'] for r in rows if r['ok']],
                     [r['NU'] for r in rows if not r['ok']])
    print(f"(compare: NU-alone AUC = {a_nu_only:.4f}; "
          f"nu_hat*sqrt(eff)-alone AUC = 0.9348 pooled per Thread 24 W5)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
