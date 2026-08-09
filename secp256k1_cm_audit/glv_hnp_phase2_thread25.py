"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  W5/W6 showed recovery = f(NU, X), X ~ mu-driven and independent of NU at
  fixed eff.  H25: within the ambiguous NU band [1.04, 2.20] (where the
  nearest-plane certificate gives no answer), AUC(-mu -> Kannan-LLL
  recovery) stays >= 0.8.
  Falsifier: if mu's apparent power is entirely mediated by NU, the band AUC
  collapses toward 0.5 and the W5 result is a stratification artifact.

Secondary: step = log2(||b*_{m+1}||) - log2(||b*_1||) (the GS profile's
head-to-tail-block jump, W1b).  Test whether step predicts the wall better
than NU or mu.

Uses the 500-row dump from glv_hnp_phase2_gsprofile_strat.py
(--dump-json), so this is pure re-analysis: no new lattice computation.

Run: python3 glv_hnp_phase2_thread25.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from glv_hnp_phase2_gsprofile import auc, spearman

DUMP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "glv_hnp_phase2_gsprofile_strat_dump.json")


def band_auc(rows, lo, hi, key, invert=False):
    sub = [r for r in rows if lo <= r["NU"] <= hi]
    pos = [r for r in sub if r["ok"]]
    neg = [r for r in sub if not r["ok"]]
    if not pos or not neg:
        return None, len(sub), len(pos)
    xs_pos = [r[key] for r in pos]
    xs_neg = [r[key] for r in neg]
    a = auc(xs_pos, xs_neg)
    if invert:
        a = 1.0 - a
    return a, len(sub), len(pos)


def logistic_fit(xs, ys, lr=0.05, iters=20000):
    """2-feature logistic regression by gradient descent, no deps."""
    n = len(xs)
    d = len(xs[0])
    w = [0.0] * d
    b = 0.0
    for _ in range(iters):
        gw = [0.0] * d
        gb = 0.0
        for x, y in zip(xs, ys):
            z = sum(wi * xi for wi, xi in zip(w, x)) + b
            p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
            err = p - y
            for j in range(d):
                gw[j] += err * x[j]
            gb += err
        w = [wi - lr * gwi / n for wi, gwi in zip(w, gw)]
        b -= lr * gb / n
    return w, b


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a genuine second coordinate?")
    print("=" * 78)

    rows = json.load(open(DUMP))
    print(f"\n{len(rows)} rows loaded from {os.path.basename(DUMP)}")

    print("\n" + "-" * 78)
    print("H25: AUC(-mu -> recovery) inside the ambiguous NU band [1.04, 2.20]")
    print("-" * 78)
    LO, HI = 1.04, 2.20
    a_mu, n_band, n_pos = band_auc(rows, LO, HI, "mu")
    a_nuhat, _, _ = band_auc(rows, LO, HI, "nuhat")
    a_eff, _, _ = band_auc(rows, LO, HI, "eff")
    a_step, _, _ = band_auc(rows, LO, HI, "step", invert=True)
    print(f"band size N={n_band}  recoveries={n_pos}")
    print(f"  AUC(-mu)        = {a_mu:.4f}" if a_mu is not None else "  (degenerate)")
    print(f"  AUC(-nu_hat)    = {a_nuhat:.4f}" if a_nuhat is not None else "")
    print(f"  AUC(-eff)       = {a_eff:.4f}" if a_eff is not None else "")
    print(f"  AUC(+step)      = {a_step:.4f}" if a_step is not None else "")

    verdict = "CONFIRMED" if (a_mu is not None and a_mu >= 0.8) else "FALSIFIED"
    print(f"\nH25 (AUC(-mu) >= 0.8 inside the band): {verdict} "
          f"(observed {a_mu:.4f})" if a_mu is not None else "\nH25: degenerate band")

    print("\n" + "-" * 78)
    print("Cross-check: same test restricted to each eff stratum separately,")
    print("since the pooled band still spans 5 bias strengths (n varies little,")
    print("K1 varies with eff, so this checks the band result isn't itself")
    print("eff-driven).")
    print("-" * 78)
    EFFS = sorted(set(r["effq"] for r in rows))
    print(f"{'eff':>5} {'N_band':>7} {'rec':>5} | {'AUC -mu':>8} {'AUC -nu_hat':>11}")
    for eff in EFFS:
        sub = [r for r in rows if r["effq"] == eff]
        a, nb, npos = band_auc(sub, LO, HI, "mu")
        a2, _, _ = band_auc(sub, LO, HI, "nuhat")
        if a is None:
            print(f"{eff:>5.2f} {nb:>7} {npos:>5} | {'(degenerate)':>8}")
        else:
            print(f"{eff:>5.2f} {nb:>7} {npos:>5} | {a:>8.4f} {a2:>11.4f}")

    print("\n" + "-" * 78)
    print("Secondary: does step = log2||b*_{m+1}|| - log2||b*_1|| beat NU/mu")
    print("pooled over ALL rows (not just the band)?")
    print("-" * 78)
    pos_all = [r for r in rows if r["ok"]]
    neg_all = [r for r in rows if not r["ok"]]
    print(f"  AUC(+step -> recovery) = "
          f"{1.0 - auc([r['step'] for r in pos_all], [r['step'] for r in neg_all]):.4f}"
          "   (larger step = more separated blocks = easier)")
    print(f"  AUC(-NU   -> recovery) = "
          f"{auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all]):.4f}")
    print(f"  AUC(-mu   -> recovery) = "
          f"{auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all]):.4f}")
    print(f"  Spearman(step, NU)     = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu)     = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    if a_mu is not None and a_mu >= 0.8:
        print("\n" + "-" * 78)
        print("H25 confirmed -> fit logistic decision boundary on (log NU, log mu)")
        print("-" * 78)
        eps = 1e-9
        xs = [[math.log(r["NU"] + eps), math.log(r["mu"] + eps)] for r in rows]
        mx0 = sum(x[0] for x in xs) / len(xs)
        mx1 = sum(x[1] for x in xs) / len(xs)
        sx0 = (sum((x[0] - mx0) ** 2 for x in xs) / len(xs)) ** 0.5
        sx1 = (sum((x[1] - mx1) ** 2 for x in xs) / len(xs)) ** 0.5
        xs_n = [[(x[0] - mx0) / sx0, (x[1] - mx1) / sx1] for x in xs]
        ys = [1.0 if r["ok"] else 0.0 for r in rows]
        w, b = logistic_fit(xs_n, ys)
        # unnormalise: w0*(lnNU-mx0)/sx0 + w1*(lnmu-mx1)/sx1 + b = 0
        a0 = w[0] / sx0
        a1 = w[1] / sx1
        c = b - w[0] * mx0 / sx0 - w[1] * mx1 / sx1
        print(f"decision boundary: {a0:.4f}*ln(NU) + {a1:.4f}*ln(mu) + {c:.4f} = 0")
        preds = []
        for x, y in zip(xs, ys):
            z = a0 * x[0] + a1 * x[1] + c
            p = 1.0 / (1.0 + math.exp(-z)) if abs(z) < 30 else (1.0 if z > 0 else 0.0)
            preds.append(1.0 if p > 0.5 else 0.0)
        acc = sum(1 for p, y in zip(preds, ys) if p == y) / len(ys)
        print(f"train accuracy (same 500 rows, no held-out split): {acc:.4f}")
    else:
        print("\nH25 falsified: skipping logistic-boundary fit "
              "(mu's power inside the band did not clear 0.8; likely mediated "
              "by NU / eff after all).")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
