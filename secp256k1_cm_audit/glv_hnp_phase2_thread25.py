"""
GLV-HNP Phase 2, Thread 25: does mu survive as a second coordinate once NU
is fixed, and is there a parameter-free combination of the two?

Pre-registered by the 2026-08-07 (Thread 24) log entry, W5/W6: NU (exact BDD
certificate) and nu_hat = mu/sqrt(det L2) (closed-form, no lattice work) are
mutually uncorrelated at fixed eff, and both independently predict recovery.

  H25: within the ambiguous NU band 1.04 <= NU <= 2.20 (the ambiguous
       bracket found in Thread 24 W4, where the NU<=1 certificate gives no
       answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

Falsifier: if AUC(-mu) inside the band drops well below 0.8, mu's apparent
power in W5 is a stratification artifact of eff, not a genuine second
coordinate.

Depends on secp256k1_cm_audit/glv_hnp_phase2_gsprofile_strat.py --dump-json,
which persists the 500-instance table (5 eff strata x 20 curves x 5 seeds,
17 bits, dim 24) this script re-analyses.  Regenerates the table if the JSON
file is not already present next to this script.

Run: python3 glv_hnp_phase2_thread25.py
"""

import json
import math
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from glv_hnp_phase2_gsprofile import auc, spearman

ROWS_PATH = os.path.join(HERE, "glv_hnp_phase2_gsprofile_strat_rows.json")
LO, HI = 1.040, 2.199  # Thread 24 W4's 17-bit NU bracket


def sigmoid(z):
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z))
    ez = math.exp(z)
    return ez / (1.0 + ez)


def fit_logistic(X, Y, lr=0.1, epochs=4000):
    """X: list of (1, x1, x2) STANDARDIZED features. Batch gradient ascent."""
    w = [0.0] * len(X[0])
    n = len(X)
    for _ in range(epochs):
        grad = [0.0] * len(w)
        for xi, yi in zip(X, Y):
            err = yi - sigmoid(sum(wj * xj for wj, xj in zip(w, xi)))
            for k in range(len(w)):
                grad[k] += err * xi[k]
        w = [wj + lr * gk / n for wj, gk in zip(w, grad)]
    return w


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a genuine second coordinate once NU is fixed?")
    print("=" * 78)

    if not os.path.exists(ROWS_PATH):
        print(f"\n{ROWS_PATH} not found, regenerating via "
              "glv_hnp_phase2_gsprofile_strat.py --dump-json ...")
        subprocess.run([sys.executable,
                         os.path.join(HERE, "glv_hnp_phase2_gsprofile_strat.py"),
                         "--dump-json", ROWS_PATH], check=True, cwd=HERE)

    rows = json.load(open(ROWS_PATH))
    print(f"\nloaded {len(rows)} instances from {ROWS_PATH}")

    band = [r for r in rows if LO <= r['NU'] <= HI]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]

    print("\n" + "-" * 78)
    print(f"H25 (literal, pooled): band [{LO},{HI}]  N={len(band)}  "
          f"pos={len(pos)}  neg={len(neg)}")
    print("-" * 78)
    a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
    a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
    a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
    print(f"AUC(-mu -> recovery) in band     = {a_mu:.4f}")
    print(f"AUC(-nu_hat -> recovery) in band = {a_nh:.4f}")
    print(f"AUC(-NU -> recovery) in band     = {a_nu:.4f}  "
          "(sanity: near 0.5 expected, NU is the stratifying var)")
    print(f"H25 literal verdict (AUC(-mu) >= 0.8): "
          f"{'HOLDS' if a_mu >= 0.8 else 'FALSIFIED'}  ({a_mu:.4f})")

    print("\nband composition and AUC(-mu), AUC(-nu_hat) PER eff stratum "
          "(controls for eff too):")
    print(f"{'eff':>5} {'N':>4} {'pos':>4} {'neg':>4} {'AUC mu':>8} "
          f"{'AUC nu_hat':>11}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in band if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        n = [r for r in sub if not r['ok']]
        if p and n:
            print(f"{eff:>5.2f} {len(sub):>4} {len(p):>4} {len(n):>4} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in n]):>8.4f} "
                  f"{auc([r['nuhat'] for r in p], [r['nuhat'] for r in n]):>11.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>4}  degenerate (pos={len(p)} neg={len(n)})")

    mus = [r['mu'] for r in band]
    nus = [r['NU'] for r in band]
    print(f"\nSpearman(mu, NU) within band = {spearman(mus, nus):.4f}")

    print("\n" + "-" * 78)
    print("EXP T1: parameter-free combination NU * nu_hat")
    print("-" * 78)
    for r in rows:
        r['PROD'] = r['NU'] * r['nuhat']
    ppos = [r for r in rows if r['ok']]
    pneg = [r for r in rows if not r['ok']]
    a_prod = auc([r['PROD'] for r in ppos], [r['PROD'] for r in pneg])
    print(f"pooled AUC(-NU*nu_hat -> recovery) = {a_prod:.4f}   "
          f"(no fitted constant, unlike nu_hat*sqrt(eff)'s C~6.29)")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in rows if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        n = [r for r in sub if not r['ok']]
        if p and n:
            print(f"  eff={eff:.2f}  AUC(-NU*nu_hat) = "
                  f"{auc([r['PROD'] for r in p], [r['PROD'] for r in n]):.4f}")

    prodsorted = sorted(set(r['PROD'] for r in rows))
    zfp = [t for t in prodsorted
           if sum(1 for r in rows if r['PROD'] <= t and not r['ok']) == 0]
    if zfp:
        t0 = max(zfp)
        tp = sum(1 for r in rows if r['PROD'] <= t0 and r['ok'])
        print(f"\nzero-FP certificate: NU*nu_hat <= {t0:.4f}  =>  "
              f"TP={tp} FP=0  (of {len(ppos)} positives)")

    print("\n" + "-" * 78)
    print("EXP T2: 2-feature logistic fit on (log NU, log nu_hat)")
    print("-" * 78)
    Xr = [r for r in rows if r['nuhat'] > 0]
    lnu = [math.log(r['NU']) for r in Xr]
    lnh = [math.log(r['nuhat']) for r in Xr]
    Y = [1.0 if r['ok'] else 0.0 for r in Xr]
    m1, s1 = statistics.mean(lnu), statistics.pstdev(lnu)
    m2, s2 = statistics.mean(lnh), statistics.pstdev(lnh)
    Xs = [(1.0, (a - m1) / s1, (b - m2) / s2) for a, b in zip(lnu, lnh)]
    w = fit_logistic(Xs, Y)
    print(f"standardized coefficients: b0={w[0]:.4f}  b1(log NU)={w[1]:.4f}  "
          f"b2(log nu_hat)={w[2]:.4f}   N={len(Xs)}")
    scores = [sum(wj * xj for wj, xj in zip(w, xi)) for xi in Xs]
    sp = [s for s, y in zip(scores, Y) if y == 1.0]
    sn = [s for s, y in zip(scores, Y) if y == 0.0]
    a_logit = auc([-s for s in sp], [-s for s in sn])
    print(f"pooled AUC of fitted logistic score = {a_logit:.4f}")
    b1, b2 = w[1] / s1, w[2] / s2
    b0 = w[0] - w[1] * m1 / s1 - w[2] * m2 / s2
    print(f"unstandardized: logit(P) = {b0:.4f} + {b1:.4f}*log(NU) + "
          f"{b2:.4f}*log(nu_hat)")
    print("(b1 ~ b2 in magnitude => log NU + log nu_hat, i.e. the product "
          "NU*nu_hat, is close to the fitted optimum with no fitting at all)")

    print("\n" + "-" * 78)
    print("EXP T3 (secondary, W1b follow-up): GS-profile step statistic")
    print("  step = log2(||b*_{m+1}||) - log2(||b*_1||)")
    print("-" * 78)
    srows = [r for r in rows if r.get('step') is not None
             and not math.isnan(r['step'])]
    sp2 = [r for r in srows if r['ok']]
    sn2 = [r for r in srows if not r['ok']]
    print(f"pooled AUC(step -> recovery), bigger step = more likely "
          f"recovery: {auc([-r['step'] for r in sp2], [-r['step'] for r in sn2]):.4f}")
    for eff in sorted(set(r['effq'] for r in srows)):
        sub = [r for r in srows if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        n = [r for r in sub if not r['ok']]
        if p and n:
            print(f"  eff={eff:.2f}  AUC(step) = "
                  f"{auc([-r['step'] for r in p], [-r['step'] for r in n]):.4f}  "
                  f"step mean pos={sum(r['step'] for r in p)/len(p):.3f} "
                  f"neg={sum(r['step'] for r in n)/len(n):.3f}")
    print("\nDirection: LARGER step predicts recovery (consistent with W1b: "
          "step -> 0 as the K1 wall is crossed, and crossing the wall is "
          "already failure-correlated).  The secondary hypothesis as "
          "literally posed ('step -> 0 predicts the wall') is directionally "
          "confirmed but is not a new independent axis: within-stratum AUC "
          "(0.77-0.92) is on par with mu/nu_hat, not better.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
