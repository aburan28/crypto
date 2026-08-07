"""
GLV-HNP Phase 2, Thread 25: is mu a second, independent coordinate of the
LLL wall, or is its apparent power mediated by NU after all?

Thread 24 (W5/W6) found that the exact BDD certificate NU and the closed-form
nu_hat = mu/sqrt(det L2) are mutually UNCORRELATED at fixed bias strength
(Spearman -0.28..+0.16 inside every eff stratum) yet both predict Kannan-LLL
recovery.  Two uncorrelated predictors of the same event means two mechanisms.
The pre-registered test:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 -- where the nearest-plane
       guarantee gives NO answer -- AUC(-mu -> recovery) stays >= 0.8.

  Falsifier: if AUC(-mu | NU band) collapses toward 0.5, mu's power in W5 was
       a stratification artifact and the closed form should be retired.

X0  build + dump the table (JSON) so it survives the run
X1  H25 proper: AUC(-mu) inside NU bands, single- and double-conditioned
X2  two-parameter logistic fit on (log NU, log mu); decision boundary
X3  W1b tertiary: the GS two-block step
        step = log2||b*_{m+1}|| - log2||b*_1||
    W1b showed the step vanishes as the K1 wall is crossed.  Does step
    predict the wall better than NU or mu?

Gram-Schmidt is float, justified by W0/W4 of glv_hnp_phase2_gsprofile.py
(max relative NU error vs exact Fractions ~1e-15 at dim 24).

Run: python3 glv_hnp_phase2_nu_conditioned.py [--dump-json FILE]
"""

import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# Thread 24 used 5 seeds; conditioning on NU splits the table further, so we
# need more instances per (curve, eff) cell to keep the bands populated.
SEEDS25 = [42, 1234, 9999, 555, 31337, 7, 20260807, 314159, 2718, 161803,
           99991, 123123, 4444, 8675309, 271828]
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
M17 = 12

# Thread 24 W4, 17 bits: sufficient NU < 1.040, necessary NU > 2.199.
BAND_LO, BAND_HI = 1.040, 2.199


def auc_p(pos, neg):
    """(AUC, normal-approximation two-sided p) for the score -x.

    Mann-Whitney U under H0: AUC = 0.5, sd = sqrt((n1+n2+1)/(12 n1 n2)).
    Ties are handled by the 0.5 credit in auc(); the sd is the untied one, so
    p is mildly conservative when ties are common.
    """
    a = auc(pos, neg)
    n1, n2 = len(pos), len(neg)
    if not n1 or not n2:
        return a, float('nan')
    sd = math.sqrt((n1 + n2 + 1.0) / (12.0 * n1 * n2))
    z = (a - 0.5) / sd if sd else 0.0
    return a, math.erfc(abs(z) / math.sqrt(2.0))


def logistic_fit(X, y, ridge=1e-6, iters=100):
    """Newton-IRLS with a ridge term.  X: rows WITHOUT intercept.

    Returns (beta, loglik) with beta[0] the intercept.
    """
    n = len(X)
    d = len(X[0]) + 1
    A = [[1.0] + list(row) for row in X]
    beta = [0.0] * d
    for _ in range(iters):
        eta = [sum(b * v for b, v in zip(beta, row)) for row in A]
        pr = [1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, e)))) for e in eta]
        g = [sum(A[i][j] * (y[i] - pr[i]) for i in range(n)) - ridge * beta[j]
             for j in range(d)]
        H = [[sum(A[i][j] * A[i][k] * pr[i] * (1 - pr[i]) for i in range(n))
              + (ridge if j == k else 0.0) for k in range(d)] for j in range(d)]
        step = _solve(H, g)
        if step is None:
            break
        beta = [b + s for b, s in zip(beta, step)]
        if max(abs(s) for s in step) < 1e-10:
            break
    eta = [sum(b * v for b, v in zip(beta, row)) for row in A]
    pr = [1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, e)))) for e in eta]
    ll = sum(y[i] * math.log(max(pr[i], 1e-15))
             + (1 - y[i]) * math.log(max(1 - pr[i], 1e-15)) for i in range(n))
    return beta, ll


def _solve(M, v):
    """Gaussian elimination with partial pivoting."""
    d = len(v)
    A = [row[:] + [v[i]] for i, row in enumerate(M)]
    for c in range(d):
        piv = max(range(c, d), key=lambda r: abs(A[r][c]))
        if abs(A[piv][c]) < 1e-14:
            return None
        A[c], A[piv] = A[piv], A[c]
        for r in range(d):
            if r == c:
                continue
            f = A[r][c] / A[c][c]
            for k in range(c, d + 1):
                A[r][k] -= f * A[c][k]
    return [A[i][d] / A[i][i] for i in range(d)]


def build_table():
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    t0 = time.time()
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in SEEDS25:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M17, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                prof = r['prof']
                # W1b two-block step: head is m copies of lambda_1(L2); the
                # second block starts at index m and collapses onto the head
                # exactly as the wall is crossed.
                step = (math.log2(prof[M17] / prof[0])
                        if prof[0] > 0 and prof[M17] > 0 else float('nan'))
                rows.append({
                    'n': n, 'K1': k1b, 'seed': seed, 'ok': bool(rk['ok']),
                    'eff': k1b * k2b / n, 'effq': eff,
                    'NU': r['NU'], 'mu': r['mu'], 'nuhat': r['nuhat'],
                    'l2': r['l2'], 'det2': r['det2'], 'argmax': r['argmax'],
                    'enorm': r['enorm'], 'step': step,
                    'lamstar': lam_star(lam, n), 'k': r['k'],
                })
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    return rows


def band_report(rows, score_keys, label, sign=None):
    """AUC table over NU bands for each score in score_keys."""
    sign = sign or {}
    bands = [(0.0, BAND_LO, "sufficient"),
             (BAND_LO, 1.5, "ambiguous-lo"),
             (1.5, BAND_HI, "ambiguous-hi"),
             (BAND_LO, BAND_HI, "AMBIGUOUS (H25)"),
             (BAND_HI, 1e9, "necessary")]
    print(f"\n{label}")
    hdr = f"{'NU band':>16} {'N':>6} {'rec':>9} |"
    for k in score_keys:
        hdr += f" {'AUC ' + k:>12} {'p':>9} |"
    print(hdr)
    for lo, hi, name in bands:
        sub = [r for r in rows if lo <= r['NU'] < hi
               and not math.isnan(r.get('step', 0.0))]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        line = (f"{name:>16} {len(sub):>6} "
                f"{str(len(pos)) + '/' + str(len(sub)):>9} |")
        if not pos or not neg:
            print(line + "   (degenerate — one class empty)")
            continue
        for k in score_keys:
            s = sign.get(k, 1.0)
            a, pv = auc_p([s * r[k] for r in pos], [s * r[k] for r in neg])
            line += f" {a:>12.4f} {pv:>9.2e} |"
        print(line)


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a second coordinate?  Conditioning on NU (H25)")
    print("=" * 78)

    dump = None
    if "--dump-json" in sys.argv:
        dump = sys.argv[sys.argv.index("--dump-json") + 1]

    rows = build_table()
    if dump:
        with open(dump, "w") as f:
            json.dump(rows, f)
        print(f"table dumped to {dump} ({os.path.getsize(dump)} bytes)")

    ok = sum(1 for r in rows if r['ok'])
    print(f"overall recovery {ok}/{len(rows)} = {ok/len(rows):.3f}")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1a: H25 proper — AUC inside NU bands.")
    print("  AUC > 0.5 means SMALLER score -> more likely to recover.")
    print("  'step' is signed -x, i.e. LARGER step -> recover (W1b direction).")
    print("-" * 78)
    band_report(rows, ['mu', 'nuhat', 'NU', 'eff', 'step', 'lamstar'],
                "all strata pooled", sign={'step': -1.0})

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1b: double-conditioned — NU band AND eff fixed.  This is the")
    print("  strict test: neither NU nor the bias strength can carry the AUC.")
    print("-" * 78)
    print(f"{'eff':>5} {'N in band':>10} {'rec':>9} | {'AUC mu':>9} {'p':>9} "
          f"| {'AUC NU':>9} | {'AUC step':>9}")
    tot_pairs = 0
    conc_mu = 0.0
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff
               and BAND_LO <= r['NU'] < BAND_HI]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>10} "
                  f"{str(len(pos)) + '/' + str(len(sub)):>9} |  (degenerate)")
            continue
        a_mu, p_mu = auc_p([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu, _ = auc_p([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st, _ = auc_p([-r['step'] for r in pos], [-r['step'] for r in neg])
        # Mantel-Haenszel style pooling of the concordance over strata.
        tot_pairs += len(pos) * len(neg)
        conc_mu += a_mu * len(pos) * len(neg)
        print(f"{eff:>5.2f} {len(sub):>10} "
              f"{str(len(pos)) + '/' + str(len(sub)):>9} | {a_mu:>9.4f} "
              f"{p_mu:>9.2e} | {a_nu:>9.4f} | {a_st:>9.4f}")
    if tot_pairs:
        print(f"\nstratum-pooled AUC(-mu | NU band, eff fixed) = "
              f"{conc_mu / tot_pairs:.4f}   ({tot_pairs} pos/neg pairs)")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1c: the converse — does NU separate inside a fixed mu?")
    print("  Curves are grouped by their (curve, eff) cell, in which mu is")
    print("  EXACTLY constant, so only NU (seed-driven) varies.")
    print("-" * 78)
    cells = {}
    for r in rows:
        cells.setdefault((r['n'], r['effq']), []).append(r)
    live = [(kk, v) for kk, v in cells.items()
            if 0 < sum(1 for r in v if r['ok']) < len(v)]
    tp = 0
    cn = 0.0
    cs = 0.0
    for _kk, v in live:
        pos = [r for r in v if r['ok']]
        neg = [r for r in v if not r['ok']]
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        a_st = auc([-r['step'] for r in pos], [-r['step'] for r in neg])
        w = len(pos) * len(neg)
        tp += w
        cn += a_nu * w
        cs += a_st * w
    print(f"{len(live)} mixed cells (of {len(cells)}), {tp} pos/neg pairs")
    if tp:
        print(f"  within-cell pooled AUC(-NU)   = {cn/tp:.4f}")
        print(f"  within-cell pooled AUC(+step) = {cs/tp:.4f}")
    print("  (mu is constant within a cell, so it has no within-cell AUC.)")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: two-parameter logistic fit on (log NU, log mu).")
    print("-" * 78)
    y = [1.0 if r['ok'] else 0.0 for r in rows]
    lNU = [math.log(max(r['NU'], 1e-12)) for r in rows]
    lMU = [math.log(max(r['mu'], 1e-12)) for r in rows]
    ll0 = None
    models = [
        ("intercept only", [[0.0] for _ in rows]),
        ("log NU", [[a] for a in lNU]),
        ("log mu", [[a] for a in lMU]),
        ("log NU + log mu", [[a, b] for a, b in zip(lNU, lMU)]),
    ]
    print(f"{'model':>18} {'loglik':>12} {'dev vs null':>12} {'AUC(fit)':>10}"
          f"   coefficients")
    for name, X in models:
        beta, ll = logistic_fit(X, y)
        if ll0 is None:
            ll0 = ll
        eta = [beta[0] + sum(b * v for b, v in zip(beta[1:], row))
               for row in X]
        pos = [e for e, yy in zip(eta, y) if yy > 0.5]
        neg = [e for e, yy in zip(eta, y) if yy < 0.5]
        a = auc([-e for e in pos], [-e for e in neg])
        cf = " ".join(f"{b:+.4f}" for b in beta)
        print(f"{name:>18} {ll:>12.2f} {2*(ll-ll0):>12.2f} {a:>10.4f}   {cf}")

    beta2, _ = logistic_fit([[a, b] for a, b in zip(lNU, lMU)], y)
    if abs(beta2[2]) > 1e-9:
        print(f"\ndecision boundary (p = 1/2):  "
              f"log mu = {-beta2[0]/beta2[2]:+.4f} "
              f"{-beta2[1]/beta2[2]:+.4f} * log NU")
        print(f"i.e.  mu * NU^{beta2[1]/beta2[2]:+.4f} = "
              f"{math.exp(-beta2[0]/beta2[2]):.4g}  is the 50% contour")

    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: the W1b two-block step as a wall predictor.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>9} | {'mean step':>10} "
          f"{'AUC step':>9} {'p':>9} | {'AUC mu':>8} {'AUC NU':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff
               and not math.isnan(r['step'])]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        ms = sum(r['step'] for r in sub) / len(sub)
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos)) + '/' + str(len(sub)):>9} | {ms:>10.4f} "
                  f"  (degenerate)")
            continue
        a_st, p_st = auc_p([-r['step'] for r in pos], [-r['step'] for r in neg])
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nu = auc([r['NU'] for r in pos], [r['NU'] for r in neg])
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(pos)) + '/' + str(len(sub)):>9} | {ms:>10.4f} "
              f"{a_st:>9.4f} {p_st:>9.2e} | {a_mu:>8.4f} {a_nu:>8.4f}")

    good = [r for r in rows if not math.isnan(r['step'])]
    gp = [r for r in good if r['ok']]
    gn = [r for r in good if not r['ok']]
    print(f"\npooled AUC(+step) = "
          f"{auc([-r['step'] for r in gp], [-r['step'] for r in gn]):.4f}   "
          f"(N={len(good)})")
    print(f"Spearman(step, log mu)  = "
          f"{spearman([r['step'] for r in good], [math.log(r['mu']) for r in good]):+.4f}")
    print(f"Spearman(step, log NU)  = "
          f"{spearman([r['step'] for r in good], [math.log(r['NU']) for r in good]):+.4f}")
    print(f"Spearman(step, eff)     = "
          f"{spearman([r['step'] for r in good], [r['eff'] for r in good]):+.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
