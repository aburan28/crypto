"""
GLV-HNP Phase 2, Thread 25: is mu a SECOND coordinate, or is its apparent
power mediated by NU?

Pre-registered by the 2026-08-07 run #2 (Thread 24) log entry:

  H25: within the ambiguous band  1.04 <= NU <= 2.20  (where the nearest-plane
       guarantee gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  If yes -> mu is a genuine second coordinate and (NU, mu) is the 2-parameter
            viability test; deliverable is a logistic fit on (log NU, log mu)
            with its decision boundary.
  If no  -> mu's power is entirely mediated by NU, the W5 result is a
            stratification artifact, and the closed form should be retired.

Why.  Thread 24 W5/W6 established that at FIXED bias strength eff = K1*K2/n
the exact BDD certificate NU and the lattice-free statistic nu_hat = mu/sqrt(det L2)
are mutually uncorrelated (Spearman -0.28..+0.16) yet both predict recovery.
Two uncorrelated predictors of one event means two mechanisms.  Conditioning
on NU is the direct way to see whether the second one is real.

Secondary (W1b follow-up, also pre-registered): the GS profile of L0 is two
blocks, the head being m exact copies of lambda_1(L2).  Define

    step = log2(||b*_{m+1}||) - log2(||b*_1||)

and test whether step -> 0 predicts the LLL wall better than NU or mu.

Numerics.  Float GS, justified by Thread 24 W0 (max relative NU error vs exact
Fractions 6.6e-16 at dim 24 / 17 bits).

Run:  python3 glv_hnp_phase2_thread25.py            # builds + caches the table
      python3 glv_hnp_phase2_thread25.py --rebuild  # ignore the cache
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

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "glv_hnp_phase2_thread25_table.json")

# Thread 24 W5 grid, with the seed list doubled so the NU bands have enough
# instances to condition on.  eff = K1*K2/n is the bias strength.
BITS = 17
M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
SEEDS25 = [42, 1234, 9999, 555, 31337, 271828, 141421, 173205, 224744, 264575]

# The Thread 24 W4 bracket at 17 bits: NU < 1.040 sufficient, NU > 2.199
# necessary.  Everything between is where nearest-plane says nothing.
BAND_LO, BAND_HI = 1.040, 2.199


# ---------------------------------------------------------------------------
# Table construction
# ---------------------------------------------------------------------------

def build_table(bits=BITS, m=M17, effs=EFFS, seeds=SEEDS25):
    t0 = time.time()
    curves = search_curves(1 << (bits - 1), 1 << bits, per_bin=2, nbins=10)
    print(f"{len(curves)} {bits}-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    M17, EFFS, SEEDS25 = m, effs, seeds
    t0 = time.time()
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves:
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
                # W1b two-block step: head is m copies of lambda_1(L2).
                step = (math.log2(prof[M17]) - math.log2(prof[0])
                        if prof[0] > 0 and prof[M17] > 0 else float('nan'))
                rows.append({
                    'n': n, 'lam': lam, 'K1': k1b, 'K2': k2b, 'seed': seed,
                    'effq': eff, 'eff': k1b * k2b / n,
                    'ok': bool(rk['ok']),
                    'NU': r['NU'], 'mu': r['mu'], 'nuhat': r['nuhat'],
                    'l2': r['l2'], 'det2': r['det2'], 'argmax': r['argmax'],
                    'lamstar': lam_star(lam, n), 'step': step,
                    'k': r['k'], 'enorm': r['enorm'],
                })
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    return rows


def load_table(rebuild=False, cache=CACHE, **kw):
    if not rebuild and os.path.exists(cache):
        with open(cache) as f:
            rows = json.load(f)
        print(f"loaded {len(rows)} cached instances from "
              f"{os.path.basename(cache)}")
        return rows
    rows = build_table(**kw)
    with open(cache, "w") as f:
        json.dump(rows, f)
    print(f"dumped table -> {os.path.basename(cache)}")
    return rows


# ---------------------------------------------------------------------------
# Logistic regression (IRLS, tiny ridge for separation stability)
# ---------------------------------------------------------------------------

def logistic_fit(X, y, ridge=1e-6, iters=60):
    """X: list of feature rows (no intercept column).  y: 0/1 list.
    Returns (beta including intercept at index 0, loglik)."""
    d = len(X[0]) + 1
    A = [[1.0] + list(r) for r in X]
    beta = [0.0] * d
    for _ in range(iters):
        g = [0.0] * d
        H = [[0.0] * d for _ in range(d)]
        for a, yi in zip(A, y):
            z = sum(bj * aj for bj, aj in zip(beta, a))
            z = max(-30.0, min(30.0, z))
            pr = 1.0 / (1.0 + math.exp(-z))
            w = max(pr * (1 - pr), 1e-9)
            for i in range(d):
                g[i] += (yi - pr) * a[i]
                for j in range(d):
                    H[i][j] += w * a[i] * a[j]
        for i in range(d):
            g[i] -= ridge * beta[i]
            H[i][i] += ridge
        delta = solve(H, g)
        if delta is None:
            break
        beta = [b + dl for b, dl in zip(beta, delta)]
        if max(abs(x) for x in delta) < 1e-10:
            break
    ll = 0.0
    for a, yi in zip(A, y):
        z = max(-30.0, min(30.0, sum(bj * aj for bj, aj in zip(beta, a))))
        pr = 1.0 / (1.0 + math.exp(-z))
        ll += yi * math.log(max(pr, 1e-12)) + (1 - yi) * math.log(max(1 - pr, 1e-12))
    return beta, ll


def solve(Amat, bvec):
    n = len(bvec)
    M = [row[:] + [bvec[i]] for i, row in enumerate(Amat)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-14:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        for j in range(c, n + 1):
            M[c][j] /= pv
        for r in range(n):
            if r == c:
                continue
            f = M[r][c]
            if f:
                for j in range(c, n + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][n] for i in range(n)]


def predict(beta, X):
    out = []
    for r in X:
        a = [1.0] + list(r)
        z = max(-30.0, min(30.0, sum(bj * aj for bj, aj in zip(beta, a))))
        out.append(1.0 / (1.0 + math.exp(-z)))
    return out


def auc_score(scores, y):
    """AUC where LARGER score -> more likely y=1."""
    pos = [s for s, t in zip(scores, y) if t]
    neg = [s for s, t in zip(scores, y) if not t]
    if not pos or not neg:
        return float('nan')
    conc = sum((1.0 if a > b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return conc / (len(pos) * len(neg))


def acc_at_half(pr, y):
    return sum(1 for p, t in zip(pr, y) if (p >= 0.5) == bool(t)) / len(y)


# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------

def exp_x1(rows):
    print("\n" + "-" * 78)
    print("EXP X1 (H25): AUC(-mu -> recovery) INSIDE NU bands.")
    print("        NU is held (approximately) fixed, so any separation left")
    print("        is a second mechanism.  AUC > 0.5 => smaller mu recovers.")
    print("-" * 78)
    bands = [
        ("NU < 1.040  (nearest-plane SUFFICIENT)", lambda r: r['NU'] < BAND_LO),
        ("1.040 <= NU <= 2.199  (AMBIGUOUS)",
         lambda r: BAND_LO <= r['NU'] <= BAND_HI),
        ("NU > 2.199  (nearest-plane says nothing)", lambda r: r['NU'] > BAND_HI),
    ]
    print(f"{'band':<42} {'N':>5} {'rec':>9} | {'AUC mu':>8} {'AUC NU':>8} "
          f"{'AUC nuhat':>10} {'AUC step':>9}")
    for label, pred in bands:
        sub = [r for r in rows if pred(r)]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        rec = f"{len(pos)}/{len(sub)}"
        if not pos or not neg:
            print(f"{label:<42} {len(sub):>5} {rec:>9} |  (degenerate)")
            continue
        print(f"{label:<42} {len(sub):>5} {rec:>9} | "
              f"{auc([r['mu'] for r in pos], [r['mu'] for r in neg]):>8.4f} "
              f"{auc([r['NU'] for r in pos], [r['NU'] for r in neg]):>8.4f} "
              f"{auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg]):>10.4f} "
              f"{auc([r['step'] for r in pos], [r['step'] for r in neg]):>9.4f}")

    # finer conditioning: NU deciles, so the band width cannot hide a trend
    print("\nfiner: NU deciles (band edges from the data, not the bracket)")
    srt = sorted(rows, key=lambda r: r['NU'])
    q = 10
    print(f"{'decile':>7} {'NU range':>17} {'N':>5} {'rec':>9} | "
          f"{'AUC mu':>8} {'AUC nuhat':>10} {'AUC step':>9}")
    for i in range(q):
        sub = srt[i * len(srt) // q:(i + 1) * len(srt) // q]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        rng = f"{sub[0]['NU']:.3f}-{sub[-1]['NU']:.3f}"
        rec = f"{len(pos)}/{len(sub)}"
        if not pos or not neg:
            print(f"{i+1:>7} {rng:>17} {len(sub):>5} {rec:>9} |  (degenerate)")
            continue
        print(f"{i+1:>7} {rng:>17} {len(sub):>5} {rec:>9} | "
              f"{auc([r['mu'] for r in pos], [r['mu'] for r in neg]):>8.4f} "
              f"{auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg]):>10.4f} "
              f"{auc([r['step'] for r in pos], [r['step'] for r in neg]):>9.4f}")


def exp_x2(rows):
    print("\n" + "-" * 78)
    print("EXP X2: logistic fit on (log NU, log mu) and the nested-model")
    print("        likelihood-ratio test.  LR = 2*(ll_full - ll_reduced),")
    print("        chi2 with 1 df: 3.84 at p=0.05, 10.83 at p=0.001.")
    print("-" * 78)
    y = [1 if r['ok'] else 0 for r in rows]
    lnu = [math.log(r['NU']) for r in rows]
    lmu = [math.log(r['mu']) for r in rows]
    lnh = [math.log(r['nuhat']) for r in rows]

    models = {
        'NU only': [[a] for a in lnu],
        'mu only': [[a] for a in lmu],
        'nuhat only': [[a] for a in lnh],
        'NU + mu': [[a, b] for a, b in zip(lnu, lmu)],
        'NU + nuhat': [[a, b] for a, b in zip(lnu, lnh)],
    }
    fits = {}
    print(f"{'model':<12} {'loglik':>10} {'AUC':>7} {'acc':>7}   coefficients")
    for name, X in models.items():
        beta, ll = logistic_fit(X, y)
        pr = predict(beta, X)
        fits[name] = (beta, ll)
        cs = "  ".join(f"{b:+.4f}" for b in beta)
        print(f"{name:<12} {ll:>10.2f} {auc_score(pr, y):>7.4f} "
              f"{acc_at_half(pr, y):>7.4f}   [{cs}]")

    lr_mu = 2 * (fits['NU + mu'][1] - fits['NU only'][1])
    lr_nu = 2 * (fits['NU + mu'][1] - fits['mu only'][1])
    print(f"\nLR(add log mu    to NU-only) = {lr_mu:8.2f}   (1 df)")
    print(f"LR(add log NU    to mu-only) = {lr_nu:8.2f}   (1 df)")
    print(f"LR(add log nuhat to NU-only) = "
          f"{2*(fits['NU + nuhat'][1] - fits['NU only'][1]):8.2f}   (1 df)")
    print(f"LR(add log NU to nuhat-only) = "
          f"{2*(fits['NU + nuhat'][1] - fits['nuhat only'][1]):8.2f}   (1 df)")

    beta = fits['NU + mu'][0]
    if abs(beta[2]) > 1e-12:
        print(f"\ndecision boundary (p=0.5):  "
              f"{beta[0]:+.4f} {beta[1]:+.4f}*log NU {beta[2]:+.4f}*log mu = 0")
        print(f"  i.e.  log mu = {-beta[0]/beta[2]:+.4f} "
              f"{-beta[1]/beta[2]:+.4f}*log NU")
        print(f"  i.e.  NU^{beta[1]/-beta[2]:+.4f} * mu  =  "
              f"{math.exp(-beta[0]/beta[2]):.4g}   (recover below)")


def exp_x3(rows):
    print("\n" + "-" * 78)
    print("EXP X3: the W1b two-block step as a predictor.")
    print("        step = log2||b*_{m+1}|| - log2||b*_1||;  W1b saw it vanish")
    print("        exactly as the K1 wall is crossed.  AUC > 0.5 => SMALLER")
    print("        step recovers, which is the sign W1b predicts.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>9} | {'AUC step':>9} {'AUC NU':>8} "
          f"{'AUC mu':>8} | {'mean step':>10} {'rho(step,NU)':>13}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        ms = sum(r['step'] for r in sub) / len(sub)
        rho = spearman([r['step'] for r in sub], [r['NU'] for r in sub])
        rec = f"{len(pos)}/{len(sub)}"
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} {rec:>9} |  (degenerate)"
                  f"{'':>18} {ms:>10.4f} {rho:>13.4f}")
            continue
        print(f"{eff:>5.2f} {len(sub):>5} {rec:>9} | "
              f"{auc([r['step'] for r in pos], [r['step'] for r in neg]):>9.4f} "
              f"{auc([r['NU'] for r in pos], [r['NU'] for r in neg]):>8.4f} "
              f"{auc([r['mu'] for r in pos], [r['mu'] for r in neg]):>8.4f} | "
              f"{ms:>10.4f} {rho:>13.4f}")
    pos = [r for r in rows if r['ok']]
    neg = [r for r in rows if not r['ok']]
    print(f"\npooled (N={len(rows)}):  AUC step = "
          f"{auc([r['step'] for r in pos], [r['step'] for r in neg]):.4f}   "
          f"AUC NU = {auc([r['NU'] for r in pos], [r['NU'] for r in neg]):.4f}   "
          f"AUC mu = {auc([r['mu'] for r in pos], [r['mu'] for r in neg]):.4f}")
    print(f"pooled Spearman(step, NU) = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):+.4f}   "
          f"Spearman(step, mu) = "
          f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):+.4f}")


def exp_x4(rows):
    print("\n" + "-" * 78)
    print("EXP X4: mediation control — is mu just a proxy for NU?")
    print("        If Spearman(mu, NU) ~ 0 at fixed eff (Thread 24 W6) AND mu")
    print("        separates inside NU bands (X1), the two are independent.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} | {'rho(mu,NU)':>11} {'rho(nuhat,NU)':>14} "
          f"{'rho(mu,nuhat)':>14}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        print(f"{eff:>5.2f} {len(sub):>5} | "
              f"{spearman([r['mu'] for r in sub], [r['NU'] for r in sub]):>11.4f} "
              f"{spearman([r['nuhat'] for r in sub], [r['NU'] for r in sub]):>14.4f} "
              f"{spearman([r['mu'] for r in sub], [r['nuhat'] for r in sub]):>14.4f}")
    print(f"{'ALL':>5} {len(rows):>5} | "
          f"{spearman([r['mu'] for r in rows], [r['NU'] for r in rows]):>11.4f} "
          f"{spearman([r['nuhat'] for r in rows], [r['NU'] for r in rows]):>14.4f} "
          f"{spearman([r['mu'] for r in rows], [r['nuhat'] for r in rows]):>14.4f}")

    # NU-band x eff cell table: the strongest form of the control.
    print("\nAUC(-mu) in each (eff, NU-band) cell.  '.' = degenerate cell.")
    print(f"{'eff':>5} | {'NU<1.04':>16} {'1.04-2.20':>16} {'NU>2.20':>16}")
    for eff in EFFS:
        cells = []
        for lo, hi in ((0.0, BAND_LO), (BAND_LO, BAND_HI), (BAND_HI, 1e18)):
            sub = [r for r in rows
                   if r['effq'] == eff and lo <= r['NU'] < hi]
            pos = [r for r in sub if r['ok']]
            neg = [r for r in sub if not r['ok']]
            if not pos or not neg:
                cells.append(f"{'.':>7} ({len(sub):>3})")
            else:
                a = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
                cells.append(f"{a:>7.4f} ({len(sub):>3})")
        print(f"{eff:>5.2f} | {cells[0]:>16} {cells[1]:>16} {cells[2]:>16}")


def exp_x5(rows):
    print("\n" + "-" * 78)
    print("EXP X5: the X2 boundary has near-equal exponents in log NU and")
    print("        log nuhat, i.e. it collapses to a PRODUCT criterion.")
    print("        Test  P = NU * nuhat  directly, unfitted.")
    print("-" * 78)
    y = [1 if r['ok'] else 0 for r in rows]
    lnu = [math.log(r['NU']) for r in rows]
    lnh = [math.log(r['nuhat']) for r in rows]
    beta, ll = logistic_fit([[a, b] for a, b in zip(lnu, lnh)], y)
    ratio = beta[2] / beta[1] if beta[1] else float('nan')
    print(f"free fit   : {beta[0]:+.4f} {beta[1]:+.4f}*log NU "
          f"{beta[2]:+.4f}*log nuhat      loglik {ll:.2f}")
    print(f"  exponent ratio  beta_nuhat/beta_NU = {ratio:.4f}   "
          f"(1.0000 => exact product)")
    thr_free = math.exp(-beta[0] / beta[1]) if beta[1] else float('nan')
    print(f"  boundary   NU * nuhat^{ratio:.4f} = {thr_free:.4f}")

    lp = [a + b for a, b in zip(lnu, lnh)]           # log(NU*nuhat)
    beta1, ll1 = logistic_fit([[a] for a in lp], y)
    thr = math.exp(-beta1[0] / beta1[1]) if beta1[1] else float('nan')
    print(f"\nconstrained (exponents forced equal):")
    print(f"  {beta1[0]:+.4f} {beta1[1]:+.4f}*log(NU*nuhat)      "
          f"loglik {ll1:.2f}   LR vs free = {2*(ll-ll1):.2f} (1 df)")
    print(f"  fitted threshold  NU*nuhat = {thr:.4f}")

    P = [math.exp(a) for a in lp]
    pos = [p for p, t in zip(P, y) if t]
    neg = [p for p, t in zip(P, y) if not t]
    print(f"\n  AUC(-NU*nuhat)          = {auc(pos, neg):.4f}")
    print(f"  AUC(-NU) alone          = "
          f"{auc([r['NU'] for r in rows if r['ok']], [r['NU'] for r in rows if not r['ok']]):.4f}")
    print(f"  AUC(-nuhat) alone       = "
          f"{auc([r['nuhat'] for r in rows if r['ok']], [r['nuhat'] for r in rows if not r['ok']]):.4f}")

    print(f"\n  confusion at the unfitted cut  NU*nuhat < 1:")
    tp = sum(1 for p, t in zip(P, y) if p < 1 and t)
    fp = sum(1 for p, t in zip(P, y) if p < 1 and not t)
    fn = sum(1 for p, t in zip(P, y) if p >= 1 and t)
    tn = sum(1 for p, t in zip(P, y) if p >= 1 and not t)
    print(f"    TP {tp:>4}  FP {fp:>4}  FN {fn:>4}  TN {tn:>4}   "
          f"acc {(tp+tn)/len(y):.4f}  precision "
          f"{tp/(tp+fp) if tp+fp else float('nan'):.4f}")
    print(f"\n  P = NU*nuhat  |  success: min {min(pos):.4f}  median "
          f"{sorted(pos)[len(pos)//2]:.4f}  max {max(pos):.4f}")
    print(f"                |  failure: min {min(neg):.4f}  median "
          f"{sorted(neg)[len(neg)//2]:.4f}  max {max(neg):.4f}")
    print(f"  bracket : sufficient P < {min(neg):.4f} , "
          f"necessary P > {max(pos):.4f}")

    print(f"\n  per-eff-stratum AUC(-P), eff held FIXED:")
    print(f"  {'eff':>5} {'N':>5} {'rec':>9} {'AUC P':>8} {'AUC NU':>8} "
          f"{'AUC nuhat':>10}")
    for eff in sorted({r['effq'] for r in rows}):
        sub = [r for r in rows if r['effq'] == eff]
        p_ = [r for r in sub if r['ok']]
        n_ = [r for r in sub if not r['ok']]
        rec = f"{len(p_)}/{len(sub)}"
        if not p_ or not n_:
            print(f"  {eff:>5.2f} {len(sub):>5} {rec:>9} {'(degenerate)':>8}")
            continue
        pr = lambda k, g: [x[k] * x['nuhat'] if k == 'NU_' else x[k] for x in g]
        print(f"  {eff:>5.2f} {len(sub):>5} {rec:>9} "
              f"{auc([x['NU']*x['nuhat'] for x in p_], [x['NU']*x['nuhat'] for x in n_]):>8.4f} "
              f"{auc([x['NU'] for x in p_], [x['NU'] for x in n_]):>8.4f} "
              f"{auc([x['nuhat'] for x in p_], [x['nuhat'] for x in n_]):>10.4f}")


def exp_x6(rows):
    print("\n" + "-" * 78)
    print("EXP X6: does `step` add anything to (NU, nuhat)?  The head of the")
    print("        profile IS lambda_1(L2) (W1b), so step = log2||b*_{m+1}||")
    print("        - log2(mu) is partly a reparametrisation of mu.")
    print("-" * 78)
    y = [1 if r['ok'] else 0 for r in rows]
    lnu = [math.log(r['NU']) for r in rows]
    lnh = [math.log(r['nuhat']) for r in rows]
    st = [r['step'] for r in rows]
    _, ll2 = logistic_fit([[a, b] for a, b in zip(lnu, lnh)], y)
    b3, ll3 = logistic_fit([[a, b, c] for a, b, c in zip(lnu, lnh, st)], y)
    print(f"  loglik (NU, nuhat)        = {ll2:.2f}")
    print(f"  loglik (NU, nuhat, step)  = {ll3:.2f}   "
          f"LR = {2*(ll3-ll2):.2f} (1 df; 3.84 at p=0.05)")
    print(f"  coefficients: [{'  '.join(f'{x:+.4f}' for x in b3)}]")
    print(f"\n  Spearman(step, log mu)  = "
          f"{spearman(st, [math.log(r['mu']) for r in rows]):+.4f}")
    print(f"  Spearman(step, log nuhat) = "
          f"{spearman(st, lnh):+.4f}")
    print(f"  step is monotone in eff (X3 mean column), and its zero crossing")
    print(f"  is NOT the wall: mean step is ~0 in the EASIEST stratum here,")
    print(f"  while W1b's 12-bit example had step ~0 in the HARDEST one.")


if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, or mediated by NU?")
    print("=" * 78)
    rows = load_table(rebuild="--rebuild" in sys.argv)
    npos = sum(1 for r in rows if r['ok'])
    print(f"\n{len(rows)} instances, {npos} recover ({npos/len(rows):.1%}), "
          f"dim {rows[0]['k']}, {BITS} bits, m={M17}")
    print(f"NU range {min(r['NU'] for r in rows):.3f} .. "
          f"{max(r['NU'] for r in rows):.3f}")

    exp_x1(rows)
    exp_x2(rows)
    exp_x3(rows)
    exp_x4(rows)
    exp_x5(rows)
    exp_x6(rows)

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
