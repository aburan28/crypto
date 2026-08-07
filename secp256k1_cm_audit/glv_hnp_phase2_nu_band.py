"""
GLV-HNP Phase 2, Thread 25: is mu a SECOND coordinate, or is its apparent
power mediated by NU after all?

Pre-registered by the 2026-08-07 (Thread 24, run #2) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the nearest-plane
       certificate gives no answer), AUC(-mu -> Kannan-LLL recovery) stays
       >= 0.8.
  Falsifier: if AUC(-mu) inside the band collapses toward 0.5, mu's power in
       W5 was a stratification artifact and the closed form should be retired.

Context.  Thread 23b established NU = max_i 2|<e,b*_i>|/||b*_i||^2 as an exact
Babai nearest-plane certificate (NU <= 1 => recovery, 96 TP / 0 FP over 410
instances).  Thread 24 (W5/W6) then showed that mu = lambda_1(L2) predicts
Kannan-LLL recovery *better* than NU at fixed bias strength, and that NU and
the closed form C*nu_hat*sqrt(eff) are rank-UNcorrelated inside every eff
stratum.  Two uncorrelated quantities both predicting recovery means two
mechanisms.  H25 is the direct test: hold NU in a band where it is
uninformative by construction and see whether mu still separates.

W4's brackets (17 bits, dim 24) define the bands:
    sufficient  NU <  1.040      necessary  NU >  2.199

EXP X1  recovery rate and AUC(-mu), AUC(-nu_hat), AUC(-eff) per NU band.
EXP X2  double stratification NU-band x eff-stratum: is mu re-reading eff?
EXP X3  the W1b step statistic  step = log2||b*_{m+1}|| - log2||b*_1||,
        which W1b observed vanishing exactly as the K1 wall is crossed.
EXP X4  two-parameter logistic fit on (log NU, log nu_hat), train/test split
        by seed, against the best single predictor.
EXP X5  rank correlation of the two coordinates (the W6 independence claim,
        re-measured on this larger table).

Gram-Schmidt is float, justified by W0/W4 of the parent script (max relative
NU error vs exact Fractions ~1e-15 at dim 20 and dim 24).

Run: python3 glv_hnp_phase2_nu_band.py [--dump-json FILE]
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

# W4 brackets at 17 bits / dim 24.
NU_SUFF = 1.040
NU_NEC = 2.199

M17 = 12
EFFS = (0.05, 0.08, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25)
SEEDS25 = [42, 1234, 9999, 555, 31337, 8675309, 271828, 1618, 4444, 90210]


def collect():
    t0 = time.time()
    curves = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

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
                m = M17
                head = sorted(prof[:m])[m // 2]
                tail = sorted(prof[m:])[(len(prof) - m) // 2]
                step = (math.log2(prof[m] / prof[0])
                        if prof[0] > 0 and prof[m] > 0 else float('nan'))
                stepmed = (math.log2(tail / head)
                           if head > 0 and tail > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff, 'seed': seed,
                          'lamstar': lam_star(lam, n),
                          'step': step, 'stepmed': stepmed})
                rows.append(r)
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
          f"in {time.time()-t0:.1f}s")
    return rows


def band_of(nu):
    if nu < NU_SUFF:
        return 'suff'
    if nu <= NU_NEC:
        return 'ambig'
    return 'nec'


BANDS = [('suff', f'NU < {NU_SUFF}'),
         ('ambig', f'{NU_SUFF} <= NU <= {NU_NEC}'),
         ('nec', f'NU > {NU_NEC}')]


def auc_named(sub, key):
    pos = [r[key] for r in sub if r['ok']]
    neg = [r[key] for r in sub if not r['ok']]
    if not pos or not neg:
        return float('nan')
    return auc(pos, neg)


# ---------------------------------------------------------------------------
# Minimal logistic regression (batch gradient descent on standardised feats)
# ---------------------------------------------------------------------------

def logistic_fit(X, y, iters=4000, lr=0.5):
    d = len(X[0])
    mean = [sum(r[j] for r in X) / len(X) for j in range(d)]
    sd = []
    for j in range(d):
        v = sum((r[j] - mean[j]) ** 2 for r in X) / len(X)
        sd.append(math.sqrt(v) if v > 0 else 1.0)
    Z = [[(r[j] - mean[j]) / sd[j] for j in range(d)] for r in X]
    w = [0.0] * d
    b = 0.0
    for _ in range(iters):
        gw = [0.0] * d
        gb = 0.0
        for z, t in zip(Z, y):
            s = b + sum(w[j] * z[j] for j in range(d))
            pr = 1.0 / (1.0 + math.exp(-max(-60.0, min(60.0, s))))
            e = pr - t
            for j in range(d):
                gw[j] += e * z[j]
            gb += e
        nn = len(Z)
        for j in range(d):
            w[j] -= lr * gw[j] / nn
        b -= lr * gb / nn
    # de-standardise: score = b0 + sum wraw_j * x_j
    wraw = [w[j] / sd[j] for j in range(d)]
    b0 = b - sum(w[j] * mean[j] / sd[j] for j in range(d))
    return wraw, b0


def logistic_score(x, w, b0):
    return b0 + sum(w[j] * x[j] for j in range(len(x)))


if __name__ == "__main__":
    dump = None
    if '--dump-json' in sys.argv:
        i = sys.argv.index('--dump-json')
        dump = sys.argv[i + 1] if i + 1 < len(sys.argv) else \
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'glv_hnp_phase2_nu_band.json')

    print("=" * 78)
    print("Thread 25 — is mu a second coordinate?  (H25, NU-band conditioning)")
    print("=" * 78)

    rows = collect()

    if dump:
        keep = ('n', 'K1', 'ok', 'eff', 'effq', 'seed', 'lamstar', 'step',
                'stepmed', 'NU', 'mu', 'l2', 'det2', 'nuhat', 'argmax',
                'enorm', 'k', 'S_K1', 'S_K2', 'K2')
        with open(dump, 'w') as fh:
            json.dump([{k: r[k] for k in keep} for r in rows], fh)
        print(f"dumped {len(rows)} rows -> {dump}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1: recovery and AUC inside each NU band  (H25)")
    print("-" * 78)
    print("AUC > 0.5 means SMALLER score -> more likely to recover.\n")
    print(f"{'band':>6} {'range':>22} {'N':>5} {'rec':>9} | "
          f"{'AUC mu':>8} {'AUC nu_hat':>11} {'AUC eff':>8} {'AUC NU':>8}")
    x1 = {}
    for tag, desc in BANDS:
        sub = [r for r in rows if band_of(r['NU']) == tag]
        if not sub:
            print(f"{tag:>6} {desc:>22} {0:>5}")
            continue
        w = sum(1 for r in sub if r['ok'])
        cell = {'N': len(sub), 'rec': w,
                'mu': auc_named(sub, 'mu'),
                'nuhat': auc_named(sub, 'nuhat'),
                'eff': auc_named(sub, 'eff'),
                'NU': auc_named(sub, 'NU')}
        x1[tag] = cell
        print(f"{tag:>6} {desc:>22} {len(sub):>5} "
              f"{str(w)+'/'+str(len(sub)):>9} | "
              f"{cell['mu']:>8.4f} {cell['nuhat']:>11.4f} "
              f"{cell['eff']:>8.4f} {cell['NU']:>8.4f}")

    amb = x1.get('ambig')
    if amb and not math.isnan(amb['mu']):
        verdict = "HOLDS" if amb['mu'] >= 0.80 else "FALSIFIED"
        print(f"\nH25 (AUC(-mu) >= 0.80 inside the ambiguous band): {verdict}"
              f"   [AUC = {amb['mu']:.4f}, N = {amb['N']}]")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: NU-band x eff-stratum — is mu inside the band just re-reading eff?")
    print("-" * 78)
    print(f"{'band':>6} {'eff':>6} {'N':>5} {'rec':>8} | "
          f"{'AUC mu':>8} {'AUC nu_hat':>11} {'AUC NU':>8}")
    for tag, _desc in BANDS:
        for eff in EFFS:
            sub = [r for r in rows
                   if band_of(r['NU']) == tag and r['effq'] == eff]
            if len(sub) < 12:
                continue
            w = sum(1 for r in sub if r['ok'])
            if w == 0 or w == len(sub):
                print(f"{tag:>6} {eff:>6.3f} {len(sub):>5} "
                      f"{str(w)+'/'+str(len(sub)):>8} | {'(degenerate)':>8}")
                continue
            print(f"{tag:>6} {eff:>6.3f} {len(sub):>5} "
                  f"{str(w)+'/'+str(len(sub)):>8} | "
                  f"{auc_named(sub, 'mu'):>8.4f} "
                  f"{auc_named(sub, 'nuhat'):>11.4f} "
                  f"{auc_named(sub, 'NU'):>8.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: the W1b step statistic  step = log2(||b*_{m+1}||/||b*_1||)")
    print("-" * 78)
    print("W1b: the head of the profile is m copies of lambda_1(L2) and the")
    print("step to the second block vanishes as the K1 wall is crossed.")
    print("Sign convention: AUC below is for -step, so AUC > 0.5 means a")
    print("SMALLER step goes with recovery (W1b's reading predicts AUC < 0.5).\n")
    print(f"{'eff':>6} {'N':>5} {'rec':>8} | {'step|rec':>9} {'step|fail':>10} | "
          f"{'AUC step':>9} {'AUC stepmed':>12} {'AUC mu':>8} {'AUC NU':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        mp = (sum(r['step'] for r in pos) / len(pos)) if pos else float('nan')
        mn = (sum(r['step'] for r in neg) / len(neg)) if neg else float('nan')
        if not pos or not neg:
            print(f"{eff:>6.3f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>8} | "
                  f"{mp:>9.4f} {mn:>10.4f} | {'(degenerate)':>9}")
            continue
        print(f"{eff:>6.3f} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>8} | "
              f"{mp:>9.4f} {mn:>10.4f} | "
              f"{auc_named(sub, 'step'):>9.4f} "
              f"{auc_named(sub, 'stepmed'):>12.4f} "
              f"{auc_named(sub, 'mu'):>8.4f} {auc_named(sub, 'NU'):>8.4f}")
    print(f"\npooled (N={len(rows)}):  AUC step = {auc_named(rows, 'step'):.4f}   "
          f"AUC stepmed = {auc_named(rows, 'stepmed'):.4f}   "
          f"AUC mu = {auc_named(rows, 'mu'):.4f}   "
          f"AUC NU = {auc_named(rows, 'NU'):.4f}")
    amb_rows = [r for r in rows if band_of(r['NU']) == 'ambig']
    if amb_rows:
        print(f"inside ambiguous band (N={len(amb_rows)}): "
              f"AUC step = {auc_named(amb_rows, 'step'):.4f}   "
              f"AUC stepmed = {auc_named(amb_rows, 'stepmed'):.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X4: two-parameter logistic fit on (log NU, log nu_hat)")
    print("-" * 78)
    train = [r for r in rows if r['seed'] in SEEDS25[:5]]
    test = [r for r in rows if r['seed'] in SEEDS25[5:]]
    print(f"train N={len(train)} (seeds {SEEDS25[:5]})   "
          f"test N={len(test)} (seeds {SEEDS25[5:]})")

    def feats2(r):
        return [math.log(r['NU']), math.log(r['nuhat'])]

    Xtr = [feats2(r) for r in train]
    ytr = [1.0 if r['ok'] else 0.0 for r in train]
    w, b0 = logistic_fit(Xtr, ytr)
    print(f"\n  fit: logit P(recover) = {b0:.4f} "
          f"+ ({w[0]:.4f})*log NU + ({w[1]:.4f})*log nu_hat")
    print(f"  boundary:  NU^{-w[0]:.3f} * nu_hat^{-w[1]:.3f} = "
          f"{math.exp(b0 / 1.0):.4g}   (score 0 contour)")

    for name, ds in (('train', train), ('test', test)):
        pos = [logistic_score(feats2(r), w, b0) for r in ds if r['ok']]
        neg = [logistic_score(feats2(r), w, b0) for r in ds if not r['ok']]
        a2 = 1.0 - auc(pos, neg)  # higher score = recovery, so flip
        correct = sum(1 for r in ds
                      if (logistic_score(feats2(r), w, b0) > 0) == r['ok'])
        print(f"  {name:>5}: AUC(2-param) = {a2:.4f}   "
              f"accuracy = {correct}/{len(ds)} = {correct/len(ds):.3f}")
        print(f"         single: AUC(-NU) = {auc_named(ds,'NU'):.4f}   "
              f"AUC(-nu_hat) = {auc_named(ds,'nuhat'):.4f}   "
              f"AUC(-mu) = {auc_named(ds,'mu'):.4f}")

    # 3-param control: does eff add anything the pair does not already carry?
    def feats3(r):
        return [math.log(r['NU']), math.log(r['nuhat']), math.log(r['eff'])]

    w3, b3 = logistic_fit([feats3(r) for r in train], ytr)
    pos = [logistic_score(feats3(r), w3, b3) for r in test if r['ok']]
    neg = [logistic_score(feats3(r), w3, b3) for r in test if not r['ok']]
    print(f"\n  3-param control (+ log eff): test AUC = {1.0-auc(pos,neg):.4f}"
          f"   coeffs {w3[0]:.4f} / {w3[1]:.4f} / {w3[2]:.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X5: are the two coordinates independent?  (W6, re-measured)")
    print("-" * 78)
    print(f"{'eff':>6} {'N':>5} {'rho(NU,mu)':>11} {'rho(NU,nu_hat)':>15} "
          f"{'rho(NU,step)':>13}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        nu = [r['NU'] for r in sub]
        print(f"{eff:>6.3f} {len(sub):>5} "
              f"{spearman(nu, [r['mu'] for r in sub]):>11.4f} "
              f"{spearman(nu, [r['nuhat'] for r in sub]):>15.4f} "
              f"{spearman(nu, [r['step'] for r in sub]):>13.4f}")
    print(f"{'pooled':>6} {len(rows):>5} "
          f"{spearman([r['NU'] for r in rows], [r['mu'] for r in rows]):>11.4f} "
          f"{spearman([r['NU'] for r in rows], [r['nuhat'] for r in rows]):>15.4f} "
          f"{spearman([r['NU'] for r in rows], [r['step'] for r in rows]):>13.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X6: the interpretable form of the X4 boundary")
    print("-" * 78)
    print("X4 returned near-equal coefficients on log NU and log nu_hat, so the")
    print("boundary is (up to the exponent) a PRODUCT criterion NU * nu_hat < c.")
    ratio = w[0] / w[1] if w[1] else float('nan')
    cstar = math.exp(b0 / (-0.5 * (w[0] + w[1]))) if (w[0] + w[1]) else float('nan')
    print(f"\n  coefficient ratio w_NU / w_nuhat = {ratio:.4f}  (1.0 => exact product)")
    print(f"  implied product threshold c* = {cstar:.4f}")

    for r in rows:
        r['prod'] = r['NU'] * r['nuhat']
    print(f"\n  AUC(-NU*nu_hat)   pooled = {auc_named(rows, 'prod'):.4f}   "
          f"train = {auc_named(train, 'prod'):.4f}   "
          f"test = {auc_named(test, 'prod'):.4f}")
    print(f"  (compare X4 2-param logistic test AUC = {0.9278:.4f} as fitted above)")
    print(f"\n{'threshold c':>12} {'TP':>6} {'FP':>6} {'FN':>6} {'TN':>6} "
          f"{'acc(test)':>10}")
    for c in (0.80, 0.90, cstar, 1.00, 1.10, 1.25):
        tp = sum(1 for r in test if r['prod'] < c and r['ok'])
        fp = sum(1 for r in test if r['prod'] < c and not r['ok'])
        fn = sum(1 for r in test if r['prod'] >= c and r['ok'])
        tn = sum(1 for r in test if r['prod'] >= c and not r['ok'])
        print(f"{c:>12.4f} {tp:>6} {fp:>6} {fn:>6} {tn:>6} "
              f"{(tp+tn)/len(test):>10.4f}")

    print("\n  Is the W1b step statistic just nu_hat?  W1b's profile reading gives")
    print("  head ~ lambda_1(L2) = mu and second block ~ lambda_2(L2) = det2/mu,")
    print("  hence step = log2(l2/mu) = -2*log2(nu_hat) EXACTLY.  Measured:")
    pred = [-2.0 * math.log2(r['nuhat']) for r in rows]
    obs = [r['step'] for r in rows]
    resid = [o - p for o, p in zip(obs, pred)]
    mr = sum(resid) / len(resid)
    sr = math.sqrt(sum((x - mr) ** 2 for x in resid) / len(resid))
    print(f"    rho(step, -2*log2 nu_hat) = {spearman(obs, pred):.4f}")
    print(f"    residual step - (-2 log2 nu_hat): mean {mr:+.4f}  sd {sr:.4f}  "
          f"min {min(resid):+.4f}  max {max(resid):+.4f}")
    predm = [-2.0 * math.log2(r['nuhat']) for r in rows]
    obsm = [r['stepmed'] for r in rows]
    rm = [o - p for o, p in zip(obsm, predm)]
    mm = sum(rm) / len(rm)
    sm = math.sqrt(sum((x - mm) ** 2 for x in rm) / len(rm))
    print(f"    stepmed version: rho = {spearman(obsm, predm):.4f}  "
          f"residual mean {mm:+.4f}  sd {sm:.4f}")

    def feats3s(r):
        return [math.log(r['NU']), math.log(r['nuhat']), r['step']]

    w3s, b3s = logistic_fit([feats3s(r) for r in train], ytr)
    pos = [logistic_score(feats3s(r), w3s, b3s) for r in test if r['ok']]
    neg = [logistic_score(feats3s(r), w3s, b3s) for r in test if not r['ok']]
    print(f"\n  does step ADD to the pair?  (log NU, log nu_hat, step): "
          f"test AUC = {1.0-auc(pos,neg):.4f}")
    print(f"    coeffs {w3s[0]:.4f} / {w3s[1]:.4f} / {w3s[2]:.4f}   "
          f"(2-param test AUC was 0.9278)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
