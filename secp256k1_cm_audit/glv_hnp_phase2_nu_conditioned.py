"""
GLV-HNP Phase 2, Thread 25: is mu a second coordinate, or is its apparent
power mediated by NU after all?

Pre-registered by the 2026-08-07 (Thread 24) log entry, verbatim:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  If yes, mu is a genuine second coordinate and the pair (NU, mu) is the
  2-parameter viability test Phase 2 has been missing; the deliverable is a
  logistic fit on (log NU, log mu) with its decision boundary.  If no -- if
  mu's apparent power is entirely mediated by NU after all -- then the W5
  result is a stratification artifact and the closed form should be retired.

  Secondary: define step = log2(||b*_{m+1}||) - log2(||b*_1||) and test
  whether step -> 0 predicts the wall better than either NU or mu.

Why the pre-registered test is not sufficient on its own.  Conditioning on
the NU band does NOT control the bias strength eff = K1*K2/n, and mu depends
on eff through S_K1 = n//K1.  A pooled in-band AUC could therefore be bought
entirely by eff, exactly the failure mode W3 fell into and W5 was built to
avoid.  X3 below repeats the test doubly conditioned on (eff, NU band), which
is the honest version; X2 reports the pre-registered pooled number for the
record and is read together with X3, not instead of it.

Input: glv_hnp_phase2_strat_rows.json, the 500-instance 17-bit table dumped
by `glv_hnp_phase2_gsprofile_strat.py --dump-json`.  No new lattice work.
The dump is ~400 KB and regenerates deterministically in ~3 s, so it is not
committed; if it is absent this script rebuilds it.

Run: python3 glv_hnp_phase2_nu_conditioned.py [rows.json]
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

# The Thread 24 W4 bracket at 17 bits, measured on the 300-instance U2 grid.
BAND_LO, BAND_HI = 1.040, 2.199
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)


# ---------------------------------------------------------------------------
# derived per-row quantities
# ---------------------------------------------------------------------------

def enrich(rows):
    for r in rows:
        prof = [x for x in r['prof']]
        lg = [math.log2(x) if x > 0 else float('nan') for x in prof]
        m = r['m']
        r['lgprof'] = lg
        # the pre-registered secondary statistic
        r['step'] = lg[m] - lg[0]
        # head-block vs tail-block drop, a smoothed version of the same idea
        r['drop'] = sum(lg[:m]) / m - sum(lg[m:]) / (len(lg) - m)
        # GSA slope: least-squares slope of log2||b*_i|| against i
        k = len(lg)
        xm = (k - 1) / 2.0
        ym = sum(lg) / k
        sxy = sum((i - xm) * (lg[i] - ym) for i in range(k))
        sxx = sum((i - xm) ** 2 for i in range(k))
        r['gsa'] = sxy / sxx
        r['head_over_mu'] = prof[0] / r['mu']
    return rows


def auc_pairs(pos, neg):
    """Concordant pairs and total pairs for the score -x."""
    if not pos or not neg:
        return 0.0, 0
    conc = sum((1.0 if a < b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return conc, len(pos) * len(neg)


def stratified_auc(cells, key):
    """Pool concordance over cells: sum(concordant) / sum(pairs)."""
    c = t = 0.0
    for sub in cells:
        pos = [r[key] for r in sub if r['ok']]
        neg = [r[key] for r in sub if not r['ok']]
        a, b = auc_pairs(pos, neg)
        c += a
        t += b
    return (c / t if t else float('nan')), int(t)


def fmt_auc(sub, key):
    pos = [r[key] for r in sub if r['ok']]
    neg = [r[key] for r in sub if not r['ok']]
    if not pos or not neg:
        return "     --"
    return f"{auc(pos, neg):>7.4f}"


# ---------------------------------------------------------------------------
# logistic regression (Newton-Raphson with ridge), pure python
# ---------------------------------------------------------------------------

def _solve(A, b):
    k = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(k):
        piv = max(range(c, k), key=lambda i: abs(M[i][c]))
        if abs(M[piv][c]) < 1e-14:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        for j in range(c, k + 1):
            M[c][j] /= pv
        for i in range(k):
            if i != c and M[i][c]:
                f = M[i][c]
                for j in range(c, k + 1):
                    M[i][j] -= f * M[c][j]
    return [M[i][k] for i in range(k)]


def logistic_fit(X, y, ridge=1e-3, iters=60):
    """X includes no intercept column; one is prepended.  Returns beta."""
    n, d = len(X), len(X[0]) + 1
    Z = [[1.0] + list(x) for x in X]
    beta = [0.0] * d
    for _ in range(iters):
        g = [0.0] * d
        H = [[0.0] * d for _ in range(d)]
        for i in range(n):
            z = Z[i]
            eta = sum(beta[j] * z[j] for j in range(d))
            eta = max(-30.0, min(30.0, eta))
            pr = 1.0 / (1.0 + math.exp(-eta))
            w = max(pr * (1 - pr), 1e-9)
            res = y[i] - pr
            for j in range(d):
                g[j] += res * z[j]
                for l in range(d):
                    H[j][l] -= w * z[j] * z[l]
        for j in range(d):
            g[j] -= ridge * beta[j]
            H[j][j] -= ridge
        step = _solve([[-h for h in row] for row in H], g)
        if step is None:
            break
        beta = [beta[j] + step[j] for j in range(d)]
        if max(abs(s) for s in step) < 1e-10:
            break
    return beta


def logistic_score(beta, x):
    eta = beta[0] + sum(beta[j + 1] * x[j] for j in range(len(x)))
    return 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, eta))))


def grouped_cv_auc(rows, feats, ngroup=5):
    """Group k-fold by curve, so no fold shares a curve with its training set.
    AUC here is of the fitted probability, so LARGER = recovery: the sign
    convention is flipped relative to auc(), hence the 1 - auc(...)."""
    curves = sorted({r['n'] for r in rows})
    fold = {n: i % ngroup for i, n in enumerate(curves)}
    preds, ys = [], []
    for f in range(ngroup):
        tr = [r for r in rows if fold[r['n']] != f]
        te = [r for r in rows if fold[r['n']] == f]
        if not te or not tr:
            continue
        beta = logistic_fit([[fe(r) for fe in feats] for r in tr],
                            [1.0 if r['ok'] else 0.0 for r in tr])
        for r in te:
            preds.append(logistic_score(beta, [fe(r) for fe in feats]))
            ys.append(r['ok'])
    pos = [p for p, o in zip(preds, ys) if o]
    neg = [p for p, o in zip(preds, ys) if not o]
    return 1.0 - auc(pos, neg), len(preds)


FEATS = {
    'logNU':    lambda r: math.log(r['NU']),
    'lognuhat': lambda r: math.log(r['nuhat']),
    'logmu':    lambda r: math.log(r['mu']),
    'logeff':   lambda r: math.log(r['eff']),
    'step':     lambda r: r['step'],
}


if __name__ == "__main__":
    path = (sys.argv[1] if len(sys.argv) > 1 else
            os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'glv_hnp_phase2_strat_rows.json'))
    if not os.path.exists(path):
        print(f"{os.path.basename(path)} absent — rebuilding "
              f"(deterministic, ~3s)")
        from glv_hnp_phase2_gsprofile_strat import build_rows, dump_rows
        dump_rows(build_rows(verbose=False), path)
    rows = enrich(json.load(open(path)))
    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, or is it mediated by NU?")
    print("=" * 78)
    print(f"\n{len(rows)} rows from {os.path.basename(path)} "
          f"(17 bits, dim {rows[0]['k']}, m={rows[0]['m']}), "
          f"{sum(1 for r in rows if r['ok'])} recoveries")

    # ----------------------------------------------------------------- X1
    print("\n" + "-" * 78)
    print("EXP X1: the NU bracket on THIS table (the W4 bracket came from the")
    print("        300-instance U2 grid; check it transfers)")
    print("-" * 78)
    pos = [r['NU'] for r in rows if r['ok']]
    neg = [r['NU'] for r in rows if not r['ok']]
    pos.sort()
    neg.sort()
    print(f"  NU | success : min {pos[0]:.3f}  median {pos[len(pos)//2]:.3f} "
          f"  max {pos[-1]:.3f}  (n={len(pos)})")
    print(f"  NU | failure : min {neg[0]:.3f}  median {neg[len(neg)//2]:.3f} "
          f"  max {neg[-1]:.3f}  (n={len(neg)})")
    tp = sum(1 for r in rows if r['ok'] and r['NU'] <= 1.0)
    fp = sum(1 for r in rows if not r['ok'] and r['NU'] <= 1.0)
    print(f"  NU <= 1 : TP {tp}  FP {fp}")
    print(f"  measured bracket : sufficient NU < {neg[0]:.3f} , "
          f"necessary NU > {pos[-1]:.3f}")
    print(f"  pre-registered   : [{BAND_LO}, {BAND_HI}]  (used below)")

    inband = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    print(f"  in-band rows: {len(inband)} "
          f"({sum(1 for r in inband if r['ok'])} recoveries)")

    # ----------------------------------------------------------------- X2
    print("\n" + "-" * 78)
    print("EXP X2: PRE-REGISTERED test — AUC inside the ambiguous NU band,")
    print("        pooled over eff.  Read with X3: eff is NOT controlled here.")
    print("-" * 78)
    bands = [("NU < %.2f" % BAND_LO, lambda r: r['NU'] < BAND_LO),
             ("%.2f-%.2f  (ambiguous)" % (BAND_LO, BAND_HI),
              lambda r: BAND_LO <= r['NU'] <= BAND_HI),
             ("NU > %.2f" % BAND_HI, lambda r: r['NU'] > BAND_HI)]
    print(f"{'NU band':>24} {'N':>5} {'rec':>8} | {'mu':>7} {'nu_hat':>7} "
          f"{'eff':>7} {'step':>7} {'NU':>7}")
    for name, pred in bands:
        sub = [r for r in rows if pred(r)]
        if not sub:
            continue
        rec = sum(1 for r in sub if r['ok'])
        print(f"{name:>24} {len(sub):>5} {str(rec)+'/'+str(len(sub)):>8} | "
              f"{fmt_auc(sub,'mu')} {fmt_auc(sub,'nuhat')} "
              f"{fmt_auc(sub,'eff')} {fmt_auc(sub,'step')} "
              f"{fmt_auc(sub,'NU')}")
    print("\n  AUC > 0.5 means SMALLER score -> more likely to recover.")

    # ----------------------------------------------------------------- X3
    print("\n" + "-" * 78)
    print("EXP X3: HONEST test — doubly conditioned on (eff stratum, NU band).")
    print("        Inside a cell both the bias strength and the BDD certificate")
    print("        are fixed, so any AUC left is pure cross-curve geometry.")
    print("-" * 78)
    print(f"{'eff':>5} {'NU band':>22} {'N':>5} {'rec':>8} | {'mu':>7} "
          f"{'nu_hat':>7} {'step':>7} {'gsa':>7}")
    cells_amb = []
    for eff in EFFS:
        for name, pred in bands:
            sub = [r for r in rows if r['effq'] == eff and pred(r)]
            if len(sub) < 4:
                continue
            rec = sum(1 for r in sub if r['ok'])
            if name.endswith("(ambiguous)"):
                cells_amb.append(sub)
            print(f"{eff:>5.2f} {name:>22} {len(sub):>5} "
                  f"{str(rec)+'/'+str(len(sub)):>8} | {fmt_auc(sub,'mu')} "
                  f"{fmt_auc(sub,'nuhat')} {fmt_auc(sub,'step')} "
                  f"{fmt_auc(sub,'gsa')}")

    print("\n  stratified pooling over the ambiguous cells "
          "(concordant pairs / all pairs):")
    for key in ('mu', 'nuhat', 'step', 'gsa', 'lamstar'):
        a, t = stratified_auc(cells_amb, key)
        print(f"    AUC(-{key:<8}) = {a:.4f}   ({t} pairs)")
    allcells = [[r for r in rows if r['effq'] == eff] for eff in EFFS]
    print("\n  same pooling over ALL eff cells (no NU conditioning), "
          "for reference:")
    for key in ('mu', 'nuhat', 'NU', 'step', 'gsa'):
        a, t = stratified_auc(allcells, key)
        print(f"    AUC(-{key:<8}) = {a:.4f}   ({t} pairs)")

    # ----------------------------------------------------------------- X4
    print("\n" + "-" * 78)
    print("EXP X4: is NU redundant given mu?  Spearman inside each eff stratum")
    print("        and inside the ambiguous band.")
    print("-" * 78)
    print(f"{'subset':>26} {'N':>5} {'rho(NU,mu)':>11} {'rho(NU,nuhat)':>14} "
          f"{'rho(mu,eff)':>12}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        print(f"{'eff='+format(eff,'.2f'):>26} {len(sub):>5} "
              f"{spearman([r['NU'] for r in sub],[r['mu'] for r in sub]):>11.4f} "
              f"{spearman([r['NU'] for r in sub],[r['nuhat'] for r in sub]):>14.4f} "
              f"{spearman([r['mu'] for r in sub],[r['eff'] for r in sub]):>12.4f}")
    sub = inband
    print(f"{'ambiguous band (pooled)':>26} {len(sub):>5} "
          f"{spearman([r['NU'] for r in sub],[r['mu'] for r in sub]):>11.4f} "
          f"{spearman([r['NU'] for r in sub],[r['nuhat'] for r in sub]):>14.4f} "
          f"{spearman([r['mu'] for r in sub],[r['eff'] for r in sub]):>12.4f}")

    # ----------------------------------------------------------------- X5
    print("\n" + "-" * 78)
    print("EXP X5: the 2-parameter viability test.  Grouped 5-fold CV by CURVE,")
    print("        so every test curve is unseen at fit time.")
    print("-" * 78)
    models = [
        ('NU alone',            ['logNU']),
        ('nu_hat alone',        ['lognuhat']),
        ('eff alone',           ['logeff']),
        ('NU + nu_hat',         ['logNU', 'lognuhat']),
        ('NU + eff',            ['logNU', 'logeff']),
        ('nu_hat + eff',        ['lognuhat', 'logeff']),
        ('NU + nu_hat + eff',   ['logNU', 'lognuhat', 'logeff']),
        ('NU + nu_hat + eff + step',
         ['logNU', 'lognuhat', 'logeff', 'step']),
    ]
    print(f"{'model':>28} {'CV AUC':>8} {'in-sample AUC':>15}")
    for name, keys in models:
        fs = [FEATS[k] for k in keys]
        cv, ntest = grouped_cv_auc(rows, fs)
        beta = logistic_fit([[f(r) for f in fs] for r in rows],
                            [1.0 if r['ok'] else 0.0 for r in rows])
        pr = [logistic_score(beta, [f(r) for f in fs]) for r in rows]
        ins = 1.0 - auc([p for p, r in zip(pr, rows) if r['ok']],
                        [p for p, r in zip(pr, rows) if not r['ok']])
        print(f"{name:>28} {cv:>8.4f} {ins:>15.4f}")

    print("\n  full-data coefficients, decision boundary at logit = 0:")
    for name, keys in (('NU + nu_hat', ['logNU', 'lognuhat']),
                       ('NU + nu_hat + eff',
                        ['logNU', 'lognuhat', 'logeff'])):
        fs = [FEATS[k] for k in keys]
        beta = logistic_fit([[f(r) for f in fs] for r in rows],
                            [1.0 if r['ok'] else 0.0 for r in rows])
        terms = "  ".join(f"{beta[i+1]:+.3f}*{k}" for i, k in enumerate(keys))
        print(f"    {name:<20} logit = {beta[0]:+.3f}  {terms}")

    # ----------------------------------------------------------------- X6
    print("\n" + "-" * 78)
    print("EXP X6: SECONDARY — the two-block step at 17 bits / dim 24.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'mean step':>10} {'sd':>7} {'min':>7} "
          f"{'max':>7} | {'mean gsa':>9} {'head/mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        st = [r['step'] for r in sub]
        mu_ = sum(st) / len(st)
        sd = math.sqrt(sum((x - mu_) ** 2 for x in st) / len(st))
        g = [r['gsa'] for r in sub]
        h = [r['head_over_mu'] for r in sub]
        print(f"{eff:>5.2f} {len(sub):>5} {mu_:>10.4f} {sd:>7.4f} "
              f"{min(st):>7.4f} {max(st):>7.4f} | {sum(g)/len(g):>9.5f} "
              f"{sum(h)/len(h):>8.4f}")
    total_range = [max(r['lgprof']) - min(r['lgprof']) for r in rows]
    print(f"\n  full profile range max(log2||b*||) - min : "
          f"mean {sum(total_range)/len(total_range):.3f} bits, "
          f"max {max(total_range):.3f} bits")
    print("  (W1b saw a ~1.6-bit step at 12 bits / dim 16; compare above)")
    print(f"\n  is step just nu_hat re-measured from the profile?")
    print(f"    rho(step, log nu_hat) = "
          f"{spearman([r['step'] for r in rows],[math.log(r['nuhat']) for r in rows]):+.4f}")
    print(f"    rho(step, log NU)     = "
          f"{spearman([r['step'] for r in rows],[math.log(r['NU']) for r in rows]):+.4f}")
    print(f"    rho(gsa,  log nu_hat) = "
          f"{spearman([r['gsa'] for r in rows],[math.log(r['nuhat']) for r in rows]):+.4f}")

    # ----------------------------------------------------------------- X7
    print("\n" + "-" * 78)
    print("EXP X7: the X5 boundary has near-equal coefficients on logNU and")
    print("        lognu_hat, i.e. it is the product rule NU*nu_hat <= t.")
    print("        Test it directly — no fitting, no free parameter.")
    print("-" * 78)
    for r in rows:
        r['prod'] = r['NU'] * r['nuhat']
    print(f"{'subset':>26} {'N':>5} {'rec':>9} {'AUC(-prod)':>11}")
    print(f"{'pooled':>26} {len(rows):>5} "
          f"{str(sum(1 for r in rows if r['ok']))+'/'+str(len(rows)):>9} "
          f"{fmt_auc(rows,'prod'):>11}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        print(f"{'eff='+format(eff,'.2f'):>26} {len(sub):>5} "
              f"{str(sum(1 for r in sub if r['ok']))+'/'+str(len(sub)):>9} "
              f"{fmt_auc(sub,'prod'):>11}")
    print(f"{'ambiguous NU band':>26} {len(inband):>5} "
          f"{str(sum(1 for r in inband if r['ok']))+'/'+str(len(inband)):>9} "
          f"{fmt_auc(inband,'prod'):>11}")
    a, t = stratified_auc(allcells, 'prod')
    print(f"\n  stratified over eff cells: AUC(-prod) = {a:.4f}  ({t} pairs)")
    a, t = stratified_auc(cells_amb, 'prod')
    print(f"  stratified over ambiguous cells: AUC(-prod) = {a:.4f}  ({t} pairs)")

    print(f"\n  threshold sweep on NU*nu_hat (pooled, N={len(rows)}):")
    print(f"{'t':>7} {'TP':>5} {'FP':>5} {'FN':>5} {'TN':>5} {'prec':>7} "
          f"{'recall':>7} {'acc':>7}")
    for t_ in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.962, 1.0, 1.1, 1.25, 1.5):
        TP = sum(1 for r in rows if r['ok'] and r['prod'] <= t_)
        FP = sum(1 for r in rows if not r['ok'] and r['prod'] <= t_)
        FN = sum(1 for r in rows if r['ok'] and r['prod'] > t_)
        TN = sum(1 for r in rows if not r['ok'] and r['prod'] > t_)
        prec = TP / (TP + FP) if TP + FP else float('nan')
        print(f"{t_:>7.3f} {TP:>5} {FP:>5} {FN:>5} {TN:>5} {prec:>7.3f} "
              f"{TP/(TP+FN):>7.3f} {(TP+TN)/len(rows):>7.3f}")
    print("\n  t = 0.962 is the logit=0 crossing of the fitted NU+nu_hat model.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
