"""
GLV-HNP Phase 2, Thread 25: is there a SECOND mechanism, independent of the
BDD certificate NU, that controls the Kannan-LLL wall?

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

Why the sign conventions matter here.  `auc(pos, neg)` scores -x, so AUC > 0.5
always means "SMALLER score -> more likely to recover".  mu = lambda_1(L2) is
NOT dimensionless: it scales like sqrt(n*S_K1*S_K2), so across NU bands (which
pool eff strata and therefore pool K1) a raw-mu AUC can be reading size rather
than geometry.  Every mu column below is therefore reported next to its
dimensionless twin nu_hat = mu/sqrt(det L2), and X2 re-runs the whole thing
with eff pinned as well so neither can hide.

Run: python3 glv_hnp_phase2_nu_conditional.py [--json PATH] [--big]
     --json  read the table dumped by glv_hnp_phase2_gsprofile_strat.py
             --dump-json instead of regenerating it
     --big   also run the 15-stratum x 10-seed extension (X5)
"""

import json
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_gsprofile import auc, spearman
from glv_hnp_phase2_gsprofile_strat import collect_rows

M17 = 12
EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
SEEDS10 = [42, 1234, 9999, 555, 31337, 7, 20260807, 99, 123456, 8675309]

# The 17-bit nearest-plane bracket measured in Thread 24 W4:
#   NU < 1.040 sufficient (71 TP / 0 FP), NU > 2.199 necessary.
BAND_LO, BAND_HI = 1.040, 2.199


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _midranks(vals):
    """Average ranks (1-based), ties averaged."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            r[order[t]] = avg
        i = j + 1
    return r


def fast_auc(vals, labels):
    """AUC of the score -vals for separating labels==True from labels==False.

    Mann-Whitney via mid-ranks, O(N log N).  Agrees with the O(N^2) `auc` of
    glv_hnp_phase2_gsprofile to within float rounding, ties included; the
    permutation loops below would be intractable with the quadratic form.
    """
    npos = sum(1 for l in labels if l)
    nneg = len(labels) - npos
    if npos == 0 or nneg == 0:
        return float('nan')
    r = _midranks(vals)
    rsum = sum(ri for ri, l in zip(r, labels) if l)
    # U for "pos has SMALLER score", which is what AUC(-x) means here.
    u = npos * nneg + npos * (npos + 1) / 2.0 - rsum
    return u / (npos * nneg)


def auc_perm(pos, neg, score, nperm=20000, seed=11):
    """Two-sided permutation p-value for AUC(-score) against 0.5."""
    obs = auc(pos, neg)
    if obs != obs:
        return obs, float('nan')
    allv = list(pos) + list(neg)
    labels = [True] * len(pos) + [False] * len(neg)
    # Ranks are permutation-invariant, so rank once and shuffle the ranks.
    ranks = _midranks(allv)
    npos, nneg = len(pos), len(neg)
    denom = npos * nneg
    base = npos * nneg + npos * (npos + 1) / 2.0
    rng = random.Random(seed)
    dev = abs(obs - 0.5)
    hit = 0
    for _ in range(nperm):
        rng.shuffle(ranks)
        a = (base - sum(ranks[:npos])) / denom
        if abs(a - 0.5) >= dev - 1e-12:
            hit += 1
    return obs, (hit + 1) / (nperm + 1)


def partial_spearman(x, y, z):
    """Spearman(x, y) with z partialled out, via rank residuals."""
    rxy, rxz, ryz = spearman(x, y), spearman(x, z), spearman(y, z)
    den = math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    return (rxy - rxz * ryz) / den if den > 1e-12 else float('nan')


def _solve(A, rhs):
    """Gaussian elimination with partial pivoting; A is small (<= 4x4)."""
    k = len(rhs)
    M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for c in range(k):
        piv = max(range(c, k), key=lambda i: abs(M[i][c]))
        if abs(M[piv][c]) < 1e-14:
            return None
        M[c], M[piv] = M[piv], M[c]
        for i in range(k):
            if i == c:
                continue
            f = M[i][c] / M[c][c]
            for j in range(c, k + 1):
                M[i][j] -= f * M[c][j]
    return [M[i][k] / M[i][i] for i in range(k)]


def logistic_fit(X, y, l2=1e-3, iters=40, tol=1e-10):
    """Ridge-regularised logistic regression by IRLS (Newton) on standardised
    features.  Converges in ~8 steps at these sizes; the earlier plain-GD
    version needed 4000 sweeps and did not fully converge on the 3-feature fit.

    Returns (w_std, b_std, mean, sd, w_raw, b_raw); the raw pair lets the
    decision boundary be quoted in the original units.
    """
    nfeat = len(X[0])
    N = len(X)
    mean = [sum(r[j] for r in X) / N for j in range(nfeat)]
    sd = [math.sqrt(sum((r[j] - mean[j]) ** 2 for r in X) / N) or 1.0
          for j in range(nfeat)]
    # design matrix with the intercept as column 0
    Z = [[1.0] + [(r[j] - mean[j]) / sd[j] for j in range(nfeat)] for r in X]
    k = nfeat + 1
    beta = [0.0] * k
    for _ in range(iters):
        H = [[0.0] * k for _ in range(k)]
        g = [0.0] * k
        for zi, yi in zip(Z, y):
            t = sum(beta[j] * zi[j] for j in range(k))
            t = max(-40.0, min(40.0, t))
            pr = 1.0 / (1.0 + math.exp(-t))
            wgt = max(pr * (1.0 - pr), 1e-9)
            d = yi - pr
            for a in range(k):
                g[a] += d * zi[a]
                za = zi[a] * wgt
                for bb in range(a, k):
                    H[a][bb] += za * zi[bb]
        for a in range(k):
            for bb in range(a):
                H[a][bb] = H[bb][a]
        for a in range(1, k):           # ridge on the slopes only
            H[a][a] += l2 * N
            g[a] -= l2 * N * beta[a]
        step = _solve(H, g)
        if step is None:
            break
        beta = [beta[j] + step[j] for j in range(k)]
        if max(abs(s) for s in step) < tol:
            break
    b, w = beta[0], beta[1:]
    w_raw = [w[j] / sd[j] for j in range(nfeat)]
    b_raw = b - sum(w[j] * mean[j] / sd[j] for j in range(nfeat))
    return w, b, mean, sd, w_raw, b_raw


def logloss(X, y, w_raw, b_raw):
    tot = 0.0
    for xi, yi in zip(X, y):
        t = b_raw + sum(w_raw[j] * xi[j] for j in range(len(xi)))
        pr = 1.0 / (1.0 + math.exp(-max(-40.0, min(40.0, t))))
        pr = min(max(pr, 1e-12), 1 - 1e-12)
        tot += -(yi * math.log(pr) + (1 - yi) * math.log(1 - pr))
    return tot / len(X)


def grouped_cv_auc(rows, feats, groupkey, nfold=5, seed=3):
    """Leave-groups-out CV.  Groups are CURVES, so a fit cannot memorise a
    curve's geometry and then be scored on the same curve."""
    groups = sorted({r[groupkey] for r in rows})
    rng = random.Random(seed)
    rng.shuffle(groups)
    folds = [groups[i::nfold] for i in range(nfold)]
    scores, labels = [], []
    for f in folds:
        held = set(f)
        tr = [r for r in rows if r[groupkey] not in held]
        te = [r for r in rows if r[groupkey] in held]
        if not te or len({r['ok'] for r in tr}) < 2:
            continue
        Xtr = [[fn(r) for fn in feats] for r in tr]
        ytr = [1.0 if r['ok'] else 0.0 for r in tr]
        _, _, _, _, wr, br = logistic_fit(Xtr, ytr)
        for r in te:
            t = br + sum(wr[j] * feats[j](r) for j in range(len(feats)))
            scores.append(-t)          # auc() scores -x, so negate the logit
            labels.append(bool(r['ok']))
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    return auc(pos, neg), len(pos), len(neg)


def step_stat(r, m=M17):
    """log2(||b*_{m+1}||) - log2(||b*_1||), the W1b two-block step."""
    p = r['prof']
    if len(p) <= m or p[0] <= 0 or p[m] <= 0:
        return float('nan')
    return math.log2(p[m]) - math.log2(p[0])


def head_flatness(r, m=M17):
    """max-min of log2||b*_i|| over the head block i < m."""
    p = [x for x in r['prof'][:m] if x > 0]
    return math.log2(max(p)) - math.log2(min(p)) if p else float('nan')


def col(rows, key):
    return ([key(r) for r in rows if r['ok']],
            [key(r) for r in rows if not r['ok']])


def aucs(rows, keys):
    out = {}
    for name, fn in keys.items():
        pos, neg = col(rows, fn)
        out[name] = auc(pos, neg)
    return out


SCORES = {
    'mu':      lambda r: r['mu'],
    'nu_hat':  lambda r: r['nuhat'],
    'NU':      lambda r: r['NU'],
    'eff':     lambda r: r['eff'],
    'lam*':    lambda r: r['lamstar'],
    'step':    step_stat,
    'n':       lambda r: r['n'],
    'PROD':    lambda r: r['NU'] * r['nuhat'],
}


# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("Thread 25 — conditioning on NU: is mu a second coordinate?")
    print("=" * 78)

    if "--json" in sys.argv:
        path = sys.argv[sys.argv.index("--json") + 1]
        rows = json.load(open(path))
        print(f"\nloaded {len(rows)} rows from {path}")
    else:
        t0 = time.time()
        curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
        rows = collect_rows(curves17, M17, EFFS)
        print(f"\n{len(rows)} instances (float GS, dim {rows[0]['k']}) "
              f"in {time.time()-t0:.1f}s")

    nrec = sum(1 for r in rows if r['ok'])
    print(f"recovery {nrec}/{len(rows)}   "
          f"NU range [{min(r['NU'] for r in rows):.3f}, "
          f"{max(r['NU'] for r in rows):.3f}]")

    # fast_auc replaces the O(N^2) auc() inside the permutation loops only;
    # check the two agree before relying on it.
    worst = 0.0
    for name, fn in SCORES.items():
        pos, neg = col(rows, fn)
        v = [fn(r) for r in rows]
        lab = [bool(r['ok']) for r in rows]
        worst = max(worst, abs(auc(pos, neg) - fast_auc(v, lab)))
    print(f"fast_auc vs auc: max |difference| over all score columns = "
          f"{worst:.2e}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1 (H25): AUC inside NU bands.  NU is ~constant within a band,")
    print("        so anything that still separates is a second mechanism.")
    print("-" * 78)
    bands = [(0.0, 1.040, "NU < 1.040  (sufficient)"),
             (1.040, 1.500, "1.040-1.500 (ambiguous, lower)"),
             (1.500, 2.199, "1.500-2.199 (ambiguous, upper)"),
             (1.040, 2.199, "1.040-2.199 (AMBIGUOUS BAND = H25)"),
             (2.199, 99.0, "NU > 2.199  (necessary violated)")]
    hdr = ('band', 'N', 'rec', 'mu', 'nu_hat', 'NU', 'eff', 'lam*', 'step',
           'n')
    print(f"{hdr[0]:<32}{hdr[1]:>5}{hdr[2]:>8} |" +
          "".join(f"{h:>9}" for h in hdr[3:]))
    band_rows = {}
    for lo, hi, name in bands:
        sub = [r for r in rows if lo <= r['NU'] <= hi] if lo == 0 else \
              [r for r in rows if lo < r['NU'] <= hi]
        band_rows[name] = sub
        w = sum(1 for r in sub if r['ok'])
        if not sub or w == 0 or w == len(sub):
            print(f"{name:<32}{len(sub):>5}{str(w)+'/'+str(len(sub)):>8} |"
                  f"{'(degenerate)':>27}")
            continue
        a = aucs(sub, SCORES)
        print(f"{name:<32}{len(sub):>5}{str(w)+'/'+str(len(sub)):>8} |" +
              "".join(f"{a[k]:>9.4f}" for k in hdr[3:]))

    key = band_rows["1.040-2.199 (AMBIGUOUS BAND = H25)"]
    pos, neg = col(key, SCORES['mu'])
    a_mu, p_mu = auc_perm(pos, neg, 'mu')
    pos, neg = col(key, SCORES['nu_hat'])
    a_nh, p_nh = auc_perm(pos, neg, 'nu_hat')
    print(f"\nH25 threshold is AUC(-mu) >= 0.80 in the ambiguous band:")
    print(f"  AUC(-mu)     = {a_mu:.4f}   permutation p = {p_mu:.2e}   "
          f"-> H25 {'HOLDS' if a_mu >= 0.80 else 'FALSIFIED'}")
    print(f"  AUC(-nu_hat) = {a_nh:.4f}   permutation p = {p_nh:.2e}")
    print(f"  N = {len(key)}, recoveries = {sum(1 for r in key if r['ok'])}")

    print("\nis mu's power inside the band mediated by NU?  (rank statistics)")
    x = [r['mu'] for r in key]
    yb = [1.0 if r['ok'] else 0.0 for r in key]
    z = [r['NU'] for r in key]
    print(f"  Spearman(mu, NU)          = {spearman(x, z):+.4f}")
    print(f"  Spearman(mu, recovery)    = {spearman(x, yb):+.4f}")
    print(f"  partial(mu, recovery | NU)= {partial_spearman(x, yb, z):+.4f}")
    xe = [r['eff'] for r in key]
    print(f"  Spearman(mu, eff)         = {spearman(x, xe):+.4f}")
    print(f"  partial(mu, recovery |eff)= {partial_spearman(x, yb, xe):+.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: double conditioning — NU band AND eff stratum at once.")
    print("        Neither the bias strength nor the BDD certificate can move")
    print("        inside a cell, so this is the strongest available control.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>8} | {'AUC mu':>8} {'AUC nu_hat':>11} "
          f"{'AUC NU':>8} {'AUC step':>9}")
    tot_pairs = {'mu': [0.0, 0], 'nu_hat': [0.0, 0], 'NU': [0.0, 0],
                 'step': [0.0, 0], 'PROD': [0.0, 0], 'n': [0.0, 0]}
    for eff in EFFS:
        sub = [r for r in key if r['effq'] == eff]
        w = sum(1 for r in sub if r['ok'])
        if not sub or w == 0 or w == len(sub):
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(w)+'/'+str(len(sub)):>8} | {'(degenerate)':>8}")
            continue
        a = aucs(sub, SCORES)
        npair = w * (len(sub) - w)
        for k in tot_pairs:
            tot_pairs[k][0] += a[k] * npair
            tot_pairs[k][1] += npair
        print(f"{eff:>5.2f} {len(sub):>5} {str(w)+'/'+str(len(sub)):>8} | "
              f"{a['mu']:>8.4f} {a['nu_hat']:>11.4f} {a['NU']:>8.4f} "
              f"{a['step']:>9.4f}")
    print("\npair-weighted mean AUC over the non-degenerate cells:")
    for k in ('mu', 'nu_hat', 'PROD', 'NU', 'step', 'n'):
        s, npair = tot_pairs[k]
        print(f"  {k:>7}: {s/npair:.4f}" if npair else f"  {k:>7}:   n/a")
    print("  ('n' is the size control: within a cell K1/K2 are pinned to n,")
    print("   so a large AUC here would mean the geometry columns are proxies")
    print("   for curve size rather than for lattice shape.)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: the 2-parameter viability test — logistic on (log NU, log mu)")
    print("-" * 78)
    y = [1.0 if r['ok'] else 0.0 for r in rows]
    variants = [
        ("NU only",            [lambda r: math.log(r['NU'])]),
        ("mu only",            [lambda r: math.log(r['mu'])]),
        ("nu_hat only",        [lambda r: math.log(r['nuhat'])]),
        ("(log NU, log mu)",   [lambda r: math.log(r['NU']),
                                lambda r: math.log(r['mu'])]),
        ("(log NU, log nu_hat)", [lambda r: math.log(r['NU']),
                                  lambda r: math.log(r['nuhat'])]),
        ("(log NU, log nu_hat, log eff)",
                               [lambda r: math.log(r['NU']),
                                lambda r: math.log(r['nuhat']),
                                lambda r: math.log(r['eff'])]),
    ]
    print(f"{'features':<32}{'in-sample AUC':>15}{'logloss':>10}"
          f"{'grouped-CV AUC':>16}")
    fits = {}
    insample = {}
    for name, feats in variants:
        X = [[fn(r) for fn in feats] for r in rows]
        _, _, _, _, wr, br = logistic_fit(X, y)
        sc = [-(br + sum(wr[j] * X[i][j] for j in range(len(feats))))
              for i in range(len(rows))]
        a_in = auc([s for s, r in zip(sc, rows) if r['ok']],
                   [s for s, r in zip(sc, rows) if not r['ok']])
        ll = logloss(X, y, wr, br)
        a_cv, _, _ = grouped_cv_auc(rows, feats, 'n')
        fits[name] = (wr, br)
        insample[name] = a_in
        print(f"{name:<32}{a_in:>15.4f}{ll:>10.4f}{a_cv:>16.4f}")

    wr, br = fits["(log NU, log mu)"]
    print(f"\ndecision boundary, all {len(rows)} rows, P(recover) = 1/2 at")
    print(f"  {br:+.4f} {wr[0]:+.4f}*log(NU) {wr[1]:+.4f}*log(mu) = 0")
    if abs(wr[1]) > 1e-9:
        print(f"  i.e.  NU^{wr[0]:.3f} * mu^{wr[1]:.3f} = "
              f"{math.exp(-br):.4g}")
    wr2, br2 = fits["(log NU, log nu_hat)"]
    print(f"dimensionless version:")
    print(f"  {br2:+.4f} {wr2[0]:+.4f}*log(NU) {wr2[1]:+.4f}*log(nu_hat) = 0")
    if abs(wr2[1]) > 1e-9:
        print(f"  i.e.  NU^{-wr2[0]:.3f} * nu_hat^{-wr2[1]:.3f} = "
              f"{math.exp(br2):.4g}   (recover BELOW this)")
        print(f"  exponent ratio wNU/wnuhat = {wr2[0]/wr2[1]:.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3b: the fitted exponents are nearly equal — does the plain")
    print("         PRODUCT NU * nu_hat, with no free parameter, do as well?")
    print("-" * 78)
    print(f"{'alpha (score = NU * nu_hat^alpha)':<36}{'AUC':>10}")
    best = (None, -1.0)
    for alpha in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0):
        pos, neg = col(rows, lambda r: r['NU'] * r['nuhat'] ** alpha)
        a = auc(pos, neg)
        if a > best[1]:
            best = (alpha, a)
        mark = "  <- the plain product" if alpha == 1.0 else ""
        print(f"{('alpha = %.2f' % alpha):<36}{a:>10.4f}{mark}")
    print(f"\nbest alpha on this grid = {best[0]:.2f} (AUC {best[1]:.4f}); "
          f"the free 2-param logistic reached "
          f"{insample['(log NU, log nu_hat)']:.4f} in-sample.")
    pos, neg = col(rows, SCORES['PROD'])
    a_pr, p_pr = auc_perm(pos, neg, 'PROD')
    print(f"AUC(-NU*nu_hat) = {a_pr:.4f}  (perm p = {p_pr:.2e}), "
          f"zero fitted parameters.")
    print("per-eff-stratum, so the product cannot score by tracking eff:")
    print(f"{'eff':>5} {'N':>5} {'rec':>8} | {'AUC PROD':>9} {'AUC NU':>8} "
          f"{'AUC nu_hat':>11}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        w = sum(1 for r in sub if r['ok'])
        if not sub or w == 0 or w == len(sub):
            print(f"{eff:>5.2f} {len(sub):>5} {str(w)+'/'+str(len(sub)):>8} | "
                  f"{'(degenerate)':>9}")
            continue
        a = aucs(sub, SCORES)
        print(f"{eff:>5.2f} {len(sub):>5} {str(w)+'/'+str(len(sub)):>8} | "
              f"{a['PROD']:>9.4f} {a['NU']:>8.4f} {a['nu_hat']:>11.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X4 (secondary): the two-block step of W1b.")
    print("        step = log2||b*_{m+1}|| - log2||b*_1||;  W1b says step -> 0")
    print("        exactly as the K1 wall is crossed.")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>8} | {'mean step':>10} "
          f"{'step|rec':>9} {'step|fail':>10} | {'AUC step':>9} "
          f"{'AUC head-flat':>14}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        st = [step_stat(r) for r in sub]
        sr = [step_stat(r) for r in sub if r['ok']]
        sf = [step_stat(r) for r in sub if not r['ok']]
        w = len(sr)
        a_st = auc(sr, sf) if sr and sf else float('nan')
        pos, neg = col(sub, head_flatness)
        a_hf = auc(pos, neg) if pos and neg else float('nan')
        print(f"{eff:>5.2f} {len(sub):>5} {str(w)+'/'+str(len(sub)):>8} | "
              f"{sum(st)/len(st):>10.4f} "
              f"{(sum(sr)/len(sr) if sr else float('nan')):>9.4f} "
              f"{(sum(sf)/len(sf) if sf else float('nan')):>10.4f} | "
              f"{a_st:>9.4f} {a_hf:>14.4f}")
    pos, neg = col(rows, step_stat)
    a_st_all, p_st = auc_perm(pos, neg, 'step')
    print(f"\npooled AUC(-step) = {a_st_all:.4f}  (perm p = {p_st:.2e})")
    print(f"pooled AUC(-NU)   = {aucs(rows, SCORES)['NU']:.4f}")
    print(f"Spearman(step, log eff)  = "
          f"{spearman([step_stat(r) for r in rows], [math.log(r['eff']) for r in rows]):+.4f}")
    print(f"Spearman(step, log NU)   = "
          f"{spearman([step_stat(r) for r in rows], [math.log(r['NU']) for r in rows]):+.4f}")
    print(f"Spearman(step, log mu)   = "
          f"{spearman([step_stat(r) for r in rows], [math.log(r['mu']) for r in rows]):+.4f}")
    print("\nhead-block check of W1b: max_i<m |log2||b*_i|| - log2 mu|")
    devs = []
    for r in rows:
        lm = math.log2(r['mu'])
        devs.append(max(abs(math.log2(p) - lm) for p in r['prof'][:M17]
                        if p > 0))
    devs.sort()
    print(f"  median {devs[len(devs)//2]:.4f} bits   "
          f"p90 {devs[int(0.9*len(devs))]:.4f}   max {devs[-1]:.4f}")

    # -----------------------------------------------------------------------
    big17 = None
    if "--big" in sys.argv:
        print("\n" + "-" * 78)
        print("EXP X5: power extension — 15 eff strata x 20 curves x 10 seeds.")
        print("-" * 78)
        t0 = time.time()
        curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
        effs = tuple(round(0.03 + 0.02 * i, 2) for i in range(15))
        big = collect_rows(curves17, M17, effs, seeds=SEEDS10)
        big17 = big
        print(f"{len(big)} instances in {time.time()-t0:.1f}s   "
              f"recovery {sum(1 for r in big if r['ok'])}/{len(big)}")
        kb = [r for r in big if BAND_LO < r['NU'] <= BAND_HI]
        w = sum(1 for r in kb if r['ok'])
        a = aucs(kb, SCORES)
        print(f"\nambiguous band: N={len(kb)}  rec={w}/{len(kb)}")
        for k in ('mu', 'nu_hat', 'NU', 'eff', 'lam*', 'step'):
            print(f"  AUC(-{k:<7}) = {a[k]:.4f}")
        pos, neg = col(kb, SCORES['mu'])
        a_mu2, p_mu2 = auc_perm(pos, neg, 'mu')
        print(f"  permutation p for mu = {p_mu2:.2e}")
        print("\n  double conditioning (band x eff), pair-weighted mean:")
        tp = {'mu': [0.0, 0], 'nu_hat': [0.0, 0], 'NU': [0.0, 0],
              'step': [0.0, 0]}
        for eff in effs:
            sub = [r for r in kb if r['effq'] == eff]
            w = sum(1 for r in sub if r['ok'])
            if not sub or w == 0 or w == len(sub):
                continue
            a = aucs(sub, SCORES)
            npair = w * (len(sub) - w)
            for k in tp:
                tp[k][0] += a[k] * npair
                tp[k][1] += npair
        for k in ('mu', 'nu_hat', 'NU', 'step'):
            s, npair = tp[k]
            print(f"    {k:>7}: {s/npair:.4f}  ({npair} pairs)" if npair
                  else f"    {k:>7}: n/a")
        ycv = [1.0 if r['ok'] else 0.0 for r in big]
        for name, feats in variants:
            a_cv, _, _ = grouped_cv_auc(big, feats, 'n')
            print(f"  grouped-CV AUC {name:<32} {a_cv:.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X6: cross-size stability.  Thread 24 W4 found NU is a")
    print("        size-DEGRADING separator (AUC 0.978 at 12 bits -> 0.860 at")
    print("        17).  Does the product NU*nu_hat hold its threshold?")
    print("-" * 78)
    from glv_hnp_common import find_generator
    from glv_hnp_phase2_projected import HIST
    from glv_hnp_phase2_gsprofile import instance
    from glv_hnp_phase2_projected import run_new as _run_new
    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    U2 = [("12-bit/2557", HIST[1], 8), ("12-bit/2677", HIST[2], 10)]
    r12 = []
    for label, h, m in U2:
        _lab, p, b, n, lam, _k1, _m = h
        G = find_generator(p, b, n)
        curve = (p, b, n, lam, G)
        for k1 in K1_GRID:
            for seed in [42, 1234, 9999, 555, 31337]:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                rr = instance(curve, m, d_trial, k1, seed, exact=False)
                if rr is None:
                    continue
                rk = _run_new(curve, m, d_trial, k1, seed)
                rr.update({'label': label, 'n': n, 'K1': k1, 'seed': seed,
                           'ok': bool(rk['ok']),
                           'eff': k1 * (math.isqrt(n) + 1) / n,
                           'lamstar': 0.0, 'effq': k1})
                r12.append(rr)
    print(f"12-bit grid: {len(r12)} instances, "
          f"recovery {sum(1 for r in r12 if r['ok'])}/{len(r12)}")

    def bracket(rs, fn):
        p = [fn(r) for r in rs if r['ok']]
        q = [fn(r) for r in rs if not r['ok']]
        if not p or not q:
            return None
        return max(p), min(q), auc(p, q)

    print(f"\n{'set':<14}{'stat':<12}{'AUC':>8}{'max|rec':>10}"
          f"{'min|fail':>10}{'overlap':>9}")
    for setname, rs in (("12-bit", r12), ("17-bit", rows),
                        ("17-bit big", big17 if big17 else rows)):
        if rs is None:
            continue
        for sname, fn in (("NU", SCORES['NU']), ("NU*nu_hat", SCORES['PROD'])):
            br = bracket(rs, fn)
            if br is None:
                continue
            hi_rec, lo_fail, a = br
            print(f"{setname:<14}{sname:<12}{a:>8.4f}{hi_rec:>10.4f}"
                  f"{lo_fail:>10.4f}{hi_rec/lo_fail:>9.2f}x")
    print("\n'overlap' = max score among recoveries / min score among failures.")
    print("1.00x would be a perfectly separating threshold; the smaller the")
    print("ratio, the tighter the sufficient/necessary bracket.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
