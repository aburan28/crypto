"""
GLV-HNP Phase 2, Thread 25: is mu a SECOND coordinate, or is its apparent
power mediated by NU after all?

Pre-registered by the 2026-08-07 run #2 (Thread 24) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the nearest-plane
       certificate gives no answer), AUC(-mu -> Kannan-LLL recovery) stays
       >= 0.8.
  If yes  -> mu is a genuine second coordinate and (NU, mu) is the
             2-parameter viability test Phase 2 has been missing; deliverable
             is a logistic fit on (log NU, log mu) with its decision boundary.
  If no   -> mu's power is mediated by NU, W5 was a stratification artifact,
             and the closed form should be retired.

Secondary (also pre-registered):
  define  step = log2||b*_{m+1}|| - log2||b*_1||   (the two-block gap of W1b)
  and test whether step -> 0 predicts the wall better than NU or mu.

Background.  Thread 23 established NU = max_i 2|<e,b*_i>|/||b*_i||^2 as an
exact, unfitted sufficient certificate for Babai nearest-plane (NU <= 1 =>
recovery, 96 TP / 0 FP over 410 instances).  Thread 24 W5/W6 found that
nu_hat = mu/sqrt(det L2) predicts *Kannan-LLL* recovery better than NU does at
fixed bias strength (AUC 0.75-0.93 vs 0.35-0.73) while being UNCORRELATED with
NU inside every eff stratum (Spearman -0.28..+0.16).  Two uncorrelated
predictors of the same event means two mechanisms; H25 asks whether they
combine.

Note on mu vs nu_hat.  Thread 24 W5 held eff and n nearly fixed inside a
stratum, where mu and nu_hat = mu/sqrt(det L2) are the same statistic up to a
constant (AUCs agreed to 3 decimals).  An NU band is NOT a fixed-size slice --
it mixes curves and K1 -- so the sqrt(det) normalisation now matters.  H25 is
stated on mu, so mu is reported as the primary; nu_hat is reported alongside
as the size-free version, and X1 says explicitly which one the verdict rests
on.

Numerics: float Gram-Schmidt, justified by W0/W4 of glv_hnp_phase2_gsprofile.py
(max relative NU error vs exact Fractions ~1e-15 at dim 20 and dim 24).

Run: python3 glv_hnp_phase2_nu_mu.py
"""

import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import lam_star, search_curves
from glv_hnp_phase2_projected import run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# Thread 23/24 bracket at 17 bits, dim 24 (log line ~6300):
#   sufficient  NU < 1.040     necessary  NU > 2.199
BAND_LO, BAND_HI = 1.040, 2.199

SEEDS10 = [42, 1234, 9999, 555, 31337, 8675309, 271828, 161803, 5772, 1618]
EFFS = (0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25)


# ---------------------------------------------------------------------------
# table construction
# ---------------------------------------------------------------------------

def build_table(curves, m, effs, seeds):
    rows = []
    for eff in effs:
        for (p, b, n, lam, G) in curves:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in seeds:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), m, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), m, d_trial, k1b, seed)
                prof = r['prof']
                step = (math.log2(prof[m]) - math.log2(prof[0])
                        if prof[0] > 0 and prof[m] > 0 else float('nan'))
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff, 'seed': seed,
                          'lamstar': lam_star(lam, n), 'step': step, 'm': m})
                rows.append(r)
    return rows


def split_auc(rows, key, sign=-1.0):
    """AUC for the score `sign*key`, house convention (auc() scores -x, so
    AUC > 0.5 means the argument is SMALL on recoveries).

    sign=-1 (default): AUC > 0.5 <=> SMALLER key -> more likely to recover.
    sign=+1          : AUC > 0.5 <=> LARGER  key -> more likely to recover.
    """
    s = -1.0 if sign > 0 else 1.0
    return auc([s * r[key] for r in rows if r['ok']],
               [s * r[key] for r in rows if not r['ok']])


# ---------------------------------------------------------------------------
# logistic regression (IRLS, pure python; no numpy in this container)
# ---------------------------------------------------------------------------

def solve(A, b):
    k = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(k):
        piv = max(range(c, k), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-300:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        for r in range(k):
            if r == c:
                continue
            f = M[r][c] / pv
            if f:
                for j in range(c, k + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][k] / M[i][i] for i in range(k)]


def logistic_fit(X, y, l2=1e-2, iters=60):
    d = len(X[0])
    w = [0.0] * d
    for _ in range(iters):
        g = [0.0] * d
        H = [[l2 if i == j else 0.0 for j in range(d)] for i in range(d)]
        for xi, yi in zip(X, y):
            z = max(-30.0, min(30.0, sum(w[k] * xi[k] for k in range(d))))
            pi = 1.0 / (1.0 + math.exp(-z))
            r = pi - yi
            wt = max(pi * (1.0 - pi), 1e-9)
            for a in range(d):
                g[a] += r * xi[a]
                for c in range(d):
                    H[a][c] += wt * xi[a] * xi[c]
        for a in range(d):
            g[a] += l2 * w[a]
        delta = solve(H, g)
        if delta is None:
            break
        w = [w[k] - delta[k] for k in range(d)]
        if max(abs(x) for x in delta) < 1e-10:
            break
    return w


def predict(w, X):
    out = []
    for xi in X:
        z = max(-30.0, min(30.0, sum(w[k] * xi[k] for k in range(len(w)))))
        out.append(1.0 / (1.0 + math.exp(-z)))
    return out


def auc_scores(scores, labels):
    """AUC with the convention: LARGER score -> label 1."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float('nan')
    conc = sum((1.0 if a > b else 0.5 if a == b else 0.0)
               for a in pos for b in neg)
    return conc / (len(pos) * len(neg))


def acc_at_half(scores, labels):
    return sum(1 for s, l in zip(scores, labels)
               if (s >= 0.5) == bool(l)) / len(labels)


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 78)
    print("Thread 25 — does mu separate INSIDE an NU band?  (H25)")
    print("=" * 78)

    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    t0 = time.time()
    rows = build_table(curves17, 12, EFFS, SEEDS10)
    dt = time.time() - t0
    print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}, "
          f"{len(EFFS)} eff x {len(curves17)} curves x {len(SEEDS10)} seeds) "
          f"in {dt:.1f}s")
    nrec = sum(1 for r in rows if r['ok'])
    print(f"recovery {nrec}/{len(rows)} = {nrec/len(rows):.3f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X0: the NU bracket on this table — does Thread 23's certificate")
    print("        still hold, and how many instances land in the band?")
    print("-" * 78)
    succ = [r['NU'] for r in rows if r['ok']]
    fail = [r['NU'] for r in rows if not r['ok']]
    succ.sort()
    fail.sort()
    print(f"NU | success : min {succ[0]:.3f}  median {succ[len(succ)//2]:.3f}"
          f"  max {succ[-1]:.3f}   (n={len(succ)})")
    print(f"NU | failure : min {fail[0]:.3f}  median {fail[len(fail)//2]:.3f}"
          f"  max {fail[-1]:.3f}   (n={len(fail)})")
    tp = sum(1 for r in rows if r['NU'] <= 1.0 and r['ok'])
    fp = sum(1 for r in rows if r['NU'] <= 1.0 and not r['ok'])
    print(f"NU <= 1 : TP {tp}  FP {fp}   "
          f"(nearest-plane guarantee: FP must be 0)")
    print(f"observed bracket here: sufficient NU < {fail[0]:.3f} , "
          f"necessary NU > {succ[-1]:.3f}")
    print(f"pre-registered band [{BAND_LO}, {BAND_HI}]")

    print(f"\n{'NU band':>16} {'N':>6} {'rec':>10} {'rate':>7}")
    EDGES = [0.0, 0.7, 1.04, 1.4, 1.8, 2.199, 3.0, 1e9]
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        sub = [r for r in rows if lo <= r['NU'] < hi]
        if not sub:
            continue
        w = sum(1 for r in sub if r['ok'])
        hs = "inf" if hi > 1e8 else f"{hi:.2f}"
        print(f"{f'[{lo:.2f},{hs})':>16} {len(sub):>6} "
              f"{str(w)+'/'+str(len(sub)):>10} {w/len(sub):>7.3f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1: H25 — AUC of mu (and nu_hat) INSIDE the ambiguous NU band.")
    print("        NU is uninformative here by construction, so any AUC > 0.5")
    print("        is signal NU does not carry.")
    print("-" * 78)
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    bw = sum(1 for r in band if r['ok'])
    print(f"band [{BAND_LO}, {BAND_HI}] : N={len(band)}  recovery "
          f"{bw}/{len(band)} = {bw/len(band):.3f}\n")
    print(f"{'score':>22} {'AUC in band':>12}")
    for key, label in (('mu', 'mu  (H25 primary)'),
                       ('nuhat', 'nu_hat = mu/sqrt(det)'),
                       ('NU', 'NU  (control)'),
                       ('lamstar', 'lam*  (control)'),
                       ('eff', 'eff')):
        print(f"{label:>22} {split_auc(band, key):>12.4f}")
    a_ns = auc([r['nuhat'] * math.sqrt(r['eff']) for r in band if r['ok']],
               [r['nuhat'] * math.sqrt(r['eff']) for r in band if not r['ok']])
    print(f"{'nu_hat*sqrt(eff)':>22} {a_ns:>12.4f}")
    print(f"{'step (larger=better)':>22} "
          f"{split_auc(band, 'step', sign=+1.0):>12.4f}")

    print(f"\nH25 threshold is AUC(mu) >= 0.8 inside the band.")
    a_mu = split_auc(band, 'mu')
    a_nh = split_auc(band, 'nuhat')
    print(f"  AUC(mu)     = {a_mu:.4f}  -> "
          f"{'HOLDS' if a_mu >= 0.8 else 'FAILS'}")
    print(f"  AUC(nu_hat) = {a_nh:.4f}  -> "
          f"{'HOLDS' if a_nh >= 0.8 else 'FAILS'}")

    # sub-bands, to check the band result is not driven by residual NU trend
    print("\nsub-bands (does NU still trend inside the band?):")
    print(f"{'sub-band':>16} {'N':>6} {'rate':>7} {'AUC mu':>8} "
          f"{'AUC nu_hat':>11} {'AUC NU':>8}")
    SUB = [(BAND_LO, 1.4), (1.4, 1.8), (1.8, BAND_HI)]
    for lo, hi in SUB:
        sub = [r for r in rows if lo <= r['NU'] < hi]
        w = sum(1 for r in sub if r['ok'])
        if not sub or w == 0 or w == len(sub):
            print(f"{f'[{lo:.2f},{hi:.2f})':>16} {len(sub):>6} "
                  f"{'(degenerate)':>7}")
            continue
        print(f"{f'[{lo:.2f},{hi:.2f})':>16} {len(sub):>6} {w/len(sub):>7.3f} "
              f"{split_auc(sub, 'mu'):>8.4f} {split_auc(sub, 'nuhat'):>11.4f} "
              f"{split_auc(sub, 'NU'):>8.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: 2-parameter logistic fit.  Cross-curve split: 10 curves")
    print("        train / 10 curves test, so the boundary must generalise to")
    print("        curves it has never seen.")
    print("-" * 78)
    ns = sorted({r['n'] for r in rows})
    train_n = set(ns[0::2])
    tr = [r for r in rows if r['n'] in train_n]
    te = [r for r in rows if r['n'] not in train_n]
    print(f"train {len(tr)} inst / {len(train_n)} curves    "
          f"test {len(te)} inst / {len(ns)-len(train_n)} curves")

    FEATS = {
        'NU only':          lambda r: [1.0, math.log(r['NU'])],
        'nu_hat only':      lambda r: [1.0, math.log(r['nuhat'])],
        'mu only':          lambda r: [1.0, math.log(r['mu'])],
        'NU + nu_hat':      lambda r: [1.0, math.log(r['NU']),
                                       math.log(r['nuhat'])],
        'NU + mu':          lambda r: [1.0, math.log(r['NU']),
                                       math.log(r['mu'])],
        'NU + nu_hat + eff': lambda r: [1.0, math.log(r['NU']),
                                        math.log(r['nuhat']),
                                        math.log(r['eff'])],
    }
    print(f"\n{'model':>20} {'test AUC':>9} {'test acc':>9} {'band AUC':>9} "
          f"  weights")
    band_te = [r for r in te if BAND_LO <= r['NU'] <= BAND_HI]
    for name, f in FEATS.items():
        Xtr = [f(r) for r in tr]
        ytr = [1.0 if r['ok'] else 0.0 for r in tr]
        w = logistic_fit(Xtr, ytr)
        Xte = [f(r) for r in te]
        yte = [1.0 if r['ok'] else 0.0 for r in te]
        pte = predict(w, Xte)
        a = auc_scores(pte, yte)
        ac = acc_at_half(pte, yte)
        pb = predict(w, [f(r) for r in band_te])
        ab = auc_scores(pb, [1.0 if r['ok'] else 0.0 for r in band_te])
        ws = " ".join(f"{x:+.3f}" for x in w)
        print(f"{name:>20} {a:>9.4f} {ac:>9.4f} {ab:>9.4f}   [{ws}]")
    print("\nbaseline: always-predict-majority accuracy on test = "
          f"{max(sum(1 for r in te if r['ok']), sum(1 for r in te if not r['ok']))/len(te):.4f}")
    print(f"band test set N = {len(band_te)}, recovery "
          f"{sum(1 for r in band_te if r['ok'])}/{len(band_te)}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: the two-block step of W1b.  step = log2||b*_{m+1}|| -")
    print("        log2||b*_1||; W1b saw it vanish as the wall is crossed.")
    print("-" * 78)
    print(f"{'eff':>6} {'N':>5} {'rec':>9} {'mean step':>10} {'AUC step':>9} "
          f"{'AUC NU':>8} {'AUC nu_hat':>11}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        w = sum(1 for r in sub if r['ok'])
        ms = sum(r['step'] for r in sub) / len(sub)
        if w == 0 or w == len(sub):
            print(f"{eff:>6.3f} {len(sub):>5} {str(w)+'/'+str(len(sub)):>9} "
                  f"{ms:>10.3f} {'(degenerate)':>9}")
            continue
        print(f"{eff:>6.3f} {len(sub):>5} {str(w)+'/'+str(len(sub)):>9} "
              f"{ms:>10.3f} {split_auc(sub, 'step', sign=+1.0):>9.4f} "
              f"{split_auc(sub, 'NU'):>8.4f} {split_auc(sub, 'nuhat'):>11.4f}")
    print(f"\npooled AUC(step, larger=better) = "
          f"{split_auc(rows, 'step', sign=+1.0):.4f}")
    print(f"pooled AUC(-NU)                 = {split_auc(rows, 'NU'):.4f}")
    print(f"pooled AUC(-nu_hat)             = {split_auc(rows, 'nuhat'):.4f}")
    print(f"Spearman(step, log nu_hat)      = "
          f"{spearman([r['step'] for r in rows], [math.log(r['nuhat']) for r in rows]):.4f}")
    print(f"Spearman(step, log NU)          = "
          f"{spearman([r['step'] for r in rows], [math.log(r['NU']) for r in rows]):.4f}")
    print(f"Spearman(step, log eff)         = "
          f"{spearman([r['step'] for r in rows], [math.log(r['eff']) for r in rows]):.4f}")
    print(f"AUC(step) inside the NU band    = "
          f"{split_auc(band, 'step', sign=+1.0):.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X5: X0 showed the recovery rate is NON-MONOTONE in NU (it dips")
    print("        at [1.4,1.8) then rises at [1.8,2.2)).  Is that a real")
    print("        effect or an eff-mixing artifact?  Condition on eff.")
    print("-" * 78)
    print(f"{'eff':>6} {'N':>5} | " + " ".join(f"{'T'+str(i+1):>13}"
                                               for i in range(3)))
    print(f"{'':>6} {'':>5} | " + " ".join(f"{'rate (NU med)':>13}"
                                           for _ in range(3)))
    for eff in EFFS:
        sub = sorted((r for r in rows if r['effq'] == eff),
                     key=lambda r: r['NU'])
        t = len(sub) // 3
        cells = []
        for i in range(3):
            g = sub[i * t:(i + 1) * t] if i < 2 else sub[2 * t:]
            w = sum(1 for r in g if r['ok'])
            med = g[len(g) // 2]['NU']
            cells.append(f"{w/len(g):>6.3f} ({med:4.2f})")
        print(f"{eff:>6.3f} {len(sub):>5} | " + " ".join(f"{c:>13}"
                                                        for c in cells))
    print("\nT1/T2/T3 = NU terciles WITHIN the eff stratum (T1 = smallest NU).")
    print("Monotone-in-NU would give rate(T1) > rate(T2) > rate(T3) in every row.")
    nonmono = 0
    for eff in EFFS:
        sub = sorted((r for r in rows if r['effq'] == eff),
                     key=lambda r: r['NU'])
        t = len(sub) // 3
        rs = []
        for i in range(3):
            g = sub[i * t:(i + 1) * t] if i < 2 else sub[2 * t:]
            rs.append(sum(1 for r in g if r['ok']) / len(g))
        if not (rs[0] >= rs[1] >= rs[2]):
            nonmono += 1
    print(f"strata violating monotonicity: {nonmono}/{len(EFFS)}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X6: if NU enters non-monotonically, the X2 logistic is")
    print("        MISSPECIFIED in log NU.  Does a quadratic NU term rescue")
    print("        the 2-parameter model?")
    print("-" * 78)
    FEATS2 = {
        'nu_hat only':          lambda r: [1.0, math.log(r['nuhat'])],
        'NU^2 + nu_hat':        lambda r: [1.0, math.log(r['NU']),
                                           math.log(r['NU']) ** 2,
                                           math.log(r['nuhat'])],
        'NU^2 only':            lambda r: [1.0, math.log(r['NU']),
                                           math.log(r['NU']) ** 2],
        'nu_hat + eff':         lambda r: [1.0, math.log(r['nuhat']),
                                           math.log(r['eff'])],
        'NU^2 + nu_hat + eff':  lambda r: [1.0, math.log(r['NU']),
                                           math.log(r['NU']) ** 2,
                                           math.log(r['nuhat']),
                                           math.log(r['eff'])],
    }
    print(f"\n{'model':>22} {'test AUC':>9} {'test acc':>9} {'band AUC':>9}")
    for name, f in FEATS2.items():
        w = logistic_fit([f(r) for r in tr],
                         [1.0 if r['ok'] else 0.0 for r in tr])
        pte = predict(w, [f(r) for r in te])
        yte = [1.0 if r['ok'] else 0.0 for r in te]
        pb = predict(w, [f(r) for r in band_te])
        ab = auc_scores(pb, [1.0 if r['ok'] else 0.0 for r in band_te])
        print(f"{name:>22} {auc_scores(pte, yte):>9.4f} "
              f"{acc_at_half(pte, yte):>9.4f} {ab:>9.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X4: size check — repeat X1 at 14 bits (dim 20).")
    print("-" * 78)
    t0 = time.time()
    curves14 = search_curves(1 << 13, 1 << 14, per_bin=2, nbins=10)
    rows14 = build_table(curves14, 10, EFFS, SEEDS10)
    print(f"{len(curves14)} curves, {len(rows14)} instances "
          f"(dim {rows14[0]['k']}) in {time.time()-t0:.1f}s")
    s14 = sorted(r['NU'] for r in rows14 if r['ok'])
    f14 = sorted(r['NU'] for r in rows14 if not r['ok'])
    print(f"14-bit bracket: sufficient NU < {f14[0]:.3f} , "
          f"necessary NU > {s14[-1]:.3f}   "
          f"(FP at NU<=1: {sum(1 for r in rows14 if r['NU'] <= 1 and not r['ok'])})")
    for lo, hi, tag in ((BAND_LO, BAND_HI, "17-bit band"),
                        (f14[0], s14[-1], "own bracket")):
        sub = [r for r in rows14 if lo <= r['NU'] <= hi]
        w = sum(1 for r in sub if r['ok'])
        if not sub or w == 0 or w == len(sub):
            print(f"  {tag} [{lo:.3f},{hi:.3f}] N={len(sub)}  (degenerate)")
            continue
        print(f"  {tag} [{lo:.3f},{hi:.3f}] N={len(sub)} rec {w}/{len(sub)}"
              f"  AUC mu {split_auc(sub, 'mu'):.4f}"
              f"  AUC nu_hat {split_auc(sub, 'nuhat'):.4f}"
              f"  AUC NU {split_auc(sub, 'NU'):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
