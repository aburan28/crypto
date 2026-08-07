"""
GLV-HNP Phase 2, Thread 25: find the second mechanism by conditioning on NU.

Pre-registered by the 2026-08-07 (Thread 24) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  Falsifier: if mu's apparent power is entirely mediated by NU — i.e. the
       in-band AUC collapses toward 0.5 — then the W5 result of Thread 24
       was a stratification artifact and the closed form should be retired.

Why it matters.  Thread 24 (W5/W6) established two facts that do not sit
together comfortably:
  * NU = max_i 2|<e,b*_i>|/||b*_i||^2 is a SOUND certificate for Babai
    nearest-plane (96 TP / 0 FP over 410 instances at two sizes), but a
    DEGRADING separator for Kannan-LLL recovery (AUC 0.978 -> 0.860 from
    12 to 17 bits), and inside a fixed-eff stratum it can be ANTI-predictive
    (AUC 0.350 at eff = 0.15).
  * nu_hat = lambda_1(L2)/sqrt(det L2) — a closed form needing NO lattice
    reduction — holds AUC 0.75-0.93 in every non-degenerate stratum.
  * The two are UNCORRELATED at fixed eff (Spearman -0.28 .. +0.16), so
    nu_hat is not NU in disguise.

Two mutually uncorrelated quantities both predicting recovery means they
govern different things.  This script asks the direct question: with NU held
(approximately) fixed, does mu / nu_hat still separate?  If yes, (NU, nu_hat)
is a genuine 2-coordinate viability test and Phase 2 finally has one.

Experiments:
  X0  build the 17-bit table.  10 seeds x 20 curves x 5 eff strata = 1000
      instances.  The first 5 seeds are Thread 24b's SEEDS verbatim, so
      rows 1-500 are a bit-exact replication of that run.
  X1  H25 proper: AUC(-mu) and AUC(-nu_hat) inside the pre-registered
      ambiguous band [1.040, 2.199], with a bootstrap CI.
  X2  the same conditioned on NU quintiles (equal-count bins) — H25 with the
      band choice removed as a free parameter.
  X3  the converse control: condition on nu_hat quintiles and ask whether NU
      still separates.  If NU dies here but nu_hat survives X2, the ordering
      of the two mechanisms is settled.
  X4  logistic fit on (log NU, log nu_hat): the 2-parameter viability test,
      with its decision boundary and pooled AUC against each single predictor.
  X5  the Thread 24 W1b secondary: step = log2||b*_m|| - log2||b*_0||, the
      height of the jump from the first GS block to the second.  W1b observed
      it vanishing exactly as the K1 wall is crossed.

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
from glv_hnp_phase2_projected import SEEDS, run_new
from glv_hnp_phase2_gsprofile import instance, auc, spearman

# Thread 24b's five seeds, then five more.  Keeping the original five first
# and unchanged makes rows[0:500] a replication of the Thread 24b table.
EXTRA_SEEDS = [7, 2718, 161803, 90210, 4242]
ALL_SEEDS = list(SEEDS) + EXTRA_SEEDS

EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
M17 = 12

# Pre-registered ambiguous band, from Thread 24 W4 at 17 bits:
#   sufficient NU < 1.040 , necessary NU > 2.199
BAND_LO, BAND_HI = 1.040, 2.199


# ---------------------------------------------------------------------------
# small statistics helpers
# ---------------------------------------------------------------------------

def auc_fast(pos, neg):
    """Same value as auc(), via the Mann-Whitney rank identity: O(N log N)
    instead of O(|pos|*|neg|).  Needed for the bootstrap at N ~ 3000."""
    if not pos or not neg:
        return float('nan')
    n1, n2 = len(pos), len(neg)
    allv = sorted([(v, 0) for v in pos] + [(v, 1) for v in neg])
    ranks = [0.0] * len(allv)
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            ranks[t] = avg
        i = j + 1
    r1 = sum(ranks[t] for t in range(len(allv)) if allv[t][1] == 0)
    # U = #{pos < neg} + 0.5*#ties, matching auc()'s "smaller score = positive"
    u = n1 * n2 + n1 * (n1 + 1) / 2.0 - r1
    return u / (n1 * n2)


def auc_boot(pos, neg, nboot=2000, seed=20260807):
    """(auc, lo95, hi95) by nonparametric bootstrap over both classes."""
    a = auc(pos, neg)
    if not pos or not neg:
        return a, float('nan'), float('nan')
    rng = random.Random(seed)
    vals = []
    for _ in range(nboot):
        bp = [pos[rng.randrange(len(pos))] for _ in range(len(pos))]
        bn = [neg[rng.randrange(len(neg))] for _ in range(len(neg))]
        vals.append(auc_fast(bp, bn))
    vals.sort()
    return a, vals[int(0.025 * nboot)], vals[int(0.975 * nboot)]


def quantile_bins(values, nbins):
    """Equal-count bin edges (interior only) for `values`."""
    s = sorted(values)
    return [s[int(round(i * len(s) / nbins))] for i in range(1, nbins)]


def bin_of(x, edges):
    i = 0
    while i < len(edges) and x >= edges[i]:
        i += 1
    return i


def solve(A, b):
    """Gaussian elimination with partial pivoting; A is n x n, b is n."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-300:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for k in range(c, n + 1):
                M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def logistic(X, y, ridge=1e-3, iters=60):
    """Ridge-penalised IRLS.  X includes its own intercept column."""
    d = len(X[0])
    w = [0.0] * d
    for _ in range(iters):
        p = []
        for xi in X:
            z = sum(a * b for a, b in zip(w, xi))
            z = max(-30.0, min(30.0, z))
            p.append(1.0 / (1.0 + math.exp(-z)))
        H = [[ridge if i == j else 0.0 for j in range(d)] for i in range(d)]
        g = [-ridge * w[i] for i in range(d)]
        for xi, pi, yi in zip(X, p, y):
            wi = max(pi * (1.0 - pi), 1e-9)
            for i in range(d):
                g[i] += (yi - pi) * xi[i]
                for j in range(d):
                    H[i][j] += wi * xi[i] * xi[j]
        step = solve(H, g)
        if step is None:
            break
        w = [a + b for a, b in zip(w, step)]
        if max(abs(s) for s in step) < 1e-10:
            break
    return w


def score_of(w, xi):
    return sum(a * b for a, b in zip(w, xi))


# ---------------------------------------------------------------------------
# X0 — build the table
# ---------------------------------------------------------------------------

def build_table(verbose=True):
    t0 = time.time()
    curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
    if verbose:
        print(f"\n{len(curves17)} 17-bit j=0 GLV curves in {time.time()-t0:.1f}s")

    t0 = time.time()
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves17:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in ALL_SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M17, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                m = r['k'] // 2
                p0, pm = r['prof'][0], r['prof'][m]
                r.update({
                    'n': n, 'K1': k1b, 'seed': seed, 'ok': bool(rk['ok']),
                    'eff': k1b * k2b / n, 'effq': eff,
                    'lamstar': lam_star(lam, n),
                    'step': (math.log2(pm) - math.log2(p0)
                             if p0 > 0 and pm > 0 else float('nan')),
                    'headflat': (max(r['prof'][:m]) / min(r['prof'][:m])
                                 if min(r['prof'][:m]) > 0 else float('nan')),
                })
                # the raw GS profile is large and unused downstream; drop it
                # from the dumped rows but keep the derived scalars.
                r.pop('prof', None)
                r.pop('nus', None)
                rows.append(r)
    if verbose:
        print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
              f"in {time.time()-t0:.1f}s")
    return rows


def conditioned_table(rows, key, nbins, scores, title, note=""):
    """AUC of each score within equal-count bins of rows[key]."""
    print("\n" + "-" * 78)
    print(title)
    if note:
        print(note)
    print("-" * 78)
    edges = quantile_bins([r[key] for r in rows], nbins)
    hdr = (f"{key+' bin':>16} {'N':>5} {'rec':>8} | " +
           " ".join(f"{'AUC '+s:>11}" for s, _ in scores))
    print(hdr)
    lo = float('-inf')
    for bi in range(nbins):
        hi = edges[bi] if bi < len(edges) else float('inf')
        sub = [r for r in rows if bin_of(r[key], edges) == bi]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        rng = (f"[{lo:.3f},{hi:.3f})" if lo != float('-inf')
               else f"(-inf,{hi:.3f})")
        if bi == nbins - 1:
            rng = f"[{lo:.3f},inf)"
        cells = []
        for _, f in scores:
            cells.append("     (degen)" if not pos or not neg
                         else f"{auc([f(r) for r in pos], [f(r) for r in neg]):>11.4f}")
        print(f"{rng:>16} {len(sub):>5} "
              f"{str(len(pos))+'/'+str(len(sub)):>8} | " + " ".join(cells))
        lo = hi
    return edges


if __name__ == "__main__":
    dump = None
    if "--dump-json" in sys.argv:
        i = sys.argv.index("--dump-json")
        dump = sys.argv[i + 1] if i + 1 < len(sys.argv) else "nu_conditioned.json"

    print("=" * 78)
    print("Thread 25 — condition on NU: is nu_hat a second coordinate? (H25)")
    print("=" * 78)

    rows = build_table()
    if dump:
        with open(dump, "w") as fh:
            json.dump(rows, fh)
        print(f"table dumped to {dump} ({os.path.getsize(dump)} bytes)")

    # -----------------------------------------------------------------------
    # X0b — replication check against Thread 24b (first 5 seeds only)
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X0b: replication of the Thread 24b W5 table (seeds 42/1234/"
          "9999/555/31337)")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>8} | {'AUC nu_hat':>11} {'AUC NU':>8} "
          f"{'AUC lam*':>9} | {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff and r['seed'] in SEEDS]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>8} | {'(degenerate)':>11}")
            continue
        cells = [auc([r[k] for r in pos], [r[k] for r in neg])
                 for k in ('nuhat', 'NU', 'lamstar', 'mu')]
        print(f"{eff:>5.2f} {len(sub):>5} {str(len(pos))+'/'+str(len(sub)):>8} | "
              f"{cells[0]:>11.4f} {cells[1]:>8.4f} {cells[2]:>9.4f} | "
              f"{cells[3]:>8.4f}")
    print("Thread 24b reported 0.4242/0.8687/0.3232/0.4242 (eff 0.05) ... "
          "0.8816/0.7277/0.1612/0.8816 (eff 0.25)")

    print(f"\nfull table, all {len(ALL_SEEDS)} seeds:")
    print(f"{'eff':>5} {'N':>5} {'rec':>8} | {'AUC nu_hat':>11} {'AUC NU':>8} "
          f"{'AUC lam*':>9} | {'AUC mu':>8}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        if not pos or not neg:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(pos))+'/'+str(len(sub)):>8} | {'(degenerate)':>11}")
            continue
        cells = [auc([r[k] for r in pos], [r[k] for r in neg])
                 for k in ('nuhat', 'NU', 'lamstar', 'mu')]
        print(f"{eff:>5.2f} {len(sub):>5} {str(len(pos))+'/'+str(len(sub)):>8} | "
              f"{cells[0]:>11.4f} {cells[1]:>8.4f} {cells[2]:>9.4f} | "
              f"{cells[3]:>8.4f}")

    # -----------------------------------------------------------------------
    # X1 — H25 proper
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print(f"EXP X1 (H25): inside the pre-registered ambiguous band "
          f"{BAND_LO} <= NU <= {BAND_HI}")
    print("=" * 78)
    below = [r for r in rows if r['NU'] < BAND_LO]
    band = [r for r in rows if BAND_LO <= r['NU'] <= BAND_HI]
    above = [r for r in rows if r['NU'] > BAND_HI]
    for name, sub in (("NU < 1.040  (sufficient)", below),
                      ("in band", band),
                      ("NU > 2.199  (necessary)", above)):
        w = sum(1 for r in sub if r['ok'])
        print(f"  {name:<26} N={len(sub):>5}  recovered {w}/{len(sub)}"
              f" = {w/len(sub) if sub else float('nan'):.3f}")

    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"\nin-band AUC (bootstrap 95% CI, 2000 resamples), N={len(band)}, "
          f"{len(pos)} recovered:")
    for k in ('mu', 'nuhat', 'NU', 'eff', 'lamstar', 'n'):
        a, lo, hi = auc_boot([r[k] for r in pos], [r[k] for r in neg])
        flag = ""
        if k in ('mu', 'nuhat'):
            flag = "   <-- H25 PASS" if lo >= 0.8 else \
                   ("   <-- H25 point-estimate pass, CI crosses 0.8"
                    if a >= 0.8 else "   <-- H25 FAIL")
        print(f"  AUC(-{k:<8}) = {a:.4f}  [{lo:.4f}, {hi:.4f}]{flag}")
    print("\n('n' and 'eff' are controls: if they match nu_hat, the in-band")
    print(" signal is a size/bias effect and not lattice geometry.)")

    # -----------------------------------------------------------------------
    # X2 / X3 — the band choice removed as a free parameter
    # -----------------------------------------------------------------------
    SCORES = [('mu', lambda r: r['mu']), ('nuhat', lambda r: r['nuhat']),
              ('NU', lambda r: r['NU']), ('eff', lambda r: r['eff']),
              ('lam*', lambda r: r['lamstar'])]
    conditioned_table(
        rows, 'NU', 5, SCORES,
        "EXP X2: AUC within NU quintiles — NU is (nearly) constant in each row",
        "AUC > 0.5 means SMALLER score -> more likely to recover.")
    conditioned_table(
        rows, 'nuhat', 5, SCORES,
        "EXP X3: the converse — AUC within nu_hat quintiles",
        "If NU dies here while nu_hat survives X2, nu_hat is the primary axis.")

    # -----------------------------------------------------------------------
    # X4 — the 2-parameter viability test
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X4: logistic fit  logit P(recover) = w0 + w1*log(NU) "
          "+ w2*log(nu_hat)")
    print("-" * 78)
    y = [1.0 if r['ok'] else 0.0 for r in rows]
    feats = {
        'NU only':      lambda r: [1.0, math.log(r['NU'])],
        'nu_hat only':  lambda r: [1.0, math.log(r['nuhat'])],
        'both':         lambda r: [1.0, math.log(r['NU']), math.log(r['nuhat'])],
        'both + eff':   lambda r: [1.0, math.log(r['NU']), math.log(r['nuhat']),
                                   math.log(r['eff'])],
    }
    fits = {}
    print(f"{'model':<14} {'weights':<46} {'AUC':>8} {'acc':>7}")
    for name, f in feats.items():
        X = [f(r) for r in rows]
        w = logistic(X, y)
        s = [score_of(w, xi) for xi in X]
        # logistic score is HIGHER for recovery; auc() wants smaller = recovery
        a = auc([-s[i] for i in range(len(rows)) if rows[i]['ok']],
                [-s[i] for i in range(len(rows)) if not rows[i]['ok']])
        acc = sum(1 for i in range(len(rows))
                  if (s[i] > 0) == rows[i]['ok']) / len(rows)
        fits[name] = w
        ws = " ".join(f"{x:+.4f}" for x in w)
        print(f"{name:<14} {ws:<46} {a:>8.4f} {acc:>7.4f}")

    w = fits['both']
    print(f"\nfitted:  logit p = {w[0]:+.4f} {w[1]:+.4f}*log(NU) "
          f"{w[2]:+.4f}*log(nu_hat)")
    if abs(w[1]) > 1e-9:
        print(f"decision boundary (p = 0.5):  log(NU) = "
              f"{-w[0]/w[1]:+.4f} {-w[2]/w[1]:+.4f}*log(nu_hat)")
        print(f"  i.e.  NU * nu_hat^{w[2]/w[1]:+.4f} = "
              f"{math.exp(-w[0]/w[1]):.4f}   is the viability curve")
    print("\nIn-band accuracy of the 2-parameter rule (the cases NU alone "
          "cannot call):")
    Xb = [feats['both'](r) for r in band]
    sb = [score_of(w, xi) for xi in Xb]
    accb = sum(1 for i in range(len(band)) if (sb[i] > 0) == band[i]['ok']) / len(band)
    base = max(sum(1 for r in band if r['ok']), sum(1 for r in band if not r['ok'])) / len(band)
    print(f"  N={len(band)}  accuracy {accb:.4f}   (majority-class baseline "
          f"{base:.4f})")

    # -----------------------------------------------------------------------
    # X5 — the W1b secondary: the GS-profile step
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X5: step = log2||b*_m|| - log2||b*_0||  (Thread 24 W1b)")
    print("-" * 78)
    print("W1b: the head of the profile is m copies of lambda_1(L2) and the")
    print("step to the second block VANISHES as the K1 wall is crossed, so")
    print("LARGER step should mean easier -> feed -step to the AUC.\n")
    print(f"{'eff':>5} {'N':>5} {'rec':>8} {'step|rec':>10} {'step|fail':>10} "
          f"{'AUC(+step)':>11} {'headflat':>9}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        pos = [r for r in sub if r['ok']]
        neg = [r for r in sub if not r['ok']]
        mp = sum(r['step'] for r in pos) / len(pos) if pos else float('nan')
        mn = sum(r['step'] for r in neg) / len(neg) if neg else float('nan')
        a = (auc([-r['step'] for r in pos], [-r['step'] for r in neg])
             if pos and neg else float('nan'))
        hf = sum(r['headflat'] for r in sub) / len(sub)
        print(f"{eff:>5.2f} {len(sub):>5} {str(len(pos))+'/'+str(len(sub)):>8} "
              f"{mp:>10.4f} {mn:>10.4f} {a:>11.4f} {hf:>9.4f}")
    p_all = [-r['step'] for r in rows if r['ok']]
    n_all = [-r['step'] for r in rows if not r['ok']]
    print(f"\npooled AUC(+step -> recovery) = {auc(p_all, n_all):.4f}   "
          f"N={len(rows)}")
    print(f"Spearman(step, log nu_hat) = "
          f"{spearman([r['step'] for r in rows], [math.log(r['nuhat']) for r in rows]):.4f}")
    print(f"Spearman(step, log NU)     = "
          f"{spearman([r['step'] for r in rows], [math.log(r['NU']) for r in rows]):.4f}")
    print(f"Spearman(step, eff)        = "
          f"{spearman([r['step'] for r in rows], [r['eff'] for r in rows]):.4f}")
    print("\nheadflat = max/min of the first m GS norms; W1b claims ~1.000 "
          "(m exact copies of lambda_1(L2)).")
    hf = [r['headflat'] for r in rows]
    print(f"headflat over all {len(rows)}: mean {sum(hf)/len(hf):.4f}  "
          f"min {min(hf):.4f}  max {max(hf):.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
