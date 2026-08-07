"""
GLV-HNP Phase 2, Thread 25: is mu a genuine SECOND coordinate, or is its
apparent predictive power mediated by NU after all?

Thread 24 (W5/W6) established two facts that cannot both be about one
mechanism:

  * NU = max_i 2|<e,b*_i>|/||b*_i||^2 is a SOUND certificate for Babai
    nearest-plane (96 TP / 0 FP over 410 instances) but a size-DEGRADING
    separator for Kannan-LLL recovery (AUC 0.978 -> 0.860, 12 -> 17 bits).
  * nu_hat = lambda_1(L2)/sqrt(det L2) — which needs no lattice reduction —
    beats it cross-curve at fixed eff (pooled AUC 0.935 vs 0.800), and the
    two are rank-UNCORRELATED inside every eff stratum (Spearman -0.28..+0.16).

Pre-registered hypothesis, verbatim from the Thread 24 next-step proposal:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

  If yes, mu is a genuine second coordinate and the pair (NU, mu) is the
  2-parameter viability test Phase 2 has been missing.  If no — if mu's
  apparent power is entirely mediated by NU after all — then the W5 result is
  a stratification artifact and the closed form should be retired.

Experiments
  X1  H25 proper: AUC of each score inside each NU band.
  X2  the strict version: NU-band x eff-stratum, so BOTH the certificate and
      the bias strength are held fixed and only cross-curve geometry is left.
  X3  logistic fit on (log NU, log nu_hat) by IRLS, with likelihood-ratio
      tests against each single-predictor nested model.  This is the
      mediation test: if the log nu_hat coefficient survives conditioning on
      log NU, mu is a real second coordinate.
  X4  the W1b secondary: step = log2||b*_{m+1}|| - log2||b*_1||, the height of
      the two-block jump in the GS profile, which W1b observed vanishing
      exactly as the K1 wall is crossed.  Does step beat NU or nu_hat?
  X5  out-of-sample validation: fit the X3 boundary on 17 bits, score 12-bit
      instances the model has never seen.  A 2-parameter law that does not
      transfer across sizes is a curve fit.

Gram-Schmidt is float, justified by W0 of the parent script (max relative NU
error vs exact Fractions 6.6e-16 at dim 24).

Run: python3 glv_hnp_phase2_conditional.py [--json PATH]
     --json reuses/creates the Thread-24b dumped table instead of regenerating.
"""

import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_common import search_curves
from glv_hnp_phase2_gsprofile import auc, spearman
from glv_hnp_phase2_gsprofile_strat import (
    EFFS, M17, build_table, dump_json, load_json,
)

# The Thread 24 W4 bracket at 17 bits: NU < 1.040 is sufficient for
# nearest-plane, NU > 2.199 was never observed to recover.  The open band
# between them is where a second coordinate would have to live.
BAND_LO, BAND_HI = 1.040, 2.199


# ---------------------------------------------------------------------------
# Logistic regression by IRLS, with a small ridge so a separating hyperplane
# cannot send the coefficients to infinity.
# ---------------------------------------------------------------------------

def logistic_irls(X, y, ridge=1e-3, iters=60):
    """X: list of feature rows (WITHOUT intercept).  y: 0/1.

    Returns (beta, loglik) with beta[0] the intercept.
    """
    n = len(X)
    d = len(X[0]) + 1
    A = [[1.0] + list(row) for row in X]
    beta = [0.0] * d
    for _ in range(iters):
        eta = [sum(b * a for b, a in zip(beta, row)) for row in A]
        mu = [1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, e)))) for e in eta]
        w = [max(m * (1.0 - m), 1e-9) for m in mu]
        # gradient and Hessian of the ridge-penalised log-likelihood
        g = [sum(A[i][j] * (y[i] - mu[i]) for i in range(n)) - ridge * beta[j]
             for j in range(d)]
        H = [[sum(w[i] * A[i][j] * A[i][k] for i in range(n))
              + (ridge if j == k else 0.0) for k in range(d)] for j in range(d)]
        step = _solve(H, g)
        if step is None:
            break
        beta = [b + s for b, s in zip(beta, step)]
        if max(abs(s) for s in step) < 1e-10:
            break
    eta = [sum(b * a for b, a in zip(beta, row)) for row in A]
    ll = 0.0
    for i in range(n):
        e = max(-30.0, min(30.0, eta[i]))
        ll += y[i] * e - math.log1p(math.exp(e))
    return beta, ll


def _solve(M, v):
    """Gaussian elimination with partial pivoting."""
    d = len(v)
    a = [list(M[i]) + [v[i]] for i in range(d)]
    for c in range(d):
        piv = max(range(c, d), key=lambda r: abs(a[r][c]))
        if abs(a[piv][c]) < 1e-14:
            return None
        a[c], a[piv] = a[piv], a[c]
        for r in range(d):
            if r == c:
                continue
            f = a[r][c] / a[c][c]
            for k in range(c, d + 1):
                a[r][k] -= f * a[c][k]
    return [a[i][d] / a[i][i] for i in range(d)]


def chi2_p_1df(stat):
    """Upper tail of chi^2_1 = erfc(sqrt(stat/2))."""
    return math.erfc(math.sqrt(max(stat, 0.0) / 2.0))


def lr_test(ll_full, ll_sub, df=1):
    stat = 2.0 * (ll_full - ll_sub)
    return stat, (chi2_p_1df(stat) if df == 1 else float('nan'))


# ---------------------------------------------------------------------------
# The GS block step (W1b secondary)
# ---------------------------------------------------------------------------

def block_step(r):
    """log2||b*_{m+1}|| - log2||b*_1||.

    W1b: the head of the LLL GS profile is m exact copies of lambda_1(L2) and
    the jump to the second block collapses to 0 exactly as the wall is
    crossed.  Positive step = two clean blocks = below the wall.
    """
    prof, m = r['prof'], r['m']
    if len(prof) <= m or prof[0] <= 0 or prof[m] <= 0:
        return float('nan')
    return math.log2(prof[m]) - math.log2(prof[0])


def band_of(nu):
    if nu < BAND_LO:
        return 'NU<1.04'
    if nu <= BAND_HI:
        return '1.04-2.20'
    return 'NU>2.20'


BANDS = ('NU<1.04', '1.04-2.20', 'NU>2.20')

SCORES = (
    ('mu', lambda r: r['mu']),
    ('nu_hat', lambda r: r['nuhat']),
    ('NU', lambda r: r['NU']),
    ('step', lambda r: -block_step(r)),   # negated: LARGER step = better
    ('eff', lambda r: r['eff']),          # confound control
    ('lam*', lambda r: r['lamstar']),     # Thread 20 null control
)


def auc_row(sub, scores=SCORES):
    pos = [r for r in sub if r['ok']]
    neg = [r for r in sub if not r['ok']]
    if not pos or not neg:
        return None
    return [auc([f(r) for r in pos], [f(r) for r in neg]) for _, f in scores]


def main():
    json_path = None
    if "--json" in sys.argv:
        json_path = sys.argv[sys.argv.index("--json") + 1]

    print("=" * 78)
    print("Thread 25 — is mu a second coordinate, or is it mediated by NU?")
    print("=" * 78)

    if json_path and os.path.exists(json_path):
        rows = load_json(json_path)
        print(f"\nloaded {len(rows)} instances from {json_path}")
    else:
        t0 = time.time()
        curves17 = search_curves(1 << 16, 1 << 17, per_bin=2, nbins=10)
        rows = build_table(curves17, M17, EFFS)
        print(f"\n{len(curves17)} 17-bit curves, {len(rows)} instances "
              f"(dim {rows[0]['k']}) in {time.time()-t0:.1f}s")
        if json_path:
            dump_json(rows, json_path)
            print(f"table dumped to {json_path}")

    hdr = f"{'band':>10} {'N':>5} {'rec':>8} |" + \
          "".join(f"{nm:>9}" for nm, _ in SCORES)

    # ---------------- X1 ----------------
    print("\n" + "-" * 78)
    print("EXP X1: H25 proper — AUC inside each NU band.")
    print("        AUC > 0.5 means SMALLER score -> more likely to recover.")
    print("-" * 78)
    print(hdr)
    for bd in BANDS:
        sub = [r for r in rows if band_of(r['NU']) == bd]
        if not sub:
            continue
        nrec = sum(1 for r in sub if r['ok'])
        a = auc_row(sub)
        cells = ("".join(f"{x:>9.4f}" for x in a) if a
                 else f"{'(degenerate)':>27}")
        print(f"{bd:>10} {len(sub):>5} {str(nrec)+'/'+str(len(sub)):>8} |{cells}")

    print(f"\nH25 threshold: AUC(-mu) >= 0.80 inside the 1.04-2.20 band.")

    # ---------------- X2 ----------------
    print("\n" + "-" * 78)
    print("EXP X2: the strict test — NU band x eff stratum.  Both the")
    print("        certificate and the bias strength are held fixed.")
    print("-" * 78)
    print(f"{'eff':>5} " + hdr[10:])
    for bd in BANDS:
        print(f"  -- band {bd} --")
        for eff in EFFS:
            sub = [r for r in rows
                   if band_of(r['NU']) == bd and r['effq'] == eff]
            if not sub:
                continue
            nrec = sum(1 for r in sub if r['ok'])
            a = auc_row(sub)
            cells = ("".join(f"{x:>9.4f}" for x in a) if a
                     else f"{'(one class only)':>30}")
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(nrec)+'/'+str(len(sub)):>8} |{cells}")

    # ---------------- X3 ----------------
    print("\n" + "-" * 78)
    print("EXP X3: mediation test — logistic fit, recovery ~ log NU + log nu_hat")
    print("-" * 78)
    y = [1 if r['ok'] else 0 for r in rows]
    lnu = [math.log(r['NU']) for r in rows]
    lnh = [math.log(r['nuhat']) for r in rows]
    lst = [block_step(r) for r in rows]

    b_full, ll_full = logistic_irls([[a, b] for a, b in zip(lnu, lnh)], y)
    b_nu, ll_nu = logistic_irls([[a] for a in lnu], y)
    b_nh, ll_nh = logistic_irls([[a] for a in lnh], y)
    _, ll_null = logistic_irls([[0.0] for _ in y], y)

    print(f"{'model':>28} {'loglik':>11} {'coefficients':>34}")
    print(f"{'intercept only':>28} {ll_null:>11.3f}")
    print(f"{'~ log NU':>28} {ll_nu:>11.3f}  "
          f"b0={b_nu[0]:+.3f}  b_NU={b_nu[1]:+.3f}")
    print(f"{'~ log nu_hat':>28} {ll_nh:>11.3f}  "
          f"b0={b_nh[0]:+.3f}  b_nh={b_nh[1]:+.3f}")
    print(f"{'~ log NU + log nu_hat':>28} {ll_full:>11.3f}  "
          f"b0={b_full[0]:+.3f}  b_NU={b_full[1]:+.3f}  b_nh={b_full[2]:+.3f}")

    s1, p1 = lr_test(ll_full, ll_nu)
    s2, p2 = lr_test(ll_full, ll_nh)
    print(f"\nLR add log nu_hat to log NU : chi2_1 = {s1:9.3f}   p = {p1:.3e}")
    print(f"LR add log NU to log nu_hat : chi2_1 = {s2:9.3f}   p = {p2:.3e}")
    print("Both significant => neither mediates the other => 2 coordinates.")

    fit = [b_full[0] + b_full[1] * a + b_full[2] * b for a, b in zip(lnu, lnh)]
    pos = [-fit[i] for i in range(len(rows)) if y[i]]
    neg = [-fit[i] for i in range(len(rows)) if not y[i]]
    print(f"\npooled AUC(fitted 2-param score) = {auc(pos, neg):.4f}")
    print(f"pooled AUC(-log NU)              = "
          f"{auc([lnu[i] for i in range(len(rows)) if y[i]], [lnu[i] for i in range(len(rows)) if not y[i]]):.4f}")
    print(f"pooled AUC(-log nu_hat)          = "
          f"{auc([lnh[i] for i in range(len(rows)) if y[i]], [lnh[i] for i in range(len(rows)) if not y[i]]):.4f}")
    if abs(b_full[2]) > 1e-9:
        print(f"\ndecision boundary (p=0.5):  log nu_hat = "
              f"{-b_full[0]/b_full[2]:+.4f} {-b_full[1]/b_full[2]:+.4f} * log NU")

    # ---------------- X4 ----------------
    print("\n" + "-" * 78)
    print("EXP X4: the W1b block step as a predictor in its own right")
    print("-" * 78)
    fin = [s for s in lst if s == s]
    print(f"step over all {len(fin)} instances: "
          f"min {min(fin):+.3f}  median {sorted(fin)[len(fin)//2]:+.3f}  "
          f"max {max(fin):+.3f}")
    okst = [lst[i] for i in range(len(rows)) if y[i] and lst[i] == lst[i]]
    nost = [lst[i] for i in range(len(rows)) if not y[i] and lst[i] == lst[i]]
    print(f"  step | success : median {sorted(okst)[len(okst)//2]:+.3f} "
          f"(n={len(okst)})")
    print(f"  step | failure : median {sorted(nost)[len(nost)//2]:+.3f} "
          f"(n={len(nost)})")
    print(f"\nSpearman(step, log nu_hat) = {spearman(lst, lnh):+.4f}")
    print(f"Spearman(step, log NU)     = {spearman(lst, lnu):+.4f}")
    b_st, ll_st = logistic_irls([[a] for a in lst], y)
    b_3, ll_3 = logistic_irls([[a, b, c] for a, b, c in zip(lnu, lnh, lst)], y)
    s3, p3 = lr_test(ll_3, ll_full)
    print(f"~ step alone           loglik {ll_st:>10.3f}  b_step={b_st[1]:+.3f}")
    print(f"~ NU + nu_hat + step   loglik {ll_3:>10.3f}")
    print(f"LR add step to the pair: chi2_1 = {s3:.3f}   p = {p3:.3e}")

    # ---------------- X5 ----------------
    print("\n" + "-" * 78)
    print("EXP X5: out-of-sample — fit on 17 bits, score UNSEEN 12-bit data")
    print("-" * 78)
    t0 = time.time()
    curves12 = search_curves(1 << 11, 1 << 12, per_bin=2, nbins=10)
    rows12 = build_table(curves12, 8, EFFS)
    print(f"{len(curves12)} 12-bit curves, {len(rows12)} instances "
          f"(dim {rows12[0]['k']}) in {time.time()-t0:.1f}s")
    y12 = [1 if r['ok'] else 0 for r in rows12]
    nrec12 = sum(y12)
    if nrec12 in (0, len(y12)):
        print("held-out set is single-class; no AUC available")
    else:
        f12 = [b_full[0] + b_full[1] * math.log(r['NU'])
               + b_full[2] * math.log(r['nuhat']) for r in rows12]
        p12 = [-f12[i] for i in range(len(rows12)) if y12[i]]
        n12 = [-f12[i] for i in range(len(rows12)) if not y12[i]]
        print(f"recovered {nrec12}/{len(rows12)}")
        print(f"  transferred 2-param score AUC = {auc(p12, n12):.4f}")
        a = auc_row(rows12)
        print("  in-sample-free single scores:  " +
              "  ".join(f"{nm} {v:.4f}" for (nm, _), v in zip(SCORES, a)))
        acc = sum(1 for i in range(len(rows12))
                  if (f12[i] > 0) == bool(y12[i])) / len(rows12)
        print(f"  transferred boundary accuracy (p=0.5 cut) = {acc:.4f}  "
              f"(base rate {max(nrec12, len(y12)-nrec12)/len(y12):.4f})")

    # ---------------- X6 ----------------
    print("\n" + "-" * 78)
    print("EXP X6: the X3 boundary has slope -1.05, i.e. it is almost exactly")
    print("        log NU + log nu_hat = const.  Test the product law directly.")
    print("-" * 78)
    print(f"fitted slope  b_NU/b_nh = {b_full[1]/b_full[2]:.4f}   (product law => 1)")
    print(f"fitted cut    NU*nu_hat = {math.exp(-b_full[0]/b_full[2]):.4f}   "
          f"(=> the natural constant 1?)")

    def prod(r):
        return r['NU'] * r['nuhat']

    for label, rr, yy in (("17 bits (in-sample)", rows, y),
                          ("12 bits (held out)", rows12, y12)):
        pp = [prod(r) for r in rr]
        po = [pp[i] for i in range(len(rr)) if yy[i]]
        ne = [pp[i] for i in range(len(rr)) if not yy[i]]
        if not po or not ne:
            continue
        print(f"\n{label}   N={len(rr)}  rec={sum(yy)}")
        print(f"  AUC(-NU*nu_hat)        = {auc(po, ne):.4f}")
        print(f"  product | success : min {min(po):.4f}  "
              f"median {sorted(po)[len(po)//2]:.4f}  max {max(po):.4f}")
        print(f"  product | failure : min {min(ne):.4f}  "
              f"median {sorted(ne)[len(ne)//2]:.4f}  max {max(ne):.4f}")
        print(f"  {'cut':>8} {'acc':>8} {'TP':>6} {'FP':>6} {'FN':>6} {'TN':>6}")
        for cut in (0.5, 0.75, 0.961, 1.0, 1.25, 1.5, 2.0):
            tp = sum(1 for i in range(len(rr)) if yy[i] and pp[i] < cut)
            fp = sum(1 for i in range(len(rr)) if not yy[i] and pp[i] < cut)
            fn = sum(yy) - tp
            tn = len(rr) - sum(yy) - fp
            print(f"  {cut:>8.3f} {(tp+tn)/len(rr):>8.4f} "
                  f"{tp:>6} {fp:>6} {fn:>6} {tn:>6}")

    print("\n" + "-" * 78)
    print("EXP X6b: certificate bracket for the product, and the fixed-eff")
    print("         control (the test NU fails and nu_hat passes).")
    print("-" * 78)
    for label, rr, yy in (("17 bits", rows, y), ("12 bits", rows12, y12)):
        pp = [prod(r) for r in rr]
        po = [pp[i] for i in range(len(rr)) if yy[i]]
        ne = [pp[i] for i in range(len(rr)) if not yy[i]]
        if not po or not ne:
            continue
        print(f"{label}: sufficient NU*nu_hat < {min(ne):.4f}   "
              f"necessary NU*nu_hat < {max(po):.4f}   "
              f"ambiguity {max(po)/min(ne):.2f}x")
    print("  (Thread 24 W4, NU alone: 17 bits [1.040, 2.199] = 2.11x, "
          "12 bits [1.188, 1.869] = 1.57x)")

    print(f"\n{'eff':>5} {'N':>5} {'rec':>8} | {'AUC product':>12} "
          f"{'AUC NU':>8} {'AUC nu_hat':>11}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        po = [r for r in sub if r['ok']]
        ne = [r for r in sub if not r['ok']]
        if not po or not ne:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(po))+'/'+str(len(sub)):>8} | {'(degenerate)':>12}")
            continue
        print(f"{eff:>5.2f} {len(sub):>5} "
              f"{str(len(po))+'/'+str(len(sub)):>8} | "
              f"{auc([prod(r) for r in po], [prod(r) for r in ne]):>12.4f} "
              f"{auc([r['NU'] for r in po], [r['NU'] for r in ne]):>8.4f} "
              f"{auc([r['nuhat'] for r in po], [r['nuhat'] for r in ne]):>11.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
