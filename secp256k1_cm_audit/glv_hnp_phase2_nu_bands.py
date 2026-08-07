"""
GLV-HNP Phase 2, Thread 25: is mu a SECOND coordinate, or is its apparent
power mediated by NU?

Pre-registered by the 2026-08-07 (autolab run #2, Thread 24) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where nearest-plane
       gives no answer), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
  If yes  -> mu is a genuine second coordinate; deliverable is a logistic fit
             on (log NU, log mu) with its decision boundary.
  If no   -> W5's result was a stratification artifact and the closed form
             should be retired.

  Secondary (also pre-registered): define
       step = log2(||b*_{m+1}||) - log2(||b*_1||)
  and test whether step -> 0 predicts the wall better than either NU or mu.
  W1b found the GS profile is m exact copies of lambda_1(L2) followed by a
  second block, and that the step between blocks vanishes right as the K1
  wall is crossed.

Why the conditioning matters.  Thread 24 (W5/W6) showed nu_hat = mu/sqrt(det
L2) out-AUCs the exact BDD certificate NU inside every eff stratum, and that
C = NU/(nu_hat*sqrt(eff)) has Spearman ~0 (or negative) against NU inside a
stratum.  Two uncorrelated quantities both predicting recovery means either
(a) two mechanisms, or (b) one of them is riding a lurking variable.  eff was
already held fixed in W5, so the remaining lurking variable is NU itself.
Conditioning on NU is the clean test.

Numerics.  Float Gram-Schmidt, justified by W0/W4 of glv_hnp_phase2_gsprofile
(max relative NU error vs exact Fractions 6.6e-16 at dim 24 / 17 bits).

Run:  python3 glv_hnp_phase2_nu_bands.py [--seeds N] [--dump-json PATH]
"""

import argparse
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

EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
M17 = 12
# dimension is held at 2*M17 = 24 for EVERY bit-length in X5, so the only
# thing varying across the size sweep is n.
SIZES = (12, 14, 17, 20)

# The 17-bit bracket measured by Thread 24 W4 on 300 instances.  X1 recomputes
# it on the present table; these are the pre-registered values.
BAND_LO_PRE, BAND_HI_PRE = 1.040, 2.199


# ---------------------------------------------------------------------------
# scores
# ---------------------------------------------------------------------------

def auc_lo(rows, key):
    """AUC for 'smaller score -> recovery'."""
    pos = [key(r) for r in rows if r['ok']]
    neg = [key(r) for r in rows if not r['ok']]
    return auc(pos, neg)


def auc_hi(rows, key):
    """AUC for 'larger score -> recovery'."""
    return auc_lo(rows, lambda r: -key(r))


# ---------------------------------------------------------------------------
# logistic regression, IRLS with ridge (pure python, 2-3 params)
# ---------------------------------------------------------------------------

def solve(A, b):
    """Gaussian elimination with partial pivoting."""
    k = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(k):
        piv = max(range(c, k), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-14:
            return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(k):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for j in range(c, k + 1):
                M[r][j] -= f * M[c][j]
    return [M[i][k] / M[i][i] for i in range(k)]


def logistic(X, y, ridge=1e-4, iters=60):
    """IRLS.  X rows already include the leading 1.0.  Returns (w, loglik)."""
    k = len(X[0])
    w = [0.0] * k
    for _ in range(iters):
        H = [[ridge if i == j else 0.0 for j in range(k)] for i in range(k)]
        g = [0.0] * k
        for xi, yi in zip(X, y):
            z = sum(a * b for a, b in zip(w, xi))
            z = max(-30.0, min(30.0, z))
            pi = 1.0 / (1.0 + math.exp(-z))
            s = pi * (1.0 - pi)
            for i in range(k):
                g[i] += (yi - pi) * xi[i]
                for j in range(k):
                    H[i][j] += s * xi[i] * xi[j]
        for i in range(k):
            g[i] -= ridge * w[i]
        step = solve(H, g)
        if step is None:
            break
        w = [a + b for a, b in zip(w, step)]
        if max(abs(s) for s in step) < 1e-10:
            break
    ll = 0.0
    for xi, yi in zip(X, y):
        z = max(-30.0, min(30.0, sum(a * b for a, b in zip(w, xi))))
        pi = 1.0 / (1.0 + math.exp(-z))
        pi = min(max(pi, 1e-12), 1 - 1e-12)
        ll += yi * math.log(pi) + (1 - yi) * math.log(1 - pi)
    return w, ll


def fit_report(rows, feats, name):
    """Fit, report coefficients, in-sample AUC and accuracy."""
    X = [[1.0] + [f(r) for f in feats] for r in rows]
    y = [1.0 if r['ok'] else 0.0 for r in rows]
    w, ll = logistic(X, y)
    scored = []
    for xi, r in zip(X, rows):
        z = sum(a * b for a, b in zip(w, xi))
        scored.append((z, r['ok']))
    pos = [z for z, o in scored if o]
    neg = [z for z, o in scored if not o]
    a = auc([-z for z in pos], [-z for z in neg])
    acc = sum(1 for z, o in scored if (z > 0) == o) / len(scored)
    base = max(sum(y), len(y) - sum(y)) / len(y)
    print(f"  {name:<26} w = [" +
          ", ".join(f"{c:+.4f}" for c in w) + f"]  loglik {ll:8.2f}  "
          f"AUC {a:.4f}  acc {acc:.4f} (base {base:.4f})")
    return w, ll, a


# ---------------------------------------------------------------------------

def collect(nseeds, bits=17, verbose=True):
    t0 = time.time()
    curves = search_curves(1 << (bits - 1), 1 << bits, per_bin=2, nbins=10)
    if verbose:
        print(f"\n{len(curves)} {bits}-bit j=0 GLV curves in "
              f"{time.time()-t0:.1f}s")
    seeds = (SEEDS * 4)[:nseeds] if nseeds <= len(SEEDS) else \
        SEEDS + [10007 * i + 3 for i in range(nseeds - len(SEEDS))]
    seeds = seeds[:nseeds]
    t0 = time.time()
    rows = []
    for eff in EFFS:
        for (p, b, n, lam, G) in curves:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            for seed in seeds:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                r = instance((p, b, n, lam, G), M17, d_trial, k1b, seed,
                             exact=False)
                if r is None:
                    continue
                rk = run_new((p, b, n, lam, G), M17, d_trial, k1b, seed)
                prof = r['prof']
                m = M17
                # pre-registered pointwise step, and its block-mean sibling
                st = (math.log2(prof[m]) - math.log2(prof[0])
                      if prof[m] > 0 and prof[0] > 0 else 0.0)
                hb = sum(prof[:m]) / m
                tb = sum(prof[m:]) / (len(prof) - m)
                stb = math.log2(tb) - math.log2(hb) if hb > 0 and tb > 0 else 0.0
                r.update({'n': n, 'K1': k1b, 'ok': bool(rk['ok']),
                          'eff': k1b * k2b / n, 'effq': eff, 'seed': seed,
                          'lamstar': lam_star(lam, n), 'step': st,
                          'stepb': stb, 'bits': bits,
                          'psi': r['NU'] * r['nuhat']})
                rows.append(r)
    if verbose:
        print(f"{len(rows)} instances (float GS, dim {rows[0]['k']}) "
              f"in {time.time()-t0:.1f}s")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--dump-json", default="")
    args = ap.parse_args()

    print("=" * 78)
    print("Thread 25 — condition on NU: is mu a second coordinate?  (H25)")
    print("=" * 78)

    rows = collect(args.seeds)

    if args.dump_json:
        # everything else (det2 = n*S_K1*S_K2, l2 = det2/mu, psi = NU*nuhat,
        # eff = K1*K2/n) is derivable from these columns.
        ints = ('n', 'K1', 'K2', 'seed', 'k', 'bits')
        flts = ('effq', 'NU', 'mu', 'nuhat', 'lamstar', 'step', 'stepb',
                'enorm')
        with open(args.dump_json, "w") as fh:
            json.dump([dict([(k, r[k]) for k in ints]
                            + [(k, float(f"{r[k]:.6g}")) for k in flts]
                            + [('ok', r['ok'])]) for r in rows], fh)
        print(f"dumped {len(rows)} rows -> {args.dump_json} "
              f"({os.path.getsize(args.dump_json)/1024:.0f} KB)")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X0: recompute the NU bracket on this table")
    print("-" * 78)
    P = [r['NU'] for r in rows if r['ok']]
    N = [r['NU'] for r in rows if not r['ok']]
    print(f"NU | success : min {min(P):.3f}  median {sorted(P)[len(P)//2]:.3f}"
          f"  max {max(P):.3f}  (n={len(P)})")
    print(f"NU | failure : min {min(N):.3f}  median {sorted(N)[len(N)//2]:.3f}"
          f"  max {max(N):.3f}  (n={len(N)})")
    tp = sum(1 for r in rows if r['NU'] <= 1.0 and r['ok'])
    fp = sum(1 for r in rows if r['NU'] <= 1.0 and not r['ok'])
    print(f"NU <= 1 certificate: TP {tp}  FP {fp}")
    lo, hi = min(N), max(P)
    print(f"observed bracket : sufficient NU < {lo:.3f} , necessary NU > {hi:.3f}")
    print(f"pre-registered   : [{BAND_LO_PRE:.3f}, {BAND_HI_PRE:.3f}]  "
          f"(Thread 24 W4, N=300)")
    print(f"AUC(-NU) overall = {auc_lo(rows, lambda r: r['NU']):.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1: THE TEST — AUC inside the ambiguous NU band")
    print("-" * 78)
    for tag, blo, bhi in (("pre-registered", BAND_LO_PRE, BAND_HI_PRE),
                          ("observed", lo, hi)):
        band = [r for r in rows if blo <= r['NU'] <= bhi]
        w = sum(1 for r in band if r['ok'])
        print(f"\nband {tag} [{blo:.3f}, {bhi:.3f}] : N={len(band)}  "
              f"recovered {w}/{len(band)}")
        if not w or w == len(band):
            print("  (degenerate — no discrimination possible)")
            continue
        print(f"  AUC(-mu)               {auc_lo(band, lambda r: r['mu']):.4f}"
              "   <-- H25 target >= 0.80")
        print(f"  AUC(-nu_hat)           "
              f"{auc_lo(band, lambda r: r['nuhat']):.4f}")
        print(f"  AUC(-nu_hat*sqrt eff)  "
              f"{auc_lo(band, lambda r: r['nuhat']*math.sqrt(r['eff'])):.4f}")
        print(f"  AUC(-NU) [residual]    {auc_lo(band, lambda r: r['NU']):.4f}"
              "   (should be ~0.5 if the band is narrow)")
        print(f"  AUC(+step)             {auc_hi(band, lambda r: r['step']):.4f}")
        print(f"  AUC(+step_block)       "
              f"{auc_hi(band, lambda r: r['stepb']):.4f}")
        print(f"  AUC(-eff)  [control]   {auc_lo(band, lambda r: r['eff']):.4f}")
        print(f"  AUC(-lam*) [control]   "
              f"{auc_lo(band, lambda r: r['lamstar']):.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X1b: DOUBLE conditioning — NU band x eff stratum")
    print("-" * 78)
    print("eff is fixed within a row and NU is bracketed within a column, so")
    print("neither can be the lurking variable.\n")
    print(f"{'eff':>5} {'band':>16} {'N':>5} {'rec':>7} {'AUC mu':>8} "
          f"{'AUC nuhat':>10} {'AUC step':>9}")
    for eff in EFFS:
        for tag, blo, bhi in (("NU<lo", 0.0, BAND_LO_PRE),
                              ("ambiguous", BAND_LO_PRE, BAND_HI_PRE),
                              ("NU>hi", BAND_HI_PRE, 1e9)):
            sub = [r for r in rows
                   if r['effq'] == eff and blo <= r['NU'] <= bhi]
            if not sub:
                continue
            w = sum(1 for r in sub if r['ok'])
            if w == 0 or w == len(sub):
                print(f"{eff:>5.2f} {tag:>16} {len(sub):>5} "
                      f"{str(w)+'/'+str(len(sub)):>7} {'--':>8} {'--':>10} "
                      f"{'--':>9}")
                continue
            print(f"{eff:>5.2f} {tag:>16} {len(sub):>5} "
                  f"{str(w)+'/'+str(len(sub)):>7} "
                  f"{auc_lo(sub, lambda r: r['mu']):>8.4f} "
                  f"{auc_lo(sub, lambda r: r['nuhat']):>10.4f} "
                  f"{auc_hi(sub, lambda r: r['step']):>9.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X2: logistic fit — does mu add anything on top of NU?")
    print("-" * 78)
    print("full table:")
    _, ll0, _ = fit_report(rows, [lambda r: math.log(r['NU'])], "log NU")
    _, ll1, _ = fit_report(rows, [lambda r: math.log(r['nuhat'])], "log nu_hat")
    w2, ll2, _ = fit_report(rows, [lambda r: math.log(r['NU']),
                                   lambda r: math.log(r['nuhat'])],
                            "log NU + log nu_hat")
    _, ll3, _ = fit_report(rows, [lambda r: math.log(r['NU']),
                                  lambda r: math.log(r['mu'])],
                           "log NU + log mu")
    _, ll4, _ = fit_report(rows, [lambda r: math.log(r['NU']),
                                  lambda r: math.log(r['nuhat']),
                                  lambda r: r['step']],
                           "log NU + log nu_hat + step")
    print(f"\n  LR test  (NU) vs (NU, nu_hat) : chi2_1 = {2*(ll2-ll0):8.2f} "
          f"-> {'SIGNIFICANT' if 2*(ll2-ll0) > 3.84 else 'n.s.'} at p=0.05")
    print(f"  LR test  (NU) vs (NU, mu)     : chi2_1 = {2*(ll3-ll0):8.2f} "
          f"-> {'SIGNIFICANT' if 2*(ll3-ll0) > 3.84 else 'n.s.'} at p=0.05")
    print(f"  LR test  (NU,nu_hat) vs +step : chi2_1 = {2*(ll4-ll2):8.2f} "
          f"-> {'SIGNIFICANT' if 2*(ll4-ll2) > 3.84 else 'n.s.'} at p=0.05")

    band = [r for r in rows if BAND_LO_PRE <= r['NU'] <= BAND_HI_PRE]
    if band and 0 < sum(1 for r in band if r['ok']) < len(band):
        print(f"\nambiguous band only (N={len(band)}):")
        _, bl0, _ = fit_report(band, [lambda r: math.log(r['NU'])], "log NU")
        _, bl1, _ = fit_report(band, [lambda r: math.log(r['nuhat'])],
                               "log nu_hat")
        wb, bl2, _ = fit_report(band, [lambda r: math.log(r['NU']),
                                       lambda r: math.log(r['nuhat'])],
                                "log NU + log nu_hat")
        print(f"\n  LR test  (NU) vs (NU, nu_hat) : chi2_1 = "
              f"{2*(bl2-bl0):8.2f} "
              f"-> {'SIGNIFICANT' if 2*(bl2-bl0) > 3.84 else 'n.s.'}")
        if abs(wb[2]) > 1e-9:
            print(f"\n  decision boundary (p=0.5) inside the band:")
            print(f"    {wb[0]:+.4f} {wb[1]:+.4f}*log NU {wb[2]:+.4f}"
                  f"*log nu_hat = 0")
            print(f"    log nu_hat = {-wb[0]/wb[2]:+.4f} "
                  f"{-wb[1]/wb[2]:+.4f}*log NU")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X3: the step statistic on its own")
    print("-" * 78)
    print(f"AUC(+step)        overall = {auc_hi(rows, lambda r: r['step']):.4f}")
    print(f"AUC(+step_block)  overall = "
          f"{auc_hi(rows, lambda r: r['stepb']):.4f}")
    print(f"AUC(-NU)          overall = {auc_lo(rows, lambda r: r['NU']):.4f}")
    print(f"AUC(-nu_hat)      overall = {auc_lo(rows, lambda r: r['nuhat']):.4f}")
    print(f"\nSpearman(step, NU)     = "
          f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):+.4f}")
    print(f"Spearman(step, nu_hat) = "
          f"{spearman([r['step'] for r in rows], [r['nuhat'] for r in rows]):+.4f}")
    print(f"Spearman(NU, nu_hat)   = "
          f"{spearman([r['NU'] for r in rows], [r['nuhat'] for r in rows]):+.4f}")

    print(f"\nstep by recovery status:")
    for tag, sel in (("success", True), ("failure", False)):
        v = sorted(r['step'] for r in rows if r['ok'] == sel)
        if v:
            print(f"  {tag:<8} N={len(v):<4} min {v[0]:+.4f}  "
                  f"median {v[len(v)//2]:+.4f}  max {v[-1]:+.4f}")
    zero = [r for r in rows if abs(r['step']) < 0.05]
    if zero:
        print(f"\n|step| < 0.05 (profile flat, W1b's 'wall crossed'): "
              f"{sum(1 for r in zero if r['ok'])}/{len(zero)} recovered "
              f"vs {sum(1 for r in rows if r['ok'])}/{len(rows)} overall")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X4: the fitted boundary is a PRODUCT LAW")
    print("-" * 78)
    print("The full-table coefficients on (log NU, log nu_hat) came out nearly")
    print("EQUAL, so the boundary collapses to  NU * nu_hat^a = const with")
    print("a ~ 1.  Test the unfitted a = 1 form directly:  PSI = NU * nu_hat.\n")
    a_fit = w2[2] / w2[1] if abs(w2[1]) > 1e-12 else float('nan')
    thr_fit = math.exp(-w2[0] / w2[1]) if abs(w2[1]) > 1e-12 else float('nan')
    print(f"  fitted exponent  a = w_nuhat / w_NU = {a_fit:.4f}   (a = 1 is the"
          " product law)")
    print(f"  fitted threshold     = {thr_fit:.4f}   (at a = a_fit)")
    print(f"  AUC(-PSI)            = {auc_lo(rows, lambda r: r['psi']):.4f}   "
          f"vs AUC(-NU) {auc_lo(rows, lambda r: r['NU']):.4f}, "
          f"AUC(-nu_hat) {auc_lo(rows, lambda r: r['nuhat']):.4f}")

    print(f"\n{'thr':>8} {'TP':>5} {'FP':>5} {'FN':>5} {'TN':>5} {'acc':>7} "
          f"{'prec':>7} {'recall':>7}")
    best = None
    for thr in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.1, 1.2, 1.5, 2.0):
        tp = sum(1 for r in rows if r['psi'] <= thr and r['ok'])
        fp = sum(1 for r in rows if r['psi'] <= thr and not r['ok'])
        fn = sum(1 for r in rows if r['psi'] > thr and r['ok'])
        tn = sum(1 for r in rows if r['psi'] > thr and not r['ok'])
        acc = (tp + tn) / len(rows)
        prec = tp / (tp + fp) if tp + fp else float('nan')
        rec = tp / (tp + fn) if tp + fn else float('nan')
        if best is None or acc > best[1]:
            best = (thr, acc)
        print(f"{thr:>8.2f} {tp:>5} {fp:>5} {fn:>5} {tn:>5} {acc:>7.4f} "
              f"{prec:>7.4f} {rec:>7.4f}")
    print(f"\nbest accuracy {best[1]:.4f} at PSI <= {best[0]:.2f}  "
          f"(NU-only best was acc 0.7800)")

    P = sorted(r['psi'] for r in rows if r['ok'])
    Nn = sorted(r['psi'] for r in rows if not r['ok'])
    print(f"\nPSI | success : min {P[0]:.4f}  median {P[len(P)//2]:.4f}  "
          f"max {P[-1]:.4f}")
    print(f"PSI | failure : min {Nn[0]:.4f}  median {Nn[len(Nn)//2]:.4f}  "
          f"max {Nn[-1]:.4f}")
    print(f"PSI bracket   : sufficient PSI < {Nn[0]:.4f} , "
          f"necessary PSI > {P[-1]:.4f}  "
          f"(overlap {P[-1]/Nn[0]:.2f}x; NU bracket overlap "
          f"{max(r['NU'] for r in rows if r['ok'])/min(r['NU'] for r in rows if not r['ok']):.2f}x)")

    # held-out seeds: the boundary must not be an in-sample artifact
    seenseeds = sorted({r['seed'] for r in rows})
    if len(seenseeds) >= 4:
        half = len(seenseeds) // 2
        tr = [r for r in rows if r['seed'] in seenseeds[:half]]
        te = [r for r in rows if r['seed'] in seenseeds[half:]]
        Xtr = [[1.0, math.log(r['NU']), math.log(r['nuhat'])] for r in tr]
        ytr = [1.0 if r['ok'] else 0.0 for r in tr]
        wtr, _ = logistic(Xtr, ytr)
        zs = [(wtr[0] + wtr[1]*math.log(r['NU']) + wtr[2]*math.log(r['nuhat']),
               r['ok']) for r in te]
        a = auc([-z for z, o in zs if o], [-z for z, o in zs if not o])
        acc = sum(1 for z, o in zs if (z > 0) == o) / len(zs)
        print(f"\nheld-out seeds  (train {len(tr)} on seeds {seenseeds[:half]},"
              f" test {len(te)}):")
        print(f"  fitted exponent a = {wtr[2]/wtr[1]:.4f}   "
              f"out-of-sample AUC {a:.4f}  acc {acc:.4f}")
        atest = auc_lo(te, lambda r: r['psi'])
        acc1 = sum(1 for r in te if (r['psi'] <= 1.0) == r['ok']) / len(te)
        print(f"  UNFITTED  PSI <= 1  on the same held-out rows: "
              f"AUC {atest:.4f}  acc {acc1:.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X5: does PSI ~ 1 survive a change of n?  (dim held at 24)")
    print("-" * 78)
    print("If PSI is a law and not a 17-bit fit, the fitted exponent and")
    print("threshold should be stable as n moves over 8 bits.\n")
    print(f"{'bits':>5} {'N':>5} {'rec':>9} {'AUC NU':>8} {'AUC nuh':>8} "
          f"{'AUC PSI':>8} | {'a_fit':>7} {'thr_fit':>8} {'acc@1':>7}")
    xsize = {}
    for bits in SIZES:
        try:
            rr = collect(args.seeds, bits=bits, verbose=False)
        except Exception as exc:                       # noqa: BLE001
            print(f"{bits:>5}  FAILED: {type(exc).__name__}: {exc}")
            continue
        if not rr:
            print(f"{bits:>5}  no instances")
            continue
        w = sum(1 for r in rr if r['ok'])
        if w == 0 or w == len(rr):
            print(f"{bits:>5} {len(rr):>5} {str(w)+'/'+str(len(rr)):>9} "
                  f"  (degenerate)")
            continue
        Xr = [[1.0, math.log(r['NU']), math.log(r['nuhat'])] for r in rr]
        yr = [1.0 if r['ok'] else 0.0 for r in rr]
        wr, _ = logistic(Xr, yr)
        af = wr[2] / wr[1] if abs(wr[1]) > 1e-12 else float('nan')
        tf = math.exp(-wr[0] / wr[1]) if abs(wr[1]) > 1e-12 else float('nan')
        acc1 = sum(1 for r in rr if (r['psi'] <= 1.0) == r['ok']) / len(rr)
        print(f"{bits:>5} {len(rr):>5} {str(w)+'/'+str(len(rr)):>9} "
              f"{auc_lo(rr, lambda r: r['NU']):>8.4f} "
              f"{auc_lo(rr, lambda r: r['nuhat']):>8.4f} "
              f"{auc_lo(rr, lambda r: r['psi']):>8.4f} | "
              f"{af:>7.4f} {tf:>8.4f} {acc1:>7.4f}")
        xsize[bits] = rr

    if len(xsize) >= 2:
        pool = [r for rr in xsize.values() for r in rr]
        print(f"\npooled across sizes (N={len(pool)}):")
        print(f"  AUC(-NU)   {auc_lo(pool, lambda r: r['NU']):.4f}    "
              f"AUC(-nu_hat) {auc_lo(pool, lambda r: r['nuhat']):.4f}    "
              f"AUC(-PSI)  {auc_lo(pool, lambda r: r['psi']):.4f}")
        acc1 = sum(1 for r in pool if (r['psi'] <= 1.0) == r['ok']) / len(pool)
        tp = sum(1 for r in pool if r['psi'] <= 1.0 and r['ok'])
        fp = sum(1 for r in pool if r['psi'] <= 1.0 and not r['ok'])
        fn = sum(1 for r in pool if r['psi'] > 1.0 and r['ok'])
        print(f"  UNFITTED PSI <= 1 : TP {tp}  FP {fp}  FN {fn}  "
              f"acc {acc1:.4f}")
        # cross-size transfer: fit on one size, score another
        print(f"\n  cross-size transfer of the fitted boundary (AUC / acc):")
        ks = sorted(xsize)
        print("    " + "train\\test".ljust(12) +
              "".join(f"{b:>16}" for b in ks))
        for btr in ks:
            Xr = [[1.0, math.log(r['NU']), math.log(r['nuhat'])]
                  for r in xsize[btr]]
            yr = [1.0 if r['ok'] else 0.0 for r in xsize[btr]]
            wr, _ = logistic(Xr, yr)
            cells = []
            for bte in ks:
                zs = [(wr[0] + wr[1]*math.log(r['NU'])
                       + wr[2]*math.log(r['nuhat']), r['ok'])
                      for r in xsize[bte]]
                p_ = [-z for z, o in zs if o]
                n_ = [-z for z, o in zs if not o]
                a = auc(p_, n_) if p_ and n_ else float('nan')
                ac = sum(1 for z, o in zs if (z > 0) == o) / len(zs)
                cells.append(f"{a:.3f}/{ac:.3f}")
            print("    " + str(btr).ljust(12) +
                  "".join(f"{c:>16}" for c in cells))

    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("EXP X6: is the threshold PSI ~ 1 a LAW or a dim-24 coincidence?")
    print("-" * 78)
    print("X5 held dim = 2m = 24 at every n, so it cannot separate 'threshold")
    print("is 1' from 'threshold is f(dim) which happens to be 1 at dim 24'.")
    print("Sweep m at fixed n (17 bits).  If thr_fit drifts with m, PSI <= 1")
    print("is a fit; if it stays at ~1, it is a law.\n")
    print(f"{'m':>4} {'dim':>5} {'N':>5} {'rec':>9} {'AUC NU':>8} "
          f"{'AUC PSI':>8} | {'a_fit':>7} {'thr_fit':>8} {'acc@1':>7}")
    global M17
    m_saved = M17
    for m in (6, 8, 10, 12, 16, 20):
        M17 = m
        try:
            rr = collect(args.seeds, bits=17, verbose=False)
        except Exception as exc:                       # noqa: BLE001
            print(f"{m:>4} {2*m:>5}  FAILED: {type(exc).__name__}: {exc}")
            continue
        if not rr:
            print(f"{m:>4} {2*m:>5}  no instances")
            continue
        w = sum(1 for r in rr if r['ok'])
        if w == 0 or w == len(rr):
            print(f"{m:>4} {2*m:>5} {len(rr):>5} "
                  f"{str(w)+'/'+str(len(rr)):>9}   (degenerate)")
            continue
        Xr = [[1.0, math.log(r['NU']), math.log(r['nuhat'])] for r in rr]
        yr = [1.0 if r['ok'] else 0.0 for r in rr]
        wr, _ = logistic(Xr, yr)
        af = wr[2] / wr[1] if abs(wr[1]) > 1e-12 else float('nan')
        tf = math.exp(-wr[0] / wr[1]) if abs(wr[1]) > 1e-12 else float('nan')
        acc1 = sum(1 for r in rr if (r['psi'] <= 1.0) == r['ok']) / len(rr)
        print(f"{m:>4} {2*m:>5} {len(rr):>5} {str(w)+'/'+str(len(rr)):>9} "
              f"{auc_lo(rr, lambda r: r['NU']):>8.4f} "
              f"{auc_lo(rr, lambda r: r['psi']):>8.4f} | "
              f"{af:>7.4f} {tf:>8.4f} {acc1:>7.4f}")
    M17 = m_saved

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
