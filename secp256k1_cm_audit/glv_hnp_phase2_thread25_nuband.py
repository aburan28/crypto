"""
GLV-HNP Phase 2, Thread 25: does mu still separate INSIDE the NU-ambiguous
band, or is its apparent power in W5 entirely mediated by NU?

H25: within the ambiguous band 1.04 <= NU <= 2.20 (Thread 24's W4 bracket,
     where the nearest-plane certificate gives no answer either way),
     AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.

Falsifier: AUC inside the band collapses towards 0.5 -> mu's W5 power was a
           NU-mediated stratification artifact.

Re-analysis only, no new lattice work: reads the 500-row table dumped by
`glv_hnp_phase2_gsprofile_strat.py --dump-json thread25_rows.json` (2026-08-08
17-bit run, 5 eff strata x 20 curves x 5 seeds).

Run: python3 glv_hnp_phase2_thread25_nuband.py thread25_rows.json
"""

import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from glv_hnp_phase2_gsprofile import auc, spearman

NU_LO, NU_HI = 1.040, 2.199  # Thread 24 W4 17-bit bracket


def logistic_fit_2d(rows, xkey, ykey, iters=4000, lr=0.05):
    """Minimal gradient-descent logistic regression on (log x, log y) ->
    ok, no external deps. Returns (w0, w1, w2) for
    p = sigmoid(w0 + w1*log(x) + w2*log(y))."""
    xs = [math.log(r[xkey]) for r in rows]
    ys = [math.log(r[ykey]) for r in rows]
    labels = [1.0 if r["ok"] else 0.0 for r in rows]
    mx, sx = sum(xs) / len(xs), (sum((v - sum(xs) / len(xs)) ** 2 for v in xs) / len(xs)) ** 0.5
    my, sy = sum(ys) / len(ys), (sum((v - sum(ys) / len(ys)) ** 2 for v in ys) / len(ys)) ** 0.5
    xs = [(v - mx) / sx for v in xs]
    ys = [(v - my) / sy for v in ys]
    w0 = w1 = w2 = 0.0
    n = len(rows)
    for _ in range(iters):
        g0 = g1 = g2 = 0.0
        for x, y, lab in zip(xs, ys, labels):
            z = w0 + w1 * x + w2 * y
            p = 1.0 / (1.0 + math.exp(-z)) if z > -30 else 0.0
            err = p - lab
            g0 += err
            g1 += err * x
            g2 += err * y
        w0 -= lr * g0 / n
        w1 -= lr * g1 / n
        w2 -= lr * g2 / n
    # unstandardize: w1*(logx - mx)/sx + w2*(logy - my)/sy + w0
    #   = (w1/sx)*logx + (w2/sy)*logy + (w0 - w1*mx/sx - w2*my/sy)
    b1, b2 = w1 / sx, w2 / sy
    b0 = w0 - w1 * mx / sx - w2 * my / sy
    return b0, b1, b2


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "thread25_rows.json"
    rows = json.load(open(path))
    print("=" * 78)
    print("Thread 25 — does mu separate inside the NU-ambiguous band?")
    print("=" * 78)
    print(f"\n{len(rows)} rows loaded from {path}")

    band = [r for r in rows if NU_LO <= r["NU"] <= NU_HI]
    below = [r for r in rows if r["NU"] < NU_LO]
    above = [r for r in rows if r["NU"] > NU_HI]
    print(f"NU band [{NU_LO}, {NU_HI}]: {len(band)} rows "
          f"({len(below)} below / sufficient, {len(above)} above / necessary-fail)")

    pos = [r for r in band if r["ok"]]
    neg = [r for r in band if not r["ok"]]
    print(f"  in-band recoveries: {len(pos)}/{len(band)}")
    if pos and neg:
        a_mu = auc([r["mu"] for r in pos], [r["mu"] for r in neg])
        a_mueff = auc([r["mu"] * math.sqrt(r["eff"]) for r in pos],
                       [r["mu"] * math.sqrt(r["eff"]) for r in neg])
        a_nh = auc([r["nuhat"] for r in pos], [r["nuhat"] for r in neg])
        a_nu = auc([r["NU"] for r in pos], [r["NU"] for r in neg])
        print(f"  AUC(-mu             -> recovery), in-band  = {a_mu:.4f}")
        print(f"  AUC(-mu*sqrt(eff)   -> recovery), in-band  = {a_mueff:.4f}")
        print(f"  AUC(-nu_hat         -> recovery), in-band  = {a_nh:.4f}")
        print(f"  AUC(-NU             -> recovery), in-band  = {a_nu:.4f}  (should be ~0.5: NU is flat here by construction)")
        print(f"\nH25 verdict (raw mu, pooled across eff): "
              f"{'HOLDS' if a_mu >= 0.8 else 'FALSIFIED'} (threshold 0.8, observed {a_mu:.4f})")
    else:
        print("  degenerate band (all-pos or all-neg), cannot compute AUC")

    print("\n" + "-" * 78)
    print("Per-eff-stratum in-band AUC(-mu), to check H25 isn't itself an")
    print("eff-mediated artifact (band membership correlates with eff)")
    print("-" * 78)
    effs = sorted(set(r["effq"] for r in rows))
    print(f"{'eff':>5} {'N_band':>7} {'rec':>7} {'AUC mu':>8}")
    for eff in effs:
        sub = [r for r in band if r["effq"] == eff]
        p = [r for r in sub if r["ok"]]
        ng = [r for r in sub if not r["ok"]]
        if p and ng:
            print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f}")
        else:
            print(f"{eff:>5.2f} {len(sub):>7} {str(len(p))+'/'+str(len(sub)):>7} {'(degenerate)':>8}")

    print("\n" + "-" * 78)
    print("Logistic fit on (log NU, log X) -> recovery, full 500-row table")
    print("-" * 78)
    for xkey, label in (("mu", "mu (raw, not eff-normalized)"),
                         ("nuhat", "nu_hat (eff-normalized)")):
        b0, bNU, bx = logistic_fit_2d(rows, "NU", xkey)
        print(f"  X = {label}")
        print(f"  logit(p) = {b0:.4f} + {bNU:.4f}*log(NU) + {bx:.4f}*log(X)")

        def score(r, b0=b0, bNU=bNU, bx=bx, xkey=xkey):
            return b0 + bNU * math.log(r["NU"]) + bx * math.log(r[xkey])
        sp = [score(r) for r in rows if r["ok"]]
        sn = [score(r) for r in rows if not r["ok"]]
        # score is a logit of P(recovery), so LARGER score -> more likely to
        # recover; auc() expects "smaller -> recovery", so negate.
        a_fit = auc([-s for s in sp], [-s for s in sn])
        print(f"  in-sample AUC of fitted (NU, {xkey}) score = {a_fit:.4f}\n")
    print(f"  (pooled single-predictor AUCs from Thread 24 W5: "
          f"nu_hat*sqrt(eff) 0.9348, NU 0.7996, nu_hat 0.6889, raw mu untested there)")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
