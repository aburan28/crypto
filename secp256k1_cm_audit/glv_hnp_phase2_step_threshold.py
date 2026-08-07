"""
GLV-HNP Phase 2, Thread 25b: is `step` an absolute THRESHOLD or only a rank
statistic?

W9 (glv_hnp_phase2_conditional.py) found that in every non-degenerate eff
stratum the mean GS step splits by sign:  step|recovered > 0 > step|failed.
Thread 24's W7 showed nu_hat has no clean cut ("a rank statistic, not a
threshold").  If sign(step) is a cut, it is strictly more useful than nu_hat:
it needs no calibration and no comparison set.

Reads the table dumped by glv_hnp_phase2_conditional.py --dump-json; no new
lattice work.

Run: python3 glv_hnp_phase2_step_threshold.py [rows.json]
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EFFS = (0.05, 0.10, 0.15, 0.20, 0.25)
DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "glv_hnp_phase2_conditional_rows.json")


def confusion(rows, pred):
    tp = sum(1 for r in rows if pred(r) and r['ok'])
    fp = sum(1 for r in rows if pred(r) and not r['ok'])
    fn = sum(1 for r in rows if not pred(r) and r['ok'])
    tn = sum(1 for r in rows if not pred(r) and not r['ok'])
    return tp, fp, fn, tn


def youden(rows, key, sign=1.0):
    """Best threshold t for the rule sign*key <= t, and its Youden J."""
    vals = sorted({sign * r[key] for r in rows})
    best = (-2.0, None)
    for t in vals:
        tp, fp, fn, tn = confusion(rows, lambda r, t=t: sign * r[key] <= t)
        if (tp + fn) == 0 or (fp + tn) == 0:
            continue
        j = tp / (tp + fn) - fp / (fp + tn)
        if j > best[0]:
            best = (j, t)
    return best[1], best[0]


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    rows = json.load(open(path))
    print("=" * 78)
    print("Thread 25b — is sign(step) an absolute cut?")
    print("=" * 78)
    print(f"\n{len(rows)} rows from {os.path.basename(path)}")

    print("\n" + "-" * 78)
    print("EXP W11: the fixed rule  `recover iff step > 0`, per stratum.")
    print("-" * 78)
    print(f"\n{'eff':>5} {'N':>5} {'rec':>8} | {'TP':>4} {'FP':>4} {'FN':>4} "
          f"{'TN':>4} | {'acc':>6} {'TPR':>6} {'FPR':>6}")
    for eff in list(EFFS) + [None]:
        sub = rows if eff is None else [r for r in rows if r['effq'] == eff]
        tp, fp, fn, tn = confusion(sub, lambda r: r['step'] > 0)
        acc = (tp + tn) / len(sub)
        tpr = tp / (tp + fn) if (tp + fn) else float('nan')
        fpr = fp / (fp + tn) if (fp + tn) else float('nan')
        lab = "ALL" if eff is None else f"{eff:.2f}"
        nrec = sum(1 for r in sub if r['ok'])
        print(f"{lab:>5} {len(sub):>5} {str(nrec)+'/'+str(len(sub)):>8} | "
              f"{tp:>4} {fp:>4} {fn:>4} {tn:>4} | {acc:>6.3f} "
              f"{tpr:>6.3f} {fpr:>6.3f}")

    print("\n" + "-" * 78)
    print("EXP W12: is the OPTIMAL cut stratum-stable?  (Youden-J threshold")
    print("         per stratum for step, nu_hat and NU.)")
    print("-" * 78)
    print(f"\n{'eff':>5} | {'t*(step)':>9} {'J':>6} | {'t*(nuhat)':>10} "
          f"{'J':>6} | {'t*(NU)':>8} {'J':>6}")
    for eff in EFFS:
        sub = [r for r in rows if r['effq'] == eff]
        ts, js = youden(sub, 'step', -1.0)
        th, jh = youden(sub, 'nuhat')
        tn_, jn = youden(sub, 'NU')
        print(f"{eff:>5.2f} | {-ts:>9.4f} {js:>6.3f} | {th:>10.4f} {jh:>6.3f} "
              f"| {tn_:>8.4f} {jn:>6.3f}")
    print("\n(t*(step) printed as the step value itself: recover iff "
          "step >= t*.)")

    print("\n" + "-" * 78)
    print("EXP W13: audit of the W1b head-block identity  ||b*_1|| == mu.")
    print("         If it held exactly, step would BE log2(||b*_{m+1}||/mu)")
    print("         and could not be an independent coordinate.")
    print("-" * 78)
    if rows[0].get('prof0') is None:
        print("\nprof0 not in the dump; re-run the parent script.")
    else:
        def stats(v):
            v = sorted(v)
            return v[0], v[len(v) // 2], v[-1]

        print(f"\n{'eff':>5} | {'min':>8} {'median':>8} {'max':>8} "
              f"{'exact (1e-9)':>13}")
        for eff in list(EFFS) + [None]:
            sub = rows if eff is None else [r for r in rows
                                            if r['effq'] == eff]
            v = [r['prof0'] / r['mu'] for r in sub if r['mu'] > 0]
            lo, md, hi = stats(v)
            ex = sum(1 for x in v if abs(x - 1.0) < 1e-9)
            lab = "ALL" if eff is None else f"{eff:.2f}"
            print(f"{lab:>5} | {lo:>8.4f} {md:>8.4f} {hi:>8.4f} "
                  f"{str(ex)+'/'+str(len(v)):>13}")
        v = [r['profm'] / r['l2'] for r in rows if r['l2'] > 0]
        lo, md, hi = stats(v)
        print(f"\n||b*_{{m+1}}|| / lambda_2(L2): min {lo:.4f} median {md:.4f} "
              f"max {hi:.4f}")
        print("(Thread 24 W2 measured the same ratio at the LAST index: "
              "mean 0.419.)")
    print("\nSpearman(step, mu) = -0.6627 (parent, W9): step is largely a")
    print("reparametrisation of mu, not an independent coordinate.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
