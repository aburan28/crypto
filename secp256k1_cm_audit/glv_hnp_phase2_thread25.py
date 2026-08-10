"""
GLV-HNP Phase 2, Thread 25: does mu survive as a second coordinate once NU
is held fixed, and does the block-transition "step" in the GS profile beat
both?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry, "Next step
proposal":

  H25: within the ambiguous NU band [1.040, 2.199] (17-bit bracket from
       Thread 23b/W4 — nearest-plane gives no answer there), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate and (NU, mu) is a 2-parameter
  viability test.  If no: W5's apparent mu power is entirely mediated by NU
  and the closed form should be retired.

Secondary (W1b-suggested): step = log2(||b*_{m+1}||) - log2(||b*_m||), the
jump from the lambda_1(L2)-repeated head block to the K1-moving tail block
of the L0 Gram-Schmidt profile (0-indexed: prof[0..m-1] is the head,
prof[m..2m-1] the tail; m=12, dim=24 here).  Test whether step -> 0 predicts
the wall better than NU or mu.

Data: reuses the 500-instance 17-bit eff-fixed table dumped by
`glv_hnp_phase2_gsprofile_strat.py --dump-json` (no new curve search, no new
lattice work — this script is pure re-analysis of existing rows, each of
which already carries 'prof', the full ||b*_i|| profile of L0).

Run: python3 glv_hnp_phase2_gsprofile_strat.py --dump-json   (if the json is
     stale or missing)
     python3 glv_hnp_phase2_thread25.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "glv_hnp_phase2_gsprofile_strat_data.json")

M = 12  # matches M17 in glv_hnp_phase2_gsprofile_strat.py
NU_BAND = (1.040, 2.199)  # 17-bit bracket, Thread 23 2026-08-07 log entry


def logistic_auc_ci_note(auc_val, npos, nneg):
    """Rough Hanley-McNeil SE for context, not a formal test."""
    q1 = auc_val / (2 - auc_val)
    q2 = 2 * auc_val ** 2 / (1 + auc_val)
    var = (auc_val * (1 - auc_val) + (npos - 1) * (q1 - auc_val ** 2)
           + (nneg - 1) * (q2 - auc_val ** 2)) / (npos * nneg)
    return math.sqrt(max(var, 0.0))


if __name__ == "__main__":
    if not os.path.exists(DATA_PATH):
        print(f"missing {DATA_PATH} -- run glv_hnp_phase2_gsprofile_strat.py "
              "--dump-json first")
        sys.exit(1)

    rows = json.load(open(DATA_PATH))
    for r in rows:
        r['step'] = math.log2(r['prof'][M]) - math.log2(r['prof'][M - 1])
    print("=" * 78)
    print(f"Thread 25 — mu vs NU inside the ambiguous band; GS-profile step")
    print(f"loaded {len(rows)} instances from {os.path.basename(DATA_PATH)}")
    print("=" * 78)

    print("\n" + "-" * 78)
    print(f"H25 primary: within NU in {NU_BAND}, does mu still separate?")
    print("-" * 78)
    band = [r for r in rows if NU_BAND[0] <= r['NU'] <= NU_BAND[1]]
    pos = [r for r in band if r['ok']]
    neg = [r for r in band if not r['ok']]
    print(f"band N = {len(band)} ({len(pos)} recover / {len(neg)} fail) "
          f"out of {len(rows)} total ({100*len(band)/len(rows):.1f}%)")
    if pos and neg:
        a_mu = auc([r['mu'] for r in pos], [r['mu'] for r in neg])
        a_nh = auc([r['nuhat'] for r in pos], [r['nuhat'] for r in neg])
        a_step_neg = auc([r['step'] for r in pos], [r['step'] for r in neg])
        a_step_pos = 1.0 - a_step_neg  # step's natural sign is "bigger=better"
        se_mu = logistic_auc_ci_note(a_mu, len(pos), len(neg))
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}  (SE~{se_mu:.3f})")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(+step   -> recovery) = {a_step_pos:.4f}   "
              f"(step = log2||b*_m|| - log2||b*_{{m-1}}||; W1b's finding is "
              f"that step -> 0 as the wall is crossed, so recovery should "
              f"correlate with a BIGGER step, not smaller)")
        verdict = "SURVIVES" if a_mu >= 0.8 else "DOES NOT SURVIVE"
        print(f"\n  H25 verdict: mu {verdict} inside the ambiguous band "
              f"(threshold 0.8, got {a_mu:.4f}).")
        print(f"  Side finding: nu_hat (mu normalised by sqrt(det L2)) "
              f"{'SURVIVES' if a_nh >= 0.8 else 'does not survive'} "
              f"the same band test ({a_nh:.4f}), unlike raw mu.")
    else:
        print("  degenerate band (all-one-class) -- cannot compute AUC")

    print("\n" + "-" * 78)
    print("Same test at coarser bands, for robustness to the exact cutoff")
    print("-" * 78)
    for lo, hi in [(1.0, 2.5), (1.1, 2.0), (1.2, 1.9), (0.9, 3.0)]:
        sub = [r for r in rows if lo <= r['NU'] <= hi]
        p = [r for r in sub if r['ok']]
        n = [r for r in sub if not r['ok']]
        if p and n:
            a_mu = auc([r['mu'] for r in p], [r['mu'] for r in n])
            a_nh = auc([r['nuhat'] for r in p], [r['nuhat'] for r in n])
            a_stp = 1.0 - auc([r['step'] for r in p], [r['step'] for r in n])
            print(f"  NU in [{lo:.1f},{hi:.1f}]  N={len(sub):>3}  "
                  f"AUC(-mu)={a_mu:.4f}  AUC(-nu_hat)={a_nh:.4f}  "
                  f"AUC(+step)={a_stp:.4f}")
        else:
            print(f"  NU in [{lo:.1f},{hi:.1f}]  N={len(sub):>3}  degenerate")

    print("\n" + "-" * 78)
    print("Pooled (all 500 instances): step vs NU vs mu vs nu_hat")
    print("-" * 78)
    P = [r for r in rows if r['ok']]
    N = [r for r in rows if not r['ok']]
    print(f"  AUC(-NU)      = {auc([r['NU'] for r in P], [r['NU'] for r in N]):.4f}")
    print(f"  AUC(-mu)      = {auc([r['mu'] for r in P], [r['mu'] for r in N]):.4f}")
    print(f"  AUC(-nu_hat)  = {auc([r['nuhat'] for r in P], [r['nuhat'] for r in N]):.4f}")
    print(f"  AUC(+step)    = {1.0 - auc([r['step'] for r in P], [r['step'] for r in N]):.4f}")
    print(f"  Spearman(step, NU)  = {spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
    print(f"  Spearman(step, mu)  = {spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

    print("\n" + "-" * 78)
    print("Per-stratum (eff fixed) AUC(+step), for comparison with the W5 table")
    print("-" * 78)
    print(f"{'eff':>5} {'N':>5} {'rec':>7} | {'AUC +step':>9} {'AUC mu':>8} "
          f"{'AUC NU':>8}")
    for eff in sorted(set(r['effq'] for r in rows)):
        sub = [r for r in rows if r['effq'] == eff]
        p = [r for r in sub if r['ok']]
        n = [r for r in sub if not r['ok']]
        if not p or not n:
            print(f"{eff:>5.2f} {len(sub):>5} "
                  f"{str(len(p))+'/'+str(len(sub)):>7} | (degenerate)")
            continue
        print(f"{eff:>5.2f} {len(sub):>5} {str(len(p))+'/'+str(len(sub)):>7} | "
              f"{1.0 - auc([r['step'] for r in p], [r['step'] for r in n]):>9.4f} "
              f"{auc([r['mu'] for r in p], [r['mu'] for r in n]):>8.4f} "
              f"{auc([r['NU'] for r in p], [r['NU'] for r in n]):>8.4f}")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
