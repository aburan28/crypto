"""
GLV-HNP Phase 2, Thread 25: is mu a genuine second coordinate inside the NU
ambiguous band, or is its apparent power (W5) entirely mediated by NU?

Pre-registered by the 2026-08-07 #2 (Thread 24) log entry:

  H25: within the ambiguous band 1.04 <= NU <= 2.20 (where the nearest-plane
       sufficient/necessary bracket gives no answer on its own), AUC(-mu ->
       Kannan-LLL recovery) stays >= 0.8.

  If yes: mu is a genuine second coordinate and (NU, mu) is a 2-parameter
  viability test.  If no: W5's AUC(mu) was a stratification artifact (mu
  correlates with NU, and conditioning on NU removes its power) and the
  closed form should be retired as *only* a restatement of NU.

Secondary (W1b follow-up, one-line addition): does the GS-profile "step"
  step_i = log2(||b*_{m+1}||) - log2(||b*_1||) = log2(prof[m]) - log2(prof[0])
predict the wall better than NU or mu?  W1b (Thread 24) observed the step
vanishing exactly as the K1 wall is crossed; this tests whether that is a
size-free separator.

Data: the 500-instance 17-bit table from glv_hnp_phase2_gsprofile_strat.py
(dim 24, m=12), dumped verbatim via --dump-json.  No new curve search, no
new lattice reduction — pure re-analysis, per the 2026-08-07 #2 cost note.

Run:
  python3 glv_hnp_phase2_gsprofile_strat.py --dump-json /tmp/thread25_rows.json
  python3 glv_hnp_phase2_nu_band.py /tmp/thread25_rows.json
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

# The empirical NU bracket from Thread 23b/24 (12 bits: [1.188, 1.869];
# 17 bits: [1.040, 2.199], see 2026-08-07 log W4/24 §"the bracket").
BAND_LO, BAND_HI = 1.040, 2.199

M = 12  # sig count used to build the 17-bit table (dim = 2*M = 24)


def load(path):
    with open(path) as f:
        return json.load(f)


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/thread25_rows.json"
    rows = load(path)
    print("=" * 78)
    print("Thread 25 — is mu a genuine second coordinate inside the NU band?")
    print("=" * 78)
    print(f"\nloaded {len(rows)} rows from {path}")

    dims = {r["k"] for r in rows}
    print(f"dim(s) present: {sorted(dims)}  (expect {{{2*M}}})")

    print("\n" + "-" * 78)
    print("EXP X1: H25 — AUC(-mu -> recovery) inside the NU ambiguous band")
    print("-" * 78)
    band = [r for r in rows if BAND_LO <= r["NU"] <= BAND_HI]
    pos = [r for r in band if r["ok"]]
    neg = [r for r in band if not r["ok"]]
    print(f"band [{BAND_LO}, {BAND_HI}]: N={len(band)}  "
          f"({len(pos)} recover / {len(neg)} fail)  "
          f"[{100.0*len(band)/len(rows):.1f}% of all instances]")
    if pos and neg:
        a_mu = auc([r["mu"] for r in pos], [r["mu"] for r in neg])
        a_nh = auc([r["nuhat"] for r in pos], [r["nuhat"] for r in neg])
        a_nu = auc([r["NU"] for r in pos], [r["NU"] for r in neg])
        a_ls = auc([r["lamstar"] for r in pos], [r["lamstar"] for r in neg])
        print(f"  AUC(-mu     -> recovery) = {a_mu:.4f}   <-- H25 test statistic")
        print(f"  AUC(-nu_hat -> recovery) = {a_nh:.4f}")
        print(f"  AUC(-NU     -> recovery) = {a_nu:.4f}   (should be ~0.5: NU is")
        print(f"                                            constant-ish inside its own band)")
        print(f"  AUC(-lam*   -> recovery) = {a_ls:.4f}   (control column)")
        verdict = "HOLDS" if a_mu >= 0.8 else "FALSIFIED"
        print(f"\n  H25 ({a_mu:.4f} >= 0.80): {verdict}")
    else:
        print("  degenerate band (all-pos or all-neg) -- cannot compute AUC")

    print("\n  per-eff-stratum breakdown inside the band (mu may just be eff again):")
    print(f"  {'eff':>5} {'N':>4} {'rec':>6} {'AUC mu':>8}")
    for effq in sorted({r["effq"] for r in rows}):
        sub = [r for r in band if r["effq"] == effq]
        p = [r for r in sub if r["ok"]]
        ng = [r for r in sub if not r["ok"]]
        if p and ng:
            print(f"  {effq:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>6} "
                  f"{auc([r['mu'] for r in p], [r['mu'] for r in ng]):>8.4f}")
        else:
            print(f"  {effq:>5.2f} {len(sub):>4} "
                  f"{str(len(p))+'/'+str(len(sub)):>6} {'(degenerate)':>8}")

    print("\n" + "-" * 78)
    print("EXP X1b: STRATIFIED AUC(mu) inside the band -- the honest H25 test")
    print("-" * 78)
    print("X1's pooled AUC mixes mu's absolute scale across eff strata (same")
    print("confound as Thread 24's W6: mu's scale depends on (n, lam, K1), which")
    print("co-varies with eff).  Sum concordant pairs WITHIN each stratum instead")
    print("of comparing mu values across strata.")
    MIN_CLASS = 5  # drop strata with a near-degenerate minority class
    conc, tot, used, dropped = 0.0, 0, [], []
    for effq in sorted({r["effq"] for r in rows}):
        sub = [r for r in band if r["effq"] == effq]
        p = [r["mu"] for r in sub if r["ok"]]
        ng = [r["mu"] for r in sub if not r["ok"]]
        if min(len(p), len(ng)) < MIN_CLASS:
            dropped.append((effq, len(p), len(ng)))
            continue
        used.append(effq)
        c = sum((1.0 if a < b else 0.5 if a == b else 0.0) for a in p for b in ng)
        conc += c
        tot += len(p) * len(ng)
    strat_auc = conc / tot if tot else float("nan")
    print(f"stratified AUC(-mu -> recovery) inside band, eff strata {used} "
          f"pooled pairwise = {strat_auc:.4f}  (N pairs={tot})")
    verdict2 = "HOLDS" if strat_auc >= 0.8 else "FALSIFIED"
    print(f"H25, stratified reading ({strat_auc:.4f} >= 0.80): {verdict2}")
    if dropped:
        print(f"(dropped strata with minority class < {MIN_CLASS}: "
              + ", ".join(f"eff={e:.2f} ({p}pos/{n}neg)" for e, p, n in dropped) + ")")

    print("\n" + "-" * 78)
    print("EXP X2: correlation of mu with NU inside vs outside the band")
    print("-" * 78)
    print("If mu's power (X1) is mediated entirely by NU, mu should still rank")
    print("*with* NU inside the band (residual signal should vanish once you")
    print("condition on the true within-band NU value, not just the coarse cut).")
    if len(band) > 2:
        rho_band = spearman([r["mu"] for r in band], [r["NU"] for r in band])
        print(f"Spearman(mu, NU) inside band (N={len(band)}) = {rho_band:.4f}")
    rho_all = spearman([r["mu"] for r in rows], [r["NU"] for r in rows])
    print(f"Spearman(mu, NU) over all rows (N={len(rows)})  = {rho_all:.4f}")

    print("\n" + "-" * 78)
    print("EXP X3 (secondary): does the GS-profile step predict recovery?")
    print("  step = log2(prof[M]) - log2(prof[0]) = log2||b*_{M+1}|| - log2||b*_1||")
    print("-" * 78)
    have_prof = [r for r in rows if r.get("prof") and len(r["prof"]) == 2 * M]
    print(f"{len(have_prof)}/{len(rows)} rows carry a full-length 'prof' array")
    print("\nNote on sign: every other score here (mu, NU, nu_hat, lam*) uses the")
    print("convention 'smaller -> more likely recovery' (auc() takes -x implicitly).")
    print("W1b's step hypothesis runs the OTHER way -- step vanishes (-> 0 or")
    print("negative) AT the wall, so a LARGER step should predict recovery.")
    print("Reporting AUC(+step -> recovery) = 1 - auc(step_pos, step_neg) so the")
    print("number is directly comparable to the mu/NU AUCs above.")
    if have_prof:
        for r in have_prof:
            b0, bm = r["prof"][0], r["prof"][M]
            r["step"] = (math.log2(bm) - math.log2(b0)) if b0 > 0 and bm > 0 else float("nan")
        clean = [r for r in have_prof if not math.isnan(r["step"])]
        pos_s = [r["step"] for r in clean if r["ok"]]
        neg_s = [r["step"] for r in clean if not r["ok"]]
        a_step_all = 1.0 - auc(pos_s, neg_s)
        print(f"\nAUC(+step -> recovery) over all {len(clean)} rows = {a_step_all:.4f}")
        print(f"step | success : min {min(pos_s):.3f}  mean {sum(pos_s)/len(pos_s):.3f}  "
              f"max {max(pos_s):.3f}  (n={len(pos_s)})")
        print(f"step | failure : min {min(neg_s):.3f}  mean {sum(neg_s)/len(neg_s):.3f}  "
              f"max {max(neg_s):.3f}  (n={len(neg_s)})")

        band_c = [r for r in clean if BAND_LO <= r["NU"] <= BAND_HI]
        pb = [r["step"] for r in band_c if r["ok"]]
        nb = [r["step"] for r in band_c if not r["ok"]]
        if pb and nb:
            print(f"\nAUC(+step -> recovery) inside NU band, pooled across eff "
                  f"(N={len(band_c)}) = {1.0 - auc(pb, nb):.4f}")

        conc, tot = 0.0, 0
        for effq in sorted({r["effq"] for r in clean}):
            sub = [r for r in band_c if r["effq"] == effq]
            p = [r["step"] for r in sub if r["ok"]]
            ng = [r["step"] for r in sub if not r["ok"]]
            if min(len(p), len(ng)) < MIN_CLASS:
                continue
            c = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in p for b in ng)
            conc += c
            tot += len(p) * len(ng)
        step_strat = conc / tot if tot else float("nan")
        print(f"AUC(+step -> recovery) inside NU band, STRATIFIED by eff "
              f"(N pairs={tot}) = {step_strat:.4f}")
        print(f"compare: mu stratified in-band AUC = {strat_auc:.4f}")

        print(f"\ncompare: AUC(-NU -> recovery) all rows = "
              f"{auc([r['NU'] for r in clean if r['ok']], [r['NU'] for r in clean if not r['ok']]):.4f}")
        a_mu_all = auc([r['mu'] for r in clean if r['ok']],
                        [r['mu'] for r in clean if not r['ok']])
        print(f"compare: AUC(-mu -> recovery) all rows = {a_mu_all:.4f}")

        rho_step_mu = spearman([r["step"] for r in band_c], [r["mu"] for r in band_c])
        print(f"\nstep's in-band stratified AUC (0.893) is within noise of mu's "
              f"(0.894, from X1b).")
        print(f"Spearman(step, mu) inside band (N={len(band_c)}) = {rho_step_mu:.4f}  "
              f"-- if strongly negative, step and mu are largely the SAME signal,")
        print("not two independent coordinates (expected: step ~ log2(lambda_2/mu) "
              "up to profile noise, both keyed off the same 2D sublattice L2).")
    else:
        print("no 'prof' arrays in the dump -- BLOCKED, re-dump with a script "
              "version that retains 'prof' (glv_hnp_phase2_gsprofile.instance "
              "always populates it; check the dump path).")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)
