"""
GLV-HNP Phase 2, Thread 25: does mu (or nu_hat) separate recovery INSIDE the
NU ambiguous band, i.e. is there a genuine second coordinate beyond NU?

H25: within the ambiguous band 1.04 <= NU <= 2.20 (17-bit bracket, Thread 24
     W4), AUC(-mu -> Kannan-LLL recovery) stays >= 0.8.
Falsifier: if AUC drops back to ~0.5 inside the band, mu's apparent power in
     W5 is entirely mediated by NU (a stratification artifact) and the
     closed form should be retired.

Secondary: step = log2(profm) - log2(prof0), the GS-profile head/tail jump
     identified in Thread 24 W1b.  Does step -> 0 predict the wall better
     than NU or mu, inside the same band?

Reads secp256k1_cm_audit/thread25_rows.json, written by
glv_hnp_phase2_gsprofile_strat.py --dump-json (500 rows, 17-bit, dim 24,
5 eff strata x 20 curves x 5 seeds).
"""

import json
import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from glv_hnp_phase2_gsprofile import auc, spearman

with open("thread25_rows.json") as f:
    rows = json.load(f)

for r in rows:
    r["step"] = math.log2(r["profm"]) - math.log2(r["prof0"])

NU_LO, NU_HI = 1.040, 2.199  # Thread 24 W4, 17-bit bracket

print("=" * 78)
print("Thread 25 — is mu a genuine second coordinate inside the NU band?")
print("=" * 78)

band = [r for r in rows if NU_LO <= r["NU"] <= NU_HI]
below = [r for r in rows if r["NU"] < NU_LO]
above = [r for r in rows if r["NU"] > NU_HI]
print(f"\nN={len(rows)} total: below bracket {len(below)}, "
      f"in band [{NU_LO},{NU_HI}] {len(band)}, above {len(above)}")
print(f"below-bracket recovery rate: "
      f"{sum(r['ok'] for r in below)}/{len(below)}")
print(f"above-bracket recovery rate: "
      f"{sum(r['ok'] for r in above)}/{len(above)}")
print(f"in-band recovery rate:       "
      f"{sum(r['ok'] for r in band)}/{len(band)}")

pos = [r for r in band if r["ok"]]
neg = [r for r in band if not r["ok"]]
print(f"\nin-band: {len(pos)} recovered, {len(neg)} failed")

if pos and neg:
    a_mu = auc([r["mu"] for r in pos], [r["mu"] for r in neg])
    a_nh = auc([r["nuhat"] for r in pos], [r["nuhat"] for r in neg])
    a_nu = auc([r["NU"] for r in pos], [r["NU"] for r in neg])
    a_st = auc([r["step"] for r in pos], [r["step"] for r in neg])
    a_ls = auc([r["lamstar"] for r in pos], [r["lamstar"] for r in neg])
    print(f"AUC(-mu     -> recovery), in-band  = {a_mu:.4f}")
    print(f"AUC(-nu_hat -> recovery), in-band  = {a_nh:.4f}")
    print(f"AUC(-NU     -> recovery), in-band  = {a_nu:.4f}  "
          f"(expect ~0.5: NU is what defines the band)")
    print(f"AUC(-step   -> recovery), in-band  = {a_st:.4f}")
    print(f"AUC(-lam*   -> recovery), in-band  = {a_ls:.4f}  (control)")
    print(f"\nH25 verdict: AUC(mu) {'>= 0.8 -- HOLDS' if a_mu >= 0.8 else '< 0.8 -- FALSIFIED'}")
else:
    print("degenerate band (all-pos or all-neg) -- cannot test H25 here")

# Per-eff-stratum breakdown of the band, since eff also varies inside it.
print("\n" + "-" * 78)
print("in-band AUC(mu), split by eff stratum (mu's power might be eff-driven)")
print("-" * 78)
for eff in sorted(set(r["effq"] for r in rows)):
    sub = [r for r in band if r["effq"] == eff]
    p = [r for r in sub if r["ok"]]
    ng = [r for r in sub if not r["ok"]]
    if not p or not ng:
        print(f"eff={eff:.2f}  N={len(sub):3d}  degenerate "
              f"({len(p)} rec / {len(sub)})")
        continue
    print(f"eff={eff:.2f}  N={len(sub):3d}  rec={len(p)}/{len(sub)}  "
          f"AUC(mu)={auc([r['mu'] for r in p], [r['mu'] for r in ng]):.4f}  "
          f"AUC(NU)={auc([r['NU'] for r in p], [r['NU'] for r in ng]):.4f}")

# Logistic-style boundary sanity check: does mu separate even after
# conditioning on eff via Spearman partial-like split (median eff)?
print("\n" + "-" * 78)
print("secondary: does 'step' (head/tail GS jump) predict recovery overall?")
print("-" * 78)
pos_all = [r for r in rows if r["ok"]]
neg_all = [r for r in rows if not r["ok"]]
print(f"AUC(-step -> recovery), pooled N={len(rows)}: "
      f"{auc([r['step'] for r in pos_all], [r['step'] for r in neg_all]):.4f}")
print(f"AUC(-NU   -> recovery), pooled: "
      f"{auc([r['NU'] for r in pos_all], [r['NU'] for r in neg_all]):.4f}")
print(f"AUC(-mu   -> recovery), pooled: "
      f"{auc([r['mu'] for r in pos_all], [r['mu'] for r in neg_all]):.4f}")
print(f"Spearman(step, NU) pooled = "
      f"{spearman([r['step'] for r in rows], [r['NU'] for r in rows]):.4f}")
print(f"Spearman(step, mu) pooled = "
      f"{spearman([r['step'] for r in rows], [r['mu'] for r in rows]):.4f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
