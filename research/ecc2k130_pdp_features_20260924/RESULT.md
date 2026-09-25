# Affine prescreen on toy two-summand PDPs

Status: **toy diagnostic**, 2026-09-24. The original 716-case producer was rerun against main source commit f0d4aceb and reproduced every base, target, feature, root count, group label, and check from the prior provisional artifact after timing fields were excluded. The saved producer_results.json is the fresh main run: 640 natural cases, 76 planted controls, and zero producer correctness disagreements. The analysis recomputed all 640 natural profiles and complete fiber-root counts, then recomputed the group labels for all 55 double-held-out cases.

The safe rule is exact for the tested S3 equation: an affine contradiction from the original equation row span cannot coexist with a root. It rejected 88 of 640 natural cases and retained all 111 natural group hits. On the confirmation intersection of held-out base and held-out target Frobenius orbit, it rejected 4 of 55 cases and retained all 13 hits. Group labels from different bases on the same target orbit were kept together for splitting.

The table reports **screen plus remaining complete fiber solves / fiber solves on every case**. All feature costs are charged. Each native counter is shown separately because a field multiplication, square, and row XOR have no measured conversion to a common curve-operation unit.

| Split | Cases | Hits retained | Certified skips | Field mul ratio | Field sqr ratio | Row XOR ratio |
|:--|--:|--:|--:|--:|--:|--:|
| Training intersection | 292 | 47/47 | 35 | 0.9457 | 0.9519 | 1.0069 |
| Held base only | 73 | 18/18 | 11 | 0.8507 | 0.8572 | 0.8838 |
| Held orbit only | 220 | 33/33 | 38 | 0.9226 | 0.9310 | 0.9909 |
| Double holdout | 55 | 13/13 | 4 | 0.9415 | 0.9504 | 0.9806 |

CPU timing sums were retained with three paired per-case repetitions. The confirmation ratio in the saved run is about 0.956; two prior local replays gave about 0.895 and 0.894. These small process-time differences are sensitive to run conditions and are not an end-to-end speed measurement. The deterministic counter ratios above are the repeatable observation.

The rank priority hypothesis is **not supported**. Among the 51 confirmation cases that survive the affine screen, the target lying in the product span has 6/21 group hits, versus 7/30 outside it. The aggregate difference reverses at n=7, and the n=13 and n=17 confirmation cells have no hits at all. A rule that discards high-rank cases would also lose real hits. The rank identity predicts the quadratic coefficient rank for nonzero targets, but does not predict solver hardness or useful relation yield.

This result supports one narrow implementation: use an affine contradiction as a certified early exit, with its setup charged. It does not select factor bases or prove a degree-131 win. The smallest next test is to cache the product span per base, then apply the same screen to naturally sampled m=3 or chained PDPs on implicit orbit bases while counting the first group-verified relation and independent row gain. Preserve failed targets and compare to the same baseline and rho accounting contract.

Reproduce from a checkout of main containing the stated source hash, using Python 3.10 or newer:

    python3 research/ecc2k130_pdp_features_20260924/profile_probe.py --source-repo . --out /tmp/pdp-profile-rerun
    python3 research/ecc2k130_pdp_features_20260924/analyse.py --source-repo . --input /tmp/pdp-profile-rerun/results.json --out /tmp/pdp-holdout-rerun.json

The commands refuse existing output paths. CONTRACT.md freezes the split and acceptance rules. producer_results.json and holdout_results.json contain the complete saved cases and counters. The original producer contract is pdp_profile_contract.md. Source and evidence hashes are embedded in the JSON.
