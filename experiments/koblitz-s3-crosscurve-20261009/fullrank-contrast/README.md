# Full-rank S3 cross-curve pilot contrast

The frozen N41 C800 and N53 C1200 full-rank S3 pipelines each showed an
exploratory paired online gain above 2x. Their relative gains are not yet
resolved: the N53/N41 ratio of geometric gains is **1.447**, with a
hierarchical target-and-repeat bootstrap 95% interval **[0.591, 3.722]**.
This compares two complete candidate/baseline pairs. The factor bases,
baseline ranks, target laws, and curve sizes differ, so the contrast does not
isolate field degree as a cause.

The [reproducible analyzer](analyze.py) verifies the original freeze, raw
ledger, analysis, and independent Sage replay hashes for both experiments.
Its [receipt](analysis.json) names all candidate IDs and input digests. It
resamples independent target blocks separately on each curve and resamples
the two repeated online timings within each arm. N41 contributes 24 frozen
targets T033–T056 and 96 verified native and Sage solves; N53 contributes
eight frozen targets T025–T032 and 32 verified native and Sage solves.

| Paired one-target online measure | N41 C244 / C800 | N53 C396 / C1200 |
| --- | ---: | ---: |
| Geometric baseline/full-rank gain | 3.514x | 5.085x |
| Hierarchical bootstrap 95% interval | [2.615, 4.663] | [2.175, 12.156] |
| A/A p95 absolute log spread, baseline / full rank | 0.885 / 0.889 | 0.958 / 1.381 |

With the observed mean log contrast and variance, a conventional two-sided
5% test with 80% power would require approximately **123 independent targets
per curve**, with two repetitions per arm: **984 one-target solves** in
total. Four repetitions per arm reduce that estimate to 112 targets per
curve but raise the total to 1,792 solves. These are pilot-based planning
numbers, conditional on the observed effect persisting, independent targets,
unchanged input laws, and a stable host. The local A/A spreads fail the 10%
noise gate, and neither pilot has a host-wide exclusive CPU receipt; a Modal
guest-fidelity pilot must first recompute the variance. Additional local
repetitions on the same contended host would not satisfy the isolation gate.

The separately frozen [matched rank-243 study](../README.md) asks a different
cross-curve question and retains 200 target indices per curve. Its candidate
and baseline IDs must not be substituted for the full-rank pairs analyzed
here.

Reproduce this analysis from the repository root:

```sh
python3 experiments/koblitz-s3-crosscurve-20261009/fullrank-contrast/analyze.py \
  --output /private/tmp/s3-fullrank-crosscurve-analysis.json
```
