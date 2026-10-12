# Cold single-target continuation

The initial tournament operation was published in [PR #372](https://github.com/aburan28/crypto/pull/372), which is merged. These rounds continue the single-target request.

## Fixed-support round 0006b

All 1,584 profiled trials and paired native runs verified. The best challenger, `combined_batch1`, reduced complete cold instructions by 21.46% and native process time by 11.54% in confirmation. Replay agreed. It did not meet the predeclared >=20% reduction in both metrics, so the incumbent was retained. Single-target rho parity remains open.

- [Plan](PLAN.md), [full report](../runs/round-0006b-single/REPORT.md), [decision](../runs/round-0006b-single/decision.json).
- [Candidate scope check](candidate-scope-check.json) and [orbit-table equivalence test](preflight.log).
- [Retained preparation failure](../runs/round-0006-single/prepare_failure.json): CPU affinity was rejected before any fixtures or trials. The successor rehashed all source and reused exact completed worker builds.

## Separate policy round

The [predeclared policy panel](POLICY_PLAN.md) tests explicit smaller factor bases on the same cold single-target ECDLP workloads. Base support and rank floors can differ across arms; their actual values must be reported. Round 0007 promoted `orbits2`: confirmation cost was 0.3753 of the incumbent in instructions and 0.7091 in native time; replay passed again. Against matched rho it still cost 1.1520 instructions and 1.1008 time, so neither parity nor strict beating was established.

The audit checked 1,560 receipts: 1,542 verified solves and 18 retained smoke rejections. `orbits1` (12/12 rejected) and `cube_root` (6/12 rejected) encountered the frozen checker's minimum of two logarithm columns. They were excluded before development. This is a protocol limitation, not evidence that their final scalars were wrong.

[Full report](../runs/round-0007-single-policy/REPORT.md), [all-candidate admission](../runs/round-0007-single-policy/admission.json), [strict rho verdict](../runs/round-0007-single-policy/strict-rho.json).

## Fixed-support implementation continuation

[Round 0008 plan](IMPLEMENTATION_PLAN.md) tests orbit-wise cofactor projection and serial collection on the promoted two-orbit support. Both preflight equivalence tests passed. The round finished and audited: **retained — incumbent**. See [combined results](RESULTS.md) and [strict rho verdict](../runs/round-0008-single-implementation/strict-rho.json).
