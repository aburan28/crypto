<!-- block-batch-results -->

## Fresh eight-target batch with shared setup

[Batch contract](batch_contract.json) · [targets](batch_targets.json) · [compressed raw](batch_raw.jsonl.gz) · [summary](batch_summary.json)

The cold screen selected n30, prefix d8 for this follow-up. Every variant gets the same four fresh uniform and four fresh supported targets, with 60 seconds for the entire eight-target enumeration batch. Field/curve/support setup is charged once per batch; S3 pair tables and hybrid support images are also charged once. Target-dependent coefficient circuits are rebuilt and charged per target.

| Variant | Add calls | Mul calls | Square calls | Field API sum | API ratio to S3 | Complete / 8 | Verified relations |
|---|---|---|---|---|---|---|---|
| bilinear-reference | 13,974,348 | 655,536 | 8,798 | 14,638,682 | 25.2091 | 8 | 4 |
| bilinear-blocks | 10,315,694 | 655,536 | 8,798 | 10,980,028 | 18.9086 | 8 | 4 |
| hybrid-filtered | 19,837,230 | 5,515,666 | 585,578 | 25,938,474 | 44.6684 | 8 | 4 |
| pair-invariants-s3 | 341,254 | 207,401 | 32,035 | 580,690 | 1.0000 | 8 | 4 |

This table is in field API calls, not a calibrated machine-operation unit. Complete workloads include empty targets and all verification. The following separate counter tracks only the explicitly instrumented binary operations; zero does not mean zero machine bit work. These components cannot be added without a measured conversion.

| Variant | Selected binary word operations | Calibrated total-cost ratio | Ratio to full-DLP floor |
|---|---|---|---|
| bilinear-reference | 8,851,713 | unmeasured | unmeasured |
| bilinear-blocks | 486,006 | unmeasured | unmeasured |
| hybrid-filtered | 1,153,470 | unmeasured | unmeasured |
| pair-invariants-s3 | 0 | unmeasured | unmeasured |

Wall time is a single observation, in seconds; no confidence interval or unqualified runtime improvement is claimed:

| Variant | Setup seconds | Queries seconds | All phases seconds |
|---|---|---|---|
| bilinear-reference | 0.0006 | 17.5258 | 17.5265 |
| bilinear-blocks | 0.0004 | 11.5901 | 11.5905 |
| hybrid-filtered | 0.0529 | 54.2583 | 54.3113 |
| pair-invariants-s3 | 0.9640 | 1.4257 | 2.3898 |

The replay checks 16 relation certificates and 12,902 block rejection certificates, with zero failures. Uniform and supported outcomes: {"known_decomposable": {"targets": 4, "targets_with_relations": 4}, "uniform": {"targets": 4, "targets_with_relations": 0}}.

The original three-size calibrated cost goal remains unmet. Both the cold regressions and this selected follow-up must accompany any claim about pruning.

The frozen constructor alone requires 530,432 multiplications, versus 207,401 measured total S3 multiplications on this batch. Even free perfect pruning cannot remove that cost. See [the new bound and next hypothesis](BOUNDS_AND_NEXT.md).
<!-- block-batch-results-end -->
