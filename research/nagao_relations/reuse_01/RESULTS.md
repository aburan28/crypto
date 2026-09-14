# Coefficient reuse: complete audit

Engineering diagnostic. The original calibrated three-size goal remains unmet. All variants use identical F4-defined curves, prefix supports, targets and arity three.

[Contract](contract.json) · [raw](raw.jsonl.gz) · [continuation](continuation_raw.jsonl.gz) · [batch raw](batch_raw.jsonl.gz) · [audit](audit.json) · [paired comparisons](comparison.json)

Cold 3-second screen, 16 frozen and eight fresh targets at n18/n30, d8/d9. The unresolved original S3 runner incident is excluded from completion counts and cost pairs, even when replay succeeds. These counts are conservative and depend on the budget.

| Variant | First / 24 | Enumeration / 24 |
| --- | --- | --- |
| bilinear-blocks | 21 | 8 |
| reuse-circuit | 21 | 8 |
| reuse-batch-inverse | 21 | 8 |
| pair-invariants-s3 | 20 | 19 |
| chained-s3 | 0 | 0 |
| s4-symmetric | 0 | 0 |

Complete enumeration of the same eight n30 d8 targets as the predecessor batch: four uniform and four supported. Setup, failed targets, extraction and exact verification are charged.

| Variant | Add | Mul | Square | Field API sum | API ratio to S3 | Complete / 8 | Relations |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bilinear-blocks | 10,315,694 | 655,536 | 8,798 | 10,980,028 | 18.9086 | 8 | 4 |
| reuse-batch-inverse | 10,405,574 | 87,552 | 2,086 | 10,495,212 | 18.0737 | 8 | 4 |
| reuse-circuit | 10,405,574 | 201,768 | 2,086 | 10,609,428 | 18.2704 | 8 | 4 |
| pair-invariants-s3 | 341,254 | 207,401 | 32,035 | 580,690 | 1.0000 | 8 | 4 |

Field API calls are an uncalibrated equal-weight component sum. The next table keeps selected binary work separate; zero is not zero machine work. Cached slots are logical storage, not RSS. Timing is one observation, with no confidence interval or runtime claim.

| Variant | Selected binary word operations | Cached field slots | Observed seconds |
| --- | --- | --- | --- |
| bilinear-blocks | 486,006 | 0 | 10.7529 |
| reuse-batch-inverse | 486,006 | 20736 | 7.0503 |
| reuse-circuit | 486,006 | 20736 | 8.2667 |
| pair-invariants-s3 | 0 | 0 | 2.2577 |

The audit replays 826 relation certificates and 929,262 rejection certificates with zero failures. Complete reuse batches preserve the predecessor rejection certificates and rejection-phase counter vector exactly.

Changing coefficient preparation alone leaves a floor of 9,794,908 rejection field API calls, 16.8677 times the measured total S3 API count. This is a bound on this traversal and this unit, not a universal or calibrated lower bound. Cached fixed column spans are the next falsifiable change. See [execution notes](EXECUTION_NOTES.md).
