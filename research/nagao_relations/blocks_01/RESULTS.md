# Block rejection: measured solver-stage outcome

[Derivation and rejection proof](README.md) · [contract](contract.json) · [raw](raw.jsonl.gz) · [summary](summary.json) · [paired comparison](comparison.json) · [audit](audit.json)

The 192 three-second cold trials cover 16 identical inputs per variant: eight frozen and eight fresh, at 18 and 30 bits with dimensions 8 and 9. All six variants are rerun. First-mode resolutions include proved empty targets. Timeouts remain unknown and over-budget completions do not count as resolved.

| Variant | First resolved | Enumeration complete | Correctness failures | Common-cost speedup |
|---|---|---|---|---|
| bilinear-reference | 15/16 | 12/16 | 0 | unmeasured |
| bilinear-blocks | 14/16 | 6/16 | 0 | unmeasured |
| hybrid-filtered | 11/16 | 0/16 | 0 | unmeasured |
| pair-invariants-s3 | 14/16 | 14/16 | 0 | unmeasured |
| chained-s3 | 0/16 | 0/16 | 0 | unmeasured |
| s4-symmetric | 0/16 | 0/16 | 0 | unmeasured |

Classification: engineering diagnostic; no advance or calibrated total-cost speedup established. One timing repetition cannot establish a runtime improvement. The separate field and binary counters do not include a measured conversion for SAT, allocation or control. The broad regression, full-DLP and three-size goal gates remain open.

## Matched cost diagnostics

Every numeric cost cell below is an equal-weight sum of field additions, multiplications and squarings, with expanded inversion internals. Ratios use the direct S3 costs on the same jointly resolved tasks. Each row can have a different subset, so compare its ratio within the row, not totals across rows. API ratios are component diagnostics. SAT has no resolved pair here unless shown. Complete primitive vectors, binary work and every regression are retained in comparison.json.

| Variant | Mode | Joint resolved tasks | Field API sum | API ratio to matched S3 | Total-cost ratio | Ratio to floor |
|---|---|---|---|---|---|---|
| bilinear-reference | first | 14 | 6,520,230 | 1.0755 | unmeasured | unmeasured |
| bilinear-reference | enumerate | 12 | 21,972,578 | 5.2209 | unmeasured | unmeasured |
| bilinear-blocks | first | 14 | 7,296,527 | 1.2036 | unmeasured | unmeasured |
| bilinear-blocks | enumerate | 6 | 8,494,561 | 3.2210 | unmeasured | unmeasured |
| hybrid-filtered | first | 11 | 1,646,516 | 0.3471 | unmeasured | unmeasured |
| hybrid-filtered | enumerate | 0 | — | — | unmeasured | unmeasured |
| pair-invariants-s3 | first | 14 | 6,062,236 | 1.0000 | unmeasured | unmeasured |
| pair-invariants-s3 | enumerate | 14 | 6,280,210 | 1.0000 | unmeasured | unmeasured |
| chained-s3 | first | 0 | — | — | unmeasured | unmeasured |
| chained-s3 | enumerate | 0 | — | — | unmeasured | unmeasured |
| s4-symmetric | first | 0 | — | — | unmeasured | unmeasured |
| s4-symmetric | enumerate | 0 | — | — | unmeasured | unmeasured |

## Predeclared pruning screen

The batch gate **passes**. There are 6 complete equal-output nonzero-target enumeration pairs between the pruned and unpruned circuits. The gate requires at least four, 20% fewer field API calls and no increase in separately counted binary word work.

```json
{
  "complete_equal_enumeration_pairs": 6,
  "reference_field_api_sum": 10981053,
  "candidate_field_api_sum": 8494561,
  "candidate_over_reference_field_api": 0.773565249161442,
  "reference_binary_words": 6637332,
  "candidate_binary_words": 378808,
  "passed": true,
  "interpretation": "Exploratory screening only; no calibrated total-cost claim."
}
```

A failed screen ends this candidate before optional matched batches. Pruned-branch counts alone cannot reverse the verdict. Profiles include work actually performed on timed-out searches and must not be read as complete-work speedups.

## Matched coverage by size

| Bits | d | Variant | First resolved | Enumeration complete |
|---|---|---|---|---|
| 18 | 8 | bilinear-reference | 6/6 | 6/6 |
| 18 | 8 | bilinear-blocks | 6/6 | 0/6 |
| 18 | 8 | hybrid-filtered | 6/6 | 0/6 |
| 18 | 8 | pair-invariants-s3 | 6/6 | 6/6 |
| 18 | 8 | chained-s3 | 0/6 | 0/6 |
| 18 | 8 | s4-symmetric | 0/6 | 0/6 |
| 18 | 9 | bilinear-reference | 2/2 | 0/2 |
| 18 | 9 | bilinear-blocks | 2/2 | 0/2 |
| 18 | 9 | hybrid-filtered | 2/2 | 0/2 |
| 18 | 9 | pair-invariants-s3 | 2/2 | 2/2 |
| 18 | 9 | chained-s3 | 0/2 | 0/2 |
| 18 | 9 | s4-symmetric | 0/2 | 0/2 |
| 30 | 8 | bilinear-reference | 6/6 | 6/6 |
| 30 | 8 | bilinear-blocks | 6/6 | 6/6 |
| 30 | 8 | hybrid-filtered | 3/6 | 0/6 |
| 30 | 8 | pair-invariants-s3 | 6/6 | 6/6 |
| 30 | 8 | chained-s3 | 0/6 | 0/6 |
| 30 | 8 | s4-symmetric | 0/6 | 0/6 |
| 30 | 9 | bilinear-reference | 1/2 | 0/2 |
| 30 | 9 | bilinear-blocks | 0/2 | 0/2 |
| 30 | 9 | hybrid-filtered | 0/2 | 0/2 |
| 30 | 9 | pair-invariants-s3 | 0/2 | 0/2 |
| 30 | 9 | chained-s3 | 0/2 | 0/2 |
| 30 | 9 | s4-symmetric | 0/2 | 0/2 |

Cohort and uniform/supported strata remain separate in summary.json. No random-target yield estimate is inferred from supported targets.

## Correctness and boundaries

The independent replay checked 512 relation certificates and 193,763 rejection certificates in the larger trials, with zero failures. Tiny validation covers all 53 finite GF64 targets on three d4 spaces: 159 target/space pairs, 318 new-solver runs, 7,018 rejection certificates and 7,168 exact circuit identity evaluations including fresh 18/30-bit samples. Both new solvers agree with exhaustive signed triples and an independent pair oracle.

The signed-triple success ceiling from bounds_01 is unchanged. On this prefix support it is 100% at 18 bits (a vacuous ceiling), and 0.333461% / 2.233163% at 30 bits for d8 / d9. These are counting ceilings over uniform affine targets in the full curve group, not measured success rates or prime-subgroup claims. The old seven-multiplication branch floor does not apply after this circuit change. No ratio to an invented new full-DLP floor is reported.

Raw SHA256: `c0d1a3cc79a0d12fd569248dcd1a99ca256ffcd4d243770a68268969cb285a5f`. All recorded source and input hashes match. Run analyze.py to reproduce the audit and these literal tables. Frozen solver source is mapped to the published Git tree in publication.json.

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
