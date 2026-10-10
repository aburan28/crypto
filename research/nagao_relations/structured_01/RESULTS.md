# Structured-base scaling: the direct S3 table changes the comparison

Executed 240 matched cold trials, including all eight frozen d8 targets and sixteen fresh holdouts, plus ten charged eight-target enumeration batches. All saved answers and complete sets pass independent replay. The pair-invariant method is a direct Semaev solver and is therefore a stronger comparator to the hybrid; a win for it cannot be claimed as a function-first advantage.

[Proof and contract](README.md) · [raw evidence](raw.jsonl) · [summary](summary.json) · [paired counter comparison](comparison.json) · [audit](audit.json)

The matched three-second cold results are:

| Variant | First-mode resolved / 24 | Exact enumeration / 24 | Common-cost speedup | Correctness failures |
|---|---|---|---|---|
| hybrid-reference | 13/24 | 0/24 | unmeasured | 0 |
| hybrid-filtered | 13/24 | 0/24 | unmeasured | 0 |
| pair-invariants-s3 | 14/24 | 14/24 | unmeasured | 0 |
| chained-s3 | 0/24 | 0/24 | unmeasured | 0 |
| s4-symmetric | 0/24 | 0/24 | unmeasured | 0 |

The table is an engineering diagnostic. First-mode resolutions include proofs that a target has no admissible relation. At 18 bits the hybrid resolves more first-mode trials (11/12) than the table (8/12); the table advantage in complete enumeration does not establish universal first-hit superiority. All d10 cold table trials time out while building the table.

## Matched cold trials

First-relation resolutions include proved empty targets; enumeration resolves only after the exact whole set is returned. Timeout is unknown and partial enumeration remains incomplete. The three-second budget includes supplied-basis expansion, curve and solver setup, search, extraction and verification. SAT soft overruns are retained.

18-bit field, first; resolved within three seconds / attempts:

| Variant | prefix d8 | prefix d9 | prefix d10 | f4-stable d8 | f4-stable d10 |
|---|---|---|---|---|---|
| hybrid-reference | 4/4 | 2/2 | 2/2 | 1/2 | 2/2 |
| hybrid-filtered | 4/4 | 2/2 | 2/2 | 1/2 | 2/2 |
| pair-invariants-s3 | 4/4 | 2/2 | 0/2 | 2/2 | 0/2 |
| chained-s3 | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| s4-symmetric | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |

18-bit field, enumerate; resolved within three seconds / attempts:

| Variant | prefix d8 | prefix d9 | prefix d10 | f4-stable d8 | f4-stable d10 |
|---|---|---|---|---|---|
| hybrid-reference | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| hybrid-filtered | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| pair-invariants-s3 | 4/4 | 2/2 | 0/2 | 2/2 | 0/2 |
| chained-s3 | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| s4-symmetric | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |

30-bit field, first; resolved within three seconds / attempts:

| Variant | prefix d8 | prefix d9 | prefix d10 | f4-stable d8 | f4-stable d10 |
|---|---|---|---|---|---|
| hybrid-reference | 2/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| hybrid-filtered | 2/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| pair-invariants-s3 | 4/4 | 0/2 | 0/2 | 2/2 | 0/2 |
| chained-s3 | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| s4-symmetric | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |

30-bit field, enumerate; resolved within three seconds / attempts:

| Variant | prefix d8 | prefix d9 | prefix d10 | f4-stable d8 | f4-stable d10 |
|---|---|---|---|---|---|
| hybrid-reference | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| hybrid-filtered | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| pair-invariants-s3 | 4/4 | 0/2 | 0/2 | 2/2 | 0/2 |
| chained-s3 | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |
| s4-symmetric | 0/4 | 0/2 | 0/2 | 0/2 | 0/2 |

## Charged pair-table batches

Each row builds a new S3 table and enumerates eight new targets: four uniform and four sampled from signed factor triples. No hybrid or SAT batch comparison was run, so these are capacity/amortization diagnostics. Setup remains quadratic. The batch wall clock also includes the independent harness oracle checks between queries; their time plus orchestration is exposed separately in summary.json. The table below reports the measured setup and sum of timed queries, which both include solver-side point verification. Counts are projected relations per target, not matrix-independent rows.

| Bits | Base | d | Signed points | S3 table entries | Setup s | Eight queries s | Batch wall s | Uniform relations | Supported relations |
|---|---|---|---|---|---|---|---|---|---|
| 18 | prefix | 8 | 254 | 16002 | 0.4646 | 0.4815 | 1.3698 | 34 | 42 |
| 18 | prefix | 9 | 512 | 65280 | 2.0733 | 1.4365 | 4.6003 | 320 | 364 |
| 18 | prefix | 10 | 1068 | 284622 | 9.0493 | 6.9073 | 20.0828 | 2943 | 3021 |
| 18 | f4-stable | 8 | 148 | 5402 | 0.1729 | 0.2755 | 0.6898 | 9 | 17 |
| 18 | f4-stable | 10 | 976 | 237656 | 7.6875 | 6.6089 | 18.0086 | 2338 | 2366 |
| 30 | prefix | 8 | 280 | 19460 | 1.0241 | 1.3052 | 3.5868 | 0 | 4 |
| 30 | prefix | 9 | 526 | 68906 | 3.5397 | 2.4337 | 8.3202 | 0 | 4 |
| 30 | prefix | 10 | 1010 | 254520 | 12.9623 | 4.6446 | 22.0892 | 0 | 5 |
| 30 | f4-stable | 8 | 304 | 22952 | 1.2467 | 1.4467 | 4.0600 | 0 | 4 |
| 30 | f4-stable | 10 | 1012 | 255530 | 13.7105 | 4.5404 | 22.5954 | 1 | 4 |

## What the experiments establish

- Both structured dimensions are real F4-linear, Frobenius4-stable subspaces; they are not contained in a proper subfield. Structured d9 is impossible because F4-linearity forces even binary dimension. Measured signed-point support varies substantially: at 18 bits and d8 the structured base has 148 points versus 254 for the prefix base. Their timing difference cannot be interpreted as a fixed-support speedup.
- The early filter preserves exact function acceptance on 10,638 candidate checks. All 1,920 parity-membership checks agree with Gaussian image recovery. All five solvers agree with exhaustive signed triples on 106 tiny target/space pairs.
- Independent replay verifies 11985 signed-point certificates and every larger complete oracle set. Frozen targets and source hashes match. No correctness failures were found.
- The saved Frobenius counterexample shows that conjugated factors remain in the structured base but sum to R^4, not the original R. Stable support does not justify fixed-target orbit collapse.
- Ordinary setup caching was overestimated in the proposed plan: previous d8 setup was only 0.24455% and 0.20975% of wall time. The new pair table is different: it deliberately pays quadratic precomputation once per named batch.

The counting ceiling min(1,8*C(M/2,3)/(#E-1)) changes when the base changes. No subquadratic cold-search theorem or exponent fit follows from these cells. The direct S3 table uses quadratic setup; the hybrid still enumerates O(2^(2d)) branches. Any lower field-API sum for the parity filter must also account for its new binary word operations. comparison.json retains both vectors and makes incomplete-pair ratios null.

## Hypothesis decisions

| Proposal | Observation or proof | Decision |
|---|---|---|
| Cache ordinary image setup for a 20% gain | Even removing all prior setup saves under 0.25% | Disproved for this measured workload |
| Early parity rejection improves completion | Same 13/24 first resolutions and 0/24 complete enumerations as the reference | Correct filter; completion hypothesis not supported |
| Stable support permits fixed-target Frobenius pruning | Saved conjugated relation sums to a different target | Disproved without a target stabilizer |
| Pair-invariant S3 strengthens enumeration | 14/24 cold completions versus 0/24 for both hybrids | Retain as a stronger Semaev baseline |
| Charged reuse reaches larger bases | All ten eight-target batches complete, including d10 at 30 bits | Capacity established; setup remains quadratic |
| Larger structured bases guarantee useful random-target yield | Only one of twenty uniform 30-bit batch targets has a relation | Yield remains an obstacle; supported targets do not estimate it |

The filter uses fewer field-API calls on ten of the thirteen matched completed first-relation searches, and more on three. Its candidate/reference API-sum ratios range from 0.5037 to 1.1799; these exclude its separately recorded binary word work and are not common-operation speedups. Regressions remain in comparison.json.

The next function-first experiment must reduce coefficient branches and compare against the direct S3 table on identical cold and named-batch workloads. Improving the early filter alone leaves O(2^(2d)) coefficient branching intact. Random-target yield and first-relation latency must remain separate from complete enumeration.

Classification: engineering diagnostic. Calibrated common-operation speedup, full-DLP S, rho/floor ratios and the original 20% goal remain unproved. These small family-specific trials include frozen replays and fresh holdouts, but do not pass the broader 60-input/repetition and full-pipeline promotion gates. No full ECDLP pipeline or scalar recovery was run in this campaign.

## Replay

Run `python research/nagao_relations/structured_01/analyze.py` to rebuild the independent group-law oracles, check every saved result and certificate, and regenerate the report. An identical existing scoreboard panel is accepted; a different panel is rejected. This analysis is outside the measured budget.

For a timing rerun, use a separate checkout of the published frozen source commit in publication.json, with Python 3.12 and pycryptosat 5.14.7. The measured runner refuses to overwrite raw.jsonl. The analysis script and reports were added after the campaign and are not part of its frozen measured source.
