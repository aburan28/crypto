# Degree-83 a=1/K600 evaluation against strong rho

**Status: all five matched pairs completed; all 22 fresh answers independently verified.**
This directory records a new native evaluation of the retained one-point degree-83 IC observation. No performance row is promoted by this diagnostic study. The result applies to five repeats of one fixed public target on one virtualized host.

The frozen paired geometric mean gives **IC/strong online cost 1.382330**, with 95% paired-bootstrap interval **[1.035774, 1.882300]**. The corresponding fresh-process wall-cost ratio is **218.034569**, interval **[172.358384, 285.573727]**. These are reciprocals of the preregistered strong/IC estimator, not ratios of arithmetic means. The strong/strong A/A online ratio was **0.997036**, interval **[0.967928, 1.023496]**. There were no excluded pairs, failed arms or resource caps in the evaluation panel.

| Secondary descriptive mean across five formal runs | Strong rho | IC K600 |
| --- | ---: | ---: |
| Complete prepared-online interval, ms | 4837.039151 | 6329.179072 |
| Reusable preparation, ms | 191.515873 | 1047924.982921 |
| Fresh-process wall, s | 5.067556 | 1054.552457 |
| Child user + system CPU, s | 5.029242 | 1054.438292 |
| Peak RSS, bytes | 116588544 | 5811525222.4 |

The arithmetic means above are explicitly secondary descriptions computed by the separate native exporter. All five raw pairs, including the fourth pair where IC was faster online, are retained in [records/panel/statistics.json](records/panel/statistics.json). The native exporter preserves failures, missing values and all charged costs without zero imputation: [report_data_secondary.json](records/panel/report_data_secondary.json).

The independent native checker replayed 16 strong answers and six IC answers/witnesses, with zero invalid rows. Native reservation monitoring recorded zero contended samples and exit status 0. Original factor-base selection/construction and OS-cold cache conditions are outside this fresh-process boundary. The 100 ms outer exit poll also affects short process-wall observations; CPU and producer timings corroborate them. See [record provenance](records/README.md) and the source-linked PDF under `reports/run-20261012T0226Z/`.

The frozen design is [PROTOCOL.md](PROTOCOL.md), with executable runs and estimator in [runner/protocol.json](runner/protocol.json). Source baseline: a8194d616ed8c04489ec46177d369f1b1670360a. Current-source measurements are separate evidence from the archived October 7 binaries.

## Exact fixture and identity correction

The task concerns the legacy **a=1, approximately 52.928-bit prime subgroup**, not the a=0 degree-83/81-bit primary sweep. Both measured methods receive the identical public point without its validation-only scalar.

| Item | Exact value |
| --- | --- |
| Curve | y² + xy = x³ + x² + 1 |
| Field | GF(2)[z]/(z^83 + z^7 + z^4 + z² + 1) |
| Subgroup order r | 8569786107849059 |
| Cofactor h | 1128547018 |
| Group order | 9671406556910881928556062 |
| Generator | (6614244651398184860292015,8868260410029249391539529) |
| Supplied public Q | (2655898345537166180532960,3881526307991673265292615) |
| Frobenius eigenvalue | 2998004823950457, exact order 83 |
| Signed-Frobenius action | 166 elements |
| Retained base | 600 signed-orbit columns, 99,600 points |
| Declared point-set digest | d4a4da26d6e7cb09367de922d9806631fe49989f70d378cae76f347f8d76a9f7 |

Inputs:
[base_n83_K600.jsonl](../../experiments/koblitz-single-target-n83-20261006/base_n83_K600.jsonl)
and
[target_points.jsonl](../../experiments/koblitz-single-target-n83-20261006/frozen/target_points.jsonl).
The coordinate words belong to the exact polynomial representation above.

**Identity erratum, without changing frozen inputs or protocol:** the strong worker computes the legacy representation's slug as icv1-f2m83-t6151469093347-9eba6cd5, with full identity:

ICV1:f2m-83-a5fd0d39:6151469093347:9671406556910881928556062:0x1:-7:unk:r:9eba6cd5d359

The frozen protocol prose and archive checker's label used historical icv1-f2m83-t6151469093347-cdcc5432. That label identifies the different z^83+z^45+z²+z+1 diagnostic representation. The actual modulus, generator, Q, subgroup, input hashes and worker assertions bind this evaluation correctly. Preserve frozen files and report the nickname correction; do not merge the models by degree or subgroup order. In particular, known_n83_k1() is not the constructor for this retained coordinate fixture.

## Historical reference correction and new adapter

The [archived runner](../../experiments/koblitz-single-target-n83-20261006/run_n83_pairs.py) used the packed backend of koblitz_rho_fixture at n=83/a=1/signed_frobenius. Its three recorded rho/IC ratios were approximately 8.181, 5.550 and 12.671. Both answers independently replayed, but this backend is now excluded by the [reference contract](../notes/index-calculus/RESEARCH_FIXED_REFERENCE_RHO_20261001.md).

Packed rho stores every step, scans the signed-Frobenius orbit in polynomial coordinates and inverts per step. Strong rho uses distinguished points, normal-basis canonicalization and batched inversion across lockstep lanes. The old ratio is therefore not evidence against the required strong control. Host contention, unmatched resource information, incomplete manifests and phase accounting also prevented controlled admission. Old measurements and correctness receipts remain unchanged.

The existing fixture CLI accepts only packed mode above degree 63, although the generic strong library supports wide fields. The new [n83_strong_point_eval.rs](../../examples/n83_strong_point_eval.rs) calls StrongRhoG<WideGf2,u64>::from_parts on KoblitzCurve::new(1,83), asserting the exact legacy modulus, r, h and generator.

The adapter accepts one JSON coordinate pair and a walk seed. It refuses scalar sidecars, multiple targets and noncanonical field words. Strong defaults remain 32 lockstep lanes within one operating-system worker, DP bits 4, step-cap factor 2000 and fruitless escape disabled. Invalid, capped and unverified outcomes exit nonzero. Small-field tests compare narrow/wide strong trajectories and charges; their separate tiny-quotient DP policy does not change production n=83 parameters.

## Frozen execution and accounting

The [native runner](runner/README.md) freezes 22 child processes: one excluded warmup per method, five same-seed strong/strong A/A pairs, then five interleaved strong/IC pairs. Strong seeds are 202610061–202610065 and IC rank seeds 21–25. Matched order alternates by round.

The common enforced envelope is one pinned CPU, CPU 3; one IC rank worker and one target worker; 12 GiB virtual-address limit through RLIMIT_AS; and 1800-second process-group timeout per arm. Children start with empty DP/index/log state and a clean whitelisted environment. Create-new raw stdout/stderr, command/environment/input/binary hashes, monotonic process wall, user/system CPU, Linux wait4 peak RSS scaled from KiB to bytes, and semantic exit/signal/cap disposition are retained. No builds, tests or independent heavy jobs run alongside a timed child.

IC reusable curve/base/normal-basis/root-index/rank/log preparation is excluded from online and separately recorded. Strong curve/normal-basis/G-only jump preparation is likewise reusable. Strong online includes Q membership/subgroup validation, Q-dependent starts, walk/collision work and a separate generic scalar replay.

The current packed jump setup constructs aG+bQ, so it is target dependent. Strong jumps are G-only. State this difference instead of assigning both algorithms the same setup boundary.

IC preserves query, PDP, relation check, descent and recovery-check costs. Its relation-check timer covers the final lift addition/equality; earlier lift work remains in PDP. Preserve the residual between their sum and target_ms. Parallel worker check-duration sums are not exclusive wall time, so target threads stay one. Exit zero alone is insufficient: the runner requires one verified scalar row on Q, full rank 600, expected base digest, one solved target and zero failed targets.

The primary estimator is the geometric mean of within-pair strong_online_ms / IC_online_ms. The prospectively fixed 95% paired percentile bootstrap uses 20,000 draws over mean log ratio, seed 2026101283. Whole-process ratios remain separate. Failed, mismatched, capped or unverified pairs remain explicit exclusions and never wins. One pair receives no interval. Repeats concern one point and one host, not a fresh-target population.

Native isolated_bench supplies reservation/condition monitoring. It does not by itself establish sealed ecbench L2, physical-host exclusivity or independent-host admission. The sealed IC route does not yet admit this exact legacy wide adapter. Timings remain diagnostic without those remaining gates. Common operation-equivalent S is unknown: IC root-index probes and rho walk steps are different units.

## Frozen hashes

The mathematical set digest differs from the SHA-256 of serialized input bytes. Per-child start receipts bind actual inputs and executable.

| Artifact | SHA-256 |
| --- | --- |
| Retained base bytes | 0d290b14a350db6271c9a49f0cd6ba246edd8c9bc6754387ad1ceebf08ec0c09 |
| Public target file | d4f9872d9d4e22d26770ecb749b29f1da944d701fb85c8bbc8f20489cd0345c5 |
| Archived validation-only fixture | f3605f0d1739e4542ced9b5c12200c5582f60ae3a4b94bc48189bf525873aea5 |
| Frozen baseline IC executable | 983f28400b7852fc6628ba97d3ea83b2ab9e32d9944658c7eea831ee89a4abf2 |
| Frozen strong executable | 98949b7127baeb1fa55420d23dc13604069fcc8fde8d92b8f2e6bebd09718542 |
| Strong adapter source | 82bc92d24e096109d4c35288d034c02f63c08241ed87dc83c93655ed668c4cbf |
| Original compiled native runner source | d6d39bcd6e9df4b477aeb4606d1073e7f0bd42a3957737e61e19ef953ed4d557 |
| Formatted delivery runner source, portable test path | 063f6424da8c8c8816742ab337099fb1791a4204d07185155bcea61c87fdf681 |
| Runner protocol JSON | ac7266f6caad2a4142ae06884ee748f980f9ee200c4e2649a0149f59e91adf29 |
| Frozen prose protocol | b98c9443fc6d3f748584863f5e971780bb32674ae84df73c32c113826de52d4f |

Executed identities are in [records/frozen_panel_inputs.sha256](records/frozen_panel_inputs.sha256) and every raw child receipt; preparation snapshots have separate provenance. A filename alone is not a source/executable binding.

## Independent native checks and validation

The standalone C++ archive checker uses separate GF(2^83) and affine-curve arithmetic. Its completed archive receipt checks irreducibility, deterministic 64-bit primality, Frobenius-recurrence group order, subgroup/cofactor product, G/Q membership and annihilation, eigenvalue order 83, signed action 166, all 600 representatives and 99,600 point/label records, and the actual four-point witnesses plus both recovered scalars for all three archived runs.

It does not reconstruct the absent 600 rank equations/factor-log database, retime old runs or qualify their hardware. Its historical hard-coded curve label has the erratum above; arithmetic uses the correct legacy polynomial and coordinates.

| Independent artifact | SHA-256 |
| --- | --- |
| Archive checker source | b3ab0742f072ae19595707afb906b89e143bc91a43dd903655cdeeb8d2bf090c |
| Archive checker executable | 4717970bb693912861ec6cfd27cc95c94440872977b2295c674c8e2334bfa3de |
| Archive checker receipt | 7d4e290524ea9b83ff75d48d20b17449fb60a6991bae787fdfee53ef28d461d5 |
| Fresh checker source prepared before final replay | 09e8162dc9157daf46777789fd09a24fca52cdd128d566fae1f0c892a7cd7115 |
| Fresh checker executable | aac8afbe40366845e069318dcda8fd0fa3a1889e2c89d627055be52bbe19c559 |

Fresh panel replay passed all 22 answers and all six IC witnesses. Receipt SHA-256: `980a3d7d52db96a559dd36ad12e317c9b12be9f1c161a2d0d86803aafcc72441`. The unexported rank equations and factor-log database remain unverified independently.

Before timing, the existing IC example's 10 tests and the new strong adapter's two tests passed. The initial tiny-curve failure is retained: default DP bits could select no distinguished class in the n=7 quotient. Only test parameters changed to DP bits zero/lanes four; n=83 production defaults did not change.

The initial full library run had 3937 passes, 112 ignored and 10 failures: nine sandbox-denied socket tests and one baseline generic-collector defect. The permitted pollard_collab rerun passed 37 tests. The unchanged failed sampler test passed after the repair below. These are scoped reruns, not a second complete-suite claim.

The final complete default-feature `cargo test --release --locked --lib` run passed **3947 tests, with zero failures and 112 ignored**, in 260.41 seconds. Network permission allowed the tests' local sockets. The delivery runner passed seven tests; the secondary exporter passed four checks, including null/censored charges and reciprocal interval handling. Initial failures and corrected attempts remain in `validation_logs/` and `evidence/runner_checks/`.

The final default-feature touched-example checks also passed: **10 compact-orbit IC tests and two exact-point strong-rho tests**, with zero failures. See `validation_logs/examples_tests_final.log`.

## Baseline discarded-draw defect

The baseline generic collector in [koblitz_index_calculus.rs](../../src/cryptanalysis/koblitz_index_calculus.rs) drew two nonzero coefficients with sample_relation_scalar, discarded them by shadowing, and then drew the intended zero-inclusive additive a and nonzero target b. The discarded calls still advanced the RNG.

The unchanged test cryptanalysis::koblitz_index_calculus::tests::n83_relation_queries_use_full_width_zero_inclusive_additive_draw expects gen_biguint_below(r), then gen_biguint_below(r−1)+1 per query. Its baseline mismatch was 899991076855513638720967 versus 1842027667052925756179639. The repair removes **only two discarded draws**, preserving the intended coefficients, seed and test assertions.

The standalone compact-orbit example uses its own guided-rank/scan code; it never calls this collector. The new adapter uses strong rho, also outside that path. The measured IC executable was frozen before repair and remains unchanged. This is a deterministic reproducibility fix, not a measured IC optimization or a reason to rewrite old receipts.

## Delivery and remaining qualification

| Deliverable | Current status |
| --- | --- |
| Protocol, native runner, strong adapter | Prepared; tests recorded |
| Native archived fixture/base/witness replay | Complete with stated limits |
| Warmup dispositions | Complete, excluded from formal estimates |
| Five strong A/A and five matched pairs | Complete; all children verified |
| Online/fresh-process statistics and fresh replay | Complete; native evidence with stated scope |
| Raw start/heartbeat/final receipts, stdout/stderr | Delivered in records/panel/ |
| Source-linked report, quantitative graph and inspected PDF | Complete; 14 pages, all pages inspected |
| Controlled vs_rho admission or ledger promotion | Not performed |

The timed panel has ended. The remaining qualification gaps are the full rank/log transcript, calibrated common operation S, canonical ecbench adapter/manifests for this exact legacy wide curve, physical-host/L2 requirements and independent-host replication. Five repeats do not establish a distribution over new targets or another curve.

## Canonical graph review

A narrow search used the retained a=1/K600 run prefix, subgroup order, 8.181 ratio, archived experiment path and degree-83 retained anchor.

- docs/index-calculus-scoreboard.html#n83-retained-gates-20261010 concerns the **a=0 81-bit** construction/storage/capacity study, explicitly engineering/L0 with complete runtime unknown. It does not report the a=1 online ratio.
- docs/ic/progress-timeline.json, docs/performance-gains/ and editable/text sources in figures/ yielded no relevant promotion of the old a=1/K600 8.181 observation under those exact terms.
- The canonical boundary ledger already labels the old observation exploratory, with controlled speedup unknown and historical verdict retained as provenance.

The canonical ledger already stores the correct legacy modulus and G/Q, without either erroneous curve slug; no identity join repair is needed there. Its controlled speedup, crossover and isolation fields remain null, and the historical label is explicitly qualified. No canonical graph or ledger is edited or promoted by this diagnostic. This scoped check does not exhaust historical branches or figures.
