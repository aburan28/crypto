# Degree-83 IC against the strong rho reference

Exact retained a=1 / K600 fixture · completed five-pair native diagnostic

12 October 2026

> Observed IC/strong cost: 1.3823× prepared online and 218.0346× fresh-process wall from retained base. These are native paired geometric means on one Q and one host. The historical 8.1813× packed-reference observation does not establish a strong-reference advantage.

## 1. Diagnostic conclusion and the historical correction

> For this fixed point and frozen five-pair panel, the native IC/strong paired geometric mean is 1.3823 [1.0358, 1.8823] online and 218.0346 [172.3584, 285.5737] for fresh-process wall. Values above 1 mean greater IC cost. No controlled vs_rho row is promoted.

The source-bound comparison does not reproduce the previous synthesis’s interpretation of a degree-83 online lead against the required reference. Its 8.1813× median recorded rho/IC ratio used the old packed reference. That worker stores every step, scans the orbit in polynomial coordinates and inverts per step. Current reference rules exclude it as the primary strong control. The old scalar answers remain valid; its reference interpretation is corrected.

The fresh reference uses distinguished points, normal-basis canonicalization and batched inversion over 32 lockstep lanes in one CPU worker. Both fresh methods solve the exact same public Q without receiving its scalar. The standalone replay verifies every one of the 22 outcomes, including the two excluded warm-ups.

The intervals describe five seed/repeat pairs on this one Q and one host. They do not estimate performance on independently drawn targets, other hosts or deployed curves. The whole-process comparison begins with the retained base file and omits original base selection/construction; it is not a construction-from-curve total.

| Recorded estimator | Value | 95% native interval / n |
| --- | --- | --- |
| IC/strong prepared-online cost | 1.382330 | [1.035774, 1.882300] · 5 pairs |
| IC/strong fresh-process wall | 218.034569 | [172.358384, 285.573727] · 5 pairs |
| Strong A / strong B online control | 0.997036 | [0.967928, 1.023496] · 5 pairs |
| Archived packed/IC median | 8.181302 | Historical only; different reference and estimator |

*Sources: [R02](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/statistics.json), [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md).*

## 2. Exact field, subgroup, point and frozen design

The correct exact-field slug is `icv1-f2m83-t6151469093347-9eba6cd5`. Full identity: `ICV1:f2m-83-a5fd0d39:6151469093347:9671406556910881928556062:0x1:-7:unk:r:9eba6cd5d359`. The frozen protocol and original archive replay copied the wrong cdcc5432 alias, which belongs to low terms [0,1,2,45]. Their actual [0,2,4,7] field and G/Q arithmetic are unchanged. The separate erratum preserves frozen evidence.

Field degree 83 is distinct from the subgroup’s approximately 52.928 bits. This is the retained a=1 instance, not the a=0 degree-83 / 81-bit ECC2K confidence gate. Exact polynomial, subgroup and G/Q coordinates are the primary binding.

The panel requests two excluded warm-ups, five same-seed strong/strong A/A pairs, and five interleaved strong/IC pairs; all 22 children completed verified. Strong seeds are 202610061–202610065 and IC rank seeds 21–25. Matched order alternates by round. One fixed Q is repeated in every child.

| Item | Frozen value |
| --- | --- |
| Curve and field | y²+xy=x³+x²+1; GF(2)[z]/(z^83+z^7+z^4+z²+1) |
| Prime subgroup r / cofactor h | 8569786107849059 / 1128547018 |
| Generator G | (6614244651398184860292015,8868260410029249391539529) |
| Public Q | (2655898345537166180532960,3881526307991673265292615) |
| Retained base | 600 signed-Frobenius representatives; 99,600 usable points |
| IC workers | Frozen compact-orbit producer; rank threads 1, target threads 1 |
| Strong policy | Wide-field exact-point adapter; 32 lockstep lanes; DP bits 4; step-cap factor 2000; fruitless escape false |
| Envelope | CPU 3, one pinned CPU; 12 GiB RLIMIT_AS; 1800 s process-group timeout |

*Sources: [R01](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/protocol.json), [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md), [R06](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/CURVE_ID_ERRATUM.md), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md).*

## 3. Three spaces and a one-dimensional subgroup

The arithmetic field has 83 coordinates over F2. The retained base is a point-defined union of 600 signed-Frobenius orbits, not a small linear subspace of those field coordinates. Its relation labels are formal vectors in (Z/rZ)^600. The cyclic prime-order subgroup itself has one scalar dimension over Z/rZ; 600 columns do not create 600 independent group directions.

A normal basis changes field representation and makes Frobenius squaring a rotation; it does not change factor-base membership. The displayed four-vector witness is retained, not newly discovered. Independent replay checks the actual four-point sum and scalar answer; the full 600-row preparation transcript remains missing.

![Exact retained signed-orbit labels and four-vector witness, from INVARIANTS.md and the independent replay. The diagram describes existing algebra and supplies no runtime claim. [R04–R06]](figures/three_spaces.png)

Exact retained signed-orbit labels and four-vector witness, from INVARIANTS.md and the independent replay. The diagram describes existing algebra and supplies no runtime claim. [R04–R06]

*Sources: [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json), [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md), [R06](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/CURVE_ID_ERRATUM.md).*

## 4. The invariant shared by IC and strong rho

On the subgroup, τ(P)=[λ]P with λ=2998004823950457, λ²−λ+2=0 mod r and exact order 83. Signed Frobenius multiplies each attached scalar coefficient by ελ^k. IC records that transport in its factor-log columns; rho records it in two state scalars for P=[a]G+[b]Q. Both are formal labels in the same one-dimensional subgroup.

At a rho distinguished-point collision, the nonzero denominator b1−b2 permits d=(a2−a1)/(b1−b2) mod r. A zero denominator yields no recovered scalar. IC instead extracts a four-point witness and applies its factor-log dot product. Each final answer requires [d]G=Q replay. The shared invariant does not make their preparation, target work or runtime equivalent.

![Existing forward signed-Frobenius transport in the retained IC labeling and strong rho state. This is a method comparison, not a fresh identity or timing claim. [R05]](figures/coefficient_transport.png)

Existing forward signed-Frobenius transport in the retained IC labeling and strong rho state. This is a method comparison, not a fresh identity or timing claim. [R05]

*Sources: [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md).*

## 5. Per-round recorded online and process values

Online values are each producer’s complete declared target interval. The primary statistic does not replace those intervals with a sum of component clocks. Every matched row solves Q; there are no exclusions, failed answers, overall caps or missing matched receipts.

Strong online includes Q membership/subgroup validation, target-dependent starts, walk/collision work and generic scalar replay. The IC target stage assumes the retained fixture is prevalidated; that declared validation asymmetry remains visible.

M4 has a lower measured IC online interval: 6204.636907 ms versus strong 6698.575519 ms, native strong/IC ratio 1.079608. M5 is near a tie, native ratio 0.999600. The conditional paired geometric mean over all five rounds still gives greater IC cost; this does not assert that strong is faster in every seed round.

| Pair / order | Strong online ms | IC online ms | Strong/IC native ratio |
| --- | --- | --- | --- |
| M1 · rho→IC | 4516.159970 | 6437.368839 | 0.701554 |
| M2 · IC→rho | 4279.021477 | 6416.446390 | 0.666883 |
| M3 · rho→IC | 2514.613847 | 6408.147497 | 0.392409 |
| M4 · IC→rho | 6698.575519 | 6204.636907 | 1.079608 |
| M5 · rho→IC | 6176.824943 | 6179.295727 | 0.999600 |

![Recorded intervals for run-20261012T0226Z, n=5 seed pairs on one Q. Right: fresh-process wall from retained base, logarithmic axis; includes measured reusable preparation. [R02,R03]](figures/paired_runtime.png)

Recorded intervals for run-20261012T0226Z, n=5 seed pairs on one Q. Right: fresh-process wall from retained base, logarithmic axis; includes measured reusable preparation. [R02,R03]

*Sources: [R02](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/statistics.json), [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json).*

## 6. Native uncertainty and same-seed A/A controls

The frozen estimator takes the arithmetic mean of within-pair log(strong/IC) ratios and exponentiates it. It uses 20,000 paired percentile-bootstrap resamples. The runner applies seed 2026101283 to A/A and seed+1=2026101284 to matched pairs. The native secondary exporter presents IC/strong as the reciprocal, with interval [1/upper,1/lower]; it performs no new resampling or estimator change.

These geometric means are not ratios of arithmetic mean times. The A/A control’s paired geometric mean is 0.997036, interval [0.967928,1.023496], with observed online ratios spanning 0.948233–1.026892. It diagnoses repeat noise at fixed seeds on this host; it does not establish physical-host exclusivity.

n=5 is a small conditional panel. Seed sensitivity and repeated measurement of one Q cannot supply a fresh-target extraction-tail distribution. Keep all raw pairs, avoid population-wide significance claims, and prospectively freeze independent unseen targets for the next question.

![Native paired geometric means and recorded 95% intervals, run-20261012T0226Z. IC/strong >1 means greater IC cost; A/A is strong-A/strong-B. Each interval covers five pairs on one Q and one host. [R02,R03]](figures/paired_intervals.png)

Native paired geometric means and recorded 95% intervals, run-20261012T0226Z. IC/strong >1 means greater IC cost; A/A is strong-A/strong-B. Each interval covers five pairs on one Q and one host. [R02,R03]

*Sources: [R01](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/protocol.json), [R02](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/statistics.json), [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R10](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/src/main.rs).*

## 7. Fresh-process cost from the retained base

The outer process starts from the retained factor-base file. IC loads its point arrays, builds the root index and acquires/solves its guided rank/log preparation. It does not newly select/build the original 600 representatives. OS page cache is not flushed. Thus cold means fresh process-local state, not OS-cold or total construction-from-curve runtime. Adding original construction cost would add IC work, but its unmeasured value stays unknown.

Strong curve/normal-basis construction and G-only jumps are measured. Supplied subgroup/catalog facts are parameters, not newly timed curve discovery or factorization. The native means below are secondary descriptive arithmetic means and ranges over five successful formal runs; they replace neither the paired primary estimator nor its interval.

| Metric: mean [min,max] | Strong | IC |
| --- | --- | --- |
| Preparation ms | 191.516 [186.924, 194.744] | 1047924.983 [1009528.197, 1139341.566] |
| Process wall s | 5.068 [2.804, 6.909] | 1054.552 [1015.960, 1145.812] |
| Child user+system CPU s | 5.029 [2.709, 6.887] | 1054.438 [1015.897, 1145.722] |
| Outer peak RSS bytes | 116588544.0 [55263232.0, 208527360.0] | 5811525222.4 [5811470336.0, 5811605504.0] |

- IC prepared-online extraction is measured after its reusable rank/log preparation. Fresh cold children reconstruct index/log state independently; no cross-process warm cache is inherited.
- The retained root table has 29,878,182 entries. Table-entry counts are stored states, not independently verified factor logs or useful-row counts.
- The root/index and guided preparation must stay charged in this workload. Amortized reuse across new targets requires an explicitly different workload count.

*Sources: [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md).*

## 8. Phase residuals and the outer wall boundary

IC PDP subtracts its nested final-lift check from extraction wall at one target thread. That check is only the final addition/equality; earlier lift work stays under PDP. All five formal IC rows have exact five-phase conformance=false. Strong’s three component timers also leave nonzero residuals; its five-phase flag is inapplicable/null. The native primary uses complete online intervals.

Outer wall starts before spawn and ends when wait4(WNOHANG) observes exit. The runner polls every 100 ms: wall includes nominal 0–100 ms detection delay plus scheduling/poll/heartbeat effects, not a strict 100 ms bound. Child CPU excludes parent sleeping. Formal M3 strong reports producer process_before_receipt_io=2.707912963 s, CPU=2.708705 s and parent wall=2.803706954 s. These boundaries do not prove an exact causal overhead decomposition. Frozen wall values remain untouched.

IC retains five component clocks. Each residual below is the native complete interval minus its reported phase sum; no phase sum is silently substituted for online.

| Pair | Query ms | PDP ms | Lift ms | Descent ms | Replay ms | Residual ms |
| --- | --- | --- | --- | --- | --- | --- |
| M1 | 0.010927 | 6437.211479 | 0.012360 | 0.000611 | 0.103506 | 0.029956 |
| M2 | 0.000030 | 6416.329133 | 0.012258 | 0.000640 | 0.103356 | 0.000973 |
| M3 | 0.000020 | 6408.030341 | 0.012148 | 0.000601 | 0.103426 | 0.000961 |
| M4 | 0.000020 | 6204.505849 | 0.025700 | 0.000630 | 0.103586 | 0.001122 |
| M5 | 0.000071 | 6179.177809 | 0.012648 | 0.000831 | 0.103577 | 0.000791 |

Strong reports Q query/validation, walk/collision and generic replay. The different component definitions are preserved.

| Pair | Q validation ms | Walk/collision ms | Generic replay ms | Residual ms |
| --- | --- | --- | --- | --- |
| M1 | 0.232440 | 4515.698654 | 0.228475 | 0.000401 |
| M2 | 0.232550 | 4278.561313 | 0.227333 | 0.000281 |
| M3 | 0.248956 | 2514.137047 | 0.227503 | 0.000341 |
| M4 | 0.261765 | 6698.084148 | 0.229306 | 0.000300 |
| M5 | 0.233803 | 6176.356336 | 0.234494 | 0.000310 |

Both native producers have nonzero timer residuals. The complete-interval primary and descriptive process metrics remain available as diagnostic observations; exact phase-accounting admission remains unmet.

The parent observer-clock caveat is retained: before-spawn→observed-exit includes polling and scheduling/heartbeat effects. CPU and producer process timestamps corroborate it without retroactively adjusting wall. Future harnesses should use event-driven pidfd/SIGCHLD or native blocking wait timing.

*Sources: [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md), [R10](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/src/main.rs).*

## 9. Counters and independent correctness at their actual scope

> Independent answer/base/witness correctness is established at the stated replay scope. It does not certify every missing preparation equation, reproducibility of wall time, physical-host isolation or a speedup across independent targets.

IC returns the same four indices [87298,27061,36601,66281] and 8,845,441 target probes in each formal repeat of Q. Probe units count partner candidates reaching canonicalization/lookup, not every S3 call or group addition. The producer records rank=600, 600 admitted rank relations and zero rank failures; the complete rank equations and factor-log vector are not exported.

The standalone C++ replay revalidates the field/curve/subgroup, all 600 representatives and 99,600 point/label records, every actual IC four-point witness and all 22 recovered answers. It reports PASS, 22 independently verified answers and zero invalid rows. The validation-only scalar sidecar is checked after answers, not passed to either solver.

Rho steps and IC probes are different counted units. The strong lane counters below describe internal lane retirement and cycle handling; nonzero local capped_walks do not mean a process failed or reached the overall cap. Common calibrated operation-equivalent S stays null.

| Pair | IC probes | Rho walk steps | Local capped lanes | Fruitless cycles |
| --- | --- | --- | --- | --- |
| M1 | 8845441 | 10203425 | 148 | 1845 |
| M2 | 8845441 | 9436758 | 147 | 1699 |
| M3 | 8845441 | 5353957 | 87 | 973 |
| M4 | 8845441 | 14229982 | 234 | 2510 |
| M5 | 8845441 | 13159092 | 187 | 2319 |

*Sources: [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json), [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md).*

## 10. Host, validation and admission boundaries

The host receipt describes virtualized KVM Linux x86_64, kernel 6.18.44, AMD EPYC 9V74 exposed model, five guest logical cores and one NUMA node. Cgroup RAM is 16 GiB with CPU quota 4. Each child has CPU-3 affinity and a 12 GiB virtual-address cap. Native isolation finishes exit 0 with contended_samples=0 and threads_moved=29. These are guest monitoring observations, not physical-host exclusivity or sealed ecbench L2.

All 22 children completed verified: two excluded warm-ups, ten A/A children and ten matched children. There are no process failures, overall caps or excluded pairs. Fresh independent arithmetic replay verifies all 22 answers. Final default-feature cargo test --release --locked --lib passes 3947 tests with zero failures and 112 ignored (260.41 s). Delivery runner tests pass 7/7, native exporter 4/4, and final combined example tests pass IC 10/10 plus strong 2/2. These validation checks are separate from timed execution.

The frozen baseline IC executable predates the isolated discarded-RNG-draw collector repair. The strong adapter was compiled after it, but neither measured solver calls that generic collector path. The repair is a reproducibility fix outside this measured algorithm, not an IC optimization or a reason to rewrite frozen timing receipts.

The exact legacy wide-field IC adapter is absent from the current sealed ecbench path. Canonical candidate/workload identities remain null, common S is uncalibrated, the full rank/log transcript is missing and there is no second-host timing/replay admission. Conductor invocation returned command-not-found, and no tracked .conductor/WORKFLOW.md was available. These gaps remain explicit; this report does not promote vs_rho.

- Correct the old packed-reference interpretation and wrong-modulus alias in new prose; preserve old frozen receipts.
- Next measure independent unseen targets with a prospectively frozen preparation/selector policy; five seed repeats of Q do not resolve target tails.
- For full algebraic preparation certification, retain every admitted rank relation, its witness and setup scalar, rank/matrix evidence, the solved factor logs and independent [L_j]G=R_j replays.
- For a construction-from-curve objective, newly measure original representative selection/construction; for admission, add exact ecbench integration, calibrated common costs and applicable isolation/second-host gates.
- Use event-driven child completion in a future runner, while retaining CPU totals and producer timestamps; leave this frozen panel’s wall unchanged.

*Sources: [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md), [R08](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/README.md), [R09](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/README.md), [R10](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/src/main.rs). Validation logs: validation_logs/full_lib_tests_final.log and validation_logs/examples_tests_final.log.*

## 11. Evidence register and reproducible report source

This is an accounting/reference correction and a source-bound diagnostic of the existing retained candidate. It introduces no selector, algebraic identity, asymptotic improvement or admitted performance row. The quantitative plots lay out native recorded values and intervals; Python performs no statistical estimation, arithmetic replay or benchmark computation. The authoritative executed identity is records/frozen_panel_inputs.sha256, checked unchanged after the panel and corroborated by every child start/execution receipt. Earlier frozen_before_execution hashes are preparation history only and are not used as executed bindings.

The report lives beside the raw records at research/n83_evaluation_20261012/reports/run-20261012T0226Z. Native JSON copies under data/ preserve the recorded source outputs; figures.py, editable SVG diagrams, report_content.json and report_source.md reproduce its layout.

Canonical graph review found the n83 retained scoreboard anchor concerns the distinct a=0 / 81-bit engineering capacity study. The boundary ledger already labels the old a=1 observation exploratory; no controlled frontier or ledger row is promoted here. The scope of graph/source checks is recorded in the study README.

| Evidence / source pin | SHA-256 or exact reference |
| --- | --- |
| Native primary statistics | 43b00e963b379aaf718a16b864a3125a866a8be9f7f83b8eaa109d23c3d1560e |
| Native secondary export | b8101e76589558c8329831e3fe4649ea010661a34c6653bdc2fec40b5063b0b3 |
| Independent fresh replay | 980a3d7d52db96a559dd36ad12e317c9b12be9f1c161a2d0d86803aafcc72441 |
| Paired-record input | e1f1b019ed5f80a739b0acc9c2e7f971d8e3400f113b47700067b0a64f8c0215 |
| Frozen IC executable | 983f28400b7852fc6628ba97d3ea83b2ab9e32d9944658c7eea831ee89a4abf2 |
| Frozen strong executable | 98949b7127baeb1fa55420d23dc13604069fcc8fde8d92b8f2e6bebd09718542 |
| Executed native runner | 8da1b77f9521fccdc5d8490a493146f0a6e2af92a2fab07c4e2e374954b9d2b4 |
| Executed input protocol | ac7266f6caad2a4142ae06884ee748f980f9ee200c4e2649a0149f59e91adf29 |

Native source hashes above identify the copied recorded outputs and executables. Links below lead to the durable study records and source explanations.

| ID | Evidence / path |
| --- | --- |
| [R01](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/protocol.json) | Frozen protocol: runner/protocol.json |
| [R02](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/statistics.json) | Native paired statistics: records/panel/statistics.json |
| [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json) | Native secondary export: records/panel/report_data_secondary.json |
| [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json) | Independent fresh replay: verification/fresh_replay_receipt.json |
| [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md) | Three spaces and signed-orbit witness: verification/INVARIANTS.md |
| [R06](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/CURVE_ID_ERRATUM.md) | Curve-identity erratum: verification/CURVE_ID_ERRATUM.md |
| [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md) | Source audit: source_audit/AUDIT.md |
| [R08](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/README.md) | Host and isolation records: records/README.md |
| [R09](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/README.md) | Study delivery and validation: README.md |
| [R10](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/src/main.rs) | Native runner source: runner/src/main.rs |
| [R11](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/descriptive_summary/src/main.rs) | Native secondary exporter: descriptive_summary/src/main.rs |

*Sources: [R01](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/protocol.json), [R02](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/statistics.json), [R03](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/panel/report_data_secondary.json), [R04](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/fresh_replay_receipt.json), [R05](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/INVARIANTS.md), [R06](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/verification/CURVE_ID_ERRATUM.md), [R07](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/source_audit/AUDIT.md), [R08](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/records/README.md), [R09](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/README.md), [R10](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/runner/src/main.rs), [R11](https://github.com/aburan28/crypto/blob/codex/n83-strong-evaluation-20261012/research/n83_evaluation_20261012/descriptive_summary/src/main.rs).*
