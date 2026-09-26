# n53 combined L384: one index and one rho table on the same 384 Q

The [pre-outcome protocol](PROTOCOL.md) committed the exact 384 distinct public Q values, one fresh 512-target/220-column compact-orbit training solve, one cold point-query index, and one fresh signed-Frobenius Kuhn–Struik rho table. The fixed Q file is the byte-exact concatenation of the three disjoint L128 blocks from [PR #754](../autolab_shared_log_n53_20260925/RESULT.md). [Run 36127210063](https://github.com/aburan28/crypto/actions/runs/36127210063) completed every stage on one GitHub-hosted Ubuntu x86_64 runner. Independent group-law audit checked all 512 training relations, all 23,320 orbit labels, 220 base logs, 384/384 IC point relations/logs and 384/384 rho scalars against those same Q coordinates and withheld labels. No final-run query failed, timed out, or crossed the 2 GiB process-group RSS limit.

For this exact stream and pre-certified base, **fully charged operational IC took 353.649 s versus 154.799 s for one-table rho**, a 2.285× IC/rho wall ratio. Even the incomplete IC two-Rust-child lower cost was 276.592 s, or 1.787× rho. This satisfies the preregistered fixed-stream no-crossover condition. It resolves the previous three-block L128 portfolio's straddling interval for this *combined one-table protocol*; it is not a distributional break-even, an n=131/ECC2K-130 result, or an attack-speed comparison in common operation units.

## Frozen inputs, source and rank

The Q input SHA-256 is `d5185187014a12516aeef29306b65d4d20864293bb2c60667a9684fa8e51ac97`; the withheld validator manifest SHA-256 is `f1843670a169d65645bf83886aeffee2e60e25627faa0002763eb1e7c13d8362`. The certified 220-orbit/23,320-point base has point-set SHA-256 `d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71` and gzip SHA-256 `23397af2ef668aed0775bcb409e1ae19555357ded635452818c9a3812f679d08`. The compact producer and rho source SHA-256 values were `c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09` and `fedacb54e441979c8c32860e7b5639e43799234741d49677010164564536d2c8`; the exact built binary SHA-256 values are in the archived panel. [SOURCE_FREEZE.json](SOURCE_FREEZE.json) pins all Python and Rust inputs from before any L384 outcome.

The fresh compact training producer completed 512/512 targets. Operational sparse elimination reached rank 220 at relation 461, predicted the remaining 51 training scalars from its solved base logs, and independently checked all 220 representative logs. A separate full-convolution group-law audit reproduced the same rank transition, base logs and all 512 group witnesses. The point and rho children received only the public Q file, never validator scalar labels; operational point recovery used only those Q values, fresh training logs and its own four-point lifts.

## Whole-process costs and native work

The first four IC rows are mandatory sequential work: producer training, operational row reconstruction/rank solve, one point-index process, and scalar recovery. Independent audit is separately timed and excluded from *both* primary operational costs. The lower bound omits mandatory rank and recovery; it is not a recovered-log cost. Compilation, creation of the fixed public-Q file, and construction of the already certified base are outside these process timings.

| Arm or phase | Cold wall s | CPU user+system s | Peak resident MiB | Result |
| --- | ---: | ---: | ---: | --- |
| IC training producer | 169.876 | 172.406 | 396.6 | 512/512 relations |
| IC mandatory rank/solve | 49.765 | 49.753 | 38.9 | rank 220; 220 logs |
| IC one-index point producer | 106.717 | 109.304 | 394.1 | 384/384 relations |
| IC mandatory scalar recovery | 27.291 | 27.281 | 38.2 | 384/384 logs |
| **IC complete operational** | **353.649** | **358.745** | **396.6** | **2.285× rho** |
| IC two-child lower, incomplete | 276.592 | — | — | 1.787× rho |
| IC audit-inclusive diagnostic | 504.545 | — | — | 3.259× rho |
| **One-table signed-Frobenius rho** | **154.799** | **154.773** | **76.8** | **384/384 logs** |
| Separate independent IC/rho audit | 150.896 | 150.803 | 43.4 | all checks pass |

IC training made 23,689,479 exact S3 calls; the single L384 point process made 14,931,236, for 38,620,715 total. Each fresh compact process scanned 2,565,200 regular states and built 5,081,560 index entries. The L384 point process returned 29,862,297 partner roots. Rho took 11,485,626 walk steps and 12,255,749 group additions, built 765,929 table entries, and recorded 382 cross-target solves. The complete archive preserves per-query S3, partner-root, indexed-hit, group-lift, hit/miss and rho charge counters. S3 and rho group additions are different units: no measured conversion to a common group-addition-equivalent `S` was frozen, so `S` and an attack-speed crossover remain unset.

The matched runner's training load average was 3.97 before its first process and 1.42 afterward; the rho/IC point stages began at 1.26/1.21. Operational CPU also exceeded rho by 2.318×, close to the 2.285× wall ratio. These are fixed-host observations, not a seed-distribution interval. The panel ended at 659.421 s, below its 2,700 s global ceiling. A review found that the runner checks that ceiling at stage entry rather than interrupting a stage already running; the observed completion is far below the ceiling, so this implementation limit did not affect this result.

## Durable replay, retained checker failure and decision

The 30-file [raw evidence manifest](evidence/archive_manifest.json) records archive SHA-256 `bef32e5d3993fe0a6cb5cc6be7c92b8b094cbaa24ff01031400d3f7e3df583b6`; all stage manifests, output SHA-256 values, CPU, RSS, load, commands, source/input/binary hashes, full JSONL and stderr remain in the archive. The measured workflow completed successfully in [run 36127210063](https://github.com/aburan28/crypto/actions/runs/36127210063).

An initial archive-only replay after merging main failed before opening any relation because its verifier applied a 64-hex SHA-256 predicate to a 40-hex Git checkout commit. [Run 36128714271](https://github.com/aburan28/crypto/actions/runs/36128714271), the byte-exact [original verifier](verifier_amendment_1/verify_archive_original.py), and the [amendment manifest](verifier_amendment_1/VALIDATION_AMENDMENT.json) preserve that failure and the checker-only correction. The pre-outcome measured-source hash still refers to the original verifier; the final replay verifier is separately hashed and pinned. [Corrected archive-only CI run 36129022158](https://github.com/aburan28/crypto/actions/runs/36129022158) rehashed the immutable archive and freshly replayed all 512 training relations, 384 IC logs and 384 rho logs (`PASS_COMPLETE_INDEPENDENT_REPLAY`). The targets, Rust/Python producer, operational solve and audit, resource limits, raw archive and measured children did not change.

The next decision is to gate further compact-index tuning on a preregistered model that charges base construction, root support, failed-query S3 work and a calibrated common operation unit at the proposed rung. This n53 fixed-stream negative makes another larger-batch run unjustified without a measured reason to expect the total-cost ratio to fall. The preregistered [rotated implicit-support gate in PR #762](https://github.com/aburan28/crypto/pull/762) is the concrete next experiment; its outcome is pending. For ECC2K-130, the separately established n=131 materialized-root coverage/memory admission failure remains the architectural obstacle; an implicit low-memory producer must first show enough supported target mass for independently verified rank and point-log recovery. No measured n=131 or Certicom speedup follows from this toy control.
