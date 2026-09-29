# Prime F4 lazy AVX-512 experiment

[CI run 36471302637](https://github.com/aburan28/crypto/actions/runs/36471302637) completed the frozen one-thread AVX-512 experiment on an Intel Xeon Platinum 8573C with AVX2, AVX-512F and BMI2. The direct AVX-512 accumulator boundary test, all prime F4 tests and all 132 benchmark processes passed. Across the frozen seed and three holdouts, every wide, current and candidate call matched the reference's basis and solve fingerprints, step count, matrix shape and basis size on all 13 cases. The runner pinned one CPU and one Rayon thread. Both release binaries are retained in the CI artifact. The [compressed full receipt](runs/36471302637-t1.json.gz) has SHA-256 `29c284e67cdbb19abaf597780c3d3996ae194460bce27fa76ca64e514a76c652`; the uncompressed JSON SHA-256 is `b2416d2c9afd2bda08edd06808623f63c8fc99f9a7e9eb784a9ccd86354202a1`. It retains every process output/status, source and binary hashes, CPU features, affinity, load and per-case timing.

The primary case is `quad_n8_p65521`. Ratios below are medians of five paired reference/candidate complete-call ratios with exact five-pair bootstrap 95% intervals. Current is the existing deferred 32-bit scalar default; wide is the original 64-bit arithmetic. On the frozen seed, marginal median times were 2161.31 ms wide, 1497.17 ms current and 1429.91 ms candidate. Paired ratios need not equal ratios of marginal medians.

| Seed workload | Current/candidate | 95% interval | Wide/candidate | 95% interval | A/A range |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen | **1.047×** | 1.030–1.048× | **1.513×** | 1.503–1.524× | 0.993–1.000× |
| Holdout A | 1.041× | 1.040–1.042× | 1.509× | 1.503–1.511× | 0.998–1.003× |
| Holdout B | 1.047× | 1.041–1.052× | 1.510× | 1.502–1.512× | 0.993–1.003× |
| Fresh holdout C | 1.044× | 1.041–1.045× | 1.511× | 1.507–1.515× | 0.993–1.004× |

All primary incremental holdouts are above their A/A maxima, and no smaller case on any seed has a current/candidate complete-call median below its own A/A minimum. The frozen incremental median misses the predeclared 1.12× promotion gate; the direct wide/candidate ratio misses the user's 2× complete-call target. The scalar kernel remains the default. `F4_FP_LAZY_AVX512=1` remains an explicit hardware-gated research option, with no four-thread control or timing-driven repeat. The one-off timing workflow was removed after archiving this receipt; the normal F4 and repository CI checks passed on the measured code head.

[CI run 36469467083, attempt 1](https://github.com/aburan28/crypto/actions/runs/36469467083/attempts/1) passed the AVX-512 boundary and prime F4 test step and built both binaries, but its AMD EPYC 7763 runner lacked AVX-512F. The benchmark recorded `unsupported_host` and zero process calls. The complete [zero-call receipt](runs/36469467083-attempt1-unsupported.json.gz) has compressed SHA-256 `85fd394c082154fff951848dd23483754bb0f7764e05c330850bfc5c70bfec28` and uncompressed SHA-256 `ca205921d1fac42ccbc3001b5709dcc943de6fd2a84797e1582295307de86089`. No timing sample was discarded.

Attempt 2 of the same run was canceled during the exactness-test step, before either binary build or any benchmark process, for the pre-measurement AVX2 dispatch-guard amendment recorded in the protocol. It produced no timing sample.

These measurements cover complete prime-field F4 solver calls only. They do not establish a 2× primary complete-call gain, one-target IC online time or DLP speedup. The 1.756× AVX2 result from [PR #903](https://github.com/aburan28/crypto/pull/903) was measured on a different AMD host and is not multiplied or directly compared with this Intel-host result.
