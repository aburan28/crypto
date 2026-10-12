# Combined F5 options: direct complete-call experiment

The frozen same-binary comparison used the current matrix-F5 default as baseline and enabled selective echelon output, fused row construction and AVX2 table XOR together as the candidate. [CI run 36467065014](https://github.com/aburan28/crypto/actions/runs/36467065014) completed all 88 processes across one frozen seed and three holdouts. Every call succeeded and matched the baseline's canonical row-space fingerprint, rank, pruned-row count and criterion word operations on all seven cases. Raw row fingerprints were stable within each arm. On the primary `f5_n24_m24_d4` case, raw fingerprint changed from `ed5234ba018bc079` to `3f659516eff553b8` as expected for echelon output, while canonical row space remained `ed5234ba018bc079`.

The pinned one-thread Linux host was AMD EPYC 9V74, with AVX2, AVX-512F and BMI2. Thus baseline used its normal AVX-512 dispatch and candidate explicitly used AVX2. The [compressed complete receipt](runs/36467065014-t1.json.gz) has SHA-256 `b94164579b76783d7b2bfdf5910c565b29bb520b814af3e38ab108e18b6b379e`; the uncompressed JSON SHA-256 is `dbe7c56a2f6b3eb9abc79f26d79fae97e52a5d974fac1c2de8f2bec3ab4a3cd3`. The CI artifact retains the release binary. Source/binary digests, host features and load, CPU affinity, every process output and status, phase timings and per-case signatures are in the receipt.

Ratios are medians of five paired baseline/candidate complete-call ratios with exact five-pair bootstrap 95% intervals. The A/A range is the five baseline/baseline control ratios.

| Primary seed | Complete-call ratio | 95% interval | A/A range | Build ratio | Reduction ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen | **1.530×** | 1.528–1.546× | 0.997–1.007× | 1.353× | 2.225× |
| Holdout A | 1.514× | 1.503–1.519× | 0.995–1.012× | 1.328× | 2.199× |
| Holdout B | 1.552× | 1.547–1.559× | 0.994–1.008× | 1.341× | 2.234× |
| Fresh holdout C | 1.547× | 1.544–1.566× | 0.999–1.008× | 1.344× | 2.266× |

On the frozen primary, the marginal median full-call times were 198.55 ms baseline and 129.61 ms candidate. Median phase times fell from 32.37 to 23.99 ms for build and 107.68 to 48.39 ms for reduction; unpacking remained 57.01 versus 55.65 ms. No smaller case on any seed had a complete-call median below its own A/A minimum. The measured combined options nonetheless miss the predeclared 2× complete-call gate, so this experiment stops without a four-thread control or timing-driven repeat. The default remains unchanged; the three options remain explicit opt-ins with their individual receipts and this combined receipt. This AVX-512-host result does not establish a ratio on an AVX2-only host.

[Run 36465188507](https://github.com/aburan28/crypto/actions/runs/36465188507) attempts 1 and 2 both passed the existing F5 and GF2 tests and built the benchmark, then stopped before any benchmark process because the hosts had AVX-512F. Their respective zero-call receipts are [attempt 1](runs/36465188507-attempt1-unsupported.json.gz) and [attempt 2](runs/36465188507-attempt2-unsupported.json.gz), with compressed SHA-256 `561d3a7b278585f6899b8a82aca6c6fbe93ae2c8686442b5e679450922631bbf` and `b9dbb8450d9fc43b9f4a7d080dd1e17cfc3aad46cf58145b479104da76687fa0`. The protocol was amended before measurement to accept any AVX2 host and compare against that host's true default dispatch. Neither attempt supplied a timing sample.

This is a complete matrix-F5 solver-call measurement only. It does not establish a 2× F5 call, a one-target IC online result or a DLP speedup. The one-off CI workflow was removed after the receipt was archived.
