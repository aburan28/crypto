# F5 fused row-build experiment

CI run [36450865831](https://github.com/aburan28/crypto/actions/runs/36450865831) at PR head `44d5d288` passed the exact F5 row and cap tests but stopped before any benchmark call. The sparse checkout used to build the unmodified reference omitted `docs/ic/calibration.json`, which `src/cryptanalysis/ic_boundary.rs` includes at compile time. No reference or candidate timing was produced; this is a setup failure, not a performance result. The workflow now includes `docs/ic` in that checkout. The original run log and failed job remain available at the linked CI run.

## Frozen one-thread result

[CI run 36451830641](https://github.com/aburan28/crypto/actions/runs/36451830641) at head `08b98155f6ed653011647956990d417c01ca0d26` completed all 88 process calls (four seed workloads, two warmups and ten A/A plus ten A/B calls each). Every call returned successfully; raw and canonical row fingerprints, rank, pruned-row count, criterion word operations and elimination word operations matched exactly across arms. The exact row/cap unit tests passed. The host was AMD EPYC 9V45 with AVX2, AVX-512F and BMI2; affinity pinned one CPU and Rayon used one thread. The unmodified reference and candidate binaries are retained in the CI artifact.

The raw receipt is [36451830641-t1.json.gz](runs/36451830641-t1.json.gz), SHA-256 `7eb75b36bb98cab6fef624c1d6cd4cc33f8e93964e375dac4e5989a1fbdc5d93`. Extract with `gzip -dc research/f5_row_build_20260928/runs/36451830641-t1.json.gz`. It contains full output, per-call status, source and binary digests, timing, load and host information. The original uncompressed receipt SHA-256 is `73180d8d0bb007eb91d01a48e66ebac51caa05d5dbc458e3d1b9b2e10b349961`.

For the primary `f5_n24_m24_d4` case on the frozen seed, median reference and candidate build times were 26.04 and 17.24 ms. Median complete-call times were 158.97 and 146.84 ms. Ratios below are the medians of the five paired reference/candidate ratios; intervals are exact five-pair bootstrap 95% intervals.

| Seed workload | Build ratio | Build interval | Complete-call ratio | Complete-call interval | A/A build range |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen | 1.522× | 1.468–1.569× | 1.074× | 1.049–1.111× | 0.996–1.023× |
| Holdout A | 1.469× | 1.464–1.496× | 1.052× | 1.039–1.068× | 0.979–1.018× |
| Holdout B | 1.522× | 1.454–1.541× | 1.073× | 1.019–1.091× | 0.999–1.029× |
| Fresh holdout C | 1.539× | 1.493–1.558× | 1.082× | 1.059–1.115× | 0.989–1.002× |

The primary build and full-call gates passed. Every holdout primary build ratio was above its A/A maximum. Across all four seeds and the six smaller cases, no complete-call median fell below its own A/A minimum. The frozen smaller-case complete-call ratios ranged from 1.076× to 1.276×. This is a solver-stage measurement; it is not a 2× complete-call result or an end-to-end IC speedup.

The four-thread control specified in the protocol is pending in the next CI run. The workflow change for that control does not alter the candidate source or benchmark harness.
