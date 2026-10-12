# X86 guarded direct F5 output unpack: rejected after qualified paired runs

The [protocol](PROTOCOL.md) was committed as `14f5e9c37` before candidate code or timing. The representative [assembly probe](ASM_PROBE.md) motivated a separate x86 test after the Apple ARM64 bounds-free screen failed; it was static evidence only. [GitHub Actions run 36679082741](https://github.com/aburan28/crypto/actions/runs/36679082741) tested the opt-in implementation from PR head `1d8f231759e2decd323c79fd3a6c62dccedaa398` (synthetic merge `b7225cd0962ce3d8b73973a3a32414cd42248e00`). Both one- and two-thread jobs passed release shared-kernel and F5 tests, as recorded in [job metadata](RUN_CANDIDATE.json), then completed the frozen paired measurements on an AMD EPYC 7763, Linux x86-64, Rust 1.98.1, with AVX2/BMI2. Both modes used the same release binary, SHA-256 `963bd09ead62003d31c05d1bc8d5180e37365f83b529a0d4b1d4b3a8a884752c`.

The candidate F5 source had SHA-256 `bbb8597f1fe174e4e287647ac0a1783875b263b13dd231eec8edf57e4f9f792a`; the shared eliminator had `6fd3f10ef8e681490784c624e703eeb8fa887c5817ec9bdf34844e3f943c7f5c`. The benchmark source had `3c6be73582637c61b6b27780a9ef60cffcbedd18e6e263391fdd68c06aa28957`. The [measured source patch](measured_candidate.patch.gz) has gzip SHA-256 `d742f9581f42db55a4de01ab75bf47f84908f0235c22827c3124900fb7e77cff` and uncompressed SHA-256 `bb43677869bbeb7c8c8d986ac58b2c0fd90d385f645b89c6755e9f1ab97e8b5d`; the archived [workflow](WORKFLOW_CANDIDATE.yml) has `76a954304fd40e7552707b8d657aca760ec792c77f38e29e584fbac9f9e13f69`. The patch is relative to protocol commit `14f5e9c37` and contains the complete runtime and benchmark-output changes. Replay it with `gzip -dc measured_candidate.patch.gz | git apply` on that commit. The harness scripts remain in this directory.

Every selected seed block had 22/22 successful processes, seven cases per process, matching raw and canonical row fingerprints, rank, term count, build/pruning counts, and reduction word operations between arms. Both thread-count jobs selected the first qualified attempt for each seed. The one-thread job preserved three zero-call first attempts rejected by the preflight PSI limit, then selected attempt two for those holdouts; the two-thread job selected attempt one for all seeds. Every selected reservation left zero other eligible user threads and recorded zero contended samples. The frozen primary output was 6,924 rows, 12,951 columns, and 13,734,979 terms, with raw fingerprint `3f659516eff553b8` and row-space fingerprint `ed5234ba018bc079` in both modes.

## Complete-call decision

Ratios are medians of five paired reference/candidate complete-call ratios on the **same host and binary**; larger is better. The intervals are the protocol's exact 3,125-resample bootstrap 95% intervals. A/A minima come from five reference/reference pairs for each seed and thread count. These are complete matrix-F5 solver-call timings, including criterion, build, elimination, and output materialization.

| Seed | One-thread full ratio (95% interval) | One-thread A/A minimum | Two-thread full ratio (95% interval) | Two-thread A/A minimum |
| --- | ---: | ---: | ---: | ---: |
| Frozen `0` | **0.9984× (0.9927–1.0226×)** | 0.9864× | 0.9983× (0.9818–1.0092×) | 0.9816× |
| Holdout `badc0de1` | 1.0030× (0.9991–1.0070×) | 0.9914× | 0.9962× (0.9813–1.0184×) | 0.9921× |
| Holdout `5eed2026` | 1.0048× (1.0003–1.0151×) | 0.9937× | 0.9962× (0.9953–1.0082×) | 0.9919× |
| Holdout `f5c02a28` | 1.0034× (0.9948–1.0110×) | 0.9946× | 1.0066× (0.9996–1.0095×) | 0.9940× |

The frozen one-thread marginal complete-call medians were 106.680 ms reference and 106.723 ms candidate. Median unpack phases were 52.038 ms and 51.549 ms; the paired unpack ratio was only 1.0085×. These marginal medians are descriptive and are not the paired decision statistic. Absolute times from this EPYC 7763 are not compared with other hosts.

The predeclared advance gate required a frozen one-thread full-call ratio of at least 1.05× and a bootstrap lower bound above 1.02×. Both fail. Two smaller-case one-thread medians also fell below their own A/A minima: `f5_n12_m12_d4` on `badc0de1` was 0.99514× versus 0.99942×, and `f5_n24_m24_d3` on `f5c02a28` was 0.99996× versus 1.00107×. No two-thread smaller-case median fell below its A/A minimum. The candidate is **rejected**; no default-on confirmation is warranted. It does not approach the further 2× full-call target. The runtime option and one-off workflow were removed after archiving the tested bytes, leaving the original F5 source and benchmark source byte-identical to reference commit `bff733dc58c5a4f6bbc011607df5266e5b450cec`.

## Raw evidence

The deterministic [one-thread archive](runs/36679082741/MANIFEST-t1.json) contains 31 raw files, including all three PSI refusals and all selected calls. Its `segments-t1.tar.gz` SHA-256 is `8e464186fe194c78728da7eb43c4f9e26803c15578073c7729ebc58450ddf3d6`. The [two-thread archive](runs/36679082741/MANIFEST-t2.json) contains 22 raw files and its `segments-t2.tar.gz` SHA-256 is `cf076e92c1fdf5923eccf4d7cde7755e605635213d9d91c9e851dbb5f0a569ef`. Each manifest records every file's digest and the first-clean selection; [archive_receipts.py](archive_receipts.py) reproduces the archives from the downloaded Actions artifacts. Extract with `tar -xzf segments-t1.tar.gz` or `tar -xzf segments-t2.tar.gz` in an empty directory.

This is a solver-stage performance result. It establishes neither a one-target index-calculus DLP online speedup nor a Pollard-rho crossover.
