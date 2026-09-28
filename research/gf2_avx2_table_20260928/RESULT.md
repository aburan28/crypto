# GF(2) AVX2 table-XOR experiment

The initial CI run [36452894497](https://github.com/aburan28/crypto/actions/runs/36452894497) was canceled while queued by the pre-measurement host amendment; it made no benchmark call. The completed [run 36454200069](https://github.com/aburan28/crypto/actions/runs/36454200069) at PR head `7d5578db` passed all four exact shared-kernel tests with AVX2 forced, built the unmodified reference and candidate, and completed all 88 process calls. Every case on all four seed workloads matched raw and canonical row fingerprints, rank, pruning, criterion word operations and elimination word operations exactly. The one-CPU pinned runner was AMD EPYC 7763 with AVX2 and BMI2 but no AVX-512F; Rust and both binary digests are in the raw receipt.

The [complete compressed receipt](runs/36454200069-t1.json.gz) has SHA-256 `65a3966ac8ff69575f7f58f3d97212fe49c8e733c1e349d12c821714b2f56405`; `gzip -dc research/gf2_avx2_table_20260928/runs/36454200069-t1.json.gz` extracts it. Its uncompressed SHA-256 is `7e26fbaaa8df28bd5c84272b76772218ebfc33ef0a64b9fd9d270e879135a8e8`. The CI artifact retains both binaries. The receipt includes every call, output, status, phase timing, host load, affinity, source hash and binary hash.

| Primary `f5_n24_m24_d4` seed | Elimination ratio, 95% interval | Complete-call ratio, 95% interval | Elimination A/A range |
| --- | ---: | ---: | ---: |
| Frozen | **1.151×, 1.132–1.160×** | 1.088×, 1.073–1.090× | 0.996–1.006× |
| Holdout A | 1.141×, 1.115–1.146× | 1.078×, 1.056–1.084× | 0.992–1.019× |
| Holdout B | 1.143×, 1.138–1.153× | 1.081×, 1.079–1.089× | 0.983–1.000× |
| Fresh holdout C | 1.155×, 1.149–1.169× | 1.092×, 1.089–1.097× | 0.992–1.001× |

Ratios are medians of five paired reference/candidate ratios. On the frozen primary case, marginal median elimination time fell from 153.75 to 133.60 ms, while complete-call time fell from 259.52 to 238.71 ms. All holdout primary elimination gains exceed their A/A maxima, and no smaller complete-call median falls below its own A/A minimum. The frozen primary elimination ratio is below the preregistered 1.20× gate, however, so the candidate **does not become the default** and no four-thread control is run. The AVX2 implementation remains available only with `KIC_GF2_FORCE_AVX2=1` (and `KIC_GF2_SIMD` not set to `0`) for research or later combinations. On an AVX2 machine without that flag, the original scalar path remains selected; AVX-512 machines retain their existing default. This is a solver-stage measurement, not a 2× complete F5 call or a one-target IC/DLP speedup.
