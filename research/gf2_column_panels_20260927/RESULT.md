# Cache-sized column panels for GF(2) elimination

The initial column-panel candidate matched all 66 matrix-F5 calls exactly, including row fingerprints, rank, pruning and criterion counts. CI run [36305543942](https://github.com/aburan28/crypto/actions/runs/36305543942) at PR head `726e7486` completed on a pinned one-thread AMD EPYC 7763 Linux runner. Its [full raw receipt](runs/36305543942/ci-result.json) has SHA-256 `4fb312a904689dddc83c0e2221ee05ce49c0139c06ff15956ad993a9c48d582a`.

| `f5_n24_m24_d4` workload | Elimination ratio, 95% interval | Complete F5 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 0.457, 0.444–0.463 | 0.572, 0.558–0.579 |
| Holdout A | 0.442, 0.438–0.448 | 0.557, 0.554–0.562 |
| Holdout B | 0.443, 0.431–0.445 | 0.553, 0.542–0.555 |

The ratios are paired reference/candidate medians. Counted elimination word XORs were unchanged on the primary cell, so the new table replay and memory access cost outweighed its reduction in full-matrix sweeps. The default remains the original per-block kernel. A distinct compact, fused replay is frozen in the protocol for one bounded follow-up measurement; the initial regression is retained.
