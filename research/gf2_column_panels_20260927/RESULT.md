# Cache-sized column panels for GF(2) elimination

The initial column-panel candidate matched all 66 matrix-F5 calls exactly, including row fingerprints, rank, pruning and criterion counts. CI run [36305543942](https://github.com/aburan28/crypto/actions/runs/36305543942) at PR head `726e7486` completed on a pinned one-thread AMD EPYC 7763 Linux runner. Its [full raw receipt](runs/36305543942/ci-result.json) has SHA-256 `4fb312a904689dddc83c0e2221ee05ce49c0139c06ff15956ad993a9c48d582a`.

| `f5_n24_m24_d4` workload | Elimination ratio, 95% interval | Complete F5 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 0.457, 0.444–0.463 | 0.572, 0.558–0.579 |
| Holdout A | 0.442, 0.438–0.448 | 0.557, 0.554–0.562 |
| Holdout B | 0.443, 0.431–0.445 | 0.553, 0.542–0.555 |

The ratios are paired reference/candidate medians. Counted elimination word XORs were unchanged on the primary cell, so the new table replay and memory access cost outweighed its reduction in full-matrix sweeps. The default remains the original per-block kernel. A distinct compact, fused replay is frozen in the protocol for one bounded follow-up measurement; the initial regression is retained.

## Compact fused replay result

CI run [36306509051](https://github.com/aburan28/crypto/actions/runs/36306509051) at PR head `980c5337` completed 66 successful one-thread calls. All seven cases on the frozen seed and both holdouts matched exact row fingerprints, rank, pruning and criterion counts. The [complete raw receipt](runs/36306509051/ci-result.json) has SHA-256 `f46a5a180f42bdc9a681726fabce808b8adb964181e81af425e300436e28a32a`.

| `f5_n24_m24_d4` workload | Elimination ratio, 95% interval | Complete F5 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 0.404, 0.401–0.427 | 0.519, 0.506–0.544 |
| Holdout A | 0.383, 0.368–0.392 | 0.486, 0.481–0.501 |
| Holdout B | 0.386, 0.374–0.395 | 0.488, 0.476–0.496 |

All 18 smaller complete-call medians also fell below their own A/A lower bounds. This revision is slower than both the original kernel and the first column-panel candidate, so it fails the predeclared primary and no-smaller-regression gates. No four-thread control or further timing-driven revisions are warranted. The tested implementation is preserved in commits `726e7486` and `980c5337` and the raw receipts; the final PR removes its runtime and workflow changes, leaving this research record only. There is no measured F5 speedup from column panels and no one-target IC or DLP speedup claim.
