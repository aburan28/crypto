# Two-block GF(2) word batching: negative stage result

The candidate matched the original matrix-F5 row fingerprints, rank, pruning and criterion work on every frozen and holdout case, and the shared GF(2) kernel unit suite matched exact RREF. It did not improve the primary elimination stage and slowed several smaller F5 calls. The production kernel remains unchanged; the experimental source is preserved in PR #883 commit `b2500466` for reproduction.

CI run [36303490393](https://github.com/aburan28/crypto/actions/runs/36303490393) completed 66 calls on a pinned one-thread Linux runner. The [full raw receipt](runs/36303490393/ci-result.json) has SHA-256 `bd2373083921c45a5f03ee642bb076f16c1860c3d22fa354a40a703de6a95096`. Ratios are per-pair reference/candidate medians, with exact five-pair bootstrap 95% intervals.

| `f5_n24_m24_d4` workload | Elimination ratio, 95% interval | Complete F5 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen | 0.997, 0.969–1.011 | 1.005, 0.979–1.017 |
| Holdout A | 1.000, 0.989–1.004 | 1.008, 0.999–1.016 |
| Holdout B | 1.011, 1.007–1.022 | 1.012, 1.008–1.016 |

The primary elimination result misses the preregistered 1.20 gate by a wide margin. Smaller `n12_m12_d4` and `n16_m16_d4` complete-call medians fell below their own A/A ranges on all three workloads. The two blocks from one 64-column word did not eliminate enough traffic to pay for extra tables and selected-row updates. Any broader cache-sized column-panel method is a separate experiment. These are internal F5 stage diagnostics, not one-target IC or DLP speedups.
