# Boolean F4 unsorted matrix products: negative paired result

The opt-in parity-hash path matched exact reduced bases, step counts, divisor tests, word XOR counts and matrix dimensions in all seven cases on the frozen seed and both holdouts. It missed the preregistered gain gate by a wide margin. Keep canonical sorted products as the production path; do not promote this candidate. Its full implementation remains reproducible at commit `f33cb0e4`; the final branch tree retains the protocol, harness and receipt without adding the slower runtime path.

CI run [36297079619](https://github.com/aburan28/crypto/actions/runs/36297079619) at PR head `f33cb0e4` passed both F4 test modes and completed 66 calls: one warmup per arm and five A/A plus five alternating A/B pairs per workload. The [complete raw receipt](runs/36297079619/ci-result.json) has SHA-256 `f56ad34366bd3dc7579b8e1e11bf817a8f17908d2595aba0552cc19bd47cc486`. The runner was Linux x86-64, AMD EPYC 7763, Rust 1.98.1, with one pinned CPU. The one-minute load average was 1.95 at start and 1.55 at end. Only ratios paired within this run are compared.

| `n20_m30` workload | Sorted build | Unsorted build | Paired build ratio, 95% interval | Sorted full F4 | Unsorted full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 290.7 ms | 397.6 ms | 0.734, 0.722–0.738 | 674.0 ms | 784.7 ms | 0.857, 0.852–0.864 |
| Holdout `1ac0ffee` | 292.4 ms | 395.9 ms | 0.738, 0.735–0.746 | 766.9 ms | 871.1 ms | 0.876, 0.875–0.887 |
| Holdout `2468ace0` | 288.0 ms | 394.7 ms | 0.728, 0.727–0.770 | 762.8 ms | 868.0 ms | 0.880, 0.873–0.922 |

Milliseconds are each arm's five A/B median; ratios are medians of paired reference/candidate samples. Frozen A/A ratios ranged 0.998–1.009 for build and 0.988–1.018 for the complete call. Elimination was unchanged (frozen paired median 1.003), so the lost time is in building matrix products. Every smaller case regressed in the complete F4 call beyond its A/A noise. This result shows that a fresh hash set per product costs more than the existing packed-key sort despite eliminating term order work. A different cancellation mechanism would need its own frozen experiment.

This is a solver-stage result. No one-target IC or DLP speedup, rho comparison, or candidate scoreboard value is inferred.
