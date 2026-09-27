# Boolean F4 bitmap membership: one-thread paired result

The opt-in bitmap passed the preregistered one-thread stage gate. Its exact output, counted work and matrix shapes matched the hash-set reference across all seven frozen cases and both new holdout workloads. This is an **engineering** improvement to Boolean F4's build phase; no IC candidate, one-target online DLP or rho speedup is measured or claimed.

CI run [36294892907](https://github.com/aburan28/crypto/actions/runs/36294892907) at PR head `3a8d3633` passed the reference and bitmap F4 tests, then completed 66 calls: one warmup per arm, five A/A reference pairs and five alternating A/B pairs for each workload. The complete receipt is [runs/36294892907/ci-result.json](runs/36294892907/ci-result.json), SHA-256 `65ab9560bb2402113fc225cd62472cb917920ecbe05cf724950ae379d2c68a92`. The pinned runner was Linux x86-64 on AMD EPYC 9V74, Rust 1.98.1, with four logical CPUs, 16.4 GB memory and AVX2/AVX-512/BMI2/POPCNT/PCLMULQDQ. Its one-minute load average was 1.98 at start and 1.30 at end; virtual-runner wall time remains hardware-specific.

| Workload, `n20_m30` | Reference build | Bitmap build | Paired build ratio, 95% bootstrap interval | Reference full F4 | Bitmap full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 248.2 ms | 150.0 ms | 1.650, 1.624–1.665 | 586.4 ms | 486.9 ms | 1.199, 1.191–1.214 |
| Holdout `1ac0ffee` | 250.7 ms | 153.0 ms | 1.639, 1.636–1.655 | 673.0 ms | 571.9 ms | 1.177, 1.172–1.224 |
| Holdout `2468ace0` | 250.7 ms | 152.8 ms | 1.633, 1.625–1.655 | 666.8 ms | 571.6 ms | 1.165, 1.149–1.170 |

The milliseconds are the median of each arm's five A/B calls; ratios are medians of paired reference/candidate ratios and need not equal the ratio of the marginal medians. The frozen A/A build ratio range was 0.995–1.009. Elimination was effectively unchanged (frozen paired median ratio 1.001). Smaller frozen cases had build ratios from 1.077 to 1.396 and full-F4 ratios from 1.040 to 1.083; no smaller case showed a median regression. All output fingerprints, basis lengths, steps, divisor-test counts, word-XOR counts and largest matrix dimensions matched in every call.

Promotion still requires a parallel control and a final default-on replay. Until then the bitmap remains opt-in. This measured gain is only a solver-stage diagnostic and does not change the repository's index-calculus scoreboard or one-target online accounting.
