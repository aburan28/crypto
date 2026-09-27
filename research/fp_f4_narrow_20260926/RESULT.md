# Prime-field F4 narrow arithmetic: first paired run

The 32-bit row and branch-free update path is correct on the frozen and two holdout suites, but **not promoted**: its primary complete-F4 gain missed the preregistered 1.2× median gate. This is an internal solver-stage engineering diagnostic, not an IC or end-to-end DLP speedup.

CI run [36293775047](https://github.com/aburan28/crypto/actions/runs/36293775047) at PR head `814b9f958` passed the `f4_fp` tests, built the release benchmark and completed all 66 warmup/A/A/A/B calls. The complete receipt is [runs/36293775047/ci-result.json](runs/36293775047/ci-result.json), SHA-256 `6d70c7f571333520b7cc8bc623b804d94f715d74dd96bef1e7027892c4f66b68`. The one-CPU pinned Ubuntu x86-64 runner reported AMD EPYC 7763, four logical CPUs, 16.4 GB memory, AVX2/BMI2/POPCNT/PCLMULQDQ and Rust 1.98.1. Its one-minute load average fell from 1.87 to 1.02 over the run.

| Workload, `quad_n8_p65521` | Reference F4 | Narrow F4 | Paired reference/candidate median, 95% bootstrap interval |
| --- | ---: | ---: | ---: |
| Frozen seed | 2229.5 ms | 2004.1 ms | 1.113, 1.111–1.121 |
| Holdout `badc0de1` | 2229.4 ms | 2003.2 ms | 1.113, 1.111–1.114 |
| Holdout `5eed2026` | 2215.2 ms | 2000.0 ms | 1.108, 1.104–1.113 |

Milliseconds are medians of each arm's five A/B calls, while ratios are medians of paired ratios. The frozen A/A range was 0.995–1.006. All 13 cases on all three workloads matched the reference basis and solution fingerprints, step counts, matrix shapes and basis sizes. The larger `quad_n9_p65521` cell improved by about 1.13× on all three workloads, while tiny cells were mostly inside the noise and `overdet_n6_m9_p101` regressed slightly on the frozen and first holdout. The complete raw receipt retains every cell, including this regression.

The first candidate eliminates a branch and halves dense row storage, but still reduces each multiply-add. The next bounded iteration will defer reductions across each row's existing pivots, with `F4_FP_NARROW=2` as a distinct arm. It must pass the same complete-F4 and exact-output gates before changing the default. The rho reference, one-target IC cost and scoreboard remain unchanged and unknown for this solver-stage experiment.

## Deferred mode-2 one-thread result

CI run [36296331241](https://github.com/aburan28/crypto/actions/runs/36296331241) at PR head `b6861df0` passed the expanded F4 exact-output tests and completed all 66 paired calls. The [complete raw receipt](runs/36296331241/ci-result.json), SHA-256 `dd93dd72f8361cb1c9835255155e1cf09275fab94e5e0af60b8821e3a209526e`, records the original mode `0` against deferred mode `2`. Its one-CPU pinned Linux x86-64 runner reported AMD EPYC 7763 and Rust 1.98.1.

| `quad_n8_p65521` workload | Reference F4 | Deferred F4 | Paired reference/candidate median, 95% bootstrap interval |
| --- | ---: | ---: | ---: |
| Frozen seed | 2264.9 ms | 1384.8 ms | 1.641, 1.631–1.649 |
| Holdout `badc0de1` | 2258.0 ms | 1374.5 ms | 1.641, 1.616–1.654 |
| Holdout `5eed2026` | 2254.3 ms | 1375.6 ms | 1.639, 1.621–1.647 |

The frozen A/A range was 0.996–1.007, and the complete-F4 1.2× gate is clear. Every basis and solve fingerprint, step count, matrix shape, basis size and verdict matched mode `0` across all 13 cases and both holdouts. The larger `overdet_n12_m24_p65521` case improved by a paired 2.101× (2.086–2.155) on the frozen workload. No smaller case had a median regression outside its A/A range. The milliseconds are each arm's five A/B median; paired ratio medians need not equal the ratio of those marginal medians.

Default promotion awaits the preregistered four-thread regression control and a final default-on replay. These are point-decomposition solver-stage timings, not one-target IC online timings or verified DLP speedups.

## Four-thread control and default decision

CI run [36298267078](https://github.com/aburan28/crypto/actions/runs/36298267078) at PR head `9341c7da` passed F4 tests and completed 66 four-thread paired calls. The [complete raw receipt](runs/36298267078/ci-result-threads4.json), SHA-256 `0f430e17634598c654bcfc9433faef1cab619a70c4c341fb610e249a46b46e82`, records a Linux x86-64 AMD EPYC 9V74 runner with four pinned CPUs. It used a different runner from the one-thread result; all ratios below are paired within this run.

| `quad_n8_p65521` workload | Four-thread A/A range | Paired mode-0/mode-2 complete-F4 ratio, 95% interval |
| --- | ---: | ---: |
| Frozen seed | 0.997–1.002 | 1.580, 1.577–1.594 |
| Holdout `badc0de1` | 0.994–1.001 | 1.588, 1.573–1.591 |
| Holdout `5eed2026` | 0.996–1.003 | 1.586, 1.575–1.591 |

Every case matched the reference basis and solve fingerprints, steps and matrix shapes. The tiny `quad_n4_p31` frozen cell's median ratio of 1.002 sits 0.002 below its A/A minimum of 1.004 but still favors mode `2`; no material smaller-case regression appeared. The primary result easily clears the four-thread regression control. Mode `2` is selected by default for supported primes, with `F4_FP_NARROW=0` retaining the original wide path and `=1` the earlier intermediate. A final one- and four-thread replay against explicit mode `0` is required before merging. These remain F4 stage timings, not online IC or verified DLP speedups.
