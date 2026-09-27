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
