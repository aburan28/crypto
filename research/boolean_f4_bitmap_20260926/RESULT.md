# Boolean F4 bitmap membership: one-thread paired result

The opt-in bitmap passed the preregistered one-thread stage gate. Its exact output, counted work and matrix shapes matched the hash-set reference across all seven frozen cases and both new holdout workloads. This is an **engineering** improvement to Boolean F4's build phase; no IC candidate, one-target online DLP or rho speedup is measured or claimed.

CI run [36294892907](https://github.com/aburan28/crypto/actions/runs/36294892907) at PR head `3a8d3633` passed the reference and bitmap F4 tests, then completed 66 calls: one warmup per arm, five A/A reference pairs and five alternating A/B pairs for each workload. The complete receipt is [runs/36294892907/ci-result.json](runs/36294892907/ci-result.json), SHA-256 `65ab9560bb2402113fc225cd62472cb917920ecbe05cf724950ae379d2c68a92`. The pinned runner was Linux x86-64 on AMD EPYC 9V74, Rust 1.98.1, with four logical CPUs, 16.4 GB memory and AVX2/AVX-512/BMI2/POPCNT/PCLMULQDQ. Its one-minute load average was 1.98 at start and 1.30 at end; virtual-runner wall time remains hardware-specific.

| Workload, `n20_m30` | Reference build | Bitmap build | Paired build ratio, 95% bootstrap interval | Reference full F4 | Bitmap full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 248.2 ms | 150.0 ms | 1.650, 1.624–1.665 | 586.4 ms | 486.9 ms | 1.199, 1.191–1.214 |
| Holdout `1ac0ffee` | 250.7 ms | 153.0 ms | 1.639, 1.636–1.655 | 673.0 ms | 571.9 ms | 1.177, 1.172–1.224 |
| Holdout `2468ace0` | 250.7 ms | 152.8 ms | 1.633, 1.625–1.655 | 666.8 ms | 571.6 ms | 1.165, 1.149–1.170 |

The milliseconds are the median of each arm's five A/B calls; ratios are medians of paired reference/candidate ratios and need not equal the ratio of the marginal medians. The frozen A/A build ratio range was 0.995–1.009. Elimination was effectively unchanged (frozen paired median ratio 1.001). Smaller frozen cases had build ratios from 1.077 to 1.396 and full-F4 ratios from 1.040 to 1.083; no smaller case showed a median regression. All output fingerprints, basis lengths, steps, divisor-test counts, word-XOR counts and largest matrix dimensions matched in every call.

## Four-thread control and default decision

CI run [36296268679](https://github.com/aburan28/crypto/actions/runs/36296268679) at PR head `25122240` repeated the one-thread comparison and ran the preregistered four-thread control: 66 calls per thread count, exact fingerprints and counters throughout. The complete receipts are [one thread](runs/36296268679/ci-result.json), SHA-256 `980bbdd30501d82bfeb26ac84e83a3a2a8486492f9481697debd60c85a71d587`, and [four threads](runs/36296268679/ci-result-threads4.json), SHA-256 `e150b0d23921ca67f0b322f6f3cab9b5964001808925fa9c7a9870edaeee0464`. This run used a different Linux x86-64 runner, AMD EPYC 9V45, Rust 1.98.1; ratios are paired within each run and thread count.

| Thread count, `n20_m30` | Frozen build ratio, 95% interval | Frozen full-F4 ratio, interval | Holdout A full ratio | Holdout B full ratio |
| --- | ---: | ---: | ---: | ---: |
| One | 1.694, 1.655–1.744 | 1.211, 1.168–1.265 | 1.208 | 1.206 |
| Four | 1.010, 0.984–1.121 | 0.979, 0.911–1.086 | 1.033 | 1.004 |

The four-thread frozen full-F4 ratio of 0.979 sits inside its A/A range of 0.973–1.110; the two holdout ratios sit inside their own A/A ranges. Every smaller four-thread cell's median was within its A/A range or improved. This clears the preregistered parallel regression control, while giving no evidence of a four-thread speed gain. The bitmap is therefore enabled by default for at most 22 variables; `F4_F2_BITMAP_SEEN=0` retains the hash-set reference. A final default-on replay of both thread counts is required before merging. The measured gains remain solver-stage diagnostics and do not change the index-calculus scoreboard or one-target online accounting.

## Final default-on confirmation

CI run [36297004647](https://github.com/aburan28/crypto/actions/runs/36297004647) at PR head `3cfaec73` passed the explicit hash-reference and default-bitmap tests, then completed 66 calls at each thread count. The receipts are [one thread](runs/36297004647/ci-result.json), SHA-256 `f42931269479e291cc2ab946566a160bb4ecfe944fcf19b06196984971694ab9`, and [four threads](runs/36297004647/ci-result-threads4.json), SHA-256 `c6cc5123377f37fa3638f1788faf7be0a1e9312f8a1da48a6df2e90401688f95`. The runner was Linux x86-64 on AMD EPYC 9V74, Rust 1.98.1; all compared runs were pinned to the declared CPU count.

| Thread count, `n20_m30` | Frozen build ratio, 95% interval | Frozen full-F4 ratio, interval | Holdout A full ratio | Holdout B full ratio |
| --- | ---: | ---: | ---: | ---: |
| One | 1.659, 1.644–1.664 | 1.199, 1.194–1.203 | 1.169 | 1.173 |
| Four | 1.032, 1.000–1.117 | 1.019, 0.996–1.045 | 1.020 | 1.026 |

All seven cases on both thread counts and all three workloads retained identical fingerprints, counted operations and matrix shapes. No full-F4 cell had a median regression beyond its own A/A noise. The one-thread build and complete-call gain passed again; the four-thread result is compatible with no material change. This completes the default-on promotion gate for at most 22 variables. The measured numbers remain Boolean F4 stage results, not a one-target IC or DLP speedup.
