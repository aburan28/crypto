# N71 single-target IC/rho boundary run

One paired run recovered the same public target on the Koblitz curve over \(\mathbb{F}_{2^{71}}\): compact-orbit index calculus took 1,699.089 ms online and Pollard rho took 28,513.091 ms online. Producer checks and independent Sage replay verify both scalars.

| Candidate | Workload | Run | Targets | IC online | rho online | rho / IC | Correctness |
|---|---|---|---:|---:|---:|---:|---|
| `IC1N71Ckb1fb85200PDP4rootRCguidedLAgaussTDdirectISO0h419f6af6f42e` | `cd9b5b764d88` | `...R1` | 1 | 1,699.089 ms | 28,513.091 ms | **16.7814×** | Both producer checks and independent Sage replay pass |

The workload freezes the single point
`Q = ["1850589826078099102726","2157895285144092536687"]` in polynomial-basis decimal encoding. IC receives only that point. The known-answer scalar `314159265358979` is stored in a validation-only sidecar and is not part of IC input or target generation during either timed interval.

The IC online clock starts immediately before target-dependent query hashing, after factor-base construction, the root index, and all factor-base logarithms are ready. It includes direct target decomposition and relation verification, scalar descent, and the in-process `[d]G=Q` check. Its exclusive online phase accounting is:

| IC phase | Time |
|---|---:|
| Target query orchestration, including the clock reconciliation residual | 0.000458 ms |
| Target PDP plus exact relation/group check (combined timer in this producer) | 1,698.943709 ms |
| Target descent | 0.026250 ms |
| Recovery check | 0.118875 ms |
| **IC online wall time** | **1,699.089292 ms** |

Rho receives the identical point. Its online clock includes target-specific jump setup (9.016 ms), the walk (28,496.694 ms), and in-process scalar/group validation (7.381 ms). Point construction is excluded for both arms.

IC target-independent preparation took 1,527,245.078 ms in the same process: factor-base load/build 240.889 ms, normal-basis setup/self-test 8.205 ms, index build 41,313.107 ms, guided rank acquisition 1,485,679.143 ms, and final Gaussian elimination 3.734 ms. The whole IC process took 1,528,946.921 ms. These costs are retained as setup data and excluded from the primary one-target online ratio.

The independent Sage replay uses the checked launcher at `/Volumes/SSD990/cryptanalysis/sage`; it verifies the exact prime order 5,513,228,015,079,457, Frobenius action used by the signed-orbit labels, and that both recovered scalars map the generator to the frozen point. The paired wall ratio is 16.7814×. Its controlled speedup remains unknown because a host-level isolation receipt is absent; one pair also provides no statistical uncertainty interval. IC probes and rho walk steps are retained in `claim_report_vs_rho.json` as separate native units, leaving `S = total_operations / sqrt(r)` unknown until calibration fixes a common boundary.

Per-process peak RSS is unknown: the IC producer returned null and macOS `time -l` emitted a sysctl permission warning before returning its wall/user/system times; both producer rows themselves report successful solves. The run used one process and one worker per arm under the declared 16 GiB envelope, but actual peak-RSS compliance could not be verified. The host was arm64 macOS 26.6; the CPU model query was denied.

The immutable candidate identity is in [candidate-manifest.json](candidate-manifest.json). The frozen workload and paired measurement are in [runs/R1/single-target-result.json](runs/R1/single-target-result.json) and [single-target-results.csv](single-target-results.csv). Raw producer rows, timing output, runtime receipt, and independent replay artifacts are retained under `runs/R1/`.

The retained `claim_report_vs_rho.json` has status **FAIL** under the current autolab `vs_rho` schema. Canonical identity and replay-digest fields, exact five-phase IC timing, the complete resource envelope, and rho policy are not present in that report; the verified raw result remains available at the paths above.
