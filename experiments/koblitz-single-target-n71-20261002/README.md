# N71 single-target IC/rho boundary run

The retained ledger series recovered one frozen public target on the Koblitz curve over \(\mathbb{F}_{2^{71}}\) in three paired runs. The median online wall ratio is **2,864.1×**: IC took 10.208 ms and rho 29,237.4 ms in the median pair. The same 14,554-probe IC relation was found in all three runs, and independent field and curve replay passes 15/15 checks for each. The [ledger record](LEDGER_RECORD.md) gives the full series and its frozen target.

| Ledger pair | IC online (ms) | rho online (ms) | rho / IC |
|---|---:|---:|---:|
| R1 | 14.457 | 70,747.3 | 4,893.5× |
| R2 | 10.208 | 29,237.4 | 2,864.1× |
| R3 | 11.874 | 18,200.9 | 1,532.8× |

The ratios are exploratory because the runs have no auditable host-level CPU isolation receipt. Repeats on this one point measure timing variation, not variation across target points. The earlier independent pilot below uses another point and remains a separate observation.

## Earlier pilot pair

The pilot recovered the public point `Q = ["1850589826078099102726","2157895285144092536687"]`. Its IC and rho online times were 1,699.089 ms and 28,513.091 ms, with producer checks and independent Sage replay of both scalars.

| Candidate | Workload | Run | Targets | IC online | rho online | rho / IC | Correctness |
|---|---|---|---:|---:|---:|---:|---|
| `IC1N71Ckb1fb85200PDP4rootRCguidedLAgaussTDdirectISO0h419f6af6f42e` | `cd9b5b764d88` | `...R1` | 1 | 1,699.089 ms | 28,513.091 ms | **16.7814×** | Both producer checks and independent Sage replay pass |

The pilot workload freezes that point in polynomial-basis decimal encoding. IC receives only the point. The known-answer scalar `314159265358979` is stored in a validation-only sidecar and is not part of IC input or target generation during either timed interval.

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

The independent Sage replay uses the checked launcher at `/Volumes/SSD990/cryptanalysis/sage`; it verifies the exact prime order 5,513,228,015,079,457, Frobenius action used by the signed-orbit labels, and that both pilot scalars map the generator to the frozen point. The pilot paired wall ratio is 16.7814×. One pilot pair provides no statistical uncertainty interval. IC probes and rho walk steps are retained in the ledger series' `claim_report_vs_rho.json` as separate native units, leaving `S = total_operations / sqrt(r)` unknown until calibration fixes a common boundary.

Per-process peak RSS is unknown: the IC producer returned null and macOS `time -l` emitted a sysctl permission warning before returning its wall/user/system times; both producer rows themselves report successful solves. The run used one process and one worker per arm under the declared 16 GiB envelope, but actual peak-RSS compliance could not be verified. The host was arm64 macOS 26.6; the CPU model query was denied.

The immutable candidate identity is in [candidate-manifest.json](candidate-manifest.json). The frozen workload and paired measurement are in [runs/R1/single-target-result.json](runs/R1/single-target-result.json) and [single-target-results.csv](single-target-results.csv). Raw producer rows, timing output, runtime receipt, and independent replay artifacts are retained under `runs/R1/`.

The retained ledger `claim_report_vs_rho.json` has status **FAIL** under the current autolab `vs_rho` schema. Canonical identity and replay-digest fields, exact five-phase IC timing, the complete resource envelope, and rho policy are not present in that report; the verified raw ledger and pilot results remain available at the paths above and in [LEDGER_RECORD.md](LEDGER_RECORD.md).
