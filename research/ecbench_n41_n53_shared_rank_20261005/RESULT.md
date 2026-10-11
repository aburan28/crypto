# n41/n53 shared-rank K8/K16 counted panel: the n37 parameters do not transfer

The [preregistered protocol](PROTOCOL.md) carried the n37 shared-rank
index-calculus arms unchanged (K8 and K16 `compact-orbit-scan` bases, same
rank and target policies, identical K16 control) to the 39-bit and 44-bit
Koblitz subgroups, beside same-point strong signed-Frobenius rho, on 16
fresh public one-target workloads per curve, five measured rounds each.
Every run is L0 on macOS; counts are the metric, wall is exploratory.

**Decision: `frozen_parameters_do_not_transfer`** (rule 3 of the
protocol). Class: **accounting** (a new size priced with the n37 method
unchanged; no admitted speedup).

| curve | log₂ r | arm | measured | verified | failures | mean cold S, lower bound |
|:--|--:|:--|--:|--:|:--|--:|
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | strong rho | 80 | 80 | — | 0.157 |
| | | IC K8 | 80 | 0 | 80 rank setup did not verify every base column | unknown |
| | | IC K16 (1,312 points, 16 columns) | 80 | 23 | 57 exhausted | 42.57 |
| | | IC K16 control | 80 | 23 | 57 exhausted | 42.57 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | strong rho | 80 | 80 | — | 0.141 |
| | | IC K8 | 80 | 0 | 80 rank setup did not verify every base column | unknown |
| | | IC K16 | 80 | 0 | 80 rank setup did not verify every base column | unknown |
| | | IC K16 control | 80 | 0 | 80 rank setup did not verify every base column | unknown |

Secondary counted ratio of sums, verified pairs only, 20,000-resample
workload-block interval at seed `202610055103`:

| curve | comparison | ratio of sums | 95% interval | pairs |
|:--|:--|--:|:--|--:|
| n41 | K16 / rho | **256.2** | [202.5, 332.8] | 23 |
| n41 | K16 / control | 1.000 | [1.000, 1.000] | 23 |
| n41 | K8 / rho, K16 / K8 | not evaluable | — | 0 |
| n53 | every comparison | not evaluable | — | 0 |

- **H1** (counted cold IC/rho above 1): *retained* for K16 at n41
  (interval wholly above 1); *not evaluable* for K8 at both sizes and for
  K16 at n53, because those arms never verified a solve.
- **H2** (K16 target-only work below K8's): *not evaluable* at both sizes;
  K8 has no verified run.
- **A/A**: the counted K16/control ratio is exactly 1.0 on all 23 verified
  pairs, as determinism requires. Its L0 wall deviation reaches 76.8%,
  which says only that the host was loaded.
- **Primary target-zero rows**: incomplete at both sizes (no paired verified
  IC round on the designated workload), so no primary quotient exists.

**What it means.** At n37 the same arms solved every target, at a counted
IC/rho lower bound of 4.55× for K8 and 5.03× for K16. At n41 the K8 base no longer
reaches full rank in 100,000 target-blind trials, and K16 reaches it but
then exhausts its 512 target attempts on 57 of 80 runs; when K16 does solve,
it costs 256× rho counted. At n53 neither base reaches full rank. The frozen
n37 parameters are tuned to n37 and do not scale, so an L2 wall-time run of
this configuration at n41/n53 would measure a method that mostly fails. The
next step on this line is a base whose size scales with the subgroup, not
an isolated host for these parameters.

**Integrity.** `ecbench verify --replay-all` reproduced all 206 deterministic
verified runs, 0 problems (`AUDIT.json`, receipt SHA-256
`8cb94fac9194ac3bf622b5a36401bb6896d3f1feea365245744ca023164e4a83`). Session
`ECBS1he98c2d9648b6`, 768 records (128 warmup, 640 measured), every failure
and exhaustion kept.

**One analyzer correction, made after the run and stated here.** The
analyzer's integrity guard rejected any unverified run that carried a cost.
The harness charges the work an `exhausted` run spent, which the protocol
requires (failed attempts are charged and kept, never counted as solves),
so the guard was wrong for this outcome, which the pilot never produced.
The guard now accepts a cost on `exhausted`, `error` and `timeout` runs
and still rejects one on any other unverified status. No decision rule,
statistic, seed or cap changed, and no such run enters any quotient.
`examples/ecbench_n41_n53_shared_rank_analyze.rs` reproduces `DECISION.json`
from the frozen files, the session and the receipt.
