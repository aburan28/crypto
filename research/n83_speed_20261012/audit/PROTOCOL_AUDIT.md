# Frozen campaign audit

The pre-launch source audit passes for the stated exploratory comparison.
This document reviews the protocol, not its eventual measured outcome.
Final completion, interval and replay status come from the native campaign
records and independent replay receipt.

| Requirement | Frozen implementation | Pre-launch status |
|---|---|---|
| Compare with already optimized rho | Current arm is Fold/reference scalar, 64 lanes, DP6, escape off, in the same executable as the candidate | Verified by source inspection |
| Tune fairly | 18 arithmetic runs and 81 policy runs on three separate training targets with three seeds each; select full nine-case process-wall mean | Verified by source inspection |
| Avoid best-seed/target selection | Complete nine-case candidate required; default/reference preferred within 2%; one policy frozen before holdout | Verified by source inspection and runner tests |
| Preserve failures | Every attempted row retained; unsuccessful attempt charged 60s; wrong answer rejects; no zero-cost imputation | Verified by source inspection and runner tests |
| Check arithmetic equivalence | All nine reference/Barrett pairs must match answers, exported outcomes, counters, preparation charges and policy | Verified by source inspection; not a full per-step trace |
| Check local timing noise | Five A/A pairs on the identical training target and seed, 10 launches, preceding the selection stages | Verified by source inspection |
| Fresh evaluation | Five new holdout targets, three fresh paired walk seeds each, 30 launches; no retuning from those results | Verified by source inspection; fixture/orbit checks have separate native receipt |
| One-target workload | Fresh child process and collision table for every run; no target batch or shared collision work | Verified by source inspection |
| Same resources | CPU3, one solver OS thread, 12 GiB address-space cap, 60s process cap; native isolation reservation and child receipts | Verified by source inspection; actual isolation receipt controls measured status |
| Complete cost | Primary statistic uses whole-process wall cost; online, CPU, RSS and walk counters are secondary records | Verified by source inspection |
| Honest interval | Per-target mean across three paired seeds, then five-target geometric ratio; 20,000 paired target-cluster bootstrap draws | Verified by source inspection and statistics tests |
| Complete optimization cost | Native analysis retains actual costs of all 139 children, including selection and A/A runs | Verified by source inspection |
| Independent correctness | Separate native C++ field/curve replay of every completed answer; sidecar not given to solver children | Required final gate; do not infer completion from producer booleans |

The full budget is 139 child launches: 10 A/A, 18 arithmetic, 81 policy and
30 holdout. The selected policy is the best eligible policy in this bounded
training grid after the stated 2% preference, not a global rho optimum.
The primary ratio compares current and selected configurations on new
targets in this campaign; it does not divide new measurements by old
timings or update the earlier IC/rho ratio using mismatched IC evidence.

The five-target cluster interval reflects the target/seed law and this VM.
It does not turn 15 seeded runs into 15 independent target clusters and
does not include selection-stage uncertainty as an unconditional search
estimator. A/A supplies an empirical local noise diagnostic, not a
certificate of physical host exclusivity or an ecbench L2 admission.
Nominal monitoring every 5ms adds scheduling and poll-body overhead;
there is no strict 5ms wall-clock bound under arbitrary scheduling.

No calibrated operation-equivalent `S`, independent-host timing audit,
full walk trace, asymptotic advantage or deployed-curve conclusion is
claimed. The final report must state all failed and invalid attempts,
actual isolation status, whole-process and online intervals, confidence
limits, arithmetic-mean comparison and the independent replay digest.
