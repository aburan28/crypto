# Cost assessment: the additional 30% goal

Outcome: **47.72% lower complete cold-solve instruction cost** against the previous verified fast-lifting batch1/window8 winner. Confirmation and fresh-process replay passed. Final campaign wall time: **23.5 minutes**.

## Measured research resources

| Work | Container CPU-hours | Peak RAM (GiB) | Wall time |
|---|---:|---:|---:|
| Candidate A development | 0.0631 | 1.43 | not separately metered |
| Final AutoLab campaign | 0.4031 | 1.69 | 23.5 min |

Total measured container use: **0.4662 CPU-hours**. The supervisor additionally recorded 64.82 host child-process CPU-seconds; the separate development audit used 2.98 host CPU-seconds. The final container resource snapshot includes its build, profiler and independent audit. No out-of-memory events were recorded.

The final campaign was capped at 2 CPUs, 8 GiB and 2 hours: at most 4 allocated CPU-hours and 16 GiB-hours. Allocation bounds are distinct from actual CPU use. The service completed and its container was removed.

## Dollar assessment

The benchmark harness used no model API calls and launched no cloud machines. Its incremental model API charge is $0. This does not price this assistant session.

The existing machine has no supplied CPU-hour or electricity rate. Total research dollar cost is therefore **unpriced**. The measured container component is `0.4662 × local CPU-hour rate`; add the separately recorded host CPU time and any priced preparation, storage, electricity and assistant usage. Host test compilation and image preparation were not fully metered, so this is not a complete invoice.

## Complete solver phase cost

The following are exclusive instruction sums across 180 matched confirmation processes per arm. They explain cost attribution; the headline uses the frozen, curve-balanced paired geometric mean. All values come from [frozen phase measurements](confirmation-phase-costs.json).

| Phase | Previous winner Ir | Candidate Ir | Matched rho Ir | Candidate / previous winner |
|---|---:|---:|---:|---:|
| collection and decomposition | 497,689,546 | 495,956,607 | 0 | 0.9965 |
| curve and targets | 1,238,475,504 | 1,220,589,063 | 1,238,475,504 | 0.9856 |
| factor base and tables | 11,087,651,858 | 8,334,604,484 | 0 | 0.7517 |
| final verification | 337,743,940 | 15,614,817 | 328,226,454 | 0.0462 |
| individual log | 1,059,430,105 | 736,941,581 | 0 | 0.6956 |
| log certification | 2,780,800,037 | 294,441,895 | 0 | 0.1059 |
| reporting and cleanup | 287,580,124 | 287,718,772 | 13,846,240 | 1.0005 |
| rho solve | 0 | 0 | 846,790,282 | not applicable |
| startup and input | 70,286,712 | 70,301,112 | 70,284,765 | 1.0002 |
| verify filter and linear algebra | 6,184,395,188 | 783,650,871 | 0 | 0.1267 |

Factor-base and table setup accounts for 68.1% of the candidate's summed confirmation cost. Verification/filtering/matrix work and log certification show the largest proportional reductions among the major previous costs. Startup, reporting and cleanup remain charged, including their small increases.

Ir includes the entire cold user-space worker: target/base construction, failed attempts, relation checks, matrix work, log certification, descent, final verification, serialization and cleanup. Kernel/device work, profiler execution, compilation and the external checker are outside solver Ir and belong to research overhead. The result is implementation engineering on the tested small Koblitz cases; native runtime and arithmetic complexity are not inferred.

[Machine-readable cost ledger](cost-assessment.json) · [Research result](README.md) · [Earlier campaign costs](../ic_autolab_harness_20260915/COSTS.md).
