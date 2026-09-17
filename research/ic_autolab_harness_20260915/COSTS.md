# AutoLab cost assessment

Final outcome: **26.9% lower complete solve instruction cost**, independently confirmed and replayed against the original batch16 baseline. The final campaign completed in 25.9 minutes.

## Measured resources

| Work | Container CPU-hours | Peak RAM (GiB) | Wall time |
|---|---:|---:|---:|
| autolab-20260915T193057Z | 1.1191 | 2.27 | 73.6 min |
| autolab-20260915T205319Z | 0.4445 | 1.70 | 25.9 min |
| Source-development probes | 0.1307 | 1.49 | not separately metered |

Total measured container use: **1.694 CPU-hours**. Host orchestration/report child processes additionally recorded 130.90 CPU-seconds; the separate search-history audit used 21.76 host CPU-seconds.

Each full campaign was capped at 2 CPUs, 8 GiB RAM and 2 hours wall time: at most 4 allocated CPU-hours and 16 GiB-hours per campaign. Both finished within those bounds; no out-of-memory events were recorded. These allocation bounds differ from actual CPU time.

## Dollar assessment

The benchmark coding agent reported **$0 inference cost** across 102 completed cost-bearing steps. Its provider catalogs listed zero prices; unfinished requests were not independently invoiced. The source-development probes and final frozen-candidate verifier made no additional model API calls. This does not price this assistant session.

No cloud machines were purchased. The existing local machine has no supplied CPU-hour, storage or electricity rate, so total research dollar cost remains **unpriced**. The measured container component is `1.6943 × local CPU-hour rate`; add the separately recorded host CPU usage and any priced preparation, storage or electricity. Initial installation/image setup, host compilation and other unmetered overhead prevent treating this as a complete invoice.

## Algorithm cost boundary

The promoted ratio is measured in Valgrind 3.22 amd64 Ir and covers the whole cold worker process, including target/base construction, relation collection, required checks, matrix work, final scalar recovery, serialization and cleanup. Compiler work, the profiler implementation, external audits and kernel/device work are research costs outside solver Ir. No arithmetic-complexity or native-runtime speedup is inferred.

[Machine-readable cost ledger](cost-assessment.json) · [Source and reproduction](source-improvement/README.md) · [Configuration-search history](search-history/README.md).
