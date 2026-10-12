# Goal28: measured GPU profiling and tuning

**28 B/s has not been reached.** G7e requests for the target RTX PRO 6000
returned insufficient capacity. These diagnostics and measurements use one
`g7.4xlarge` RTX PRO 4500 Blackwell Server Edition, driver 595.91.07, native
sm_120 code compiled by CUDA 13.3.73. They cannot establish a speedup over the
historical RTX PRO 6000 rates of 14.637530 B/s benchmark and 14.106673 B/s DP34.

This is **engineering** of the Pollard-walk kernel. The generic-group boundary
remains Ω(√(n/262)); field products remain 5.3125 per scalar update, ratio 1.0
to control. No full-DLP speedup or exponent change was measured. The frozen
index-calculus/WDSat solver suite is outside this kernel experiment's scope.

## Equal-work grid screen

The [contract](benchmarks/goal28/contract.json) preceded allocation. Each of
18 samples completed **201,863,462,912 scalar updates**, with zero dropped
reports. There are three repetitions per grid and workload. Each row below
uses median B completed updates/s; ratios compare with the original grid.

| Workers | Field products/update | Product ratio | Benchmark B/s | Rate ratio | DP34 B/s | Rate ratio | Correctness / class |
|---:|---:|---:|---:|---:|---:|---:|---|
| 385,024 control | 5.3125 | 1.0 | 4.649429 | 1.000000 | 4.602291 | 1.000000 | passed / engineering |
| 192,512 | 5.3125 | 1.0 | 4.534123 | 0.975200 | 4.483364 | 0.974159 | passed / engineering |
| 96,256 | 5.3125 | 1.0 | 4.537225 | 0.975867 | 4.473096 | 0.971928 | passed / engineering |

Both smaller grids were slower. Grids use different seed panels, with more
launches on smaller grids to preserve completed work. The middle repetition
reverses grid order. First rows are retained without asserting thermal steady
state. This is a screen, not a matched-corpus candidate promotion.

All GPU arithmetic, storage, shared-sigma and checkpoint/restart integration
checks passed. An offline archive audit revalidated all 18 raw rows and the
compiled binary hash. It also verified identical sorted DP multisets across
the three repetitions **within each grid**; multisets differ across grids.

## Nsight findings

The initial Nsight Systems 2026.3.2 jobs produced two usable `.nsys-rep` captures. In the
screening job's short trace, four walk kernels totaled **628.592 ms** and
initialization took 310.023 ms. All five `cudaLaunchKernel` API calls totaled
**0.113 ms**. This supports investigating work inside the kernel before
launch overhead. The six `cudaDeviceSynchronize` calls wait for GPU work;
their 938.636 ms overlaps GPU execution and must not be added as CPU overhead.
Initialization occupies a large share of this deliberately short trace.
**No profiler duration is used as a throughput measurement.**

Nsight Compute 2026.3.0.0, build 38525999, failed before producing usable
hardware counters. The full section set with kernel replay failed on both
385,024 and 96,256 workers. A separate VM retried six hardware sections
(`SpeedOfLight`, `Occupancy`, `SchedulerStats`, `WarpStateStats`,
`MemoryWorkloadAnalysis`, `ComputeWorkloadAnalysis`) with application replay;
both grids failed again. All four commands returned code 9 with
`Failed to prepare kernel for profiling` and `Unknown error on device 0.`
The cause remains unresolved. No memory, scheduler-stall or pipeline-counter
bottleneck is claimed.

The control's static resource receipt reports 104 registers/thread and zero
stack/local bytes. The runtime distinguishes 1,792 bytes of kernel shared
storage and 1,024 bytes reserved by the driver; the static total is 2,816.
These observations motivate the separate
[launch-bound experiment](benchmarks/goal28/occupancy-contract.json), not a
claim that reducing registers will improve performance.

Two further diagnostic jobs tested a client with `ECC_PROFILE_RANGE=1` and
an independent CUDA smoke workload: 4,096 unsigned LCG results, each with
1,024 updates, checked against the CPU. Both instrumented-client integration
and ordinary smoke execution passed. Nsight Systems produced two additional
timelines. These instrumented durations are excluded from throughput results.

The smoke workload reproduced the kernel-preparation error with code 9. After
correcting an unsupported `--profile-from-start` option in the first range
attempt, smoke **application range replay also failed**, with code 9 and
`Failed to prepare range for profiling` / `Unknown error on device 0.`
The failure therefore is not confined to ECC arithmetic. Its underlying
profiler/platform cause remains unresolved; a specific driver or permission
cause has not been established.

The corrected ECC range commands returned code 0 but warned **“No ranges were
profiled.”** No `.ncu-rep` or numeric counters were produced. The runner's
report-existence check correctly kept the overall result invalid. The retained
second-range selection did not produce a usable capture; this is not evidence
that ECC range profiling succeeded or reached the same preparation failure.
[range-summary.json](benchmarks/goal28/range-summary.json) and
[range-evidence.tgz](benchmarks/goal28/range-evidence.tgz) preserve both attempts,
their exact sources, full logs and timelines. Range-capture selection and the
independent smoke profiling failure remain open diagnostic issues.

## Paired launch-bound experiment

The [separately frozen contract](benchmarks/goal28/occupancy-contract.json)
changed only `ECC_MINBLOCKS` from 2 to 3. Both clients passed GPU probes,
reference-report/checkpoint integration and byte-identical cross-binary
checkpoints at 8, 128 and 257 workers. CUDA's occupancy API confirmed an
increase from two to three resident blocks per SM. Registers fell from 104
to 80, with zero stack/local bytes in both builds. Static walk instructions
increased from 4,513 to 4,569; this count is not a dynamic throughput measure.

Four warmups were excluded. Five alternating pairs per workload each completed
201,863,462,912 scalar updates per sample. Every one of the ten pairs favored
the control. All DP34 warmups and timed samples produced the same multiset of
5,149 records, with zero drops.

| Launch bound | Field products/update | Product ratio | Benchmark B/s | Rate ratio | DP34 B/s | Rate ratio | Correctness / class |
|---|---:|---:|---:|---:|---:|---:|---|
| MINBLOCKS=2 control | 5.3125 | 1.0 | 4.692014 | 1.000000 | 4.632696 | 1.000000 | passed / engineering |
| MINBLOCKS=3 | 5.3125 | 1.0 | 4.593011 | 0.978900 | 4.542281 | 0.980483 | passed / engineering, rejected |

The candidate is **2.11% slower in the benchmark and 1.95% slower in DP34**.
Its paired geometric speed ratios are 0.982453 and 0.980758. The two-sided
95% Student-t intervals for paired log ratios, transformed back to speed
ratios, are [0.974553, 0.990417] and [0.973004, 0.988573], respectively.
These exploratory intervals assume independent, approximately normal log
ratios; neither supports improvement. More resident warps did not make this
binary faster. No runtime stall cause is inferred from this experiment.

[occupancy-summary.json](benchmarks/goal28/occupancy-summary.json) preserves
all pair rates/counts, resource receipts and corpus hashes.
[occupancy-evidence.tgz](benchmarks/goal28/occupancy-evidence.tgz) preserves
the exact tested source, full raw runner JSON and launch receipt. Both archived
binaries were checked against the hashes recorded before and after measurement.

## Retained profiling evidence and failures

[profiling-summary.json](benchmarks/goal28/profiling-summary.json) contains
individual rows, ratios, hashes, hardware identities and profiler status.
[profiling-evidence.tgz](benchmarks/goal28/profiling-evidence.tgz) contains
unaltered raw JSON from all three attempts, launch receipts, both timelines,
profiler logs and exact source archives for the two profiling attempts.
Complete binaries, corpora and bootstrap logs remain in the approved private
benchmark bucket at the receipt's `s3ResultPrefix` under `results.tgz`; archive
SHA-256 values are recorded. No credentials or signed URLs are included.

The first attempt completed correctness and one timing sample, then exposed
a parser bug with padded columns. That sample is not part of the screen.
The second completed the entire screen and timeline, then failed a postrun
tuple/dictionary bookkeeping check. Its overall `valid:false` is preserved;
the archive checks described above independently establish the retained
screen evidence. Both harness bugs were fixed and the actual timing log was
used to verify the parser correction. The third attempt completed its checks
and timeline but retained `valid:false` because counter collection failed.

No arithmetic flag, launch-bound default or production preset is promoted.

Allocation cleanup is recorded separately in
[cleanup.json](benchmarks/goal28/cleanup.json). Every allocated ID received an
exact-ID termination request. The receipt distinguishes currently returned
terminated instances from an older ID no longer returned by EC2 and checks
for any pending, running, stopped, stopping or shutting-down IDs.
The final check found none: all five currently returned IDs were terminated;
the oldest allocation was no longer returned by EC2.
