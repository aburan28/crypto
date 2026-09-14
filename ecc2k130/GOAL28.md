# Goal: 28 billion scalar updates per second

The target is one RTX PRO 6000 Blackwell Server Edition at 28 B completed
scalar updates/s, with correct DP collection and no regression in the DP34
workload. The historical measured references are 14.637530 B/s for the
benchmark and 14.106673 B/s for DP34. These are hardware measurements, not
algorithmic lower bounds. A new same-device paired control is required.

The [frozen contract](benchmarks/goal28/contract.json) names the goal, equal-work
screen and promotion requirements before measurement. This is Pollard-walk
engineering. The generic-group square-root boundary and the 5.3125 field
products/update accounting remain unchanged. Full-DLP speedup is unmeasured.
The WDSat index-calculus solver regression is inapplicable to this CUDA kernel.

The [measured profiling and grid-screen results](GOAL28-RESULTS.md) are retained.
The previous native-square and Karatsuba candidates failed acceptance on G7.
Fewer static instructions did not establish a speedup, so the next round
collects counters before selecting another arithmetic change.

## Profiling and tuning screen

`codegen/profile_goal28.py` builds the frozen control and its GPU probes. It
checks device arithmetic/storage/shared masks and checkpoint/restart integration,
then uses both profilers:

- Nsight Systems (`nsys`) records CUDA API calls, kernels and transfers. Its
  timeline separates startup, launch/synchronization overhead and GPU execution.
- Nsight Compute (`ncu`) attempts the packed walk at the original grid and a
  quarter-size grid, including the full section set. A counter-only recovery
  mode selects six hardware sections with application replay. Both collection
  modes failed at kernel preparation on the measured G7; no usable counters
  were obtained. Replay uses explicit cache flushing and no clock changes.
  Profiler durations never enter throughput comparisons.

An unprofiled screen measures 385,024 / 192,512 / 96,256 workers with 32 / 64 /
128 launches respectively. Each sample still completes 201,863,462,912 scalar
updates. Three repetitions cover benchmark and DP34, with reversed grid order
on the middle repetition. Grid sizes use distinct seed panels, and the first
rows are retained without assuming thermal steady state. The screen does not
claim a matched-corpus gain or promote a preset. A promising setting needs a
fresh paired confirmation with at least five pairs, a 95% paired confidence
interval excluding no improvement and independent correctness checks.

The smaller-state-footprint hypothesis failed this throughput screen: both
smaller grids were slower. The timeline shows tiny kernel-launch API time
relative to walk execution. Counter collection failed, so neither a memory
limit nor an issue/stall limit is established.

The next [frozen experiment](benchmarks/goal28/occupancy-contract.json) compares
`ECC_MINBLOCKS=2` with `3` while leaving all arithmetic flags and workload sizes
identical. It measures CUDA's actual resident-block limit, abandons the change
if it fails to increase residency or spills, and otherwise runs five alternating
pairs per workload after excluded warmups. Cross-binary checkpoints and every
DP34 multiset must agree. Its 95% paired log-ratio t interval assumes independent,
approximately normal log ratios; five pairs remain exploratory evidence.

## Run

From `ecc2k130/`, using the existing AWS credential chain:

```sh
python3 aws/native_benchmark.py --profile --presigned-transfer \
  --launch-template YOUR_EXISTING_TEMPLATE --automatic-placement \
  --instance-type g7e.4xlarge --out /tmp/goal28-new-run
```

Use `--profile-counters-only` in place of `--profile` for the reduced-section
application-replay retry without a repeated grid screen. Use `--occupancy` for
the unprofiled register-budget comparison; it does not grant SYS_ADMIN.

The existing template supplies the AMI, security group and disposable root
volume. `g7.4xlarge` can collect a separate RTX PRO 4500 diagnostic when G7e
has no capacity; its rates cannot establish the RTX PRO 6000 target. The job
uses the existing approved benchmark prefix, adds a progress object, and
retains its 45-minute workload, 50-minute controller and independent 55-minute
shutdown limits. The profiler container receives SYS_ADMIN for GPU performance
counter access. No production template, IAM policy or campaign worker is changed.

Profiler versions, raw logs, `.nsys-rep` / `.ncu-rep` files, CSV summaries,
source and binary hashes, and unprofiled timing rows are retained. Package
versions are selected from NVIDIA's configured CUDA repository and recorded;
CUDA compilation remains pinned to 13.3.73. The controller prints progress
stages as they change and retrieves the final archive before termination.

Sources for profiler commands and their semantics:
[Nsight Systems User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/),
[Nsight Systems installation](https://docs.nvidia.com/nsight-systems/InstallationGuide/index.html),
[Nsight Compute CLI](https://docs.nvidia.com/nsight-compute/NsightComputeCli/index.html).
