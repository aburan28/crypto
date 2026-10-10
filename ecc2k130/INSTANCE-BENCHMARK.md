# Measure complete instance throughput

`benchmark_instance.py` launches one packed CUDA client per selected physical
GPU and counts completed scalar rho updates over a shared host time window.
Use this on the intended instance to measure the **40 billion updates/s**
target. GPU model and instance type are recorded separately.

## Run

From `ecc2k130/`, using a validated packed binary built for the host GPU:

```bash
python3 benchmark_instance.py --devices all \
  --out benchmarks/instance-40b-benchmark \
  --run-id-base 54000

python3 benchmark_instance.py --devices all \
  --out benchmarks/instance-40b-dp34 \
  --run-id-base 54100 --collect-dp34
```

Every output directory must be new. Select a subset with `--devices 0,1`, or
use comma-separated physical GPU UUIDs. `--binary` selects another validated
packed executable. `--batch` must match its compiled logical batch size
(default 16); it does not rebuild or alter the executable. Linux process
ownership and physical NVIDIA GPUs are required; MIG is not supported.

Defaults are 524,288 logical workers per GPU, batch 16, 1,024 steps per launch,
32 launches and three independent timed cohorts. Each timing client performs
274,877,906,944 complete updates. Before timing, each GPU performs a separate
small DP52 admission with 16 independent CPU replays. Use `--timeout` to bound
each cohort (default 900 seconds). Run IDs are distinct across all admissions
and repetitions within one invocation. Choose disjoint ID ranges for separate
invocations; IDs must fit 16 bits. Walk populations cannot exceed 2^32.

## What the rate means

```text
instance updates/s = sum(completed updates from every concurrent client)
                     / common host elapsed seconds
```

The window begins before the first process launch and ends after all processes
exit and observation/cleanup completes. It includes each timed client's setup,
launch skew and teardown. Independent admission checks precede that window.
Per-client rates remain diagnostics because their individual timing windows
differ. Hardware identity, physical GPU UUIDs, process ownership, exact work,
source/binary hashes, logs, DP files and elapsed endpoints remain in the output.

The driver rejects duplicate GPU assignments, repeated or overflowing run IDs,
geometry/work mismatches, missing CPU replays, dropped/truncated DP records,
worker errors, timeouts and observed foreign GPU processes. Finish compilers
before starting. Every child PID remains owned and unreaped through GPU
observations, preventing PID reuse from confusing interference checks. GPU
process sampling cannot exclude activity between observations. Failed cohorts
are retained and stop the experiment without replacement samples.

`summary.json` marks the absolute target met only if at least three complete
cohorts each reach 40 B/s. This is an absolute observation for that recorded
hardware and workload; relative speedup claims still require matched baseline
and candidate comparisons. Generic rho work is unchanged and full-DLP S is null.

## Current validation

[The saved validation](benchmarks/g7-40b-instance-driver/RESULTS.md) passes eight
simulated-worker tests and the available one-GPU checks. On g7.2xlarge / RTX PRO
4500, single setup-inclusive cohorts measured **5.329455 B/s** benchmark and
**5.323868 B/s** DP34. These method checks leave the accepted kernel records
unchanged. **Multi-GPU hardware and G7e throughput remain unmeasured.**
