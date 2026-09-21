> Archived experiment: the source inputs matching the saved provenance are
> preserved in `sources/`. The reproduction commands below describe the benchmarked
> candidate checkout; this UI PR does not install its runtime API or Make targets.

# Bounded CPU/GPU scheduling: local stage measurements

Class: **engineering, stage diagnostic**. The two-slot pipeline overlaps CPU
preparation and result handling with GPU dense RREF and host/device transfers.
The existing reduction arithmetic is unchanged. This is local pinned-memory
staging, not network RDMA; no RDMA device is present on the validation host.

## Boundary and falsification target

The [frozen contract](contract.json), established before running the comparison,
requires the same matrices, IDs, ranks, pivot lists and canonical-result digests
as the independently CPU-verified serial reducer. Any missing, duplicate or
incorrect result rejects the run. Scheduling alone cannot reduce the underlying
algebraic work or move the generic-group floor. The reference is the one-slot
serialized version of the same API; the candidate uses two slots. This isolates
overlap from changes in kernels, arithmetic, pinned-memory policy or CPU cutoff.
The old `bench throughput` is not a valid transport reference because it omits
transfers, output retrieval, CPU work and setup from its timer.

This is the equivalent-suite exception to the WDSat protocol: dense RREF over a
small prime field does not consume binary ANF solver inputs. It follows the
[parent accounting contract](../../../../research/index_calculus_baseline_20260914/ec_index_calculus_contract.json)
by retaining unavailable whole-DLP operation counts, rho ratios and floor ratios
as null. No full-DLP runtime or operation improvement is claimed. The API is not
yet integrated into Rust residual decomposition.

## Frozen comparison

RTX PRO 4500 Blackwell Server Edition, driver 595.91.07, `sm_120`.
Prime 2083, shape 226×286; 65 jobs per run, batch 16, CPU cutoff 4.
Each run reduces 64 matrices on GPU and the one-job tail on CPU. Three training
seeds and three fresh holdout seeds, three repetitions each, two variants:
36 successful runs and 2,340 verified matrix results (390 distinct inputs, each
replayed six times). The suite includes zero, full-rank, rank-deficient and
leading-zero-column matrices. Every matrix digest, rank, full pivot list and
exactly-once ID check passed. No failed runs were omitted.

One unit throughout this supplementary table: **seconds per cold stage
invocation**, including allocation/context setup, preparation, copies, reduction,
result checking and teardown. Ratios are medians of *paired* candidate/reference
times, not ratios of the displayed marginal medians.

| Variant | Split | Median seconds | Paired time / serialized | Common ops / rho | Common ops / floor | Correctness | Class |
|---|---|---:|---:|---|---|---|---|
| Serialized, one slot | Training | 0.323409 | 1.000000 | unmeasured | unmeasured | 585/585 | engineering reference |
| Overlap, two slots | Training | 0.306742 | 0.940133 | unmeasured | unmeasured | 585/585 | engineering |
| Serialized, one slot | Holdout | 0.321892 | 1.000000 | unmeasured | unmeasured | 585/585 | engineering reference |
| Overlap, two slots | Holdout | 0.307124 | 0.949346 | unmeasured | unmeasured | 585/585 | engineering |

Evidence: [raw paired runs](run-01/raw.jsonl), [summary](run-01/summary.json),
[source and executable hashes](run-01/provenance.json).

The holdout paired median is descriptively 5.07% lower. This is **not a statistical
speedup claim**: the host was not exclusive and no paired 95% confidence interval
was established. Cold setup dominates (~0.225 seconds median for overlap).
The serial CPU tail takes ~0.052 seconds; its work overlaps GPU activity in the
candidate. `min_gpu_batch` is configurable and was held fixed, not calibrated
as an optimal routing policy. Setting it to 1 avoids CPU tails when appropriate.
This comparison does not establish the best cutoff or that a CPU tail is faster
than sending the tail to an already initialized GPU.

CPU preparation and checking each take ~0.004 seconds. GPU reduction remains
~0.045 seconds: the kernel itself was not accelerated. Phase times overlap and
must not be summed. CPU oracle construction occurs before the timed invocation;
`process_wall_s` retains the full process time including that oracle and startup.
There is no measured operation-unit conversion, matched rho, real-residual yield,
full-DLP result, or scaling exponent. Existing attack verdicts remain unchanged.

## Reproduce and validate

From `gpu/macaulay`:

```sh
make test bench pipeline test_pipeline ARCH=sm_120
./bench selftest 8
./test_pipeline
CUDA_VISIBLE_DEVICES='' ./pipeline 3 16 2 12345 4
compute-sanitizer --tool memcheck --error-exitcode 99 ./pipeline 35 16 2 12345 4
compute-sanitizer --tool memcheck --error-exitcode 99 ./test_pipeline
compute-sanitizer --tool synccheck --error-exitcode 99 ./pipeline 35 16 2 20261 4
python3 benchmark_pipeline.py --output /tmp/macaulay-fresh-run
```

The runner refuses an existing output directory. All 36 raw records include
commands, stdout, stderr, exit status, phase times and job placement. Failures and
timeouts are retained and reject the summary. CPU/Python arithmetic and serial
RREF tests passed at both the default shape and p=271, 64×80. GPU selftests and
the mixed CPU/GPU pipeline passed at both shapes. CUDA memory and synchronization
checks reported zero errors, including callback exceptions with work in flight.
