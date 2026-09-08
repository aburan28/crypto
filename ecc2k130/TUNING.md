# Measuring ECC2K-130 throughput changes

The reference measurement supplied on 2026-09-08 was **857.163 M iterations/s**
on an NVIDIA RTX PRO 6000 Blackwell Server Edition (`sm_120`), with
`leaf=0`, batch 32, 128 threads/block and two requested resident blocks/SM.
The current generator resolves leaf 0 to 66 for GF(2^131).

The experimental options below have **no claimed GPU throughput gain**. They
preserve the walk function and remain disabled by default. CPU arithmetic and
state-equivalence checks do not replace CUDA correctness or performance tests.

## What this implementation lets you compare

| Option | Default | Purpose |
|---|---|---|
| `--stream-karat` / `STREAM_KARAT=1` | off | Compute and consume one Karatsuba subproduct at a time |
| `--smem-spill` / `SMEM_SPILL=1` | off | CUDA 13+ opt-in shared-memory register spilling |
| `--global-cg` / `GLOBAL_CG=1` | off | Compile global loads with the L2-only cache policy |
| `--prefer-l1` | off | Request more L1 cache for the walk kernel |
| `--workers N` | automatic | Set total CUDA threads independently of block size/register budget |
| `--repeats N` | 3 | Rank by median completed-run throughput; retain min/max and every sample |
| `ECC_CUDA_VERSION` | `12.8.1` | Select the Modal toolkit image before import/build |

The three compile flags are accepted by `make gpu`. The CLI flags are accepted
by `bench`, `autotune` and, for the arithmetic/cache options, `validate`.
`--threads` still controls **block size** in Modal; `--workers` controls **total
threads**. The C++ binary uses `--threads` for total threads.

At the outer 131-to-66 split, the stream schedule changes the source scratch
inventory from `4*66 + 3*131 = 657` words to `2*66 + 131 = 263` words. This
saves 1,576 bytes of source buffers per thread at that split. CUDA can change
lifetimes and reuse storage, so these are not measured register/stack savings.
The schedule also changes the output accumulation order and may add memory
accesses. Keep the old schedule as the control until a device selects a winner.

One small common-path change is unconditional: the last cumulative-product
store is omitted because inversion consumes the product directly and reverse
traversal only reads earlier prefixes. Batch size 1 needs no prefix stores.

## Baseline, arithmetic candidate, and leaf/batch sweep

From this directory:

```bash
ECC_GPU=RTX-PRO-6000 modal run modal_app.py::bench --repeats 3
ECC_GPU=RTX-PRO-6000 modal run modal_app.py::validate --stream-karat
ECC_GPU=RTX-PRO-6000 modal run modal_app.py::bench --stream-karat --repeats 3

ECC_GPU=RTX-PRO-6000 modal run modal_app.py::autotune \
  --configs '0:16:128:2,0:32:128:2,0:64:128:2,33:16:128:2,33:32:128:2,33:64:128:2,131:16:128:2,131:32:128:2,131:64:128:2' \
  --repeats 3
```

The baseline and candidate commands can land on different physical devices;
their returned GPU UUID/clock metadata exposes that. Repeat finalist controls
in alternating order on the same device allocation when practical. The
autotuner repeats each build consecutively to avoid recompiling for every
sample. Its default three repetitions increase benchmark time versus the old
one-sample behavior; `--repeats 1` is available for screening.

Leaf 131 removes the outer C++ recursion but still uses a generated Karatsuba
circuit internally. It trades call/local-array traffic for instruction count
and a larger generated body. It is not a recommended default.

For grid tuning, use the SM count printed by the binary. Compare total worker
counts of `SM_count*128` and `SM_count*256`, both with `--threads 128
--min-blocks 2`. This changes the launch footprint without changing the register
budget. Do not assume maximum occupancy gives maximum throughput.

The arithmetic cost at batch B is `(5B+5)/B = 5+5/B` field multiplications per
scalar iteration: B−1 forward multiplications, eight for inversion, and 4B−2
reverse multiplications. Moving batch 32 to 64 saves only 1.52% of multiplication
work. Larger effects reflect storage, cache, scheduling, or overhead.

## Compiler and cache experiments

Compile without renting a GPU:

```bash
ECC_GPU=RTX-PRO-6000 modal run modal_app.py::compile_check --stream-karat
ECC_GPU=RTX-PRO-6000 ECC_CUDA_VERSION=13.0.2 \
  modal run modal_app.py::compile_check --smem-spill
```

Then compare CUDA 13 with and without the feature:

```bash
ECC_GPU=RTX-PRO-6000 ECC_CUDA_VERSION=13.0.2 \
  modal run modal_app.py::bench --repeats 3
ECC_GPU=RTX-PRO-6000 ECC_CUDA_VERSION=13.0.2 \
  modal run modal_app.py::validate --smem-spill
ECC_GPU=RTX-PRO-6000 ECC_CUDA_VERSION=13.0.2 \
  modal run modal_app.py::bench --smem-spill --repeats 3
```

The pragma targets register spills, not every addressable local array or device
function frame. Check shared bytes/block, local bytes/thread, registers/thread,
and resulting occupancy. Whole-program compilation and explicit launch bounds
are retained; dynamic shared allocation is not used. Older toolkits are rejected
when this option is requested. [NVIDIA shared-memory spilling documentation](https://developer.nvidia.com/blog/how-to-improve-cuda-kernel-performance-with-shared-memory-register-spilling/).

Test `--prefer-l1` and `--global-cg` individually before combining them. The
hypothesis is that streaming global field state evicts reused local arithmetic
intermediates. However, the reverse pass can reuse field state too, so L1 bypass
can lose. Shared-memory spilling and cache carveout also interact.

## Profile the exact binary

On a GPU machine with Nsight Compute and access to performance counters:

```bash
ncu --kernel-name 'regex:eccWalkKernel' \
  --launch-skip 2 --launch-count 1 \
  --section SpeedOfLight --section LaunchStats \
  --section Occupancy --section MemoryWorkloadAnalysis \
  --section SchedulerStats --section WarpStateStats \
  --section SourceCounters --section InstructionStats \
  -o ecc2k130-baseline \
  ./ecc2k130 --curve 131 --bench --steps 64 --launches 4 --verify 0
```

This profiles a later walk launch rather than initialization. Save the profiler
version and clock/cache-control settings; replay changes the execution
environment. Use ordinary, unprofiled runs for throughput comparisons. Available
sections can be checked with `ncu --list-sections`. Some GPU providers restrict
counter access; compile success does not establish profiling access.
[Nsight Compute CLI documentation](https://docs.nvidia.com/nsight-compute/NsightComputeCli/index.html).

Look for actual DRAM saturation, local-memory traffic at L1/L2, load/store issue
pressure, dependency stalls, and instruction-fetch stalls. Static instruction
counts, source array sizes and ptxas spill reports do not alone distinguish
these bottlenecks.

In the current common path, x/y/pchain perform `7−2/B` logical field transfers
per slot after removal of the unused last prefix store. Before that change the
count was `7−1/B`, or 114.113 bytes/scalar iteration at B=32. That was about
97.814 GB/s of logical field-array requests at 857.163 M it/s, excluding local
arrays, spills and rare paths. It does not establish DRAM saturation. The old
README estimate of roughly 98 bytes/iteration omitted one field-state read.

## Evidence and correctness

Benchmark JSON contains the source digest, generated leaf, binary digest,
compiler version, ptxas build log when rebuilt, GPU metadata, command, return
code and full stdout/stderr for every repetition. The kernel prints resource
attributes even when worker count is explicit. `rate` is the median, in M it/s.

Nonzero exits, missing/invalid final rates, interrupted runs and mismatches
invalidate a repetition. A candidate with a failed repetition is not eligible
for `best`; an entirely failed sweep returns `best: null` and the CLI fails.
Every sweep receives a timestamped JSON file in `/data/autotune`, while the
legacy per-GPU filename remains a pointer in practice to the latest contents.

Initialization precedes the throughput timer. Fresh benchmark rates count
scalar bit lanes. For resumed searches, both throughput numerators subtract
the restored iteration base, so historical work no longer inflates the rate.
`--bench` suppresses almost all reports via cutoff zero; it does not measure
full search/report/reseed/persistence throughput or remove the reporting code.

`validate` builds the requested candidate and retains the small-curve GPU
end-to-end tests. It also requires actual report replay on both GF(2^97) and
GF(2^131), including at least 16 verified reports, no mismatches, zero drops and
a successful process exit. Failure raises an error rather than printing a
plausible validation summary and returning success.

Local checks:

```bash
make check-cli
make test-schedule
make test
# Build baseline/candidate binaries with matching batch and lane settings first:
python3 codegen/testwalkvariants.py BASELINE_BINARY CANDIDATE_BINARY
```

The convolution checks cover odd/even splits, 32/64-bit words and guarded output
boundaries. The paired-client check requires identical resumed 131-bit walk
state, correct current-run iteration counts, 131-bit report replay, and planted
discrete logs on both binaries. These tests can run on CPU or GPU clients; the
tested backend must be reported with the result.
