# Shared walk-permutation masks

`PACKED_SHARED_SIGMA=1` is an experimental packed-kernel option. Its default
is 0, and it requires weighted-prefix mode 2 and the walk permutation network.
The measured RTX preset continues to use the global-memory helper while this
candidate undergoes device validation and timing.

The paired Frobenius helper uses a 56-by-8 table of 32-bit masks. The candidate
copies all 448 words, or 1,792 bytes, into a separate shared-memory table once
per block. Every block thread participates in the copy and barrier before
inactive workers return. The arithmetic helper contains no barrier. Its
generated operations and mask indices match the original helper; the scalar,
global paired and inverse helpers remain separately available.

Rows remain contiguous. A fixed-row lookup with eight possible indices maps
to eight different banks, while equal indices read the same address. This is
an address model; actual performance still requires measurement.

The [host/source review](benchmarks/shared-sigma/source-host-probe-review.json)
passed 30 actual host kernel runs and 15 comparisons, including coordinates,
reports, seed/start/dead data and raw denominator/prefix scratch. It records
9,590 reference steps, 130 zero-batch steps and 632 report replays. Paired
helper checks cover 6,240 pairs in each of three ordinary/emulated compiler
branches. The actual copy loop was also executed sequentially for 12 block
widths with a counted dummy barrier. Host checks model post-barrier values;
they do not establish CUDA concurrency or memory ordering.

The [independent compiler review](benchmarks/shared-sigma/compiler-review.json)
checked two CUDA 13.3.73 sm_120 wrapper builds:

| Compiler observation | Global control | Shared candidate |
|---|---:|---:|
| Walk registers | 108 | 104 |
| Init registers | 98 | 98 |
| Compiled walk shared extent | 0 bytes | 2,816 bytes |
| Paired-helper mask loads | 56 global read-only | 56 shared |
| Paired-helper non-NOP instructions | 505 | 502 |

Reported stack/local/spill values are zero. The control cubin matches the
previous compact implementation exactly, and init code is identical. The
candidate's 1,792-byte table begins at offset `0x400`, after a 1,024-byte
system reservation, explaining the 2,816-byte compiled extent. The actual
runtime `cudaFuncAttributes.sharedSizeBytes` value has not yet been measured;
runtime validation must record it separately from the device reservation.

The dedicated `test-shared-sigma-cuda` target exercises the original and
selected helpers against independent routing. For example:

```sh
make test-shared-sigma-cuda BATCH=16 THREADS=256 \
  ARCH='-gencode arch=compute_120,code=sm_120' \
  PACKED_WEIGHTED_PREFIX=2 PACKED_PERM_SIGMA=3 PACKED_SHARED_SIGMA=1
```

Use mode 0 for the global control. The planned test covers 21 scenarios,
21,036 input pairs, 114 complete block-mask snapshots and 51,072 mask words.
It poisons shared storage before initialization, checks output guards and
fully inactive blocks, and interleaves all eight powers within warps in one
repeat. The [probe resolution](benchmarks/shared-sigma/probe-resolution.json)
records that additive coverage change. This probe has been reviewed but has
not yet been compiled or run on the GPU.

No throughput improvement is established by these host and compiler results.
The current objective is 26 B complete scalar iterations/s on one RTX PRO
6000; the validated public median remains 14.472716 B/s. Public Modal option
wiring and preset selection will follow only if the candidate qualifies in
the complete-walk comparison.
