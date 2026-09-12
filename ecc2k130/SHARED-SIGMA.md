# Shared walk-permutation masks

`PACKED_SHARED_SIGMA=1` is an experimental packed-kernel option. Its default
is 0, and it requires weighted-prefix mode 2 and the walk permutation network.
The RTX PRO 6000 preset selects it after a finite matched comparison measured
gains of 0.801% in the complete scalar benchmark and 0.954% with DP34 collection.
`RTX_PRO6000_SHARED_SIGMA=0` retains the global-memory control. The new public
Make/Modal integration still requires its first public-command GPU audit.

The paired Frobenius helper uses a 56-by-8 table of 32-bit masks. The candidate
copies all 448 words, or 1,792 bytes, into a separate shared-memory table once
per block. Every block thread participates in the copy and barrier before
inactive workers return. The arithmetic helper contains no barrier. Its
generated operations and mask indices match the original helper; the scalar,
global paired and inverse helpers remain separately available.

Rows remain contiguous. A fixed-row lookup with eight possible indices maps
to eight different banks, while equal indices read the same address. This is
an address model. The measured gain below does not establish bank-conflict
behavior or a performance-counter explanation.

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
system reservation, explaining the 2,816-byte compiled extent. The first
untimed runtime calibration observed `cudaFuncAttributes.sharedSizeBytes=1792`
and a separate device reservation of 1,024 bytes on device 0. The function
attribute excludes the reservation on this device/runtime. Control attributes
were 108 registers, zero local bytes and zero function shared bytes; candidate
attributes were 104 registers, zero local bytes and 1,792 function shared bytes.
Both reported two resident 256-thread blocks per SM across 188 SMs.

The [calibration payload](benchmarks/shared-sigma/runtime-calibration.json) is
hash-bound to the first 128-update client observation for each mode, the source,
and both executables. Later recorded client markers match those exact values.
The compiled extent, function attribute and reservation remain distinct fields;
the generic ranking helper does not assume platform-specific shared-memory sizes.

The dedicated `test-shared-sigma-cuda` target exercises the original and
selected helpers against independent routing. For example:

```sh
make test-shared-sigma-cuda BATCH=16 THREADS=256 \
  ARCH='-gencode arch=compute_120,code=sm_120' \
  PACKED_WEIGHTED_PREFIX=2 PACKED_PERM_SIGMA=3 PACKED_SHARED_SIGMA=1
```

Use mode 0 for the global control. The test passed on the GPU in both modes and covers 21 scenarios,
21,036 input pairs, 114 complete block-mask snapshots and 51,072 mask words.
It poisons shared storage before initialization, checks output guards and
fully inactive blocks, and interleaves all eight powers within warps in one
repeat. The [probe resolution](benchmarks/shared-sigma/probe-resolution.json)
records that additive coverage change. Both probe modes completed in the
frozen full-client comparison. The public audit forwards
the selected mode and all arithmetic flags to this target before integration
or timing. It requires the exact scenario, pair, snapshot and word counts and
one matching mode marker. This probe applies to WP2 with the walk network;
WP0/WP1 audits retain their existing arithmetic/integration paths.

The [complete comparison](benchmarks/shared-sigma/comparison.json) and
[independent review](benchmarks/shared-sigma/comparison-review.json) bind the
same kernel, arithmetic and probe source now used here. One RTX PRO 6000
Blackwell Server Edition ran native sm_120 code built with CUDA 13.3.73, driver
580.95.05. Both modes used B16/T256/min2/COMPACT1/WP2/CLMAD1/TILE256, with
385,024 workers and 6,160,384 scalar walks. Each timed row completed exactly
201,863,462,912 updates (1,024 steps and 32 launches), with zero dropped reports.

| Three-repeat measurement | Global median B/s | Shared median B/s | Median gain |
|---|---:|---:|---:|
| Complete scalar benchmark | 14.296584 | 14.411102 | 0.8010165% |
| DP34 collection | 13.960757 | 14.093912 | 0.9537807% |

Benchmark ranges were 14.284151–14.327940 B/s for the control and
14.407711–14.478495 B/s for shared sigma. Collection ranges were
13.960680–13.974133 B/s and 14.092497–14.103560 B/s. All three alternating pairs
favored shared sigma for each measurement. These are finite measurements on
one allocation, not a cross-device result or a population-level estimate.

The fixed panel completed 17 rows: two excluded warmups, a three-row
control/candidate/control screen, six benchmark confirmations and six collection
rows. The candidate narrowly exceeded the preset screen threshold of
14.5556964 B/s; control drift across the screen was −0.7297518%. Every fresh
collection produced 5,149 records (164,768 bytes), and all six sorted corpus
hashes matched. Timed rows used `--verify 0` after separate reference checks;
the archive retains collection counts and hashes rather than temporary corpus bytes.

Before timing, both modes passed all six arithmetic suites, 128 storage cases
covering 297,344 records, the shared-table probe, and full client replay/resume/
guard integration. Normalized state comparisons covered 64 and 16,384 scalar
slots. The checkpoint gate completed 26 successful children, two expected
worker-geometry rejections with preserved inputs, and 12 comparisons.

The [evidence index](benchmarks/shared-sigma/evidence-index.json) records all
archive hashes. The raw comparison SHA-256 is
`ce331c88a096637a23f9f5a4d632a535f1c5d89e2d112d31e41599f7446df406`;
the review SHA-256 is
`3073c6ee74924851293bc410e8cd05892743bbe51a59453d8549d5f269f5bc7c`.
The exact CPU prebuild receipt is embedded in the raw comparison, with its
original hash and reconstruction recorded in the index. Driver, source,
executable, complete native-code and calibration bindings stayed unchanged.

The preset and control are selected with:

```sh
make audit-rtx-pro6000
make audit-rtx-pro6000 RTX_PRO6000_SHARED_SIGMA=0
```

These public commands have not yet supplied a GPU audit of this integration.
Missing, duplicated or different shared-mode markers invalidate packed benchmark
and collection rows; applicable audits stop before integration/timing if the
shared-table probe fails. The previously validated public benchmark median is
14.472716 B/s. This paired gain does not establish a new public absolute-rate
record. The objective remains 26 B complete scalar iterations/s on one GPU,
approximately 1.804 times the candidate's measured benchmark median.
