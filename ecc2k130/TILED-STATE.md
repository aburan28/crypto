# Tiled packed walk state

`PACKED_STATE_TILE=256` groups the packed coordinate and batch scratch fields
by blocks of 256 workers. The RTX PRO 6000 preset enables it; the general
default remains `0`. Tiling requires `THREADS=256`, polynomial coordinate
storage, denominator caching and polynomial chains.

The layout changes addressing and storage only. The generated-product
arithmetic, walk rule, logical worker IDs, scalar iteration count,
distinguished-point records and checkpoint version remain the same.

## Build and run

From `ecc2k130`, with Modal installed and authenticated:

```sh
make bench-rtx-pro6000
make audit-rtx-pro6000
```

The preset keeps batch 32, 256 threads per block, minBlocks 2, 192,512
workers and CUDA 13.3.1. The audit checks GPU arithmetic and full client
integration before three benchmarks and three DP34 collection runs.
`MODAL=/path/to/modal` selects a specific client installation.

For a custom build, the Make setting is `PACKED_STATE_TILE=0|256`; the Modal
environment setting is `ECC_PACKED_STATE_TILE=0|256`. An enabled Modal image
bakes a compatible 256-thread binary. Build reuse includes the tile mode,
and incompatible block sizes fail before compilation. Benchmark/audit
results require exactly one matching runtime tile marker; missing, wrong
or duplicated markers invalidate the rate.

## Physical layout and checkpoints

The untiled field index is `(slot * 5 + word) * workers + tid`. The tiled
index is `((tid / 256) * batch * 5 + slot * 5 + word) * 256 + tid % 256`.
Each warp still accesses adjacent words. X, Y, denominator and prefix
allocations round the physical worker count up to a multiple of 256.
Logical walk identity remains `slot * workers + tid`; seeds, dead flags and
iteration metadata retain their logical indexing.

Checkpoints remain version 2 with logical-size coordinates in the normal
basis and untiled order. Export removes physical padding before basis
conversion. Import converts the logical field and tiles it into zero-padded
storage. X and Y reset their staging sizes independently. The common engine
distinguishes physical transfer size from logical checkpoint size, with the
original behavior as the default for other engines.

Resuming across tile modes requires the checkpoint's original worker and
batch counts. Supply that worker count explicitly; changing block size can
change the automatically selected count. No checkpoint migration or format
change is required for a compatible packed checkpoint.

## Controlled GPU measurements

The [combined comparison](benchmarks/tiled-state/comparison.json) uses
generated-product arithmetic in both modes on one RTX PRO 6000 Blackwell
Server Edition with 188 SMs, driver 580.95.05 and CUDA 13.3.73. Only the tile
mode changes. Both use B32/T256/minBlocks2 and 192,512 workers. Every timed
sample completes **201,863,462,912 scalar updates**, through 32 launches of
1,024 steps.

| Workload | Untiled median B/s | Tiled median B/s | Gain |
|---|---:|---:|---:|
| Complete scalar benchmark | 6.852888 | **6.994745** | **2.070032%** |
| DP34 collection | 6.766691 | **6.894434** | **1.887821%** |

Benchmark ranges across three confirmations are 6.847417–6.855230 B/s
untiled and 6.991535–6.995327 B/s tiled. Collection ranges are
6.764705–6.769657 and 6.894273–6.897658 B/s. Every paired repetition favors
the tiled candidate. Two warmups and the initial control/candidate/control
screen are excluded from these confirmation medians.

All six collections contain 5,149 records, 164,768 bytes and zero drops,
with sorted record hash
`ab237b6352380547fd37fcdc9e2aa83f6ae1a718b842590b9a19e4224e5d336b`.
The [independent artifact audit](benchmarks/tiled-state/comparison-review.json)
verifies all 17 timed rows, source and binary bindings, complete walk/init
instruction encodings and resources, arithmetic/client checks, 64- and
16,384-slot normalized states, and both checkpoint resume directions at
8 and 257 logical workers. Temporary checkpoint/corpus contents are assessed
through retained producer receipts and hashes; those temporary files were
not retained for separate local decoding.

The earlier [native-arithmetic comparison](benchmarks/tiled-state/native-arithmetic-comparison.json)
measured 6.891133 to 7.028159 B/s benchmark (+1.988439%) and 6.790145 to
6.923774 B/s collection (+1.967984%). It also passed the full checks, including
matching collection hashes. That comparison used a separate allocation and
the handwritten multiplier. Do not combine its absolute rates with the
generated-product comparison to estimate another gain.

These measurements do not establish the 15 B/s target or performance on
other GPUs/configurations.

## Compiler and validation evidence

On the frozen generated-product source, the untiled control reproduces the
previous public G1 cubin exactly. The tiled walk uses 128 registers, zero
stack bytes, no shared memory and no caller/callee spills. Its complete
common-path instruction count is 4,204.5625 per scalar update, versus
4,249.5 untiled. Field-product multiplicities and active global transfer
bytes are unchanged. These counts are static code observations, not cycle
counts or measurements of DRAM traffic, cache hits or TLB behavior.

The [host interaction review](benchmarks/tiled-state/host-interaction-review.json)
binds the exact composed source. Both modes passed actual-kernel/scalar
reference cases covering reports, reseeds, guards and zero inputs under
Clang/UBSan. Focused checkpoint tests include 12 byte comparisons,
33,920 resumed updates, 84 init/reseed boundary calls and four comparisons
against retained native checkpoint goldens. The public kernel, arithmetic
and generator files preserve those tested bytes; wrapper/preset changes are
validated separately by `make check-cli` and the native audit command.

Combined raw result SHA256:
`87e60bd3da41f80066bfed5057b2830e16c24743c40572c5172a3ca473d12b31`.
