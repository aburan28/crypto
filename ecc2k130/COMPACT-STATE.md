# Compact physical field storage

`PACKED_COMPACT_STATE=1` enables an experimental storage layout for the packed
CUDA backend. The default is 0. It requires `PACKED_STATE_TILE=256`, polynomial
state, polynomial chains and the denominator cache. The tested arithmetic
configuration uses weighted-prefix mode 2 and batch 16.

This prototype has passed host validation and a two-build CUDA compiler
comparison. Full device validation and complete-walk timing are pending.
The current measured RTX preset is documented in [RTX-PRO6000.md](RTX-PRO6000.md).

## Representation

Each field uses a 256-worker tile per batch slot:

| Region | Bytes per tile | Access per worker |
|---|---:|---|
| Four low limbs | 4,096 | One aligned 16-byte record |
| High limb and optional jump tag | 256 | One unsigned byte |
| Total | **4,352** | **17 bytes** |

The previous tile occupies 5,120 bytes. The new layout reduces each field
allocation by 15%. At 385,024 workers and batch 16, the four field arrays
occupy 399.5 MiB instead of 470 MiB. The separate dead/seed/start arrays
occupy 117.5 MiB. These are allocation sizes, not measured cache or memory
traffic statistics.

For `q = (tid / 256) * BATCH + slot`, a worker's low limbs start at byte
`4352*q + 16*(tid % 256)` and its high byte is at
`4352*q + 4096 + tid % 256`. The regions and worker records are disjoint.
Coordinates and weighted prefixes use high values 0 through 7. Cached
denominators include three jump bits and use values 0 through 63; storage
preserves all six bits.

Device access uses `uint4` for low limbs and an unsigned byte for the tail.
Host checkpoint staging uses `memcpy` for low limbs, so it does not depend
on vector alignment or vector object lifetimes in a host `vector<unsigned>`.

Logical fields contain `5*BATCH*workers` words. Physical transfers contain
`1088*BATCH*ceil(workers/256)` opaque words. Checkpoint version 2 stores
logical normal-basis coordinates and the existing metadata. Export decodes
physical storage before the basis conversion; import converts the logical
coordinates before encoding into a zero-filled physical buffer. The generic
checkpoint methods resize staging for each X/Y field, including geometries
where the physical buffer is smaller than the logical payload.

## Validation so far

The [host receipt](benchmarks/compact-state/host-review.json) records:

- 400 storage-codec cases and 1,493,328 records, including all 64 tail values,
  padding, canaries and shifted host storage.
- 30 actual host kernel runs and 15 comparisons between layouts across
  batches 1/4/16/32 and partial-tile boundaries. These include 9,590 scalar
  reference steps, 130 zero-denominator batch steps and 632 report replays.
- Actual checkpoint methods and field hooks under checked host transfer
  shims: 18 byte comparisons, 54 malformed-input rejections before writes,
  and resumed steps at 8, 256 and 257 workers.
- 16 tests of the actual automatic-worker method with injected memory and
  occupancy responses, 14 default preprocessing comparisons and seven
  invalid-configuration checks.

An [independent source and retained-evidence review](benchmarks/compact-state/host-independent-review.json)
passed. These host checks use serial CUDA emulation and software field
arithmetic; they do not establish device concurrency or device correctness.
Size arithmetic was checked for 64-bit `size_t` and the selected batches.

The [CUDA compiler comparison](benchmarks/compact-state/compiler-review.json)
used CUDA 13.3.73, sm_120, B16/T256/minBlocks2, weighted mode 2 and native
carryless arithmetic. The mode-0 whole cubin equals the previous weighted
control. Walk registers were 100/108 for modes 0/1, with 98 initialization
registers in both. Reported stack, local, shared and spill counts were zero;
the inspected caller/helper code contains no local-memory instructions.

The compact walk caller contains eight 128-bit loads, eight unsigned-byte
loads, five 128-bit stores and five unsigned-byte stores. Initialization
contains two of each store width. These are actual native-code observations.
Static caller instruction counts and register counts do not predict a speedup.

## Synthetic device storage test

The new `test-packed-storage-cuda` Make target compiles and runs the actual
storage accessors. For example, on a CUDA GPU with a matching `ARCH` setting:

```sh
make test-packed-storage-cuda BATCH=16 THREADS=256 \
  ARCH='-gencode arch=compute_120,code=sm_120' \
  PACKED_STATE_TILE=256 PACKED_POLY_STATE=1 PACKED_POLY_CHAIN=1 \
  PACKED_CACHE_DENOM=1 PACKED_COMPACT_STATE=1
```

Use `PACKED_COMPACT_STATE=0` for its control. At batch 16, each version
checks 128 cases and 297,344 records. The independently encoded physical
image includes untouched padding and surrounding guards; separate logical
reads include output guards. Tail patterns identify both worker and slot
across rounds and sweep all 64 legal values. The
[static test review](benchmarks/compact-state/storage-test-static-review.json)
records the address-signature checks. The device test has not yet run.

The next measurement uses the same logical population and complete scalar
iteration budget as the current weighted preset, after device arithmetic,
storage, client and checkpoint checks. This prototype is currently exposed
through direct Make builds; Modal option wiring and a preset change await
a successful complete-walk comparison. The 15 B/s target remains unachieved.
