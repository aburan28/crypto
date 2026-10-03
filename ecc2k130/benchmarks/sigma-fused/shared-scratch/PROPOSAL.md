# Fused-sigma shared scratch: static source proposal

Status: **static feasibility only; implementation and GPU measurement are not
part of this change.**  The source boundary is current `main` commit
`18ee85c1075010ce83be60a4a02636e9aa4ac83a`.  The native result is
[`model.json`](model.json).

## Question and conclusion

The selected counter-free sigma walk is B16/T256/minBlocks2 with
`SIGMA_FUSED=1`.  Its fused step writes a polynomial prefix value to `pchain`
and a tagged polynomial denominator to `denominators`, then the following
reverse pass reloads both values.  Can the first two, three or four slots owned
by each CUDA thread keep those two ephemeral fields in per-block shared memory
without changing the walk or losing the selected two-block geometry?

The source and native ledger say **yes, as an implementation candidate**:

- four cached slots remove exactly **17 logical compact-field bytes per scalar
  update** from the hot global path;
- the replacement is **20 logical shared-memory bytes, or five 32-bit shared
  operations, per scalar update**;
- the four-slot array adds 40,960 bytes per block, for 42,752 bytes of reported
  function shared memory and a 43,776-byte compiled extent when the previously
  observed 1,024 driver-reserved bytes are charged;
- two such extents use 87,552 of the recorded 102,400 bytes/SM, while the
  retained 126-register report remains below the 128-register two-block launch
  ceiling; and
- five full-word slots would use 108,032 shared bytes for two blocks and is the
  frozen capacity stop.

This establishes layout and capacity feasibility.  It does **not** establish a
speedup, a bandwidth limit, actual cache transactions, retained occupancy after
compilation, or progress to 26 B/s.

## Source audit and launch boundaries

The selected path is the `ECC_SIGMA_FUSED` branch in
`include/packedkernels.cuh`:

1. The launch prologue visits slots `0..15`.  `sigmaFusedSelect` writes one
   `pchain` value and one tagged denominator for every slot.
2. Step zero visits the slots in descending order and reads those values.  A
   non-final reverse pass immediately calls `sigmaFusedSelect` on the new point
   and writes the next generation while the point remains in registers.
3. Slot order alternates.  The next step therefore consumes the prefix chain
   in the reverse of the order in which it was built.
4. The final reverse pass reads the last scratch generation and deliberately
   performs no next selection.  Nothing consumes `pchain` or `denominators`
   after the kernel returns.

For a launch of `S > 0` steps, each scratch field and slot consequently has
exactly `S` producers and `S` consumers: the prologue replaces the producer
omitted by the final step.  The native checker exercises `S =
1,2,3,4,7,16,95,1024` and rejects a read without the same launch's preceding
producer.

Checkpoint persistence does not create another boundary.  `CudaEngine::save`
and `restore` serialize `x`, `y`, `dead`, optional witness counters and the
backend's lane arrays.  They do not serialize `pchain` or `denominators`.
Keeping the existing full global allocations while leaving the cached slots
stale at kernel exit is therefore safe: the next launch overwrites every cached
slot in shared memory before reading it.  Allocation shrinking and L2-window
changes are outside this proposal.

## Logical traffic and capacity ledger

Compact global fields use a 16-byte low record and one tail byte.  The shared
proposal deliberately uses five full `uint32_t` words so the denominator's
three jump bits in the top word are preserved without packing.  The model
counts source-visible logical field bytes, not cache-line or DRAM bytes.

At batch `B = 16` and `S` steps:

```
coordinate global bytes = 34*B + 68*B*S
scratch global bytes    = 68*(B-C)*S
shared scratch bytes    = 80*C*S
updates                 = B*S
```

The `34*B` term is the launch prologue's `x/y` read.  There is no scratch
boundary term.  At the registered `S = 1024`, the exact table is:

| cached slots `C` | scratch bytes/block | function shared | compiled extent | two-block extent | static block bound | global field B/update | saved B/update | replacement shared B/update | global / control |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 1,792 | 2,816 | 5,632 | 2 (registers) | 136.033203125 | 0 | 0 | 1.000000 |
| 2 | 20,480 | 22,272 | 23,296 | 46,592 | 2 | 127.533203125 | 8.50 | 10 | 0.937515 |
| 3 | 30,720 | 32,512 | 33,536 | 67,072 | 2 | 123.283203125 | 12.75 | 15 | 0.906273 |
| 4 | 40,960 | 42,752 | 43,776 | 87,552 | 2 | 119.033203125 | **17.00** | 20 | 0.875031 |
| 5, stop only | 51,200 | 52,992 | 54,016 | **108,032** | **1** | 114.783203125 | 21.25 | 25 | 0.843788 |

The active-block entries are static upper bounds using the retained 126
registers/thread, 65,536 registers/SM, 1,536 threads/SM, 102,400 shared
bytes/SM and the previously observed 1,024 driver-reserved bytes/block.  The
implementation gate must query the occupancy API on the actual RTX PRO 6000;
this table cannot substitute for that measurement.

The baseline 136.033 bytes cover only `x`, `y`, `pchain` and the denominator.
The steady full-step value is exactly 136 bytes/update; 0.033203125 is the
1,024-step launch's amortized `x/y` prologue.  The four-slot values are 119 and
0.033203125 respectively.  These values omit `dead`, seeds, iteration
metadata, reports, instruction fetches, cache-line amplification and writeback
policy.  A saved logical byte is not a saved DRAM byte.  The selected compact
layout is coalesced and may already hit L1 or L2.

## Proposed source shape

The follow-up implementation should add one closed knob,
`SIGMA_FUSED_SHARED_SLOTS`, with values `0`, `2`, `3` and `4`.  Nonzero values
must be rejected unless the exact B16/T256/minBlocks2 compact, shared-sigma,
one-chain fused path is compiled.

The kernel should declare a full-word SoA:

```c++
__shared__ uint32_t scratch[2][SIGMA_FUSED_SHARED_SLOTS][5][ECC_THREADS];
```

The last dimension is the thread index.  For a fixed field, slot and word, a
warp accesses 32 adjacent words and covers the 32 shared banks exactly once.
The model enumerates every owner and proves there is no alias.  One field is
`pchain`; the other is the tagged denominator.

Two small device helpers should route loads and stores for `slot < C` to that
thread's shared words and leave later slots on the existing compact global
path.  `sigmaFusedSelect` and the reverse pass must use the same helpers.  The
condition is block-uniform for a given loop iteration.  No scratch
initialization or inter-thread barrier is required: the producing and consuming
instructions have the same CUDA thread owner.  The existing shared-sigma mask
initialization barrier remains before the partial-block early return, and
inactive threads never access scratch.

The implementation should retain the global `pchain` and denominator
allocations, engine setup, checkpoint version, point arithmetic, prefix order,
reporting and all launch boundaries.  It should add a runtime marker for the
cached-slot count and exact static shared bytes.  Removing or resizing global
allocations would be a separate experiment.

## Risks and existing evidence

The exchange is not free.  The four-slot arm replaces compact vector-plus-tail
global operations with five scalar shared loads/stores per update.  Those
operations share the MIO path with the existing 1,792-byte sigma-mask table.
Forty-one additional KiB of static shared memory per block can also reduce the
unified L1 capacity available to `x`, `y` and metadata even while two blocks
remain resident.  The kernel is already at 126 registers/thread under a
two-block launch bound.  There are at most two nominal registers/thread before
the 128-register capacity edge, and hardware allocation granularity may already
round the reported 126 to that edge; new address state may create spills.

Repository history contains no exact test of fused `pchain` plus denominator
scratch.  Related evidence is mixed:

- the selected sigma fusion removed one `x/y` load pass and measured a 1.028257
  paired median, so field traffic has mattered on this path;
- the fused one-knob star's persisting-L2 arm measured a 0.989916 geometric
  ratio, so cache placement alone did not pay;
- the table walk's tag-denominator route removed 34 logical bytes/update and
  gained about 1.5% on a different kernel and geometry; and
- historical commit `98efc958a88b171523c92619e7e054a74ad16cc9`
  tested a compact shared-X cache in a different split/block-inverse G7 path.
  Its B32-local gain was 1.47--1.75%, but B16 stayed faster.  That experiment
  used a legacy Python orchestration layer, is not on the current source path,
  and is context only under the repository's native-only rule.

None of these results predicts this arm.  The proposed panel in
[`PROTOCOL.md`](PROTOCOL.md) treats a speedup as unknown until matched native
correctness, resource and timing gates pass.

## Stop and next gate

This work stops at `STATIC_FEASIBLE_IMPLEMENTATION_NOT_MEASURED`.  No CUDA
candidate, binary, GPU process or throughput sample is produced here.

The next gate is a separate implementation PR.  Before it may allocate a GPU,
all four arms must compile to native `sm_120`, report their exact registered
shared extents, retain two active blocks/SM, use at most 128 registers/thread,
and have zero stack, local memory and ptxas spill loads/stores.  A native device
test must then prove per-thread scratch ownership, top-word preservation,
partial-block safety and equality with the compact global helpers.  Any failure
stops the experiment before timing.  Only after those gates and independent
source review may the frozen one-GPU panel run.
