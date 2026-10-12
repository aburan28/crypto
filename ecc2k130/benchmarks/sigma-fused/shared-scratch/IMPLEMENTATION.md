# Fused-sigma shared scratch: implementation status

Status: **draft implementation; compile gate passed, actual-device gates
pending.**  This branch
starts from `main` merge `bc217d318dde444014cde4e81dea02c70b126995`, which
contains the frozen [proposal](PROPOSAL.md) and [protocol](PROTOCOL.md).  No GPU
has been allocated for this implementation.

## Implemented scope

- `SIGMA_FUSED_SHARED_SLOTS` defaults to zero and accepts only `0`, `2`, `3`
  or `4`.
- A nonzero value is compile-time restricted to the exact counter-free
  B16/T256/minBlocks2 compact, shared-sigma, inline-product fused path.
- The first `C` `pchain` and tagged-denominator slots use the registered
  full-word `[field][slot][word][thread]` shared SoA.  Every other slot retains
  the existing compact global path.
- Producer and consumer call the same typed routing helpers.  The scratch adds
  no barrier and uses `threadIdx.x` only for the same thread's shared lane.
- Global `pchain` and denominator allocations, `x/y` storage, arithmetic,
  alternating traversal, checkpoint version and checkpoint payload are
  unchanged.
- Runtime output records the selected shared-slot count.

The implementation deliberately does not include compact shared packing,
`x/y` caching, global-allocation shrinking, a square-table combination or a
different block geometry.

## Native and device controls

`test_native.cpp` independently constructs compact fields, copies canonical
prefixes and six-bit tagged denominators through the full-word scratch helpers,
and compares the exact shared image and compact round trip.  It covers cache
sizes 2/3/4, worker counts 1/3/31/32/255/256, inactive lanes and both leading
and trailing canaries.

`testsigmasharedscratchcuda.cu` uses the production routing helpers.  For each
compiled cache arm it covers worker counts 1/3/255/256/257/511/512/513 plus an
entirely inactive block, verifies cached and global fallback slots, preserves
all canaries and six-bit tags, and queries the production walk's registers,
local bytes, static shared bytes and active blocks/SM.  It requires exactly two
active blocks before returning success.

The compile-only producer builds four production clients and three device
helpers with CUDA 13.3.73 for native `sm_120`.  Its independent C++ audit
requires the exact baseline flags, exact static shared ledger, no more than 128
registers, zero stack/local/spills and a static two-block capacity bound.  The
hosted workflow retains every build log, binary, resource dump and audit result.
Compile feasibility is not an occupancy measurement; the device test remains
mandatory on the actual RTX PRO 6000.

## Current evidence and stop

Local native results:

```
PASS: shared SoA ownership/banks, launch boundaries, traffic and capacity stop
PASS: 144 native compact/shared cases, 83232 field records,
      36414 tagged denominators; canaries and partial blocks intact
```

Strict compilation, UBSAN, the existing fused-star native suite and default
CPU compilation pass.  Exact-head hosted run
[`37053839217`](https://github.com/aburan28/crypto/actions/runs/37053839217)
compiled source `fc3941078db32d1a81f2732fd3622c13866208e6` with CUDA
13.3.73 for native `sm_120`.  All four arms report 126 registers/thread, zero
walk stack/local/spills, the registered ptxas static shared bytes, an exact
1,024-byte compiled reserve, and the registered compiled extents.  The static
capacity calculation remains two blocks/SM for every arm.  The retained
receipt is [`results/compile-fc394107/`](results/compile-fc394107/).

An independent read-only source and producer review passes the same-thread
routing, full tag, array bounds, partial-block barrier, launch boundary,
checkpoint and exact-work contracts at that source.  GPU helper execution,
runtime occupancy, replay/corpus/checkpoint equivalence and timing remain
unmeasured.  The branch stops before GPU dispatch until the root coordinator
releases the single frozen allocation.
