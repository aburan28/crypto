# GPU-wide hint resolver block geometry: frozen preregistration

Source parent: `7e5281e6b276b89428dbb74a1c05a50ca5b61a02`.

Status: **implementation and safe-compile only.**  The active five-pair
GPU-wide-hint comparison must reach a terminal audited result before this
follow-up may allocate a GPU.  A fresh matched screen also requires explicit
root authorization.  This PR does not select a resolver geometry.

Class: compatible scheduling engineering.  The exact-v3 map, queue ownership,
cold cycle-anchor code, select grid, hot update kernel, shared table, worker
population, queue capacity, checkpoint state and all arithmetic remain fixed.
Only the resolver kernel's compile-time block width changes.

## Motivation and fixed costs

The production preflight at the source parent established identical
1,710,384-record corpora and passing full, partial and bidirectional checkpoint
gates.  Its resource tuple is:

| stage | registers/thread | local bytes/thread | shared bytes/block | launch geometry |
|---|---:|---:|---:|---|
| hot walk | 126 | 0 | 0 static | existing B16/T512/min1 grid |
| select | 56 | 0 | 57,052 dynamic | existing worker grid |
| resolver | 124 | 400 | 57,052 dynamic | 188 blocks x 128 threads |

The resolver table allocation admits only one block per SM, so the current
128-thread geometry exposes four resident warps per SM.  The same kernel at
256 or 512 threads would expose eight or sixteen warps while keeping one block
per SM.  At 124 registers/thread, the simple register footprints are 15,872,
31,744 and 63,488 registers/block respectively; the 512-thread form remains
below a 65,536-register SM budget but must be confirmed by ptxas.

At the frozen timing population, the queue remains exactly 6,160,384 entries
or 24,641,536 bytes plus one four-byte counter.  The resolver grid remains 188
blocks.  Every block still loads the same 57,052-byte table and every queued
owner still executes the same `tableResolveHintSlot` once.  The block widths
change only how those owners and table-copy words are distributed among lanes:

| resolver threads | grid lanes | warps/SM | approximate table words/lane |
|---:|---:|---:|---:|
| 128 (reference/default) | 24,064 | 4 | 112 |
| 256 | 48,128 | 8 | 56 |
| 512 | 96,256 | 16 | 28 |

`64` is deliberately not admitted: the shared-table limit would still allow
only one block per SM, leaving two resident warps while preserving all table
copies and launch overhead.  It has no occupancy mechanism that can improve on
128 and would expand the later panel without a positive static hypothesis.

Owner decoding (`owner % workers`, `owner / workers`) is a separate possible
hypothesis.  It must not be combined with this geometry experiment.

## Patch contract

`TABLE_GLOBAL_HINT_THREADS` is a compile-time Make knob with default `128`.
Only `128`, `256` and `512` are accepted.  It controls:

1. `resolveGlobalHints` launch bounds;
2. the resolver launch block width;
3. runtime identity and occupancy markers; and
4. native/device queue-control geometry.

It does not change `ECC_THREADS`, the 188-block resolver grid, queue memory,
select or hot-kernel geometry, dynamic shared bytes, map logic or checkpoint
format.  Default builds must remain byte-for-byte source-compatible in behavior
and continue to print 128 resolver threads.

## Required pre-GPU gates

Before any timing authorization:

- host ownership/phase controls pass for all admitted block widths, empty/full/
  sparse queues, partial worker blocks, production population and the 32-bit
  stride boundary;
- CUDA device queue controls compile separately for 128, 256 and 512 threads;
- default-128 CUDA controls retain their existing execution target;
- each client geometry safely compiles for native `sm_120` with the exact
  resolver launch bounds and 57,052 dynamic shared bytes;
- ptxas/cuobjdump receipts bind resolver registers, local bytes, static/dynamic
  shared bytes and maximum active blocks/SM;
- runtime markers bind queue entries/bytes, counter bytes, grid blocks, block
  threads, warps/SM and launch bounds; and
- no test changes owner decode, queue order requirements, table contents,
  selected/resolved history, report behavior or hot update work.

This draft PR remains unmerged while CUDA device gates are pending.  It may
merge as a preregistered plan/patch only after all applicable compile and
device-control checks pass; merging it still does not authorize a performance
run.

## Later matched screen, if authorized

The later screen holds the complete GPU-wide candidate fixed and compares
resolver widths 128, 256 and 512 on one RTX PRO 6000, CUDA 13.3.73 and native
`sm_120`.  It retains the existing odd-launch replay, spread replay, full and
partial corpus identity, bidirectional checkpoints, exact worker/update budget,
all cold costs, five A/A pairs and five alternating matched pairs.

A non-default geometry qualifies only if its median complete-update ratio to
128 is at least 1.05, every paired ratio is above one, A/A maximum symmetric
drift is below 1%, and every correctness/resource gate passes.  The original
GPU-wide candidate/control engineering threshold of 1.10 remains unchanged;
an isolated resolver improvement does not select the global map or establish
progress toward 26 B/s.  Otherwise retain 128 and preserve the negative result.

No search, solver, collision recovery or key recovery is in scope.
