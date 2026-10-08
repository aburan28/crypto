# Resolver geometry: static cost ledger

This ledger prices only the scheduling quantity changed by the preregistered
patch.  It is not a throughput result.  `analyze.cpp` regenerates every row and
fails if the production population, shared-table size, queue size, SM count or
observed 124-register resolver baseline drifts.

| threads | admitted | grid lanes | warps/SM | registers/block | register headroom | table words/lane | table bytes copied/grid | queue bytes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 64 | no | 12,032 | 2 | 7,936 | 57,600 | 223 | 10,725,776 | 24,641,536 |
| 128 | yes, default | 24,064 | 4 | 15,872 | 49,664 | 112 | 10,725,776 | 24,641,536 |
| 256 | yes | 48,128 | 8 | 31,744 | 33,792 | 56 | 10,725,776 | 24,641,536 |
| 512 | yes | 96,256 | 16 | 63,488 | 2,048 | 28 | 10,725,776 | 24,641,536 |

All rows keep 188 resolver blocks, 57,052 dynamic shared bytes per block,
6,160,384 queue entries, 24,641,536 queue bytes and one four-byte counter.
The table-copy byte total is identical; larger blocks distribute it among more
lanes.  Queue owners, cold calls and hot work are identical.

The 64-thread form is a static no-go.  The table already limits the resolver to
one block per SM, so 64 threads expose only two warps while retaining every
table copy, queue byte and launch.  The admitted forms test the only positive
mechanism: increasing resident resolver warps from four to eight or sixteen.

The 512-thread row has only 2,048 simple register slots of headroom at the
observed 124 registers/thread.  This ledger establishes feasibility, not ptxas
allocation; native-sm120 compilation and device resources remain required
before any GPU performance screen.
