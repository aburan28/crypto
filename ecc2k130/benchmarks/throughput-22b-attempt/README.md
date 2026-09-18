# RTX PRO 6000 22 B/s attempt

Canonical report: [`../../THROUGHPUT-22B.md`](../../THROUGHPUT-22B.md).

The target was not met. Six equal-work alternating pairs confirm
**20.134326 B/s** against **18.164748 B/s**, a 1.1084 ratio. Both binaries
replayed 300/300 reports with zero dropped.

Files:

- `result.json`: frozen configuration, table rows, paired statistics and
  profile counters.
- `balanced-comparison.log`: six alternating equal-update pairs.
- `build-{baseline,candidate}.log`: CUDA 13.3.73 compiler output.
- `verify-{baseline,candidate}.log`: independent scalar replay.
- `profile-candidate.log`: Nsight Compute pipeline, scheduler and memory
  counters.
- `scout-*.log` and `rejected-*.log`: retained intermediate and negative rows.

The candidate differs from the 18 B/s composed reference by:

```text
BATCH=17 THREADS=512 MINBLOCKS=1 UNROLL_SLOTS=1
TABLE_RECOMPUTE_DENOM=1 TABLE_BANK_PAD=1
```

Both arms otherwise use the inline-both table walk, ALU polynomial square,
paired CLMAD, reduced-input conversion and L2 persistence settings frozen in
`result.json`.
