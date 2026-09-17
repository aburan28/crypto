# G7 larger-batch X-cache experiment

Status: complete; selected B16/cache4/min3 retained. See [RESULTS.md](RESULTS.md).

## Objective

Test whether retaining half of each split worker's X coordinates in shared
memory recovers the coordinate traffic lost by the B24/B32 experiment while
preserving its larger block-inversion batch.  The selected B16 kernel caches
four of eight local X slots.  The matched candidates therefore cache six of
twelve local slots for B24 and eight of sixteen local slots for B32.

## Frozen controls and candidates

All variants use 256 threads, split2, compact polynomial state, block
inversion, weighted prefix 2, four-byte sigma tables, and the CUDA 13.3
front-end / CUDA 13.4 assembler toolchain.  The experiment changes only batch,
the number of shared X slots, and the minimum-block launch bound.

| Label | Batch | Shared X slots per thread | Minimum blocks/SM | Role |
|---|---:|---:|---:|---|
| selected | 16 | 4 | 3 | production reference |
| batch24-cache4-min3 | 24 | 4 | 3 | prior larger-batch reference |
| batch32-cache4-min3 | 32 | 4 | 3 | prior larger-batch reference |
| batch24-cache4-min2 | 24 | 4 | 2 | launch-bound control |
| batch24-cache6-min2 | 24 | 6 | 2 | B24 candidate |
| batch32-cache4-min2 | 32 | 4 | 2 | launch-bound control |
| batch32-cache8-min2 | 32 | 8 | 2 | B32 candidate |

The candidate source broadens `ECC_PACKED_SHARED_X_SLOTS` from the historical
set 0/2/4 to the closed range 0..8.  Existing runtime layout checks remain.
No production source is modified.

## Gates

1. Build all four new variants and record binary hashes and ptxas resources.
2. Require byte-identical built-in test output.
3. Compare normalized full and ragged CUDA client state, DP multisets and CPU
   replay against the selected B16 implementation for each candidate.
4. Require actual-walk shared-X unit checks and memcheck, initcheck and
   synccheck with zero errors for timing-eligible candidates.
5. Time every eligible arm in three interleaved equal-work repetitions on the
   same local g7.2xlarge / RTX PRO 4500 at 165 W.  Each sample performs exactly
   34,359,607,296 complete scalar updates.
6. A candidate qualifies only when its paired geometric speedup over selected
   has a two-sided 95% confidence interval wholly above 1.  Promotion also
   requires a fresh held-out confirmation.

Rates are billions of complete scalar rho updates per second.  Generic work is
reported at `sqrt(n/262)` with ratio 1; full-DLP S remains null.

## Result

B32/cache8 produces a reproducible local gain over B32/cache4: 1.75% with the
min2 schedule and 1.47% with the min3 schedule. Both larger-batch candidates
remain slower than selected B16. The best follow-up candidate reaches 6.050341
B/s versus 6.208908 B/s selected, with paired speedup 0.963305 and 95% CI
[0.935844, 0.991572]. No candidate qualifies for promotion and the 12 B/s goal
remains unmet.
