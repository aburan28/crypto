# Amendment 1 (2026-10-14, before any data for these arms)

Recorded while the frozen grid (`PROTOCOL.md`) was running; it adds arms
and changes none of the frozen ones.

## Why

The frozen engine comparison hands CADO `bwc` our *unmerged* core, because
our merge produces coefficients that are arbitrary residues mod p and the
`bwc` format carries 32-bit integers. That isolates the Krylov engines but
is not how CADO-NFS runs: it filters with its own `purge` and `merge-dl`,
whose merges keep small integer coefficients. The pipeline-against-pipeline
comparison needs CADO's filter too.

## Added arms (m = 3, D = 1,024 … 262,144 in powers of two, seeds 1–3,
four threads, same timeout and stop rule as the frozen grid)

| arm | what |
|:--|:--|
| `cado-full` | `purge -keep 160` → `merge-dl -target_excess 0` (default `target_density` 100) → `replay-dl` → native rhs assembly → `bwc.pl` m = n = 4, `thr=2x2` → native propagation reconstruction (CADO's `reconstructlog-dl` does this step for NFS; ours is the same iteration over the original rows without the NFS parts) |
| `full-bw-pos` | `full-bw` on the same matrices |
| `full-bl-pos` | `full-bl` on the same matrices |

These run on **`--coeffs positive`** matrices: every summand coefficient +1
(a repeated column gives +2). CADO's relation format stores an ideal once
per unit of its nonnegative valuation, so ±1 rows cannot be written. The
generator draws the same random stream as `unit`, so the column pattern is
identical to the frozen grid's matrices; only the signs differ.

Timed for `cado-full`: each CADO binary from start to exit, plus the native
rhs assembly and reconstruction (in-process timers). Not timed: the export of
the relation file (input preparation, as for our arms' generation).

Verification: all D logs reconstructed and equal to the planted vector.

## Contention

Two builds of the instrument and two 1-thread smoke runs overlapped the
frozen grid (windows in `raw/contention.log`, all during D ≤ 8,192 runs).
Runs whose interval overlaps a window are flagged in the report.
