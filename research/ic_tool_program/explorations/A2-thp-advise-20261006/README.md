# A2 exploration: huge pages advised on the pair tables, on v3

**What this is.** A timed exploration of the plan's A2, huge pages for the
scan's tables (`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8). It
ran on 2026-10-06 from 05:57 to 06:05 UTC, on v3, before any round was
declared on it. It is a record, not a round's measurement. Its outcome is
that **no round is declared**: the advice takes effect, and cold time does
not move.

## The candidate

`b2e9d484`, one commit on v3 (`995ea207`), frozen here as
[`thp-advise.patch`](thp-advise.patch):
- **Zeroed allocation and advice:** each pair-table builder allocates its
  tables zeroed (the presence filter, the bucket starts and the stored
  words) and advises `MADV_HUGEPAGE` on their 2 MiB-aligned interior
  before the first write.
- **The filter is filled in place** through an atomic view, which also
  drops the copy out of the atomic filter.
- **The tables are unchanged,** bit for bit.

**Why it was worth a run.** The earlier test, through glibc's tunable,
reads 0.96–1.00× in cold time at the three largest sizes (the
scoreboard's "Huge pages are not a lever"). That test advised every
allocation; this one advises only the tables the scan reads at random.

**This host honours the advice.** `transparent_hugepage/enabled` reads
`[madvise]`. The candidate's processes took far fewer minor page faults,
the mean over each arm's six processes (from each run's isolation
record):

| curve | base | candidate |
|:--|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 56,215 | 39,902 |
| `icv1-f2m59-tm943548413-98844ecc` | 87,907 | 52,496 |
| `icv1-f2m61-t158598901-ab42b6c5` | 127,212 | 69,010 |

## What ran

- **The rows:** `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
  `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`.
- **The schedule:** three rounds, the order alternating (base first,
  candidate first, base first), every process `ic price --single-target`,
  isolated on CPU 2, `RAYON_NUM_THREADS=1`: 36 processes
  ([`explore.sh`](explore.sh)).
- **The binaries** are hashed in [`binaries.sha256`](binaries.sha256).
- **Three runs were marked contended,** all at `2^44.3`. That leaves four
  clean pairs there and six at each other size.
- **Every pair recovered the same logarithm** in both arms
  ([`pairs.tsv`](pairs.tsv)).

## Results

Cold time, base over candidate (above 1, the candidate is faster), 95% t
intervals of the geometric mean:

| curve | `log₂ r` | clean pairs | cold time | 95% interval |
|:--|--:|--:|--:|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 4 | 0.960 | [0.798, 1.156] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 6 | 1.068 | [0.979, 1.165] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 6 | 1.007 | [0.950, 1.067] |

**No size gains.** Every interval contains 1. At `2^47.2`, the size with
the largest tables and the most faults saved, the point estimate is
1.007.

## Reading

- **Faults were not the scan's cost.** Halving the minor faults saves
  their kernel time, about 60k faults, a few tens of milliseconds at
  most, inside a cold time of 3.4 s.
- **The scan's misses are cache misses, not TLB misses.** If the TLB
  misses had mattered, huge pages would have shown in collection. R04's
  probes on v3 put the filter and the admitted keys at 10–12 ns a summand
  ([record](../scan-stages-v3-20261006/README.md)): the cost of reaching
  DRAM at all, which a larger page does not change.
- **A2 is set aside.** The patch stays on record. No round is declared
  on it.

## Files

- [`runs.tar.xz`](runs.tar.xz) holds each process's report, isolation
  record and logs, under the suite's directory names, in `runs2/`.
- **The tables come from it:**
  - [`pairs-tsv.sh`](pairs-tsv.sh) gives [`pairs.tsv`](pairs.tsv), as
    `pairs-tsv.sh runs2 base cand`;
  - [`pairs2.sh`](pairs2.sh) gives the intervals.
- [`SHA256SUMS`](SHA256SUMS) checks every file here.
