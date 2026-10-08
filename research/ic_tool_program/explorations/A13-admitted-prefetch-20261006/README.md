# A13: each admitted key's run prefetched, on the folding kernel

**What this is.** A timed exploration of A13, the admitted keys' cost in
the collection scan (`research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§8). It ran on 2026-10-06 from 08:13 to 08:21 UTC on the folding-kernel
candidate (R08's, `a4c19ed7`), before any round was declared on it. It is
a record, not a round's measurement. None of its runs is pooled with any
round's.

**Why.** On the folding-kernel candidate the scan's stages put the
admitted keys at 4.4–6.6 ns of a scanned summand's 21–22.5 ns, and each
admitted key at 210–225 ns
([record](../R08-fold-kernel-20261006/README.md)). That is the cost of one
random access into the largest table, about 221 ns on this host at that
size. Nearly every admitted key is a filter false positive, so nearly all
of that time is spent waiting.

**The outcome.** Collection is 1.12–1.29× faster, and the build and every
logarithm are unchanged. Cold time falls by 7–11%. **A round is to be
declared on it,** on the baseline that R08 decides, on a host of the
class it ran on.

## The candidate

`73799778`, one commit on `a4c19ed7`, frozen here as a plain diff,
[`admitted-prefetch.patch`](admitted-prefetch.patch) (38 lines):
- **Before:** the scan prefetched the presence filter ahead of its keys
  and, for an admitted key, the bucket entry. The run the entry points at
  was then read cold.
- **After:** between the filter pass and the lookups, a second pass
  (`prefetch_run`) reads each admitted key's bucket entry, cached by then,
  and prefetches the run's first and last words. The runs' misses overlap.
- **Nothing a prefetch reads changes a value:** the lookups, the matches
  and their order are what they were.

**The tests:** `cargo test --release --lib -- koblitz_fast
koblitz_index_calculus` on `73799778`'s tree: 173 passed, none failed, 3
ignored ([`tests.log`](tests.log), the build's log).

## What ran

- **The rows:** `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
  `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`.
- **The arms:** `sub` (`ic-r09sub-on-34ad98f9-a4c19ed7`, the folding-kernel
  candidate) and `admit` (`ic-r11admit-on-a4c19ed7-73799778`). Both are
  hashed in [`binaries.sha256`](binaries.sha256), and every run's report
  carries its arm's one BLAKE3.
- **The schedule:** four rounds, the order alternating, every process
  `ic price --single-target` with the suite's rho seed,
  `RAYON_NUM_THREADS=1`, isolated on CPU 2: 48 processes
  ([`explore.sh`](explore.sh)).
- **Three of the 48 were marked contended,** one `sub` and two `admit`,
  leaving six clean pairs at `2^44.3`, seven at `2^44.5` and eight at
  `2^47.2`.
- **Every pair recovered the same logarithm** in both arms, all 24
  ([`pairs.tsv`](pairs.tsv)).

## Results

Cold time, `sub` over `admit` (above 1, the prefetch is faster), 95% t
intervals of the geometric mean, on the clean pairs:

| curve | `log₂ r` | clean pairs | cold time | 95% interval | all pairs |
|:--|--:|--:|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 6 | **1.075** | [0.984, 1.174] | 1.060 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 7 | **1.125** | [1.018, 1.243] | 1.114 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 8 | **1.086** | [0.975, 1.209] | 1.086 |

The phases, as the geometric mean of the clean pairs' ratios, `sub` over
`admit`:

| curve | collection | lowest pair | build |
|:--|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.149 | 1.054 | 0.992 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.288 | 1.094 | 0.996 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.117 | 0.923 | 1.007 |

The median over each arm's clean processes, in ms
([`phases.tsv`](phases.tsv)):

| curve | arm | processes | cold | build | collect |
|:--|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | `sub` | 7 | 523.7 | 192.8 | 268.1 |
| | `admit` | 7 | 480.7 | 190.0 | 230.0 |
| `icv1-f2m59-tm943548413-98844ecc` | `sub` | 8 | 590.1 | 203.0 | 299.2 |
| | `admit` | 7 | 527.4 | 199.6 | 245.7 |
| `icv1-f2m61-t158598901-ab42b6c5` | `sub` | 8 | 2,059.9 | 384.2 | 1,567.2 |
| | `admit` | 8 | 1,839.9 | 370.4 | 1,347.6 |

## Reading

- **Collection falls by 38, 54 and 220 ms** at the medians, and cold time
  by 43, 63 and 220 ms. The build, which the commit does not touch, does
  not move.
- **About half to two-thirds of the admitted keys' waiting is gone.**
  This is an estimate from two sets of processes, not a measurement: the
  probes on the candidate priced the admitted stage at about 82, 90 and
  327 ms a cold run, and collection fell by 47%, 59% and 67% of that.
- **Cold time is noisy here.** The per-pair spread of the log ratio is
  0.08, 0.11 and 0.13, the largest of this day's explorations. One pair at
  `2^47.2` reads 0.92 in collection. Two rows in four rounds do not pin a
  7–11% gain; a round's forty pairs a set would, at a 95% half-width of
  about 4%.
- **The lever is portable,** unlike the folding kernel. A prefetch hint
  runs on any CPU, so the round that measures it can name the class it
  ran on. On a host without the vector kernels the scan's other stages are
  slower, so the share it saves is smaller there.

## Its limits

- **Two rows a size, in four rounds:** enough to declare, not to accept.
- **One host class,** this host's Xeon, with AVX-512F, GFNI, VBMI2 and
  VPCLMULQDQ. No other class was measured.
- **Wall time on this container** carries the ±5–10% residual noise that
  AGENTS.md §10 describes.

## Files

- [`runs.tar.xz`](runs.tar.xz) holds each process's report, isolation
  record and logs, under the suite's directory names, in `runs/`, and the
  exploration's log in `logs/`.
- **The tables come from it:**
  - [`pairs-tsv.sh`](pairs-tsv.sh) gives [`pairs.tsv`](pairs.tsv), as
    `pairs-tsv.sh runs sub admit`;
  - [`pairs2.sh`](pairs2.sh) gives the intervals;
  - [`phase-medians.sh`](phase-medians.sh) gives [`phases.tsv`](phases.tsv).
- [`SHA256SUMS`](SHA256SUMS) checks every file here.
