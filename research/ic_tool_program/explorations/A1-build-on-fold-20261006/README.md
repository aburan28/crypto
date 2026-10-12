# A1 on the folding kernel: the build's partition streams and per-run filters

**What this is.** A timed exploration of the plan's A1, the pair-table
build (`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8), run on the
folding-kernel candidate (R08's, `a4c19ed7`). It ran on 2026-10-06 from
08:04 to 08:13 UTC, before any round was declared on it. It is a record,
not a round's measurement. None of its runs is pooled with any round's.

**Why again.** On v3 the same lever cut the build by 10–25%, but no
variant moved cold time reliably, because the build was 14–30% of cold
time there ([that record](../A1-build-partitions-20261006/README.md)).
On the folding-kernel candidate the build is a larger share of a shorter
cold time: about 36%, 34% and 17% at the three sizes below.

**The outcome.** The build phase is 1.40–1.47× faster, on every pair.
Collection is unchanged, and so are the table and every logarithm. Cold
time falls by 7–12%. **A round is to be declared on it,** on the
baseline that R08 decides.

## The candidate

Two commits on `a4c19ed7`, frozen here as plain diffs:

| commit | what changes | patch |
|:--|:--|:--|
| `54ccf235` | A1's v1 streams, ported: each pair's `(hash, word)` goes straight into one of 256 coarse partitions, appended to its run's stream, instead of a per-row list that is counted and scattered. The rows are cut into contiguous runs, one a thread, so the table is word for word the same. | [`a1-streams-on-fold.patch`](a1-streams-on-fold.patch) |
| `4c71a00e` | The filter is set a run at a time: each run sets a filter of its own with plain ORs, prefetching the word sixteen pairs ahead, and the runs' filters are ORed at the end. At one thread there is one run, and its filter is the filter. | [`a1-run-filters-on-fold.patch`](a1-run-filters-on-fold.patch) |

[`build-candidate.patch`](build-candidate.patch) is the two together, the
diff from `a4c19ed7` to `4c71a00e`.

**The port.** `54ccf235` is A1's v1 (`41d80f7b` on v3) carried onto
`a4c19ed7`, whose build's rows already take the folding kernel past their
own orbit. The two diffs differ only where they meet that change: a row
now hands its keys over a slice at a time, not whole.

**What replaced the locked OR.** On v3, v2's tagged streams set the
filter with an atomic OR as each pair was keyed. That waited out every
store before it. `4c71a00e` takes v2b's remedy: per-run filters, ORed
after the sort. It adds a prefetch, so each run's misses overlap.

**The tests.**
- `the_partitioned_fold_build_is_one_table_word_for_word_on_any_thread_count`
  compares the table and the filter word for word at one, two, three and
  five threads, and now covers `2^59` and `2^61` too.
- `cargo test --release --lib -- koblitz_fast koblitz_index_calculus` on
  `4c71a00e`'s tree: 174 passed, none failed, 3 ignored
  ([`tests.log`](tests.log)). They ran after the timed runs, under the
  isolation tool's lock.

## What ran

- **The rows:** `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
  `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`.
- **The arms:** `sub` (`ic-r09sub-on-34ad98f9-a4c19ed7`, the folding-kernel
  candidate) and `build` (`ic-r10build-on-a4c19ed7-4c71a00e`). Both are
  hashed in [`binaries.sha256`](binaries.sha256), and every run's report
  carries its arm's one BLAKE3.
- **The schedule:** four rounds, the order alternating, every process
  `ic price --single-target` with the suite's rho seed,
  `RAYON_NUM_THREADS=1`, isolated on CPU 2: 48 processes
  ([`explore.sh`](explore.sh)).
- **Four of the 48 were marked contended,** one `sub` and three `build`,
  leaving seven clean pairs at `2^44.3`, five at `2^44.5` and eight at
  `2^47.2`.
- **Every pair recovered the same logarithm** in both arms, all 24
  ([`pairs.tsv`](pairs.tsv)).

## Results

Cold time, `sub` over `build` (above 1, the build lever is faster), 95% t
intervals of the geometric mean, on the clean pairs:

| curve | `log₂ r` | clean pairs | cold time | 95% interval | all pairs |
|:--|--:|--:|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 7 | **1.119** | [1.011, 1.237] | 1.110 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 5 | **1.101** | [1.016, 1.192] | 1.106 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 8 | **1.071** | [0.986, 1.164] | 1.071 |

The phases, as the geometric mean of the clean pairs' ratios, `sub` over
`build`:

| curve | build phase | lowest pair | collection |
|:--|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.395 | 1.172 | 1.000 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.430 | 1.310 | 0.976 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.473 | 1.358 | 1.005 |

The median over each arm's clean processes, in ms
([`phases.tsv`](phases.tsv)):

| curve | arm | processes | cold | build | collect |
|:--|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | `sub` | 8 | 569.8 | 203.8 | 299.9 |
| | `build` | 7 | 519.7 | 148.2 | 306.6 |
| `icv1-f2m59-tm943548413-98844ecc` | `sub` | 7 | 588.5 | 201.0 | 295.6 |
| | `build` | 6 | 536.4 | 143.2 | 313.3 |
| `icv1-f2m61-t158598901-ab42b6c5` | `sub` | 8 | 2,106.5 | 362.8 | 1,630.1 |
| | `build` | 8 | 1,937.3 | 252.1 | 1,569.7 |

## Reading

- **The build gain is steady, and cold time follows it.** The build falls
  by 56, 58 and 111 ms, and the cold medians by 50, 52 and 169 ms.
  Collection, which reads the same table, does not move.
- **Cold time is noisy here.** The per-pair spread of the log ratio is
  0.11, 0.06 and 0.10, against 0.02–0.09 in the folding kernel's
  exploration. One pair at `2^44.3` reads 0.90 and one reads 1.30. Two
  rows in four rounds cannot pin a 7–12% gain. A round with forty pairs
  a set would: at a spread of 0.11 its 95% half-width is about 3.5%.
- **The per-run filters matter more here than on v3.** v1's streams alone
  cut v3's build by 10–18%. Here the two commits together cut it by
  28–32% at the medians. No arm ran v1 alone on the folding kernel, so
  the split between the two commits is not measured.

## Its limits

- **Two rows a size, in four rounds:** enough to declare, not to accept.
- **One host class,** this host's Xeon. The build's code has no
  CPU-specific path, but no other class was measured.
- **Wall time on this container** carries the ±5–10% residual noise that
  AGENTS.md §10 describes.

## Files

- [`runs.tar.xz`](runs.tar.xz) holds each process's report, isolation
  record and logs, under the suite's directory names, in `runs/`, and the
  exploration's log in `logs/`.
- **The tables come from it:**
  - [`pairs-tsv.sh`](pairs-tsv.sh) gives [`pairs.tsv`](pairs.tsv), as
    `pairs-tsv.sh runs sub build`;
  - [`pairs2.sh`](pairs2.sh) gives the intervals;
  - [`phase-medians.sh`](phase-medians.sh) gives [`phases.tsv`](phases.tsv).
- [`SHA256SUMS`](SHA256SUMS) checks every file here.
