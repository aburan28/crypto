# A1 explorations: the folded build keyed straight into its partitions, on v3

**What this is.** Timed explorations of the plan's A1, the pair-table
build (`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8). They ran on
2026-10-06 on v3, main's head `995ea207`, before any round was declared
on them. They are a record, not a round's measurement, and none of their
runs is pooled with any round's.

**The outcome.** The build phase falls by 10–25%, but no variant moves
cold time reliably on v3, so **no round is declared on these
candidates.** The lever continues on the folding-kernel candidate (the
plan's §8), where the build's rows take the same kernel. That is a
separate exploration ([record](../A1-build-on-fold-20261006/README.md)).

## Where the build's time goes (the stage probes)

`build-probes.patch` adds clock reads at the folded build's stages (a
build of its own, `ic-a1-probes-on-995ea207`). It ran on `M1`'s two rows
at the three largest sizes, two rounds, isolated
([`probe.sh`](probe.sh)). At `icv1-f2m61-t158598901-ab42b6c5`, with 6.29 M
stored pairs, the fill took 410–450 ms:
- **Keying the rows, 200–230 ms:** the additions 58–69 ms, the canonical
  keys 68–82 ms, and the per-row lists' growth (`emit`) 70–80 ms.
- **Counting and scattering the lists into partitions,** 14–15 and
  73–76 ms.
- **Sorting each partition and filling the table and filter,** 103–118 ms.

So the lists' growth, count and scatter were about 26 ns of the 61–68 ns
a stored pair costs at `n = 53`–`61`.

## The candidates

Three commits, each on the one before, frozen here as plain diffs:

| candidate | commit | what changes | patch |
|:--|:--|:--|:--|
| v1 | `41d80f7b` | each pair's `(hash, word)` goes straight into one of 256 coarse partitions, appended to its run's stream, instead of a per-row list that is counted and scattered; the rows are cut into contiguous runs, one a thread, so the table is word for word the same | [`a1-v1-streams.patch`](a1-v1-streams.patch) |
| v2 | `6563a33b` | a tagged build streams four bytes a pair (the bucket's place in its partition and sixteen bits of the hash); the orbit comes back from per-row counts; the filter is set as each pair is keyed | [`a1-v2-tagged.patch`](a1-v2-tagged.patch) |
| v2b | `76a8de3e` | each run sets the bits in a filter of its own with plain ORs, and the runs' filters are ORed together after the sort: the same filter, word for word | [`a1-v2b-run-filters.patch`](a1-v2b-run-filters.patch) |

**Every candidate's table is the base's, word for word** (tested,
including a forced spill, at one and three threads). Every pair below
recovered the same logarithm in both arms.

## What ran

`M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
`icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`,
three rounds, every process `ic price --single-target`, isolated on CPU 2,
`RAYON_NUM_THREADS=1`:
1. **v1 against v3** ([`explore.sh`](explore.sh)), the order alternating;
2. **v2 against v3** ([`explore2.sh`](explore2.sh));
3. **v3, v1 and v2b together** ([`explore3.sh`](explore3.sh)), the order
   rotating by round.

The binaries are hashed in [`binaries.sha256`](binaries.sha256), and every
pair is in the `pairs-*.tsv` files.

## Results

Cold time, v3 over the candidate (above 1, the candidate is faster), 95%
t intervals of the geometric mean, on the clean pairs:

| exploration | arm | `2^44.3` | `2^44.5` | `2^47.2` |
|:--|:--|:--|:--|:--|
| 1 | v1 | 1.049 [0.975, 1.128] (4) | 1.036 [0.967, 1.110] | 1.020 [0.987, 1.054] |
| 2 | v2 | 0.974 [0.847, 1.119] | 1.012 [0.897, 1.142] | 0.977 [0.949, 1.005] |
| 3 | v1 | 1.010 [0.962, 1.061] | 1.038 [0.974, 1.105] (5) | 1.016 [0.956, 1.079] |
| 3 | v2b | 1.092 [1.026, 1.163] | 1.135 [1.007, 1.279] | 0.963 [0.898, 1.034] |

Six pairs a size unless marked. The build phase, the median over each
arm's processes, in ms:

| exploration | arm | `2^44.3` | `2^44.5` | `2^47.2` |
|:--|:--|--:|--:|--:|
| 1 | v3 | 238.8 | 254.3 | 444.6 |
| 1 | v1 | 213.2 | 227.8 | 364.8 |
| 2 | v3 | 230.8 | 263.1 | 437.2 |
| 2 | v2 | 225.3 | 245.1 | 451.2 |
| 3 | v3 | 267.6 | 283.6 | 468.5 |
| 3 | v1 | 222.3 | 239.3 | 396.9 |
| 3 | v2b | 203.8 | 213.3 | 429.0 |

## Reading

- **v1 is the steady gain.** Its build is 10–18% faster in both
  explorations that ran it, at every size. That is 25–80 ms, a few per
  cent of cold time, inside the noise of six pairs.
- **v2 alone is no faster.** The tagged streams save bytes, but the
  build does not fall. v2b's commit names the likely reason: a locked
  filter OR in the keying loop waits for every store before it, the
  streams' included.
- **v2b removes that wait, and gains at two sizes.** Its build is the
  fastest at `2^44.3` and `2^44.5`, 22–25% below v3. At `2^47.2` it is
  8% below v3 and slower than v1. No probe split v2b's build, so why is
  not established here.
- **None of them is a round on v3.** The build is 14–30% of cold time
  there, so even the best build gain is a few per cent of cold time, and
  three rounds cannot resolve it. On the folding-kernel candidate the
  build's additions are cheaper and its share of cold time is larger.
  The lever is explored again there, with v1's streams and per-run
  filters that are prefetched, not locked.

## Files

- [`runs.tar.xz`](runs.tar.xz) holds the run trees, each process's report,
  isolation record and logs under the suite's directory names, and their
  logs:
  - `runs/` holds exploration 1 and the probes (`runs/probe`);
  - `runs2/` holds exploration 2;
  - `runs3/` holds exploration 3.
- **The tables come from it:**
  - [`pairs-tsv.sh`](pairs-tsv.sh) gives the `pairs-*.tsv` files, as
    `pairs-tsv.sh runs base cand`;
  - [`pairs2.sh`](pairs2.sh) gives the intervals;
  - [`phase-medians.sh`](phase-medians.sh) gives the phases.
- [`SHA256SUMS`](SHA256SUMS) checks every file here.
