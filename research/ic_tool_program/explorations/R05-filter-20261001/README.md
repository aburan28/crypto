# R05 exploration: a sharper presence filter, before declaration

**What this is.** An exploratory timing of the presence filter that R05
then declared ([`../../rounds/R05-presence-filter/PROTOCOL.md`](../../rounds/R05-presence-filter/PROTOCOL.md)).
It ran on 2026-10-01, 17:06–17:22 UTC, between R03's timed processes,
under the benchmark lock. It tested whether the lever was worth a round.
It is not a round: it claims nothing, its arms are not R05's, and R05
neither reuses nor pools its runs.

## Why the filter

R04's check after its run (`../../rounds/R04-scan-probes/filter_check.json`)
found that from `2^36.6` up, the admitted fraction of scanned keys equals
the presence filter's false-positive rate to within 2%. The filter sets
one bit a stored pair, at 4–8 bits a pair. So the admitted stage, 20–42%
of the scan at the top sizes, is spent almost entirely on keys that are
not there.

## The prototype

[`prototype.patch`](prototype.patch) is commit `0318f016` against its
parent `eaa842aa`. `eaa842aa` is the tip of Track B's local stack (B3b,
work in progress), which is on record in Track B's bundle
([`stack-20261001.txt`](../../track-b/stack-20261001.txt) lists it).
The prototype sat there because the stack's two-word table builds the
filter too.

What it changes:
- **Three bits a key, in one 64-bit word.** `filter_slot` sets bit
  `h & 63` of word `(h & mask) / 64`, as before, and two more from `h`'s
  top twelve bits, which the word index never reads
  (`h = pair_filter_hash(key)`). A probe still reads one word.
- **Eight bits a stored key**, before rounding up to a power of two (it
  was four).
- **Everywhere the filter is built or read:** the five one-word builders,
  the probe, the two-word table, and the GPU fold kernel in
  `gpu/ecc2k/pairtable.cuh`. `FoldedParts` carries `filter_probes`.

R05's candidate is this change ported to the newest baseline, with the
GPU table file's version bumped. The protocol lists the differences.

## How it ran

[`explore.py`](explore.py) `run`, from a scratch copy whose harness is
`research/ic_tool_program/harness` as on `main`:
- **Rows:** `M1`'s two rows at three sizes: `2^39.0`, `2^44.3` and
  `2^47.2`.
- **Rounds:** three, ABAB, each process isolated by the harness.
- **Processes:** 37, of which one was contended and re-run, so 36 clean
  pairs (`summary.json`, `accounting`). None failed.

The binaries, in [`binaries.sha256`](binaries.sha256), were built in
release from `eaa842aa` (the base) and `0318f016` (the prototype). The
run tree is [`runs.tar.xz`](runs.tar.xz) (`runs.tar.xz.sha256`), and
`explore.py summarise` makes [`summary.json`](summary.json) from it.

## What it found

Base over prototype, so a ratio above 1 means the prototype is faster.
Six pairs a size; the interval is the 95% `t` interval of the geometric
mean.

| curve | `log₂ r` | cold time | collection | outputs |
|:--|--:|--:|--:|:--|
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 1.125 [0.976, 1.296] | 1.204 [1.034, 1.402] | identical |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.274 [1.211, 1.340] | 1.498 [1.440, 1.558] | identical |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.374 [1.315, 1.435] | 1.445 [1.377, 1.518] | identical |

- **The outputs are identical at every pair:** the counts, both arms'
  scalars, rho's counts and the verification flags.
- **Collection runs 1.45–1.50× faster at the two largest sizes,** and
  cold time 1.27–1.37×.
- **The spread is small where the lever is large.** The per-pair log
  standard deviation is 0.048 and 0.041 at `2^44.3` and `2^47.2`, against
  0.135 at `2^39.0`.

## What it does and does not establish

- **It set R05's predictions and its power estimate.** The rule R05
  applies is R02b's, unchanged.
- **It is not R05's comparison.** The base here is Track B's stack, not
  an accepted baseline, and the prototype lacks the candidate's file
  version and comment fixes.
- **`2^44.5` was not run.** R05's prediction there is the model's alone.
- **Its rows are among R05's suite rows** (`M1`). R05's fresh holdouts
  follow suite v1's construction at the next unused seeds, 210–213, and
  no process has run them.
