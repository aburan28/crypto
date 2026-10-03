# R05: a sharper presence filter

**Declared 2026-10-01, before any R05 timed run.** This is Track A, the
collection scan (A2) in `research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§8, at its admitted stage. An exploration ran first and is disclosed
below. Nothing below changes after the first timed run, except by a
dated amendment appended at the end.

## Hypothesis

**At the top sizes, a fifth to two fifths of the scan is spent on keys
the presence filter admits wrongly.**
- R04 priced the scan's stages ([`../R04-scan-probes/README.md`](../R04-scan-probes/README.md)).
  From `2^36.6` up, the admitted stage is 20–42% of the scan.
- R04's check after the run (`filter_check.json`) found the admitted
  fraction equal, within 2%, to the false-positive rate of the filter as
  built: one bit a key under one hash, at 4–8 bits a stored pair. True
  hits are one summand in 7,200 to 221,000.
- So the admitted stage is spent almost entirely on false positives, at
  118–159 ns an admitted key at the three largest sizes.

**A filter with three bits a key in one word, at eight to sixteen bits a
stored pair, passes far fewer absent keys at the same probe cost.**
- A probe still reads one 64-bit word, the same word as now. The filter
  stage's memory traffic stays at one access a key, in a filter twice as
  large.
- At the suite's stored-pair counts, the blocked-filter formula
  ([`predict.py`](predict.py)) gives these false-positive rates. The
  formula with one bit a key reproduces R04's measured rate.

  | curve | `log₂ r` | stored pairs | filter now → candidate | false positives now → candidate |
  |:--|--:|--:|:--|:--|
  | `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 3.71 M | 2 → 4 MiB | 0.198 → 0.028 |
  | `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 3.88 M | 2 → 4 MiB | 0.207 → 0.031 |
  | `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 6.27 M | 4 → 8 MiB | 0.170 → 0.019 |

## The phase and its share

- **Collection** is 66%, 72% and 88% of cold time at the three target
  sizes, in v0's profile (R01's medians, as `prediction.json` reads
  them).
- **The admitted stage** is 42%, 33% and 32% of the scan there (R04).

## Class

Engineering, if it gains. The answer is unchanged by construction, and
on the suite so is every count.

## The candidate

**Frozen here as [`candidate.patch`](candidate.patch)** (SHA-256 in
`candidate.patch.sha256`), against v0′ `c1a2e5f8`. It is the
exploration's prototype ported to the baseline. It adds the GPU table
file's version, a test, and documentation. It also corrects the GPU host
test's own copy of the byte estimate, which the prototype missed and the
GPU host tests caught.
- **The filter.** `filter_slot` sets three bits of word
  `(h & mask) / 64`, where `h = pair_filter_hash(key)`:
  - bit `h & 63`, as now;
  - two from `h`'s top twelve bits, which the word index never reaches.

  The filter is sized at eight bits a stored key before rounding up to a
  power of two (it was four), within `2^6`–`2^32` bits.
- **Everywhere the filter is built or read:**
  - the five builders in `koblitz_index_calculus.rs`, and the probe;
  - `gpu/ecc2k/pairtable.cuh`'s fold kernel, so a GPU-built table still
    matches the CPU's word for word.
- **External tables keep their meaning.** `FoldedParts` carries
  `filter_probes` (1 or 3), and the probe reads as many bits as were set.
  `gpu/ecc2k/fold_io.hpp` now writes `PTFOLD2`, and
  `examples/load_fold_table.rs` reads `PTFOLD1` as one bit a key.
  `from_folded_parts` refuses any other count.
- **The answer is unchanged by construction.** Neither filter has false
  negatives, and no pinned count reads the filter.
- **Three things change besides time:**
  - the filter's memory doubles, as in the table above;
  - the tables' byte estimates count the filter at eight bits a stored
    pair, not four. Under a memory budget an estimate chooses the tier,
    so a budget within half a byte a pair of an estimate can now get a
    narrower tier, or a refusal. Every suite size uses the folded tier,
    the narrowest, tens of megabytes inside the default 4 GiB, so no
    suite choice moves;
  - `perfbench`'s pair-table fingerprint, which hashes the filter's
    words, changes by construction. Nothing pins it.

The two-word table on Track B's local stack (B3b, not on `main`) takes
the same change when the stack is rebased onto R05's decision, by a
dated B3b amendment.

## Arms

- **The base**: the newest accepted baseline when R05 runs, as a commit
  recorded in the manifest. That is R03's candidate `30f6c153` if R03 is
  accepted, and v0′ `c1a2e5f8` otherwise.
- **The candidate**: the base plus `candidate.patch`, applied with
  `git apply` and committed with nothing else. It applies to both
  possible bases. R03 touches only `factorise_u64`, far from the
  filter's code.
- **The manifest** (`run.py manifest`) records both binaries' SHA-256
  and the commits they were built from.

## Before any timed step

1. **The tests,** on the candidate's tree:
   - `cargo test --release --lib -- koblitz_index_calculus`, which
     includes the patch's two tests (the filter's layout, and a table
     read as one bit a key) and `from_folded_parts`' refusals;
   - `make -C gpu/ecc2k test`. It runs the emulated fold kernel against
     the CPU's own table, loads the result back through `PTFOLD2`, and
     ends under ThreadSanitizer.
2. **The pin, untimed:** the candidate on all 90 suite rows. Every output
   must equal v0's from R01's profile pass: the counts, both arms'
   scalars, rho's counts and the verification flags. The names must
   agree, as in R02's pin.

## Rows

- **The suite rows (40 pairs at each target size):**
  - `M1`'s 22 rows: two targets at each of the eleven sizes;
  - at the three target sizes, the other six suite rows each (`M2`–`M4`).

  That is 40 rows. Five rounds, ABAB, isolated: 400 processes.
- **The fresh holdouts (40 pairs at each target size):** at the three
  target sizes only, eight rows each, by suite v1's own construction
  (`make_suite.py`), as R02b's were built:
  - recipe seeds 210, 211, 212 and 213;
  - two targets per seed, `T111` to `T118`;
  - `public_hash_seed` 23111 to 23118, and rho seeds `0x230000 + 111` to
    `+ 118`.

  They are committed in [`holdouts/`](holdouts/) (`SHA256SUMS`). Neither
  the seeds nor the targets have been run by any round or by the
  exploration. Five rounds, ABAB, isolated: 240 processes.
- **The extension.** At a target size, a set whose interval half-width
  (`hi / geomean − 1`) exceeds 3% after five rounds gets rounds 6–10, and
  its figure pools all ten. This depends only on the interval's width,
  never on its position.

## Power

The exploration's per-pair spread gives a 95% half-width of about 1.5%
and 1.3% at forty pairs at `2^44.3` and `2^47.2`. `2^44.5` was not
explored; R02's spread there gives 2.6–3.4%.

## Prediction

| curve | `log₂ r` | model | exploration, 6 pairs | predicted, both sets |
|:--|--:|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.31 | 1.274 [1.211, 1.340] | **1.27×** |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.25 | not run | **1.25×** |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.34 | 1.374 [1.315, 1.435] | **1.37×** |

- **The model** ([`prediction.json`](prediction.json)) scales R04's
  admitted stage by the predicted admitted keys a summand and keeps
  every other stage. It then weights collection by its share of v0's
  cold time. It leaves out the larger filter's own cost.
- **The other eight sizes: 1.02–1.13×** by the model, all above 1. None
  is predicted to regress.

## Success and stop

**Accepted** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At each of the three target sizes, the paired cold-time ratio's 95%
   interval lies above 1.10, on the suite rows and on the fresh holdouts
   separately.
3. No size regresses beyond its A/A band. A regression is a size whose
   interval's upper end lies below the lower end of R01's A/A interval,
   or the round's own if the host has changed.

**Rejected** otherwise.

**Stopped** if a test fails, an output differs, or a verification fails.

**If accepted**, the candidate becomes the next baseline, with its
binary hash, host, and `S` from this round's runs. It also gets the
rule's comparison: §23's protocol at its six sizes, with 64 targets, on
the candidate (384 processes).

## Inadmissible

- Reusing the exploration's runs, pooling them with R05's, or choosing
  between the two.
- Changing the filter's layout or sizing, the suite, the recipes or the
  targets, or drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting collection, the scan or the admitted stage as the speedup: the
  speedup is cold time.

## Cost

- The tests and the pin: a few minutes, and 90 untimed processes.
- The suite rows: 400 timed processes. The holdouts: 240, and up to 240
  more if they are extended. About three hours together.
- If accepted, the rule's comparison: 384 processes, about 40 minutes.

## Order, and plan §11

- **R05 runs after R03's decision and before R02b.** R02b's amendment 1,
  in the same pull request, moves R02b after R05's decision.
- **Why this order.** Both rounds are on the collection scan, after R02's
  rejection there. Plan §11 sets a phase aside after two consecutive
  failed rounds on it, so the second of the two to run would not run if
  the first failed. The exploration measured this lever at 1.27–1.37× in
  cold time at the two sizes it ran, against R02's 1.16–1.27× for the
  kernel. The larger lever goes first.
- **If R05 is rejected**, the scan has failed twice running (R02, R05).
  It is set aside, and R02b does not run.
- **If R05 is accepted**, R02b runs on R05's candidate as its base. A
  later failure there is then not a second consecutive one.

## Exploration before this declaration (disclosed)

[`../../explorations/R05-filter-20261001/README.md`](../../explorations/R05-filter-20261001/README.md)
records it. In short:
- **When:** 2026-10-01, 17:06–17:22 UTC, under the benchmark lock.
- **What:** `M1`'s two rows at `2^39.0`, `2^44.3` and `2^47.2`, three
  rounds, ABAB, isolated. That is 37 processes, one of them contended and
  re-run.
- **The arms:**
  - the base was `eaa842aa`, the tip of Track B's local stack;
  - the prototype was `0318f016`, that plus the filter change.
- **What it found:** identical outputs everywhere, and cold time 1.125
  [0.976, 1.296], 1.274 [1.211, 1.340] and 1.374 [1.315, 1.435].
- **What it set here:** the predictions and the power estimate. The rule
  is R02b's, unchanged.
