# R11: the pair-table build by partition streams and per-run filters, on v6

**Declared 2026-10-16, before any R11 timed run.** This is Track A, the
pair-table build (A1) in `research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§8. Its exploration on R08's candidate said a round was to be declared on
it "on the baseline that R08 decides". R08 decided v6. Nothing below
changes after the first timed run, except by a dated amendment appended at
the end.

## Why now

**At v6 the build is a third of the set-up at the 2^44 sizes.** R08's
candidate arm is v6. Over every repetition of its timed runs, the set-up
divides as follows:

| curve | `log₂ r` | set-up, median | collection | build | select |
|:--|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 560 ms | 55.2% | 34.6% | 4.5% |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 578 ms | 52.8% | 35.1% | 5.2% |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 2,002 ms | 78.0% | 17.6% | 2.0% |

The linear algebra is under 0.5% at every size. The rule at v6 left rho
7.90×, 7.40× and 12.06× faster cold at these sizes. Collection and the
build are what is left to cut.

## Hypothesis

**The build's partition streams and per-run filters make the build about
1.4× faster, with the same table word for word.**
- **Streams:** each pair's `(hash, word)` goes straight into one of 256
  coarse partitions, appended to its run's stream. This replaces a per-row
  list that is counted and then scattered. The rows are cut into contiguous
  runs, one a thread, so the table is word for word the same.
- **Per-run filters:** each run sets a presence filter of its own with
  plain ORs, prefetching the word sixteen pairs ahead, and the runs'
  filters are ORed at the end. At one thread there is one run, and its
  filter is the filter.

**The outputs are identical:**
- `the_partitioned_fold_build_is_one_table_word_for_word_on_any_thread_count`
  compares the table and the filter word for word at one, two, three and
  five threads, including at `2^59` and `2^61`;
- the pin checks every output against v0's.

## Class

Engineering, if it gains: the table, and with it every count and every
answer, is unchanged by construction.

## The candidate

**Frozen here as [`candidate-v6.patch`](candidate-v6.patch)** (SHA-256 in
`candidate-v6.patch.sha256`): v6 to local commit
`c2af9b7ccd92469830048dc2deaf97e1d194b27a`, one file,
`koblitz_index_calculus.rs`.
- It is the exploration's
  [`build-candidate.patch`](../../explorations/A1-build-on-fold-20261006/build-candidate.patch),
  applied to its declared base and cherry-picked onto v6. The base is
  R08's candidate, rebuilt from the frozen patches.
- **It merges without a conflict and adds and removes exactly the same
  lines.**
- **No run-time detection:** it changes the build's data flow, not a
  kernel, so it runs on every host. The claim is made for the reference
  class only, where it is measured.

## Arms

- **The base: v6,** the binary R08 measured, SHA-256
  `4ccb1951134ff31b5183fc00033933b9add2e3c0c52003f9041fd04bac1f4012`
  (local commit `3d4c5c0e`, tree `ab57f3fc`).
- **The candidate:** v6 plus `candidate-v6.patch`, built with
  `rustc 1.97.0` and main's lock (`--locked`), and `IC_BUILD_COMMIT` set.
  Its SHA-256 is recorded in an amendment, and in the manifest, before the
  first step runs.

## Before any timed step

1. **The tests,** on the candidate's tree:
   `cargo test --release --locked --lib -- koblitz_fast koblitz_index_calculus`,
   logged in `tests.log`.
2. **The pin, untimed:** the candidate on all 90 suite rows. Every output
   must equal v0's.
3. **The A/A:** v6 against a byte-identical copy, on `M1`'s 22 rows, five
   rounds.

## Rows

- **The suite rows:** `M1`'s 22 rows, and the other six suite rows at
  each target size: 40 rows, five rounds, ABAB, isolated, 400 processes.
- **The fresh holdouts:** eight rows at each target size, drawn by
  `icprog holdouts` (which first reproduced R10's holdouts byte for byte on
  this host):
  - recipe seeds 234 to 237, two targets a seed, `T159` to `T166`;
  - `public_hash_seed` 23159 to 23166, and rho seeds `0x230000 + 159` to
    `+ 166`.

  They are committed in [`holdouts/`](holdouts/) (`SHA256SUMS`). No round,
  exploration or diagnostic has run them. 24 rows, five rounds, ABAB,
  isolated: 240 processes.
- **The extension.** A set whose interval half-width exceeds 3% after five
  rounds gets rounds 6–10, pooled.

## Prediction

**From the build's share at v6 and the exploration's build gain**
(1.395×, 1.430× and 1.473×), with collection unchanged:

| curve | `log₂ r` | exploration, cold [95%] | predicted |
|:--|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.119 [1.011, 1.237] | **1.08–1.12×** |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.101 [1.016, 1.192] | **1.08–1.12×** |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.071 [0.986, 1.164] | **1.04–1.07×** |

Elsewhere the build is a smaller or similar share. No size is predicted
to regress.

## Success and stop

**Accepted, as baseline v7,** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At each target size, the paired cold-time ratio's 95% interval lies
   above 1.03, on the suite rows and on the fresh holdouts separately.
3. No size regresses beyond its A/A band.

**Rejected** otherwise. **Stopped** if a test fails, an output differs,
or a verification fails.

**If accepted:**
- the candidate becomes v7, with its binary hash, host and `S` from this
  round's runs;
- the patch lands in the tool, byte-identical to the file measured.

## Inadmissible

- Reusing the exploration's runs, pooling them with R11's, or choosing
  between them.
- Changing the candidate, the suite, the recipes or the targets, or
  drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting the build or collection as the speedup: the speedup is cold
  time.
- Crediting the gain to one of the two parts: the round measures them
  together.

## Cost

- **The tests, the pin and the A/A:** about 90 untimed and 220 timed
  processes, about an hour.
- **The suite rows and holdouts:** 640 timed processes, about two and a
  half hours, and up to 640 more if extended.

## Order and restarts

- **R11 runs on v6, on this host,** next.
- **The commands** are R08's, with `r11` for the round, through
  [`run.sh`](run.sh) (`run.sh <scratch dir> <candidate commit>`), as a
  session-tracked task.
- **After a container restart,** a resume runs
  `icprog run r11 manifest-resumed` first.
