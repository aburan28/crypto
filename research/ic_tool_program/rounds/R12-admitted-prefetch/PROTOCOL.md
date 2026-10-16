# R12: the collection scan's admitted keys, their runs prefetched, on v7

**Declared 2026-10-16, before any R12 timed run.** This is Track A, the
admitted keys' cost in the collection scan (A13) in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8. Its exploration on
R08's candidate said a round was to be declared on it "on the baseline
that R08 decides". R08 decided v6, and R11 then moved the baseline to v7,
so R12 runs on v7. Nothing below changes after the first timed run, except
by a dated amendment appended at the end.

## Why now

**At v7 collection is the larger share of the set-up at every target
size.** R11 made the build 1.50–1.51× faster and left collection where it
was. At v6 collection was 55%, 53% and 78% of the set-up and the build 35%,
35% and 18%, so at v7 collection is about 62%, 59% and 83%.

**Inside collection, the admitted keys are waiting.** On the folding-kernel
candidate, the admitted keys cost 4.4–6.6 ns of a scanned summand's
21–22.5 ns, and each admitted key 210–225 ns: about one random access into
the largest table. Nearly every admitted key is a filter false positive,
so nearly all of that time is spent waiting
([record](../../explorations/R08-fold-kernel-20261006/README.md)).

## Hypothesis

**A second pass that prefetches each admitted key's run makes collection
1.1–1.3× faster, with every lookup, match and order unchanged.**
- **Before:** the scan prefetches the presence filter ahead of its keys
  and, for an admitted key, the bucket entry. The run the entry points at
  is then read cold.
- **After:** between the filter pass and the lookups, a second pass
  (`prefetch_run`) reads each admitted key's bucket entry, cached by then,
  and prefetches the run's first and last words. The runs' misses overlap.
- **Nothing a prefetch reads changes a value.** The lookups, the matches
  and their order are what they were, and the pin checks every output
  against v0's.

## Class

Engineering, if it gains: no value the scan computes changes.

## The candidate

**Frozen here as [`candidate-v7.patch`](candidate-v7.patch)** (SHA-256 in
`candidate-v7.patch.sha256`): v7 to local commit `23138d20b1d64cbb7f54a27e41cd08e1d1a38b0c`, 38 lines in
`koblitz_index_calculus.rs`.
- It is the exploration's
  [`admitted-prefetch.patch`](../../explorations/A13-admitted-prefetch-20261006/admitted-prefetch.patch),
  applied to its declared base and cherry-picked onto v7. The base is
  R08's candidate, rebuilt from the frozen patches.
- **It merges without a conflict and adds exactly the same lines.**
- **No run-time detection:** a prefetch is an x86 hint. The claim is made
  for the reference class only.

## Arms

- **The base: v7,** the binary R11 measured, SHA-256
  `85116c2a0da19e5e2e31d7ec1ff8d58a28406bbf9ffe30d330476a67a182315a`
  (local commit `c2af9b7c`, tree `bdd975b7`).
- **The candidate:** v7 plus `candidate-v7.patch`, built with
  `rustc 1.97.0` and main's lock (`--locked`), and `IC_BUILD_COMMIT` set.
  Its SHA-256 is recorded in an amendment, and in the manifest, before the
  first step runs.

## Before any timed step

1. **The tests,** on the candidate's tree:
   `cargo test --release --locked --lib -- koblitz_fast koblitz_index_calculus`,
   logged in `tests.log`.
2. **The pin, untimed:** the candidate on all 90 suite rows. Every output
   must equal v0's.
3. **The A/A:** v7 against a byte-identical copy, on `M1`'s 22 rows, five
   rounds.

## Rows

- **The suite rows:** `M1`'s 22 rows, and the other six suite rows at each
  target size: 40 rows, five rounds, ABAB, isolated, 400 processes.
- **The fresh holdouts:** eight rows at each target size, drawn by
  `icprog holdouts`:
  - recipe seeds 238 to 241, two targets a seed, `T167` to `T174`;
  - `public_hash_seed` 23167 to 23174, and rho seeds `0x230000 + 167` to
    `+ 174`.

  They are committed in `holdouts/` (`SHA256SUMS`) before any run. No
  round, exploration or diagnostic has run them. 24 rows, five rounds,
  ABAB, isolated: 240 processes.
- **The extension.** A set whose interval half-width exceeds 3% after five
  rounds gets rounds 6–10, pooled.

## Prediction

**From collection's share at v7 and the exploration's collection gain**
(1.149×, 1.288× and 1.117×), with the build unchanged:

| curve | `log₂ r` | exploration, cold [95%] | predicted |
|:--|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.075 [0.984, 1.174] | **1.06–1.11×** |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.125 [1.018, 1.243] | **1.10–1.17×** |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.086 [0.975, 1.209] | **1.06–1.12×** |

The exploration's pairs are few and its intervals wide. No size is
predicted to regress.

## Success and stop

**Accepted, as baseline v8,** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At each target size, the paired cold-time ratio's 95% interval lies
   above 1.03, on the suite rows and on the fresh holdouts separately.
3. No size regresses beyond its A/A band.

**Rejected** otherwise. **Stopped** if a test fails, an output differs,
or a verification fails.

**If accepted:**
- the candidate becomes v8;
- the patch lands in the tool, byte-identical to the file measured.

## Inadmissible

- Reusing the exploration's runs, pooling them with R12's, or choosing
  between them.
- Changing the candidate, the suite, the recipes or the targets, or
  drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting collection as the speedup: the speedup is cold time.

## Order and restarts

- **R12 runs on v7, on this host,** after the rule at v7.
- **The commands** are R11's, with `r12` for the round, through
  [`run.sh`](run.sh) (`run.sh <scratch dir> <candidate commit>`), as a
  session-tracked task.
- **After a container restart,** a resume runs
  `icprog run r12 manifest-resumed` first.
