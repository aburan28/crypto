# R07: the programme re-based on main's head

**Declared 2026-10-05, before any R07 timed run.** This is Track A in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md`. Plan §5 builds the
baseline "from the unmodified head", and main's head has moved on from
the programme's baseline. A diagnostic ran first and is disclosed below.
Nothing below changes after the first timed run, except by a dated
amendment appended at the end.

## Why a round, and why now

**The programme's baselines have been a lineage of their own.**
- v0 was main at `46ae2014`.
- v1 (R03, `30f6c153`) and v2 (R05, `edcb0bec`) are candidates built on
  it. Their code reached main through their results pull requests
  (#1166, #1187).
- Every round since R01 has compared two binaries of that lineage.

**Main has moved on.** Since #1195 it has gained about 700 commits.
- **#1242** (merged 2026-10-02 as `809a7318`) changed the field
  arithmetic every phase uses. `semaev_decomp::Gf2` now reduces a product
  by two carry-less folds with the modulus's sparse tail, instead of a
  byte table, whenever the tail is short enough. Every suite field's tail
  is. Its own review measured, on perfbench at `n = 53`:
  - 2.44× on multiply and square, 2.49× on inversion and 1.66× on
    batched inversion;
  - 1.64× on the pair-table build kernel.
- The suite's two wide-tail sizes gain most. There the scan's
  subtraction runs the scalar path (R02, R04), and #1242 inlines its
  products.
- **Main's `ic` is what anyone building this repository runs.** A round
  measured on v2 now measures a tool nobody builds.

**What this round is.** It asks whether main's head computes v2's answers,
and what it costs, size by size. It claims no lever of its own: any gain
is main's, chiefly #1242's.

## Hypothesis

- **Main's head (`995ea207`) gives v0's outputs on every pinned row:**
  - the counts, relations and logarithms;
  - rho's walks and both arms' scalars;
  - the verification flags.
- **It is faster than v2 at every size:**
  - most at `icv1-f2m59-tm943548413-98844ecc` and
    `icv1-f2m61-t158598901-ab42b6c5`, whose scan runs the scalar
    subtraction;
  - least at `icv1-f2m53-tm56619371-dac20a85`, whose scan runs the 8-lane
    kernel, which #1242 does not touch.

## Class

**Engineering, if it gains.**
- The pin makes the counts identical, so the ratio to the floor is flat.
- The ledger credits the gain to main (#1242 and the rest of the drift),
  not to a lever of the programme's.

## Arms

- **The base: v2.**
  - Commit `edcb0bec`, binary SHA-256 `76a2a2fd05acd47d…` (`baselines.json`).
  - It is the same binary R02b's stopped runs used (R02b amendment 3).
- **The candidate: main at `995ea2071cc7453877a503d30eb561ae82cddab9`.**
  - Built with `cargo build --release --bin ic` and v2's `Cargo.lock`
    (SHA-256 `4f17b356fa7bac39…`, the lock R02b's manifest records), so
    no dependency moves between the arms.
  - Binary SHA-256
    `9c3320390d74ae48b74382d2ec3f8b4b5b1585e5e8aa2448c0d78369807a26b4`,
    with `rustc 1.94.1`, as v2 was built.
- **The manifest** (`icprog run r07 manifest`) records both binaries'
  SHA-256 and the commits they were built from.

## Before any timed step

1. **The tests.** On the candidate's tree:
   `cargo test --release --lib -- koblitz_fast koblitz_index_calculus semaev_decomp`.
   These are main's own tests of every kernel on the timed path, among
   them #1242's tests pinning each field operation to the previous
   code's.
2. **The pin, untimed.** The candidate runs on all 90 suite rows
   (`icprog run r07 pin`). Every output must equal v0's from R01's
   profile pass, and the names must agree, as in R02b and R05.
3. **The A/A.**
   - **Why.** This container's kernel build is `6.18.44-fc-v70`, and
     R01's was `-fc-v50`, so the bands are the round's own.
   - **What.** The base against a byte-identical copy, on `M1`'s 22 rows,
     five rounds, ABAB, isolated (`icprog run r07 aa`): 220 processes.

## Rows

- **The suite rows (40 pairs at each target size).**
  - `M1`'s 22 rows: two targets at each of the eleven sizes.
  - At the three target sizes, the other six suite rows each (`M2`–`M4`).

  That is 40 rows. Five rounds, ABAB, isolated: 400 processes.
- **The fresh holdouts (40 pairs at each target size).** These are eight
  rows at each of the three target sizes:
  - recipe seeds 214, 215, 216 and 217, two targets per seed, `T119` to
    `T126`;
  - `public_hash_seed` 23119 to 23126, and rho seeds `0x230000 + 119` to
    `+ 126`.

  They were drawn by `icprog holdouts`, the native port of suite v1's
  construction. Before drawing them, it reproduced byte for byte:
  - suite v1's 88 S rows;
  - R02b's and R05's holdouts, with their `SHA256SUMS`.

  They are committed in [`holdouts/`](holdouts/) (`SHA256SUMS`). Neither
  their seeds nor their targets have been run by any round or by the
  diagnostic. Five rounds, ABAB, isolated: 240 processes.
- **The extension.** At a target size, a set whose interval half-width
  (`hi / geomean − 1`) exceeds 3% after five rounds gets rounds 6–10, and
  its figure pools all ten. This depends only on the interval's width,
  never on its position.
- **The callgrind cross-check (plan §5).**
  - **What runs.** `ic workflow` under callgrind on `M1`'s first target
    at `icv1-f2m61-t158598901-ab42b6c5` and
    `icv1-f2m53-tm56619371-dac20a85`, both arms, untimed
    (`icprog run r07 callgrind`).
  - **What valgrind sees.** It hides AVX-512 from both arms. It does run
    their PCLMULQDQ paths, so #1242's reduction is in the counts.
  - **What it reports.** The instruction ratio and the functions whose
    counts differ.
  - **What it decides.** Only one thing: both arms must recover the same
    logarithm. It is a deterministic cross-check of the timed figures,
    not a control.

## Prediction

**From the diagnostic below**, set-up time, base over candidate, one
process each:

| curve | `log₂ r` | v2 set-up | main set-up | ratio |
|:--|--:|--:|--:|--:|
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 36.8 ms | 28.0 ms | 1.32 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 852 ms | 809 ms | 1.05 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1225 ms | 888 ms | 1.38 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 4908 ms | 3190 ms | 1.54 |

- **Cold time, both sets:** about **1.05×**, **1.38×** and **1.54×** at
  the three target sizes. These are the diagnostic's ratios: one process
  per arm, so they are a prediction, not an estimate.
- **Elsewhere:** above 1 at every size, about 1.1–1.35×, from the build
  and the field operations every phase shares.
- **The diagnostic's one surprise:** at
  `icv1-f2m53-tm56619371-dac20a85` collection read 459 ms on v2 and
  509 ms on main, while build fell from 296 to 241 ms. One process per
  arm cannot tell a regression from noise. The A/A rule below decides
  it.

## Success and stop

**Accepted, as baseline v3** (main's head), if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. No size regresses beyond its A/A band. A regression is a size whose
   cold interval's upper end lies below the lower end of the round's own
   A/A interval.
3. Under callgrind, both arms recover the same logarithm at both sizes.

There is **no minimum gain**:
- the programme moves to the code its users build whatever the gain;
- the figures report, size by size, what that move did.

**Rejected** if a size regresses beyond its A/A band. Then:
- v2 stays the baseline;
- the regression is located (callgrind phase split, R04's probes if
  needed) and fixed on main by its own pull request;
- R07 is re-declared on a head that carries the fix (R07b).

No Track A round runs on v2 meanwhile.

**Stopped** if a test fails, an output differs, or a verification fails.
- **A pin failure means main's head changed an output.** The programme
  cannot re-base on main until the change is identified, and either:
  - declared as an algorithm change, with its own measurement; or
  - fixed.
- The differing rows are reported in R07's results.

**If accepted:**
- **v3 is entered.** The candidate becomes v3 in the ledger and in
  `baselines.json`, with its binary hash, host, and `S` from this
  round's runs, beside v2's `S` from the same runs.
- **The rule's comparison runs at v3** (plan §5).
  - It uses rule v2's declared protocol with v3's binary, declared in
    R07's results pull request as `rule/v3`.
  - Rule v2 was declared for v2 and has not run. It stays on record
    unrun, superseded by v3's.
  - If R07 is rejected, rule v2 runs as declared.

## Inadmissible

- Pooling the diagnostic with R07's runs, or choosing between them.
- Changing either binary, the commit or the lock after a run.
- Changing the suite, the recipes, the holdouts or the targets, or
  drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting collection, build, rho or a perfbench kernel as the speedup:
  the speedup is cold time.
- Crediting the gain to the programme: it is main's.

## Cost

- **The tests and the pin:** a few minutes, and 90 untimed processes.
- **The A/A:** 220 timed processes, about an hour.
- **The suite rows and holdouts:** 400 and 240 timed processes, and up
  to 240 more if extended. About three hours together.
- **The callgrind cross-check:** four runs under valgrind, about an hour.
- **If accepted, the rule's comparison at v3:** about 40 minutes, and
  the reference check after it.

## Order

**R07 runs first.**
- Before R02b, which amendment 3 suspends in this pull request.
- Before R06, the funnel-shift key. R06 is declared after R07's decision,
  on the newest baseline.

**The container has restarted mid-run three times since 2026-10-01**
(R05's holdouts, R02b's two attempts).
- Every step resumes what its run tree holds.
- Each resume runs `icprog run r07 manifest-resumed`. The first writes
  `host-resumed.json`, and each later one writes
  `host-resumed-<k>.json`. The analysis lists them all.
- A pair that straddles a resume is paired as R05's was, and the results
  say so.
- A resume onto a host whose CPU model, flags, cores, memory or kernel
  build differs from `host.json`'s stops the round. The rest then
  re-runs from the A/A on the new host, as a dated amendment.

## The diagnostic before this declaration (disclosed)

**When and how.** 2026-10-05, about 19:20 UTC.
- One process per arm on `M1`'s first target at four sizes.
- Each ran under `taskset -c 2` and the benchmark lock (`isolated_bench
  busy`), on an otherwise idle machine, but not through `isolated_bench
  run`. It is a diagnostic, not a measurement.
- The binaries are the arms above, built from the same commits and
  lock.

**What it found.**
- **The answers matched.** At all four sizes the counts, rho's counts
  and every certificate were identical between the arms, and so was the
  recovered logarithm.
- **The set-up times** are the prediction table above.
- **Collection per scanned summand on main** was 35–38 ns at all three
  target sizes: 2615 ms over 74.4 M summands at
  `icv1-f2m61-t158598901-ab42b6c5`. On v2 it was 33–58 ns, with the
  scalar wide-tail sizes at the top. That is why R02b's amendment 3
  re-asks R02b's question on the new base before running it.
- **What the diagnostic set:** the prediction above. The rows, the rule
  and the A/A requirement were not chosen from it.

## Commands

```
# the candidate
git worktree add main-wt 995ea2071cc7453877a503d30eb561ae82cddab9
cp <v2's Cargo.lock> main-wt/Cargo.lock
cargo build --release --bin ic          # in main-wt; SHA-256 9c332039…

# the round, from a checkout of this pull request's merge
for step in manifest pin aa compare holdout extend callgrind; do
  icprog run r07 $step --base <v2's ic> --cand <main's ic> \
    --isolate <isolated_bench> \
    --base-commit edcb0bec948ffab566a3c8998dc3c92d52375da2 \
    --cand-commit 995ea2071cc7453877a503d30eb561ae82cddab9
done
icprog analyse r07 > research/ic_tool_program/rounds/R07-main-head/analysis.json
```
