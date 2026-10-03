# B0: refuse what is wrong today

**Declared 2026-10-01, before any B0 measurement.** This is Track B's
first round in `research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9,
B0). Nothing below changes after the first measurement, except by a
dated amendment appended at the end.

**Disclosure: the code came first.** B0's code was written before this
protocol, as local commit `a1187008` (2026-10-01, 01:55 UTC). It was
ported unchanged onto `main` at `c1a2e5f8`, without conflicts. What this
protocol fixes in advance is everything that decides the round: the
measurements, the acceptance rule and the stop rules. No B0 measurement
ran before this protocol was committed. The unit tests ran while the
code was written, as tests do.

## What B0 does

Each defect the plan's survey lists (§9) becomes a precise refusal or a
fix. Each has a case in conformance suite v1
(`research/ic_tool_program/conformance/v1/`) or, where no input can
reach it, a unit test.

| defect (plan §9) | fix | checked by |
|:--|:--|:--|
| `r ≥ 2^64` is read as its low word | the single-word paths take `r` whole or not at all: `subgroup_order_word` refuses `r ≥ 2^63` at all four sites | unit test `subgroup_order_word_refuses_what_it_cannot_hold`; no input reaches it while `n ≤ 63` |
| the synthetic-degree refusal says `3..=62` while 63 passes | the message names the range it enforces | C001; unit test `degree_refusal_names_the_range_it_enforces` |
| a foreign `state.json` with a short digest panics | the digest is shown as far as it goes | C002 |
| a huge `points` overflows in debug builds and draws for hours first | `points` outside `1..=2·MAX_ABSCISSAE` is refused before drawing; the budget saturates | C003 |
| large random-binary degrees hang in point counting | `ic boundary`, `rho` and `bench` refuse degrees outside `5..=32` (`MAX_ENUMERATED_DEGREE`), with the reason | C004, C005, C006 |
| a too-small `pair_table_bytes` is reported as "field too wide" | the refusal names the budget, or the field when the field is the cause | C007 |
| `ic search`'s validation child loses `--subfield`, `--curve-b`, `--wdsat-binary` and `--wdsat-timeout-ms` | one function builds the child's command, with every flag that decides the curve, the base, the solver or the target | unit test `validation_child_carries_the_curve_and_the_solver` |
| a report's `commit` is the working directory's, not the binary's | a binary built with `IC_BUILD_COMMIT` reports it as `build_commit`, beside `commit` | C008 |

C008 is also the suite's positive control: a valid file runs to a
verified scalar.

## Class

**Robustness.** B0 changes what the tool refuses and what it reports,
not what it computes. It claims no speedup.

## Arms

- **The base**: the newest accepted baseline when B0 runs, as a commit,
  recorded in B0's manifest. It is v1 if R02 is accepted, and the
  baseline after R03 if R03 is accepted first.
- **B0**: that commit plus B0's change, built with
  `IC_BUILD_COMMIT=$(git rev-parse HEAD)`.

B0 touches neither `koblitz_fast.rs` nor `factorise_u64`, so it merges
with R02's and R03's changes without interaction.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   the two modules B0 touches (`ic_boundary`,
   `koblitz_index_calculus::tests`), on B0's tree.
2. **Conformance.** Suite v1 on B0, with `--build-commit` set to B0's
   commit. Suite v1 also runs on the base, without `--build-commit`; its
   failures are recorded as the defects as they were.
3. **The pin, untimed.** B0 runs on all 90 suite rows, and every output
   must equal v0's from R01's profile pass. These are the counts, both
   arms' scalars, rho's counts and the verification flags.
4. **No slowdown, timed.** The base against B0 on `M1`'s 22 rows, five rounds
   ABAB, isolated: 220 processes, with the programme's runner. The
   figure is the paired cold-time ratio per size.

## Acceptance

B0 is **accepted** when all of the following hold:
- every test passes;
- the conformance suite passes 8 of 8 on B0;
- the pin holds on every row;
- no size regresses beyond its A/A band. A regression is a size whose
  interval's upper end lies below the lower end of that size's A/A
  interval: R01's, or the current round's own when the host has
  changed.

It is **stopped** if any output differs, since a refusal round must not
change an answer. It is **rejected** if a test or case fails, or if a
size regresses. A rejected B0 is fixed and declared again.

## Inadmissible

- Loosening a conformance case's expectation after a run.
- Dropping a defect from the table.
- Counting a contended run.

## Cost

- The pin: 90 untimed processes.
- The timing check: 220 timed processes, about 40 minutes.
- The tests and the conformance suite: a few minutes.

## Amendment 1 (2026-10-01, before any measurement)

**Measurement 4, the timing check, runs in the chain of Track B steps.**
- **What it was:** the base against B0 on `M1`'s 22 rows, five rounds
  ABAB, 220 processes, for this step alone.
- **What it is now:** one interleave over the newest accepted baseline
  and every Track B step's arm, in the queue's order: the baseline, B0,
  B1, B3, B2, B2b, B7a and B3b (`bround.py chain`). It runs on `M1`'s 22
  rows, five rounds, isolated, and the order reverses every other
  round.
- **B0's figure** is the paired cold-time ratio of the arm before it
  in the chain over B0's own, per size. It is read against R01's A/A
  bands, as before.
- **Why the chain is a valid base.** Each step's arm is built on the one
  before, so the arm before B0 in the chain is its base. The two run
  back to back in every round, with ten pairs a size, as in its own
  ABAB.
- **What it saves.** The seven steps' checks take 880 processes
  together, against 1,540 one by one.
- **If an earlier step is rejected,** the steps after it wait. The chain
  runs again from that step, once it is fixed or removed, since the
  later arms carry its change.
- **Nothing else changes:** the acceptance rule, the A/A bands and the
  other measurements.
