# Round 0011 pre-registration: the general arithmetic in words

Written before any measured stage of this round was run.

## Where the job's cost sat after round 0010

With the evaluator's fork gone and the executable static, the archived
round-0010 receipts put the largest phase of the IC arm, on every cell, in
the **general-arithmetic final check** that both arms run on every recovered
scalar: 0.78M of its 1.87M instructions on `n13a0`, 2.83M of 6.82M on
`n23a0`. The **curve construction and target lift** came second (0.52M and
1.39M). By function, the final check spent 53% of its instructions in the
general field's squaring, 18% in its schoolbook multiplication and 20% in
heap allocation: `F2mElement::square` set one output bit per input bit and
`reduce` walked every bit position from the top of the product down, one
low term at a time. The construction rebuilt the single-word curve's
reduction tables for every candidate abscissa of the generator search and
again for every target lift, re-ran the sparse irreducible search (a sort
of about 1,500 candidate masks and the irreducibility tests) on every job,
took the eigenvalue's square root in big integers, and probed the CPU
feature cache of the standard library inside the general `Gf2` constructor
(73–97 µs on this virtual machine, for both arms, which is why round
0010's ablation of the IC arm's own probe measured nothing: the shared one
had already paid it).

## What changes, and for whom

**Baseline source** (`--source-root`, [patch](round11-wordfield.patch) on
the round-0010 incumbent source), shared by every arm including rho:

- `binary_ecc::f2m`: squaring by bit spreading a half-word at a time,
  multiplication as the carry-less product word by word (`pclmulqdq` where
  present, a word loop otherwise), and reduction a word of high bits at a
  time. The public API and every result are unchanged: the previous
  bit-serial implementation is retained under `#[cfg(test)]` as the
  reference, and a test compares square, product and repeated squaring on
  31 widths from 1 to 571 bits, random elements and random moduli, on both
  carry-less paths.
- The construction builds one single-word curve per job and passes it to
  the point count, the generator search and the eigenvalue check; the
  worker's target lift uses one for all its attempts. The sparse irreducible
  polynomial is tabulated for every degree below 64, with a test that
  re-runs the search for each degree and pins the table to it. The
  eigenvalue's square root and its two candidates are computed in one word
  when the subgroup order fits, step for step the same Tonelli–Shanks as the
  big-integer routine (the same root, not merely a root), tested equal on
  every worker cell and on assorted primes of both residue classes.
- One cached `CPUID` probe (`f2m::pclmul_available`) serves the general
  field, the single-word `Gf2` and the tiny IC.

Rho is the shipped solver, untouched; it runs on the single-word curve and
uses the general arithmetic where it always did, so it also gets faster.

**Challengers:** `wordfield_fastio` ([patch](round11-fastio.patch)), the
IC-only change: the report's decimal numbers are written two digits at a
time from a table and appended without a UTF-8 validation pass; a test
compares the text with the standard formatting on thousands of values. In
round 0010 that serialisation was 65% of the `log_certification` phase
(0.108M of 1.87M instructions on `n13a0`, 0.265M of 6.82M on `n23a0`).
`wordfield_batch4`, the configuration control (`batch_trials: 4`).

Development evidence before freezing (not a claim): the baseline executable
reproduces, byte for byte apart from `elapsed_seconds`, all 304 frozen
fixture outputs of round 0010 (152 IC and 152 rho jobs) and every output
passes the independent checker. Under Callgrind on one confirmation case per
cell the final check fell to 0.20M–1.09M instructions (from 0.54M–3.26M), the
construction to 0.19M–0.54M (from 0.51M–1.39M), the IC arm's whole job to
0.50–0.60 of its round-0010 instructions and rho's to 0.59–0.66 (rho's own
solve also shrank, 1.47M→1.14M on `n13a0`, 7.81M→5.68M on `n23a0`, where it
uses the general arithmetic). Native, on two confirmation fixtures per cell
with ten repetitions under the evaluator's spawn, the IC arm took 1.55–1.97 ms
(from 1.74–2.55) and rho 1.75–2.35 ms (from 1.93–3.23); the IC/rho ratio was
0.83–0.89 on every cell. The final check remains the largest IC phase on the
largest cell (its cost is now the per-operation heap allocation of the
general element type), which is the natural target of a further round.

The whole library test suite was run on the baseline tree: 2,425 tests
pass and 5 fail, all in the general pair-table tests of
`koblitz_index_calculus` (`every_stored_pair_sum_is_found_by_lookup`,
`pair_table_agrees_with_enumeration_for_two_three_and_four_summands`,
`the_compact_table_answers_exactly_as_the_full_one`,
`the_compact_table_decomposes_at_every_summand_count`,
`the_folded_table_answers_exactly_as_the_full_one`). The same five fail
identically on the unmodified round-0010 source (the same assertion values),
so they are a pre-existing condition of the frozen research snapshot, which
descends from the round-0005 lineage rather than `main`, and not a regression
of this round; they exercise the general pair table the tiny IC does not use.

## Objective and gates

`--objective rho`, as rounds 0007–0010. Incumbent: the round-0010 winner
configuration (`batch_trials: 1`) on the baseline above. `wordfield_fastio`
is promoted only if it passes the no-regression gate against the incumbent
(instruction ratio at most 0.98 with the paired upper limit below one, native
upper limit below one, no cell more than 10% worse) and the strict rho gate
on confirmation and replay; its expected gain is 2–6% of instructions, so a
promotion is plausible on the small cells' weight and a *retained* verdict
is a fully reported outcome. The control is expected to be retained. Either
way the decision records the winner's `winner_over_rho`, `beats_rho_strict`
and `rho_parity`.

## Parent, seed, budget

Fresh seed 2026091611; target count 1; pilot profile; 1,800 paired-job
budget; one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0 `Ir`;
native progress recorded; the round-0010 evaluator (spawn without a fork,
caps before job delivery).

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps and report fields
are unchanged; the final check is still `[d]G = Q` in the general
arithmetic on every recovered scalar, only the arithmetic's implementation
moved. Class: engineering. Absolute times and instruction counts are
comparable with round 0010 (same evaluator, same host, same compiler) and
not with earlier rounds. No arithmetic-complexity, family-wide or
cryptographic-size claim follows. Fresh fixtures from the new seed; every
failure retained; panels stay separate.
