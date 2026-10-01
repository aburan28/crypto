# R03: the curve's construction at composite degrees

**Declared 2026-10-01, before any R03 run.** This is Track A, A9, in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md`. Nothing below
changes after the first run except by a dated amendment appended at the
end.

## Hypothesis

**At a composite degree, building the curve is mostly factorising its
order, and the factorisation tests primality at every trial divisor.**

- `factorise_u64` (`src/cryptanalysis/koblitz_index_calculus.rs`)
  divides by `d = 2, 3, 4, …`. Before each division it runs a
  deterministic Miller–Rabin test, twelve bases with 128-bit modular
  products, on what is left.
- At a prime degree the order is a small cofactor times the prime `r`,
  so what is left is prime almost at once.
- At a composite degree the cofactor carries the orders of the subfield
  curves. The loop then runs up to their largest prime, with a full test
  at every step.
- **R01's profile shows the cost.** At `icv1-f2m57-tm747311035-c1f545af`
  (`n = 57`, `2^38.0`), the set-up phase that builds the curve is 45% of
  the set-up: about 53 ms of 118 ms. At `icv1-f2m45-tm6236725-40939294`
  (`n = 45`, `2^24.8`) it is 25%: about 0.8 ms of 3.1 ms. The suite's
  other nine sizes are prime degrees.
- **An exploratory instruction count, run before this protocol.**
  Callgrind of `ic workflow --stop-after select` on
  `icv1-f2m57-tm747311035-c1f545af` `M1-T01`, untimed:
  - `factorise` is 71% of the 430 M instructions up to factor-base
    selection;
  - `__umodti3`, the 128-bit remainder inside the test's products, is
    130 M of them.
  
  This informed the prediction. It is not evidence for the round, and no
  candidate was timed.

## The candidate

- Test primality when a division changes what is left, not at every
  divisor.
- Skip even divisors after 2.

The factors are identical, divisor for divisor. A unit test compares the
new function with the old one, kept verbatim, on every value it covers,
including the orders of `K_0` and `K_1` for `n = 3..=63`.

**Disclosure: the code came first.** The candidate was written before
this protocol, as local commit `30f6c153`, and its unit test ran then.
No R03 measurement has run.

## Class

Engineering, if it gains. Nothing the pipeline decides changes, so the
ratio to the floor is flat.

## Arms

- **Base:** R02's outcome. That is v1, R02's candidate commit
  (`a1645ac6`), if R02 is accepted, and v0′ (`c1a2e5f8`) otherwise. The
  choice is fixed by R02's recorded decision, before any R03 run, so
  that R03 need not wait for R02's pull request to merge.
- **Candidate:** the base plus this change, which touches only
  `src/cryptanalysis/koblitz_index_calculus.rs`. R02's change touches
  only `koblitz_fast.rs`, so the two merge without interaction.

## Pinned outputs

Before any timed step, the candidate runs all 90 rows untimed. Every
output must equal v0's from R01's profile pass:
- the counts;
- both arms' scalars;
- rho's counts;
- the verification flags.

R02's pin already ties the base to v0.

## Timed comparison

All runs use the programme's runner, isolated, on one thread, and alternate the order ABAB.

- **The two composite sizes, every row.** That is
  `icv1-f2m57-tm747311035-c1f545af` and
  `icv1-f2m45-tm6236725-40939294`: 16 rows, five rounds, 160 processes.
- **The nine prime sizes, `M1`'s rows**, as the check that nothing else
  moved: 18 rows, five rounds, 180 processes.
- **The holdouts,** at the two composite sizes: R02's frozen files for
  seed 205, targets `T101` and `T102`. That is 4 rows, five rounds, 40
  processes.
- **The A/A:** R01's, if the host manifest matches; otherwise R03's own,
  on `M1`'s rows.

## Prediction

- **`icv1-f2m57-tm747311035-c1f545af`: 1.75×** [1.5, 2.0] in cold time,
  base over candidate. The construction phase falls from about 53 ms to
  under 2 ms.
- **`icv1-f2m45-tm6236725-40939294`: 1.3×** [1.1, 1.5].
- **The nine prime sizes: 1.00×.** Their construction never reaches the
  loop.

## Success and stop

**Accepted** if all of the following hold:
1. Every pinned output is identical.
2. At `icv1-f2m57-tm747311035-c1f545af`, the paired cold-time ratio's
   95% interval lies above 1.3. This must hold on the suite rows and on
   the holdouts separately.
3. No size regresses beyond its A/A band. A regression is a candidate
   interval whose upper end lies below the lower end of the size's A/A
   interval.

**Rejected** otherwise. The code is kept on record and the numbers go in
the ledger.

**Stopped early** if any output differs, or if a verification fails.

## Inadmissible

- Quoting the construction phase's own ratio as the speedup. The speedup
  is cold time.
- Claiming a gain at a prime degree.
- Pooling contended runs.

## Cost

- The pin: 90 untimed processes.
- The comparison: 380 timed processes, about an hour, half of it the
  isolation tool's settling between processes.
