# P-256 Dickson branched quadratic-S3 protocol, round 4

Date frozen: 2026-10-04

Round 3 lowered the input degree to 2 but retained solving degree 4.  This
round uses the triangular Dickson structure rather than treating the full
fibre as one ideal: it branches on one or two terminal preimages, removes the
fixed tail variables, and solves every resulting component.

## Hypothesis

For a depth-5 Dickson fibre, fixing the boundary value after the first
`5-s` recurrence steps splits membership into `2^s` depth-`5-s` fibres.  Two
summands therefore require `4^s` component systems per target.  With the
round-3 grouped quadratic `S3` formulation, either `s = 1` or `s = 2` has:

1. maximum certified component solving degree at most 3;
2. correct verdicts for every positive and negative component;
3. total field operations, summed over all components, explicitly charged.

The primary success condition is the degree reduction.  A reduction obtained
by omitting negative components, using knowledge of the witness branch, or
failing to decide a component is inadmissible.

## Frozen inputs

- dependency: round 3 ending at local commit `c2ff3ce6f`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy field and curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- depth: 5;
- terminals: `c = 0` and round-2 winner `c = 369`;
- branch depths: `s = 0, 1, 2`, where `s = 0` is the unbranched reference;
- targets: the first three common-positive affine targets from round 3,
  recomputed exhaustively;
- equations: round-3 grouped direct quadratic `S3`;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 12;
- per-component stop: 30 seconds.

For branch depth `s`, enumerate every `r in F_p` satisfying
`D_(2^s)(r) = c`.  Replace the last `s` chain variables by the chosen boundary
and retain a depth-`5-s` chain ending at `r`.  All boundary pairs are run in
increasing `(r1,r2)` order.

## Correctness and accounting

For each boundary, enumerate its signed liftable points.  Exhaustive addition
classifies each `(r1,r2,target)` component as positive or negative.  A
completed positive F4 component must be consistent; a completed negative
component must be inconsistent.  Any certified mismatch falsifies the round.
An incomplete component makes that branch depth inconclusive.

The result table has one row per `(terminal, branch depth)`:

| terminal | branch depth | components/target | targets | positive components | negative components | correct | complete | max solving degree | median solving degree | max columns | total field ops | ops / unbranched | classification |
|---:|---:|---:|---:|---:|---:|:--:|:--:|---:|---:|---:|---:|---:|:--|

The degree boundary is the unbranched grouped quadratic system's degree 4.
Total operations include every positive and negative component for all three
targets.  Wall time is secondary only.

## Stop conditions and scope

Stop after the two terminals, three branch depths, and three targets.  Preserve
all extension-field/non-refutation cases as certified mismatches or
inconclusive cells according to the F4 report; do not reinterpret them as
finite-field proofs.  Commit raw JSON, the table, commands, and hashes.

This native Rust experiment measures only small-prime component ideals.  It
does not run a P-256 Gröbner basis or establish an ECDLP speedup.  A degree win
would justify a separately priced higher-depth branching study; it would not
by itself transfer to depth 18.
