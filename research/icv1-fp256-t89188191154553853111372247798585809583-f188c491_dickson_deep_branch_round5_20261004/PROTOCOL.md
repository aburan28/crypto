# P-256 Dickson deep-branch degree protocol, round 5

Date frozen: 2026-10-04

Round 4 retained solving degree 4 through two terminal branch levels, while
reducing total field operations by about 90%.  This final bounded rung fixes
three and four terminal levels, leaving respectively two and one Dickson chain
variables per summand.

## Hypothesis

On at least one of branch depths 3 or 4, every certified component of the
grouped quadratic `S3` system has solving degree at most 3 for both terminal
zero and terminal 369.  Every component must be run and charged; the maximum
component degree, not the degree of a selected witness branch, decides the
hypothesis.

## Frozen inputs and procedure

- dependency: round 4 ending at local commit `297975e51`;
- registered curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- toy curve: `p = 1151`,
  `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- original Dickson depth: 5;
- terminals: 0 and 369;
- branch depths: 3 and 4;
- components per target: 64 and 256;
- targets: the same frozen first three common-positive points as round 4;
- equations and variable layout: round-4 grouped direct quadratic `S3`;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 12;
- per-component stop: 30 seconds.

Boundary enumeration, signed-point component classification, and finite-field
correctness checks are identical to round 4.  Positive components must remain
consistent; negative components must be refuted.  Any completed mismatch
falsifies the round and any incomplete component makes its row inconclusive.

## Boundary and table

The degree boundary is round 4's unbranched degree 4.  The result table is:

| terminal | branch depth | components/target | targets | positive | negative | correct | complete | max solving degree | median solving degree | max columns | total field ops | ops / unbranched | classification |
|---:|---:|---:|---:|---:|---:|:--:|:--:|---:|---:|---:|---:|---:|:--|

Total field operations include every component for all targets.  The fixed
unbranched denominators are 55,498,248 operations for terminal zero and
52,441,153 for terminal 369.

## Stop conditions and scope

Stop after both branch depths and terminals.  Preserve all failures and raw
component records.  Commit native Rust output, commands, and hashes.

This is the final small-prime branching rung in this sequence.  Even a success
would require a separately preregistered scaling/cost study before any P-256
or ECDLP claim.  No P-256 Gröbner basis is run here.
