# Result 1: nonlinear factor-base label gauges

The exact 30-class census is a **negative degree result and a useful
relabelling control**.  No nonlinear label class refuted either discovery
system through degree 6.  The identity encoding refuted both at degree 6, as
the frozen reference requires.  The selected nonlinear class made the visible
matrices much smaller but raised the input degree from 3 to 5 and still did not
produce a refutation at degree 6.  That is work moved outside the measured
window, not a solving improvement.

## One table, one unit

The unit is the exact bounded-Macaulay outcome and column count on one frozen
system.  `—` means no refutation through degree 6, not a lower degree.

| split / draw | variant | input degree | degree-5 columns | degree-6 columns | refutation degree through 6 | correct |
|---|---|---:|---:|---:|---:|---|
| discovery / 3 | identity reference | 3 | 6,513 | 14,642 | 6 | yes |
| discovery / 3 | nonlinear winner `[0,1,2,3,4,5,7,6]` | 5 | 3,577 | 9,893 | — | yes |
| discovery / 4 | identity reference | 3 | 6,276 | 14,396 | 6 | yes |
| discovery / 4 | nonlinear winner | 5 | 2,077 | 6,739 | — | yes |
| holdout / 9 | identity reference | 3 | 6,274 | 14,396 | 6 | yes |
| holdout / 9 | nonlinear winner | 5 | 2,041 | 6,690 | — | yes |
| holdout / 11 | identity reference | 3 | 6,276 | 14,396 | 6 | yes |
| holdout / 11 | nonlinear winner | 5 | 2,077 | 6,739 | — | yes |

The winner was serialized before either holdout was evaluated.  Its degree-6
column counts are 53.1% and 53.2% below the identity on the holdouts, but the
identity has already refuted the systems at that degree while the nonlinear
encoding has not.  Counting the smaller unresolved matrices as a gain would
be exactly the repository's **relabelling** failure mode.

## Search and correctness

The program constructed all 1,344 elements of `AGL(3,2)`, partitioned all
40,320 label permutations into exactly 30 right cosets, and evaluated the
lexicographically least representative of every class.  All 29 non-affine
representatives had quadratic coordinate ANFs.  Every representative was
checked on every one of the 65,536 Boolean assignments of both discovery
systems.  Baseline and winner received the same exhaustive check on both
holdouts.  All checks commuted with the original equations and every frozen
system retained solution count zero.

Run 01 is preserved as an excluded resource failure: the default matrix caps
stopped the known baseline before degree 6.  It exposed no holdout result.
Run 02 used the registered source with explicit limits
`F4_F2_MAX_ROWS=2000000`, `F4_F2_MAX_COLS=200000`, and
`RAYON_NUM_THREADS=1`; it completed with exit status zero.  The isolation
record reports 0.03 other CPU-seconds during 23.96 wall-seconds, but Docker
threads could not be moved off the reserved CPU.  Wall time is therefore only
a practicality note and supports no speed claim.

## Verdict and scope

- Primary gate (both holdouts refute by degree 5): **fail**.
- Secondary gate (same degree 6, fewer columns on both holdouts): **fail**,
  because the nonlinear winner does not refute at degree 6.
- Classification: **relabelling / negative stage diagnostic**.

The finite-cube idea appears distinct in the focused literature and repository
search recorded in the protocol, but novelty is not proven.  More importantly,
the bounded experiment falsifies this version of it.  This is `GF(2^7)` data,
not the required `m = 83` confidence gate and not evidence at `GF(2^131)`.
Full IC cost, `S`, calibrated operations and the rho ratio remain null.

## Result 2: affine follow-up

The exact affine tuner is also negative.  All 1,344 gauges produced exactly
the same summed discovery matrix dimensions at degree 5: 12,789 columns and
13,496 rows.  Input support did move—summed generator terms ranged from 950 to
1,221—so the census was capable of distinguishing presentations.  The 12
preregistered finalists all produced the same summed degree-6 dimensions:
29,038 columns and 49,364 rows, and all refuted at degree 6.

The frozen winner was `[0,4,3,7,1,5,2,6]`.  On the holdouts its matrix sizes
were identical to identity:

| split / draw | variant | input terms | degree-6 rows | degree-6 columns | refutation degree | correct |
|---|---|---:|---:|---:|---:|---|
| holdout / 9 | identity | 565 | 22,498 | 14,396 | 6 | yes |
| holdout / 9 | affine winner | 520 | 22,498 | 14,396 | 6 | yes |
| holdout / 11 | identity | 489 | 22,498 | 14,396 | 6 | yes |
| holdout / 11 | affine winner | 524 | 22,498 | 14,396 | 6 | yes |

The term-count change is mixed: lower on draw 9, higher on draw 11.  It did not
move the registered exact-work metrics.  The 5% gate and the weaker strict-
reduction gate both fail.  Classification: **accounting / negative stage
diagnostic**.  The result also gives a useful invariant for this frozen panel:
affine gauges move generator sparsity but not the actual Macaulay dimensions
observed at degrees 5 and 6.

Run 03 completed with exit status zero after serializing the winner before the
holdouts.  The isolation record reports 0.46 other CPU-seconds during 495.75
wall-seconds and no contention, but Docker threads remained on the reserved
CPU; the time is not used as evidence.  Every one of the 1,344 discovery
gauges and both holdout arms passed exhaustive assignment equivalence and kept
solution count zero.
