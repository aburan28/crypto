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

