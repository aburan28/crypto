# Full-cost rho parity: frozen engineering iterations

Reference: `a2963f0d6fb5cb74a20e678cf4def0c7cbafa1d9` (PR347 including PR348).
This experiment pursues parity; it does not assume parity is attainable.

## Boundaries and acceptance, frozen before production changes

Retain the ten cases, every factor base, two-summand decomposition, 256
relation trials, batch size one, node budget 20,000, direct-relation shortcut
disabled, and cache-off conditions from the Weil composition contract. The
counting boundary remains at most `B(B+1)/2` unordered sums from `B` points,
and hence uniform-target yield at most `min(1,B(B+1)/(2(N-1)))`. Rho is the
existing signed-Frobenius implementation on exactly the same curve, subgroup,
target and seed. Cold cost includes curve/base/plan setup, target generation,
all failed attempts, extraction, verification and relation-matrix work.
The prior cold chart/rho ratios were approximately 3–17 on useful toy cases.
The two zero-column controls remain reported and cannot count as successes.

Use all four prior seed/log pairs `(101,5),(202,5),(503,7),(607,11)` as frozen
inputs, plus fresh pairs `(1201,13),(1429,17),(1601,19),(1873,23)`. Three paired
repetitions in alternating order, one CPU, one thread, 2 GiB address-space
limit, 10 second per-process limit. Keep the chart baseline, enumeration,
pair table, and rho. Independently validate every oracle call of every
completed full DLP and final `[k]G=Q`, including all failures and controls.
Run the seven native S3/S4 regression cases on both frozen and fresh inputs,
all three repetitions, with reference and candidate builds. The external
WDSat S4 binary does not call this native encoding; the existing equivalent
native suite remains the regression gate. Compare complete chart root sets
on stage seeds 17, 937, 1201 and 1873, eight targets each.

An iteration passes the engineering runtime gate only with unchanged verified
workloads and a fresh paired 95% interval entirely below 1. Its numerical
target is at least 20% lower cold runtime. Rho runtime parity requires the
upper paired 95% interval of candidate/rho to be at most 1 on **every useful
case**, with every target solved; an aggregate cannot hide a losing case.
This is a toy-suite practicality gate, not degree-131 parity or an exponent
claim. Total operations, S and cost/floor ratios stay null unless the whole
pipeline is measured in a documented common calibrated unit. No isolated
kernel timing or partial counter is promoted to an attack advance.

## Iteration 1: single-word arithmetic and existing point lifts

Hypothesis: repeated generic arithmetic consumes avoidable cold and per-trial
cost. Use the existing exact single-word field/curve operations where valid,
compute abscissa lifts by the same Artin–Schreier root convention, and retain
point ordering. Reuse the factor-base points when testing algebraic roots.
Validate cached point content and ordering as well as curve/domain metadata;
the driver validates an immutable plan once per run. Keep the F4 algorithm,
complete component cover, relation draws and witness order unchanged.
Avoid launching Rayon for fewer than 256 cofactor projections. For a nonempty
negation-closed factor base, even arity passes cofactor admission by the exact
witness `P+(-P)+...=O`; no class enumeration is necessary. This necessary
condition says nothing about the decomposition yield of a prescribed target.

Abandon or correct on any root loss, false refutation, mismatched point
indices, changed input draws, invalid scalar, skipped cost or dropped failure.
Retain the measured comparison even if the hypothesis fails.
