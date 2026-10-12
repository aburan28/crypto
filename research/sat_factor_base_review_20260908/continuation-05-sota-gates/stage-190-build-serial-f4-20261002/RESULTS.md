# Stage 190: serial F4 build phases with parallel elimination

## Decision

`SELECTED_FOR_REPOSITORY_DEFAULT`. The selected current F4 path now serializes
monomial-product construction, symbolic preprocessing, frontier sorting, and
row packing inside each F4 call. Inner five-column `BlockTables` construction
and row reduction remain parallel, and all 242 independent fixed-X1 systems
remain in the twelve-worker outer batch.

Unset and `F4_F2_DISABLE_INNER_BUILD_PARALLEL=1` select the new path.
`F4_F2_DISABLE_INNER_BUILD_PARALLEL=0` retains the prior nested-build control.

## Frozen screen and confirmation

The one-pair screen passed at build-serial/current ratios `0.638598` wall,
`0.927132` total core, and `1.302465` RSS. The frozen three-pair confirmation
then produced:

| pair | order | wall | total core | peak RSS |
|---:|:---|---:|---:|---:|
| 1 | current, build-serial | 0.908681 | 0.927715 | 0.847874 |
| 2 | build-serial, current | 0.913292 | 0.910502 | 1.026266 |
| 3 | current, build-serial | 0.948125 | 0.926973 | 1.564208 |
| **median** | frozen interleaving | **0.913292** | **0.926973** | **1.026266** |

The candidate clears the pre-registered strict-below-`0.97` wall and total-core
gates. Median wall falls 8.67 percent and CPU falls 7.30 percent. Median RSS
rises 2.63 percent and is variable across pairs; this is reported rather than
hidden or converted into a post-hoc gate.

The exact selected-commit unset replay contains no scheduling environment
override and takes `17.519346` wall-seconds, `181.812510` total core-seconds,
and `4,176,412,672` bytes peak RSS. Its native seven-check replay authenticates
the meter artifacts, proves the unset route, and exactly matches the selected
confirmation structure.

## Correctness

Every screen, confirmation, and default-replay record:

- authenticates public source instance
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`;
- uses algebraic factor base `span_F2(1,z,...,z^8)` without target-subgroup
  enumeration or known factor-base log labels;
- reproduces equation fingerprint
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`;
- visits all 512 masks, skips 270 non-rational masks, constructs and completes
  all 242 systems, finds zero roots, and returns exhaustive `UNSAT`; and
- agrees exactly on equations, terms, pair and field-pair counts, dense-pair
  counts, divisor tests, matrices, basis, extraction, degree, logical XORs,
  performed XORs, table/matrix/scratch memory, and full-M4RI counters.

The phase screen passes 8/8 checks, confirmation passes 16/16, default replay
passes 7/7, and final composition passes 29/29. Final result SHA-256 is
`ea84b500f977acd66b5a36b4632620c357cc2b301a5c9987be39a096b0d24384`.

## Accounting and boundary

Stage 190 contributes a measured lower bound of 24 components,
`485.938466` wall-seconds, `3,027.193941` total core-seconds, and
`5,897,338,880` bytes peak RSS. The cumulative measured campaign lower bound is
564 components, `23,476.882993` wall-seconds, `59,943.278380` core-seconds,
and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null` because development compilation before
the exact native builds and the final writes after the charged composition
control are not all outer-metered. No modelled values replace them.

This is a selected implementation improvement on one opened decomposition
target. Same-binary direct MITM remains decisively faster, full-cost
automorphism-aware rho gate 6 remains false, and independent reproduction and
novelty review remain absent. It is not relation-yield evidence, an unknown-
scalar run, a full-DLP result, or Koblitz index-calculus SOTA.
