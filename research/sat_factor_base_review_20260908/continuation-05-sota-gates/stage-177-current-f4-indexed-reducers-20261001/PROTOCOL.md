# Stage 177 protocol: exact indexed reducers in current F4

## Hypothesis

The current repository F4 performs `4,190,633,182` active-leading-monomial
divisibility tests on the frozen single target. Each symbolic monomial normally
has low degree while the active basis has thousands of elements. Enumerating
the nonzero submasks of the monomial and probing an exact leading-monomial map
should return the identical shortest reducer with far fewer lookups.

This ports only the adaptive exact-submask lookup into the current repository
F4. It retains the current pair queue, bitmap/hash symbolic sets, adaptive
`BlockTables`, tiled elimination, and current correctness fixes. It does not
restore the older Phase B F4 fork.

## Algorithm and control

For each symbolic monomial `m`:

1. Build one exact map from active leading monomial to the active polynomial
   with the fewest terms, breaking ties by the existing lowest index.
2. If `2^deg(m) - 1 <= active.len()`, enumerate every nonzero submask of `m`
   and probe the map.
3. Otherwise scan the active list exactly as the current engine does.
4. Among all divisors found, select the fewest-term polynomial and then the
   lowest index. Thus the selected reducer is byte-for-byte identical to the
   linear reference.

The candidate is default-on. `F4_F2_INDEXED_REDUCERS=0` selects the same-binary
linear reference. Lookup counts are exported as
`divisor_submask_lookups` and `divisor_linear_tests`; their sum must equal
`divisor_tests`.

## Correctness gates

- An exhaustive differential unit test compares indexed and forced-linear
  lookup over varied active bases and every monomial in their small domains.
- The existing eight Boolean-F4 tests pass, including agreement with certified
  Buchberger bases and exact `BlockTables`/row-by-row equivalence.
- The fixed-X1 specialization identity test and three backend tests pass.
- Candidate and control reproduce the frozen equation fingerprint
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Both return exhaustive UNSAT, visit 512 masks, complete 242 systems, and
  retain exact F4 matrix, basis, row-equivalent XOR, actually-performed XOR,
  pair, and extraction counts. Only divisor accounting and timing may differ.

## Frozen benchmark

- Same opened `n=59, ell=9, m=3` target and algebraic factor base as Stages
  175–176; no subgroup enumeration and no discrete-log labels.
- `RAYON_NUM_THREADS=12`, `PQ_F4_X1_BATCH=512`, all BLAS/OpenMP controls one.
- Three interleaved pairs in fixed order:
  `linear, indexed, indexed, linear, linear, indexed`.
- Each process has the same 300-second internal and 360-second watchdog budget.
- Fresh exact-commit release build is charged separately.

## Decision rule

Select the indexed default only if every correctness gate passes and the median
paired indexed/linear ratios are both below `0.97` for wall and total
core-seconds. Peak RSS is reported and cannot be omitted. Also report absolute
ratios to the Stage 174 medians and same-binary direct MITM; a lookup-stage win
is not a full-method speedup.

All build, candidate, control, failed, and validation processes are charged.
Single-core time remains null in this twelve-worker experiment. This one-target
engineering result cannot change any SOTA gate by itself.
