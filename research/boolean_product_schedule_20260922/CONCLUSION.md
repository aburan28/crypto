# Decision: reject exact-support schedules as a new cross-instance optimization

The experiment requested in the previous iteration is implemented and has run.
It is a standalone Boolean matrix-construction experiment; the production F4
solver is unchanged. `run_01` and `run_02` retain the fixed protocols, exact sources,
compiler and executable identities, raw worker output, process resource
receipts, complete-grid verifier and manifest.

## Mathematical result

For canonical squarefree polynomials over F2, support equality is polynomial
equality. With the same ordered generators, variable count, multiplier mask
and degree bound, an exact full-support key specifies exactly the same matrix.
Thus there is no cross-instance symbolic reuse under the proposed key. Changed
constants, removed terms, and degree drops all require a fallback. This is a
property of the key, independent of the implementation's timing.

The measured candidate stores each surviving product row as a list of column
indices. The stronger control stores the completed packed matrix and uses
exactly the same signature. Both pay initial construction, key checks,
fallbacks, output allocation, validation and destruction in each cold batch.
Neither populates additional cache entries after a miss.

## Fixed-run observations

All four Rust correctness tests pass, including exhaustive supports through
three variables, evaluation checks, signature changes and changing resource
caps. The run completed all **256 cells**, **8,192 batch-arm samples** and
**174,080 independently verified matrix outputs**. Exact schedule hits on
changed canonical systems: **zero**.

Initial run (`run_01`), retained as the before measurement. One timing unit:
cold batch construction and verification in
milliseconds for 64 identical-input repetitions. Each cell is the median
across two holdout seeds and eight balanced-order repetitions. These are
matrix-construction diagnostics, not full solver measurements.

| Variant | 6 variables | 8 variables | 10 variables | 12 variables |
|---|---:|---:|---:|---:|
| Direct construction | 0.556230 | 1.138104 | 1.953042 | 3.113875 |
| Verified layout reuse | 0.485042 | 0.988375 | 1.710438 | 2.754792 |
| Exact product schedule | 0.084417 | 0.148959 | 0.243229 | 0.367833 |
| Exact packed-matrix cache | 0.070375 | 0.122396 | 0.201166 | 0.308166 |

The preregistered numerical gate was a 95% paired-bootstrap lower bound greater
than 1.05 for control total time divided by schedule total time, at both batch
16 and batch 64 for every size, against **both** layout reuse and matrix caching.
The schedule passed all eight layout comparisons and failed all eight packed
matrix comparisons. All eight upper bounds for the packed-matrix ratio are
below 1. These intervals describe repeated measurements of these fixed
fixtures, not performance across an unknown workload population.

In the initial repeated-input holdout cells, the packed-matrix cache retained less
storage at every size. At 12 variables the median retained structures were
17,546 bytes for the schedule versus 11,324 bytes for the matrix cache. These
counts include Rust structure and allocation capacities, exclude allocator
metadata, and are distinct from fresh-worker RSS and output payload size.

## Compact-storage amendment, retained as run_02

Review identified spare column-vector capacity in the schedule that the
packed-matrix cache did not retain. The amendment compacts the schedule's
column storage and charges that operation to compilation. It changes no
fixture, seed, acceptance threshold, cap or repetition. The initial evidence
is preserved. All four correctness tests and the same 256-cell grid completed
again, with 174,080 verified outputs and zero changed-input schedule hits.

| Variant | 6 variables | 8 variables | 10 variables | 12 variables |
|---|---:|---:|---:|---:|
| Direct construction | 0.571125 | 1.122459 | 1.995958 | 8.434125 |
| Verified layout reuse | 0.502396 | 0.968792 | 2.798437 | 2.801292 |
| Exact product schedule, compact columns | 0.082334 | 0.145250 | 0.244792 | 0.385563 |
| Exact packed-matrix cache | 0.070646 | 0.117395 | 0.198708 | 0.310896 |

Unit and sample selection are unchanged: milliseconds per cold batch of 64,
on the held-out repeat fixtures. Runtime varies across runs; retain the paired
within-cell comparisons rather than interpreting cross-run differences as a
gain. The schedule again passes all eight layout comparisons and fails all
eight packed-matrix comparisons, with every packed-matrix ratio upper bound
below 1. No further timing rerun was selected to improve those results.

Compaction reveals a space/time tradeoff at the largest toy size: 9,942 retained
bytes for the schedule versus 11,324 for the packed matrix (12.2% less), while
schedule replay remains slower. At the other three sizes the packed matrix
still retains less. This replaces the initial memory ranking as the current
assessment without erasing its measurement.

## Disposition

**Reject performance promotion and cross-instance reuse under the exact key.**
The apparent improvement over rebuilding products disappears as a new idea
when compared with the simpler exact matrix cache. The compact schedule's
largest-size memory saving is a bounded space/time tradeoff, not a runtime win.
Changed-input families
retain the negative controls: every changed system falls back in both exact
cache arms, while the layout-only arm can sometimes reuse a column layout
despite changing generator supports.

A parameterized support-envelope schedule would be a different experiment:
it must represent coefficients that can change, preserve parity under
colliding Boolean products, add multipliers when a generator's degree falls,
and rebuild when products leave its admissible support. None of those claims
is established by an exact full-support hit.

This finding concerns an algebraic cache contract. It supplies no ECDLP
exponent, factor-base improvement, relation-yield measurement or rho ratio.
The previous solver timings remain historical stage measurements.
