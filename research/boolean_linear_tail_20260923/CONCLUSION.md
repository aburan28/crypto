# Large-active linear-tail gains are measured; promotion remains guarded

The latest retained run measures **3.59x** lower cold cost on the declared
large-size workload mixture against the pointwise faster of the two retained
kernel strategies, with a 95% paired-bootstrap interval **[3.57, 3.62]**. It is
**1.64x** faster than the pointwise best of the four prior experimental arms,
interval **[1.62, 1.66]**. The identical-control timing gate passes.

**The guarded portfolio verdict is still REJECTED.** Four small-case
non-regression intervals fail the predeclared 0.95 lower-bound condition. Their
paired median reference/candidate ratios are 0.969, 0.975, 0.969 and 0.964. The
universal per-cell 2x gate also remains rejected: the hybrid deliberately reuses
the flat control on small active spaces and does not double its speed there.
No threshold or failed case has been removed to turn these measurements into
an accepted default replacement.

These are conditional **specialized linear-tail kernel** results. They do not
measure complete polynomial solving, recursive search or index calculus. The
broader thread goal remains open.

## What produces the gain

The source-matched contract is the canonical basis of W intersect L: all linear
and constant vectors in the complete capped degree-D product row space. The
controls retain the existing cached/fused flat method and minimum-weight sparse
bucket strategy; linear-tail restriction and sparse high-column elimination
are not claimed as new ideas.

The candidate streams degree-ordered nonlinear monomial keys while carrying
the affine component in a single word. Lighter same-leading pivots may replace
earlier pivots, with the old XOR new residual retained for further elimination.
An exact binomial-index bitmap counts original columns without a per-term
support hash set or a global matrix layout. It preserves exact source dimensions
and all row/column-cap semantics. A fixed density rule retains the identical
flat path when the degree-bounded active monomial universe has at most 512
columns. The later integer cutoffs are proven equivalent to that same rule.

In the latest run, every large-active quadratic, linear-drop and cross-cancel
cell at nominal n=20/28/36 passes the 2x gate. Paired median gains against the
better retained control range **3.17–4.45x**; the smallest lower confidence
bound is **3.13x**. The restricted-cycle family has eight occurring variables
and correctly takes the flat path. Its remaining dispatch cost matters to the
non-regression criterion.

## Latest cold costs

The milliseconds below are copied from `run_09/RESULT.md`: one cold batch of
eight inputs at nominal n=36, including schedule/layout work, product parity,
exact source accounting, high elimination, affine canonicalization, validation
and destruction. Values are medians over two holdout seeds and 42 observations.
Each observation contains sixteen independently reset cold batches; the table
normalizes by sixteen. Cache state is never retained between those batches.
Ratios here are descriptive ratios of pooled medians against the flat control;
the acceptance gates use paired ratios against both retained controls.
All measured rows have exact oracle equality PASS and class engineering.

| Method | Quadratic (ms) | Linear drop (ms) | Eight-active cycle (ms) | Cross cancellation (ms) | Flat / method, quadratic |
|---|---:|---:|---:|---:|---:|
| Cached flat tail | 1.481716 | 2.955561 | 0.074004 | 1.590747 | 1.000 |
| Existing sparse bucket | 0.943720 | 1.195530 | 0.591923 | 0.943421 | 1.570 |
| Plain sparse streaming | 0.441522 | 0.564499 | 0.579945 | 0.435589 | 3.356 |
| Pivot-exchange streaming | 0.437848 | 0.565529 | 0.466715 | 0.423594 | 3.384 |
| Exact bitmap census | 0.213036 | 0.302470 | 0.444639 | 0.226766 | 6.955 |
| **Fixed hybrid dispatcher** | **0.213638** | **0.301923** | **0.082290** | **0.223602** | **6.936** |
| Identical flat alias | 1.484120 | 2.994671 | 0.076327 | 1.607216 | 0.998 |

The mixture gives each of the twelve n>=20/family cells one cold batch per
paired seed/repetition. These are declared benchmark weights, not measured
solver-call frequencies. An overall solver claim would require actual call
traces and the rest of the solver's work.

## Experimental history and measurement limits

All complete runs and the failed launch remain intact. `RUN_LEDGER.json` binds
their protocols, workers, manifests and verdicts. The sequence is:

| Run | Change or finding | Status relevant to promotion |
|---|---|---|
| 01 | Inherited numeric term/multiplier order | Correct row spaces; not source-order timing evidence |
| 02 | Source term order and degree-layered multipliers, fresh holdouts | Exchange passes 7/12 universal gates |
| 03 | Duplicate launcher mode argument | Producer failure before any measurements; hash-bound failure retained |
| 04 | Exact bitmap census and fixed density dispatch | Aggregate ratios pass; small-path guards fail |
| 05 | Remove recursive dispatch, fresh holdouts | Aggregate ratios pass; small-path guards still fail |
| 06 | Balance predecessor labels and add an identical flat alias | Alias and guard gates fail |
| 07 | Typed timed kernels and predecessor/cycle pairing | Alias and guard gates still fail |
| 08 | Sixteen cold batches per observation, identical non-inlined flat leaf | Both aggregate ratios and all alias checks pass; one guard fails |
| 09 | Proven-equivalent integer dispatch cutoffs, fresh holdouts | Both aggregate ratios and all alias checks pass; four guards fail |

The earlier cyclic schedules balanced position but preserved neighbors. The
later schedule uses Euler circuits so every method follows every method equally
often, including itself. One fully validated unmeasured primer observation
establishes the first predecessor; its cost is in process receipts. Each paired
unit has the same predecessor and cycle. This addresses first-order carryover,
not every possible hardware-history effect. Typed labels enter the same flat
leaf; label parsing is harness work outside timing, while actual density
selection remains timed. Earlier timings are not silently replaced.

Run 08 measured a mixture ratio of 3.66x [3.64, 3.71]; run 09 measured
3.59x [3.57, 3.62]. They are fresh-seed runs on the same host, with a dispatcher
arithmetic change between them, not unaffiliated reproduction or an unchanged-
binary confirmation. Both retain their failed guards. In run 09 those lower
bounds are 0.94543 (n12 quadratic), 0.94945 (n12 restricted cycle), 0.93995
(n20 restricted cycle) and 0.92144 (n36 restricted cycle), below the fixed 0.95.

## Correctness, costs and provenance

Eight complete grids contain 1,536 cells, 271,872 measured observations and
**8,516,352 checked measured affine-tail outputs**. These output counts include
repeated cold batches; they are not counts of distinct random systems. The
failed launch contributes zero measurements. Its worker, protocol and analyzer
match the recovered run exactly; only the launcher argument construction was
fixed. `LAUNCH_FAILURE.json` retains the failure evidence.

Thirteen current Rust tests cover ordering, exact bitmap cardinalities, dispatch
equivalence, canonical affine intersections, 24,456 small matrix/split cases by
explicit span enumeration, 6,144 Boolean coefficient/degree/mask cases, larger
families, hidden affine consequences, source caps, cache guards, schedule
balance and alias identity. The evidence tests independently reconstruct exact
source support, cache hits, dispatch counts, source dimensions and cold-batch
normalization, and replay the frozen results. Every complete manifest has 15
files; all 120 hashes are retained, along with actual executable hashes and
process receipts.

All source generation, failed cache guards, index preparation, sparse merges,
low-block reduction, exact validation and destruction are charged. Reference
preparation and the timing primer are supplied-input/measurement work outside
arm totals and inside process receipts. Retained schedule masks, layout columns
and census entries are counts, not byte estimates. Packed XOR counts and sparse
merge counts are different diagnostics, not interchangeable total operations.
The latest whole-worker peak RSS is 7,634,944 bytes, including all methods and
reference data; it is not candidate-specific peak memory.

No production solver path was modified or executed. No curve input, scalar
recovery or cryptanalytic comparison occurs. Full polynomial-solving, full IC
and rho costs remain null. The next relevant decision requires complete generic
solver traces and their real density/call distribution, together with the
remaining small-case guard issue; a default promotion is not justified here.
