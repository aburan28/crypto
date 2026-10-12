# Charge Boolean construction and complete linear reduction together

The preceding turn made verified progress: it removed the exponential lookup,
tested through 36 variables and established the density-dependent cost of sparse
intermediates. Construction gains alone left a missing measurement: what happens
when all forward and backward linear reduction is charged? This experiment
measures that combined task with fresh holdouts and an exact canonical output.

This is **not a complete polynomial solver**. It constructs one bounded-degree
product matrix and returns its reduced row-echelon form (RREF). It does not
iterate Gröbner closure to completion, find roots at large n, collect elliptic
curve relations or recover scalars. Full-solving and cryptanalytic costs remain
unknown. The broad speedup goal cannot be closed by this narrower measurement.

## Exact contract and streaming proof

Return the original compact column labels in ascending monomial order and the
unique nonzero RREF rows of the complete capped product matrix. The source-row
budget counts **every nonzero generated row before elimination**, including
dependent rows. Reaching full rank or obtaining a constant polynomial never
short-circuits construction or changes resource accounting.

The common reducer inserts each row into a pivot-indexed basis, XORing away its
leftmost pivot until the row is zero or gets a new pivot. A final backward pass
clears pivots above their rows, producing canonical RREF. In the initial comparison, all staged arms use
this exact reducer on their compact matrices. The streaming arm instead builds
each product in ambient ranked coordinates and inserts it immediately, retaining
the original support before elimination. It finishes the same backward pass and
then compacts the output. Context changes rebuild and reduce directly.

Every insertion preserves the span through row XORs or zero-row removal. Sorting
the final pivots and backward elimination give RREF. Ambient columns absent from
all source rows are zero throughout elimination; deleting them is an order-
preserving projection and commutes with the row operations. Thus the streamed
answer must equal the staged RREF, not merely have the same rank.

The independent oracle uses recursive/set-parity product construction and a
column-oriented Gauss-Jordan algorithm that clears both above and below each
selected pivot immediately. Every measured output equals that oracle. Tests
also check explicit row spans of all 5,050 binary matrices with 1–3 rows and
1–4 columns, all 6,144 small Boolean coefficient/degree/mask cases, larger family
fixtures, dependent-row caps, word boundaries and context guards. A separate
test enumerates all eight Boolean points for 400 small polynomial systems to
check that the original and reduced equations have the same zero set.

## Frozen comparison

`protocol.json` fixes n=12/20/28/36, three unscreened fixture families, batches
1/4/8, discovery seeds 17/937, new holdouts 20260929/196613, and twenty balanced
repetitions. The families are quadratic coefficients, linear degree drops and
restricted-mask quadratic/linear/constant/zero cycles. There are 144 cells.

The staged controls are sorted construction, ranked packed construction and
sparse-intermediate construction, each followed by the common reducer. A dense
lookup control is included at n=12 only. The candidate streams ranked products
directly into that reducer. Streaming may lose: it eliminates using a wider
ambient coordinate system than a precompacted matrix.

The **dramatic combined-workload gate** requires a 95% paired-bootstrap lower
bound above **2.0** for the pointwise fastest staged arm / streamed total on
every family and each of n=20/28/36 at batch8. All nine comparisons must pass.
Each optimized arm is also reported against sorted construction. No seed,
family or failing gate may be dropped after observation, and no tuning uses
these holdouts. Intervals describe fixed seeds and repeated measurements on
this host, not a population or independent reproduction.

Cold totals include setup, construction, complete reduction, compaction, output
allocation, exact validation and destruction. Staged construction and reduction
have separate timers. The streaming arm has one fused timer: its unobservable
separate components are **null**, not zero. Fixture/oracle generation is common
supplied-input work outside arm totals and inside fresh-worker receipts.

Logical row XORs and executed 64-bit XORs are counted in reduction, including
the backward pass. Logical counts must agree across matched arms; physical
word counts can change with coordinate width. They are reduction diagnostics,
not a complete calibrated operation unit. Retained context capacities and the
maximum end-of-insertion basis storage are separate measurements. The latter
is not total peak memory; incoming rows and final-compaction temporaries are
excluded. Process RSS includes all arms and common reference matrices.

The reference boundary is the fastest correct staged pipeline on each pair.
Canonical output still must be materialized. A diagnostic construction-only
projection sets the measured construction phase to zero while holding all
others fixed: `total / (total - construction)`. This is a conditional projection,
not a bound on implementations that change multiple phases or on full solving.

## Execution and artifacts

```sh
python3 research/boolean_construction_reduction_20260922/run.py --out research/boolean_construction_reduction_20260922/run_01
python3 -m unittest discover -s research/boolean_construction_reduction_20260922 -p 'test_*.py'
```

Every run uses a new directory with frozen sources/protocol, executable hashes,
raw fixtures and measurements, process receipts, test output and a manifest.
Replay uses temporary copies. The previous construction studies are untouched.
No production solver path changes; the WDSat/full IC suite is inapplicable to
this standalone generic matrix task. The canonical scoreboard records the
combined cost without upgrading it to an end-to-end polynomial-solving claim.

## Additive nonzero-word and incidence experiment

The first run rejected streaming in all nine comparisons with the fastest staged
path. Its wider coordinates increased physical word XORs without changing the
logical row eliminations. The frozen run is retained unchanged.

`incidence_protocol.json` fixes a successor with new holdouts 20260930/262147 and
thirty balanced repetitions of five arms, or six at the n=12 bridge. It keeps
all prior pipelines as fresh controls and adds `incidence_reduce`: sparse
construction followed by a reducer that caches nonzero pivot-word positions
and precomputes the target rows for backward elimination.

Forward pivot rows stay immutable until backward elimination, so their cached
nonzero-word lists are exact. The forward-echelon pivot pattern determines all
backward targets: removing a higher pivot cannot change a lower pivot column,
because every higher-pivot row is zero in that column. Each backward pivot's
word list is refreshed after all higher pivots have been eliminated. Complete
canonical RREF and the logical row-XOR count must match the original reducer.

The successor retains the same strict **2x against the fastest prior pipeline**
gate on all nine larger-size/family cells. Cache construction, pivot masks,
incidence lists and list refreshes are timed inside reduction. Auxiliary storage
is recorded separately from end-of-insertion basis storage; neither number is
misrepresented as full peak memory. The current source uses `with-incidence`
only when the runner is given the successor protocol.

```sh
python3 research/boolean_construction_reduction_20260922/run.py --protocol research/boolean_construction_reduction_20260922/incidence_protocol.json --out research/boolean_construction_reduction_20260922/run_02
```
