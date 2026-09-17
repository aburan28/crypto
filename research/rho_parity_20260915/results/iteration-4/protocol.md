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

## Iteration 2: packed residual solving after linear projection

Frozen after iteration 1: its fresh cold chart/reference ratios were
0.179–0.756, but chart/rho remained 1.70–16.89. The complete root sets and
every matched relation draw/witness were unchanged.

Hypothesis: for the small residual affine spaces left by row projection,
constructing polynomials and invoking F4 costs more than exact enumeration.
When at most twelve coordinates remain, enumerate that affine space and
check the original quadratic coefficient rows by packed parity. Otherwise
retain F4. Preserve every component pair, every algebraic root and point
verification. Charge each residual assignment against the same whole-cover
node budget as F4 reductions, with a separate assignment counter; these
different work items are not a common cost unit. A budget stop is unknown.
The root visitation order can change, so independently cross-check all
oracle inputs and compare identical random-draw prefixes and final logs.

Numerical gate: fresh cold runtime at least 20% below iteration 1, with the
same complete verified workloads and the same all-cases rho parity gate.
Keep all original controls and the entire matched suite. This hypothesis
does not enlarge the factor base, query a precomputed pair table or make
an asymptotic claim.

## Iteration 3: rational support and packed factor-base construction

Reference `c71d4a8c705dfc597017e476a93312a109f80e86`. Iteration 2's fresh
chart/rho ratios were 1.35–4.21; its 2,004-process comparison passed every
check. Independent `support.py` census confirms that the difficult ordinary
degree-eleven base has rational abscissae in a **one-dimensional** subspace
of each nominal three-dimensional component. Most other useful families
retain their nominal span, so no universal rank reduction is assumed.

For each component, take the linear span of coordinate assignments whose
abscissae occur in the factor base. Its annihilator gives homogeneous linear
constraints that every valid point decomposition must satisfy. Add these to
the group-witness search before projection. The public algebraic root-set
visitor keeps the original domain and every S3 root; the point-witness search
can exclude nonlifting abscissae. This is reusable, target-independent curve
preprocessing, with no pair-sum table and no target-specific cached answers.

Also build Frobenius unions and orbit maps with one packed field context,
and reuse the already materialized factor-base lifts directly in point-array
order. A reordered base can choose a different valid witness; stale cached
indices are still rejected. Freeze the same 20% cold-runtime engineering
target and all-cases rho gates, full matched suite and correctness controls.
Classify this as engineering; the counting boundary and factor base do not
change. No degree-131 extrapolation follows from the observed support ranks.

When the support span has smaller positive dimension, also construct its
Frobenius closure as a compact internal point-search chart. Accept this
compression only when the new cover contains every original factor-base
point with unchanged indices. Its preprocessing is charged. Retain the
original charts for unrestricted algebraic queries; the full point base and
counting bound are unchanged. This avoids carrying nominal coordinates that
can never occur in a rational witness.

## Iteration 4: reuse relation preprocessing

Iteration 3 cut the degree-eleven ordinary case by another 57%; its fresh
ratio to rho was 1.81. The full suite and the exhaustive support-compression
test passed. Parity is still absent; the degree-thirteen case remains the
largest fresh runtime ratio (3.26 on that shared-host run).

A detailed instruction profile exposes unneeded `FieldStructure::new`,
repeated modular exponentiation of the same Frobenius eigenvalue, and scalar
verification through general arithmetic. Freeze a fourth engineering change:
create ambient field encoding only for a strategy that uses it, compute
`1,lambda,...,lambda^(n-1)` once by recurrence, use the existing packed curve
for production verification, and precompute binary powers of `G` and `Q` for
the unchanged independent `(a,b)` draws. Charge every precomputed point and
coefficient to the cold run. Keep external verification in general arithmetic.
The expected improvement is >=20% in fresh cold runtime and lower complete
instruction count, with byte-identical attempt records and relation outcomes.
The existing per-case rho gates and all controls remain unchanged.
