# Batched prime pair-table construction: frozen comparison

## Question and prediction

For the two-summand negation-folded pair table at 64 prime abscissae,
can one Montgomery batch inversion per table row replace the individual
affine inversions while retaining the same table keys, relations, verified
answers and counted group-addition equivalents (GAE)? The table stores
only the x coordinate of each sum; the candidate computes that coordinate
directly and leaves the probe's exact group re-addition in place.

Prediction: table-build field inversions fall from roughly one per pair
to one per base abscissa. The group-addition ledger remains unchanged.
The primary correctness gate is identical output, relation counts, table
entry counts, and GAE on paired executions. The native field-work counters
and table-build wall observations are secondary diagnostics. No wall-time
advantage is claimed from this macOS L0 run.

## Inputs and controls

The sole spec is `specs/prime20-b64-batchtable.json`, using registered
`icv1-fp20-t727-cd198a38`, subgroup order 562427, `prime_search`
bits 20 and seed 59297. Eight new public targets use target seed 20261015;
no previous cohort uses them. Each arm has the same base size 64, checked
full-rank linear algebra, algorithm seed, `max_trials=5000000`, and
120-second timeout. Run one warmup and three measured rounds per arm,
alternating order, one thread. `rho.negation` is the matched reference
and repeats as A/A. Planned executions: 128.

Baseline `mitm:m=2,negation_folded=1` uses the existing table builder.
Candidate `mitm-prime-batch:m=2` uses the batched, x-only builder with
the same generic pair-table probe and relation loop. Preserve all failed
executions. Stop after one spec run and full replay. The comparison is
accepted as a correctness-preserving engineering change only if all
measured runs verify and paired target outputs, relation trials, table
entries, and total GAE agree exactly. Report any mismatch as a failure.

Prime factor-base lifts, Legendre tests, square roots, target probes,
hash lookups, row operations and some rho native work are unpriced in
ecbench.gae. The field-operation tally is separate from GAE. Wall time
requires L2 and an A/A noise check; L0 values remain descriptive. Cold
single-target setup is charged in each arm. No extrapolation to larger
fields or full-cost ECDLP follows.
