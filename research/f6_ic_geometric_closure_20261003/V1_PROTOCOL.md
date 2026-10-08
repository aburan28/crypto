# F6-IC v1: preregistration before implementation or execution

The v0 pilot is a negative receipt: 11 support checks, zero residual lookups,
and exactly the same 16 reductions and 9 splits as inherited F4 on its frozen
one-target workload. This version changes the exact closure trigger. Its
source and candidate IDs must differ from v0; no v0 file is rewritten.

For `m=3`, when exactly one summand x-code is fixed and the enumerated base
has at most 256 points, enumerate its exact point lifts. For each lift,
enumerate every exact usable base point at a second summand, compute the
target residual for the third, and look it up in the exact base. Apply all
currently fixed bits of both remaining x-codes, ignoring affine-definition
placeholders. Return a relation only after full group replay. If no choice
works, refute this branch as containing no usable IC relation. This is a
complete finite enumeration, so a negative result is exact. At bases above
256 points or any other summand count, retain v0's two-fixed residual rule
and inherited F4 continuation. Every group addition, lookup and replay is
charged to target PDP or relation check, never omitted.

First verify the new one-fixed trigger on n9 complete enumerations, including
positive and negative cases, permuted variables, missing support and
affine-defined placeholders. The existing F4 versus F6 exact-enumeration
controls must still pass. A bounded prepared-worker pilot then uses the same
certified n17 mathematical state, 62 actual usable points, 29 folded columns,
`m=3`, degree 3, 8192 nodes, one Rayon worker, one previously unseen public
point from fixture seed `20261003020`, algorithm seed `20261003035`, one
target-query batch trial, and at most 32 queries. Fixture construction,
source and binary hashes, candidate/workload IDs, and exact inputs are
committed before timing. Run six fresh processes under a 600-second cap in
the order F4/F6, F6/F4, F4/F6. Do not change the seed or replace a run based
on its result. Retain timeouts, failures, zero-yield attempts and every
five-phase ledger.

The decision gate is all six processes verified with identical scalar and
exact phase closure, plus **fewer F4 reductions or splits** in at least one
F6 run without exceeding the point-work budget. A favorable local wall time
is exploratory only under the CPU isolation rule. There is no IC/rho claim
from this pilot. A new cross-method claim requires `ecbench`, the same public
target, strong rho, operation accounting, a sealed replay and an independent
isolated-host receipt. If the new branch enumeration increases complete
online work, preserve it as a negative result and stop this variant.
