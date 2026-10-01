# Boundary online schema repair discovered during preparation CI

The [first final-head integration run](https://github.com/aburan28/crypto/actions/runs/36769633539)
on `21222aad24ed1c3dc39828447ecdd08ec13a8347` passed all 398 main IC harness
tests, then failed two of the 16 boundary-autolab controls. The local boundary
suite reproduced the same two failures; the current main branch's ledger had
the same schema mismatch.

The accepted single-target validator/tests arrived between PR #1064's checked
head and its merge. Its live ledger schema still required the old `ic_cost`,
`rho_cost`, automorphism and whole-series fields, enumerated only legacy timing
classes, and omitted the five online phase keys. Consequently a complete
single-target report failed both timing-class validation and phase closure
(the empty declared phase list summed to zero). A legacy report also lacked
`target_count`, but the ledger did not report it among missing required fields.

The repair changes only `measurement_schema.vs_rho` in
`docs/ic/boundary_targets.json`. It declares the canonical candidate/workload/run
bindings, exactly one shared point, online times and ratio, five exclusive
phases, interval events/stages, independent replay certificates, matching
resource envelopes, verified scalars and rho policy required by the accepted
validator. The timing class is `single_target_online_wall`. Every historical
row and every other ledger field is unchanged; parsed equality outside this
schema section was checked before writing.

The validator additionally requires an integer `target_count`; Python equality
must not let `true` or `1.0` masquerade as an integer one-target record. Existing
adversarial controls now include both cases. The corrected boundary suite passes
all 16 controls locally, including complete online acceptance and rejection of
batching, mismatched targets/resources, false scalar/replay records, wrong
ratios/phase sums and legacy whole-process timing.

This is contract/validation repair, with no native benchmark, ledger promotion
or new performance result. Boundary `PASS` is schema validation; full independent
scientific/source replay remains required. The preparation certificates and
their mathematical/source seals are unchanged. Integration must pass on the
repair commit before the PR is accepted.
