# Native endomorphism search validation

This integration replaces the delivered Python execution path with Rust. The
original ZIP is immutable historical evidence, with SHA-256
`89d4f459ef440dd238864d4a39a52c9b59d2c6e086d6f2566a81aad7f82c83d6`.
Its timings describe variable-time Python on small synthetic prime-field
models. They are not native measurements and establish no ECDLP improvement.

Before native execution, freeze these inputs: valid order discriminants through
absolute value 1000, degree bound 1000; the four model triples in the archived
demo; and the archived Sage reference fixture. The source archive pins the
original commands, seeds, raw timings, schema, tests, and run summary.

Hypothesis: integer-exact native enumeration agrees with an independent
rectangular enumeration and with all 500 archived order summaries. Explicit
map constructions reproduce the archived geometric degrees and subgroup
actions; lattice decomposition and complete scalar multiplication agree with
independent binary arithmetic for every scalar in each selected subgroup.

Acceptance: zero mismatches, no failed map or group checks, all CLI rejection
controls pass, incomplete scans cannot certify absence, and altered evidence
must be rejected. Stop on the first correctness failure. Search limits are
explicit; unsupported fields, kernels and paths remain unsearched. Verification
is bounded separately by the laboratory field cap.

Native `--count-ops` is a scalar-arithmetic stage diagnostic. It records field
multiplications, field additions, inversions, group calls and integer
decomposition work separately, including map evaluation, per-call tables and
one lattice setup per batch. No calibration combines these units. There is no
wall-clock claim, no ratio to rho, no end-to-end cost, and no reported speedup.
New timing comparisons require the repository's isolation and A/A protocol.

This is not an IC experiment or cross-method ecbench session. No relation
model, collector, rho walker or recovered public logarithm is produced. The IC
leaderboard and progress chart receive no new performance figure. Curve
identity metadata must be registered before catalog or UI publication; the
runner reports registration as unresolved and does not mutate the registry.
