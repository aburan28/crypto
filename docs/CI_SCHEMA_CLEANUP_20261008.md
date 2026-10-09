# CI and schema maintenance — 2026-10-08

The cleanup starts from `c70c32d486a3ac7531fe27f7193d9f09caa58344`.
The existing developer checkout is preserved; changes were prepared on
`codex/ci-schema-cleanup-20261008` in an isolated worktree.

## Source integrity

The `fe07f6088`, `ec13b3b87`, and `92600ae6d` union merges combined conflicting
function bodies and declarations. The resulting main branch failed compilation,
including the library, IC command, curve cover command, and several examples.
The cleanup reconstructs complete implementations from the recorded merge
parents and preserves the later intended changes. The rank example retains
the complete rank and shared-log implementations in separate modules behind
their existing options.

The curve registry had a missing record boundary; its JSON delimiters are
repaired. The duplicate `sect113r1` records were consolidated with both
representations retained, leaving 321 distinct curves. The ECC2K130 Makefile had
duplicated continuation fragments. The shared cover checker now receives its
model-normalization module in both commands that include it.

The repairs also restore wide-field kernel coordinates, preserve architecture
capability checks and portable arithmetic, and route models without a
Frobenius endomorphism through the appropriate quotient paths. Endomorphism
information remains unknown for those representatives. The staged prime-field
solvers require a determined target column and independent scalar replay;
unused factor-base columns may remain free. Exhausted relation queries are
included in the recorded query counts.

## Validation boundaries

Strict Clippy checks apply to the maintained library, commands, tests, benches,
and examples rooted in `examples/`. The two hash-bound S3 comparison snapshots
remain byte-identical and remain covered by the all-target compilation gate.
`scripts/clippy-maintained.sh` obtains the maintained target list from Cargo
metadata and fails if metadata or selection is unavailable.

Correctness checks use the existing fixtures and verification requirements.
The Gröbner controls distinguish a reported timeout from a completed empty
solution set and retain their original completed-case requirements. Test
concurrency is bounded during local verification to prevent simultaneous
controls from starving each other's wall-clock budgets.

## Evidence and delivery

Validation commands and raw success/failure logs are retained under the
cleanup workspace's `logs/` directory. Rust validation uses toolchain 1.98
with matching Cargo, rustc, rustfmt, and Clippy binaries. Local runs are native
ARM64 correctness and compilation checks. GitHub Actions provides the
repository's Linux checks and release build matrix.

The pull request's current check results and final commit identify the
delivered revision. Performance promotion continues to require the repository's
single-target accounting, independent verification, and isolation receipts.
