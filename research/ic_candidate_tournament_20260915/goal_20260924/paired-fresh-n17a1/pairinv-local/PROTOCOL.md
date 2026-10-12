# Local pair-table incumbent and rho on the frozen n17a1 point

This is a diagnostic macOS ARM64 port of the accepted instrumented pair-table
source, not a reuse of its qualified Linux binary or an instruction-count
qualification. The source manifest is
`240b8daa0478aacb1f6fc248de9fd5ed9b0773c59d3ae869ced0bab8be438378`.
The checked prepared source, complete registry dependency contents, compiler,
Cargo flags, build log and resulting binary are bound in the local build
receipt. The archived Linux target configuration is overridden explicitly for
`aarch64-apple-darwin`; no source file is edited for the local build.

Before either measured arm runs, an inventory-only worker invocation resolves
the actual factor base and its sign/Frobenius columns. Registration checks the
inventory with the independent curve arithmetic, fixes the complete candidate
and workload records, and seals the job and runner bytes. The factor base is
the historical sampled subgroup-orbit policy with seed 43 and point bound 102.
Its actual usable point count and folded column count come from the inventory,
not the point bound. This differs from the standard-subspace factor base used
by the F5 and SAT arms; any comparison is of complete pipelines and base
policies, not an isolated solver substitution.

The IC arm uses the optimized three-summand pair table, one collection trial
per batch, at most 65,536 trials, and the final tiny Gaussian relation solve.
The rho arm uses the same public point, same checked source and binary, signed
Frobenius walks, the historically selected four requested walks and at most
65,536 iterations per restart. The source's small-group affordability rule
reduces this n17 arm to one effective interleaved walk; registration and audit
retain both requested and effective widths.
Both receive the target point `[52411,72106]` directly, with seed
`2026092948` recorded only as provenance. No known scalar is supplied. They
run on the same physical local host, with one worker thread and no hard outer
wall or memory limit. Child high-water RSS and whole-process wall are
diagnostics; the headline candidate interval is the producer's exclusive
target-dependent interval after reusable preparation through scalar replay.

Each arm executes once and retains its raw stdout, stderr, exit code,
process clock, inventory, worker and audit. A failure or bounded incomplete
result remains a row and has no verified online cost. Independent replay
checks the curve/subgroup, every relation scalar, matrix/rank path, recovered
logarithm and exclusive phase closure for a complete IC result. The rho audit
checks the exact walk width, scalar replay and exclusive online interval.

The [target allocation](../target-panel.json) predeclared F5, SAT, incumbent
and rho in that order. SAT was actually launched before F5. Preserve and
report this scheduling deviation when interpreting the four rows; it cannot
be repaired by relabelling their execution order. The four-walk setting
matches the historical selected source/configuration, while the local binary
and hardware are separately bound and have no Linux calibration. These
one-shot diagnostics do not establish a global winner, calibrated operation
speedup, or confirmation-stage claim. All four rows and an explicit host/order
review are required before even a conditional same-point online comparison.
