# Amendment 1: reduce the certified `k=16` core

Registered after freezing the first result in commit
`071d5ef53` and **before** changing the probe or running this follow-on.
The first protocol conditioned core reduction on the full `k=90`
certificate. That full set left 15,822 unresolved rows and correctly
stopped. The exact `k=16` certificate unexpectedly left only 1,928
unresolved rows, below the same 2,000-row core threshold. This
amendment asks whether the smaller certified core can yield a
source-only affine consequence within the existing resource cap.

Keep the same registered K0 curve, standard dimension-18 base, public
T001 at torsion offset zero, original 332-equation system, first 16
frozen source variables, 32 candidates per row, two-pass exact
certificate, and planted `[0,2,4,6,8]` control. The new binary must
recompute the `k=16` certificate and require exactly 1,928 unresolved
row IDs and the frozen witness digest
`ba12cc6bb138b135f514f6bb6198db21286101b4a55eeb6f42a5af12e4a49d74`
before core construction. This guard prevents an unnoticed algorithm
or input change from being read as a new core result.

Construct the core from all original equations and exactly those
1,928 unresolved products. Because every removed prolonged row has a
private degree-four monomial against the full `k=16` row set, removing
it preserves **all possible degree-three-or-lower consequences**.
Run the existing deterministic root reducer with its unchanged
1,500,000-column cap, observed 7-GiB RSS gate, one thread, and a
300-second process limit. Preserve a column cap, timeout, OOM or
nonzero exit as an inconclusive result.

Structural success requires a completed reduction yielding a
contradiction or a nonconstant affine row involving only source bits.
If successful, repeat all four torsion offsets and verify any resulting
point decomposition in the full group before proposing a solver.
If the core completes without either, reject this `k=16` selective
row set as a pruning path. The `k=90` unresolved core remains a
separate unsolved question. Timings on this contended Mac remain
feasibility diagnostics, not CPU speedup claims. No complete F6, IC
online or same-target rho comparison is implied.

Build one new native Rust release probe, freeze its source and binary
hashes, run planted and ordinary `k=16` once each, and retain raw
stdout/stderr, exit codes, memory, phase timings and SHA-256 manifest
without overwriting the first panel.
