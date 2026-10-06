# Amendment 2: complete the exact k=16 core below a memory gate

Registered after the Amendment 1 cap was frozen in commit `7f7a2d398`
and before changing the reducer or running this follow-on. The exact
2,260-row core contains 6,376,371 term occurrences and hit the
1,500,000-column cap at 1,500,001; no rank or affine conclusion was
obtained. Its total number of distinct columns cannot exceed its term
occurrences, so a 6,500,000-column cap admits this *same* core if memory
allows. This is a resource-bound continuation, not a changed algebraic
candidate or a speed benchmark.

Keep the registered K0 curve `icv1-f2m83-tm6151469093347-debefd74`,
standard dimension-18 base, public T001 at torsion offset zero, 332
original equations, first 16 source multipliers, 32 sampled degree-four
candidates per row, and exact two-pass certificate. Recompute and require
the frozen 1,928 unresolved rows and BLAKE3 witness digest
`ba12cc6bb138b135f514f6bb6198db21286101b4a55eeb6f42a5af12e4a49d74`
before constructing the core. Retain the planted `[0,2,4,6,8]` control
and full group replay.

Add an explicit column-cap parameter to the existing deterministic root
reducer, preserving its 1,500,000-column default for all existing callers.
For this one amended probe, set the cap to 6,500,000. Run one native Rust
release process for the ordinary T001 offset-zero core, after one planted
control process, with one thread, 7-GiB **observed** RSS gate and
300-second process timeout. A separate process wrapper samples RSS and
kills the process if it exceeds 7 GiB; record the samples and exit status.
Keep the exact same reduction ordering and no back-substitution.

Structural success means a completed reduction yields a contradiction
or a nonconstant affine row on source bits only. If that happens, repeat
all four torsion offsets with independent full-group witness verification
before proposing a solver. A completed reduction without either rejects
this k=16 prolongation as a source-pruning route. Timeout, OOM, RSS kill,
column cap, or nonzero exit is inconclusive. Do not turn the previous
column-cap row into a negative result. Preserve every raw stdout/stderr,
build log, status, observed RSS, source/binary hashes, and SHA-256 receipt.
Timings remain exploratory on the contended host; no complete F6 or
one-target IC-vs-rho speedup can follow from this structural probe.
