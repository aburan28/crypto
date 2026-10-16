# Larger logical batches on the selected split/block-inverse G7 kernel

Preregistered on the local AWS g7.2xlarge / RTX PRO 4500, GPU UUID
GPU-827e82b5-4739-d339-d95e-b7214337553e, at 165 W. The goal is 12 billion
complete ECC2K-130 scalar rho updates/s on one GPU. The immutable selected
reference is `build/ecc2k130-local-packed`, SHA-256
`8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02`.
Generic work is `sqrt(n/262)`, ratio 1; full-DLP S is null.

## Hypothesis

The selected logical batch has 16 walks, split across two CUDA threads that
process eight denominators each. Keep that two-thread split, 256-thread block,
whole-block inverse, polynomial arithmetic, eight-way jump function and all
state semantics, but test logical batches 24 and 32. Their physical threads
process 12 and 16 denominators. A larger local prefix amortizes one block-root
inversion and finite-prefix endpoints over more complete scalar updates.

Use the selected four-slot shared-X cache in every row. This retains the same
per-block shared allocation and three-block residency, but caches a smaller
fraction of X coordinates in the larger batches. Extra global coordinate
traffic, longer dependency chains or compiler scheduling can exceed the saved
inversion work. Historical batch-32 measurements predate the current physical
split, block inverse, polynomial state, weighted prefix, shared-X cache,
seven-stage conversion and CUDA 13.4 assembler, so they do not measure this
combination.

Compile from an isolated source that reproduces the validated B16 control.
Only broaden the compile-time batch guards to admit even batches from 16
through 32; loop bounds, state sizing and record headers already derive from
`ECC_BATCH`. Reject unsupported sizes and configurations.

## Gates and timing

Before timing require exact built-in arithmetic checks; field, block-inverse,
storage and shared-X units at each new batch; a normalization oracle showing
that B16/B24/B32 produce the same per-walk coordinate, seed, iteration and dead
state from the same global walk IDs; matching sorted DP records; bidirectional
resume within each batch; ragged populations; zero drops; and walk memcheck,
initcheck and synccheck. A batch-size checkpoint header intentionally prevents
cross-batch resume and must be rejected.

Use 8,388,576 total scalar walks in every timing arm: 524,286 logical workers
for B16, 349,524 for B24 and 262,143 for B32. With 1,024 steps and four
launches, each sample performs 34,359,607,296 complete updates. Run three
interleaved paired repetitions, reversing the middle order. Advance only if
the paired Student-t 95% log-ratio interval versus selected is wholly above
one. Promotion requires five fresh benchmark and DP34 pairs, both intervals
above one, exact normalized records, zero drops and the full validation suite.
Preserve regressions. No hardware, power, driver, service, cloud-resource or
publication changes.

## Result

Both candidates regress and do not advance. B24 measures 6.121832 B/s versus
6.274186 B/s selected, paired ratio 0.971864 with 95% CI
[0.961912, 0.981919]. B32 measures 5.994967 B/s, ratio 0.948253
[0.930625, 0.966215]. Every sample contains 34,359,607,296 complete updates.

The correctness gate passes in full: built-in arithmetic, native control
identity, storage, block inverse, actual-walk shared X, normalized complete
states and DPs, CPU replay, ragged and aligned populations, same-format
resume, cross-format rejection, memcheck, initcheck and synccheck. Telemetry
shows the four-slot cache no longer covers enough of the larger physical
batches: memory activity rises to 36%/50% and sustained clocks fall under the
fixed power cap. The selected B16 runtime remains unchanged.
