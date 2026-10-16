# Nine-arm held-out fiber comparison (V4)

**Hypothesis.** Direct inverse-Vélu and full source-image construction
produce the same signed eight-column destination factor base and relation
support as the four polynomial/direct-pair arms in the frozen V3 panel.
Inverse construction may reduce setup relative to the full census, but its
cold complete cost must still be compared with same-target BSGS and rho.

**Frozen inputs.** Reuse the public four order-229 target positions 64–67,
checker-only scalars, signed eight-column projected factor base, relation
random seeds, 1,832-trial cap, three rounds, and rho controls from
`HELDOUT_SEVEN_ARM_SPEC_V3.json`. The source and degree-three quotient
have different `j` invariants. The inverse and full-image executables
must construct their base from curve data, then assert exact equality
with the published eight columns. Public preimage witnesses are for
verification only. Freeze both solver binaries and the checker before
execution. Do not overwrite V3 records.

**Run.** Interleave nine arms in round-dependent cyclic order:
χ resultant, direct pair, source S3, quotient S3, direct inverse-Vélu,
full-image census, BSGS, signed rho, rho A/A. Launch every target as a
cold process with the same native whole-child meter, 300-second hard
wall cap, and 512 MiB observed Darwin RSS ceiling. Preserve exit,
stderr, process CPU/wall/RSS, and all raw scalar/row receipts.

**Acceptance.** Stop and preserve the failed row on a factor-base
or relation-set mismatch. A separate Rust checker must replay 108
scalars on both curves, all relation random streams and witnesses,
row ranks, BSGS and rho outputs. The native summary asserts exact
trial and relation-count support for all six relation arms in each
paired round-target cell. Report one cold counted-group `S_group`
table and descriptive process CPU separately. These do not establish
a portable primitive `S` or an isolated wall crossover.

The next gate after this matched comparison is nonoverlapping primitive
instrumentation with calibrated conversion intervals, plus independent
Linux replay and physically isolated timing.
