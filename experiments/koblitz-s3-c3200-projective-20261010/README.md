# C3200 projective target checks

This candidate keeps the exact C3200 factor base, S3 root index, relation
collector, rank solve, and target decomposition. It uses the repository's
portable López–Dahab projective scalar engine for the target subgroup check
and the final recovered-scalar replay. Each multiplication performs one
field inversion at its end, replacing the affine double-and-add path that
inverts at every point operation. Both checks remain inside the one-target
online interval and retain the same subgroup and group-equality conditions.

The source includes a cross-engine exact-output test for N41 and N53 on
zero, small, middle, and subgroup-boundary scalars. After it passes, the
frozen T025 control will require the same scalar, first relation witness,
factor-base digest, full rank, root-table layout, and every non-timing report
field as the verified C3200 reference. It will record the five charged target
phases, peak RSS, compressed raw report, and a separate Sage replay. The
validation and recovery phase costs will be compared first; a fresh-target
paired online screen is needed for the owner-set 2× threshold.

After writing the primary T025 result, the control runs 256 balanced
affine/projective check pairs against that same public point and scalar.
Each row records subgroup-check and replay wall time, main-thread CPU time,
and exact outcomes. The repetitions share setup and a known point, so they
measure only this check stage. Their raw ledger is archived separately from
the one-target online result.
