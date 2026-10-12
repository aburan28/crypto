# Frozen 2,048-target power extension, 2026-09-25

The 512-target masked pilot in PROTOCOL_MASK.md is preserved as a
separate, negative confirmation result: the held base had no group hit.
An exact, target-independent census of its complete three-sum set
found 7,658 distinct group points per 42-point base out of 2,099,948
source group points. Under an approximately uniform full-group target
draw, the expected count is only 1.87 in 512 attempts. Zero held hits
therefore does not diagnose the prescreen. This extension's sample
size and analysis are fixed before measuring any attempts 512..2047.

Use the **same** two geometry-selected rank-5 raw bases, degree-7 map,
secret, SHA-256 coefficient/mask domains, uniform-affine mask
rejection sampler, target-orbit hash split, complete pair lookup,
generalized affine S3 screen, projected-log row, and rank/verification
rules of PROTOCOL_MASK.md. Extend the frozen stream from 512 to exactly
2,048 attempts by continuing the attempt index. Keep the first 512
records and all outcomes as an unchanged prefix; verify their target
coordinates, masks, first witnesses, certified skip counts and rank
gains against the saved 512 run. No base, target or failure is chosen
based on observed relation yield.

The exact 7,658/2,099,948 coverage gives expected 7.47 hits per base
at 2,048 approximately uniform full-group targets. The binomial
probability of fewer than two hits is approximately 0.49%; this is a
power calculation, not a guaranteed full-rank claim. Record every
target, including duplicates/infinity and no-hit attempts. If either
base still lacks verified full rank, report a null cold-to-rank cost
and do not extend this run post hoc.

Independently replay every group hit, cofactor-projected row, skipped
S3 residual, rank gain and scalar recovery. Report double-held-out
hit retention, certified skips, full-stream and cold-to-first-rank
cost, preserving source/transported parity. The decision gate remains
the one in PROTOCOL_MASK.md: a lookup reduction is insufficient if
screen setup and evaluation raise charged solve cost. This toy
extension cannot imply an ECC2K-130 or matched-rho crossover.
