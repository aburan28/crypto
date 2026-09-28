# Frozen cofactor-masked raw-orbit PDP pilot, 2026-09-25

The prior predeclared projected-base and raw-base runs are retained. The
projected-base run saturated the tiny order-421 target subgroup; the raw
rank-8 orbit run had no group hits because subgroup-only targets did not
meet its observed triple-sum set. The following distribution and
morphology selection were fixed **before measuring masked target
relations**. This is a third, distinct toy diagnostic, not a repair to
either negative receipt.

Retain the F_(2^21) source/codomain pair, exact degree-7 orientation,
full-point map, order-421 challenge subgroup, m=3 residual S3 rule,
complete unordered pair table, S3 constant `b` on each curve,
cofactor-projected row equation, and separate native-unit accounting
specified in PROTOCOL.md and PROTOCOL_RAW.md.

## Bases

For trial indices 0..199999, generate nonzero x with Python 3.12
`random.Random(2026092504).randrange(1, 1 << 21)`. Compute the rank
of its 21 repeated-square abscissae. Only if rank <= 7, evaluate the
curve lifts and select the lexicographically first. Require a full
42-point signed-Frobenius orbit and nonzero order-421 projection
`[4988]P`. Deduplicate complete signed orbits. Select the two
eligible orbits with smallest x-span rank, breaking ties by SHA-256
`pdp-chain-mask-orbit-v1|x|y`; assign train/held base by that hash
order. Charge all 200,000 rank probes, selective lifts, projections,
and selected orbit/map construction. Neither relation hits nor screen
outcomes may affect this selection. The two 42-point raw bases must
be disjoint; their 42-point transported images have identical index
order and relation trajectories. On the codomain, coordinate
squaring is not a same-curve Frobenius endomorphism.

## Targets and split

Use the same seed 2026092503, secret, `u_i,v_i` hash domains and
512 attempt indices as PROTOCOL_RAW.md. For each attempt independently
sample a full-curve affine point `W_i` by rejection: SHA-256
`pdp-chain-mask-v1|2026092503|i|retry|x` gives an x in
`[0,2^21)`; SHA-256 with final `|sign` gives one sign bit.
Reject x without a curve point, and reject sign 1 at x=0. Accept
the corresponding full point otherwise. This maps each accepted
coordinate/sign pair one-to-one to an affine point, so `W_i` is
uniform over affine points, excluding only infinity. Record retries.
Set the published torsion mask `M_i=[421]W_i` and target
`T_i=M_i+[u_i]G+[v_i]Q`. Require `[4988]M_i=O`.
The mask is chosen without knowledge of any relation and makes
`T_i` approximately uniform on the full curve; the exact deviation
from a uniform full-group draw comes solely from excluding infinity
from `W_i`. Keep mask, infinity, duplicate and failed-target records.

Map each target and mask to the codomain and compare to direct
codomain construction with the same coefficients; a mismatch fails.
The relation after multiplying by cofactor remains
`(w1+w2+w3)d - 4988 v k = 4988 u (mod 421)`, because
`[4988]M_i=O`. Verify every accepted triple both in the full group
and after cofactor projection before admitting its row to rank.
Recover both `d` and `k` only from verified full rank, then
check `[d]G=[4988]P` and `[k]G=Q`.

Partition target attempts by their canonical source target-x
Frobenius orbit; hold an orbit if the first 64 SHA-256 bits of
`pdp-chain-mask-target-orbit-v1|key` are congruent to zero modulo 3.
This fixed hash rule is independent of hits and groups repeated
targets/signs. Train and confirmation are the train/held base
crossed with unheld/held target orbits; retain both mixed cells and
infinity separately. All 512 attempts are retained in every arm.

## Decision

The complete pair lookup is the comparator. The screen may skip a
lookup only when its cached affine row span certifies a contradiction
for the fixed-third residual. Independently verify every skipped
residual against all finite base-x S3 pairs and the group pair table.
Baseline and screened first witnesses, group-hit trajectories and
independent rank gains must match exactly. Original/transported group
trajectories must match as an isogeny-invariance check.

Report setup, 512-attempt full-stream cost, and cold cost to the first
**verified** full rank separately, including masked-target generation,
complete pair-table construction, all screened feature work, failed
attempts, map construction, and rank/verification. Field operations,
curve additions, row XORs, lookups, CPU and wall time remain separate
units; do not sum overlapping counters. If no full rank, preserve the
failure and leave cold-to-rank null.

Call the screen promising for this toy policy only if the double
holdout has certified skips, no lost algebraic/group roots, identical
rank, and lower charged solve cost after setup. A cheaper lookup
count alone is insufficient. Even a positive toy result cannot
override the n=53 matched-rho no-crossover or n=131 four-sum coverage
bound without a new fully charged full-DLP gate.
