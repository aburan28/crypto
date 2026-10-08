# Preregistered sparse-support four-policy n37 PDP gate

Status: **protocol only**. Commit and open this PR before generating the new
public targets or running the producer. This is a factor-base/PDP diagnostic,
not a cold IC/rho comparison or a claim about ECC2K-130 attack speed.

## Question and fixed input

The earlier equal-useful-size K42 comparison reached m≤3 on every one of its
2,048 held-out points, so it could not distinguish the two seed-selection
rules. Fix **K=16** by taking the first sixteen complete signed-Frobenius
columns of each of the four immutable policies in
`n37_four_policy_support_20261003/RESULT.json.gz`. Each resulting base has
16 log columns, 592 signed classes and 1,184 physical points. The archived
support gzip SHA-256 is
`8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5`;
the expanded SHA-256 is
`05664103a6dab090a8b4298f34c0964e6f253f625cd73d661562d85a2fce87af`.
No column or point may be changed after the targets are frozen. Source and
transported, and descendant-native and pullback, must have identical exact
hit decisions, witness labels, relation rows and rank trajectories after the
degree-73 map. Only original versus pullback compares different selection.

The uniform-target counting ceiling for m≤2 is
`(1+1184+1184·1185/2)/230603167 = 702705/230603167 ≈ 0.305%`.
For m≤3 the formal multiset capacity is
`C(1186,3)=277334240 ≈ 1.203r`. That capacity permits either hits or misses;
it is not an expected yield, so the fresh point experiment is necessary.
The attack reference remains same-Q signed-Frobenius rho in a later cold
process. This gate leaves complete `S`, online speedup and n131 transfer null.

## New public target freeze

After the protocol and native generator are committed, generate two ordered
blocks of 1,024 source-group public points. For block `b∈{0,1}` and candidate
index `j≥0`, SHA-256 the UTF-8 string
`n37-four-policy-sparse16-20261004|b|j`, interpret the digest as a big-endian
integer, and reduce modulo `r=230603167`. Reject zero. Compute `[d]G` with
the registered source generator. Reject a signed-Frobenius orbit already in
the new blocks or in the earlier b03/b04 point-only blocks used by the K42
gate. For this curve the orbit key is the smallest integer among the 37
Frobenius conjugates of x; sign preserves x. Accept the first 1,024 points
per block, recording candidate indices and rejection counts. Commit the two
point-only files, verifier-only scalar labels, source/input hashes and an
independent generator replay **before** running any PDP arm. Neither PDP
producer nor its replay may read the labels.

## Exact search, rank and decisions

For each policy, enumerate the identity, all 1,184 singletons and every
nondecreasing pair into a complete full-point table. There are 702,705 raw
entries. Keep the first enumeration witness for each sum. For every public Q,
scan identity then all 1,184 factors in frozen order, seeking `Q−P_i` in the
pair table; recompute each accepted group sum. A proved miss requires the
complete table and all 1,185 residuals. Also record exact m≤2 as a secondary
control. Failures, caps, malformed points and incomplete scans are unknown,
never misses.

Use SplitMix64 seed `0x6e33_375f_7331_3664` to generate the first 256
distinct nonzero scalars modulo r. On every policy, query the corresponding
`[a]G` or mapped leaf point in that fixed order, verify each witness, and
insert its coefficient row into an exact rank-16 solver. Stop at rank 16 or
the cap. If full rank, solve base logs and verify each by scalar replay; for
every public-Q hit, recover and verify its logarithm on both curves. A miss
or incomplete rank leaves that target's recovered logarithm null. Preserve
every rank miss, dependency, row, and attempted target.

The independent native replay must rebuild source support tables by a
different lookup structure and use the general binary-curve group law for
every accepted witness, rank row and scalar. Check the leaf results through
the explicit isogeny map and paired-policy identities. Mutation of a target,
witness, rank row, hit flag or decision must fail replay. Record source,
support, target, result, replay and binary hashes, host facts, operation
counts, elapsed phase times, memory and any failures. Local wall times are
descriptive L0 diagnostics. Update the canonical scoreboard with the stage
result and leave whole-method speed fields unset.

The fixed-block selection-lead rule is an absolute m≤3 hit gap of at least
21/2,048 **and** a two-sided exact conditional McNemar p-value below 0.01,
with complete paired replay. If both policies hit all 2,048, classify the
experiment as saturated again. If neither rule holds, classify it as no
fixed-block lead. Incomplete rank is reported separately, without imputing
logs. A lead licenses a new preregistered cold single-target and matched-rho
panel; it is not itself an attack improvement. A negative outcome redirects
the held-out yield comparison to a larger subgroup rather than retuning
these points.
