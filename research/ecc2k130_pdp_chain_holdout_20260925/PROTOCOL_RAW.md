# Frozen raw-orbit follow-up, 2026-09-25

The first preregistered pilot in PROTOCOL.md used signed-Frobenius orbits
**after** projecting into the tiny order-421 subgroup. Its complete pair
table saturated the subgroup: every target found a three-summand witness
at the first third-point candidate, so no affine contradiction could be
observed. Preserve that run as an exploratory negative control. This
follow-up was specified before measuring any raw-base target outcomes.

Keep the field, source curve, order-421 challenge subgroup, degree-7 map,
target-orbit split, S3 residual rule, complete pair-table baseline,
group-witness/rank verifier, and separate-unit cost accounting from
PROTOCOL.md. Change only the factor-base construction, relation label,
and sample size as follows.

Generate exactly 20,000 nonzero x candidates by Python 3.12
`random.Random(2026092503).randrange(1, 1 << 21)`, retaining their trial
indices and all rejection counts. For each candidate calculate the binary
rank of its 21 repeated-square x orbit. Evaluate a curve lift only when
that rank is at most 9; take the lexicographically first lift. Require
the full signed-Frobenius point orbit to contain exactly 42 distinct
points and the cofactor-4988 projection of the representative to be a
nonzero order-421 point. Deduplicate signed point orbits, then choose
the two eligible orbits with lowest x-span rank, using SHA-256 of
`pdp-chain-raw-orbit-v1|x|y` as tie break. Assign train and held base
by hash order. Their 42 raw points must be disjoint. Transport all 42
points of each through the fixed oriented map, preserving index order;
the codomain sets are transported source orbits, not native Frobenius
orbits. Charge every candidate rank, lift, projection, and map operation.

The raw base points need not lie in the challenge subgroup. Let
`R=[4988]P` be the representative's subgroup projection. If the
signed-Frobenius weight of a raw base point is `w`, its projected
subgroup log is `w log_G(R)`. A **group-verified** decomposition
`P1+P2+P3 = [u]G + [v]Q` therefore yields
`(w1+w2+w3) log_G(R) - 4988 v log_G(Q) = 4988 u (mod 421)`.
Use this equation for rank and independently check both recovered logs
by scalar multiplication. The same weights and coefficient equation
hold for transported points under the subgroup isomorphism.

Generate 512 natural synthetic subgroup target attempts using the
unmodified SHA-256 domains and fixed secret of PROTOCOL.md, extending
the attempt index to 511. Group the target attempts by the same
source target-x Frobenius orbits; use the same geometry-fixed held
orbit set. Train/confirmation are the two disjoint base/target cells;
both mixed cells and infinity remain in the ledger.

For every attempted target, the baseline and screened policies scan
all 42 candidate third points until the first exact pair-table witness
or exhaustion. The screen computes its exact x-span per raw base,
precomputes the base products once, and lazily memoizes each residual
abscissa. It skips a pair lookup only on a certified affine contradiction
for `S3(x1,x2,x(T-P3))` with the appropriate curve `b`.
Independently check every skipped residual against all finite base-x
pairs by direct S3 evaluation; also check that no skipped group pair
exists. A lost algebraic root or group hit fails the run. Record
per-target scans, misses, skips, first witness, independent row gain,
and verified full-rank stop, or an explicit no-rank failure.

Primary decision: on the **double-held-out** confirmation cell, does
screen setup plus all profile evaluations cost less than the complete
pair-table baseline after charging saved lookups and preserving all
group hits and rank gains? Report field multiplications, squarings,
curve additions, row XORs, lookup counts, CPU, and full cold-to-rank
cost separately; no synthetic scalar conversion between units. A
zero-skip or higher-cost result rejects this prescreen placement for
this raw-orbit toy policy. A positive result remains a diagnostic,
requiring the n=37/41/53 full-rank and matched-rho gates before any
ECC2K-130 claim.
