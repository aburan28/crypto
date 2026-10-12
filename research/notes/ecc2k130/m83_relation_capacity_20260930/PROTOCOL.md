# Exact relation-support capacity screen, 2026-09-30

This is a deterministic counting audit of the already published m=83
full-orbit control and the ECC2K-130 subgroup order. It is not a PDP run, a
new target sample, a solver comparison, or a discrete-log timing experiment.
The hypothesis is that the m=83 control bases are too small for their zero
natural m=3 relations to distinguish solver quality or native factor-base
quality. The reference is the exact cardinality of the prime-order subgroup.

## Frozen inputs and question

- m=83 source: `../m83_solver_matrix_20260927/results/validated_summary.json`,
  SHA-256 `ff5b525fd188fb30f13b1959eec823ae3a239424c8c493503382ccc751c00261`.
  Read both `signed_orbit_size` values, 332 and 498, and check their recorded
  unordered pair counts against `B(B+1)/2`.
- m=83 subgroup order and four natural targets per seed: the adjacent
  `RESULTS.md`, SHA-256
  `00a1e65642713dc04a71e158386e6d0048806fe87d73233766ea041afcc0e687`.
  Use `r83=2417851639230796216685689`.
- ECC2K-130 challenge subgroup order: `research/ecc2k130_relations/relations.py`,
  SHA-256 `0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f`.
  Use `r131=680564733841876926932320129493409985129`.
- Treat B as the number of **physical distinct subgroup points** in a fixed
  factor base, after any sign/Frobenius expansion. Count unordered m-tuples
  with repetition. For a uniform subgroup target, distinct supported sums
  are at most `min(r, C(B+m-1,m))`; collisions only lower support.

Compute the one-target upper bound for m=3,4,10 at each archived m=83 B, and
the four-target union bound per seed. Compute the smallest B whose counting
ceiling reaches 1% of the group for m=3,4,7,10 at degrees 83 and 131. The
1% line is an **illustrative planning threshold**, not a sufficient condition
for rank, a measured yield, or a necessary condition for every algorithm.
For each threshold also show `C(B+1,2)` pairs and a hypothetical 16-byte per
pair materialization size. That size is an explicit layout scenario, not a
lower bound on an implicit solver or on a compressed representation.

## Decision and validation

Use exact integer combinatorics and binary search. Independently replay the
output with `math.comb`, verify each threshold and its predecessor straddle
the 1% inequality, check the frozen input hashes and the two recorded pair
counts, and check the generic multiset bound by exhaustive sums in a small
cyclic group. Commit the script, independent verifier, machine JSON and a
short interpretation. Do not change the prior m=83 result or impute natural
hits from planted controls.

If the m=83 m=3 upper bound for four natural targets on either base exceeds
1%, the hypothesis is falsified and this audit stops without a no-go claim.
Otherwise classify those archived zero-hit observations as underpowered for
natural-yield discrimination. This is a **conditional counting no-go for those
fixed bases and m=3**, not a no-go for high-arity PDP, larger bases, chosen
targets, precomputation, or ECC2K-130 as a whole. Full-cost S and matched-rho
ratios remain null. The next empirical gate must choose B, arity, target law,
solver and memory budget before sampling, then measure verified useful rank
and full costs; the already open review-gated m10 capacity attempt is separate.
