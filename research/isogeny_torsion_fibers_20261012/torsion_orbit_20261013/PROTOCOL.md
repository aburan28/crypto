# Frozen three-torsion orbit coverage screen

2026-10-13, Conductor T-2611. This protocol precedes the new run. The
[four-size cold solve](../full_cost_20261013/REPORT.md) supplies the exact
four non-isomorphic pairs and the reference solver. This is a combinatorial
relation-yield screen, not a measured complete DLP method.

## Question and hypothesis

Does replacing each of twelve subgroup factor-base points by its complete
three-point rational torsion fiber create new, independently useful two-sum
relations after projecting to the prime-order target subgroup? The tempting
prediction is a threefold or greater yield at the same twelve unknown logs.
The exact null prediction is that all three lifts collapse to one point under
the degree-three quotient and under the subgroup projection, so the set of
target subgroup points reached by two sums is unchanged.

## Fixed input and controls

Use `(p,r)=(617,103),(2633,439),(9677,1613),(38561,6427)`, source
`E/F_p:y²=x³+1`, destination `E'/F_p:y²=x³−27`, and
`K={O,(0,1),(0,−1)}`. Select `G` and the first twelve distinct projected
factor-base points exactly as in `full_cost_20261013/runner-v3.rs`: scan
affine points by increasing `x`, take the smaller square root, multiply by
six, and deduplicate by point. Select the first 36 distinct projected
points by the same scan for the size-matched independent-base control.

For each size, enumerate pair sums with repetition from (a) the twelve
projected subgroup points, (b) all 36 torsion translates of those twelve
points, allowing the known torsion correction in the recovered relation,
and (c) the 36 independently selected subgroup points. Count distinct
reachable subgroup elements including identity and divide by `r` to get the
exact hit probability for a uniform subgroup probe. Also report raw pair
entries, quotient classes, and source group additions. The orbit arm must
pay for its 36 representatives; multiplicity is not independent rank.

Verify `T` has exact order three, every orbit has three distinct rational
points, `phi(P+T)=phi(P)` for every selected point and translate, and every
factor-base point has order `r`. Cross-check group-sum support against an
independent scalar-index calculation over all `r` subgroup elements. The
runner is native Rust; retain source/binary hashes, command, stdout, stderr,
and failures. Any failed identity or disagreeing count stops the screen.

The existing generic two-class floor is `sqrt(pi/4)=0.886227` in whole-method
`S`, and matched source negation rho is the reference. This screen does not
time or solve a DLP, so it cannot produce an `S` or speed ratio. It tests only
the proposed three-lift yield mechanism. A null result leaves open the
published torsion-symmetric decomposition-polynomial mechanism over extension
fields, which needs a separately frozen, fully priced oracle experiment.

## Decision

- If orbit support exceeds the twelve-point control after quotient and
  independent scalar replay, inspect the exact new rows and implement full
  relation recovery with charged torsion arithmetic.
- If orbit support equals the twelve-point control, reject a free three-lift
  yield multiplier for this collector. Compare the independent 36-point arm
  to show what extra *distinct* factor-base elements buy, then move the lead
  to extension-field decomposition-oracle cost rather than repeating lifts.
- Preserve every failed check and leave the full-cost question open.
