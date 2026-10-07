# P-256 Dickson-aware Riemann--Roch selector, round 296: protocol

Date frozen: 2026-10-07

## Question and boundary

Round 295 constructs the minimum-width independent-log Dickson factor base,
`FB1hc72514a2a8d3`, with 164 negation-folded columns.  Fixed arities 109 and
110 cover the P-256 subgroup, but the degree and storage of a global relation
selector remain unset.  This round tests the odd arity `m=109` through an
exact divisor/Riemann--Roch formulation instead of one Boolean selector per
point or one branch per subset.

The frozen factor-base collection boundary is 13.920747397073491 times the
rho reference, even with a free perfect relation oracle.  The selector can be
promoted only if its measured structured residual degree is at most five and
all original end-to-end gates remain satisfied.  Local quadratic generators
alone do not pass that gate.

## Frozen dependency and identity

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- Round-295 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/dickson-union-result.json`;
- required result SHA-256:
  `9550f2bbaceb9e297fb35c12e79480581ca80bd58b03a8493e86ea7c1eda91a5`;
- factor-base dump:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json`;
- required dump SHA-256:
  `d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8`;
- FB1: `FB1hc72514a2a8d3`, 164 columns / 328 signed points;
- selected components: depth 8 terminal
  `0xd280a56fdc5f82bb38bb4b9f0646a63f916aec63586aaac640427ea9e380d0e8`
  and depth 6 terminal
  `0x42ab26df6d850f09c48b015c21c81fc22c4d0f567a7619fb169a1b6a13bbacd1`.

Hash-check and rebuild both dependencies.  Abort on an identity, point,
column, or curve mismatch.

## Exact common-chain certificate

Let `D_2(z)=z^2-2` and let `u_i=D_{2^i}(x)`, represented by the recurrence
`u_1=x^2-2`, `u_{i+1}=u_i^2-2`.  The union membership condition is exactly

```text
(u_8 - t_8) (u_6 - t_6) = 0.
```

This is a single quadratic terminal constraint on one depth-8 chain, not a
generic disjunction and not `(u_8-t_8)(u_8-D_4(t_6))`, which would admit three
unwanted sibling depth-6 fibres.  For every stored column, evaluate the chain
over P-256 and require exactly the recorded depth-8 or depth-6 membership.
Rebuild both complete component inventories and require that the liftable
roots of the product condition are exactly the 164 stored abscissae with no
overlap or omission.  Record the tempting but overbroad terminal-only
alternative and count its extra liftable columns.

## Compact odd-arity selector

Form the monic, square-free support polynomial

```text
F(X) = product over the 164 stored columns of (X - x_j).
```

Evaluate `F` and `F'` at every stored abscissa, hash its canonical P-256
coefficient stream, and independently rebuild it.  For a target `R` whose
abscissa is not in the factor base, an unordered set of 109 distinct columns
with signs summing to `R` or `-R` is represented by

```text
F = g q
A^2 - (X^3 + aX + b) B^2 = c (X - x_R) g,
```

where `g` is monic of degree 109, `q` is monic of degree 55,
`deg A <= 55`, `deg B <= 53`, and `c` is nonzero.  Saturate by `c` in any
solver that claims exactness.  The coefficient system has 275 unknowns and
275 quadratic equations.  It folds permutations, enforces distinct columns
through square-free `F`, and lets the norm choose signs.  Report its exact
input monomial count and generic Macaulay dimensions at degrees 2 through 7,
but do not call a dense generic matrix a necessary storage cost for a
specialized solver.

Use the public P-256 target

```text
k = SHA256(curve || "/rr-selector-round296/public") mod n,
R = [k]G,
```

replacing zero by one.  Emit `k`, `R`, the support hash, and the complete
system-shape receipt; do not run full P-256 F4 in this round.

## Complete toy transfer

Use `p=1151`, `a=-3`, and P-256 `b mod p`.  The deterministic support pool is
the sorted union of all curve-liftable roots of
`D_32(x)=0` and `D_32(x)=369`, ordered by
`SHA256("rr-selector-round296/support/" || x)` with `x` as canonical decimal.
For `(m,B)=(3,5),(5,8),(7,11),(9,14)`, use the first `B` columns.  These are
support-prefix transfer controls, not claims about P-256 yield.

For each size create:

1. a planted target by hash-selecting `m` distinct columns and signs, retrying
   only if the sum is infinity or its abscissa is in the support; and
2. an unplanted target by hash-ranking all affine toy-curve points, taking the
   first whose abscissa is outside the support.  Preserve its observed class.

Exhaustively enumerate every `2^m * C(B,m)` signed distinct-column tuple,
deduplicate only the target classification (not candidate accounting), replay
every reported relation in the group, and hash the canonical witness stream.
Run the exact  `F=gq`, norm-form system with the repository prime-field F4
engine under grevlex, maximum degree seven, and a 20-second cooperative budget
per target.  A complete algebraic classification must equal exhaustive group
replay.  Timeouts and pairs above the bound are censored, never negatives.
For planted instances, construct an independent `(A,B,g,q,c)` witness by
linear algebra and require direct evaluation of every coefficient equation.

Record variables, equations, input degree and monomials, candidate count,
relations, replayed relations, solving degree, matrix widths, field
operations, wall time, timeout/completeness, and false positives/negatives.
Fit degree, time, field-operation, and memory-width slopes only over at least
three comparable completed sizes; otherwise leave the exponent unset.

## Hypotheses, gates, and stop

- **H1 (chain exactness):** the depth-8/depth-6 terminal product has exactly
  the selected 164 liftable columns; the depth-8-only lifted alternative has
  extra columns.
- **H2 (selector exactness):** every complete toy F4 classification agrees
  with exhaustive replay and every planted witness satisfies all equations.
- **H3 (degree gate):** the largest three complete toy sizes have structured
  solving degree at most five with a non-increasing fitted degree slope.
- **H4 (credible projection):** extrapolation to 275 variables stays below
  `2^50` materialized bytes and prices all branches.  No result based only on
  input degree or a censored run passes.

Promotion additionally requires zero false positives and false negatives on
all complete checked instances, projected cost below `2^103` per usable
relation and below `2^120` for 138,031 independent rows including duplicate,
rank, and sparse-linear-algebra allowance.  Unset is failure.

Stop after the eight toy targets, the exact P-256 identity/shape certificate,
deterministic JSON, and a byte-identical replay.  Attempt no unplanted
full-depth P-256 relation unless every gate passes.  If the degree grows,
times out, or the projection cannot be made without exponential discarded
branches, publish that obstruction as a negative result.

## Deliverables

Implement the certificate and experiment in Rust with unit tests.  Preserve
the deterministic result, isolation receipt, hashes, fitted or censored
exponents, and gate decision.  Update the progress data, scoreboard boundary
facts, comparison dashboard, and generated panel index.  Publish protocol,
implementation, tests, evidence, and result in a PR stacked on Round 295.
