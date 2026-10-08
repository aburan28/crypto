# P-256 algebraic factor-base quotient rigidity, round 304: preregistration

Date preregistered: 2026-10-08

## Objective and fixed boundary

Round 303 shows that a homomorphic same-target cover route needs an
unconstructed genus-165 Jacobian with 164 P-256 factors.  This round evaluates
the other requested parity route: a fundamentally different factor base with
a global algebraic symmetry that reduces its independent logarithm variables.

The fixed comparison is the 164-column base `FB1hc72514a2a8d3` on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`, with rho fixed
at `S=1.3`.  Its current free-perfect-oracle boundary is
`13.920747397073491` times rho.  The registered 17-term base
`FB1h2f8621cda105` remains `394.425280` times rho.

The hypothesis under test is whether a rational algebraic point symmetry can
identify factor-base logarithms without either:

1. becoming an affine scalar transport whose unknown anchor remains to be
   solved; or
2. embedding known scalar labels, reducing the proposed method to a generic
   claw/rho-style search plus its actual construction cost.

This is an exact quotient and identifiability audit, not a new relation
attempt.

## Frozen dependencies

- Round 25 scalar-orbit result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json`,
  SHA-256 `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`.
- Round 297 scalar stabilizer:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/scalar-stabilizer-result.json`,
  SHA-256 `53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982`.
- Round 298 parity boundary:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json`,
  SHA-256 `8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965`.
- Round 303 Jacobian capacity result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/jacobian-row-capacity-result.json`,
  SHA-256 `6b023f47fbb2c57eb79f870f1c25be48e0b4230d974c3823caf76784a095caa0`.

Reject every hash, schema, curve, factor-base, width, scalar-action, stabilizer,
boundary, or imported gate mismatch.

## Rigidity boundary

Every rational map between smooth projective curves extends to a morphism.
For elliptic curves, every morphism `f:E->E` has the form

```text
f(P) = phi(P) + T,
```

where `phi` is a group endomorphism and `T=f(O)`.  On the prime-order P-256
rational group, Round 25 certifies that `phi` acts as a known scalar `[a]`.
Writing `T=[b]G` gives the induced log transport

```text
log_G f(P) = a*log_G(P) + b mod n.
```

If `b` is unknown, it is an anchor variable rather than a quotient relation.
If `b` and every edge coefficient are known, the factor base was supplied
with replayable affine scalar labels.  Round 25 is the registered complete
instance of that favorable case; it must remain classified as a generic
control, not a non-generic index-calculus advance.

A coordinate permutation that does not satisfy an exact affine log identity
does not identify logarithm variables.  It may still define a different
relation system, but it receives no quotient credit without replayable group
equations.

This boundary applies to algebraic self-map transport.  It does not claim that
all non-homomorphic relation mechanisms are impossible.

## Exhaustive affine-function control

Enumerate all `7^7=823,543` functions `f:Z/7Z -> Z/7Z` in lexicographic
value-table order.  Count and retain examples of:

- exact additive functions `f(x+y)=f(x)+f(y)`;
- exact affine-difference functions
  `f(x+y)-f(0)=f(x)+f(y)-2*f(0)`;
- bijective affine-difference functions;
- additive or affine functions not equal to `a*x` or `a*x+b` respectively.

The exact expected counts are seven additive functions, 49 affine functions,
and 42 bijective affine functions, with zero non-scalar/non-affine exceptions.
Replay every admitted table against its recovered coefficients.

This finite control validates the algebra used by the quotient graph; the
elliptic-curve rigidity theorem supplies the algebraic-map classification.

## Exact 164-column transport graphs

Over both `2^61-1` and the P-256 subgroup order, build and rank these frozen
systems in the 164 unknown factor-base logs `x_i`:

1. **independent:** no transport edges, quotient dimension 164;
2. **connected chain:** 163 replayed affine edges
   `x_(i+1)=a_i*x_i+b_i`, rank 163 and one unknown anchor;
3. **dependent cycle:** add a closing edge whose affine composition is the
   identity, rank remains 163 and the anchor remains unknown;
4. **anchor-solving cycle:** add a nondegenerate closing edge, rank 164 and
   quotient dimension zero.

Derive a deterministic witness vector first, then derive every `b_i` from that
vector and the frozen nonzero `a_i`; never choose a right-hand side after
ranking.  Verify all equations, dual-modulus ranks, solution replays, cycle
compositions, and forward/reversed row-order agreement.  Record rows,
operations, materialized bytes, false positives, and false negatives.

The graph result must distinguish:

- **unknown anchor:** connected transport reduces 164 variables to one, but
  Round 25's free-oracle ratio is `1.4465047156952717`, already above rho;
- **known/solved anchor:** the ideal oracle constant is
  `0.964336477130181`, but Round 25's measured direct construction is
  `7.714691817041447` times rho and its labels are scalar-defined;
- **unquotiented geometry:** 164 independent logs and the unchanged
  `13.920747397073491` boundary.

Do not extrapolate Round 25's implementation cost to every possible
scalar-labelled base.  The supported conclusion is that algebraic symmetry
alone supplies no non-affine log quotient and that the known-anchor case still
requires a new complete below-rho algorithm.

## Hypotheses

- **H1:** the complete `Z/7Z` census finds exactly the affine functions and no
  exceptional exact transports.
- **H2:** a connected 164-node affine transport graph has quotient dimension
  one until an independent cycle solves the anchor.
- **H3:** algebraic P-256 self-map transport is affine-scalar on logs; it is
  either an unknown-anchor system or a scalar-labelled generic control.
- **H4:** no tested algebraic quotient changes the executable rho boundary;
  non-homomorphic relation structure remains open and uncredited.

## Promotion gates

Promotion requires all of:

1. exact dependency, function-census, graph-rank, composition, and replay
   checks with zero false positives and false negatives;
2. a P-256 factor base with a replayable quotient not reducible to affine
   scalar labels or a known-anchor generic search;
3. a complete relation/decomposition/recovery implementation for that base;
4. structured residual degree of regularity at most five;
5. complete cost below rho and below `2^120` field/group-equivalent
   operations;
6. cost per usable relation below `2^103`;
7. projected materialized storage below `2^50` bytes; and
8. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.  If
the exact transports are all affine-scalar and neither anchor case passes,
publish the scoped negative result and retain genuinely non-homomorphic
relation mechanisms as the next obligation.

## Deliverables and stop condition

Implement the dependency checker, exhaustive function census, dual-modulus
affine graph certificate, operation accounting, deterministic JSON, unit
tests, isolated canonical/independent replay, typed transfer assessment,
result report, and dashboard update in Rust.  Stop after every exact control
replays or on the first correctness failure.

The transfer workflow is available, but its referenced methodology and
assessment-template resources remain unavailable.  Record that limitation in
the assessment.
