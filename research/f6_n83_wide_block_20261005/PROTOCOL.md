# n=83 F6 algebraic block admission

Registered before source changes or measurement on 2026-10-05. The
current `wide_groebner` Boolean monomial mask holds 128 variables, but
its field structure, symbolic constants, subspace coordinates, and
decoded summands truncate at bit 63. A three-summand chained S3 block
for the 12-dimensional n=83 subspace has `3*12+83=119` variables and
fits that mask. The candidate adds an 83-bit field-structure path and
keeps the existing <=64-bit behavior. This is an **algebraic block
admission**, not a complete eight-summand F6 solver or a speedup claim.

Baseline: main `8323fd352aff97cc75ad10e6089d9e4a33523376`;
`wide_groebner.rs` SHA-256
`168352ca7ac877b7addcb4dec0ba366124690efdda3138b16190569fdca18b3b`.
Freeze K0 curve `icv1-f2m83-tm6151469093347-debefd74` over
`F2[z]/(z^83+z^45+z^2+z+1)`, its public 81-bit subgroup of order
`2417851639230796216685689`, and the standard polynomial subspace
of dimension 12. This source subspace yields 4,057 curve points;
cofactor-four projection yields 4,054 distinct nonidentity subgroup
points. Use deterministic source-point indices for planted block
controls. Preserve all group and field representations and input
hashes in the result.

Acceptance checks, before any performance experiment:

1. A new two-word field table must agree with the existing one at
   degrees 9, 17, and 31, and each degree-83 table product and square
   must agree with the reference `F2mElement` arithmetic. Check high
   coordinate bits explicitly.
2. The m=3 wide system must build with exactly 119 variables. For a
   planted triple from the source subspace, assign all 36 summand bits
   and the 83 intermediate bits from the exact group sum; every
   Boolean equation must evaluate to zero. Independently replay the
   three point sum in the curve group. Any remaining signs must be
   checked by point arithmetic, not assumed from the x equations.
3. The source-to-subgroup bridge must verify `R=4*S` for the planted
   sum `S`, enumerate and check the four rational 4-torsion offsets
   for this exact K0 model, and show the subgroup preimage `U` obeys
   `S=U+T` for exactly one of them. No projected point's x coordinate
   is assumed to remain in the source subspace.
4. Existing wide-backend small-curve correctness tests must pass.
   Unsupported field widths and spent search budgets must remain
   distinguishable from a proved no-decomposition result.

Run only exactness and bounded construction diagnostics here. Preserve
failures, timeouts, memory, source and executable hashes, commands, and
raw logs. Stop the n83 construction if it exceeds 120 seconds; such a
timeout is an inconclusive admission failure. No natural-target yield,
completed n83 relation, full F6 versus F4/F5 comparison, or one-target
IC/rho speedup can follow from a planted three-summand block. The
next gate remains a compositional high-arity solver that returns one
ordinary verified relation with every failed attempt charged. The
incomplete pipeline has no `IC1` candidate ID or measured online cost.

The chained auxiliary-variable approach has prior art in
[Karabina's point-decomposition analysis](https://arxiv.org/abs/1504.02347).
The specific cofactor bridge and wide implementation here are an
engineering admission test; novelty and performance require separate
evidence.
