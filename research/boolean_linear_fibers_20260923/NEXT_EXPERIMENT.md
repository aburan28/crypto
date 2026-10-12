# Prospective partial-syndrome filtering

No implementation or timing claim for this successor is present. Its selection
rule and controls must be frozen before any new holdout timing. Current and
historical failed fixtures remain regressions.

## Structural result from discovery only

The support census exhausts every independent coordinate set on the same twelve
n16/n24 discovery systems, under an explicit 100,000-set cap. All twelve searches
complete after 32–232 sets. The current maximum-cardinality choice leaves no
equation independent of the selected variables on every n24 input and on five
of the six n16 inputs.

Thus a static untouched-row check provides no rejection for the current n24
choices. Alternative sets at the same cardinality sometimes leave one untouched
equation. At n24, sets of size two can leave five or six untouched equations,
and sets of size three can leave two to four. Changing the size increases the
outside domain, so these counts alone cannot establish a net improvement.
`support_census_01` retains every maximizing set and its unaffected-equation mask;
the full setup and runtime of a different selection policy remain unmeasured.

## A different representation to test

A necessary-condition filter can evaluate a fixed subset of at most eight
original equations using an eight-bit syndrome, then evaluate all remaining
equations only for surviving points. Rejecting a nonzero partial syndrome is
exact. A zero partial syndrome never licenses acceptance without the full check.

This representation permits sixteen byte lanes per 128-bit vector instead of
four u32 lanes. That is a representation fact, not a measured fourfold gain.
The cost of constructing coefficient images, maintaining the full-equation
state for later checks, extracting survivors and performing those checks must
be included. An assumed 1/256 survival rate would require independent uniform
equation values; it is only a model and must not replace measured counts.

One bounded formulation uses six low variables, storing 64 partial syndromes
per high-coordinate step. The full syndrome at a surviving low assignment is
reconstructed from the exact high constant, six changing low linear coefficients,
and a precomputed full-word quadratic-offset table. The existing Gray derivative
contract must be re-established for both representations, including cancellations,
constant equations, partial domains and caps.

## Preserve ordering and an honest reference

To match the retained low-four-variable ordering, write its high step as
t=4T+r with 0<=r<4. The two additional low bits must be visited in order

    (r XOR (r>>1)) XOR (2*(T mod 2)),

while coordinates above them follow Gray(T). A scalar independent oracle must
check this identity and the first complete model. All lanes computed in a block
are charged even if an earlier lane produces a model. Partial evaluations,
complete checks and rejected assignments require separate counters.

The successor must include a full-u32 64-point control and a retained SIMD control
without diagnostic hashing. Otherwise a changed batching/checksum schedule could
be misreported as a mathematical gain. Test direct evaluation first; a separate
extension to retained search leaves must preserve that policy's prefix and recovery.
Any dispatcher is itself a new policy, and selecting between existing methods
cannot beat their pointwise cost minimum in an ideal cost model.

Freeze the equation subset, limits, input cohorts, reference roster and complete-
cost gate in advance. Retain all current methods, use unused holdouts, and keep
production, calibrated-operation, full-index-calculus and rho costs null unless
they are actually measured. No asymptotic or cryptanalytic claim is proposed.
