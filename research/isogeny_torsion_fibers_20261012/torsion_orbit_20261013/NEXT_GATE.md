# Next test: torsion-symmetric decomposition oracle

Pending protocol, dependent on the exact [orbit support screen](REPORT.md).
The preceding screen rejected a free three-lift hit multiplier in a finite
prime-field pair table. The surviving lead is the explicit nonfactorable
coordinate `chi=(y+1)/(y−1)` from the [equivariance check](REPORT.md),
with `chi(P+T)=1/(1−chi(P))`. Its order-three action may change the
decomposition-polynomial representation and solver cost over extension
fields. This is an implementation and experiment obligation, not a
claimed advantage.

The one-target generic two-class counting floor is
`S_floor=sqrt(pi/4)=0.886227`, with `S=all counted operations/sqrt(r)`.
The reference is a measured, same-target negation rho run, with a
six-class automorphism-aware rho control where applicable. These
boundaries are fixed before oracle tuning; a polynomial-solver or
relation-stage ratio alone cannot cross them.

## Paired mechanism test

Construct an explicit extension-field curve `E/F_(p^n)` with rational
order-three `T`, a verified degree-three quotient `phi:E→E'`, and a public
prime-order subgroup with `gcd(3,r)=1`. Freeze `p`, `n`, irreducible field
polynomial, curve and subgroup identities, map, target coordinates, factor
base, seeds and solver caps *before* running. Verify `chi` and its
order-three Möbius action over the full extension field, including poles;
the coordinate must remain nonfactorable through `phi`. Start with a bounded small
field that allows exhaustive decomposition-oracle answer checking, then add
at least four growing sizes before any scaling interpretation.

On each same source target compare (1) direct source Semaev or equivalent
decomposition over an independently chosen factor base, (2) source
torsion-orbit canonicalisation, (3) the nonfactorable `chi` coordinate
with its order-three action, and (4) invariant `x∘phi` decomposition on
the quotient curve as a control. Match the number of *distinct projected
subgroup elements* and separately record equal-memory and equal-oracle-call
comparisons. Include an isomorphism-only model arm if destination coordinates
change. Each oracle must return independently verified source points or a
checked absence certificate within its finite search domain; count
decomposition hits, independent row rank, relation recovery and target
descent, not only polynomial monomials or Gröbner time.

Price field operations, polynomial-system construction, SAT/F4 work,
factor-base setup, map evaluations, torsion corrections, linear algebra,
descent, verification, failed attempts, peak memory and amortization by the
exact reuse count. Run same-target matched rho, a source A/A control, an
isolated L2-or-better session and independent host replay before claiming
any wall or portable `S` advantage. A cheaper oracle with unchanged or worse
whole-solve cost is a stage result. A bounded solver failure is an
implementation or resource outcome, not a mathematical obstruction.

The current repository has native finite-field and Semaev components, but
no completed three-torsion extension-field equivariant oracle wired into
ecbench. Selecting the field, constructing and checking that oracle, and
running the paired solve remain open work. The earlier [full-cost next gate](../full_cost_20261013/NEXT_GATE.md)
for four-size prime-field iso-only attribution and independent-host pricing
also remains pending; this test does not replace it.
