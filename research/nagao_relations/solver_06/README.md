# Matched symmetric S4 controls

These two controls use exactly the x-subspace factor bases, targets,
three-summand restrictions, first/enumerate modes, solver and all-phase
budgets of solver_04. They do not use the repository's transformed-coordinate
factor base or admit target sums shifted by two-torsion.

## Elementary-symmetric S4

For the four abscissas including r, let e1,e3,e4 be their elementary
symmetric polynomials. On the binary curve with b=1, the repository's
24-term S4 is exactly

```
e3^4 + (e4+1)*e3² + e4*(e4+1)*e1² + e1^4.
```

The validation expands this expression over F2 and checks exact equality
of monomial sets with src/cryptanalysis/koblitz_symmetrised.rs. Thus the
circuit uses a small expression without changing the polynomial. Factor
abscissas remain explicit, constrained and ordered as in the other controls.

## Transformed S4 on the same base

For each explicit x, introduce u and impose `(x+1)*u=1`. Compute `w=u²+u`
and the full sum of u coordinates including the target. Evaluate the
repository's 18-term transformed S4. The validation substitutes
`u=1/(x+1)` and clears denominators, proving exact equality with the same
24-term x polynomial. This is a polynomial identity, not an interpolation
check on a few field elements.

The adapter retains the cost of mapping from the x-subspace base. In these
proper first-d-coordinate subspaces, x=1 is absent. A target with r=1 uses
the elementary chart, preserving that target rather than excluding it.
All emitted x-tuples undergo exact signed lifting to R; failed lifts are
blocked and charged. No R+T result is admitted as a relation to R.

## Outcome

All 96 trials had zero validation errors. Elementary S4 resolved all
first-relation cases at n=5,9,11 and all enumerations at n=5,9; one of eight
n=11 enumerations remained incomplete, despite finding its verified tuples.
The transformed adapter timed out on all larger eleven-bit cases, reflecting
this implementation and its mapping constraints rather than disproving
transformed-coordinate techniques in general.

Elementary S4 is the strongest measured matched Semaev control in this
panel. The direct quadratic function solver beats its eleven-bit diagnostic
cost, but S4 wins at smaller sizes. The stated three-size 20% criterion is
not met. Neither a globally strongest Semaev implementation nor calibrated
operation ratios are established by these adapters; solver/ordering tuning
and additional backends remain possible.

See ../goal_round_20260913.md for all variants in one set of tables. The
raw records, frozen contract, provenance hashes and source publication
mapping are preserved here. This is an engineering comparison with an
unchanged counting boundary, not an attack advance.
