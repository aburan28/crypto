# Exact image-space support and linear root recovery

This adaptive successor replaces the dominant cubic remainder evaluation
with the already-proved quadratic image condition. It runs after solver_07
has completed, on its exact targets, modes and budgets. These targets were
held out from earlier rounds but are now development data for this successor.

For each nonzero u in V, construct the F2-linear map T_u(w)=w²+uw on V.
Its kernel is {0,u}; hence its image has dimension d-1. Gaussian elimination
stores both a basis for this image and preimages of the basis elements.
Every field-vector elimination update and all image setup are counted.

For a generated function with H(z)=0, define u=h2+z and v=h1+zu. Then

```
H(X)=(X+z)*(X²+uX+v).
```

The remaining quadratic splits into distinct roots inside V exactly when
v lies in the image of T_u. Reducing v by the stored basis simultaneously
recovers w in V with T_u(w)=v. Its two roots are w and w+u. Thus the support
test also supplies the remaining abscissas, with no cubic modular squaring,
root scan, or additional inverse. Distinctness from z, nonzero/target
exclusions, on-curve recovery and the exact signed sum are still checked.

This is an exact replacement for H divides L_V in the restricted hybrid
chart, not a relaxed heuristic. It does not change the factor base, arity,
set of genuine relations, or candidate function generation. It conditions
on one abscissa and is not the pure root-free SAT formulation. The current
search still has O(2^(2d)) conditioned branches before field-operation costs;
keeping d fixed while increasing n does not establish an exponent gain.

Validation compares the returned support roots with the previous full cubic
predicate on every generated candidate for all 43 affine five-bit targets,
including rejections. Complete sets match the independent group oracle.
Every larger complete result is separately checked using group-law pair
lookup. See the combined ../scaling_23_29.md and frozen JSON records for
counts, timing limits and the unresolved full-cost goal.

A next scaling test must increase d as well as n, and use new held-out
seeds. Precomputation remains per-instance in these experiments; changing
to amortized setup requires equivalent treatment of the Semaev controls.
