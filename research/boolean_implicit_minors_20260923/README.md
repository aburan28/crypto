# Exact implicit minors: bounded construction diagnostic

This implements the prospective determinant-filter experiment from
`../boolean_projected_fibers_20260923/NEXT_EXPERIMENT.md`. It uses only generated
12- and 16-variable Boolean quadratic systems at the retained discovery seeds
17 and 937. No complete new solver, curve input, imported target or scalar-recovery
path is implemented. Full-IC, production, operation-calibrated and rho costs stay
null. This is an engineering diagnostic, not a full-method performance claim.

## Mathematical contract

An exact equation-coordinate projection removes quadratic monomials internal to
the first four variables y. For every assignment z to the outside variables,
the projected necessary subsystem is `A(z)y=b(z)`. The four columns of A are
affine in z; b is quadratic. Every original solution therefore satisfies

    det([A_R(z) | b_R(z)]) = 0

for any five selected equation coordinates R. In characteristic two determinant
signs coincide. Subset dynamic programming sums every permutation product, with
Boolean multiplication implemented by OR of monomial masks and XOR parity.
The ordinary degree is at most six; Boolean reduction may lower it or produce
the zero polynomial. Changing input coefficients requires a fresh compilation.

The projection basis determines the non-pivot equation coordinates. Four cyclic
windows of five coordinates, starting at offsets 0,1,2,3 in their ascending order,
are fixed before measurement. There is no yield-based row selection. Zero and
duplicate minors are retained. A product cap is a censored construction, not a
small or zero polynomial.

These minors are only necessary conditions. If A loses rank, all chosen minors
can vanish even when `Ay=b` is inconsistent. The tests retain the explicit case
A=0, b nonzero. A full method must still recover complete affine fibers and check
every candidate against the original equations.

## Controls and accounting

Three new filter arms accompany all 52 retained complete solvers in the same binary:

- `numeric_minor4` evaluates the same four determinants directly at every outside
  assignment, using the existing Gray coefficient updates.
- `symbolic_minor4` compiles their algebraic normal forms and converts them into
  complete truth masks using the Boolean Mobius transform.
- `affine_filter` constructs the stronger exact necessary affine-consistency mask,
  including rank-deficient systems.

The different predicates are labelled explicitly. Filter-only timings are not
complete-solver costs. Each filter clock includes projection, coefficient
construction, allocations, complete mask construction and temporary-workspace
release. Total time additionally includes checking every returned mask against
the oracle. Fixture generation, independent oracle construction, exhaustive
original-system checking, diagnostic serialization, repeat-signature checking
and returned-output destruction are outside arm clocks and inside worker receipts.
Instrumented product counters and cap checks remain in the symbolic compiler's
measured cost; no observer or phase time is subtracted.

There are 12 fixtures, 55 arms and eight random/reverse repetitions: 5,280
observations. This discovery-only cost diagnostic does not integrate a new solver
or launch holdouts. Its frozen decision determines whether such a candidate is
justified. It does not promote a performance change on the basis of discovery data.

For a full solver that retains this measured complete-mask prepass, even free
recovery cannot reduce its cost below the prepass. The optimistic observed ratio
is the fastest retained complete-reference total divided by the symbolic arm's
build time. A 95% upper interval below two blocks the proposed universal dramatic
claim for this implementation. A lower interval above two in all six groups would
only justify a complete candidate with new holdouts. The bound does not apply to
different compilers, representations, instrumentation or early-termination designs.

No hardware-independent operation conversion or applicable cryptanalytic floor is
measured. The WDSat/full-IC suite is inapplicable to this standalone mathematical
preprocessor. Its matched Boolean controls price every implemented stage; they do
not supply the missing recovery cost.

## Validation and execution

The Rust oracle evaluates the original systems directly at all affine corner
points, computes numeric determinants by row elimination, and checks consistency
with a separate affine row solver. The Python replay independently constructs
original equation truth tables as integers, checks determinants by column rank,
and verifies every filter bit. Both check that every original solution survives.
The full original truth tables cover 417,792 assignments across the 12 fixtures.

The Rust suite covers symbolic/numeric determinant equality, explicit permutation
parity, Boolean cancellation and degree drops, zero and duplicate rows, packed
truth transforms through word boundaries, rank-deficient counterexamples, caps,
and all retained solver tests. Joint rejection and pairwise overlap are recorded
exactly; marginal rejection rates are never added or treated as independent.

Run into a new directory only, with Python 3.10+ and stable Rust:

```sh
python3 research/boolean_implicit_minors_20260923/run.py --out research/boolean_implicit_minors_20260923/run_02
python3 -m unittest discover -s research/boolean_implicit_minors_20260923 -p 'test_*.py' -v
```

Source and protocol are frozen before compilation. Executable hashes, raw output,
per-worker receipts, correctness receipts and the result are sealed by a manifest,
including a failed execution if one occurs. Use the read-only verifier for replay;
never run the writing analyzer inside an already sealed directory.
