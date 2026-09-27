# Exact necessary affine fibers with original-equation recovery

This bounded experiment implements the mathematical contract proposed in
`../boolean_schedule_construction_20260923/NEXT_EXPERIMENT.md`. It uses generated
Boolean quadratic systems at 12, 16, 20 and 24 variables. It has no external-input
interface and does not invoke a curve, scalar-recovery or production solver path.

The earlier parameterized support-envelope constructor is already implemented in
`../boolean_support_envelope_20260922/`. It correctly represents changing
coefficients, cancellations, generator degree drops and newly required multipliers,
but its eager representation failed all 48 performance gates. That construction
result is preserved. This follow-up changes the mathematical scan, not its cache
identity or the meaning of those earlier measurements.

## Exact contract

Write the quadratic system as an equation-coordinate vector `F(y,z)`, where the
first k variables form y and the remaining variables form z. Each monomial has
an equation coefficient in `F2^m`. Let W be the span of the coefficients of
quadratic monomials internal to y. Construct an exact linear map P with kernel W.
Then

    P F(y,z) = P b(z) + sum_i P a_i(z) y_i.

This identity holds at every outside assignment. All original solutions survive.
The reverse implication is false: for `xy+1`, the quotient is zero, but the
original equation requires x=y=1. Quotient dimension is not a count of independent
conditional constraints or a prediction of rejection probability.

The implementation uses XOR elimination on equation-coordinate words. It evaluates
outside assignments with the retained Gray cursor and a necessary SIMD screen.
Current affine columns are updated for every assignment; column rank is recomputed
for every surviving screen query. The affine solve returns a particular solution
and a complete independent kernel basis. All kernel combinations are considered
until an original solution is verified or the fiber is exhausted. Zero columns,
rank drops, constants and a zero quotient retain their ordinary exact meanings.
Caps return UNKNOWN, never UNSAT. Inputs outside the bounded quadratic domain also
return UNKNOWN.

Completeness follows from the column tags: each retained row equals A times its
tag, and each dependent column produces a tag in the kernel of A. A dependent
tag has its own previously unused highest column, so these tags are independent.
There are exactly k-rank(A) of them. Reducing the right-hand side yields either
inconsistency or one particular solution; its translations by all kernel tags
therefore enumerate the entire affine fiber. The necessary screen uses only
linear images of coefficients and cannot reject a consistent affine system.

The three policies fix k=4,5,6 and choose coordinates 0 through k-1. The choice
is frozen before timing and is not selected from measurements. Projection setup,
coefficient updates, screening, elimination, free-variable recovery, original
equation checks and temporary-workspace release all occur inside cold solve time.
Result validation is additionally included in total time. Common fixture generation,
reference-status preparation, serialization and destruction of returned diagnostic
records are outside arm clocks and inside worker process receipts.

## Fixed comparison

`protocol.json` retains all 49 prior solver arms in the same binary. The six
construction diagnostic wrappers from the predecessor do not define different
mathematical policies and are not repeated. The previous 192 fixtures are frozen
in `REFERENCE_FIXTURES.json`; 24 new fixtures use seeds 20261024 and 3145729.
There are 216 systems, 52 methods and eight random/reverse repetitions, for
89,856 observations. No source tuning follows holdout timing.

The reference is the pointwise fastest retained complete solver on the same
fixture and repetition. The dramatic-gain criterion is a paired 95% median-ratio
lower bound above 2.0 in all 18 n16/n20/n24 × family × regression/holdout groups,
with all cells complete and verified. An incremental threshold above 1.0 is
separate. The bootstrap uses 4,000 draws with its seed frozen in the protocol.
Intervals describe this finite grid. A positive result would still require an
unchanged-source confirmation on unused holdouts. Failures and censored outcomes
remain evidence and cannot be removed to create a passing comparison.

This is a generic mathematical solver experiment. Calibrated operation conversion,
cryptanalytic floors, production cost, full index-calculus cost and rho comparisons
are unavailable and stay null. The WDSat/full-IC suite is inapplicable to this
standalone generated-input program; the complete retained Boolean suite and new
holdouts provide the matched control here. Wall-time ratios are finite generic
solver diagnostics, not a full-method or asymptotic claim.

## Run and verify

From the repository root, with Python 3.10+ and stable Rust:

```sh
python3 research/boolean_projected_fibers_20260923/run.py --out research/boolean_projected_fibers_20260923/run_02
python3 -m unittest discover -s research/boolean_projected_fibers_20260923 -p 'test_*.py' -v
```

Every output directory must be new. The runner freezes source and protocol before
compilation, records executable and input hashes, retains per-worker receipts,
and seals the directory after success or failure. Replay reads immutable evidence;
do not run the writing analyzer in an already sealed directory.

The Rust suite includes exhaustive affine-fiber membership and complete recovery,
all pairs of three-variable quadratic polynomials, direct projection identities,
changing rank and constants, capped outcomes, and generated-system truth tables.
The independent Python fixture generator matches all 192 retained fixtures.
