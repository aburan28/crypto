# Packed residual coefficients at the current Boolean-search frontier

This standalone experiment measures complete generated Boolean solves. The
retained reference includes the fixed-quadratic method from PR #664, not only
the older monomial-list methods. Production solver code and curve workflows
are outside this experiment. Full index-calculus and rho costs remain null.

## Diagnostic and hypothesis

The discovery-only `profile_01` instruments n=16/24 and seeds 17/937. Models,
logical counters and traces match the uninstrumented quadratic reference.
Median phase fractions attribute about 30% to canonical tracing, 21–26% to
affine propagation, and 29–30% to specialization. Compilation is small on these
discovery fixtures. Timers perturb execution; the figures are mechanism
diagnostics and are never used as promotion evidence.

The candidate targets all three substantial terms: local coefficient words
avoid reconstructing quadratic support during tracing, a fixed stack basis
removes affine allocations and scans only occupied pivots, and hash buckets
reduce the number of exact duplicate-equation comparisons.

## Exact representation contract

An equation with q original quadratic terms stores their coefficients in its
lowest q bits, in the original canonical term order. The next n bits are the
linear coefficients in the common variable coordinates; the final bit is the
constant. Require q<=32 and q+n<64. Unsupported inputs take the retained
fixed-quadratic path. The declared generated fixtures have twelve quadratic
terms and n<=24, so every measured candidate stays inside its domain.

Original quadratic coefficients between unassigned variables remain fixed.
For an assignment x_j=b, clear every quadratic coefficient incident to j and
the linear coefficient of j. When b=1, first toggle the constant by the current
linear coefficient of j, and toggle the linear coefficients of its remaining
original quadratic neighbors. Apply multiple assignments in increasing variable
order. Contributions with both endpoints assigned therefore reach the constant
through the same two updates. All coefficient arithmetic uses XOR.

Precompiled touch masks implement coefficient removal. A word's population
count is exactly the number of residual monomials, including a present constant.
Local quadratic bits are enumerated through the original ordered monomial
table; linear bits and the constant follow. This reproduces canonical traces
without scanning every active variable for absent quadratic terms.

Equations have different local quadratic coordinates. Their raw packed words
must never be used as cross-equation equality keys. Deduplication buckets use
the common affine coefficients and the number of residual quadratic terms.
Every bucket match still compares the actual ordered quadratic monomials.
Equal equations necessarily select the same bucket; unequal equations in one
bucket remain distinct. Keep the first equation in original generator order.

The stack affine basis uses lowest-set-bit pivots and eliminates occupied
higher pivots in descending order. It returns the same canonical row space as
the retained reducer. Direct and derived contradictions preserve their separate
trace events; forced assignments, free-variable choices and branching remain
unchanged. No inference is added or removed.

## Reference and rejection conditions fixed before measurement

`protocol.json` fixes n=12/16/20/24, three families and nine balanced arm-position
repetitions. The two discovery seeds recur. The entire prior confirmation grid,
including the two failed n16 groups, becomes a regression split. Two previously
unused holdout seeds supply a fresh split. This is 72 cells, each with all nine
methods, for 5,832 complete-solve observations if every cell finishes.

The dramatic current-frontier gate requires a 95% paired-bootstrap lower bound
above 2.0 against the pointwise fastest of **all eight** retained methods on
every family at n=16/20/24 in both regression and holdout splits. All 18 comparisons
and all grid completions/verification must pass. An additional >1.0 lower-bound
gate reports incremental improvement. A separate cumulative 2x comparison
against the historical seven methods is retained; it cannot replace the stronger
reference. No parameters change after holdout measurements.

Cold arm timing includes compilation of all source tables, state construction,
every recursive decision, propagation, tracing, context destruction and result
validation. Fixtures, completed status-reference preparation and output formatting
are outside arm timing and inside worker receipts. Logical counters stay equal;
they are not calibrated CPU-operation totals. Paired intervals describe repeated
timings on fixed fixtures, not a confidence interval over all Boolean systems.

The mutable packed state occupies 304 bytes at MAX=36, compared with 384 bytes
for the prior state. The candidate adds immutable tables, so this is not a claim
of lower total memory. Whole-worker peak RSS includes all arms. Missing results
or UNKNOWN remain censored with null completion costs and block promotion.

The WDSat/full-curve protocol is inapplicable to this standalone generated-system
driver. Any gain here is scoped to complete generic Boolean solving, with the
same logical search; it is an engineering result and does not change an exponent
or establish a cryptanalytic crossover.

## Reproduction

```sh
python3 research/boolean_quadratic_frontier_20260923/run.py \
  --out research/boolean_quadratic_frontier_20260923/run_01
python3 -m unittest discover \
  -s research/boolean_quadratic_frontier_20260923 -p 'test_*.py'
```

Use a fresh output directory. Each run freezes its worker, old and new state
implementations, algebra kernel, protocol, runner and analyzer, then runs the
correctness tests before launching the grid. Earlier artifacts remain immutable.
