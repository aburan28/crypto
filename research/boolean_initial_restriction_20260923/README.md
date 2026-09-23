# One initial affine restriction and transported syndrome blocks

This is a standalone mathematical successor to PR #675. The input remains a
generated public system of Boolean polynomials. No curve, imported target, scalar
recovery, production solver integration or index-calculus/rho campaign is present.
All corresponding unmeasured costs remain null. The predecessor's rejections
remain immutable.

The primary is complete and **all three dramatic gates are rejected (0/18 each)**.
Initial restriction passes 10/18 incremental comparisons, direct block transport
5/18, and transported leaves 0/18. None passes its universal incremental gate.
See [CONCLUSION.md](CONCLUSION.md) for complete costs, retained failures, validation
and the separately measured discovery structure behind the next hypothesis.

## Mechanisms and exact contracts

First reduce the original coefficient row space with all quadratic coordinates
before affine coordinates. The affine tail is exactly the intersection of that
row space with affine polynomials. It includes combinations of arbitrarily many
input rows, not just pairs with matching quadratic support. It is not closure
under multiplying ideal generators.

A fixed linear map from quadratic coefficient vectors into F_2^64 supplies a
cheap sufficient certificate. Each monomial contributes a deterministic one-hot
word and additions are XOR. If the projected rows have full row rank, then the
original quadratic rows also have full row rank: a dependence before projection
would remain a dependence after projection. Thus their span has no nonzero
affine member. Hash collisions can cause false negatives only. Every failed
certificate invokes complete affine-tail extraction; it never authorizes reuse
or deletion by hash equality. Projection rows and reduction XORs, failed probes,
full reductions and all setup are explicitly charged.

For a consistent affine tail of rank r, RREF defines an injective map x=A y+b
with n-r free variables. Every original equation is transformed once. The list
reference multiplies affine polynomial images with explicit Boolean parity
cancellation. The optimized form packs up to 32 equation coefficients into u32
words and applies the same affine map to every equation at once. It includes
constant-linear terms, the Boolean diagonal y_i^2=y_i, cancellations, degree
drops and all newly appearing quadratic coefficients. A recovered model is
always checked against the original equations. Rank zero has no domain-size
benefit; inconsistent affine rows prove UNSAT.
The 32-term source limit applies to original fixtures. Transformed equations
may grow within the complete degree-at-most-two envelope; their terms are never
truncated to the source limit, and transformation growth is charged.

The second mechanism maintains the sixteen syndrome values for the four low
variables. In reflected high-coordinate Gray order, flipping coordinate j
changes a block by a fixed low-affine vector plus a changing scalar word D_j.
Only the scalar parts of lower scheduled differences change. The first discovery
implementation instead updated full difference vectors and was usually slower;
its original source and negative measurements are retained in resource_probe_01.
The factored implementation is in resource_probe_02.

Scalar and native NEON/SSE2 block transport preserve the old assignment order,
first model, assignment/block counts and checksum exactly. The same transport
also runs within the retained 16-variable-leaf search policy. Its prefix and
recovery remain unchanged. Initial affine restriction changes the assignment
order relative to the original coordinates, so comparisons across that policy
check outcomes and recovered solutions, not identical trees. Its two backends
must still match each other in models, logical counters and traces.

## Evidence and frozen boundary

The two discovery probes use only n16/n24, seeds 17 and 937, three families and
three repetitions. Each has 1,008 observations and all 28 arms. They are design
evidence only. In the second probe initial restriction eliminates one coordinate
on both n24 cross-planted inputs, with ratios about 2.42 and 2.53 versus the fastest
retained arm. Small cases and several other methods regress. No holdout claim
follows from these selected discovery results.

The primary protocol fixes 144 systems: all 120 predecessor fixtures and 24
unused holdouts at 12/16/20/24 variables. All 22 predecessor methods, including
their three SIMD treatments, remain in the same binary. Three new scalar/list
controls and three new SIMD treatments bring the roster to 28.

Seven repetitions per fixture give 28,224 observations. The first arm's position
rotates by 7 times the fixture index, then by one each repetition. Across the 36
fixtures at each size, each method appears nine times in each position. This is
balance across a size cohort, not within every fixture. This successor's repetition
count is a new fixed protocol and does not change any predecessor result.

The dramatic gate is unchanged in strength: every n16/20/24, family and regression/
holdout combination must have a paired 95% lower bound above 2.0 against the
pointwise fastest fixed reference. The reference contains all retained methods
and new scalar/list controls. All 18 groups and all 144 cells must complete with
verified results. Separate >1.0 incremental and >1.05 matched-policy gates cannot
replace it. No tuning follows holdout timing. A selected passing treatment would
require unchanged-source confirmation on unused fixtures.

Cold totals include all construction, failed certificates, transformations,
enumeration, recovery, destruction and result checks. Phase timers are exclusive
subsets; uninstrumented phases stay null. Fixture generation and preparation of
the independent search-status reference are outside arm timers but within process
receipts. Whole-worker RSS includes every arm; candidate-specific allocation and
calibrated-operation ratios remain unmeasured. UNKNOWN is censored, never UNSAT.

## Correctness and replay

The 42 Rust tests include the retained 34 checks plus exact affine-tail/recovery
checks on all 16,384 pairs of three-variable quadratic polynomials; all 524,288
quadratic-polynomial/affine-map combinations; dense dimension-changing maps with
all 32 equation bits; every transported block point through n12; every hit lane;
wrapping sums; tiny domains; caps; and full-policy agreement. A deliberate
projection collision verifies the full-reduction fallback. These are producer
checks, not an independent audit or a formal proof.

Use a fresh output directory:

```sh
python3 research/boolean_initial_restriction_20260923/run.py --out /tmp/boolean-initial-replay
```

Every executed run retains sources, protocol, compiler/host data, receipts,
raw measurements, checks and a manifest. Frozen directories must not be changed.
This general generated-system benchmark is not the WDSat/full-curve suite;
its results cannot establish an index-calculus or rho crossover.
