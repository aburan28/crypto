# Affine substitution and fast exact Boolean evaluation

This standalone study continues PR #671 on generated public Boolean systems.
It implements the previously prospective affine-substitution contract, then
tests a different exact evaluation strategy after discovery measurements show
that substitution work can exceed the saved search cost. Production solver,
full index-calculus and rho costs remain null.

The primary run is now complete: all four new treatments fail the universal
strongest-reference gate. Direct SIMD passes 16/18 comparisons and the
16-variable-leaf hybrid passes 15/18. Confirmation was stopped before launch.
See [CONCLUSION.md](CONCLUSION.md) for the result, all failed intervals and costs;
the recorded confirmation plan remains unexecuted.

## Affine and Boolean-product contracts

Consistent affine rows give an injective parametrization x=A y+b. Substitution
updates every coefficient with GF(2) parity, including y_i squared=y_i. Recovery
maps compose so returned assignments are checked against the original equations.
No uniqueness assumption is needed. Direct polynomial substitution is compared
with global packed rows, compact free-coordinate rows and whole-word affine
product formulas. The last formulas retain the Boolean diagonal contribution.

Coordinate layouts are keyed by the exact remaining-variable mask. Point
operators use abstract local dimension, local variable position and assigned
value; tests cover reuse under different original-variable label embeddings.
General affine images are rebuilt for the actual map and shared only within
that transformation. Changing coefficients are never accepted using a stale
matrix or hash-only equality. Setup, coefficient growth and recovery are timed.

An additional policy admits f*x_i only when its canonical Boolean degree remains
at most two. For a quadratic row this means every quadratic term contains i.
The direct list reference constructs and checks the product; the packed version
uses the exact common-variable criterion. New monomials within the complete
degree-two envelope are allowed. Nonzero, nonduplicate products are reduced with
the existing rows until an affine consequence appears or rank stops increasing.
This runs only at at most ten remaining coordinates, bounding rank by 56. It is
sound inference, not a claim of complete ideal closure. Every admitted product
and additional reduction is charged and counted.

The implementation also tests branching on the basis alone, removing the second
original-equation representation. List and optimized implementations must match
models, logical counters and traces within each policy. These mechanisms were
selected and optimized only on discovery inputs before the final protocol.

## Syndrome evaluation contract

For at most 32 equations, pack their values into a u32 word F(x). F(x)=0 is
equivalent to satisfying every equation. Split off four low variables y, and
write the remaining variables as h:

\[
F(y,h)=C(h)+\sum_{i<4}L_i(h)y_i+\sum_{i<j<4}Q_{ij}y_iy_j.
\]

High assignments follow reflected Gray order. When step t flips high variable
j=v_2(t), maintain a difference D_j containing its linear coefficient, the fixed
lower-neighbor quadratic contribution (for j>0), and current higher-variable
contributions. Update C by D_j, update only D_i for i<j, and update the four low
linear coefficients by their cross terms with j. Direct-evaluation tests check
every resulting point through twelve variables, not only the returned solution.

Sixteen low assignments are evaluated per high assignment. Scalar, NEON and SSE2
paths must agree on the first satisfying assignment and a wrapping syndrome-sum
checksum. All sixteen points are charged even when an early lane satisfies the
system. The checksum is diagnostic; independent model checks, complete reference
outcomes and direct-evaluation tests provide separate evidence. This follows the
established finite-difference exhaustive-search direction described in
[Fast Exhaustive Search for Polynomial Systems in F2](https://eprint.iacr.org/2010/313.pdf).
No new asymptotic or cryptanalytic claim is made.

Two hybrids retain packed search and cheap affine propagation, then use exact
enumeration when at most 12 or 16 occurring variables remain. Leaf encoding,
mapping and evaluation are charged. A disabled leaf cutoff must preserve the
original packed model, counters and trace. Nonoccurring free variables can be
set to zero. Whole enumeration is bounded to 24 variables; unsupported inputs
take the retained fallback. All selected fixtures fit the enumeration domains.

## Frozen measurement and reference

Seven resource probes retain each discovery step. A separate native CryptoMiniSat
feasibility probe uses the same public equations and checks its answers, but its
Python-wrapper timings are not a same-binary comparison or promotion evidence.

`protocol.json` freezes 120 distinct generated systems: all 96 previously tested
inputs and 24 fresh holdouts, at 12/16/20/24 variables in three families. Twenty-two
methods run in the same binary with twenty-two balanced arm-position repetitions,
giving 58,080 observations if every cell finishes. Unselected discovery prototypes
remain in their immutable probe bundles and have no claimed holdout gain.

The fixed reference includes all thirteen previously retained methods, the packed
method without diagnostic hashing, the affine finalist pair and all scalar
enumeration controls. Each SIMD treatment must exceed the pointwise fastest
reference by 2x with a 95% paired-bootstrap lower bound on every family at
n16/20/24, in both the regression and holdout splits. All eighteen comparisons
and every completed/verified cell are required. The affine finalist is evaluated
against that same roster with itself removed. Matched-backend gains are reported
separately with a 1.05 lower-bound threshold. Comparisons among the three new SIMD
treatments do not move this preregistered reference boundary.

The untraced packed arm is an accounting control. It must preserve the traced
algorithm's model and logical work, and reports a null trace. Including it in the
reference prevents reduced diagnostic hashing from manufacturing the SIMD result.
Enumeration arms report evaluated points, blocks and leaf calls; zero tree nodes
is not zero work. Unmeasured phase costs are null. UNKNOWN is censored and blocks
promotion. A passing selected treatment requires unchanged-source confirmation
with fresh holdout seeds before a robust claim. No tuning follows holdout timing.

All setup, solving, checksum, recovery, destruction and result validation costs
are inside cold arm totals. Fixture generation, preparation of the completed
status reference and record formatting are outside arm timings but inside process
receipts. Whole-worker RSS includes every arm; candidate-specific allocation is
unmeasured. Paired intervals concern these fixed fixtures, not a population-wide
guarantee. SIMD performance is measured on the recorded host; cross-platform CI
validates correctness rather than timing.

## Reproduction

```sh
python3 research/boolean_affine_substitution_20260923/run.py \
  --out /tmp/boolean-affine-replay
python3 -m unittest discover \
  -s research/boolean_affine_substitution_20260923 -p 'test_*.py'
```

After a complete primary pass, `python3 research/boolean_affine_substitution_20260923/confirm.py`
validates the primary manifest and unchanged source, binds the previously recorded
confirmation plan, and runs it into the new `run_02` directory. The plan replays
the two primary holdout seeds as regressions and adds two unused holdout seeds.
A confirmation pass cannot supersede a primary rejection. `summarize.py` binds
completed runs and reports both per-group intervals and individual case ratios;
those case ratios expose regressions that a pooled summary might hide.

Use a new output directory. Sources, protocol, compiler/host metadata, raw samples
and analysis are frozen for each run. Earlier records must not be overwritten.
The WDSat/full-curve suite is inapplicable to this generic generated-system driver;
passing these gates does not establish a full index-calculus crossover.

The native SAT bundle retains its producer exactly as executed. For a new native
feasibility replay, use the separate `native_probe.py --out <NEW_DIRECTORY>`
entry point with the installed matching Python binding. It refuses an existing
output directory; do not rerun a producer inside an immutable evidence bundle.
