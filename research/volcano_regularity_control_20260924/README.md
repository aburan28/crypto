# Neighboring-curve regularity: a fixed tiny-field control

**Plain English.** We moved the same factor-base points and targets to a
neighboring curve. Every decomposition count stayed the same, and the measured
regularity stayed the same in all 240 paired cases. Rewriting the equations
using the factor-base constraints lowered the regularity number in 160 of the
480 individual systems. This is a useful control: a lower reported degree
can come from the encoding, without changing the point problem.

Classification: **accounting / correctness**. No solver optimization or
ECDLP performance iteration is performed. Full-pipeline cost, rho comparison,
and speedup remain null. This is a new F_103 measurement, not an n=53/n=83
neighbor benchmark and not a replay of the older GF(256) experiment.

## The terminology, with short translations

| Technical term | Plain-English meaning in this report |
|---|---|
| Factor base | The small collection of points we allow as building blocks. |
| Relation yield | How often a target can be expressed using those building blocks. Count distributions and the fraction with at least one solution are distinct. |
| Transported factor base | The same points mapped through the isogeny, so the comparison preserves the combinatorial problem. |
| Polynomial encoding | The equations we choose to express the point problem. Equivalent encodings may behave differently. |
| `d_reg_top` | The first degree at which the highest-degree parts of the chosen equations generate every homogeneous polynomial of that degree. It is an exact encoding-dependent algebraic metric. |
| Solving degree | The largest relevant degree reached by a specified polynomial-solving algorithm under a specified convention. We did not measure it here. |
| Endomorphism stratum | Curves sharing a particular endomorphism order. A lower stratum is not automatically easier for index calculus. |

## Frozen setup and checks

Source: `y²=x³+1` over F_103. Neighbor: `y²=x³+88x+22` over F_103.
The fixed 2-isogeny has kernel `{O,(102,0)}` and sends an affine nonkernel
point `(x,y)` to `(x+3/(x+1), y*(1-3/(x+1)²))`. Both curves have 84 rational
points. We use a cyclic subgroup of order 21, on which transport is injective.
The experiment verifies every subgroup point image and all 441 subgroup
addition pairs against the homomorphism identity. No discrete log is solved.

The source subgroup has ten distinct nonidentity x coordinates. For each
support size `B=2,3,4`, Python's fixed `random.Random(seed).sample` selects
supports for seeds 0 through 7. All 24 resulting supports are distinct.
Including both signs yields factor bases of 4, 6 or 8 points. The exact
supports are frozen in `results.json` so reproduction is not dependent on
an unrecorded sampling choice.

For each support, the full relation census covers all 21 subgroup targets,
including the identity: 504 transported comparisons. Polynomial measurements
use one target per distinct nonidentity x coordinate: 10 targets per support,
240 source-neighbor pairs, and 480 systems in each encoding. Opposite-sign
targets share the x-only polynomial system and are not counted twice.

The raw encoding is `[C_B(X), C_B(Y), S3(X,Y,x_R)]`, where `C_B` is the monic
support polynomial. The reduced encoding replaces `S3` with its remainder
modulo the two support polynomials, preserving the generated ideal. Each
encoding retains its own named regularity value.

Every ordered support pair is checked before and after reduction, and the
polynomial root set is compared with direct rational point addition. The
regularity calculation uses exact finite-field Gaussian elimination on
homogeneous Macaulay matrices, preserving every degree/rank certificate up
to the first full rank. It does not run F4 or F5, so there is no runtime claim.
For the raw systems, the top ideal is `(X^B,Y^B,X²Y²)`. Its regularity threshold
`B+1` is an independent monomial-count check for the implemented sizes.

## Results

| Factor-base points | Paired target cases | Raw `d_reg_top`, both curves | Reduced `d_reg_top`, both curves | Source systems with no relation |
|---:|---:|---:|---:|---:|
| 4 | 80 | 3 | 2 | 49 |
| 6 | 80 | 4 | 4 | 23 |
| 8 | 80 | 5 | 5 | 3 |

| Comparison | Lower on neighbor | Equal | Higher on neighbor |
|---|---:|---:|---:|
| Raw encoding | 0 | 240 | 0 |
| Reduced encoding | 0 | 240 | 0 |

All **504/504** transported decomposition counts matched. The comparison
includes all 75 zero-relation paired polynomial cases; failures to find a
relation were not removed. The 240 pairs share curves, supports, and targets;
they are not 240 statistically independent samples of a curve population.

The preprocessing change affected **160/480** individual systems, exactly
the four-point factor-base cases. In these systems the reduced equation has
a nonzero mixed quadratic highest part, so the threshold becomes 2. That
change does not prove a runtime improvement or a change in relation yield.

## What this says about the earlier GF(256) numbers

The earlier report recorded lower/equal/higher `d_reg_top` in 12/52/32 of
96 transported pairs. It used four Boolean lookup-index variables, while
this control uses two prime-field x-coordinate variables. Its measured
quantity was also tied to its particular generators. The two studies do
not contradict one another: they examine different fields and encodings.
Neither establishes that descending the volcano generally lowers regularity.
The prior GF(256) computations were reviewed from the saved report and were
not independently replayed in this run.

## Reproduction and next evidence needed

Run `python3 research/volcano_regularity_control_20260924/measure.py`.
Outputs: `results.json` (supports, targets, equation coefficients, complete
rank certificates and census counts) and captured `run.log`.

At n=53 and n=83, actual neighbor curves and certified isogenies have not
been constructed by this thread. Their factor-base relation yields and
regularity remain unmeasured. A valid comparison must specify the same
metric, encoding, support cardinality and target sampling, and report
solver work independently of algebraic degree. This control adds a verified
measurement definition; it does not extrapolate its numbers to those sizes.
