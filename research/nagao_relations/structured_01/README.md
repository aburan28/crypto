# Structured support and pair-invariant scaling experiment

This is the next executable test of the scaling proposals, with a frozen
contract and matched reference/candidate trials. It tests dimensions 8, 9, 10
on the prefix bases and dimensions 8, 10 on F4-linear, Frobenius4-stable bases.
The ambient fields are GF(2^18), GF(2^30); curve coefficients are non-F2
elements of F4. All code and evidence are under `research/`.

## Corrections to the proposed plan

The earlier claim that caching ordinary image setup should substantially
improve these runs was too optimistic. The frozen subfield_01 complete d8
records spend only 0.24455% (n18) and 0.20975% (n30) of total time in setup.
Even removing all of that setup cannot save 20%. Early candidate verification
is the more promising target: it consumes about 71% and 74% respectively.

F4-linear spaces have even F2 dimension, so the structured d9 cell is
mathematically impossible. Fourth-power Frobenius stability and F4 linearity
are separate properties. Both are checked explicitly. The map w↦w²+uw is
F2-linear, and generally is not F4-linear even on an F4-linear V.

If a triple sums to R, conjugating its points by fourth powers gives a
triple summing to R^4. Stable support alone does not allow one to discard
those branches for an arbitrary fixed R. Only the stabilizer of R can be
used for that quotient. Validation saves a concrete counterexample to the
invalid fixed-target orbit collapse; no such pruning is implemented.

## Early support rejection

The baseline constructs a function for a chosen h2∈V and z∈V, then computes
the full residual norm, exclusions, and support before extracting points.
For fixed b, the generator already knows h2 and h1 at a=0. Therefore

```
h1 = a² + a*b + h1_at_zero
u = h2 + z
v = h1 + z*u
```

can be computed first. With I_u=image(w↦w²+uw on V), a necessary support
condition is v∈I_u. We compile the image into independent parity-check masks
of rank n-d+1, and reject at the first failed check. Only survivors undergo
the original full norm, root, exclusion, curve and signed-sum checks.
All mask construction and parity operations are recorded separately from
field operations. This does not change the number of (h2,z) branches.

`validate()` compares this predicate with the original recovery on every
generated candidate for every tiny target, including rejected functions.
It checks membership over all ambient vectors for every tiny image space.
The arbitrary-space reference is independently compared with the frozen
prefix implementation on all tiny targets.

## Structured bases

Let F denote fourth-power Frobenius and ell=n/2, which is odd in this panel.
Factor T^ell+1 over F4 and compute kernels of the resulting polynomials in F.
The implementation chooses a deterministic direct sum of components of the
requested dimension and checks linear independence, closure under F and
multiplication by the F4 generator, and noncontainment in every proper
subfield. It rejects invalid or unavailable dimensions.

This avoids accidentally using only a smaller ambient field, but does not
prove good relation yield or avoid every subgroup concentration. The actual
curve-point support and exact oracle sets are recorded. Every algorithm in
a panel uses the same explicit basis and targets. Comparing different bases
moves the support boundary, so it cannot establish a fixed-support speedup.

## Pair-invariant table: a stronger Semaev baseline

For two distinct nonzero factor abscissas x,y, put u=x+y and v=xy. The
intermediate abscissa q satisfies

```
S3(x,y,q) = u²*q² + v*q + v² + B = 0.
q = (v/u²)*w
w²+w = u²*(1+B/v²).
```

Since both abscissas lift to the curve, the two roots q are the abscissas
of the two pair sums up to negation. All nonzero V inverses are computed
once by batch inversion. We build q→{unordered factor pairs} using these
quadratic S3 equations. To decompose R, loop over the two signs of each
factor P3, compute Q=R-P3, look up x(Q), and recover signs of the pair by
the actual group law. Every returned triple is checked on the curve and
against the exact target. Excluded/repeated abscissas remain excluded.

This is a direct Semaev solver, not a new root-free Nagao formulation.
Its setup remains quadratic in the number of factor abscissas. Query work
is linear in the signed base plus returned table candidates and verification.
A saved table can amortize that setup over a declared batch. The experiment
charges cold setup once plus all eight queries; it does not present this as
a subquadratic cold algorithm or compare it to cold SAT timeouts as a gain.

The independent oracle enumerates signed pair sums using affine curve
formulas and batched denominator inverses. It uses neither S3 equations nor
function coefficients, and is validated against exhaustive signed triples.
It is available only to the harness after solver execution.

## Measurement contract and limits

The 240 cold trials run all five variants in both first and enumeration
modes with identical three-second budgets: the unfiltered hybrid, filtered
hybrid, direct S3 pair table, chained S3 SAT and symmetric S4 SAT. They replay
all eight frozen d8 inputs and add sixteen fresh holdouts on larger or
structured bases. Ten additional cold batches use eight new targets each
and a 180-second all-work budget, with complete enumeration.

The applicable counting ceiling remains
min(1,8*C(M/2,3)/(#E-1)) for M signed factor points. A larger support moves
that ceiling. No algorithmic runtime lower bound or full-DLP floor is
inferred from it. All common-operation speedups, full-DLP S and rho/floor
ratios remain null. Field primitives and newly introduced binary word work
are separate vectors; an API sum alone is not a calibrated comparison.

This follows the parent accounting contract's matching, retention and
verification invariants. The frozen WDSat protocol is for a different curve,
representation and ANF interface. This exploratory family-specific suite
does not pass the full 60-input/repetition or end-to-end promotion gate.
There is one timing repetition, a small holdout sample, no full ECDLP run,
and no exponent fit. Any correctness failure invalidates the variant.

Source and contract are frozen before the campaign. Old sources/results
remain unchanged. Run with Python 3.12 and pycryptosat 5.14.7:

```sh
python research/nagao_relations/structured_01/run.py --validate-only
python research/nagao_relations/structured_01/run.py
```

The campaign refuses to overwrite `raw.jsonl`; reruns require a new sibling
folder. Results, comparisons and the canonical scoreboard accompany the
completed campaign.
