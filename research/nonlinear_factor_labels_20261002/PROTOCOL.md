# Preregistered experiment: nonlinear factor-base label gauges

Registered 2026-10-02 before computing any candidate Macaulay profile.

## Question and proposed novelty

The ordinary binary-subspace encoding labels each of eight factor-base
abscissae by its three linear coordinates.  This experiment keeps those eight
field elements unchanged but relabels them with an arbitrary bijection of the
Boolean cube.  Each old coordinate is then the algebraic-normal-form (ANF) of
the new three-bit label.  The same relabeling is applied independently to all
three summands of the chained `S3` system; intermediate-point variables are
unchanged.

This is a target-independent **nonlinear gauge of factor-base labels**, not a
new factor base.  It preserves every decomposition and relation multiplicity
by construction, while a nonlinear Boolean automorphism need not preserve the
degree filtration used by Macaulay/F4 methods.  The hypothesis is that one
gauge class lowers the measured bounded-Macaulay refutation degree, or lowers
the exact matrix work at the same degree, on held-out systems.

The repository contains linear basis changes, selector permutations, symmetry
quotients, and finite invariant encodings, but no search over all nonlinear
label bijections.  A focused Crossref title/abstract search on 2026-10-02 for
combinations of Semaev, Weil descent, Boolean change of variables, factor-base
relabeling, and reversible Boolean circuits found the standard Semaev/Weil-
descent and Groebner studies but no report of this exact search space.  That is
a literature-search result, not proof that the idea has never appeared.

## Frozen search space

There are `8! = 40,320` bijections of three-bit labels.  The affine group
`AGL(3,2)` has `8 * |GL(3,2)| = 1,344` elements.  Right-composition by an
affine relabeling is a degree-preserving change of Boolean variables, so the
bounded-Macaulay resolution degree is invariant within each right coset.  The
primary sweep therefore evaluates the lexicographically least representative
of each of the exactly `40,320 / 1,344 = 30` cosets.  The program must verify
the group and coset counts rather than assume them.

The identity is the baseline.  Candidate selection is lexicographic over the
two discovery systems:

1. fewest unresolved discovery systems through degree 5;
2. lowest maximum resolving degree through degree 6, with unresolved worse
   than any resolution;
3. lowest sum of degree-5 Macaulay columns, then rows;
4. lowest total input terms;
5. lexicographically smallest permutation.

If a candidate has input degree above the requested Macaulay degree, that cell
is explicitly ineligible at that degree.  No lower-degree subsystem is scored.
The selected nonlinear coset representative and the identity are then run on
the two untouched holdouts.  No post-holdout candidate replacement is allowed.

## Frozen systems

All systems are the `m = 3`, `b = 1` chained `S3` construction over
`GF(2^7)`, with `ell = 3`, 16 Boolean unknowns and the repository's verified
irreducible-field constructor.  They replay four committed unsatisfiable rows
from `research/dreg_fixed_surplus_20260923/runs/cell-7-3-7.jsonl`, whose
baseline bounded-Macaulay refutation degree is exactly 6.

| split | draw | basis words | target x | committed baseline |
|---|---:|---|---:|---:|
| discovery | 3 | `[22,30,88]` | 9 | 6 |
| discovery | 4 | `[35,87,26]` | 24 | 6 |
| holdout | 9 | `[120,70,119]` | 34 | 6 |
| holdout | 11 | `[12,69,104]` | 7 | 6 |

The holdout rows may be used by the verifier and final evaluation only after
the discovery winner is serialized.  Their candidate metrics must not affect
selection.

## Correctness and boundaries

For every evaluated representative, the native verifier checks that the map is
a permutation and that ANF substitution commutes with evaluation on every
assignment for both discovery systems.  For baseline, selected candidate and
both holdouts it exhausts all `2^16` assignments.  All four systems must remain
unsatisfiable, matching the committed solution count zero.  It also verifies
that applying the inverse label map restores every original leaf assignment.

The reference is the identity labeling on the same frozen system.  The hard
floor for the primary degree metric is the input-system degree: a resolution
cannot be reported below it by this harness.  A degree reduction counts only if
the candidate resolves correctly at a lower degree than baseline on both
holdouts.  At equal degree, matrix columns/rows are an engineering diagnostic,
not a degree advance.  Wall time is recorded only as a practicality note.

## Success, stop and reporting rules

Primary success requires one preselected nonlinear class to refute both
holdouts at degree at most 5, versus baseline degree 6, with all exhaustive
equivalence checks passing.  Secondary success is a strictly smaller
degree-6 column count on both holdouts with the same resolving degree and
correctness.  Mixed, capped, or unresolved outcomes fail the corresponding
gate and remain in the results.

The sweep stops after all 30 discovery representatives, one frozen winner,
and the baseline/winner holdouts.  It does not tune against the holdouts, run a
full relation collection, or extrapolate a speedup.  This is a bounded algebra
stage diagnostic.  Full IC cost, calibrated operation cost, `S`, the rho ratio,
the required `m = 83` transfer gate, and any result at `GF(2^131)` remain null.

