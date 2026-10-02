# Preregistered follow-up: exact affine gauge tuning

Registered 2026-10-02 after nonlinear run 02 and before evaluating any affine
gauge other than the identity.

The nonlinear census falsified a degree reduction but exposed a clean bounded
follow-up.  Affine relabelings preserve the Boolean degree filtration, so they
cannot change the resolution degree.  They can change generator support,
actual Macaulay columns and elimination layout.  This follow-up asks whether a
single target-independent affine gauge reduces exact degree-6 matrix width on
both held-out systems while retaining the baseline degree-6 refutation.

## Frozen two-stage search

Enumerate all 1,344 elements of `AGL(3,2)` in the same native binary.  The
identity remains the reference.

1. On discovery draws 3 and 4 only, rank every affine gauge by summed
   degree-5 columns, then rows, input terms and permutation code.
2. Freeze the top 12 from that proxy without reading holdouts.
3. Run exact degree 6 on those 12.  Discard any gauge that does not refute both
   discovery systems at degree 6.  Select by summed degree-6 columns, then
   rows, degree-5 columns, input terms and permutation code.
4. Serialize that winner, then evaluate identity and winner on untouched
   holdout draws 9 and 11.  No replacement is allowed.

All exhaustive assignment-equivalence and zero-solution checks from the parent
protocol remain mandatory.  The explicit two-million-row, 200,000-column caps
and one Rayon thread are fixed.  No wall-time claim is made.

Primary success is the same degree-6 refutation with at least 5% fewer
degree-6 columns on **each** holdout.  A strict reduction below 5% is a weaker
engineering diagnostic.  A regression on either holdout fails.  The floor is
the invariant resolution degree 6; this study cannot be an algebraic advance.
Full IC cost, the rho ratio, `m = 83`, and `GF(2^131)` remain unmeasured.

Status: completed as run 03.  Both the 5% primary gate and strict-reduction
diagnostic failed; see [RESULTS.md](RESULTS.md#result-2-affine-follow-up).
