# Protocol D-1: summation polynomials in non-even coordinates compatible with the Artin–Schreier form

Frozen 2026-10-07, before any instrument is built.  **Stage diagnostic.**
`S`, end-to-end cost and speedup are **unset**.  No work below `2^61` at
`n = 131` is claimed.  Status: **PENDING**.

## Derivation (stated before measuring)

Every coordinate the exotic-coordinates thread searched is a function on
the `x`-line, and
[`../notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md`](../notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md)
§1.1 proves that a degree-2 coordinate is a Möbius frame on that line.  Two
consequences follow for every such coordinate: the summation polynomial
`S_{m+1}` has degree `2^{m−2}` in each variable, and it vanishes on all
`2^{m−1}` sign variants `±P₁ ± … ± P_m = ∓R` of a decomposition, so a root
is an orbit of decompositions and the lift must try them.

On `y² + xy = x³ + ax² + b` over `F_{2^n}` put `u = y/x`.  Then

    u² + u = x + a + b/x²,

so the curve is an Artin–Schreier cover of the `x`-line, `u(−P) = u(P) + 1`,
and the map `u ↦ u² + u` is `F₂`-linear with image the trace-zero
hyperplane.  Two coordinates fall out, neither of them even:

- `u = y/x`, of degree 3, which separates `P` from `−P`;
- `z = x + b/x²`, of degree 3 on the `x`-line, with three abscissae per
  value, whose Weil descent is linear in the trace-zero part because `z`
  is `u² + u` up to the constant `a`.

What might be gained: a summation polynomial in `u` carries no sign fibre,
so the root count per target drops by `2^{m−1}` and the lift is one check;
a base defined by `z ∈ W` for a subspace `W` has about `3·2^{dim W}` points
at the same descent dimension, 1.58 extra bits per summand.  What is
expected to be lost: a degree-3 coordinate has summation polynomials of
degree about `2·3 = 6` per variable instead of 2, so the Weil-descended
system should have per-block Boolean degree 3 where the `x`-system has
degree 2 (bilinear) or 1 (the `(A, P)` form).  Whether the sign collapse
and the base density pay for that degree is the question.  The thread's
own §13 touched the degree-3 function `y` on `j = 0` curves and nowhere
else; neither `u` nor `z` was interpolated.

Filter 2 of the index says why this is the candidate: the Artin–Schreier
map is the only other `F₂`-linear structure on the curve equation, and it
has not been used in a factor base.

## Instrument (Rust, to build in the follow-on PR)

1. Extend `src/cryptanalysis/coordinate_search.rs` so that interpolation
   of summation polynomials from random relation tuples accepts an
   arbitrary rational function `φ` on `E`, not only a Möbius frame in `x`.
   Interpolate `S^φ_3` and `S^φ_4` for `φ ∈ {u, z}` on the registered toy
   curves `icv1-f2m13-t181-515ee569`, `icv1-f2m19-t797-b6cf2467` and
   `icv1-f2m23-t5197-69e76b73`, and on their `a = 1` twins once registered.
   Verify each interpolant on 1,000 fresh relation tuples and 1,000 fresh
   non-relation tuples.
2. Weil-descend `S^φ_3` over subspace bases `W` (random and geometric) in
   `koblitz_symmetrised`'s machinery, count unknowns, monomials and the
   Boolean degree per summand block.
3. Measure first-fall and solving degree with the existing F4 harness
   (`koblitz_groebner`, `matrix_f4_f2`) on 16 target draws per cell, found
   and refuted targets separately, at `n ∈ {13, 17, 19, 23}`.
4. **Control at equal base size.**  Every comparison with the `x`-system
   is at equal `|F|`, not equal subspace dimension: a `z`-base of dimension
   `l` is compared with an `x`-base of dimension `l + 2` or the nearest
   available, and the yield per target is recorded for both so that the
   counting law can be checked rather than assumed.

## Smoke test (to be disclosed, not cited)

One curve, `icv1-f2m13-t181-515ee569`, `φ = u`, 200 relation tuples, to
check that the interpolation converges and that `S^u_3` rejects the sign
variants.  Its numbers are not cited in the result.

## Predictions (pass/fail)

- **A1 (degree law).**  The interpolated `S^φ_3` has degree in `[4, 6]` in
  each variable for both `φ`, and `S^φ_4` degree in `[8, 12]`.
- **A2 (sign).**  `S^u_3` vanishes on every tuple with `P₁ + P₂ + P₃ = O`
  and on none of the sign variants that are not relations; the measured
  collapse factor is 1 against the `x`-line's 4.
- **A3 (descent cost).**  After Weil descent over geometric `W`, the
  per-block Boolean degree of `S^z_3` is at most 3 and the unknown count at
  equal `|F|` is at most 1.5× the `x`-system's.
- **A4 (solve).**  At `n ∈ {17, 19, 23}` and equal `|F|`, the F4 solving
  degree of the `φ`-system is at most the `x`-system's plus one, and its
  median refutation time is within 4× of the `x`-system's.
- **A5 (yield).**  Relations per target at equal `|F|` agree with the
  `x`-system within the binomial interval of 512 targets.

The thread continues only if A1, A2, A3 and A4 all pass; A5 is a control
that must pass for the comparison to be valid.  If A3 or A4 fails the
degree cost exceeds the sign and density gains and the coordinate family
is closed with the measured degrees recorded.

## What a pass would and would not mean

A pass is **engineering**: a smaller constant in the oracle for the same
`|F|`.  It becomes an **advance** only if the solving degree at fixed ratio
stops growing with `n`, which no measurement in the repository has shown
for any coordinate, and which would have to be tested separately on the
ladder of `../notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md`.  The
product law is untouched by any outcome here.

## Stop condition and inadmissible moves

The run is bounded: four toy sizes, two coordinates, 16 draws per cell.
It stops at the first failure of A1 to A4 with the failing table recorded.

Inadmissible: comparing at equal subspace dimension instead of equal
`|F|`; reporting found targets without refuted ones; omitting the group
re-check of every lifted relation; carrying a toy degree to `n = 131` as
anything but an extrapolation marked as such.
