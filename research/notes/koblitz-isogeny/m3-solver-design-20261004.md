# m = 3 decomposition: solver design (2026-10-04)

Status: design only. Nothing here is implemented or measured.

## Current state

`build_decomposition_system` already splits m = 3 into two S₃ links through an
intermediate point e (n free Boolean variables):

S₃(x₁, x₂, e) = 0 and S₃(e, x₃, x_R) = 0.

That gives 3l + n variables and 2n equations, with multidegree (1, …, 1). The pilot
measured about 295 s/call at n = 17 (`presentation-vs-curve-pilot-20261003.md`).
Splitting is therefore not the missing ingredient. What costs is the n intermediate
variables.

## Proposal: symmetrised S₄ over a geometric V

From `subspace-structure-20261004.md`: when V = θ·⟨1, g, …, g^{l−1}⟩, the elementary
symmetric functions of three summands lie in small subspaces:

- e₁ = x₁+x₂+x₃ ∈ θ·⟨g^0..g^{l−1}⟩, of dimension l;
- e₂ ∈ θ²·⟨g^0..g^{2l−2}⟩, of dimension 2l−1;
- e₃ ∈ θ³·⟨g^0..g^{3l−3}⟩, of dimension 3l−2.

This works for any geometric V, not just monomial V; that is what makes the
FGHR-style symmetrisation available beyond ⟨1..z^{l−1}⟩.

1. Precompute S₄(x₁, x₂, x₃, X) = Res_T(S₃(x₁, x₂, T), S₃(x₃, X, T)) over F₂[b]. Then
   rewrite it in e₁, e₂, e₃ (it is symmetric in x₁, x₂, x₃). This is a one-off
   symbolic computation; the char-2 forms are in the literature (Huot thesis), and we
   would re-derive them and check numerically against S₃ chains.
2. Weil descent of S₄(e₁, e₂, e₃, x_R) = 0 with e_k in the subspaces above. That gives
   6l − 3 variables and n equations, with **no intermediate variables**.
3. Solve. Then, for each root, factor T³ + e₁T² + e₂T + e₃ and keep the roots lying
   in V. The existing `lift` stays as the certificate check.

Variable count (n = 17): splitting, 3l + n = 26–32 for l = 3–5; symmetric,
6l − 3 = 15–27. Fewer variables and half the equations. The symmetric system has
higher Boolean degree, though: S₄ is degree 2 in e₃ after symmetrisation, but terms
like e₁e₂e₃ are present. So the win has to be measured, not assumed.

## Tests before any timing claim

- S₄ identity: random (x₁, x₂, x₃, x_R) with x₁+x₂+x₃ = R on the curve make S₄ vanish,
  and random non-decompositions do not.
- Agreement: on n = 11–13, the set of decompositions found by the symmetric solver
  equals the set found by the S₃ chain, and by brute force over V³.
- Null-object control: random V (not geometric). The symmetric descent must refuse
  it, because e₂ and e₃ do not lie in small subspaces.

## Cost to build

Estimated at 2–4 days of engineering: the symbolic S₄, the descent with three
subspace blocks, root-filtering, and the tests. Then a pilot at n = 17, 19 against
the 295 s/call baseline. Go / no-go: ≥ 10× per-call speed-up at n = 17 with identical
decomposition sets.

## Results (2026-10-11): built, tested, measured

This section is an additive amendment; the design above is unchanged. Code is in
`src/cryptanalysis/koblitz_sym_s4.rs`. It is selected with `IcCostOptions::m3_solver`,
or with `KOBLITZ_PILOT_M3=chain|sym|e3` in `koblitz_ic_pilot`.

### The symmetric S₄

`Res_T(S₃(x₁,x₂,T), S₃(x₃,x₄,T))`, rewritten in `e₁, e₂, e₃` with `x = x_R`, is

    x⁴b² + e₃x³b + e₃²(x⁴+b) + e₃³x + e₃⁴ + e₂²x²b + e₂²e₃x³ + e₂⁴x⁴
         + e₁²e₃xb + e₁²e₃²x² + e₁⁴b²

That is 11 terms. The Rust test `symmetric_form_equals_the_resultant` checks the
formula against the resultant at 500 random points over F_{2^17}. `e₁` and `e₂` occur
only as squares and fourth powers, which are F₂-linear. So the Weil descent is
**quadratic**, where the S₃ chain is cubic. The design note had expected degree-3
terms; there are none.

### Two solvers from it

1. **Symmetric Gröbner** (`SymmetricGroebner`). It has 6l − 3 unknowns and n quadratic
   equations. The design's go/no-go missed one thing: in the useful range,
   6l − 3 > n. The system is then underdetermined, with about 2^{6l−3−n} roots in
   (e₁, e₂, e₃). Almost none of those split over V, and the splitting solver pays a
   Gröbner node for each root.
2. **e₃ enumeration** (`SymmetricE3`). With e₃ fixed, S₄ is F₂-affine in the bits of
   (e₁, e₂). So each e₃ ∈ W₃ costs one n × (3l − 1) elimination; Gray-code order
   makes stepping between e₃ values cheap. It is exact and finds every root. Its
   cost is `2^{3l−2}·poly(n)`, which is the exponent class of brute force over V³.
   It is **not** a sub-exponential step.

### Correctness

Test `symmetric_solver_agrees_with_chain_and_brute_force`: n = 11, l = 3, monomial and
geometric V, 40 random targets each. Brute force over multisets of V, the S₃ chain,
the symmetric Gröbner solver and e₃ enumeration all return the same decomposition
sets.

Test `symmetric_e3_m3_solves_soundly_and_refuses_random_v`: n = 13, l = 4, full-group
probes, monomial and geometric V. The IC log equals the BSGS log. A random
(non-geometric) V is refused (`SkipReason::SymmetricNeedsGeometricM3`).

### Per-call cost at n = 17, m = 3, a₂ = 1, a₆ = 1, monomial V, 4 probes

Data: `experiments/koblitz_m3/solver_compare_n17.jsonl`. The table is in the run report.

### Full m = 3 solves

These are closed by the e₃ engine and checked against BSGS
(`experiments/koblitz_m3/full_solves_e3.jsonl`):
- 4 curves × {monomial, geometric g1} at n = 17, l = 5;
- 4 curves × {monomial, geometric g1} at n = 19, l = 6.

All 16 solved, and every log equals both the planted secret and the BSGS log. Wall
time per solve is 3.5–14.5 s at n = 17 and 74–287 s at n = 19. Geometric V costs
about 2–2.5× monomial per call at both sizes; this is observed, not explained here.

### Reading against the go/no-go

The criterion was a ≥ 10× per-call speed-up at n = 17 with identical decomposition
sets. The e₃ engine meets it by a wide margin, and the decomposition sets are
identical on the tested cases. The symmetric Gröbner solver does not meet it at
l ≥ 4.

The speed-up is against this repository's chain solver. It comes from replacing an
algebraic solve with an exhaustive, `2^{3l}`-class enumeration, so it says nothing
about whether a summation-polynomial approach beats generic methods asymptotically.

Next useful step: m = 3 curve-vs-presentation at n = 19–23, which this makes
affordable. At n ≥ 29 the enumeration is 2^{3l−2} with l ≥ 9, which is out of reach.
