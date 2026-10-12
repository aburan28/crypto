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

## Addendum (2026-10-11, later): the geometric-V cost gap, and pinning x₃

### Why geometric V cost 2–2.5× more per call (diagnosed)

Data: `experiments/koblitz_m3/geom_cost_diag.jsonl`; curve a₆ = 1, 20 probes per cell.
- **The algorithmic work is identical.** Monomial and geometric V visit the same
  number of e₃ values and find the same number of S₄ roots per call:
  - n = 17, l = 5: 8192 vs 7790 e₃ values (the geometric probes stop early on a
    relation), 1019 vs 966 roots.
  - n = 19, l = 6: 65536 vs 65536 e₃ values, 16410 vs 16346 roots.
- **The gap came from `WordField::mul`.** Its bit-serial loop runs over the set bits of
  one operand and the high bits of the product. Monomial-basis elements are
  low-degree, so those loops exit early.
- **With PCLMULQDQ the gap closes.** Monomial vs geometric ms/call: 9.6 vs 9.3 at
  n = 17, and 106 vs 123 at n = 19, against 19.4 vs 28.7 and 211 vs 503 with the
  serial loop.
- PCLMULQDQ is now the default; `KOBLITZ_WORDFIELD_SERIAL=1` restores the serial path.
  Both paths are tested against `F2mElement` multiplication.

### `M3Solver::PinX3`: m = 3 for 2^{l+1} linear solves per probe

With s = x₁ + x₂ and p = x₁x₂, S₃(x₁, x₂, X) = s²X² + pX + p² + b. That is F₂-linear in
the bits of (s, p) ∈ W₁ × W₂, so m = 2 decomposition is one elimination plus a check
that t² + st + p splits over V (`decompose_pair_linear`).

For m = 3: for every x₃ ∈ V carrying a point P₃, and both signs ±P₃, decompose the
point R ∓ P₃ with m = 2. That is 2^{l+1} eliminations per probe, against 2^{3l−2} for
e₃ enumeration. Every m = 3 decomposition is reached this way: fix its third summand.

**Tests.**
- `pair_linear_matches_brute_force`: n = 13, l = 4, monomial and geometric V, on targets
  that decompose. The linear solve returns exactly the brute-force pairs.
- `symmetric_e3_m3_solves_soundly_and_refuses_random_v` (extended): on the same seed
  and probes, PinX3 and e₃ enumeration produce the **same number of relations**, and
  both IC solves agree with BSGS.

**Per call, geometric V, 40 probes** (`experiments/koblitz_m3/pin_vs_e3.jsonl`). The
relations found were identical: 2/40 at n = 17 and 1/40 at n = 19.

| n, l | e₃ enumeration | PinX3 |
|---|---|---|
| 17, 5 | 9.7 ms | 0.052 ms |
| 19, 6 | 124 ms | 0.154 ms |

**Full m = 3 solves with PinX3** (`experiments/koblitz_m3/full_solves_pin.jsonl`, curve
a₆ = 1). All four solved, with every log equal to the planted secret and the BSGS log.

| n, a₂, l | monomial V | geometric V |
|---|---|---|
| 23, 1, 8 | 0.5 s, 418 probes, 118 unknowns | 0.7 s, 435 probes, 118 unknowns |
| 31, 0, 10 | 33.5 s, 5612 probes, 525 unknowns | 45.9 s, 5876 probes, 507 unknowns |

**What this is and is not.** The exponent drops from 3l to l per probe: in general,
pinning m − 2 summands costs 2^{(m−2)l} linear solves. This is the standard
"fix all but two summands" reduction, made cheap by the linearity of S₃ in
symmetric coordinates in characteristic 2.

It does not threaten the curves. With l chosen so that the yield is at most about 1
per probe, the whole index calculus costs about 2^{2n/3}·poly(n), which is worse than
Pollard rho's 2^{n/2}. At n = 31, rho needs about 2^{15.5} group operations; the
34–46 s IC solves above are a measurement of this solver, not a competitive attack.

**What it enables:** m = 3 curve-vs-presentation at n = 31, at about 6–8 ms per call.
