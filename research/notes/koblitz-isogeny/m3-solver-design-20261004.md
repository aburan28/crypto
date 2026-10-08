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
