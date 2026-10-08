# Novel directions for Koblitz index calculus: relation collection, Gröbner solving, and non-generic structure

**Date:** 2026-10-07
**Status:** hypotheses with falsifiers. Nothing here is measured. Novelty is
**unverified** for every entry (no literature search was run from this
environment); the closest prior art known to the author is named so a reviewer
can check.
**Companion:** `PLAN_IC_ACCOUNTING_FIXES_20261007.md` (accounting fixes, which
should land first so these ideas are measured against the right baseline).
**Conductor task:** T-15.

## 0. Framing: what the structure offers and what is already used

The ECC2K-130 curve `K_0: y² + xy = x³ + 1` over `F_{2^131}` has:

- Frobenius `τ` with `τ² + τ + 2 = 0`, so `End(E) = Z[τ] = O_K`, `K = Q(√−7)`,
  class number 1. `2 = τ·τ̄`; the dual Frobenius `τ̄ = −1 − τ` is the
  Verschiebung, a rational degree-2 isogeny with kernel `⟨T₂⟩ = ⟨(0,1)⟩`, and
  `τ̄² = τ − 1` has kernel `E(F₂) = ⟨T₄⟩`. On x-coordinates `x(τ̄P) = x + 1/x`.
- `E(F_{2^131}) ≅ Z/4 × Z/r`, `r ≈ 2^129` prime, `μ := τ mod r` of order 131.
- The field is `F₂(γ)` with `γ` a primitive 263rd root of unity (2 is a
  quadratic residue mod 263, so `ord₂₆₃(2) = 131` and `γ ∈ F_{2^131}`). The
  type-II optimal normal basis is `β^{2^i}` with `β = γ + γ⁻¹`; squaring is a
  cyclic shift of ONB coordinates; the multiplicative group contains `μ₂₆₃`.
- No nontrivial Frobenius-stable `F₂`-subspace exists (2 is primitive mod 131).

Already used in the repos: signed-Frobenius orbit columns (factor `2n` in
relation count, `(2n)²` in LA), the `T₂` symmetrization (`u = 1/(x+1)`,
`w = u² + u`), the `T₄` symmetry as the `π − 1 = τ̄²` transport, window and
Hamming-weight ONB bases with Frobenius-barrel selectors, S3 chains, F4/M4RI,
CryptoMiniSat, WDSat, msolve, hybrid seed pinning. Not used: everything below.

The only honest certificate of non-genericity is the preprocessing-frontier
quantity from catalogue entry A1-4: an index calculus with base size `B` and
per-target cost `D` is non-generic exactly when `B·D² ≪ N`. Every idea below
should be reported against that quantity, not against wall clock.

---

## 1. Gröbner-basis solving

### G1. F4 trace reuse across same-shape decomposition systems

**Claim.** Every per-target PDP system in a fixed configuration (same `n`, `m`,
window or subspace, same chain) differs only in the constants carrying `x_R`.
The set of useful polynomial multiples F4 selects at each degree is therefore
(almost) target-independent. Record it once, replay it for every target, and
skip symbolic preprocessing and the discovery of reductions-to-zero.

**Mechanism.** Joux–Vitse, *A variant of the F4 algorithm* (CT-RSA 2011),
used exactly this for Semaev systems on `E(F_{p^5})` and reported one to two
orders of magnitude per solve. The autoresearcher kb holds the paper
(`KN-LIT-275`, `KN-LIT-c137bd`); neither pipeline implements it. The repo's
`matrix_f4_f2` already builds explicit Macaulay rows, so the trace is the list
of `(multiplier monomial, equation index, degree step)` triples.

**Falsifier.** At `n=31, dim 16, m=2` (the landed F4 cell, 637 ms median)
and `n=31, m=3` (145 s refutation), record the trace on target A and replay
on 16 fresh targets. Fail if replay produces a different Gröbner basis on any
target, or if the per-target wall drops by less than 3×.

**Expected gain.** 10–100× per solve at fixed shape. Constant factor; moves
the F4 frontier by one or two rungs of `n`, not to 131.

**Prior art.** Joux–Vitse 2011; "F4Remake". Not novel; missing.

### G2. Frobenius-equivariant F4: block-diagonalize the Macaulay matrix by the cyclic shift

**Claim.** In the barrel formulation (summand `i` lies in `τ^{k_i}(W)` with a
one-hot selector for `k_i`), if the target's own Frobenius shift is also made
a free one-hot variable, the whole Boolean system is invariant under the
cyclic group `C₁₃₁` acting simultaneously on every selector. F4 on an ideal
invariant under an abelian group whose characters live in the coefficient
field splits into one Macaulay reduction per character, each `1/|G|` the size.
Here the characters live in `F_{2^130}` and form one Galois orbit plus the
trivial character, so only **two** reductions are needed: a trivial-weight
block over `F₂` and one nontrivial-weight block over `F_{2^130}` (its 129
conjugates come for free).

**Mechanism.** Faugère–Svartz, *Gröbner bases of ideals invariant under a
commutative group: the non-modular case* (ISSAC 2013). The group order 131 is
odd, so the non-modular case applies over `F₂` after extending scalars to
`F_{2^130}`. The discrete Fourier transform over `Z/131` maps the one-hot
selector blocks to weight spaces; monomials acquire weights in `Z/131`; the
ideal's weight-`w` component is handled by its own matrix.

**Falsifier.** Toy window system at `n=19` or `n=23`, `m=2`, barrel selectors,
free target shift. Implement the DFT change of variables over `F_{2^{n−1}}`,
run F4 per weight, and compare the recovered solution set with the plain F4
run. Fail if the solution sets differ, or if the equivariant run is not at
least `n^{ω−1}/c_mul` faster once the Macaulay matrix exceeds `10^4` columns
(`c_mul` is the measured cost ratio of `F_{2^{n−1}}` to `F₂` word
operations; expect `n^{1.8}/c_mul ≈ 10–300×` at `n=131` with `ω ≈ 2.8`).

**Expected gain.** Constant factor on the F4 linear algebra, `10–300×` at
`n=131`, less with sparse elimination. Only useful if the barrel formulation
is the one being solved (it is, in Q1109/Q1111 and the N41 screens).

**Prior art.** Faugère–Svartz 2013 for the method; Faugère–Gaudry–Huot–Renault
2014 use the symmetric group and 2-torsion but treat Frobenius only through
invariant factor bases. Application to the Frobenius shift of Weil-descent
systems: not known to the author.

### G3. Sparse Macaulay matrices and block Wiedemann (BooleanSolve) using ONB sparsity

**Claim.** In the type-II ONB the product of two basis elements has exactly
two nonzero coordinates (complexity `2n − 1`), so each Boolean component of a
product of window elements of dimension `d` has about `2d²/n` quadratic
monomials, not `d²`. The Macaulay matrices of the descended S3 chain are
therefore sparse enough for Wiedemann-style rank and kernel computation
instead of dense M4RI, which is what the 4 GiB cap and the 222-variable wall
at `n=41` are hitting.

**Mechanism.** Bardet–Faugère–Salvy–Spaenlehauer, *On the complexity of
solving quadratic Boolean systems* (2013): hybrid guessing plus a sparse
Macaulay matrix at the solving degree `D`, tested for inconsistency with
block Wiedemann. The repo's `RESEARCH_FFD_PROOF_COMPLEXITY.md` cites the
method; nothing runs it.

**Falsifier.** Build the degree-4 Macaulay matrix of the `n=41` low-power W8
S6 system (222 variables) and report rows, columns, nonzeros per row, and
the Wiedemann cost `cols × nnz`. Fail if `nnz/row > 10⁴` or if the solving
degree needed is 6 or more (then `C(222, 6)` columns are out of reach).

**Expected gain.** Converts the memory wall into a time cost; feasible only
if the solving degree is 4 or 5 at `n=41`. Does not change the exponent.

**Prior art.** BFSS 2013. Not novel; missing.

### G4. Mixed-field F4: Boolean leaves, field-valued intermediates (expected to fail; worth one day)

**Claim.** The chain intermediates cost `(m − 2)·n` Boolean unknowns, the
repo's primary metric. Keep each intermediate `y_k` as a single variable over
`F_{2^n}` and descend only the leaves (`u² = u` for leaf bits). Unknowns drop
from `m·d + (m − 2)·n` to `m·d + (m − 2)`.

**Why it probably fails.** The Boolean descent's degree falls use all
Frobenius conjugates `y^{2^j}` of the intermediates at degree one. Over
`F_{2^n}` with `y` a single variable, `y^{2^j}` is a degree-`2^j` monomial, so
those falls reappear only at degree `2^j`. The solving degree should explode.

**Falsifier.** `n=19, m=3` in Sage over `GF(2^19)` (Singular `std`) versus the
repo's Boolean F4. Fail (as expected) if the mixed system does not finish at
degree ≤ 6 where the Boolean one finishes at 3. Pass would be a surprise and
would say the chaining term is an artifact of the descent, which is worth
knowing either way.

**Prior art.** None known. Novelty unverified.

---

## 2. Relation collection

### A1. `Z[τ]`-multiplier closure of the compact-orbit base

**Claim.** The compact-orbit design uses the multipliers `{±τ^k}` (2n per
column). Every `α ∈ Z[τ]` gives a point `α(P)` with `log α(P) = α(μ)·log P`,
so the multiplier set can be any set `A ⊂ Z[τ]` of small norm without adding
columns. `α(P)` costs one τ-adic expansion (Frobenius is free in the ONB) plus
`wt_τ(α)` additions. With `L = |A|·2n` effective points per column, 4-sum
yield per target scales as `L⁴` and the pair table as `L²`. At fixed table
memory, `K` shrinks by `L/(2n)`, relations needed shrink by the same factor,
LA shrinks by its square, and total rank work drops by `L/(2n)`.

**Mechanism.** Elements of norm `≤ 2^j` in `Z[τ]` number about
`(π/√7)·2^j`; norm `≤ 16` already gives about 19 multipliers up to sign, so
`L/(2n) ≈ 19`. The Verschiebung powers `τ̄^k = 2^k τ^{−k}` are among them and
have explicit degree-`2^k` x-maps (`x ↦ x + 1/x` for `k = 1`), which is how
the `T₂`/`T₄` symmetrizations arise; for `k ≥ 3` the kernel is not rational
so there is no symmetrized polynomial, but the multiplier is still usable as
a base extension.

**Falsifier.** At `n=53, K=600`: add the norm-`≤ 16` multipliers, rebuild
the pair table, and measure probes/relation. Predicted drop: `(L/2n)² ≈ 360×`
at `360×` the memory, or `19×` fewer probes at equal memory with `K = 32`.
Fail if the measured yield is below half the prediction (which would mean
`α(P)` points collide or are not uniformly distributed), or if the relation
matrix with coefficients `α(μ)` loses rank.

**Expected gain.** Constant factor `10–100×` on precompute at fixed memory.
Catalogue A1-2 still applies: the exponent is unchanged.

**Prior art.** Frobenius orbits (GLV; GGMP 2020) are the `{±τ^k}` case.
General small-norm endomorphism closure of a factor base: not known to the
author for binary curves.

### A2. Cyclotomic-sparse factor bases: the second Frobenius-stable family at prime `n`

**Claim.** `x = Σ_{k∈S} γ^k` with `|S| = w` defines a Frobenius-stable set
(Frobenius maps `S ↦ 2S mod 263`) of size up to `C(263, w)`, with an algebraic
membership test of polynomial degree: for `w = 1`, `x²⁶³ = 1`, and since
`263 = 256 + 4 + 2 + 1` this is `x·x²·x⁴·x²⁵⁶ = 1`, a Boolean system of degree
4 in the ONB coordinates. For general `w`, introduce `z_j ∈ μ₂₆₃` with
`x = Σ z_j`. Compare the Hamming-weight ONB base, whose membership is a
cardinality constraint (friendly to SAT, hostile to Gröbner).

**Mechanism.** The field is `F₂(γ)`, so `μ₂₆₃ ⊂ F_{2^131}^×`. The same holds
at `n=11, 23, 83` (`2n+1 ≡ 7 mod 8`, `ord_{2n+1}(2) = n`), which gives a toy
ladder that includes the landed `n=83` rung. The set carries a `263·131`-element
symmetry group (exponent shifts and Frobenius); exponent shift `x ↦ γx` is
not an endomorphism of `E`, so it saves no columns, but it is a symmetry of
the membership set that a solver can exploit the way the barrel exploits
Frobenius.

**Falsifier.** At `n=23`: build the `γ`-weight-`w` base for `w = 2, 3`,
report `|F|`, the subgroup-usable count, and natural 2-sum and 3-sum yields
against a Hamming-weight ONB base of matched size; then compare F4 and
WDSat solving time and solving degree on the two systems for 16 targets.
Fail if the cyclotomic system's solving degree is not lower, or its yield
per column is lower.

**Expected gain.** Unknown. It is a new design axis rather than a predicted
constant; the point is that it is the only other natural Frobenius-stable
family with an algebraic description when no stable subspace exists.

**Prior art.** Redundant cyclotomic representations of `F_{2^n}` are standard
in multiplier design (Wu–Hasan–Blake–Gao 2002; Gauss periods). Using
`μ_{2n+1}`-sparse elements as an index-calculus factor base: not known to the
author. Novelty unverified.

---

## 3. A non-generic approach to test, and the bar it must clear

### N1. Group-algebra lift of the decomposition problem

**Claim.** `F_{2^131}` is a quotient of the group algebra `F₂[Z/263]`
(the "redundant representation"): elements are 263-bit vectors modulo a
132-dimensional kernel, multiplication is cyclic convolution, Frobenius is the
decimation `k ↦ 2k`. For a cyclotomic-sparse base (A2) the unknowns are
exponent sets, and the S3 equation lifts to a convolution identity that must
hold modulo the kernel. If a useful fraction of true solutions satisfied the
lifted identity **exactly** in the group algebra, the decomposition problem
would become an additive-combinatorics problem on `Z/263` (multiset sums of
exponents cancelling in pairs), open to list-sum and sieve methods that do
not exist for the field equation.

**Why it is a long shot.** Exact cancellation imposes 132 extra bits of
constraint; the target `x_R` is dense in every representation. The honest
expectation is that essentially no true solution lifts exactly.

**Falsifier (one afternoon).** At `n=23` (`F₂[Z/47]`), enumerate all S3
solutions with `γ`-weight-2 `x₁, x₂` for 256 targets and count how many
satisfy the lifted identity exactly. If the count is zero or negligible, the
approach is dead and should be recorded as such.

**Prior art.** None known.

### N2. The bar

Catalogue A1-4 gives the generic preprocessing bound `S·T² = Ω̃(N)`. For any
index calculus with factor base `B` and per-target online cost `D`, report
`B·D²/N`. Every rung in the current ledger has `B·D² ≫ N` once `D` is the
mean-target operation count. A value below 1 on any rung, with both arms
measured in operations, would be the first non-generic signal in this
program. None of G1–G4 or A1–A2 is expected to produce one; they are
constant-factor levers whose value is in moving the algebraic frontier by a
few rungs of `n` and in producing measured solving-degree data at sizes
nobody has published.

---

## 4. Suggested order

1. G1 (trace reuse) and A1 (multiplier closure): both reuse landed code and
   have same-day falsifiers.
2. G3 (sparse Macaulay) on the existing `n=41` S6 instance: a measurement,
   not an implementation.
3. A2 (cyclotomic bases) at `n=23` with the `n=83` ladder in view.
4. G2 (equivariant F4) on a toy, then on the barrel systems if G3 shows the
   solving degree is reachable.
5. G4 and N1 as one-day negative controls, recorded either way.
