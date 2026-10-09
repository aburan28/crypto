# Scope: G2 — Frobenius-Equivariant F4 for the Barrel Decomposition System

**Date.** 2026-10-08.  **Status.** Scoping with one measured falsifier.
**Origin.** G2 in `docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md`.
**Code.** `examples/g2_equivariant_orbit_probe.rs`; raw output in
`research/g2_equivariant_f4_20261008/g2_orbits_n7_n11.md`.

---

## 1. The claim, restated precisely

Barrel formulation: summand `i` lies in `τ^{k_i}(W)` with a one-hot selector
`s_i[k]`; giving the target its own selector `s_0[k]` (`X_R = Σ_k s_0[k] R^{2^k}`)
makes the whole system invariant under the cyclic shift `σ: k ↦ k + 1` applied
to every block at once, because `S₃` has `F₂` coefficients for Koblitz curves
and `σ` sends `(x₁, x₂, X_R) ↦ (x₁², x₂², X_R²)`.

Faugère–Svartz (ISSAC 2013): for an ideal invariant under an abelian group `G`
of order prime to the characteristic, after extending scalars to a field
containing the characters, F4 splits into one Macaulay reduction per
character, each of dimension `≈ 1/|G|` of the original. For `G = C_n` over
`F₂` the characters live in `F_{2^d}`, `d = ord_n(2)`, and form `(n−1)/d`
Galois orbits plus the trivial character; conjugate blocks need not be
recomputed.

| n | `d = ord_n(2)` | blocks |
|---:|---:|---|
| 31 | 5 | 6 over `F_{2^5}` + trivial |
| 73 | 9 | 8 over `F_{2^9}` + trivial |
| 89 | 11 | 8 over `F_{2^11}` + trivial |
| 127 | 7 | 18 over `F_{2^7}` + trivial |
| **131** | **130** | **1 over `F_{2^130}` + trivial** |

## 2. The structural objection: gauge-fixing

`σ` acts **freely** on the selector blocks and has an obvious slice:
`k₀ = 0`. Every solution `(k₀, k₁, k₂, c)` of the symmetric system is the
`σ^{k₀}`-image of a solution with `k₀ = 0`, so the gauge-fixed system
(`s₀ = e₀`, i.e. `X_R = R`) has exactly the same decompositions, `1/n` the
solutions, and `n` fewer variables. This is what any pipeline solves (the
target is fixed; the Frobenius multiplier is taken at the *relation* level by
orbits), and it is what the ladder's barrel screens solve.

The equivariant method reduces the **symmetric** system, whose Macaulay
matrix is — up to the pure-`c` monomials fixed by `σ` — `n` glued copies of
the gauge-fixed one. Its blocks have dimension `≈ (#orbits) ≈ (symmetric
dimension)/n ≈ gauge-fixed dimension`. So the expected relation is

    dim(equivariant block)  ≈  dim(gauge-fixed Macaulay matrix),

and the equivariant blocks live over `F_{2^d}` instead of `F₂`. If that holds,
G2 **cannot** beat gauge-fixing: it recovers the `n`-fold saving that fixing
`k₀ = 0` gives for free and then pays extension-field arithmetic on top
(`≈ d²/64` word operations per entry; at `n = 131`, `d = 130`).

The only way G2 could still win is if the symmetric formulation admits a
*lower solving degree* or a *smaller monomial set* than the gauge-fixed one —
which adding `n` variables cannot do — or if a second symmetry (summand
permutations, already exploited by symmetrisation) is combined with it.

## 3. Falsifier run today

`g2_equivariant_orbit_probe` builds both systems on toy `K_0` instances in a
normal basis (so `σ` permutes variables *and* equations), with the degree
kept at 3 by auxiliary coordinates `u_i = Σ_k s_i[k] τ^k(y_i)`:

- variables: `s₀ (n, symmetric only), s₁, s₂ (n each), c₁, c₂ (l each), u₁, u₂ (n each)`;
- equations: `n` normal-basis coordinates of `S₃(u₁, u₂, X_R)`, `n` per summand
  for the `u`-constraints, and the one-hot relations per block;
- it checks the symmetric equation set is exactly `σ`-invariant, enumerates
  the Macaulay rows and columns at degrees 2–4 the way `build_macaulay` does,
  counts `σ`-orbits of rows and of columns for the symmetric system, and
  prints them next to the gauge-fixed row and column counts, with ranks.

**Prediction.** `row σ-orbits ≈ gauge-fixed rows` and
`col σ-orbits ≈ gauge-fixed cols` at every degree, within the pure-`c`
monomial count. If instead the orbit counts are several times smaller than
the gauge-fixed counts, G2 has a real dimensional advantage and the
extension-field implementation is worth building.

| n | l | system | vars | eqs | degree | rows | cols | row σ-orbits | col σ-orbits | rank | seconds |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 3 | symmetric | 41 | 87 | 2 | 203 | 672 | 47 | 102 | 179 | 0.0 |
| 7 | 3 | symmetric | 41 | 87 | 3 | 5827 | 10382 | 889 | 1502 | 4406 | 0.1 |
| 7 | 3 | symmetric | 41 | 87 | 4 | 101234 | 107947 | 14570 | 15457 | — | 0.1 |
| 7 | 3 | gauge-fixed | 34 | 65 | 2 | 133 | 455 | 133 | 455 | 118 | 0.0 |
| 7 | 3 | gauge-fixed | 34 | 65 | 3 | 3397 | 6028 | 3397 | 6028 | 2589 | 0.1 |
| 7 | 3 | gauge-fixed | 34 | 65 | 4 | 50708 | 51541 | 50708 | 51541 | — | 0.0 |
| 11 | 3 | symmetric | 61 | 201 | 2 | 373 | 1514 | 53 | 144 | 337 | 0.0 |
| 11 | 3 | symmetric | 61 | 201 | 3 | 17281 | 34606 | 1631 | 3166 | 13595 | 3.3 |
| 11 | 3 | symmetric | 61 | 201 | 4 | 468132 | 539262 | 42672 | 49062 | — | 0.7 |
| 11 | 3 | gauge-fixed | 50 | 145 | 2 | 245 | 1019 | 245 | 1019 | 222 | 0.0 |
| 11 | 3 | gauge-fixed | 50 | 145 | 3 | 9845 | 19536 | 9845 | 19536 | 7764 | 0.9 |
| 11 | 3 | gauge-fixed | 50 | 145 | 4 | 224220 | 246431 | 224220 | 246431 | — | 0.1 |

**Result: the prediction of §2 is refuted, in G2's favour.** The symmetric
equation set is exactly `σ`-invariant (checked), and the symmetric Macaulay
matrix is only 1.5–2× larger than the gauge-fixed one — adding the target's
`n` selector variables to a system of `≈ 5n + 2l` variables multiplies the
multiplier count by `≈ 1.25^{D−2}`, not by `n` — while its `σ`-orbits are
`≈ 1/n` of it. So the equivariant blocks are **4–6× smaller per dimension
than the gauge-fixed matrix** already at `n = 7–11`:

| n | degree | gauge-fixed rows × cols | equivariant block rows × cols (orbits) | per-dimension ratio |
|---:|---:|---|---|---:|
| 7 | 3 | 3 397 × 6 028 | 889 × 1 502 | 3.8–4.0× |
| 7 | 4 | 50 708 × 51 541 | 14 570 × 15 457 | 3.3–3.5× |
| 11 | 3 | 9 845 × 19 536 | 1 631 × 3 166 | 6.0–6.2× |
| 11 | 4 | 224 220 × 246 431 | 42 672 × 49 062 | 5.0–5.3× |

The gauge-fixing argument fails because gauge-fixing removes `n` *solutions*
but only `n` *variables*; the Macaulay matrix is polynomial in the variable
count, so the symmetric system is nowhere near `n` copies of the gauge-fixed
one, and dividing it by `n` wins. With `ρ ≈ 1.5–2` for the symmetric/gauge
size ratio and `pen(d)` the cycles per `F_{2^d}` entry operation against
`1/64` for `F₂`, the dense elimination gain over the gauge-fixed solve is

    gain ≈ n²·d / (ρ³ · 64 · pen(d))      (summing the (n−1)/d conjugate blocks),

which is ≈ 7× at `n = 31` (`d = 5`, table arithmetic), ≈ 70× at `n = 73`
(`d = 9`), ≈ 300× at `n = 131` (`d = 130`, three-word carry-less
multiplies). These are dense-model estimates; the ladder's M4RI kernel
narrows the `F₂` side, and the exponent `ω` of the elimination changes the
powers, so the figures are a target to measure, not a result. **Decision:
build it** (§4), starting with the numerical correctness milestone
`rank(M) = rank(M₀) + d·Σ_orbits rank_{F_{2^d}}(M_j)` on the same toy systems.

### 3.1 Milestone: the block decomposition is exact

`g2_equivariant_orbit_probe --blocks` builds the trivial block over `F₂`
(orbit sums) and one character block per Galois orbit of characters over
`F_{2^d}` (`N_j[r̄, c̄] = Σ_t P[r̄, σ^t c̄] ζ^{−jt}`, fixed rows and columns
dropped), and checks `rank(M) = rank(M₀) + Σ_orbits |orbit| · rank(N_j)`:

| n | d | degree | full rank | trivial block rank | character blocks | identity |
|---:|---:|---:|---:|---:|---|---|
| 7 | 3 | 2 | 179 | 41 | χ₁: 23, χ₃: 23 (orbits of 3) | 41 + 3·23 + 3·23 = 179 ✓ |
| 7 | 3 | 3 | 4 406 | 671 | χ₁: 624, χ₃: 621 | 671 + 3·624 + 3·621 = 4 406 ✓ |
| 11 | 10 | 2 | 337 | 47 | χ₁: 29 (orbit of 10) | 47 + 10·29 = 337 ✓ |
| 11 | 10 | 3 | 13 595 | 1 275 | χ₁: 1 232 | 1 275 + 10·1 232 = 13 595 ✓ |

So the Frobenius-equivariant reduction of the barrel Macaulay matrix is
exact, with the conjugate blocks genuinely free. Timing of the two paths on
the same matrices is in §3.2.

### 3.2 First timing: equivariant blocks vs the full `F₂` rank

Same matrices, same process, two repeats; the control is the full Macaulay
rank through `koblitz_groebner::macaulay_profile` (the crate's `rref_f2`
kernel); the equivariant side is the trivial block (`F₂`, word-parallel) plus
the character blocks with shift-and-add `F_{2^d}` multiplication in the first three rows and log/antilog tables in the last (the block side is then dominated by orbit bookkeeping and matrix assembly, not arithmetic).

| n | d | degree | full `F₂` rank (s) | equivariant blocks (s) | ratio |
|---:|---:|---:|---:|---:|---:|
| 7 | 3 | 3 | 0.103 / 0.105 | 0.020 / 0.016 | 5.2× / 6.6× |
| 11 | 10 | 2 | 0.002 / 0.001 | < 0.001 | ≈ 10× |
| 11 | 10 | 3 | 1.770 / 1.740 | 0.077 / 0.047 | **22.9× / 36.9×** |
| 11 | 10 | 3 (log/antilog tables for `F_{2^10}`) | 2.428 / 2.172 | 0.060 / 0.054 | **40.5× / 40.6×** |

The ratio grows with `n` as the model predicts (blocks shrink like `n/ρ` per
dimension) and this is before any of the obvious block-side speedups
(log/antilog tables for `d ≤ 16`, carry-less multiplication for large `d`,
blocked elimination). Three honest caveats: it is a *rank* comparison, not a
full F4 pass (F4 also needs the reduced rows mapped back — the same
elimination plus an inverse DFT of cost `O(dim²)`); the control kernel is the
plain one, which the ledger's M4RI-blocked kernel beats by ≈ 1.6×; and
`n = 11` is the largest the 64-variable monomial engine admits for this
formulation (`5n + 2l ≤ 64`). The next milestone is a wide-mask engine so the
same measurement runs at `n = 17–31`, where the model predicts 10–100×, and
then integration as a `reduce_system` engine with full-F4 agreement checks.

### 3.3 Wide-mask engine: n = 13–23

`g2_wide_probe` re-implements the formulation over 192-bit monomial masks
(own Boolean polynomials, Weil restriction, Macaulay rows, orbits, blocks,
and a dense `F₂` control) so that `5n + 2l ≤ 192`. It reproduces the
64-variable engine at `n = 11` exactly (`13 595 = 1 275 + 10·1 232`) and
extends the measurement; `l = 3`, degree 3, same process for both sides:

| n | d | vars | eqs | symmetric rows × cols | block rows × cols | trivial rank | character block ranks | blocks (s) | full `F₂` rank | control (s) | ratio | identity |
|---:|---:|---:|---:|---|---|---:|---|---:|---:|---:|---:|---|
| 11 | 10 | 61 | 201 | 17 281 × 34 606 | 1 631 × 3 144 | 1 275 | χ₁×10: 1 232 | 0.05 | 13 595 | 1.4 | 31× | ✓ |
| 13 | 12 | 71 | 276 | 26 404 × 54 752 | 2 092 × 4 210 | 1 655 | χ₁×12: 1 612 | 0.1 | 20 999 | 3.9 | 39× | ✓ |
| 17 | 8 | 91 | 462 | 53 242 × 115 792 | 3 194 × 6 810 | 2 571 | χ₁×8: 2 528, χ₃×8: 2 528 | 0.4 | 43 019 | 37.9 | **105×** | ✓ |
| 19 | 18 | 101 | 573 | 71 677 × 158 558 | 3 835 × 8 344 | 3 107 | χ₁×18: 3 064 | 0.7 | 58 259 | 78.3 | **105×** | ✓ |
| 23 | 11 | 121 | 831 | 120 379 × 273 262 | 5 297 × 11 880 | 4 335 | χ₁×11: 4 292, χ₅×11: 4 292 | 1.8 | 98 759 | 232.5 | **128×** | ✓ |

Reading: the identity `rank(M) = rank(M₀) + Σ|orbit|·rank(N_j)` holds at every
size, the block dimensions are `≈ 1/n` of the symmetric matrix (and 3–6×
below the gauge-fixed matrix), and the time ratio climbs with `n` — 31×,
39×, 105×, 105×, 128× — exactly the direction the cost model predicts, with the
`F₂` control being plain dense elimination of the symmetric matrix in the
same process. The block side (`0.4–0.7 s` including assembly) is still
naive dense elimination over `F_{2^d}` with table arithmetic; nothing has
been optimised on it. Against the *gauge-fixed* matrix (the fair baseline,
`≈ 3–6×` smaller per dimension than the symmetric one) the implied gain is
roughly `ratio / (3–6)³ … ratio / (3–6)²`, i.e. still ≥ 5–10× at `n = 17–19`
and growing.

## 4. If it were to be built anyway

1. Normal-basis Weil restriction so `σ` is a permutation (done in the probe).
2. Orbit representatives for rows and columns; the trivial block = orbit
   sums over `F₂`; the nontrivial block = DFT at a primitive `n`-th root of
   unity `ζ ∈ F_{2^d}`: entry `M̂[r̄, c̄] = Σ_t M[σ^t r, c] ζ^t`.
3. Elimination per block (`F₂` and `F_{2^d}`), rank check
   `rank(M) = rank(M_triv) + d · rank_{F_{2^d}}(M_nontriv)`.
4. Inverse transform of the reduced rows back to Boolean polynomials.
5. Integrate as an alternative `reduce_system` engine in
   `koblitz_groebner`, with the current engine as the control.

Estimated effort: 3–5 days. §3 refuted the §2 prediction, so the build is on;
the next milestone is step 3 (block ranks agree with the full rank) on the
toy systems, then wall-clock against `macaulay_profile` at `n = 31`, dim 16.
