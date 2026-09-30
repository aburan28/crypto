# Hamming ideals and an ISD-like Gröbner oracle for point decomposition

**Experiment:** [`research/hamming_ideal_pdp_20260930/`](../../hamming_ideal_pdp_20260930/)
(`PROTOCOL.md`, the solver crate `hamming_pdp/`, frozen runs under
`results/main/`, `RESULT.md`).
**Related:** [`RESEARCH_ECC2K130_DECOMPOSITION.md`](../ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)
(the cost surface this oracle would have to move; its §6 is what the
construction below touches),
[`RESEARCH_SEMAEV_DECOMPOSITION.md`](RESEARCH_SEMAEV_DECOMPOSITION.md)
(the subspace oracle), [`RESEARCH_GROEBNER_F4.md`](RESEARCH_GROEBNER_F4.md)
(the repository's own engines, which this thread deliberately does not use:
their monomials are 64-bit masks and the lifted Hamming ideals need more
variables than that).

**Source paper.** R. La Scala, M. Marchesin, S. K. Tiwari, *Hamming ideals
and Gröbner bases for ISD-like syndrome decoding* (2026). Its objects:

- the **Hamming variety** `H_t = {v ∈ F_2^n : wt(v) = t}` is cut out, modulo
  the field equations, by `e_{2^k}(v) = t_k` for the binary digits `t_k` of
  `t`, because Lucas' theorem makes the Boolean elementary symmetric function
  `e_{2^k}` the `k`-th digit of the weight (their Theorem 2.3); only digits
  up to `L = ⌈log₂(max(t, n − t) + 1)⌉ − 1` are needed (Theorem 2.5);
- since `e_{2^k}` has `C(n, 2^k)` monomials, three **lifted** presentations
  with auxiliary variables along a balanced binary split of the coordinates:
  **C-Hamming** (the convolution `e_d^{(a,b)} = Σ_k e_k^{(a,m)} e_{d−k}^{(m+1,b)}`,
  `O(n L)` variables, quadratic), **FC-Hamming** (only power-of-two ESFs kept,
  every other `e_j` replaced by its Lucas factorisation `Π_{h ∈ bits(j)} e_{2^h}`,
  `O(n)` variables, degree up to `L + 1`), and **QFC-Hamming** (the Lucas
  products given their own variables, `O(n L)` variables, quadratic);
- **GBDecode**, an ISD-like decoder that fixes only `r ≤ k` coordinates of
  an information set and solves the rest algebraically, and **MultiSolve**,
  which calls a degree-truncated Gröbner computation (`GroebnerSafe`, with a
  timeout) at every node of a binary assignment tree and branches when the
  basis is not linear (`OracleT`: every node is tried).

Their finding on the Classic McEliece Category 1 parameters is negative
for the algebra: fixing an entire information set (`r = k`, Prange) is the
cheapest configuration, and each coordinate left free costs Gröbner calls
that outweigh the `1.24`-bit combinatorial gain of `r = k − 10`.

## 0. Why this touches the ECDLP thread at all

The decomposition note's §6 records that at prime `n` no Frobenius-stable
factor base of subspace shape exists between dimension 1 and `n − 1`, and
that a Frobenius-stable **set** — a union of `π`-orbits — has "no low-degree
membership polynomial", so it must be **materialised** (`2^l` stored points)
and decomposed by a lookup instead of by algebra. The subspace oracle's
module documentation says the same thing from its side: restricting each
summand to a subspace keeps the Weil-descended `S₃` quadratic, "the property
the linearised-polynomial factor base has and a weight-bounded one does
not".

The paper's construction bears on exactly that sentence. In a **normal
basis** `{α^{2^i}}` of `F_{2^n}`, Frobenius is a cyclic shift of
coordinates, so

```text
    F_w = { P ∈ E(F_{2^n}) : wt_NB(x(P)) ≤ w }
```

is `π`-stable at every `n`, prime or not, at every size
`Σ_{j ≤ w} C(n, j)`. Membership is a popcount, so the `2^l` storage of §6
is gone. And by the paper, `wt_NB(x) ≤ w` is a system of `O(n)` (FC) or
`O(n log n)` (C, QFC) equations of bounded degree in the coordinates plus
auxiliary variables: a Frobenius-stable base with an algebraic membership
ideal, which §6 said prime `n` forbids. Two of §6's three claims survive
unchanged (the operation count of the orbit-collapsed enumeration, and the
product law); the storage claim does not, and whether the algebraic
oracle is *useful* is the question this thread measures.

Semaev's `S₃` in normal-basis coordinates stays quadratic for `m = 2`
(squaring is a cyclic shift of the coordinate variables, the product
`x₁x₂` is bilinear), so the whole decomposition system is `S₃` (`n`
quadratic equations in `2n` coordinate variables) plus two copies of the
Hamming ideal.

## 1. Boundary, unit, target

Stated in `PROTOCOL.md` before any cell ran; repeated here.

- **Reference:** the exhaustive oracle, the Prange analogue of this
  setting: for each `P₁ ∈ F` test whether `x(R − P₁) ∈ F`. `|F|` group
  subtractions and `|F|` membership tests per target. Both oracles may
  quotient by Frobenius, so that common factor `n` is applied to neither.
- **Floor:** none new. The free-oracle floor (§5.3 of the decomposition
  note) and the orbit-collapse cost (§6) hold for any base of this shape.
- **Unit:** `GroebnerSafe` calls per target, the paper's unit, next to the
  candidate count `|F|`; elimination row-word XORs as the second unit,
  uncalibrated against curve operations. Wall time is a practicality note.
  This is a **stage diagnostic** throughout: no `S`, no rho ratio, no
  end-to-end figure.
- **Falsification target (H1):** correctness on every target; mean tame
  depth at most `n − 3` on the weight cells at `n ≥ 13`; and the exponent
  of mean calls per target in `|F_w|`, fitted over the four prime sizes
  `11, 13, 17, 19`, below `0.8` for at least one of C, FC, QFC. Otherwise
  the paper's negative result transfers (H0).
- **Expected class:** accounting — a storage term of §6 is removed, no
  operation count moves.

## 2. What was built

`research/hamming_ideal_pdp_20260930/hamming_pdp/`, a standalone crate
with no dependencies, so that the encodings are not constrained by the
repository engines' 64-variable monomials:

| module | contents |
|:--|:--|
| `gf2n.rs` | `F_{2^n}` on `u64`, the `find_irreducible` reduction-polynomial convention of `scripts/ecc2k130_point_decomposition.py`, the smallest normal element, both change-of-basis matrices, a self-check that Frobenius is a cyclic shift of normal coordinates |
| `curve.rs` | `K_0 : y² + xy = x³ + 1`, point count checked against the Koblitz recurrence at startup, lifting by half-trace, addition, negation, cofactor clearing |
| `semaev.rs` | symbolic field elements with Boolean-polynomial coordinates; Weil descent of `S₃(x₁, x₂, x(R))` in any coordinate basis |
| `hamming.rs` | the balanced tree `T_n`, the C-, FC- and QFC-Hamming ideals exactly as in Theorems 3.5, 4.5 and Section 5 of the paper (leaf variables identified with the coordinates, product variables introduced on demand), the digit constraints for `wt ≤ w`, and a monomial-ideal control (`MONO`: every product of `w + 1` distinct coordinates) |
| `f4.rs` | a degree-truncated Boolean F4: pairs of lcm-degree `≤ d` only, Gebauer–Möller lcm criteria (the product criterion is false in the Boolean quotient and is not used), the field pairs `(x_v + 1)·tail(f)`, retirement of elements a newer leading monomial divides, eager propagation of linear elements, a budget in row-word XORs and a matrix-size cap in place of the paper's wall-clock timeout; `GroebnerSafe`'s verdict is `1` (no solution), linear (an affine solution set, enumerated and checked against the original equations), or wild |
| `multisolve.rs` | `MultiSolve` with `OracleT`: a `GroebnerSafe` call at every node, branching on the summand coordinates in order (`b = 0` first), the residual-weight bound, forced zeros once a summand reaches weight `w`, stop at the first verified decomposition; a direct check at a fully assigned leaf so that no budget failure can drop a solution |
| `selftest.rs` | 300 random systems per seed in 4–8 variables: every F4 verdict is checked against brute force (an inconsistent verdict must have no solutions, a linear verdict must enumerate exactly the solutions, a wild verdict at full degree must have a non-affine solution set) |
| `main.rs` | factor bases, targets, the exhaustive oracle, the JSONL records |

Every verdict of the solver is sound independently of what the truncation
achieves: `1` is reported only when it is in the ideal, a linear verdict is
enumerated and every point checked against the original system and lifted
to points on the curve, and a wild verdict only branches. Completeness is
guaranteed by the leaf check. The self-test passes on 1,200 systems.

The subspace baseline (`SUB`) is the repository's oracle shape — a random
subspace of the matched dimension, polynomial-basis coordinates, the same
`S₃` — run through the same `MultiSolve` so that every row is in one unit.

## 3. The systems

Sizes at the root, before any assignment, from the frozen manifests
(`S₃` equations included, field equations implicit in the Boolean ring).
`L = ⌊log₂ n⌋` digits are carried, the minimum that determines any weight
in `[0, n]`.

SYSTEMS_TABLE

## 4. The table

RESULTS_TABLE

## 5. Reading it

RESULTS_TEXT

## 6. What this does not settle

- **`m ≥ 3`.** Only `m = 2` was run. With intermediate points the chained
  `S₃` is cubic and the Hamming ideal is unchanged, so the algebra can only
  get harder; the exhaustive oracle's `|F|^{m−1}` gets worse faster, which
  is the one direction in which the comparison could move.
- **The engine.** The tame depth depends on how much of the truncated basis
  the engine reaches within its budget: a better F4 (or Magma's) resolves
  higher in the tree. The budget-sensitivity cell in `RESULT.md` bounds
  that dependence at `n = 11`; it does not remove it. A faster run of this
  same engine, with the same `xor_words` and the same tame/wild verdict, is
  recorded in
  [`ENGINE_SPEED.md`](../../hamming_ideal_pdp_20260930/ENGINE_SPEED.md).
  That change is engineering of the stage diagnostic's wall time. It does
  not resolve a higher basis, and it is not a scoreboard row. An `n = 53`
  root call outside the frozen set is in that note: the matrix cap stops it,
  and on that call the widened engine is slower than the one it replaced.
  `n = 83` does not fit the `u64` field.
- **The `r`-sweep of `GBDecode`.** The paper fixes `r` coordinates and
  enumerates `u` of weight `t̄` outside the algebra. Here `MultiSolve` does
  the fixing itself and the tame depth reports the `r` it needed; the
  outer combinatorial loop of `GBDecode` (`C(t̄)` of their Section 7) was
  not run separately because it would multiply, not change, the per-node
  figures.
- **Larger `n`.** Nothing here is a measurement above `n = 19`; the
  exponent is a fit over four prime sizes and is marked as such.
