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

<!-- TABLE:systems -->
| n | encoding | base | \|F\| (points) | variables | of which auxiliary | generators | max degree |
|--:|:--|:--|--:|--:|--:|--:|--:|
| 7 | SUB | SUB5 | 27 | 10 | 0 | 7 | 2 |
| 7 | MONO | WT2 | 15 | 14 | 0 | 77 | 3 |
| 7 | C | WT2 | 15 | 48 | 34 | 45 | 2 |
| 7 | FC | WT2 | 15 | 42 | 28 | 39 | 3 |
| 7 | QFC | WT2 | 15 | 46 | 32 | 43 | 2 |
| 11 | SUB | SUB7 | 137 | 14 | 0 | 11 | 2 |
| 11 | MONO | WT2 | 67 | 22 | 0 | 341 | 3 |
| 11 | C | WT2 | 67 | 94 | 72 | 89 | 2 |
| 11 | FC | WT2 | 67 | 70 | 48 | 65 | 4 |
| 11 | QFC | WT2 | 67 | 86 | 64 | 81 | 2 |
| 13 | SUB | SUB7 | 125 | 14 | 0 | 13 | 2 |
| 13 | MONO | WT2 | 79 | 26 | 0 | 585 | 3 |
| 13 | C | WT2 | 79 | 114 | 88 | 107 | 2 |
| 13 | FC | WT2 | 79 | 84 | 58 | 77 | 4 |
| 13 | QFC | WT2 | 79 | 106 | 80 | 99 | 2 |
| 17 | SUB | SUB8 | 269 | 16 | 0 | 17 | 2 |
| 17 | MONO | WT2 | 137 | 34 | 0 | 1377 | 3 |
| 17 | C | WT2 | 137 | 172 | 138 | 163 | 2 |
| 17 | FC | WT2 | 137 | 120 | 86 | 111 | 5 |
| 17 | QFC | WT2 | 137 | 150 | 116 | 141 | 2 |
| 19 | SUB | SUB8 | 243 | 16 | 0 | 19 | 2 |
| 19 | MONO | WT2 | 229 | 38 | 0 | 1957 | 3 |
| 19 | C | WT2 | 229 | 196 | 158 | 185 | 2 |
| 19 | FC | WT2 | 229 | 132 | 94 | 121 | 5 |
| 19 | QFC | WT2 | 229 | 174 | 136 | 163 | 2 |
<!-- /TABLE -->

## 4. The table

One table, one unit: `GroebnerSafe` calls per target, every variant a
row, correctness first. `|F|` is the exhaustive oracle's candidate count.
"budget-wild" counts wild verdicts that came from the budget or the
matrix cap rather than from a completed non-linear basis. All targets of
all cells agree with the exhaustive oracle.

<!-- TABLE:results -->
| n | encoding | targets | decomposable | correct | timeouts | calls / target | calls, no-targets | tame | wild | budget-wild | tame depth | XOR words / call | wall / target (s) | \|F\| | calls / \|F\| |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 7 | SUB | 24 | 24 | 24/24 | 0 | 3.62 | — | 1.00 | 2.62 | 0.00 | 2.62 | 6699 | 0.01 | 27 | 0.134 |
| 7 | MONO | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 7.75e+05 | 0.19 | 15 | 0.311 |
| 7 | C | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 7.07e+06 | 1.23 | 15 | 0.311 |
| 7 | FC | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 6.06e+06 | 0.97 | 15 | 0.311 |
| 7 | QFC | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 6.57e+06 | 0.99 | 15 | 0.311 |
| 11 | SUB | 24 | 24 | 24/24 | 0 | 3.25 | — | 1.08 | 2.17 | 0.00 | 2.19 | 1.43e+06 | 0.12 | 137 | 0.024 |
| 11 | MONO | 24 | 24 | 24/24 | 0 | 4.58 | — | 1.00 | 3.58 | 0.00 | 3.58 | 1.14e+08 | 3.98 | 67 | 0.068 |
| 11 | C | 24 | 24 | 24/24 | 0 | 6.75 | — | 1.54 | 5.21 | 4.08 | 4.62 | 2.01e+08 | 35.62 | 67 | 0.101 |
| 11 | FC | 24 | 24 | 24/24 | 0 | 8.25 | — | 2.00 | 6.25 | 5.67 | 5.31 | 1.25e+08 | 19.16 | 67 | 0.123 |
| 11 | QFC | 24 | 24 | 24/24 | 0 | 6.75 | — | 1.54 | 5.21 | 4.08 | 4.62 | 2.02e+08 | 30.88 | 67 | 0.101 |
| 13 | SUB | 24 | 12 | 24/24 | 0 | 1.17 | 1.00 | 1.00 | 0.17 | 0.00 | 0.17 | 1.94e+05 | 0.02 | 125 | 0.009 |
| 13 | MONO | 24 | 14 | 24/24 | 0 | 10.12 | 13.00 | 3.62 | 6.50 | 6.00 | 4.38 | 1.96e+08 | 15.29 | 79 | 0.128 |
| 13 | C | 24 | 14 | 24/24 | 0 | 16.38 | 22.60 | 6.12 | 10.25 | 9.58 | 7.39 | 1.33e+08 | 183 | 79 | 0.207 |
| 13 | FC | 24 | 14 | 24/24 | 0 | 27.96 | 50.20 | 12.00 | 15.96 | 15.92 | 10.38 | 1.03e+08 | 55.30 | 79 | 0.354 |
| 13 | QFC | 24 | 14 | 24/24 | 0 | 12.21 | 16.40 | 4.50 | 7.71 | 7.42 | 5.57 | 1.59e+08 | 81.80 | 79 | 0.155 |
| 17 | SUB | 8 | 1 | 8/8 | 0 | 1.00 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 1.96e+05 | 0.08 | 269 | 0.004 |
| 17 | MONO | 8 | 5 | 8/8 | 0 | 22.25 | 27.00 | 9.25 | 13.00 | 13.00 | 8.76 | 1.69e+08 | 114 | 137 | 0.162 |
| 17 | C | 8 | 5 | 8/8 | 0 | 93.50 | 156 | 40.38 | 53.12 | 44.88 | 16.15 | 1.21e+08 | 633 | 137 | 0.682 |
| 17 | FC | 8 | 5 | 8/8 | 0 | 56.88 | 102 | 25.50 | 31.38 | 30.00 | 14.74 | 5.67e+07 | 231 | 137 | 0.415 |
| 17 | QFC | 8 | 5 | 8/8 | 0 | 68.75 | 118 | 29.75 | 39.00 | 34.12 | 15.35 | 1.52e+08 | 624 | 137 | 0.502 |
| 19 | SUB | 8 | 0 | 8/8 | 0 | 1.00 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 2.75e+05 | 0.02 | 243 | 0.004 |
| 19 | MONO | 8 | 4 | 8/8 | 0 | 33.75 | 41.00 | 14.00 | 19.75 | 17.50 | 11.80 | 1.29e+08 | 284 | 229 | 0.147 |
| 19 | C | 8 | 4 | 8/8 | 0 | 135 | 212 | 62.25 | 73.12 | 66.75 | 18.21 | 1.37e+08 | 891 | 229 | 0.591 |
| 19 | FC | 8 | 4 | 8/8 | 0 | 113 | 175 | 51.50 | 61.38 | 55.75 | 17.81 | 5.18e+07 | 418 | 229 | 0.493 |
| 19 | QFC | 8 | 4 | 8/8 | 0 | 135 | 212 | 62.25 | 73.12 | 66.62 | 18.21 | 1.08e+08 | 1152 | 229 | 0.591 |
<!-- /TABLE -->

Fitted exponents over the prime sizes (least squares of `log₂` mean per
target against `log₂ |F|`; the pre-registered success threshold for the
calls exponent was 0.8):

<!-- TABLE:fits -->
| encoding | sizes fitted | exponent of calls / target in \|F\| | exponent of XOR words / target in \|F\| |
|:--|:--|--:|--:|
| C | 11, 13, 17, 19 | 2.43 | 2.28 |
| FC | 11, 13, 17, 19 | 1.87 | 1.12 |
| QFC | 11, 13, 17, 19 | 2.48 | 2.09 |
| MONO | 11, 13, 17, 19 | 1.51 | 1.64 |
| SUB | 11, 13, 17, 19 | -0.87 | -1.74 |
<!-- /TABLE -->

Budget sensitivity, the same four `n = 11` targets under three per-call
budgets (matrix cap `2^27` words):

<!-- TABLE:tau -->
| encoding | budget (row-word XORs) | targets | correct | calls / target | tame depth | budget-wild calls / target | XOR words / target | wall / target (s) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| FC | 2^26 | 4 | 4/4 | 9.25 | 5.36 | 6.50 | 3.44e+08 | 21.1 |
| FC | 2^28 | 4 | 4/4 | 8.50 | 5.10 | 6.00 | 9.9e+08 | 21.9 |
| FC | 2^30 | 4 | 4/4 | 7.25 | 4.38 | 5.25 | 2.62e+09 | 24.8 |
| MONO | 2^26 | 4 | 4/4 | 6.00 | 3.71 | 4.00 | 3.31e+08 | 3.9 |
| MONO | 2^28 | 4 | 4/4 | 3.00 | 2.00 | 0.00 | 3.3e+08 | 3.5 |
| MONO | 2^30 | 4 | 4/4 | 3.00 | 2.00 | 0.00 | 3.3e+08 | 3.4 |
| QFC | 2^26 | 4 | 4/4 | 6.75 | 4.12 | 4.75 | 3.77e+08 | 22.6 |
| QFC | 2^28 | 4 | 4/4 | 6.00 | 3.71 | 4.00 | 1.26e+09 | 32.9 |
| QFC | 2^30 | 4 | 4/4 | 3.50 | 2.50 | 1.00 | 2.54e+09 | 28.9 |

| encoding | seed | calls at 2^26 | calls at 2^28 | calls at 2^30 | tame depths at 2^26 / 2^28 / 2^30 |
|:--|--:|--:|--:|--:|:--|
| FC | 0 | 7 | 6 | 6 | 6 / 5 / 5 |
| FC | 1 | 13 | 11 | 11 | 6,6,5,4,3,3 / 5,5,4,3,3 / 5,5,4,3,3 |
| FC | 2 | 10 | 10 | 6 | 7,7,6 / 7,7,6 / 5 |
| FC | 3 | 7 | 7 | 6 | 6 / 6 / 5 |
| QFC | 0 | 5 | 5 | 2 | 4 / 4 / 1 |
| QFC | 1 | 10 | 8 | 2 | 5,5,4,3,2 / 4,4,3,2 / 1 |
| QFC | 2 | 6 | 5 | 4 | 5 / 4 / 3 |
| QFC | 3 | 6 | 6 | 6 | 5 / 5 / 5 |
| MONO | 0 | 5 | 1 | 1 | 4 / 0 / 0 |
| MONO | 1 | 8 | 1 | 1 | 4,4,3,2 / 0 / 0 |
| MONO | 2 | 5 | 4 | 4 | 4 / 3 / 3 |
| MONO | 3 | 6 | 6 | 6 | 5 / 5 / 5 |
<!-- /TABLE -->

## 5. Reading it

**H0 stands.** Every target is answered correctly, and the algebra
resolves only after nearly all of one summand's coordinates are fixed:
the mean tame depth of the paper's encodings is within 1–3 of the full
summand at `n ≥ 17`. The oracle is therefore the exhaustive oracle with a
truncated Gröbner computation attached to each candidate, and its call
count grows faster than the candidate count (exponents 1.9–2.5 in `|F_w|` over
`n = 11, 13, 17, 19` against the threshold 0.8; `calls / |F_w|` rises from about 0.1 at
`n = 11` to 0.4–0.7 at `n = 17` and 0.5–0.6 at `n = 19`). The subspace baseline through the same
solver is tame at the root or one level down at every size, which is the
contrast the solver module's documentation predicts: restricting summands
to a subspace keeps everything linear in the unknowns the way a weight
bound does not, however the weight bound is presented.

Three things the table separates:

- **Presentation matters, but not in the useful direction.** Among the
  paper's three ideals FC (fewest variables, highest degree) needs the
  fewest calls at `n = 17` and `n = 19`, where C and QFC (quadratic, more
  variables) fall behind it by 20–65%; the control MONO, `C(n, w+1)` monomial
  generators and no auxiliary variables, resolves 6–8 levels higher than
  any of them. The lifted ideals trade the degree the paper wants to avoid
  for auxiliary variables whose values the truncated basis does not reach
  within a budget. At the paper's own scale (`n = 778` free coordinates)
  the monomial control is not an option, which is exactly why the paper
  builds the lifts; here it only shows what the lifts cost.
- **The budget moves the depth, not the conclusion.** Sixteen times the
  per-call budget lowers the tame depth by one to two levels at `n = 11`.
  A stronger F4 (the paper's Magma runs, 20-minute timeouts, solving
  degree 10 on FC) would resolve higher in the tree; nothing in the sweep
  suggests it would resolve at a depth that does not grow with `n`.
- **The Prange point.** Read against the paper's Proposition 7.2 and
  Table 3, this is the same finding: `C(t̄)` is minimised by fixing
  everything, and each coordinate handed to the algebra costs a call that
  the combinatorial gain of freeing it does not cover. Here the "gain" is
  the difference between `|F_w|` candidates and the tree the solver
  actually walks, and the tree is the candidate set plus overhead.

What the thread establishes positively is structural and belongs in the
decomposition note's §6, where it now is: a normal-basis weight base is a
Frobenius-stable set at prime `n` with a popcount membership test, no
storage, and (by the paper) a bounded-degree membership ideal. The
storage term of §6 was an artefact of thinking of Frobenius-stable sets
as orbit unions to be enumerated; the operation count is unchanged, and
the algebraic oracle that the ideal makes possible does not, on this
evidence, beat the enumeration it replaces. Class: accounting.

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
