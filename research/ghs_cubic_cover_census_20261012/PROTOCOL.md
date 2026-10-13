# Protocol: minimal genus of a GHS-type cover over `F_q` for elliptic curves over `F_{q³}`

Frozen 2026-10-12, before `cover_census.rs` was run. Conductor T-2261.
**Structural stage diagnostic on toy fields.** It counts covers by their
Riemann–Hurwitz genus. It constructs no curve over `F_q`, evaluates no
transfer map, solves no discrete logarithm, and claims no speedup.

## Question

For `E` over `K = F_{q³}`, what is the smallest genus of a curve `C` over
`F_q` obtained from `E` by the odd-characteristic GHS construction (Diem;
Momose–Chao), over all choices of the degree-2 map `E → P¹` up to
`PGL2(K)`? How does that change within one step of a `K`-rational
2-isogeny, and over the whole `K`-isogeny class?

Why genus 3 is the only value that matters (derived from published
exponents, not measured here): rho on `E(K)` costs about `q^{3/2}`;
index calculus on a genus-`g` Jacobian over `F_q` costs about
`q^{2−2/g}`, which is `q^{4/3}` at `g = 3` (about `q` for a plane
quartic), `q^{3/2}` at `g = 4`, and more above.

## Construction counted

The map `E → P¹` has a branch set `B` of 4 points. An `F_q`-structure on
`P¹_K` is a semilinear map `τ = A∘σ`, `σ` the `q`-Frobenius,
`A ∈ PGL2(K)` with `A·A^σ·A^{σ²} = 1`. Let `B_i = τ^i(B)`, `i = 0,1,2`.
The compositum of the three conjugate double covers is a `(Z/2)^m` cover
of `P¹`, with

- `m` = rank of `{B_0, B_1, B_2}` under symmetric difference,
- `r = |B_0 ∪ B_1 ∪ B_2|`,
- genus `g = 2^{m−2}(r − 4) + 1` for `m ≥ 2`.

`m = 1` means `B` is `τ`-stable (`E` is a base change or twist on that
structure); it is recorded and excluded, since the norm map then has a
large kernel. For `m ≥ 2` the kernel of conorm-then-norm lies in
`E[2^{m−1}]` (derived; standard argument, not verified here).

For `m = 3` the cover is hyperelliptic if some 2-dimensional subspace `W`
of the span has `|∪ supp(W)| = 3` (a genus-0 quotient by an involution).
For `m = 2` hyperellipticity is not decided by branch data and is
recorded as undetermined.

Cross-check in the tool: the genus formula must equal the sum of the
genera `(s − 2)/2` of all quadratic subcovers (`s` = support size).

## Pre-registered predictions

- **P1 (derived).** Every recorded genus is in `{3, 5, 7, …, 17}`.
- **P2 (derived).** A curve with `#E(K)` odd has minimal genus at least 5,
  and so does its whole isogeny class: both genus-3 types force a
  `K`-rational point of order 2.
- **P3 (open, measured).** The fraction of even-order curves with a
  genus-3 cover, directly and after 2-isogenies. No prediction.

## Falsifiers

P1: any genus outside the set, or a formula/subcover-sum mismatch.
P2: an odd-order curve with a genus-3 cover.

## Frozen inputs

- `q ∈ {5, 7}`: every `K`-isomorphism class of `y² = x³ + Ax + B`
  (non-singular), every cocycle `A` (all `(Q³−Q)/(q³−q)` structures).
- `q = 11`: every cocycle, on the first 120 curves drawn by the fixed-seed
  stream (`seed = 20261012`, tag 11), each with its `K`-rational
  2-isogeny neighbours.
- Fields as in `weak_curve_lab.rs` (first primitive monic modulus).

## Measures

Per curve: `#E(K)`, its largest prime factor, minimal genus, the `(m, r)`
type and hyperelliptic flag at the minimum, the number of structures
attaining it; minimal genus over the curve and its 2-isogeny neighbours;
for `q ∈ {5, 7}` the minimum over the isogeny class (same `#E(K)`).
Per field: histograms, split by parity of `#E(K)` and by cofactor.

## Not measured

The `F_q`-model of `C`, the map `E(K) → Jac(C)(F_q)` and its kernel on
the prime-order subgroup, any index-calculus run, isogenies of degree
other than 2 as explicit steps (the class minimum covers them only as an
existence bound), `n ≠ 3`, characteristic 2 and 3.

## Language

Rust (`cover_census.rs`, dependency-free, `rustc -O`). Python only plots.

## Amendment 2026-10-13 (post-run; adds references and status, changes no frozen input)

**Preregistration status.** The text above was written before the first run
but is a version-controlled file, not an externally timestamped
registration: the first run (`20261013T0305Z`) recorded the Rust source hash
but not this file's hash, and the protocol and the run land in the same pull
request. It is therefore *protocol-versioned at commit*, not independently
preregistered. Run `20261013T0403Z` records the sha256 of this file in its
`env.txt`; that hash covers this amendment too, and the amendment adds no
prediction, frozen input or measure.

**References for the exponents quoted under Question** (published
asymptotics, none measured here):

- Rho on `E(K)`, `#E(K) ≈ q³`, cost `Θ(√#E(K)) = q^{3/2}`: J. M. Pollard,
  "Monte Carlo methods for index computation (mod p)", Math. Comp. 32
  (1978) 918–924, [doi:10.2307/2006496](https://doi.org/10.2307/2006496).
- Genus-`g` Jacobians over `F_q`, `Õ(q^{2−2/g})` (so `q^{4/3}` at `g = 3`,
  `q^{3/2}` at `g = 4`): C. Diem, "On the discrete logarithm problem in
  class groups of curves", Math. Comp. 80 (2011) 443–475,
  [doi:10.1090/S0025-5718-2010-02281-1](https://doi.org/10.1090/S0025-5718-2010-02281-1)
  (main theorem, every fixed `g ≥ 2`); for hyperelliptic curves earlier,
  with the double-large-prime variation: P. Gaudry, E. Thomé,
  N. Thériault, C. Diem, Math. Comp. 76 (2007) 475–492,
  [doi:10.1090/S0025-5718-06-01900-4](https://doi.org/10.1090/S0025-5718-06-01900-4).
- Plane quartics (non-hyperelliptic genus 3), `Õ(q)`: this is a
  *heuristic* running time, not a theorem. C. Diem and E. Thomé, "Index
  calculus in class groups of non-hyperelliptic curves of genus three",
  J. Cryptology 21 (2008) 593–611,
  [doi:10.1007/s00145-007-9014-6](https://doi.org/10.1007/s00145-007-9014-6);
  general degree `d`, `Õ(q^{2−2/(d−2)})`, `d = 4` giving `q`, heuristic:
  C. Diem, "Index calculus in class groups of plane curves of small
  degree", [IACR ePrint 2005/119](https://eprint.iacr.org/2005/119)
  (a proven variant for plane models with a totally split line is in
  Diem's Leipzig preprint "On the discrete logarithm problem for plane
  curves").
- The covering construction counted here: P. Gaudry, F. Hess, N. P. Smart,
  "Constructive and destructive facets of Weil descent on elliptic
  curves", J. Cryptology 15 (2002) 19–46; odd characteristic: C. Diem,
  "The GHS attack in odd characteristic", J. Ramanujan Math. Soc. 18
  (2003) 1–32. (The genus formula and branch-set counting above are this
  protocol's own Riemann–Hurwitz derivation, not quoted from these.)

Titles, venues and exponents were checked on 2026-10-13 against publisher
and ePrint records; page ranges for Pollard, Gaudry–Hess–Smart and the
Diem odd-characteristic paper are from memory and not re-checked.
