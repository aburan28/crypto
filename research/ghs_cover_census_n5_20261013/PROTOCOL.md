# Protocol: GHS cover census at larger `q` for `n = 3`, and for `n = 5`

Written 2026-10-13 before any panel run. Conductor T-2467. Extends
`research/ghs_cubic_cover_census_20261012/` (PR #1759, "the merged census").
**Structural stage diagnostic on toy fields.** It counts covers by their
Riemann–Hurwitz genus. It constructs no curve over `F_q`, evaluates no
transfer map, solves no discrete logarithm and claims no speedup.

Status: protocol-versioned at commit (this file's sha256 is recorded in every
run's `env.txt`); not independently preregistered. What was run before this
file was frozen is listed under "Runs before freezing".

## Question

For `E: y² = x³ + a₂x² + a₄x + a₆` over `K = F_{q^n}` and the degree-2 map
`x: E → P¹` with branch set `B = x(E[2]) = {∞} ∪ {roots of the cubic}`,
what is the smallest genus of the `(Z/2)^m` compositum of the `n` conjugate
double covers, over every `F_q`-structure `τ` on the `x`-line?

1. `n = 3`: does the merged census's pattern (class minimum 13 / 5 / 3 for
   `#E(K)` odd / 2 mod 4 / 0 mod 4) hold at `q = 11, 13, 17`?
2. `n = 5`: what is the minimal genus and what governs it?

## Construction counted (as in the merged census, generalised to `n`)

`τ = A∘σ`, `σ` the `q`-Frobenius, `A ∈ PGL2(K)` with
`A·A^σ···A^{σ^{n−1}}` scalar. `B_i = τ^i B` for `i = 0..n−1`;
`m` = rank of `{B_i}` under symmetric difference, `r = |∪ B_i|`,
`g = 2^{m−2}(r − 4) + 1` for `m ≥ 2`; `m = 1` (`B` is `τ`-stable) is recorded
and excluded. Cross-check on every structure the tool evaluates in full: `g`
equals the sum of `(s − 2)/2` over all nonzero vectors of the span
(`s` = support size).

Only odd prime `n` (3 and 5). **`n = 4` is left out.** For odd `n` the
Frobenius extends to the compositum with order `n` (Diem, "The GHS attack in
odd characteristic", 2003; as in the merged census, cited, not re-verified).
For even `n` the extension of `τ` to `K(x)(√f_0, …, √f_{n−1})` sends
`√f_{n−1} ↦ c·√f_0` and the composite `τ^n` acts on `√f_0` by a scalar whose
being 1 depends on the square class of a norm of the leading coefficient of
`f` — data the census does not carry (it works with branch sets only, which
are blind to quadratic twists). Deciding descent for `n = 4` therefore needs
the forms, not the sets; it is not attempted here.

## Enumeration (new) and why it is the same set of structures

By Lang's theorem every structure is `τ = MσM^{−1}`, `M ∈ PGL2(K)`, unique up
to `M ↦ MN`, `N ∈ PGL2(F_q)`; `A = M(M^σ)^{−1}`. Then
`τ^i B = M σ^i (M^{−1}B)`, so `(m, r, g)` and every support size of
`(B, τ)` equal those of `(B' = M^{−1}B, σ)`. The tool enumerates
`N = M^{−1}` over a transversal of `PGL2(F_q)\PGL2(K)` given by normal forms
of the ordered triple `(N(∞), N(0), N(1))` (for `n` an odd prime,
`PGL2(F_q)` acts freely on `P¹(K) ∖ P¹(F_q)`):

- `N(∞) = ∞`, `N(0) = 0`, `N(1)` a representative of `K*/F_q*`;
- `N(∞) = ∞`, `N(0)` a representative of `(K∖F_q)/Aff(F_q)`, `N(1)` any other finite point;
- `N(∞)` a representative of `(K∖F_q)/PGL2(F_q)`, `N(0), N(1)` any distinct other points.

The tool asserts the total equals `(Q³ − Q)/(q³ − q)`, `Q = q^n`, and that
every curve's census counted exactly that many structures. Then
`B'_i = σ^i B'` is the coefficient-wise Frobenius of the quartic form.
Fast path: the `σ^i B'` are pairwise disjoint iff `gcd(B', σ^k B') = 1` for
`1 ≤ k ≤ (n−1)/2` (since `σ^i B' ∩ σ^j B' = σ^i(B' ∩ σ^{j−i}B')` and
`σ^n B' = B'`), and then `m = n`, `r = 4n`. Every other structure is
evaluated in full (span by subset sums, union, subcover-sum check).

The merged census's brute cocycle search (cost `Q³` per field) is kept,
generalised to `n`, as a second method (`cocycles`).

## Reproduction checks (run first; abort and fix on any difference)

- **R1.** `n = 3`, `q = 5, 7`, exhaustive, `cosets`: for every curve the
  fields `A_idx, B_idx, order, maxprime, min_genus, types_at_min,
  structures_by_genus, m1_structures, min_genus_2isog, min_genus_class`
  equal the merged run `20261013T0403Z` record for the same curve, and the
  merged `mod4` histograms (`curves_by_min_genus` per scope and class) are
  reproduced by the new `group` records with `by = order_mod4`.
- **R2.** `n = 3`, `q = 11`, `sample:120`, `cosets` (same seed and stream
  as the merged tool): the same fields (minus `min_genus_class`, not
  computed) equal the merged `q = 11` sample records.
- **R3.** `n = 5`, `q = 3`, exhaustive: `cosets` and `cocycles` give
  identical curve records (all fields except `mismatch`, whose scope differs
  between the methods).

## Panels (frozen)

Run in this order, one process at a time, 6 threads, after R1–R3 pass.
Fields as in the merged tool (first primitive monic modulus). Curves: every
`K`-isomorphism class (short Weierstrass for `p > 3`, same representatives
as the merged tool; general `y² = x³ + a₂x² + a₄x + a₆` for `p = 3`).
Samples: the merged tool's fixed-seed stream (`seed = 20261012`, tag `q` for
`n = 3`, tag `100q + n` for `n = 5`), each sampled curve with its
`K`-rational 2-isogeny neighbours; every structure of every model.

| panel | `n` | `q` | `K` | curves | structures per curve | coverage |
| --- | --- | --- | --- | --- | --- | --- |
| A | 3 | 11 | `F_{11³}` | all classes | 1 786 323 | exhaustive, class minima |
| B | 3 | 13 | `F_{13³}` | all classes | 4 855 539 | exhaustive, class minima |
| C | 3 | 17 | `F_{17³}` | sample of 200 + neighbours | 24 221 379 | every structure of each sampled model; curves sampled; no class minima |
| D | 5 | 3 | `F_{3⁵}` | all classes | 597 861 | exhaustive, class minima (also R3) |
| E | 5 | 5 | `F_{5⁵}` | sample of 60 + neighbours | 254 313 125 | every structure of each sampled model; curves sampled; no class minima |

Exhaustive vs sampled: in every panel each model is censused over **all**
of its structures, so its minimal genus is exact. Panels A, B, D cover every
curve; C and E are samples of curves, so their fractions are sample
estimates (reported with the count `n`, no interval claimed beyond the
binomial count). Stop rule: if a panel does not finish within 6 hours wall
time it is stopped and reported as incomplete; nothing is rerun with a
different sample after seeing results.

## Derived statements (proved in this file before the runs; a violation means a bug or an error in the proof)

Notation: `U = ∪_i σ^i B'` is `σ`-stable (because `σ^n B' = B'`), so it is
a union of `σ`-orbits; `r = |U|` is the sum of the sizes of the orbits that
meet `B'`. A point of `B'` of degree `d` over `K` (`d ∈ {1,2,3}`) has
`σ`-orbit size `e` with `e | nd`, `e ∤ n·d'` for `d' < d`, so for `n` prime
`e ∈ {d, nd}`. If `e` is prime to `n`, `σ^n` generates the same group as
`σ` on that orbit, so `B'` contains the whole orbit; if `B'` met only such
orbits it would be `σ`-stable (`m = 1`).

**Lemma M (the span is a cyclic module).** The span of `{σ^i v}`,
`v = 1_{B'}`, is the `F_2[x]`-module generated by `v` with `x` acting as
`σ`, and `(x^n − 1)v = 0`. So its dimension is the degree of a monic divisor
of `x^n − 1` over `F_2`. `x³ − 1 = (x+1)(x²+x+1)` and
`x⁵ − 1 = (x+1)(x⁴+x³+x²+x+1)`, both factors irreducible (2 has order 2 mod
3 and 4 mod 5). Hence `m ∈ {1, 2, 3}` for `n = 3` and **`m ∈ {1, 4, 5}` for
`n = 5`**, with `m = n − 1` iff `Σ_i σ^i v = 0`, i.e. iff every point of `U`
lies in an even number of the `σ^i B'`.

**Lemma D0.** For `n` prime and `m ≥ 2` the `n` sets `σ^i B'` are distinct,
each gives an elliptic quotient, so `g ≥ n`.

**Lemma E (the even configurations).** If `m = n − 1`, each orbit meeting
`B'` has size `nd` and meets `B'` in `|S|·d` points with `|S|` even
(`S ⊂ Z/n` the cosets of `⟨n⟩` it contains). As `B` has a `K`-rational
point, two conjugate pairs (`d = 2`) cannot fill `B'`, and `d = 3` needs 6
points. So `m = n − 1` forces **all four points of `B'` to be
`K`-rational** (`E(K)[2] = (Z/2)²`), in `n`-orbits meeting `B'` in 2 or 4
points: for `n = 5`, `r ∈ {5, 10}` and `g ∈ {5, 25}`; for `n = 3`,
`r = 6`, `g = 3`.

**Corollary (`n = 5` spectrum).** `m = 5` and `r = 5` is impossible (an
`r = 5` set `U` is a single 5-orbit containing `B'`, whose points then lie
in 4 of the 5 sets). So every `n = 5` genus is 5, 25 (`m = 4`) or
`8(r − 4) + 1 ≥ 17` (`m = 5`, `r ≥ 6`): **no genus in 6–16, genus 5 only
for `E(K)[2] = (Z/2)²`**, attained exactly by `B' = O(x) ∖ {σ^j x}` with
`x ∈ K ∖ F_q`.

**Proposition O5 (odd order, `n = 5`).** If `#E(K)` is odd (cubic
irreducible over `K`) then every structure has `m = 1` or
`g ∈ {33, 97, 129}`, and exactly one structure has `m = 1` or `g = 33`.
*Proof.* A cubic root `θ'` of `B'` has orbit size `e ∈ {3, 15}`, the
rational point `a'` has orbit size 1 or 5. `(e, |O(a')|) = (3,1)`: `m = 1`.
`(3,5)`: `r = 8`; `(15,1)`: `r = 16`; `(15,5)`: `r = 20`; in all three
`Σσ^i v ≠ 0` (the orbit of `a'` resp. the 15-orbit is covered an odd number
of times), so `m = 5` and `g = 33, 97, 129`. `PGL2(K)` acts simply
transitively on points of degree 3 over `K` (orbit counting:
`|PGL2(K)| = Q³ − Q`, the stabiliser fixes 3 points), and `PGL2(F_q)` acts
simply transitively on `F_{q³} ∖ F_q`, which lies in that set because
`gcd(3, 5) = 1`; so the `N` with `N(θ) ∈ F_{q³}` form exactly one coset,
i.e. one structure. ∎ Moreover `m = 1` occurs iff `B` is `PGL2(K)`-equivalent
to a `σ`-stable set, which forces `j(E) ∈ F_q`; for `j ∉ {0, 1728}`,
`j ∈ F_q` also suffices (the base change of a curve over `F_q` with that `j`
has a branch set of the same type and invariant). So for odd order and
`j ∉ {0,1728}`: **minimal genus 33 iff `j(E) ∉ F_q`, otherwise ≥ 97.**
(Odd-order analogue of the merged `OrbitNine`; there `3 | n` forces `e = 9`.)

**Proposition T (the `(1,1,2)` obstruction, any odd `n`).** For a set of
type (two `K`-points `a₁, a₂`, one conjugate pair `β, β̄`) let
`μ = cr(a₁, a₂; β, β̄)`. Frobenius of `K` inverts `μ`, so `μ` lies in the
norm-one group `U₁ ⊂ F_{Q²}*` (cyclic of order `Q + 1`, even), and its class
`χ ∈ U₁/U₁² ≅ Z/2` is a `PGL2(K)`-invariant independent of the orderings.
For `B = {∞, e, β, β̄}`: `μ = (β−e)^{Q−1}`, and `χ = 1` iff `β − e` is a
square in `F_{Q²}` iff `g(e) = N(β − e) ∈ K*²` (`f = (x−e)g`) iff
`T = (e,0) ∈ 2E(K)` iff `4 | #E(K)` (the 2-part is cyclic). For the
configuration `B' = {β, β^q} ∪ {x, σ^k x}`, `β ∈ F_{q²} ∖ F_q`,
`x ∈ K ∖ F_q`: with `h(t) = (β − t)/(β^q − t)`, `h(x) ∈ U₁` and
`h(σt) = h(t)^{−q}`, so `μ = h(x)^{1 − (−q)^k}` with an even exponent:
`χ = 1`. Hence **that configuration forces `4 | #E(K)`.**

**Corollary T3 (`n = 3`).** Genus 3 at `n = 3` is `(m,r) = (2,6)` (Lemma E:
split) or `(3,5)`: `U` = a 3-orbit plus two `F_q`-points (all of `B'`
rational: split) or a 3-orbit plus a 2-orbit (`B' = {β,β^q} ∪ {x, σx}`:
Proposition T). Every case has `4 | #E(K)`. So **`#E(K) ≡ 2 mod 4` ⇒ every
structure has `g ≥ 5`** at `n = 3` (the merged census observed this; this is
its proof).

**Corollary T5 (`n = 5`).** For `E(K)[2] = Z/2`, `m = 5` (Lemma E) and `U`
must contain a 5- or 10-orbit and the orbit of the pair, so `r ≥ 7`,
**`g ≥ 25`**; `r = 7` only for `B' = {β,β^q} ∪ {x, σ^k x}`, so
**`#E(K) ≡ 2 mod 4` ⇒ `g ≥ 33`.**

**Lemma H (hyperellipticity flag).** For `m ≥ 3` and `g ≥ 2` the cover is
hyperelliptic iff some involution of the cover group has a genus-0 quotient
(the hyperelliptic involution is central, and `(Z/2)³` does not embed in
`PGL2` in odd characteristic), i.e. iff a hyperplane of the span has only
support-2 vectors. A space of support-2 vectors has dimension ≤ 2, so for
`m ≥ 4` the cover is never hyperelliptic. `m = 2` stays undetermined.

## Predictions and falsifiers

- **n = 3, derived (any violation is a bug):** odd ⇒ every structure has
  `g ∈ {13, 17}` and `m ≠ 1`; `2 mod 4` ⇒ `g ≥ 5` on every structure.
- **n = 3, measured prediction P3:** at `q = 11, 13` (exhaustive) the class
  minimum is 13 for every odd curve, 3 for every `0 mod 4` curve and 5 for
  every `2 mod 4` curve. Falsifier: any curve with another class minimum.
  (The merged census already had two `2 mod 4` curves at `q = 5` with class
  minimum 9, so the `2 mod 4` row is known to have exceptions at small `q`;
  an exception here is reported as a falsification of P3, not explained
  away.) At `q = 17` (sample) class minima are not computed; the
  `≤ 1 step` minimum is reported without a prediction.
- **n = 5, derived (any violation is a bug):** D0 `g ≥ 5`; D1
  `m ∈ {1,4,5}`; D2 `m = 4` only for `E(K)[2] = (Z/2)²` with `r ∈ {5,10}`,
  no genus in 6–16, genus 5 only for `E(K)[2] = (Z/2)²`; D4 odd order:
  genera in `{33, 97, 129}`, `m1 + #(g = 33) = 1`, and for `j ∉ {0, 1728}` `m1 = 1` iff `j ∈ F_q`; D5 `E(K)[2] = Z/2` ⇒
  `g ≥ 25`, and `2 mod 4` ⇒ `g ≥ 33`. The tool checks each on every curve
  (`violations` field, `derived_checks` record).
- **n = 5, measured:** no prediction beyond the derived bounds. In
  particular no prediction for how often genus 5 (split curves), 25 (type
  `(1,1,2)` with `4 | #E(K)`) or 33 is attained, or for any dependence on
  `E(K)[4]` or `v₂(#E(K))`. Note: the genus data depend on `B` only, which
  is shared by `E` and its quadratic twist; `#E^t(K) = 2Q + 2 − #E(K)` is
  congruent to `#E(K)` mod 4 but not mod 8, so no law in `v₂ ≥ 3` or
  `E(K)[4]` that separates twists can hold.

## Measures

Per curve: `#E(K)`, `v₂`, `#E(K)[2]`, `#E(K)[4]`, `j ∈ F_q`, minimal genus,
types `(m, r, hyperelliptic)` at the minimum with structure counts, the
histogram of all structures by genus and by `(m, r)`, `m = 1` count,
minimum over the curve and its `K`-rational 2-isogeny neighbours, and (for
exhaustive panels) over the isogeny class (same `#E(K)`). Per field: the
same histograms grouped by `#E(K) mod 4`, `#E(K)[2]`,
`(#E(K)[2], #E(K)[4])`, `v₂`, `(mod 4, #E(K)[2])`, and (odd order) `j ∈ F_q`;
`(m, r, #E(K)[2])` structure totals; derived-check violation counts.

## Not measured

The `F_q`-model of `C`, any transfer map or kernel, any index calculus, `n = 4`
(reason above), `n ≥ 7`, characteristic 2, degree-2 maps other than `x`
(a `K`-rational degree-2 class `O + S` with `S ∉ 2E(K)` gives a branch set
of a different Galois type with the same `j`; not counted, as in the merged
census), isogenies other than 2 as explicit steps.

## Language

Rust (`cover_census_n.rs`, dependency-free, `rustc --edition 2021 -O`).
Comparisons with `jq` and `cmp`. Python only renders recorded values.
Lean 4 (no Mathlib) for the feasible lemmas.

## Runs before freezing

- Development pilot `n = 3`, `q = 3` (not a panel; char 3 was outside the
  merged census): both methods, identical curve records; seen before this
  file was written. It showed two `0 mod 4` curves (`j = 0`, supersingular in
  characteristic 3) whose isogeny-class minimum is 7, so the `0 mod 4 → 3`
  pattern is already known to fail at `q = 3`. Kept as a labelled pilot.
- Timing pilot `n = 5`, `q = 5`, `sample:1`, output sent to `/dev/null`
  unread. The session running it ended before it finished (about 8 minutes wall, 3.4 CPU-minutes, at load average about 90); no timing line survived and the scratch output was lost. Panel sizes above were set from that partial observation (roughly 0.3–0.4 µs of CPU per structure) and the stop rule. The `n = 3`, `q = 3` pilot records were lost with the same scratch directory; the pilot is rerun once, unchanged, and committed as `results/pilot-n3-q3.jsonl`, labelled pilot. (As run, the rerun pilot landed as
  `results/run-20261014T0050Z/pilot-n3-q3-{cosets,cocycles}.jsonl`.)

## Amendment 2026-10-14 (after panel D, before panels A, B, C, E produced output)

Adds one hypothesis; changes no panel, input or derived statement.
Panel D (`n = 5`, `q = 3`, exhaustive) showed: among `E(K)[2] = Z/2`,
`4 | #E(K)` curves, genus 25 on the curve itself occurs exactly for
`v₂(#E(K)) ≥ 3` (60 of 60), never for `v₂ = 2` (61 curves, minimum 33 or
65); genus 5 occurs only with `#E(K)[4] = 8`. A candidate reason, not a
proof: in Proposition T the configuration has `μ = h^{1 − (−q)^k}`; at
`q ≡ 3 mod 4` both exponents (`1 + q`, `1 − q²`) are divisible by 4 and
`K* ⊂ F_{Q²}*⁴`, so `μ` is a fourth power in `U₁`, a finer 2-adic condition
than Proposition T; at `q ≡ 1 mod 4` the `k = 1` exponent `1 + q` is
`2 mod 4` and no such refinement follows.

**Hypothesis H (prediction for panel E, `q = 5 ≡ 1 mod 4`):** the `v₂ ≥ 3`
restriction does not persist: some sampled `E(K)[2] = Z/2` curve with
`v₂(#E(K)) = 2` has genus 25 on itself. Falsifier: in panel E every
genus-25 `E(K)[2] = Z/2` curve has `v₂ ≥ 3`. With fewer than 5 sampled
`E(K)[2] = Z/2`, `v₂ = 2` models the test is reported as underpowered, not
as a pass.
