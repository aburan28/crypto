# Protocol: at which prime degree does each representation of an ordinary `ℓ`-isogeny become the cheapest, and does the class-group route reproduce the walker's certified edges

Frozen 2026-10-08, before any instrument is run at the registered sizes.
**Cost study, toy and asymptotic.**  `S`, end-to-end cost and ECDLP
speedup are **unset** and are not touched by any outcome below.  The
subject is the cost of one certified degree-`ℓ` edge in the `F_p`-isogeny
class of an ordinary curve, for prime `ℓ` beyond the walker's cap of 61.
The candidates are C1 to C7 of [`README.md`](README.md).

## What is already known, and the gap

- The walker (`src/cryptanalysis/isogeny_walk/`) builds `Φ_ℓ mod p` from
  `q`-series at about `ℓ⁵` multiplications and is capped at `ℓ = 61`
  (`research/p256_isogeny_multigraph_64k_20261006/`).  Its 302,653
  certified P-256 edges, replayed by an independent verifier, are the
  correctness oracle for this protocol.
- For registered P-256: `t = 89188191154553853111372247798585809583`,
  `D_π = t² − 4p` has 258 bits, `D_π ≡ 1 (mod 4)`, and `−D_π = 3 · 5 · q₁
  · q₂ · q₃` with `q₁ = 456597257999`, `q₂ = 1428624589419343516204097`,
  `q₃ = 46523541035814968339936406074986559003387` distinct primes (PARI
  `factorint`, 2026-10-08, before freezing).  `D_π` is fundamental, `Z[π]`
  is maximal and `End(E) = Z[π]` across the class.  The walker's
  `ClassInfo` records the same `D_π` with the cofactor `q₁q₂q₃` as
  `composite_unfactored`, since its trial bound is below `q₁`.
  `(D_π / ℓ) = +1` for `ℓ ∈ {11, 13, 17, 23, 29, 37, 41, 43, 47, 59}` and
  `−1` for `{7, 19, 31, 53}`, matching the degrees the walker found
  walkable.
- The literature costs in the README table are asymptotic.  None has
  been measured against this walker, and no crossover degree between
  the `Φ_ℓ` route, the rational-torsion route and the class-group route
  has been measured for any registered curve.

The gap is a measured crossover table, and a check that the class-group
route lands on the curves the walker certified.

## Shared definitions

- `λ, μ` are the Frobenius eigenvalues mod `ℓ`, roots of `x² − tx + p`
  in `F_ℓ`; they exist exactly when `(D_π / ℓ) = +1`.  `r(λ)` is the
  multiplicative order of `λ` in `F_ℓ^*`.  The `λ`-eigenspace of `E[ℓ]`
  is a cyclic subgroup of order `ℓ`, stable under Frobenius, with every
  point defined over `F_{p^{r(λ)}}` and none over a smaller field.
- `r_min(ℓ) = min(r(λ), r(μ))`.  THEOREM: `r_min(ℓ)` is the least `k ≥ 1`
  with `ℓ | #E(F_{p^k})`, where `#E(F_{p^k}) = p^k + 1 − s_k`, `s_0 = 2`,
  `s_1 = t`, `s_k = t s_{k−1} − p s_{k−2}`.  The E1 instrument self-tests
  on this identity.
- The eigenvalue model: `λ` behaves as a uniform element of the cyclic
  group `F_ℓ^*`, so `Pr[r(λ) ≤ R] = S_R(ℓ) / (ℓ − 1)` with
  `S_R(ℓ) = Σ_{d | ℓ−1, d ≤ R} φ(d)`, and
  `Pr[r_min ≤ R] ≈ 2 S_R(ℓ) / (ℓ − 1)` for `R` small against `ℓ`.
  HEURISTIC; E1 tests it.
- Cost unit: `F_p` multiplications.  A multiplication in `F_{p^r}` is
  charged `r²`.  A scalar multiplication by a `b`-bit integer over
  `F_{p^r}` is charged `12 b r²`.  A √élu evaluation of a degree-`ℓ`
  isogeny over `F_{p^r}` is charged `c_v · √ℓ · log₂ ℓ · r²` with
  `c_v = 60` (UNTESTED constant; E1 reports the crossover for `c_v ∈
  {20, 60, 200}`).  One `Φ_ℓ`-route edge at degree `ℓ` with `Φ_ℓ` already
  built is charged `2 ℓ²` (evaluation) `+ 4 ℓ log₂ p` (root) `+ 20 ℓ²`
  (kernel) `+ 20 ℓ²` (certificate); the `ℓ⁵` construction is amortised
  over the run's node count and reported separately.

## E1: eigenvalue orders and the rational-torsion window (candidate C4)

Instrument: [`eigenvalue_orders.py`](eigenvalue_orders.py), pure Python,
no dependencies.  Self-test passed 2026-10-08 on 4,286 `(p, t, ℓ)`
triples; the output path was smoke-tested at `ℓ ≤ 100` only, and no
registered size has been run.  Inputs: `(p, n)` for P-192, P-224, P-256 as in the
walker's `StartCurve`.  Sizes: every prime `ℓ ≤ 2²⁰` (82,025 primes).
Deterministic; no seed.

For every `ℓ`: the kind (`divides`, `atkin`, `elkies`), `r(λ)`, `r(μ)`,
`r_min`, the model probability `2 S_R(ℓ)/(ℓ − 1)` for `R ∈ {2, 3, 4, 6, 8,
12, 24, 48}`, and the charged cost of the C4 route against the `Φ_ℓ`
route.  Outputs: `results/e1_{curve}.jsonl` (one line per Elkies `ℓ`) and
`results/e1_summary.json`.

Self-test, run before the registered sizes and recorded in the summary:
for 200 random `(p, t)` with `p` a prime below `2²⁰` and `|t| ≤ 2√p`,
and every Elkies `ℓ < 200`, the computed `r_min` equals the least `k ≤
500` with `ℓ | p^k + 1 − s_k`.  One mismatch fails the instrument.

Predictions:

- **E1.1 (model).**  For each curve and each `R ∈ {4, 8, 12, 24}`, the
  count of Elkies `ℓ ≤ 2²⁰` with `r_min ≤ R` lies within a factor 2 of
  `Σ_ℓ 2 S_R(ℓ)/(ℓ − 1)` summed over the same `ℓ`.  *Falsified by one
  `(curve, R)` outside `[0.5, 2]`.*
- **E1.2 (twist consistency).**  `r_min = 1` occurs only at `ℓ = n`
  (never, for `ℓ ≤ 2²⁰`, since `n` is a 256-bit prime), and the set of
  `ℓ` with `r_min = 2` equals the set of primes `ℓ ≤ 2²⁰` dividing the
  quadratic-twist order `2p + 2 − n`.  *Falsified by any difference.*
  This ties E1 to the walker's `twist_factors`.
- **E1.3 (window).**  For P-256 at least 20 Elkies `ℓ` in `[61, 2²⁰]`
  have `r_min ≤ 12`, and at least one has `r_min ≤ 4` and `ℓ > 10³`.
  *Falsified if either count is short.*
- **E1.4 (crossover, charged).**  With `c_v = 60`, the C4 route is charged
  below the `Φ_ℓ` route for every `ℓ > 200` with `r_min ≤ 6`, and above it
  for every `ℓ < 100` with `r_min ≥ 8`.  *Falsified by one `ℓ` on the
  wrong side.*  This is a statement about the registered charge sheet,
  not a timing.

Decision rule: E1.3 passes and E1.4 passes: C4 is admitted as an edge
kind for the listed `ℓ`, and the walker extension in E2b carries an
`F_{p^r}` kernel-point witness.  E1.3 fails: C4 stays a curiosity and the
large-`ℓ` question is decided between C1 and C5.  E1.1 fails: the
eigenvalue model is wrong for these curves and the failing direction is
reported before any use of the model elsewhere.

## E2: the class-group route (candidate C5)

### E2 step 0: the order

Settled before freezing: `−D_π = 3 · 5 · q₁ q₂ q₃` with three distinct
primes, so `D_π` is fundamental, `Z[π] = O_K` is maximal, and `End(E) =
Z[π]` on the whole class.  E2a therefore computes in `Cl(O_K)` with
PARI's `bnfinit`, and the relation it finds is valid for every curve the
walker can reach.  The route is invalid only if `ℓ | f_π`, and `f_π = 1`.

Kept for other classes, THEOREM: for `q ∤ f`, a relation `∏ 𝔮_i^{e_i} =
(α)` in `Z[π]` extends to the same relation in any intermediate order
`O`, so a relation valid in `Cl(Z[π])` yields the correct codomain under
the `Cl(End(E))` action whatever `End(E)` is.  The instrument refuses a
non-fundamental `D` rather than silently using the maximal order.

### E2a: class group and relation norms (PARI)

Instrument: [`relation_lattice.gp`](relation_lattice.gp), PARI/GP 2.17.
Self-tested 2026-10-08 on a toy class (`p = 1000003`, 22-bit `D`, `h =
344`): every target reduced to `|e|₁ ≤ 7` and every conjugate check
passed.  Not run at the registered size.  Steps:

1. `quadclassunit(D_π)` under GRH; record `h`, the cyclic structure, and
   the generators.  `bnfcertify` is not attempted at 258 bits.
2. Factor base `B`: the Elkies primes `q ≤ 200`.  For each, the
   coordinates of `𝔮` in the cyclic decomposition (discrete logarithm
   against the generators).
3. Relation lattice `L ⊂ Z^{|B|}`, LLL-reduced.
4. For each target `ℓ` in `T = {11, 13, 17, 23, 29, 37, 41, 43, 47, 59}
   ∪ {the first Elkies prime above each of 10³, 10⁴, 10⁵, 10⁶, 2³², 2⁶⁴}`:
   one solution `e` of `∏ 𝔮_i^{e_i} ~ 𝔩`, Babai-reduced against `L`; its
   `ℓ₁` norm `|e|₁` and its charge `Σ |e_i| · (2 q_i² + 4 q_i log₂ p + 40
   q_i²)`.  A target that is itself in `B` is solved with its own column
   removed, so the relation is never the trivial one-step vector; the
   record carries an `in_base` flag.
5. Also for each `ℓ` in `T` the second solution for `𝔩̄`, which must be
   `−e` modulo `L` (a consistency check, not a prediction).

Outputs: `results/e2a_classgroup.json`, `results/e2a_relations.jsonl`.
Wall-time cap 6 hours for step 1; a step 1 that does not finish is
reported as partial and E2b is not run.

Predictions:

- **E2.1 (size).**  `log₂ h ∈ [120, 136]`.  *Falsified outside.*  (The
  class number of a discriminant of 258 bits is `√|D| · L(1, χ) / π`
  with `L(1, χ)` within a factor 20 of 1 under GRH.)
- **E2.2 (norms).**  For every `ℓ ∈ T`, `|e|₁ ≤ 400`.  *Falsified by one
  target above 400.*  The expectation under a uniform-lattice model with
  `|B| ≈ 20` and `h ≈ 2¹²⁸` is `|e|₁` of a few hundred.
- **E2.3 (independence of ℓ).**  The maximum of `|e|₁` over the six large
  targets is at most 1.5 times the maximum over the ten small targets.
  *Falsified otherwise.*
- **E2.4 (charged crossover).**  The charged cost of the class-group
  route is below the charged cost of the `Φ_ℓ` route, construction
  included at 65,536 nodes, for every `ℓ ∈ T` with `ℓ > 10³`, and above it
  for every `ℓ ≤ 61`.  *Falsified by one `ℓ` on the wrong side.*

### E2b: route verification against certified edges

Requires a walker extension not in this PR: an **eigenvalue tag** on each
certified edge, computed inside `verify_kernel` as the `λ ∈ F_ℓ` with
`x([λ]P) ≡ x^p (mod h)` on the kernel polynomial `h`, which orients the
edge as `𝔮` or `𝔮̄`.  The tag is part of the certificate and is replayed.

With the tag: for each `ℓ ∈ {11, …, 59}` of `T`, start at registered
P-256, walk the relation vector `e` step by step choosing at each step the
edge whose tag matches the ideal's eigenvalue, and compare the final
`j` with the `j′` of the walker's certified degree-`ℓ` edge of the same
orientation.  Then repeat from 32 random nodes of the 65,536-node run.

- **E2.5 (reproduction).**  All `10 × 33` routes end at the certified
  `ℓ`-neighbour, up to `F_p`-isomorphism.  *Falsified by one route that
  does not.*  One failure means either a wrong orientation convention,
  which the tag data will show as a global sign, or an `End(E)` smaller
  than assumed, which step 0 will show.
- **E2.6 (beyond the cap).**  For the six large `ℓ`, the routes for `e`
  and for `−e` from the root end at two distinct curves, and a second,
  independently reduced relation `e′ ≠ e` for the same `𝔩` ends at the
  same curve as `e`.  *Falsified by a collision or a disagreement.*  This
  is the only check available where no `Φ_ℓ` edge exists.

Decision rule: E2.2, E2.3 and E2.5 pass: C5 is admitted; a `class-group`
edge kind with the relation vector as witness is specified for the
walker, and the `2³²`-population goal of the 2026-10-06 run is restated
as navigation in `Cl(Z[π])`.  E2.2 fails: the route is reserved for
`ℓ` above the measured crossover and C1 is the walk's main line.  E2.5
fails with step 0 showing `C` non-squarefree: the experiment is rerun in
the class group of the order step 0 identifies as `End(E)`, determined
by the volcano depth at each prime of `f`, and the failure is not read
as a failure of C5.

## E3: BMSS kernel recovery (candidate C3)

Design only; needs Rust work.  Implement the [BMSS08] solve behind a
flag in `kernel.rs`, same inputs `(a, b, a′, b′, ℓ, Σ x(Q))`.  Replay all
302,653 certified P-256 edges with both derivations.

- **E3.1 (agreement).**  Identical kernel polynomial on every edge.
  *Falsified by one difference.*
- **E3.2 (scaling).**  Fitted exponent of kernel-derivation time in `ℓ`
  over `ℓ ∈ {11, …, 59}`: current `≥ 1.7`, BMSS `≤ 1.3`.

## E4: multipoint evaluation across the frontier (candidate C2)

Design only.  Evaluate the `ℓ + 1` column polynomials of `Φ_ℓ` at the
frontier in blocks of `ℓ` points by a product tree.

- **E4.1 (crossover).**  Measured over `ℓ ∈ {11, …, 59, 127, 251}` with
  `Φ_ℓ` precomputed, the degree at which the block tree beats Horner is
  at most 200.  *Falsified if Horner still wins at 251.*
- **E4.2 (exactness).**  Identical roots, hence identical edges, on the
  full replay.

## E5: construction scaling (candidate C1)

Design only.  Time `modular_polynomial(f, ℓ)` for `ℓ ∈ {11, 13, 17, 23,
29, 37, 41, 43, 47, 59, 127, 251}` as is, then with Kronecker-substitution
products, then with the [BLS12] CRT construction.

- **E5.1 (exponents).**  Fitted exponents: current `≥ 4.5`; Kronecker
  `≤ 4.2`; BLS `≤ 3.5`.  *Falsified by any fitted exponent on the wrong
  side.*
- **E5.2 (agreement).**  All three constructions give the same
  coefficient table mod `p` at every `ℓ`.
- **E5.3 (budget).**  With BLS, `Φ_ℓ mod p` for `ℓ = 1,009` is built in
  under 1 hour on one core.  *Falsified otherwise.*

## Combined reading

One table, one row per `ℓ` in `T ∪ {E1's window primes}`: charged cost
of C1-route, C4-route, C5-route; the cheapest; and whether the cheapest
is also a verified edge (E2.5, E3.1, E5.2).  The claim allowed if the
predictions hold is: *for registered P-256, the representation that
produces a certified `ℓ`-edge at least charged cost changes from the
`Φ_ℓ` route to the class-group route between `ℓ = 61` and `ℓ = 10³`,
with the rational-torsion route cheapest on the sparse set E1.3
identifies.*  No claim about any other class, any other field, or any
DLP cost is allowed.

## Stop conditions and inadmissible moves

Wall-time caps: E1 30 minutes per curve; E2a step 1 6 hours; every
E2b route 10 minutes.  A size that does not finish is reported as
partial, never extrapolated.

Inadmissible: changing `c_v` or the charge sheet after seeing a
crossover; adding or removing targets from `T`; widening the E2.2 bound;
reading a passed E2.5 as evidence about `End(E)` beyond what step 0
established; reading any row of the combined table as a statement about
`S` or about ECDLP.

[BLS12], [Sut13], [BMSS08], [BDLS20], [BKV19]: see README.
