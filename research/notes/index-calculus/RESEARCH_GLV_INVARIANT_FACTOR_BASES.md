# Endomorphism-invariant factor bases across curve families: plan, pilot, experiments E1–E9, and the road to the state of the art

**Modules:** `src/cryptanalysis/glv_invariant_base.rs` (the fold, prime-field automorphisms, Vélu degree-2 and degree-3 endomorphisms, CM instance generators, the folded rho classes), `src/cryptanalysis/ext_curve.rs` (`ExtField` over `F_{p²}` and `F_{p³}`, the generic `ExtCurve` counted group, diagonal automorphisms and Frobenius-type maps on it), `src/cryptanalysis/gls_fp2.rs` (GLS `ψ`, the `ψ`-stable line, the `j = 0` and `j = 1728` twists with their lifted automorphisms), `src/cryptanalysis/subfield_fp3.rs` (`E/F_p` on `E(F_{p³})`, the Frobenius eigenline), `src/cryptanalysis/line_oracle.rs` (the Weil-descent resultant oracle for a line, E2b), `src/cryptanalysis/glv_invariant_experiments.rs` (one relation stream feeding both arms to full rank), `src/cryptanalysis/ic_framework/plugins.rs` (`glv-orbit`, `gls-line`), `src/cryptanalysis/ic_boundary.rs` (`FactorBase::from_column_map`)
**CLI:** `ic bench --bits 20 --family j0 --factor-base glv-orbit:size=64 --oracle mitm:negation_folded=1` (control: `glv-orbit:size=64,no_fold=1`)
**Bench (pilot):** `cargo run --release --example glv_invariant_bench -- --families j0,j1728,generic,d7,d8 --bits 16,20,24 --seeds 2 --oracles subtract,mitm --json experiments/23_glv_invariant_pilot.json`; `--families gls --bits 8,10,12 --oracles subtract --json experiments/23_glv_invariant_gls_pilot.json`
**Runner (E1–E7):** `cargo run --release --example glv_invariant_experiments -- --exp e1 --bits 16,20,24,28 --seeds 6 --json experiments/23_glv_invariant_e1.json` (the exact command of every file is its `command` field)
**Data:** `experiments/23_glv_invariant_pilot.{json,log}`, `experiments/23_glv_invariant_gls_pilot.{json,log}` (2026-09-28), `experiments/23_glv_invariant_e{1,1_32,2,3,4,5,6,7}.{json,log}` (2026-09-29), `experiments/23_glv_invariant_e{8,9,11,11_12,11_13,11_14}.{json,log}` (2026-10-01); this host: Linux x86-64, 4 threads; wall time is recorded and is not a result
**Tables:** `python3 scripts/glv_invariant_tables.py experiments/23_glv_invariant_pilot.json experiments/23_glv_invariant_gls_pilot.json` (§5) and `python3 scripts/glv_invariant_experiment_tables.py experiments/23_glv_invariant_e*.json` (§6, §8); every number in §5, §6 and §8 is printed by them from the frozen files

> **Status.**  Implementation, pilot, and the seven experiments of §4
> run and read (§6).  The fold is one function over the framework's
> `CountedGroup`, verified end to end on `F_p`, `F_{p²}` and `F_{p³}`
> groups; the pilot (§5) is at `2^13`–`2^24` with two seeds per size,
> the experiments (§6) at `2^13`–`2^32` with four to six seeds, both
> arms fed one relation stream to full rank.  What they establish: the
> fold is exact (every row recovers its planted logarithm); the column
> count is the orbit count; after `columns` relations both arms reach
> the same rank fraction, so the relation count moves with its floor
> and nothing else moves; the fold by **any** finite group of
> endomorphisms is the order of the subgroup of `(Z/rZ)^*` its
> eigenvalues generate (`12` on the `j = 0` GLS twist, `4` — not `8` —
> on the `j = 1728` twist, `3` — not `9` — on the `j = 0` subfield
> curve); the degree-2 and degree-3 CM endomorphisms keep **zero** base
> points in the base at eigenvalue orders of `10²`–`10⁸`; and the rho
> folded by the same group takes `√(w/2)` fewer steps, verified.
> Class, where anything moved: **engineering** (§6.8); nothing moved a
> ratio to a floor, and no scoreboard row claims a speed.  §8 places
> this against the literature (Galbraith–Granger–Merz–Petit's invariant
> bases, Faugère–Gaudry–Huot–Renault's system symmetries), ranks what
> is next, and runs the first three: three summands on the Frobenius
> line (E8: the fold's relation count falls with the columns once the
> two-summand degeneracy of §6.5 is gone, and every phase is priced in
> one unit beside the matched rho), the folded rho on the `F_{p²}` and
> `F_{p³}` groups (E9), the fitted exponents (E10), and the algebraic
> `S₄` oracle on the line (E11: flat in `p`, the first arm of this note
> whose fitted exponent sits below rho's `1/2`, and still not faster than
> rho at any size run).

## 1. What is new, against what exists

The repository already folds a factor base by an endomorphism in three
places, each written for one family:

| where | family | fold | what it established |
|:--|:--|:--|:--|
| `koblitz_index_calculus.rs`, `ic_boundary::koblitz_factor_base`, `koblitz-orbit` | Koblitz over `F_{2^n}` | signed Frobenius orbit, `2n` points a column, coefficient `±λ^k` | the ledger's one *advance* row: an orbit base needs `2n` fewer relations (`RESEARCH_IC_BOUNDARY_LEDGER.md`) |
| `glv_gaudry.rs`, `OrbitBase::fold_terms` | `j = 0` over `F_{p³}` with the subspace base | `⟨ψ⟩` orbit, `3` abscissae a column | `S ÷ 3.0` at every size with `residuals / floor ≈ 1.0`: engineering (`RESEARCH_GLV_INDEX_CALCULUS.md` §3) |
| `residual_walk.rs`, `Fold::Automorphism` | `j = 0` over `F_p` | sixfold residual fold | the prime-field sixfold of `RESEARCH_RESIDUAL_WALKS.md` §10.3 |
| `ec_index_calculus_j0.rs` | `j = 0` over `F_p` | orbit representatives only | documents that `ψ`-shifted rows are `λ^k` multiples of one row (no rank) |
| `cryptanalysis` repository, `prime_orbit_index_calculus.rs` | `j = 0`, `j = 1728`, generic over `F_p` | `Aut(E)` orbit, gain-graph solver | the `w²/4` coverage argument for pair sums |

What none of them is: **one fold, over any group the framework
counts, for any finite list of endomorphisms**, with the control the
same call.  That is what `fold_by_endomorphisms` is.  A family adds an
`Endomorphism` (a name, a degree, the eigenvalue `λ` modulo `r`, the
map on points) and a seed set; the fold closes the seed under the maps
or refuses a set that is not closed, walks the orbits, checks that two
paths to one point carry the same `Π λ`, and returns a
`FactorBase` with one column per orbit.  `verify_endomorphism` checks
every constructor's output on random points before it is used.  With
that in place the families this note covers are:

| type | family | endomorphism | `deg` | `ord_r(λ)` | invariant base | points a column | status |
|:--|:--|:--|--:|--:|:--|--:|:--|
| A | `j = 0`, `p ≡ 1 (3)` | `ψ(x, y) = (ζx, y)` | 1 | 3 | closure of any abscissa set under `x ↦ ζx` | 6 | implemented, `glv-orbit` |
| A | `j = 1728`, `p ≡ 1 (4)` | `ι(x, y) = (−x, iy)`, `ι² = −1` | 1 | 4 | closure under `x ↦ −x` | 4 | implemented, `glv-orbit` |
| A | generic | negation | 1 | 2 | any | 2 | implemented (the control family) |
| B | Koblitz, `F_{2^n}` | `τ` | 2 | `n` | `τ`-stable `F_2`-subspaces | `2n` | existing, `koblitz-orbit` |
| B | GLS twist over `F_{p²}` | `ψ = τ_u π_p τ_u⁻¹`, `ψ² = −1` | `p` | 4 | the line `x ∈ u·s·F_p` (`s² = ν`) | 4 | implemented, `gls-line` |
| B | subfield curve `E/F_p` on `E(F_{p^n})` | `π_p` | `p` | `n` | `π`-stable `F_p`-subspaces of `F_{p^n}` not inside `F_p` | `2n` | **pending** (§4, E5) |
| A×B | `j = 1728` twisted over `F_{p²}` (FourQ-shaped) | `⟨ι, ψ⟩` | — | 8 | a line stable under both | 8 | **pending** (E3); the fold call takes both generators as is |
| C | `D = −7` (`j = −3375`), `D = −8` (`j = 8000`) | degree-2 by Vélu, `φ² − tφ + 2 = 0`, `t ∈ {±1}`, `0` | 2 | `10³`–`10⁶` measured | **none** inside `⟨G⟩` (§2) | — | implemented as a measurement: `velu_degree2_endomorphisms`, `endomorphism_overlap` |
| C | `1 + i` on `j = 1728` | Vélu through `(0, 0)`, trace `±2` | 2 | large | none | — | implemented, verified in tests |
| C | degree 3 (`D = −11`; `√−3` on `j = 0`), Q-curves (Smith) | 3-isogeny by Vélu; `ψ² = ±d` | 3, `dp` | large | none | — | **pending** (E4b) |

The one derivation the module rests on, stated once: relations are
written over `[h]P`, and `[h]φ(P) = φ([h]P) = [λ][h]P`, so
`log [h]φ(P) = λ · log [h]P` for every base point, in or out of `⟨G⟩`.
The framework's cofactor curves (`h` up to 8 in the pilot) fold by the
same code path as the prime-order ones, and the test
`the_fold_is_sound_on_a_cofactor_curve` holds it to that.

## 2. Boundaries, stated before measuring

- **The fold is bounded by the orbit length, and the orbit length by
  `ord_r(λ)`.**  A base inside `⟨G⟩` that is stable under `φ = [λ]` is
  a union of orbits of multiplication by `λ` on `Z/rZ`, each of length
  `ord_r(λ)`.  For an automorphism of order `w` that is `w` (up to the
  fixed points, which are torsion); for the Frobenius on `E(F_{q^n})`
  it is `n`; for a degree-`d ≥ 2` CM endomorphism `λ` is a root of
  `λ² − tλ + d` modulo `r` and its order divides `r − 1` with nothing
  forcing it small.  So an invariant base under a type-C map is either
  essentially all of `⟨G⟩` or empty, and the most a type-C map can
  offer a base of `|F|` points is the chance overlap `|F|² / #E` of
  free relations `log φ(P) = λ log P`.  That is the boundary for E4,
  and the pilot measures it (§5.3).
- **Relation count.**  A square system needs about `columns` relations,
  hence `columns / ρ` targets at decomposition rate `ρ`; both arms see
  the same `ρ` because they hold the same points.  The fold divides
  `columns` by `w/2` (`w` the automorphism group's order), so the floor
  divides with it.  An *advance* would be `relations / columns < 0.9`
  on the folded arm with the control at `≈ 1.0`; the same ratio on
  both is the count moving with its bound (AGENTS.md §3).  One
  qualification the pilot forced into the open: the framework's
  relation loop stops when the incremental matrix first determines the
  target's logarithm, which on a two-summand base is the first cycle of
  the relation graph, at about half the columns; so what the pilot
  compares is relations-to-first-determination, not relations-to-full
  rank, and the ratio scatters around `w/2` from seed to seed.  E1's
  protocol fixes the stopping rule (§4).
- **Reference.**  A counted Pollard rho with the negation map
  (`rho_reference_negation`, `A = 2`) on the same instance and target,
  eight walks per row in the pilot.  The matched reference for a curve
  with `Aut` of order `w` is a rho that folds by all of `Aut`
  (`√(w/2)` cheaper; `aut_folded_rho.rs` has the `j = 0` walk in
  BigInt arithmetic, the framework does not have it in counted form):
  every "vs rho" figure in this note is therefore against the
  `A = 2` walk and understates the gap by up to `√3` on `j = 0` and
  `√2` on `j = 1728` and GLS.  Marked as such wherever it appears.
- **Falsification targets.**  E1/E2: an advance requires
  `relations / columns` to fall by more than the fold on the folded
  arm; the expectation, and the pilot's reading, is that it does not.
  E4: an overlap above three times the chance fraction on any family
  at any size would contradict §2's first boundary and would be the
  interesting result.  E3: fewer than `8` points a column on the
  `⟨ι, ψ⟩` line, or a fold that does not verify, falsifies the
  composite construction.  Inadmissible: changing the point set between
  the arms (the control **is** the folded arm's points), changing the
  oracle or the stopping rule between arms, dropping the base build or
  the table from `S`, counting `λ^k`-multiples of one row as
  independent, or a row whose answer is not the planted logarithm.

## 3. What the implementation checks so that a wrong number cannot get through

- Every constructor's map goes through `verify_endomorphism`: image
  on the curve, additive on random pairs, `φ([k]G) = [λ][k]G` on
  random `k`.  A perturbed eigenvalue is caught
  (`a_wrong_eigenvalue_is_caught_by_verification`).
- Instance orders are certified, never assumed: the CM candidate
  orders come from Cornacchia (`4p = t² + |D|v²`), the one that kills
  three random points is taken, and `[r]G = O` is checked on the
  generator.  Generic instances count points in `O(p)`.
- The fold refuses what it cannot fold, by name: a base a generator
  maps outside (`Strict`), a base point in a kernel, an orbit past
  `max_orbit` (an infinite-order generator), and two paths to one
  point that disagree on `Π λ`.  The type-C test asks for exactly that
  refusal.
- The GLS derivation is checked, not trusted: `ψ² = −1` on random
  points of the whole twist, `ψ_x(u s t) = −u s t` for every `t`, and
  the other `ψ`-stable line `u · F_p` carrying at most the 2-torsion.
- The shared relation loop recovers the planted logarithm over the
  folded base and over the control on `j = 0`, `j = 1728`, generic,
  `D = −7` (cofactor 4) and GLS
  (`a_folded_base_recovers_the_planted_logarithm`,
  `the_line_base_folds_four_to_one_and_recovers_the_logarithm`).

## 4. The plan

Every experiment below is a framework sweep: the same instance, the
same planted targets, the same seeds, the folded base and its control,
every phase inside `S`, a counted rho beside it.  Sizes are at least
four per family so an exponent can be fitted (AGENTS.md §5), seeds at
least six so the rho reference has its own spread
(`RESEARCH_GLV_INDEX_CALCULUS.md` §3 needed 48 walks to settle it).

| id | question | families | sizes | arms | measure | prediction | status |
|:--|:--|:--|:--|:--|:--|:--|:--|
| **E1** | Does the automorphism fold move anything but the column count on `F_p`? | `j0`, `j1728`, `generic` | `2^16`–`2^32`, six seeds | `glv-orbit` vs `no_fold`, oracles `subtract` and `mitm` | relations to **full rank** on one target stream (a second matrix fed the same rows until every column is pinned; the loop's first-determination count kept beside it), `S`, `S / rho` | `relations / columns ≈ 1.0` on both arms; `S` ratio `→ w/2` only where the relation phase dominates and `→ 1.0` where the pair table does (§5.2) | **done** (§6.1): column ratio exactly `3`, `2`; same rank fraction after `columns` relations on both arms at every size to `2^32`; engineering |
| **E2** | Does the GLS line fold match the Koblitz orbit fold in the same unit? | `gls` at `p = 2^8`–`2^16`; Koblitz `n = 13`–`31` from the ledger | six seeds | `gls-line` vs `no_fold`; `koblitz-orbit` vs `no_fold` | column ratio, relations, `S`, `S / rho`; the `descent-algebraic` analogue for `F_{p²}` (E2b) | `2×` on GLS against `n×` on Koblitz: the fold is the eigenvalue order, `4` against `2n`, and nothing else | **done** (§6.2): `4` and `2n` points a column, column ratios `2` and `n`, one driver; engineering |
| **E2b** | An `O(1)` decomposition oracle for the line: Weil-descend `S₃(u s t₁, u s t₂, x_R)` to two `F_p` equations in `(t₁, t₂)`, resultant, roots by `gcd(t^p − t, ·)` | `gls` | as E2 | `subtract` vs the resultant oracle on the same base | oracle cost per target, hit rate agreement target by target (AGENTS.md §6 cross-check) | same hits, `O(p)` fewer group operations per target; the fold ratio unchanged | **done** (§6.2): agrees with `subtract` on every target, `O(10³)` `F_p` multiplications a target against `|F|` additions; engineering |
| **E3** | Does a composite group fold as its order says? | `j = 1728` twisted over `F_{p²}` (`⟨ι, ψ⟩`, order 8); `j = 0` twisted over `F_{p²}` (`⟨ψ₃, ψ⟩`, order 12) | `p = 2^8`–`2^14` | the line stable under both vs negation | points a column, verification, `S` | `8` and `12` points a column; engineering | **done** (§6.3): `12` on `j = 0`; `4`, not `8`, on `j = 1728` (`ι = ±ψ` on `⟨G⟩`); the fold is the order of the eigenvalue subgroup of `(Z/rZ)^*` |
| **E4** | How much of a base does a type-C map keep in the base? | `d7`, `d8`, `1 + i` on `j1728` | `2^16`–`2^28` | one base per size, every rational degree-2 map | `ord_r(λ)`, `images_in_base / base_points` against `|F| / #E` | chance level; zero at pilot sizes | **done** (§6.4): `0` of `43,684` base points at `2^16`–`2^28`; accounting |
| **E4b** | The same for degree 3 (`D = −11`, `√−3` on `j = 0`) and for a Q-curve of degree 2 or 3 over `F_{p²}` | as named | as E4 | as E4 | chance level | **done for degree 3** (§6.4): `0` of `47,520` base points on `D = −11` and `j = 0`; Q-curves **not done** (§7) |
| **E5** | Subfield curves on `E(F_{p^n})`: a `π`-stable subspace not inside `F_p` | `E/F_p` on `E(F_{p³})`, with and without `j = 0` | `p = 2^6`–`2^11` | `⟨π⟩`, `⟨ψ⟩`, `⟨π, ψ⟩` folds vs negation on one base | points a column (`3`, `3`, `9`), relations, `S` | multiplicative; `glv_gaudry`'s `3.0×` is the `⟨ψ⟩` row | **done** (§6.5): `6` a column for `⟨−1, π⟩` (generic and `j = 0`) and for `⟨−1, π, ζ⟩` alike — `ζ` is a power of `π` on `⟨G⟩`, so not `9`; the `j = 0` line is the subgroup itself and its descent degenerates |
| **E6** | The matched reference | all folding families | as E1 | rho folded by the same `Aut` (counted) beside the `A = 2` walk | `S / rho_folded` | the fold's `w/2` in the relation count against rho's `√(w/2)` in steps: the gap widens by `√(w/2)` | **done** (§6.6): steps ratio `1.72` on `j = 0` (expected `1.73`), `1.51` on `j = 1728` (expected `1.41`), 384 walks, every answer verified |
| **E7** | Does the fold interact with a three-summand oracle? | `j0` | `2^18`–`2^30` | `mitm` at `m = 3` and the `S₄` oracle on folded vs control | relations, solver calls, rows that fold to `0 = 0` (the pair-generator finding of `RESEARCH_GLV_INDEX_CALCULUS.md` §4) | canonicalisation saves nothing on a uniform stream and everything on a pair sieve | **done** (§6.7): fold unchanged at `m = 3`; `0` rows fold to `0 = 0`; orbit duplicates at the birthday count |

What every experiment owes on completion, per AGENTS.md: the frozen
JSON under `experiments/23_glv_invariant_<id>.json`, the table printed
by `scripts/glv_invariant_tables.py`, the class of every row, the
scoreboard row, and — for anything claimed as a speedup — the frozen
regression suite of §8.

## 5. Pilot

Two seeds per size, `bits ∈ {16, 20, 24}` for the prime families
(`r` between `2^12.8` and `2^23.9`, cofactors up to 8) and
`p ∈ {2^8, 2^10, 2^12}` for GLS (`r` up to `2^22.9`); base seeded with
`max(8, r^{1/3})` abscissae and closed under the group (so `3×` as
many abscissae on `j = 0`, `2×` on `j = 1728`); GLS base the whole
line (`p − 1` abscissae); targets from the guarded r-adding walk;
oracles `subtract` (no table, `|F|` additions a target) and `mitm`
(pair table, negation-folded); eight `A = 2` rho walks per row.  Square
roots and Legendre symbols in the base build are counted and unpriced
(no pinned ratio for a generated curve, as in `ic bench`).  Every row
recovered its planted logarithm.

### 5.1 The fold, by family and oracle

| family | oracle | instances | column ratio | trial ratio (min–max) | S ratio (min–max) | S fold / rho S (min–max) |
|:--|:--|--:|--:|--:|--:|--:|
| d7 | mitm | 6 | 1.0 | 1.00–1.00 | 1.00–1.00 | 7×–24× |
| d7 | subtract | 6 | 1.0 | 1.00–1.00 | 1.00–1.00 | 36×–2718× |
| d8 | mitm | 6 | 1.0 | 1.00–1.00 | 1.00–1.00 | 6×–24× |
| d8 | subtract | 6 | 1.0 | 1.00–1.00 | 1.00–1.00 | 64×–3227× |
| generic | mitm | 6 | 1.0 | 1.00–1.00 | 1.00–1.00 | 5×–24× |
| generic | subtract | 6 | 1.0 | 1.00–1.00 | 1.00–1.00 | 54×–4250× |
| gls | subtract | 6 | 2.0 | 1.31–2.04 | 1.31–1.98 | 40×–796× |
| j0 | mitm | 6 | 3.0 | 1.00–2.66 | 1.00–1.02 | 21×–81× |
| j0 | subtract | 6 | 3.0 | 1.00–2.66 | 1.00–2.65 | 18×–1328× |
| j1728 | mitm | 6 | 2.0 | 1.00–4.74 | 1.00–1.13 | 13×–94× |
| j1728 | subtract | 6 | 2.0 | 1.00–4.74 | 1.00–3.82 | 16×–2834× |

The control families (`generic`, `d7`, `d8`) fold by negation on both
arms, so their two arms are the same run: a ratio of exactly `1.00` on
every row is the check that the harness changes nothing but the
column map.

The folding families, row by row:

| family | log2 r | h | oracle | signed points | columns fold | columns control | column ratio | relations fold | relations control | trials fold | trials control | S fold | S control | S control / S fold | rho S (A = 2) | S fold / rho S | correct |
|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| gls | 13.4 | 5 | subtract | 248 | 62 | 124 | 2.00 | 44 | 69 | 97 | 157 | 169.9 | 267.0 | 1.57 | 1.92 | 88× | yes |
| gls | 15.3 | 1 | subtract | 180 | 45 | 90 | 2.00 | 29 | 57 | 89 | 182 | 69.6 | 137.9 | 1.98 | 1.75 | 40× | yes |
| gls | 16.1 | 13 | subtract | 952 | 238 | 476 | 2.00 | 93 | 114 | 222 | 290 | 594.8 | 780.6 | 1.31 | 1.62 | 368× | yes |
| gls | 17.0 | 4 | subtract | 704 | 176 | 352 | 2.00 | 97 | 172 | 253 | 445 | 366.1 | 635.5 | 1.74 | 2.15 | 170× | yes |
| gls | 22.9 | 1 | subtract | 2748 | 687 | 1374 | 2.00 | 403 | 783 | 1112 | 2072 | 821.5 | 1512.6 | 1.84 | 1.03 | 796× | yes |
| gls | 22.9 | 1 | subtract | 2816 | 704 | 1408 | 2.00 | 428 | 864 | 1094 | 2169 | 802.4 | 1577.3 | 1.97 | 1.17 | 683× | yes |
| j0 | 13.8 | 3 | mitm | 144 | 24 | 72 | 3.00 | 16 | 34 | 72 | 156 | 52.0 | 52.9 | 1.02 | 2.48 | 21× | yes |
| j0 | 15.7 | 1 | mitm | 222 | 37 | 111 | 3.00 | 17 | 17 | 49 | 49 | 56.9 | 56.9 | 1.00 | 2.27 | 25× | yes |
| j0 | 18.0 | 3 | mitm | 378 | 63 | 189 | 3.00 | 32 | 87 | 354 | 941 | 75.0 | 76.4 | 1.02 | 1.56 | 48× | yes |
| j0 | 19.7 | 1 | mitm | 564 | 94 | 282 | 3.00 | 66 | 157 | 340 | 845 | 88.3 | 89.0 | 1.01 | 1.24 | 71× | yes |
| j0 | 20.5 | 7 | mitm | 678 | 113 | 339 | 3.00 | 79 | 132 | 3560 | 5597 | 101.8 | 104.0 | 1.02 | 1.49 | 69× | yes |
| j0 | 22.3 | 3 | mitm | 1032 | 172 | 516 | 3.00 | 86 | 224 | 2611 | 6245 | 119.9 | 121.6 | 1.01 | 1.47 | 81× | yes |
| j0 | 13.8 | 3 | subtract | 144 | 24 | 72 | 3.00 | 16 | 34 | 72 | 156 | 81.4 | 168.5 | 2.07 | 2.48 | 33× | yes |
| j0 | 15.7 | 1 | subtract | 222 | 37 | 111 | 3.00 | 17 | 17 | 49 | 49 | 40.0 | 40.0 | 1.00 | 2.27 | 18× | yes |
| j0 | 18.0 | 3 | subtract | 378 | 63 | 189 | 3.00 | 32 | 87 | 354 | 941 | 250.3 | 663.9 | 2.65 | 1.56 | 160× | yes |
| j0 | 19.7 | 1 | subtract | 564 | 94 | 282 | 3.00 | 66 | 157 | 340 | 845 | 182.0 | 456.5 | 2.51 | 1.24 | 147× | yes |
| j0 | 20.5 | 7 | subtract | 678 | 113 | 339 | 3.00 | 79 | 132 | 3560 | 5597 | 1972.6 | 3095.1 | 1.57 | 1.49 | 1328× | yes |
| j0 | 22.3 | 3 | subtract | 1032 | 172 | 516 | 3.00 | 86 | 224 | 2611 | 6245 | 1159.7 | 2769.0 | 2.39 | 1.47 | 788× | yes |
| j1728 | 13.3 | 4 | mitm | 84 | 21 | 42 | 2.00 | 12 | 12 | 153 | 153 | 27.5 | 27.5 | 1.00 | 2.07 | 13× | yes |
| j1728 | 13.8 | 4 | mitm | 96 | 24 | 48 | 2.00 | 3 | 13 | 31 | 147 | 27.1 | 29.2 | 1.07 | 1.89 | 14× | yes |
| j1728 | 16.6 | 8 | mitm | 180 | 45 | 90 | 2.00 | 26 | 50 | 1584 | 2481 | 37.5 | 41.1 | 1.10 | 1.33 | 28× | yes |
| j1728 | 17.9 | 4 | mitm | 244 | 61 | 122 | 2.00 | 16 | 61 | 604 | 2227 | 35.0 | 39.6 | 1.13 | 1.92 | 18× | yes |
| j1728 | 20.9 | 8 | mitm | 496 | 124 | 248 | 2.00 | 79 | 120 | 7944 | 11651 | 52.1 | 55.2 | 1.06 | 1.05 | 50× | yes |
| j1728 | 23.0 | 2 | mitm | 808 | 202 | 404 | 2.00 | 126 | 207 | 6483 | 10828 | 59.8 | 61.4 | 1.03 | 0.63 | 94× | yes |
| j1728 | 13.3 | 4 | subtract | 84 | 21 | 42 | 2.00 | 12 | 12 | 153 | 153 | 130.2 | 130.2 | 1.00 | 2.07 | 63× | yes |
| j1728 | 13.8 | 4 | subtract | 96 | 24 | 48 | 2.00 | 3 | 13 | 31 | 147 | 30.7 | 117.3 | 3.82 | 1.89 | 16× | yes |
| j1728 | 16.6 | 8 | subtract | 180 | 45 | 90 | 2.00 | 26 | 50 | 1584 | 2481 | 913.0 | 1417.7 | 1.55 | 1.33 | 685× | yes |
| j1728 | 17.9 | 4 | subtract | 244 | 61 | 122 | 2.00 | 16 | 61 | 604 | 2227 | 302.1 | 1100.6 | 3.64 | 1.92 | 158× | yes |
| j1728 | 20.9 | 8 | subtract | 496 | 124 | 248 | 2.00 | 79 | 120 | 7944 | 11651 | 2811.1 | 4116.2 | 1.46 | 1.05 | 2675× | yes |
| j1728 | 23.0 | 2 | subtract | 808 | 202 | 404 | 2.00 | 126 | 207 | 6483 | 10828 | 1792.9 | 2994.5 | 1.67 | 0.63 | 2834× | yes |

### 5.2 Reading the pilot

- **The column fold is exact and the points are the same.**  `6.0`,
  `4.0`, `4.0` and `2.0` signed points a column on `j = 0`, `j = 1728`,
  GLS and the controls, on every row, with the folded and control arms
  holding the same point set (the harness's `signed points` column is
  one number per row).
- **Relations track columns, with the stopping rule's scatter.**  The
  loop stops at the first determination of the target's logarithm,
  which on a two-summand base arrives at about half the columns
  (relations `16`–`428` against columns `24`–`704` on the folded arms).
  The ratio of relations control-to-fold is `1.0`–`3.8` on `j = 1728`
  and `1.0`–`2.7` on `j = 0`, around the column ratios `2` and `3`;
  the rows at `1.00` (`j0 15.7`, `j1728 13.3`) are runs where the same
  early target determined the logarithm on both arms.  This is why E1
  asks for relations to full rank on one stream: the first-cycle count
  is what the framework measures for a whole pipeline, and the
  full-rank count is what the fold's `w/2` is a statement about.
- **`S` moves only where the relation phase is the cost.**  With
  `subtract` (`|F|` additions a target) the relation phase is
  `> 97 %` of `S` and the `S` ratio is the trial ratio, `1.5`–`2.7` on
  `j = 0` and `1.5`–`3.8` on `j = 1728`.  With `mitm` the pair table
  is `85`–`99 %` of `S` at these sizes, both arms build the same table
  from the same points, and the `S` ratio is `1.00`–`1.13`: the fold
  saves trials that cost one lookup each.  Neither is a speedup claim;
  both are the same engineering fold seen through two oracles, and the
  phase table in `scripts/glv_invariant_tables.py` shows where each
  `S` sits.
- **Against rho, everything is far behind**, `5×`–`4,250×` on the
  `A = 2` walk, and by `√(w/2)` more against the matched folded walk
  (E6, extrapolation).  `S` grows with `r` on every family as the
  two-summand method must; no exponent is fitted on two seeds at three
  sizes, and none is claimed.

### 5.3 Type C, measured

| family | log2 r | rational degree-2 maps | ord_r(λ) | images in base | base points | chance fraction | verified |
|:--|--:|--:|--:|--:|--:|--:|:--|
| d7 | 12.8 | 4 | 6910 | 0 | 38 | 0.00069 | yes |
| d7 | 12.8 | 4 | 3554 | 0 | 38 | 0.00067 | yes |
| d7 | 17.5 | 4 | 30135 | 0 | 112 | 0.00015 | yes |
| d7 | 17.6 | 4 | 64002 | 0 | 114 | 0.00015 | yes |
| d7 | 21.5 | 4 | 243126 | 0 | 284 | 0.00002 | yes |
| d7 | 21.7 | 4 | 3424060 | 0 | 300 | 0.00002 | yes |
| d8 | 13.2 | 2 | 4613 | 0 | 40 | 0.00072 | yes |
| d8 | 14.9 | 2 | 14861 | 0 | 60 | 0.00101 | yes |
| d8 | 17.1 | 2 | 137218 | 0 | 102 | 0.00012 | yes |
| d8 | 17.9 | 2 | 33840 | 0 | 122 | 0.00013 | yes |
| d8 | 21.3 | 2 | 2545240 | 0 | 272 | 0.00002 | yes |
| d8 | 22.5 | 2 | 2969588 | 0 | 362 | 0.00003 | yes |

Vélu through the rational 2-torsion finds the degree-2 endomorphism on
every `D = −7` and `D = −8` instance (four maps on `D = −7`: two
kernels, each with `±φ`; two on `D = −8`), each verified as
`[λ]` on `⟨G⟩` with `λ² − tλ + 2 ≡ 0`, `t = ±1` and `0`.  The
eigenvalue order is `10³`–`10⁶` — comparable to `r` — and **no** base
point maps into the base at any size, against an expectation of
`0.00002`–`0.001` of them.  Both readings are §2's first boundary
holding: a type-C map has no invariant base to offer and no free
relations worth harvesting, and `fold_by_endomorphisms` refuses it by
name (`degree_two_cm_endomorphisms_exist_verify_and_do_not_fold`).
The `1 + i` map on `j = 1728` is the same object next to a type-A
automorphism: the automorphism folds, its degree-2 relative does not
(`one_plus_i_is_a_degree_two_endomorphism_of_a_j1728_curve`).

### 5.4 Verdict, for the pilot

| lever | what moved | class | test against the boundary |
|:--|:--|:--|:--|
| automorphism fold, `F_p` (`j = 0`, `j = 1728`) | columns `÷ 3`, `÷ 2` on the same points; relations and trials with them, with first-determination scatter; `S ÷ 1.5–3.8` under `subtract`, `÷ 1.0–1.1` under `mitm` | engineering | count moves with its floor; `18×`–`2,834×` the `A = 2` rho |
| GLS line fold, `F_{p²}` | columns `÷ 2`; relations `÷ 1.2–2.0`; `S ÷ 1.3–2.0` | engineering | the eigenvalue order `4` against Koblitz's `2n`: the same fold, a smaller group |
| degree-2 CM maps | `0` base points kept in the base, `ord_r(λ) ∈ [3.5·10³, 3.4·10⁶]` | accounting | the first boundary of §2 holds at every size |
| the unified fold itself | one function for five families; the control is the same call | engineering | — |

## 6. Experiments E1–E7: what was run and what it says

**Runner:** `cargo run --release --example glv_invariant_experiments -- --exp e<k> …` (the exact
commands are in each file's `command` field).  **Data:** `experiments/23_glv_invariant_e1.json`
(`2^16`–`2^28`, six seeds), `…_e1_32.json` (`2^32`, two seeds), `…_e2.json`, `…_e3.json`,
`…_e4.json`, `…_e5.json`, `…_e6.json`, `…_e7.json`, each with its `.log`.  **Tables:**
`python3 scripts/glv_invariant_experiment_tables.py experiments/23_glv_invariant_e*.json`; every
number below is printed by it.  Host: this container (Linux x86-64, 4 threads); wall time is
recorded in the files and is not a result.

### 6.0 Two things the runs taught the driver before they taught anything else

- **Full rank has a structural ceiling on cofactor curves.**  A two-summand relation
  `R = εP + ε'Q` with `R ∈ ⟨G⟩` forces the `E[h]`-components of `P` and `Q` to be negatives
  (`εu + ε'u' = 0`), so every odd function `χ` on `E[h]` — with `χ(gu) = λ_g χ(u)` on a folded
  arm — is a functional `Σ χ(u_i)x_i` no row can ever determine.  The rank saturates `D` short
  of `columns + 1`, `D` the number of classes of components that are not their own negatives
  (`glv_invariant_experiments.rs`, header).  The logarithm is pinned all the same, which is why
  the pilot never saw it.  The driver computes `D` from `[r]P` for every base point and calls a
  matrix **full** at `rank = columns + 1 − D`; the tables carry `D` for both arms.  Three-summand
  relations carry no such constraint.
- **Full rank is a coupon-collector count.**  The last column has to be touched and its component
  connected, so "relations to full rank" carries a `ln(columns)` factor on both arms and the ratio
  between arms overshoots the column ratio by the ratio of those logarithms.  The size-matched
  reading is the **rank fraction reached after exactly `columns` relations**; it is reported
  beside the full-rank count and is the number to read the floor against.  The framework's own
  first-determination count (the pilot's) is kept as well.
- **On a subfield curve the cofactor is fixed by the Frobenius.**  `E(F_p) ⊂ E(F_{p³})` is
  cofactor, `π` fixes it, so two-summand relations on the Frobenius line pair a point with its
  own orbit and the matrix is block-diagonal by `E(F_p)`-component on both arms (§6.5).  The
  driver reports single-column rows and `D` per arm so this reads off the table.

### 6.1 E1 — the automorphism fold at full rank

**Summary by family, oracle and size (means over seeds; per-row table: the script, 106 rows)**

| family | oracle | log2 r (mean) | seeds | column ratio | full-rank ratio mean (min–max) | full-rank rel / cols, fold | full-rank rel / cols, control | rank fraction at k = cols, fold | control | first-pin ratio mean | all correct |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| generic | mitm | 14.4 | 6 | 1.0 | 1.00 (1.00–1.00) | 1.94 | 1.94 | 0.84 | 0.84 | 1.00 | yes |
| generic | mitm | 18.4 | 6 | 1.0 | 1.00 (1.00–1.00) | 2.54 | 2.54 | 0.85 | 0.85 | 1.00 | yes |
| generic | mitm | 22.6 | 6 | 1.0 | 1.00 (1.00–1.00) | 3.14 | 3.14 | 0.83 | 0.83 | 1.00 | yes |
| generic | subtract | 14.4 | 6 | 1.0 | 1.00 (1.00–1.00) | 1.94 | 1.94 | 0.84 | 0.84 | 1.00 | yes |
| generic | subtract | 18.4 | 6 | 1.0 | 1.00 (1.00–1.00) | 2.54 | 2.54 | 0.85 | 0.85 | 1.00 | yes |
| j0 | mitm | 14.2 | 6 | 3.0 | 5.48 (2.23–8.88) | 2.17 | 3.47 | 0.79 | 0.84 | 2.63 | yes |
| j0 | mitm | 18.0 | 6 | 3.0 | 4.24 (3.05–5.77) | 2.10 | 2.83 | 0.82 | 0.84 | 4.20 | yes |
| j0 | mitm | 21.9 | 6 | 3.0 | 2.83 (2.35–3.76) | 3.30 | 3.03 | 0.83 | 0.85 | 2.95 | yes |
| j0 | mitm | 25.5 | 6 | 3.0 | 3.34 (2.73–4.12) | 3.48 | 3.83 | 0.84 | 0.84 | 2.43 | yes |
| j0 | mitm | 29.5 | 2 | 3.0 | 3.46 (3.15–3.77) | 3.53 | 4.07 | 0.84 | 0.84 | 2.50 | yes |
| j0 | subtract | 14.2 | 6 | 3.0 | 5.48 (2.23–8.88) | 2.17 | 3.47 | 0.79 | 0.84 | 2.63 | yes |
| j0 | subtract | 18.0 | 6 | 3.0 | 4.65 (3.39–5.77) | 2.00 | 2.98 | 0.82 | 0.84 | 4.17 | yes |
| j1728 | mitm | 13.4 | 6 | 2.0 | 2.22 (1.63–4.03) | 2.73 | 2.79 | 0.79 | 0.79 | 2.94 | yes |
| j1728 | mitm | 17.7 | 6 | 2.0 | 2.96 (1.25–4.07) | 2.73 | 3.51 | 0.81 | 0.81 | 2.27 | yes |
| j1728 | mitm | 21.8 | 6 | 2.0 | 2.41 (1.45–3.14) | 2.77 | 3.23 | 0.84 | 0.83 | 1.97 | yes |
| j1728 | mitm | 26.0 | 6 | 2.0 | 2.76 (2.32–3.09) | 3.12 | 4.30 | 0.84 | 0.84 | 1.98 | yes |
| j1728 | mitm | 30.6 | 2 | 2.0 | 2.63 (1.97–3.30) | 3.17 | 4.04 | 0.83 | 0.83 | 1.95 | yes |
| j1728 | subtract | 13.4 | 6 | 2.0 | 2.22 (1.63–4.03) | 2.73 | 2.79 | 0.79 | 0.79 | 2.94 | yes |
| j1728 | subtract | 17.7 | 6 | 2.0 | 2.96 (1.25–4.07) | 2.73 | 3.51 | 0.81 | 0.81 | 2.27 | yes |


Reading.  The column ratio is exactly `3` on `j = 0` and `2` on `j = 1728` at every size and
seed, and `1` on the generic control where the two arms are the same run.  After `columns`
relations both arms sit at the same rank fraction (folded and control within a few percent of
each other at every size), so the fold changes the number of unknowns and nothing about how fast
relations fill them: the count moves with its floor.  Relations to full rank fall by more than the
column ratio (the `ln(columns)` factor), and the first-determination count by about the column
ratio with the scatter the pilot showed.  `subtract` and `mitm` see the same relation stream
except where a target has more than one decomposition and the two oracles return different ones
(`2^18` on `j = 0`: `4.65` against `4.24`), so the oracle changes the cost of a target and little
else.
**Class: engineering**, as predicted; no `ratio-to-floor` fell.

### 6.2 E2 — GLS line against Koblitz orbit, in one driver

**Summary by family and eigenvalue order (per-row table: the script)**

| family | eigenvalue order | instances | log2 r | points per column | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, fold / control | oracle F_p muls / call | agreement with subtract | all correct |
|:--|--:|--:|:--|--:|--:|--:|--:|--:|:--|:--|
| gls | 4 | 18 | 9.3–24.0 | 4.0 | 2.0 | 2.20 (1.31–3.41) | 0.82 / 0.84 | 2606 | 3600/3600 (0 disagree) | yes |
| koblitz | 13 | 6 | 11.0–11.0 | 25.8 | 12.92 | — (—–—) | 0.45 / — | — | — | yes |
| koblitz | 15 | 6 | 7.7–7.7 | 19.4 | 9.8 | 12.00 (9.50–13.50) | 1.00 / — | — | — | yes |
| koblitz | 17 | 6 | 16.0–16.0 | 29.3 | 14.71 | 25.80 (17.11–36.22) | 0.90 / 0.85 | — | — | yes |
| koblitz | 23 | 6 | 22.0–22.0 | 45.0 | 22.52 | 45.69 (25.45–88.58) | 0.88 / 0.84 | — | — | yes |
| koblitz | 31 | 6 | 20.5–20.5 | 59.0 | 29.5 | 54.33 (28.09–83.70) | 0.70 / 0.58 | — | — | yes |


Reading.  On the GLS twist the fold is `4` points a column and the column ratio `2`; on the
Koblitz curves it is up to `2n` points a column (`19`–`59` on average, the abscissae in a proper
subfield having shorter orbits) and the column ratio `n` less that shortfall — `ord_r(λ)` in both
cases, which is the whole content of "Frobenius-type".  Both arms sit at the same rank fraction
after `columns` relations on the GLS twist (`0.82` against `0.84`) and on Koblitz `n = 17, 23`
(`0.90` / `0.85`, `0.88` / `0.84`); on `n = 31` the folded arm is ahead (`0.70` / `0.58`).  The
line oracle (E2b) agrees with `subtract` on every one of `3,600` targets across the small
instances and costs about `2,600` `F_p` multiplications a target where `subtract` costs
`|F| ≈ p` group additions, so the GLS arm runs at `p = 2^12` (`r ≈ 2^24`) in seconds.  On the
Koblitz arm only degrees whose `x^n − 1` has a factor of intermediate degree have a base of the
right size (`n = 15, 17, 23, 31`); `n = 13`'s only divisor gives a `2,003`-column base on a
group of order `2^11`, and its stream exhausts the group's targets before either arm reaches
full rank, so its row carries the fold (`25.8` a column) and no ratio.  The Koblitz controls
carry the coupon-collector factor of a `500`–`1,000`-column base, which is why their full-rank
ratios (`12`–`54`) exceed `n`.  **Class: engineering** on both arms.

### 6.3 E3 — composite groups on twisted CM curves

**Every instance**

| family | p | log2 r | h | ord λ_ψ | ord λ_aut | aut = ±ψ on ⟨G⟩ | pts/col negation | pts/col ψ | pts/col aut | pts/col ψ+aut | cols ψ+aut | square ratio ψ+aut vs negation | square ratio ψ vs negation | correct |
|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|
| j0 | 2^8 | 13.8 | 4 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 20 | 10.33 | 2.10 | yes |
| j0 | 2^8 | 15.2 | 1 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 20 | 10.63 | 2.30 | yes |
| j0 | 2^8 | 12.7 | 4 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 15 | 5.27 | 2.30 | yes |
| j0 | 2^8 | 11.0 | 25 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 16 | 13.18 | 1.07 | yes |
| j0 | 2^10 | 16.5 | 4 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 54 | 10.82 | 2.77 | yes |
| j0 | 2^10 | 17.1 | 4 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 63 | 3.99 | 1.24 | yes |
| j0 | 2^10 | 17.6 | 4 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 73 | 12.72 | 1.90 | yes |
| j0 | 2^10 | 18.0 | 4 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 84 | 5.70 | 2.60 | yes |
| j0 | 2^12 | 19.2 | 13 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 240 | 7.14 | 1.76 | yes |
| j0 | 2^12 | 20.0 | 13 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 309 | 8.77 | 2.20 | yes |
| j0 | 2^12 | 19.8 | 13 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 273 | 7.95 | 1.85 | yes |
| j0 | 2^12 | 23.0 | 1 | 4 | 3 | False | 2.0 | 4.0 | 6.0 | 12.0 | 236 | 6.34 | 1.99 | yes |
| j1728 | 2^8 | 6.8 | 340 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 56 | — | — | folds only (r < 16·columns) |
| j1728 | 2^8 | 6.7 | 324 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 50 | — | — | folds only (r < 16·columns) |
| j1728 | 2^8 | 6.7 | 324 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 50 | — | — | folds only (r < 16·columns) |
| j1728 | 2^8 | 6.2 | 340 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 36 | — | — | folds only (r < 16·columns) |
| j1728 | 2^10 | 8.4 | 1220 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 174 | — | — | folds only (r < 16·columns) |
| j1728 | 2^10 | 6.6 | 9860 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 240 | — | — | folds only (r < 16·columns) |
| j1728 | 2^10 | 8.8 | 1780 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 216 | — | — | folds only (r < 16·columns) |
| j1728 | 2^10 | 8.8 | 1508 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 216 | — | — | folds only (r < 16·columns) |
| j1728 | 2^12 | 8.1 | 28260 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 690 | — | — | folds only (r < 16·columns) |
| j1728 | 2^12 | 10.3 | 5380 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 644 | — | — | folds only (r < 16·columns) |
| j1728 | 2^12 | 7.1 | 34112 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 544 | — | — | folds only (r < 16·columns) |
| j1728 | 2^12 | 11.0 | 7780 | 4 | 4 | True | 2.0 | 4.0 | 4.0 | 4.0 | 1034 | — | — | folds only (r < 16·columns) |


Reading.  The composite fold is not the product of the generators' orders; it is the **order of
the subgroup of `(Z/rZ)^*` their eigenvalues generate**.  On a `j = 0` twist `λ_ζ` (order 3) and
`λ_ψ` (order 4) generate the twelfth roots of unity and the line folds `12` to a column, exactly.
On a `j = 1728` twist `ψ² = −1 = ι²` and there are only two square roots of `−1` modulo `r`, so
`ι = ±ψ` on `⟨G⟩` (measured on every instance) and `⟨−1, ψ, ι⟩` folds `4`, as `⟨−1, ψ⟩` alone
does.  A `j = 1728` GLS twist also always carries a cofactor of about `p` (its `ψ` and `ι`
coincide up to sign as maps of the subgroup, and the twist is isogenous to a curve over `F_p`),
so its usable `r` is only about `p`; the rows say so.  The `j = 0` composite arm recovers the
logarithm over the line oracle with the control on the same points.  This corrects the plan's
prediction of `8` for the `j = 1728` composite.  **Class: engineering** on `j = 0`;
**accounting** on `j = 1728` (the composite adds nothing).

### 6.4 E4 — type C, degree 2 and degree 3

**Summary (per-map table: the script, 80 rows)**

| family | degree | instances | maps per instance | ord_r(λ) min–max | images in base, total | base points, total | all verified |
|:--|--:|--:|--:|--:|--:|--:|:--|
| d11 | 3 | 16 | 4 | 33–245977410 | 0 | 22056 | yes |
| d7 | 2 | 16 | 4 | 337–37072938 | 0 | 15568 | yes |
| d8 | 2 | 16 | 2 | 4613–28560484 | 0 | 9484 | yes |
| j0 | 3 | 16 | 6 | 263–66152880 | 0 | 25464 | yes |
| j1728 | 2 | 16 | 4 | 41–124755636 | 0 | 18632 | yes |


Reading.  Vélu finds every rational degree-2 map on `D = −7`, `D = −8` and `1 + i` on
`j = 1728`, and every rational degree-3 map on `D = −11` and `√−3` on `j = 0` (six on `j = 0`:
three cube roots of the isomorphism scaling, two signs), each verified as `[λ]` on `⟨G⟩`.  The
eigenvalue orders run from `10²` to `2·10⁸` — comparable to `r` at every size — and **no** base
point maps into the base on any instance, `2^13`–`2^28`, against a chance expectation of
`10⁻⁵`–`10⁻³` per point.  §2's boundary holds at every size for both degrees: no type-C map
folds, and none offers free relations worth harvesting.  **Class: accounting.**

### 6.5 E5 — subfield curves on E(F_{p³})

**Summary by family (per-instance table: the script)**

| family | instances | log2 r | pts/col negation | pts/col π | pts/col ζ | pts/col π+ζ | ζ ∈ ⟨π⟩ on ⟨G⟩ | streams | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | oracle muls/call | degenerate | all correct |
|:--|--:|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|
| subfield | 20 | 10.7–21.7 | 2.0 | 6.0 | — | — | — | 20 | 1.61 (1.19–3.31) | 0.80 / 0.83 | 1997 | 0 | yes |
| subfield-j0 | 20 | 6.2–10.8 | 2.0 | 6.0 | 6.0 | 6.0 | all | 0 | — (—–—) | — / — | — | — | folds only |


Reading.  On a generic subfield curve the Frobenius eigenline `x ∈ s^k·F_p` (the line whose
points survive `[h]`; the other line carries the other cubic twist) folds `6` to a column,
column ratio `3`, and the line oracle over `F_{p³}` recovers the logarithm on both arms — the
odd-characteristic twin of the Koblitz fold, at `ord_r(λ_π) = 3`.  But the relation count does
**not** fall with the columns: relations to full rank fall by `1.2`–`1.6` against a column ratio
of `3`, and after `columns` relations the folded arm is at `0.60`–`0.68` of its full rank against
the control's `0.81`–`0.92` (the rows with `h = #E(F_p)`).  The files say why.  `E(F_p)` is a
subgroup of `E(F_{p³})` and never of `⟨G⟩`, so it is always in the cofactor, and `π` fixes it:
every point of a `π`-orbit has the **same** `E(F_p)`-component.  A two-summand relation
`R = εP + ε'Q` forces those components to be negatives (§6.0), so — up to the birthday
coincidences between orbits — `Q = −π^i P` and the relation is `(1 − λ^i)·x_c = log R`: a
**single-column row**.  The tables carry the count: on the folded arm `1,338` of the `1,339`
relations of the `2^21.3` instance pin one column each (`D = 0`, nothing is void), and on the
control the matrix is block-diagonal with a block of three columns and one undetermined
functional per block (`D = 232` against `264` orbits: a block per distinct `±` component, a few blocks holding two orbits that share one).  Both arms are then coupon collectors over
the **same** set of blocks — the fold merges each block's three unknowns into one and lowers the
rank a block needs from two to one, which is what the `1.2`–`1.6` is — and the rank fractions
are those of a coupon collector at one row per coupon (`1 − 1/e ≈ 0.63`) against two rows per
coupon at three per block (`≈ 0.80`).  With a cofactor beyond `E(F_p)` (`h / #E(F_p) = 4, 7`)
the extra component is not fixed by `π`, `D` grows on both arms, and the fold reaches full rank
before `columns` relations.  This is a property of two-summand decomposition on any subfield
curve, not of the fold: Gaudry's `E(F_{p³})` index calculus uses three summands for exactly the
size reason and inherits none of it; the three-summand line oracle was not run here (§7).

On a `j = 0` subfield curve the picture collapses further: `N(π² + π + 1)` splits as
`N(π − ω)·N(π − ω²)` in `Z[ω]`, so the "new" part of `E(F_{p³})` is two subgroups of about
`p` each, `r ≈ p`, and the eigenline carrying `⟨G⟩` is the `F_p`-points of a cubic twist — the
base **is** the subgroup, every target decomposes in about `p/2` ways, `x_R` lies on the line
and the descent degenerates.  The fold still measures `6` for `⟨−1, π⟩`, `⟨−1, ζ⟩` and
`⟨−1, π, ζ⟩` alike, because `λ_ζ ∈ {λ_π, λ_π²}` on every instance: `ζ` is a power of `π` on
`⟨G⟩`, the same rule as E3.  This corrects the plan's prediction of `9` for `⟨π, ζ⟩`.
**Class: engineering** for the columns on the generic curve, **accounting** for the relation
count (the block structure, not the fold, sets it) and for `j = 0`.

### 6.6 E6 — the matched folded rho

**Summary (per-instance table: the script, 48 rows)**

| family | instances | walks | S ratio mean (min–max) | steps ratio mean | expected | all verified |
|:--|--:|--:|--:|--:|--:|:--|
| j0 | 24 | 192 | 1.58 (1.06–3.01) | 1.72 | 1.73 | True |
| j1728 | 24 | 192 | 1.41 (0.96–2.23) | 1.51 | 1.41 | True |


Reading.  The rho folded by the same group verifies its answer on every walk and takes about
`√3` (`j = 0`) and `√2` (`j = 1728`) fewer steps than the negation walk, with the spread eight
walks per instance give.  That is the `√(A/2)` a generic algorithm already takes from the
group, and it is what every "vs rho" figure of §5 understated: the fold's `w/2` in the relation
count against rho's `√(w/2)` in steps widens the gap by `√(w/2)`.  **Class: accounting** for the
reference.

### 6.7 E7 — three summands and orbit duplicates

**Every instance**

| log2 r | seed | m | group order | seed abscissae | points | cols fold | cols control | square rel fold | square rel control | ratio | zero-support rows fold / control | single-column rows fold / control | orbit duplicates | targets | hit rate | correct |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|
| 17.2 | 3 | 2 | 6 | 19 | 114 | 19 | 57 | 35 | 124 | 3.54 | 0 / 0 | 5 / 1 | 148 | 3101 | 0.0400 | yes |
| 17.2 | 3 | 3 | 6 | 19 | 114 | 19 | 57 | 30 | 58 | 1.93 | 0 / 0 | 0 / 0 | 0 | 70 | 0.8286 | yes |
| 17.4 | 4 | 2 | 6 | 20 | 120 | 20 | 60 | 24 | 185 | 7.71 | 0 / 0 | 7 / 2 | 289 | 4800 | 0.0385 | yes |
| 17.4 | 4 | 3 | 6 | 20 | 120 | 20 | 60 | 46 | 85 | 1.85 | 0 / 0 | 0 / 0 | 0 | 103 | 0.8252 | yes |
| 17.7 | 2 | 2 | 6 | 21 | 126 | 21 | 63 | 52 | 159 | 3.06 | 0 / 0 | 4 / 0 | 239 | 4679 | 0.0340 | yes |
| 17.7 | 2 | 3 | 6 | 21 | 126 | 21 | 63 | 41 | 104 | 2.54 | 0 / 0 | 0 / 0 | 1 | 133 | 0.7820 | yes |
| 17.7 | 1 | 2 | 6 | 21 | 126 | 21 | 63 | 32 | 145 | 4.53 | 0 / 0 | 10 / 3 | 193 | 3966 | 0.0366 | yes |
| 17.7 | 1 | 3 | 6 | 21 | 126 | 21 | 63 | 25 | 143 | 5.72 | 0 / 0 | 0 / 0 | 0 | 181 | 0.7901 | yes |
| 18.9 | 4 | 2 | 6 | 26 | 156 | 26 | 78 | 74 | 288 | 3.89 | 0 / 0 | 14 / 4 | 24149 | 82085 | 0.0035 | yes |
| 18.9 | 4 | 3 | 6 | 26 | 156 | 26 | 78 | 59 | 80 | 1.36 | 0 / 0 | 0 / 0 | 0 | 453 | 0.1766 | yes |
| 20.0 | 3 | 2 | 6 | 32 | 192 | 32 | 96 | 84 | 193 | 2.30 | 0 / 0 | 14 / 3 | 2946 | 36381 | 0.0053 | yes |
| 20.0 | 3 | 3 | 6 | 32 | 192 | 32 | 96 | 34 | 155 | 4.56 | 0 / 0 | 0 / 0 | 0 | 449 | 0.3452 | yes |
| 20.1 | 1 | 2 | 6 | 32 | 192 | 32 | 96 | 71 | 267 | 3.76 | 0 / 0 | 15 / 1 | 5310 | 52344 | 0.0051 | yes |
| 20.1 | 1 | 3 | 6 | 32 | 192 | 32 | 96 | 33 | 165 | 5.00 | 0 / 0 | 1 / 1 | 1 | 581 | 0.2840 | yes |
| 20.3 | 2 | 2 | 6 | 33 | 198 | 33 | 99 | 75 | 312 | 4.16 | 0 / 0 | 12 / 1 | 6583 | 62178 | 0.0050 | yes |
| 20.3 | 2 | 3 | 6 | 33 | 198 | 33 | 99 | 42 | 165 | 3.93 | 0 / 0 | 0 / 0 | 0 | 590 | 0.2797 | yes |
| 22.8 | 4 | 2 | 6 | 51 | 306 | 51 | 153 | 170 | 602 | 3.54 | 0 / 0 | 18 / 9 | 140858 | 714665 | 0.0008 | yes |
| 22.8 | 4 | 3 | 6 | 51 | 306 | 51 | 153 | 60 | 227 | 3.78 | 0 / 0 | 0 / 0 | 6 | 2385 | 0.0952 | yes |
| 23.1 | 2 | 2 | 6 | 54 | 324 | 54 | 162 | 80 | 303 | 3.79 | 0 / 0 | 10 / 2 | 35128 | 373187 | 0.0008 | yes |
| 23.1 | 2 | 3 | 6 | 54 | 324 | 54 | 162 | 102 | 368 | 3.61 | 0 / 0 | 0 / 0 | 5 | 4134 | 0.0890 | yes |
| 24.3 | 1 | 2 | 6 | 67 | 402 | 67 | 201 | 115 | 701 | 6.10 | 0 / 0 | 16 / 5 | 36476 | 567338 | 0.0012 | yes |
| 24.3 | 1 | 3 | 6 | 67 | 402 | 67 | 201 | 132 | 425 | 3.22 | 0 / 0 | 0 / 0 | 0 | 2580 | 0.1647 | yes |
| 25.7 | 3 | 2 | 6 | 85 | 510 | 85 | 255 | 188 | 664 | 3.53 | 0 / 0 | 5 / 0 | 3524 | 275060 | 0.0024 | yes |
| 25.7 | 3 | 3 | 6 | 85 | 510 | 85 | 255 | 124 | 641 | 5.17 | 0 / 0 | 0 / 0 | 1 | 1839 | 0.3486 | yes |


Reading.  With three summands the fold on `j = 0` divides the columns by `3` as before and the
full-rank relation count with them; the deficiency of §6.0 does not apply (`h = 1` here in any
case).  Rows whose base part folds to nothing (`0 = ha + hbd`) do not occur on a uniform target
stream; single-column rows do (a target that is a sum of points of one orbit), a few percent on
the folded arm.  Orbit duplicates among the targets — a target whose `⟨−1, ζ⟩`-orbit repeats an
earlier one's, the solver call canonicalisation would save — are the birthday count
`3T²/r` and nothing more, as `RESEARCH_GLV_INDEX_CALCULUS.md` §4 found on the `F_{p³}` harness.
**Class: accounting.**

### 6.8 Verdict after E1–E7

| lever | measured | class |
|:--|:--|:--|
| automorphism fold, `F_p` (E1) | columns `÷ 3`, `÷ 2`; same rank fraction after `columns` relations on both arms at every size to `2^32` | engineering |
| GLS line (E2) | `4` points a column, column ratio `2`, same rank fraction after `columns` relations on both arms; Koblitz up to `2n` on the same driver | engineering |
| subfield line (E5) | `6` points a column, column ratio `3`; relations to full rank `÷ 1.2–1.6` only, every folded row a single-column row: the cofactor `E(F_p)` is `π`-fixed and two-summand relations stay inside one orbit | engineering (columns) / accounting (relations) |
| composite groups (E3, E5) | the fold is the order of the eigenvalue subgroup of `(Z/rZ)^*`: `12` on the `j = 0` twist, `4` on the `j = 1728` twist, `3` on the `j = 0` subfield curve | engineering / accounting |
| type C, degree 2 and 3 (E4) | `0` base points kept in the base at every size; eigenvalue orders `10²`–`10⁸` | accounting |
| matched folded rho (E6) | `√3`, `√2` fewer steps, verified | accounting (reference) |
| three summands, orbit duplicates (E7) | fold unchanged; duplicates at the birthday count | accounting |
| the line oracle (E2b) | about `2,600` `F_p` multiplications a target on `F_{p²}`, `2,000` on `F_{p³}`; agrees with `subtract` on `3,600` of `3,600` targets | engineering |

Nothing here moves a ratio to a floor.  The fold by any finite group of endomorphisms is the
order of the subgroup of `(Z/rZ)^*` that the group's eigenvalues generate — no more, whatever
the group of maps looks like — and rho takes the square root of the same number.  Where the
group has a cofactor the endomorphisms fix, the fold buys columns and not relations (E5).  That is the
closed form the plan's E3 and E5 predictions lacked, and it is the reason the search for a
larger fold on a fixed prime-order subgroup ends here: the roots of unity modulo `r` that
endomorphisms of a curve can realise are the automorphism group's, the Frobenius's, and their
products.

## 7. What was not done

- **Q-curves (E4b, second half).**  No degree-2 or degree-3 Q-curve over
  `F_{p²}` was built; Smith's construction is not in the repository and
  the type-C boundary is only measured on CM curves over `F_p`.
- **E2's Koblitz arm covers four degrees.**  The ledger's
  `koblitz_instance` roster has usable divisors of `x^n − 1` of
  intermediate degree at `n = 15, 17, 23, 31` only; `n = 13` has a hit
  rate too low to reach full rank within the trial cap and `n = 19` a
  subspace too large for the pair table.  Four sizes is the minimum the
  plan asked for and no exponent is fitted to them.
- **E5 on `j = 0` has no relation stream.**  The eigenline carrying
  `⟨G⟩` is the subgroup itself, so the descent is degenerate and only
  the fold is measured; a `j = 0` subfield curve with `r` larger than
  `p` would need `E(F_{p^n})` for `n` prime to `3` or a different
  eigenline, neither of which the driver builds.
- **The algebraic `S₄` oracle (E11, §8.5) is run to `p = 2^13` only**,
  and its Frobenius symmetry breaking is measured as the orbit-duplicate
  count rather than reproduced from the paper, which this session could
  not retrieve.
- **`S` at full rank is measured for E8 only.**  §8.2 prices every
  phase of the three-summand pair-table pipeline and the matched rho in
  one unit on the subfield curves; E1–E7 count relations, rank, trials,
  oracle cost and rho steps, and their end-to-end `S` columns stay the
  pilot's (§5.1).  No row anywhere claims a speed.
- **The line oracle is a library and example component only.**  It is
  not a framework plugin and `ic bench` cannot select it; the `subtract`
  and `mitm` oracles remain what the CLI offers on the line bases.
- **The `F_{p²}` and `F_{p³}` "vs rho" figures of §5 are against the
  `A = 2` walk**; E9 (§8.3) measures the matched folded walk on those
  groups, so the reader can rescale them, but the §5 rows are not
  re-run.
- The pilot's base build prices no square roots or Legendre symbols
  (no pinned ratio for generated curves); at these sizes the build is
  under `3 %` of `S` on every row, but a larger sweep should pin them
  as `ic bench` does for the roster.
- The scoreboard carries the pilot and E1–E9 as two panels with their
  column, rank and relation ratios; no scoreboard row claims a speed,
  and the exponent panel is untouched.
- Nothing here bears on any deployed curve: `r ≤ 2^32`, certified toy
  instances, and a fold that rho already takes as `√(w/2)`.
- E12–E15 of §8.1 are proposed, not run.

## 8. Toward the state of the art: what the literature does, what is next, and E8–E9

### 8.0 Where E1–E7 stand against the literature

Three published lines of work bear on an endomorphism-invariant factor base, and
§6 reproduces two of them at framework level and adds two things they do not
report.

- **Galbraith, Granger, Merz, Petit, "On Index Calculus Algorithms for Subfield
  Curves"** (SAC 2020; ePrint 2020/1315).  Factor bases invariant under the
  `q`-power Frobenius on subfield (Koblitz) curves: the number of decomposition
  systems to solve falls by `1/n`, the linear algebra by `n²`, and the Frobenius
  is used for *symmetry breaking* in the polynomial systems.  E2 (Koblitz orbit,
  GLS line) and E5 (the `F_{p³}` line) are the odd-characteristic, counted
  reproduction of the first two claims; §6.8's closed form — the fold is the
  order of the eigenvalue subgroup of `(Z/rZ)^*`, so `⟨π, ζ⟩` folds `3` and not
  `9` — and §6.5's block degeneracy of two-summand relations under a
  `π`-fixed cofactor are not in that paper.  Its third lever, symmetry breaking
  inside the decomposition system, is **not** implemented here (E11 below).
- **Faugère, Gaudry, Huot, Renault, "Using Symmetries in the Index Calculus for
  Elliptic Curves Discrete Logarithm"** (ePrint 2012/199; J. Cryptology 2014),
  and the follow-up on symmetrised summation polynomials from small-order
  torsion.  These act on the *decomposition system*: the symmetric group `S_m`
  on the summands and, from a rational `2`-torsion point, `(Z/2)^{m−1}`, give an
  invariant-ring presentation whose Gröbner cost falls by an exponential factor
  in `m`; on binary Koblitz curves the Frobenius is added.  An endomorphism `φ`
  of `⟨G⟩` is **not** a symmetry of that system: `R = ΣP_i` does not give
  `R = Σ φ^{a_i} P_i` unless every `a_i` is equal, in which case the target
  moves.  So the fold of this note and the FGHR symmetries are different levers
  on different objects — the base and the system — and multiply (E13).
- **Chi-Domínguez, Rodríguez-Henríquez, Smith, "Extending the GLS endomorphism
  to speed up GHS Weil descent using Magma"** (2021; `F_{2^{155}}` solved).  The
  GLS endomorphism induces an endomorphism of the GHS Jacobian and gives a
  factor `n` there; it is the same `ord_r(λ)` lever on a different group, and
  reads as further evidence that the lever is a constant.
- In this repository, `RESEARCH_GLV_INDEX_CALCULUS.md` §3 measured the `⟨ψ⟩`
  quotient on `E(F_{p³})` with the algebraic `S₄` oracle: `S ÷ 3.0` with the
  count on its floor — the same reading as §6 with the other oracle.

The honest position: no published or measured use of an endomorphism moves the
exponent.  The state of the art is the **combination** — an invariant base (GGMP,
this note) under an FGHR-symmetrised decomposition system with Frobenius symmetry
breaking — and nothing here or in the literature has measured that combination
on an odd-characteristic subfield curve at full rank with every phase priced.
That is the ordering below.

### 8.1 Ranked next steps

| id | question | prediction | falsification target (in advance) | class if it holds | status |
|:--|:--|:--|:--|:--|:--|
| **E8** | Three summands on the `π`-line of a subfield curve: does the fold's relation count fall with the columns once the §6.5 block degeneracy is gone? | relations to full rank `÷ 3` (times the coupon factor), `0` single-column rows, both arms at one rank fraction after `columns` relations | full-rank ratio `< 2.0` on a majority of instances, or single-column rows `> 1 %`, falsifies | engineering | **done** (§8.2) |
| **E9** | The matched folded rho on the `F_{p²}` and `F_{p³}` groups, so every "vs rho" on those groups has its reference | steps ratio `√(A/2)`: `√2` (GLS), `√6` (`j = 0` twist), `√3` (subfield) | a mean steps ratio below `0.8·√(A/2)` or any unverified walk falsifies | accounting (reference) | **done** (§8.3) |
| **E10** | One `S` column in one unit (`F_p` multiplications; inversions by a measured factor; linear-algebra row operations as one each) for both arms and both walks on one instance, five sizes, exponent fitted | the fold moves `S` by its column ratio where the relation phase dominates and by its square where the linear algebra does; the exponent of neither arm falls below rho's `1/2` with the pair-table oracle | an arm whose fitted exponent is below the control's by more than the fit's scatter would be an advance; expected: none | engineering | **done**, inside E8 (§8.2) |
| **E11** | The algebraic three-summand oracle on the line: `S₄(s t₁, s t₂, s t₃, x_R)` Weil-descended, `S₃`-symmetrised in the elementary symmetric functions of `t_i`, solved as `gaudry_cubic::solve_s4_subspace` solves the `x_i ∈ F_p` base; then GGMP's Frobenius symmetry breaking — only canonical `⟨−1, π⟩`-orbit representatives admitted as solutions | oracle cost `O(1)` in `p` (the Macaulay solve at degree `10`–`13`); the fold ratio unchanged at `3`; symmetry breaking removes the `3`-fold redundancy among a system's solutions, not systems | the fold's `S` ratio below `2.7` at any size with this oracle, or an oracle cost growing with `p`, falsifies | engineering | **done** (§8.5): `9.42e+05` multiplications a call, flat in `p`; agrees with the pair table on `3600` of `3600` targets; fitted exponent `0.53` on the fold; symmetry breaking reduces to the orbit duplicates |
| **E12** | The pair table over orbit representatives: `P + φ^c P'` for representatives `P, P'` and `c < w/2` | table `w/2` smaller and `w/2` cheaper to build, one probe per target unchanged; `S` moves only by the table's share | a probe cost above `1.2×` the full table's falsifies "unchanged" | engineering (memory) | pending |
| **E13** | FGHR's `2`-torsion symmetry **and** the fold together: close the line base under translation by `T ∈ E(F_p)[2]` and under `π`, fold by `⟨−1, π⟩`, symmetrise the system by `(Z/2)^{m−1} ⋊ S_m` | the two levers multiply — columns `÷ 3`, system degree `÷ 2^{m−1}` — because one acts on the base and the other on the system | a combined `S` ratio below the product of the separate ratios by more than `20 %` falsifies "multiply" | engineering | pending; needs E11's oracle and a `+T`-closed base |
| **E14** | Q-curves of degree 2 and 3 over `F_{p²}` (the second half of E4b) | type C: `0` base points kept, no fold | any image in the base beyond chance falsifies | accounting | pending; Smith's construction not in the repository |
| **E15** | Transfer to the binary Koblitz program (ECC2K-130, AGENTS.md §8a–8b): is the §6.5 degeneracy present there? | no: `E(F_2) ⊂ E(F_{2^n})` has order `2` or `4`, so at most two component classes and no block structure; the fold there is the known `2n` | a measured deficiency `D > 2` on a prime-order-times-`4` Koblitz subgroup with two summands falsifies | accounting | pending; one `m = 31` stream with the deficiency column suffices |

What is **not** on the list, and why: a larger fold on a fixed prime-order
subgroup.  §6.8's closed form bounds it by the roots of unity the curve's
endomorphisms realise modulo `r` — `lcm(|Aut E|, ord λ_π)` at most — and every
family with `|Aut E| > 2` in odd characteristic is `j = 0` or `j = 1728`; there
is no fourth lever on the base.

### 8.2 E8 — three summands on the Frobenius line, every phase priced

**Relations, every instance**

| p | log2 r | h / #E(F_p) | cols π | cols negation | full-rank rel π | full-rank rel negation | ratio | rank fraction at k = cols, π / negation | single-column rows π / negation | orbit duplicates | targets | hit rate | correct |
|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|
| 2^7 | 10.7 | 4 | 11 | 33 | — | — | — | 0.75 / 0.74 | 0 / 0 | 1410 | 10180 | 0.0094 | yes |
| 2^8 | 11.8 | 7 | 26 | 78 | 50 | 181 | 3.62 | 0.85 / 0.90 | 0 / 0 | 598 | 1325 | 0.1366 | yes |
| 2^7 | 11.9 | 4 | 19 | 57 | — | — | — | 0.75 / 0.72 | 0 / 0 | 3110 | 31286 | 0.0084 | yes |
| 2^7 | 12.1 | 1 | 14 | 42 | 19 | 94 | 4.95 | 0.80 / 0.91 | 0 / 0 | 59 | 353 | 0.2663 | yes |
| 2^8 | 12.4 | 7 | 33 | 99 | 42 | 226 | 5.38 | 0.88 / 0.90 | 0 / 0 | 647 | 1637 | 0.1381 | yes |
| 2^7 | 13.4 | 1 | 18 | 54 | 46 | 116 | 2.52 | 0.95 / 0.89 | 0 / 0 | 110 | 755 | 0.1536 | yes |
| 2^9 | 14.8 | 7 | 87 | 261 | 100 | 672 | 6.72 | 0.98 / 0.91 | 0 / 0 | 658 | 3218 | 0.2088 | yes |
| 2^9 | 15.1 | 4 | 62 | 186 | 99 | 506 | 5.11 | 0.95 / 0.91 | 0 / 0 | 604 | 3204 | 0.1579 | yes |
| 2^8 | 15.5 | 1 | 44 | 132 | 76 | 321 | 4.22 | 0.89 / 0.90 | 0 / 0 | 117 | 1522 | 0.2109 | yes |
| 2^8 | 15.5 | 1 | 48 | 144 | 110 | 383 | 3.48 | 0.90 / 0.95 | 0 / 0 | 116 | 1579 | 0.2426 | yes |
| 2^9 | 15.7 | 4 | 81 | 243 | 225 | 527 | 2.34 | 0.88 / 0.91 | 0 / 0 | 332 | 2760 | 0.1909 | yes |
| 2^10 | 16.2 | 4 | 84 | 252 | 152 | 508 | 3.34 | 0.93 / 0.92 | 0 / 0 | 630 | 4691 | 0.1083 | yes |
| 2^10 | 16.5 | 4 | 104 | 312 | 134 | 629 | 4.69 | 0.96 / 0.93 | 0 / 0 | 369 | 3965 | 0.1586 | yes |
| 2^9 | 17.1 | 1 | 58 | 174 | 67 | 451 | 6.73 | 0.97 / 0.95 | 0 / 0 | 193 | 3319 | 0.1359 | yes |
| 2^10 | 17.2 | 4 | 133 | 399 | 276 | 676 | 2.45 | 0.93 / 0.92 | 0 / 0 | 242 | 3889 | 0.1738 | yes |
| 2^11 | 17.9 | 7 | 222 | 666 | 462 | 1430 | 3.10 | 0.91 / 0.94 | 0 / 0 | 782 | 8956 | 0.1597 | yes |
| 2^11 | 19.0 | 4 | 264 | 792 | 612 | 2197 | 3.59 | 0.92 / 0.94 | 0 / 0 | 642 | 11860 | 0.1852 | yes |
| 2^10 | 19.2 | 1 | 141 | 423 | 196 | 1091 | 5.57 | 0.95 / 0.96 | 0 / 0 | 146 | 5818 | 0.1875 | yes |
| 2^11 | 21.3 | 1 | 264 | 792 | 545 | 2569 | 4.71 | 0.95 / 0.93 | 0 / 0 | 314 | 17765 | 0.1446 | yes |
| 2^11 | 21.7 | 1 | 307 | 921 | 776 | 2370 | 3.05 | 0.94 / 0.93 | 0 / 0 | 122 | 14351 | 0.1651 | yes |

**Every phase priced (rows at full rank on both arms)**

| p | log2 r | inv / mul (measured) | muls per addition | base build π / neg | pair table | stream π / neg | LA π / neg | total π | total negation | S π | S negation | S neg / S π | rho S negation | rho S folded | rho muls per group op, negation / folded | rho steps ratio (expected) | S π / rho S folded | all verified |
|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|:--|
| 2^7 | 10.7 | — | — | not full rank on both arms: the stream exhausted the group's targets | | | | | | | | | | | | | | |
| 2^8 | 11.8 | 2.4 | 83.9 | 3.37e+05 / 3.36e+05 | 5.09e+05 | 4.09e+06 / 1.66e+07 | 1.63e+03 / 3.42e+03 | 4.94e+06 | 1.75e+07 | 83240.0 | 294518.5 | 3.54 | 221.1 | 202.4 | 89 / 104 | 1.58 (1.73) | 411 | yes |
| 2^7 | 11.9 | — | — | not full rank on both arms: the stream exhausted the group's targets | | | | | | | | | | | | | | |
| 2^7 | 12.1 | 2.4 | 85.7 | 1.21e+05 / 1.2e+05 | 1.48e+05 | 4.58e+05 / 3e+06 | 786 / 1.52e+03 | 7.28e+05 | 3.27e+06 | 10945.5 | 49116.3 | 4.49 | 189.6 | 202.5 | 90 / 106 | 1.20 (1.73) | 54 | yes |
| 2^8 | 12.4 | 2.6 | 83.5 | 4.52e+05 / 4.51e+05 | 8.21e+05 | 3.51e+06 / 2.54e+07 | 2.07e+03 / 6.54e+03 | 4.78e+06 | 2.67e+07 | 65028.0 | 362761.0 | 5.58 | 200.4 | 168.5 | 89 / 107 | 1.87 (1.73) | 386 | yes |
| 2^7 | 13.4 | 2.8 | 84.9 | 1.86e+05 / 1.86e+05 | 2.45e+05 | 3.77e+06 / 8.44e+06 | 1.31e+03 / 1.87e+03 | 4.2e+06 | 8.87e+06 | 39768.1 | 84004.3 | 2.11 | 202.6 | 161.4 | 88 / 110 | 1.77 (1.73) | 246 | yes |
| 2^9 | 14.8 | 3.1 | 82.2 | 1.29e+06 / 1.29e+06 | 5.73e+06 | 1.7e+07 / 1.24e+08 | 7.79e+03 / 4.03e+04 | 2.41e+07 | 1.31e+08 | 140331.9 | 762621.8 | 5.43 | 145.9 | 104.1 | 88 / 109 | 2.26 (1.73) | 1348 | yes |
| 2^9 | 15.1 | 3.6 | 82.6 | 9.32e+05 / 9.3e+05 | 2.93e+06 | 1.82e+07 / 9.62e+07 | 5.14e+03 / 1.45e+04 | 2.2e+07 | 1e+08 | 118711.7 | 539163.0 | 4.54 | 173.3 | 145.0 | 90 / 106 | 1.93 (1.73) | 819 | yes |
| 2^8 | 15.5 | 2.6 | 83.3 | 5.5e+05 / 5.49e+05 | 1.46e+06 | 8.59e+06 / 3.37e+07 | 3.14e+03 / 8.58e+03 | 1.06e+07 | 3.57e+07 | 49135.3 | 165494.7 | 3.37 | 203.2 | 151.3 | 87 / 107 | 2.24 (1.73) | 325 | yes |
| 2^8 | 15.5 | 2.6 | 83.1 | 6.26e+05 / 6.24e+05 | 1.74e+06 | 9.82e+06 / 3.7e+07 | 3.62e+03 / 1.56e+04 | 1.22e+07 | 3.94e+07 | 56054.2 | 181145.9 | 3.23 | 167.7 | 153.2 | 88 / 107 | 1.63 (1.73) | 366 | yes |
| 2^9 | 15.7 | 3.4 | 82.3 | 1.2e+06 / 1.2e+06 | 4.99e+06 | 4.41e+07 / 1.06e+08 | 7.15e+03 / 2.79e+04 | 5.03e+07 | 1.12e+08 | 221091.8 | 494196.1 | 2.24 | 175.5 | 136.4 | 88 / 108 | 2.15 (1.73) | 1621 | yes |
| 2^10 | 16.2 | 3.5 | 82.2 | 1.29e+06 / 1.28e+06 | 5.37e+06 | 6.1e+07 / 1.99e+08 | 6.06e+03 / 4.51e+04 | 6.76e+07 | 2.05e+08 | 243220.5 | 738352.7 | 3.04 | 152.4 | 143.8 | 88 / 110 | 1.50 (1.73) | 1691 | yes |
| 2^10 | 16.5 | 3.6 | 82.1 | 1.65e+06 / 1.65e+06 | 8.24e+06 | 4.15e+07 / 1.97e+08 | 9e+03 / 4.99e+04 | 5.14e+07 | 2.07e+08 | 167825.5 | 673631.6 | 4.01 | 157.9 | 137.8 | 88 / 111 | 1.67 (1.73) | 1218 | yes |
| 2^9 | 17.1 | 3.0 | 82.8 | 7.77e+05 / 7.75e+05 | 2.55e+06 | 1.64e+07 / 1.01e+08 | 5.15e+03 / 1.75e+04 | 1.97e+07 | 1.05e+08 | 52125.2 | 276570.0 | 5.31 | 150.3 | 95.9 | 87 / 109 | 2.67 (1.73) | 544 | yes |
| 2^10 | 17.2 | 3.6 | 81.9 | 2.19e+06 / 2.19e+06 | 1.35e+07 | 9.7e+07 / 2.43e+08 | 9.9e+03 / 1.54e+05 | 1.13e+08 | 2.59e+08 | 290849.5 | 667415.4 | 2.29 | 112.8 | 116.7 | 89 / 111 | 1.29 (1.73) | 2492 | yes |
| 2^11 | 17.9 | 4.4 | 81.6 | 3.58e+06 / 3.57e+06 | 3.79e+07 | 2.93e+08 / 9.16e+08 | 3.6e+04 / 7.48e+05 | 3.35e+08 | 9.58e+08 | 688620.5 | 1970035.7 | 2.86 | 139.1 | 122.1 | 88 / 114 | 1.59 (1.73) | 5641 | yes |
| 2^11 | 19.0 | 4.0 | 81.5 | 5.07e+06 / 5.06e+06 | 5.33e+07 | 3.86e+08 / 1.42e+09 | 5.92e+04 / 8.85e+05 | 4.45e+08 | 1.48e+09 | 615720.3 | 2044461.1 | 3.32 | 107.5 | 91.1 | 87 / 115 | 1.69 (1.73) | 6759 | yes |
| 2^10 | 19.2 | 3.7 | 81.9 | 2.11e+06 / 2.1e+06 | 1.52e+07 | 6.61e+07 / 3.85e+08 | 1.9e+04 / 1.7e+05 | 8.34e+07 | 4.03e+08 | 107666.4 | 519523.6 | 4.83 | 124.4 | 59.5 | 87 / 112 | 3.57 (1.73) | 1809 | yes |
| 2^11 | 21.3 | 4.9 | 81.6 | 4.42e+06 / 4.41e+06 | 5.39e+07 | 4.48e+08 / 2.25e+09 | 5.98e+04 / 8.87e+05 | 5.06e+08 | 2.31e+09 | 318217.3 | 1451200.0 | 4.56 | 106.8 | 85.4 | 87 / 118 | 1.73 (1.73) | 3727 | yes |
| 2^11 | 21.7 | 5.2 | 81.5 | 5.09e+06 / 5.08e+06 | 7.31e+07 | 6.79e+08 / 2.07e+09 | 7.37e+04 / 9.61e+05 | 7.57e+08 | 2.15e+09 | 412185.9 | 1171277.2 | 2.84 | 110.0 | 97.4 | 87 / 119 | 1.53 (1.73) | 4232 | yes |

**Summary**

| arms | instances (of run) | log2 r | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | S neg / S π mean (min–max) | S π / rho S folded, min–max | rho steps ratio mean (expected) | rho muls per group op, negation / folded | rho S negation / rho S folded, mean | fitted exponent of total: π, negation, rho folded, rho negation (rho: 0.50) | all verified |
|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| π-line fold vs negation, m = 3, pair table | 18 (20) | 11.8–21.7 | 3.0 | 4.20 (2.34–6.73) | 0.92 / 0.92 | 3.75 (2.11–5.58) | 54–6759 | 1.88 (1.73) | 88 / 110 | 1.24 | 0.88, 0.86, 0.37, 0.39 | yes |


Reading.  With three summands the §6.5 degeneracy is gone: **zero** single-column
rows on either arm at every size, the fold's relation count falls by `4.20`
(`2.34`–`6.73`) against a column ratio of `3` — the coupon-collector
factor on top of `3`, as in E1 — and both arms sit at the same rank fraction after
`columns` relations (`0.92` against `0.92`).  The prediction of E8
holds; the two `2^7` instances with `r < 2^12` exhausted the group's targets before
full rank and carry no ratio.  The cost column is the first end-to-end `S` of this
note on a group where index calculus is asymptotically competitive (`n = 3`): with
the pair-table oracle every phase is counted in `F_p` multiplications — the base
build (square roots on the line), the pair table (`|F|²/4` additions, shared), the
stream (`|F|` additions a target), the linear algebra (row operations, one `Z/rZ`
multiplication each) — and the matched rho, negation and `π`-folded, is counted in
the same unit on the same instance, with inversions priced by the measured
`2.2`–`5.5` multiplications each.  The fold divides `S` by `3.75`
(`2.11`–`5.58`): the stream is `63`–`90 %`
of the folded arm's total and falls with the relations, the pair table
(`6`–`24 %`) does not fold (E12), and the linear
algebra is below `1 %` on both arms at these sizes.  The folded arm costs
`54`–`6759×` the `π`-folded rho.  Fitted exponents of the total
against `r` over the `18` instances at full rank: fold `0.88`, control
`0.86`, folded rho `0.37`, negation rho `0.39` — the
pair-table arms grow faster than rho's `1/2`, as the pair table (`p² ∝ r^{2/3}`) and
the stream (`p` targets at `p` additions) say they must; the algebraic oracle (E11)
is what would bring the relation phase to `p^{1+o(1)}`.  One number the unit
exposes that the step count hides: the `π`-folded walk takes `1.88×` fewer
steps (expected `√3 = 1.73`) but canonicalises on every step — the `⟨−1, π⟩`-orbit
by breadth-first search, three Frobenius maps and their keys — and that costs
`22` multiplications a group operation on top of the addition's
`88` (`110` against `88`), so in `F_p` multiplications the
negation rho costs `1.24×` the folded one on `E(F_{p³})`, not `1.73×`.  A first
version of the map that multiplied by the unit twist constants paid `230` a step
and made the folded walk the dearer of the two in this unit; that is why the
canonicalisation's price is reported beside the step ratio, and why the same
accounting applies to a folded base's canonicalisations in a relation search.
**Class: engineering** on the fold; the `S`
column is a stage-complete measurement of a method that is not faster than rho
here, and claims nothing else.

### 8.3 E9 — the matched folded rho on `F_{p²}` and `F_{p³}`

**Summary by family (per-instance table: the script, 60 rows)**

| family | A | instances | log2 r | walks | S ratio mean (min–max) | steps ratio mean | expected | all verified |
|:--|--:|--:|:--|--:|--:|--:|--:|:--|
| gls-generic | 4 | 12 | 9.3–22.9 | 192 | 1.28 (0.90–1.87) | 1.40 | 1.41 | True |
| gls-j0 | 12 | 12 | 11.0–23.0 | 192 | 1.71 (1.21–2.17) | 2.16 | 2.45 | True |
| gls-j1728 | 4 | 12 | 6.2–11.0 | 192 | 1.06 (0.93–1.28) | 1.15 | 1.41 | True |
| subfield | 6 | 12 | 11.8–23.8 | 192 | 1.48 (1.16–2.70) | 1.74 | 1.73 | True |
| subfield-j0 | 6 | 12 | 6.2–11.9 | 192 | 1.23 (0.96–1.64) | 1.72 | 1.73 | True |


Reading.  Every one of the `960` walks verifies its answer.  The steps ratios
track `√(A/2)` with the spread eight walks per instance give: `1.40` on the GLS
twist (`A = 4`, expected `1.41`), `1.74` on the subfield curve (`A = 6`,
expected `1.73`), and `2.16` on the `j = 0` GLS twist (`A = 12`, expected
`2.45`) — the four-dimensional GLV–GLS curves of Longa–Sica, the largest fold any
ordinary curve over `F_{p²}` admits, where the walks at `r ≤ 2^23` are too short
(`2^5`–`2^11` steps) to reach their asymptote.  The `j = 1728` twist and the `j = 0`
subfield curve are measured at `r ≈ p` (§6.3, §6.5) and show `1.15` and
`1.72` at `r ≤ 2^12`, where a walk is a few hundred steps and its set-up
dominates.  Every "vs rho" figure on these groups in §5 and §6 can now be read
against the matched walk instead of the `A = 2` one, with E8's caveat that the
matched walk's cost in field multiplications depends on what a canonicalisation
costs.  **Class: accounting** (reference).

### 8.5 E11 — the algebraic `S₄` oracle on the Frobenius line

**Relations, every instance**

| p | log2 r | h / #E(F_p) | cols π | cols negation | full-rank rel π | full-rank rel negation | ratio | rank fraction at k = cols, π / negation | single-column rows π / negation | orbit duplicates | targets | hit rate | correct |
|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|
| 2^7 | 10.7 | 4 | 11 | 33 | — | — | — | 0.75 / 0.74 | 0 / 0 | 1410 | 10180 | 0.0094 | yes |
| 2^8 | 11.8 | 7 | 26 | 78 | 50 | 181 | 3.62 | 0.85 / 0.90 | 0 / 0 | 598 | 1325 | 0.1366 | yes |
| 2^7 | 11.9 | 4 | 19 | 57 | — | — | — | 0.75 / 0.72 | 0 / 0 | 3110 | 31286 | 0.0084 | yes |
| 2^7 | 12.1 | 1 | 14 | 42 | 19 | 94 | 4.95 | 0.80 / 0.91 | 0 / 0 | 59 | 353 | 0.2663 | yes |
| 2^8 | 12.4 | 7 | 33 | 99 | 42 | 226 | 5.38 | 0.88 / 0.90 | 0 / 0 | 647 | 1637 | 0.1381 | yes |
| 2^7 | 13.4 | 1 | 18 | 54 | 46 | 116 | 2.52 | 0.95 / 0.89 | 0 / 0 | 110 | 755 | 0.1536 | yes |
| 2^9 | 14.8 | 7 | 87 | 261 | 100 | 672 | 6.72 | 0.98 / 0.91 | 0 / 0 | 658 | 3218 | 0.2088 | yes |
| 2^9 | 15.1 | 4 | 62 | 186 | 99 | 505 | 5.10 | 0.95 / 0.91 | 0 / 0 | 604 | 3204 | 0.1576 | yes |
| 2^8 | 15.5 | 1 | 44 | 132 | 75 | 317 | 4.23 | 0.91 / 0.91 | 0 / 0 | 117 | 1522 | 0.2083 | yes |
| 2^8 | 15.5 | 1 | 48 | 144 | 109 | 380 | 3.49 | 0.92 / 0.95 | 0 / 0 | 116 | 1579 | 0.2407 | yes |
| 2^9 | 15.7 | 4 | 81 | 243 | 225 | 527 | 2.34 | 0.88 / 0.91 | 0 / 0 | 332 | 2760 | 0.1909 | yes |
| 2^10 | 16.2 | 4 | 84 | 252 | 152 | 508 | 3.34 | 0.93 / 0.92 | 0 / 0 | 630 | 4691 | 0.1083 | yes |
| 2^10 | 16.5 | 4 | 104 | 312 | 133 | 626 | 4.71 | 0.96 / 0.93 | 0 / 0 | 369 | 3965 | 0.1579 | yes |
| 2^9 | 17.1 | 1 | 58 | 174 | 67 | 451 | 6.73 | 0.97 / 0.95 | 0 / 0 | 193 | 3319 | 0.1359 | yes |
| 2^10 | 17.2 | 4 | 133 | 399 | 274 | 671 | 2.45 | 0.93 / 0.92 | 0 / 0 | 242 | 3889 | 0.1725 | yes |
| 2^11 | 17.9 | 7 | 222 | 666 | 462 | 1429 | 3.09 | 0.91 / 0.94 | 0 / 0 | 782 | 8956 | 0.1596 | yes |
| 2^11 | 19.0 | 4 | 264 | 792 | 612 | 2193 | 3.58 | 0.92 / 0.94 | 0 / 0 | 642 | 11860 | 0.1849 | yes |
| 2^10 | 19.2 | 1 | 141 | 423 | 195 | 1086 | 5.57 | 0.95 / 0.96 | 0 / 0 | 146 | 5818 | 0.1867 | yes |
| 2^12 | 19.2 | 7 | 350 | 1050 | 768 | 2848 | 3.71 | 0.94 / 0.94 | 0 / 0 | 1503 | 20176 | 0.1412 | yes |
| 2^12 | 20.2 | 4 | 361 | 1083 | 1003 | 2892 | 2.88 | 0.94 / 0.94 | 0 / 0 | 877 | 20684 | 0.1398 | yes |
| 2^11 | 21.3 | 1 | 264 | 792 | 544 | 2564 | 4.71 | 0.95 / 0.93 | 0 / 0 | 314 | 17765 | 0.1443 | yes |
| 2^11 | 21.7 | 1 | 307 | 921 | 776 | 2368 | 3.05 | 0.94 / 0.93 | 0 / 0 | 122 | 14351 | 0.1650 | yes |
| 2^12 | 23.7 | 1 | 614 | 1842 | 1199 | 3979 | 3.32 | 0.93 / 0.93 | 0 / 0 | 113 | 26283 | 0.1514 | yes |
| 2^12 | 23.8 | 1 | 659 | 1977 | 1589 | 5156 | 3.24 | 0.95 / 0.93 | 0 / 0 | 184 | 33288 | 0.1549 | yes |
| 2^14 | 24.2 | 4 | 1409 | 4227 | 4086 | 12670 | 3.10 | 0.93 / 0.94 | 0 / 0 | 968 | 87047 | 0.1456 | yes |
| 2^13 | 25.0 | 1 | 977 | 2931 | 2349 | 9548 | 4.06 | 0.93 / 0.94 | 0 / 0 | 260 | 56958 | 0.1676 | yes |
| 2^13 | 25.0 | 1 | 956 | 2868 | 1875 | 7993 | 4.26 | 0.93 / 0.94 | 0 / 0 | 212 | 53004 | 0.1508 | yes |
| 2^14 | 27.0 | 1 | 1935 | 5805 | 5690 | 18731 | 3.29 | 0.94 / 0.94 | 0 / 0 | 258 | 123811 | 0.1513 | yes |

**Every phase priced (rows at full rank on both arms)**

| p | log2 r | inv / mul (measured) | muls per addition | base build π / neg | pair table | stream (incl. oracle) π / neg | LA π / neg | total π | total negation | S π | S negation | S neg / S π | rho S negation | rho S folded | rho muls per group op, negation / folded | rho steps ratio (expected) | S π / rho S folded | all verified |
|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|:--|
| 2^7 | 10.7 | — | — | not full rank on both arms: the stream exhausted the group's targets | | | | | | | | | | | | | | |
| 2^8 | 11.8 | 2.5 | 92.7 | 3.38e+05 / 3.37e+05 | 0 | 2.4e+08 / 9.74e+08 | 1.63e+03 / 3.42e+03 | 2.4e+08 | 9.74e+08 | 4046821.5 | 16430454.0 | 4.06 | 221.4 | 202.6 | 89 / 105 | 1.58 (1.73) | 19972 | yes |
| 2^7 | 11.9 | — | — | not full rank on both arms: the stream exhausted the group's targets | | | | | | | | | | | | | | |
| 2^7 | 12.1 | 2.2 | 91.4 | 1.21e+05 / 1.2e+05 | 0 | 4.48e+07 / 2.93e+08 | 786 / 1.52e+03 | 4.49e+07 | 2.93e+08 | 675457.9 | 4405395.8 | 6.52 | 189.0 | 202.0 | 89 / 106 | 1.20 (1.73) | 3344 | yes |
| 2^8 | 12.4 | 2.7 | 92.9 | 4.53e+05 / 4.51e+05 | 0 | 1.69e+08 / 1.23e+09 | 2.07e+03 / 6.54e+03 | 1.7e+08 | 1.23e+09 | 2307959.1 | 16678827.1 | 7.23 | 200.7 | 168.7 | 89 / 107 | 1.87 (1.73) | 13683 | yes |
| 2^7 | 13.4 | 2.4 | 93.0 | 1.86e+05 / 1.85e+05 | 0 | 2.84e+08 / 6.35e+08 | 1.31e+03 / 1.87e+03 | 2.84e+08 | 6.35e+08 | 2686978.2 | 6017586.0 | 2.24 | 201.6 | 160.8 | 87 / 110 | 1.77 (1.73) | 16714 | yes |
| 2^9 | 14.8 | 3.2 | 92.2 | 1.29e+06 / 1.29e+06 | 0 | 3.83e+08 / 2.78e+09 | 7.79e+03 / 4.03e+04 | 3.85e+08 | 2.79e+09 | 2244087.6 | 16254073.7 | 7.24 | 146.0 | 104.1 | 88 / 109 | 2.26 (1.73) | 21550 | yes |
| 2^9 | 15.1 | 2.8 | 92.9 | 9.27e+05 / 9.25e+05 | 0 | 5.15e+08 / 2.73e+09 | 5.13e+03 / 1.45e+04 | 5.16e+08 | 2.73e+09 | 2783147.6 | 14717631.2 | 5.29 | 171.7 | 143.9 | 89 / 105 | 1.93 (1.73) | 19338 | yes |
| 2^8 | 15.5 | 3.0 | 92.5 | 5.52e+05 / 5.5e+05 | 0 | 3.4e+08 / 1.33e+09 | 3.07e+03 / 8.97e+03 | 3.41e+08 | 1.33e+09 | 1579001.2 | 6186425.2 | 3.92 | 204.1 | 151.8 | 88 / 107 | 2.24 (1.73) | 10400 | yes |
| 2^8 | 15.5 | 2.6 | 92.1 | 6.26e+05 / 6.24e+05 | 0 | 3.74e+08 / 1.41e+09 | 3.73e+03 / 1.64e+04 | 3.74e+08 | 1.41e+09 | 1722359.5 | 6482738.5 | 3.76 | 167.7 | 153.1 | 88 / 107 | 1.63 (1.73) | 11247 | yes |
| 2^9 | 15.7 | 3.4 | 92.7 | 1.2e+06 / 1.2e+06 | 0 | 1.02e+09 / 2.46e+09 | 7.15e+03 / 2.79e+04 | 1.02e+09 | 2.46e+09 | 4495401.1 | 10819235.9 | 2.41 | 175.4 | 136.4 | 88 / 108 | 2.15 (1.73) | 32962 | yes |
| 2^10 | 16.2 | 3.3 | 93.7 | 1.28e+06 / 1.28e+06 | 0 | 1.25e+09 / 4.08e+09 | 6.06e+03 / 4.51e+04 | 1.25e+09 | 4.08e+09 | 4505118.5 | 14665701.5 | 3.26 | 152.1 | 143.5 | 88 / 110 | 1.50 (1.73) | 31393 | yes |
| 2^10 | 16.5 | 3.6 | 93.1 | 1.65e+06 / 1.65e+06 | 0 | 7.46e+08 / 3.53e+09 | 8.98e+03 / 5.19e+04 | 7.47e+08 | 3.53e+09 | 2437562.0 | 11513227.3 | 4.72 | 157.9 | 137.8 | 88 / 111 | 1.67 (1.73) | 17689 | yes |
| 2^9 | 17.1 | 3.5 | 93.4 | 7.79e+05 / 7.77e+05 | 0 | 4.77e+08 / 2.95e+09 | 5.15e+03 / 1.75e+04 | 4.77e+08 | 2.95e+09 | 1262964.5 | 7795194.7 | 6.17 | 151.2 | 96.3 | 87 / 109 | 2.67 (1.73) | 13111 | yes |
| 2^10 | 17.2 | 3.9 | 93.1 | 2.2e+06 / 2.19e+06 | 0 | 1.39e+09 / 3.48e+09 | 9.86e+03 / 1.16e+05 | 1.39e+09 | 3.49e+09 | 3598552.7 | 8997342.6 | 2.50 | 113.2 | 117.0 | 89 / 111 | 1.29 (1.73) | 30746 | yes |
| 2^11 | 17.9 | 4.3 | 93.1 | 3.58e+06 / 3.57e+06 | 0 | 2.58e+09 / 8.05e+09 | 3.6e+04 / 7.72e+05 | 2.58e+09 | 8.06e+09 | 5313399.7 | 16572286.9 | 3.12 | 138.9 | 122.0 | 88 / 114 | 1.59 (1.73) | 43562 | yes |
| 2^11 | 19.0 | 4.5 | 92.9 | 5.09e+06 / 5.08e+06 | 0 | 2.99e+09 / 1.1e+10 | 5.91e+04 / 8.88e+05 | 3e+09 | 1.1e+10 | 4151761.7 | 15217225.4 | 3.67 | 108.1 | 91.5 | 88 / 115 | 1.69 (1.73) | 45364 | yes |
| 2^10 | 19.2 | 4.0 | 93.0 | 2.11e+06 / 2.11e+06 | 0 | 9.06e+08 / 5.27e+09 | 1.86e+04 / 1.62e+05 | 9.08e+08 | 5.28e+09 | 1171590.0 | 6810054.5 | 5.81 | 124.8 | 59.6 | 87 / 112 | 3.57 (1.73) | 19641 | yes |
| 2^12 | 19.2 | 5.2 | 93.5 | 6.92e+06 / 6.91e+06 | 0 | 4.73e+09 / 1.83e+10 | 1e+05 / 1.86e+06 | 4.74e+09 | 1.83e+10 | 6034763.6 | 23337516.2 | 3.87 | 122.8 | 98.5 | 88 / 116 | 1.77 (1.73) | 61252 | yes |
| 2^12 | 20.2 | 5.6 | 93.5 | 6.99e+06 / 6.98e+06 | 0 | 6.53e+09 / 1.91e+10 | 1.36e+05 / 1.96e+06 | 6.54e+09 | 1.91e+10 | 5893919.5 | 17239831.3 | 2.93 | 105.5 | 102.3 | 88 / 118 | 1.40 (1.73) | 57626 | yes |
| 2^11 | 21.3 | 5.0 | 93.5 | 4.42e+06 / 4.41e+06 | 0 | 3.26e+09 / 1.64e+10 | 5.92e+04 / 8.84e+05 | 3.26e+09 | 1.64e+10 | 2052760.4 | 10299420.6 | 5.02 | 106.9 | 85.5 | 87 / 118 | 1.73 (1.73) | 24022 | yes |
| 2^11 | 21.7 | 5.0 | 93.4 | 5.09e+06 / 5.08e+06 | 0 | 4.35e+09 / 1.33e+10 | 7.37e+04 / 9.61e+05 | 4.36e+09 | 1.33e+10 | 2372929.0 | 7238675.9 | 3.05 | 109.8 | 97.3 | 87 / 119 | 1.53 (1.73) | 24398 | yes |
| 2^12 | 23.7 | 6.3 | 93.5 | 1.19e+07 / 1.19e+07 | 0 | 7.64e+09 / 2.5e+10 | 4.44e+05 / 1.08e+07 | 7.65e+09 | 2.5e+10 | 2056739.3 | 6721324.6 | 3.27 | 93.3 | 108.0 | 88 / 121 | 1.17 (1.73) | 19049 | yes |
| 2^12 | 23.8 | 6.7 | 93.5 | 1.29e+07 / 1.28e+07 | 0 | 9.59e+09 / 3.16e+10 | 5.73e+05 / 1.19e+07 | 9.6e+09 | 3.17e+10 | 2472468.6 | 8153526.3 | 3.30 | 79.9 | 93.8 | 89 / 121 | 1.15 (1.73) | 26350 | yes |
| 2^14 | 24.2 | 7.3 | 93.6 | 3.18e+07 / 3.18e+07 | 0 | 2.66e+10 / 8.32e+10 | 4.44e+06 / 1.2e+08 | 2.67e+10 | 8.34e+10 | 6125678.4 | 19151026.7 | 3.13 | 98.9 | 60.4 | 89 / 120 | 2.41 (1.73) | 101441 | yes |
| 2^13 | 25.0 | 7.2 | 93.4 | 2.05e+07 / 2.04e+07 | 0 | 1.26e+10 / 5.43e+10 | 1.3e+06 / 4.24e+07 | 1.26e+10 | 5.43e+10 | 2200628.7 | 9465212.7 | 4.30 | 111.6 | 72.5 | 89 / 122 | 2.19 (1.73) | 30373 | yes |
| 2^13 | 25.0 | 7.2 | 93.5 | 1.95e+07 / 1.95e+07 | 0 | 1.14e+10 / 5.05e+10 | 1.92e+06 / 3.65e+07 | 1.15e+10 | 5.05e+10 | 1984973.6 | 8743451.2 | 4.40 | 106.5 | 74.5 | 89 / 122 | 2.03 (1.73) | 26639 | yes |
| 2^14 | 27.0 | 7.8 | 93.6 | 4.22e+07 / 4.21e+07 | 0 | 3.68e+10 / 1.2e+11 | 1.38e+07 / 2.6e+08 | 3.69e+10 | 1.21e+11 | 3145232.1 | 10284379.7 | 3.27 | 126.3 | 80.3 | 89 / 124 | 2.19 (1.73) | 39150 | yes |

**Summary**

| arms | instances (of run) | log2 r | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | S neg / S π mean (min–max) | S π / rho S folded, min–max | rho steps ratio mean (expected) | rho muls per group op, negation / folded | rho S negation / rho S folded, mean | fitted exponent of total: π, negation, rho folded, rho negation (rho: 0.50) | all verified |
|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| π-line fold vs negation, m = 3, the algebraic S₄ oracle | 26 (28) | 11.8–27.0 | 3.0 | 3.98 (2.34–6.73) | 0.93 / 0.93 | 4.17 (2.24–7.24) | 3344–101441 | 1.85 (1.73) | 88 / 113 | 1.25 | 0.53, 0.50, 0.41, 0.42 | yes |

**The solver, per instance**

| p | log2 r | solver calls | F_p muls per call | Macaulay muls per call | unsolved (border unreachable) | retried at 11 / 12 / 13 | root triples | unliftable systems | agreement with the pair table (S₄ only / pair table only) | wall s per call |
|--:|--:|--:|--:|--:|:--|:--|--:|--:|:--|--:|
| 2^7 | 10.7 | 1692 | 830626 | 488503 | 30 (30) | 0 / 0 / 0 | 126 | 24 | 300/300 (0 / 0) | 0.0037 |
| 2^8 | 11.8 | 1118 | 867628 | 499124 | 7 (7) | 0 / 0 / 0 | 181 | 0 | 300/300 (0 / 0) | 0.0038 |
| 2^7 | 11.9 | 3732 | 856630 | 496050 | 42 (42) | 0 / 0 / 0 | 312 | 30 | 300/300 (0 / 0) | 0.0038 |
| 2^7 | 12.1 | 344 | 848029 | 483948 | 9 (9) | 0 / 0 / 0 | 94 | 0 | 300/300 (0 / 0) | 0.0031 |
| 2^8 | 12.4 | 1393 | 876293 | 501507 | 5 (5) | 0 / 0 / 0 | 237 | 8 | 300/300 (0 / 0) | 0.0037 |
| 2^7 | 13.4 | 733 | 862917 | 495715 | 6 (6) | 0 / 0 / 0 | 116 | 0 | 300/300 (0 / 0) | 0.0043 |
| 2^9 | 14.8 | 3041 | 911512 | 504912 | 9 (9) | 0 / 0 / 0 | 672 | 0 | 300/300 (0 / 0) | 0.0041 |
| 2^9 | 15.1 | 3043 | 892869 | 502908 | 22 (22) | 0 / 0 / 0 | 509 | 4 | 300/300 (0 / 0) | 0.0040 |
| 2^8 | 15.5 | 1488 | 892496 | 500635 | 13 (13) | 0 / 0 / 0 | 324 | 7 | 300/300 (0 / 0) | 0.0043 |
| 2^8 | 15.5 | 1563 | 897067 | 499185 | 21 (21) | 0 / 0 / 0 | 393 | 6 | 300/300 (0 / 0) | 0.0035 |
| 2^9 | 15.7 | 2711 | 902830 | 504796 | 8 (8) | 0 / 0 / 0 | 530 | 3 | 300/300 (0 / 0) | 0.0040 |
| 2^10 | 16.2 | 4566 | 888582 | 505412 | 12 (12) | 0 / 0 / 0 | 515 | 6 | — | 0.0038 |
| 2^10 | 16.5 | 3878 | 905095 | 505666 | 9 (9) | 0 / 0 / 0 | 632 | 5 | — | 0.0039 |
| 2^9 | 17.1 | 3290 | 890752 | 503764 | 16 (16) | 0 / 0 / 0 | 451 | 0 | 300/300 (0 / 0) | 0.0041 |
| 2^10 | 17.2 | 3846 | 901266 | 505857 | 11 (11) | 0 / 0 / 0 | 671 | 0 | — | 0.0038 |
| 2^11 | 17.9 | 8777 | 912673 | 506741 | 17 (17) | 0 / 0 / 0 | 1429 | 0 | — | 0.0040 |
| 2^11 | 19.0 | 11733 | 930989 | 506745 | 26 (26) | 0 / 0 / 0 | 2193 | 0 | — | 0.0040 |
| 2^10 | 19.2 | 5784 | 906446 | 505544 | 23 (23) | 0 / 0 / 0 | 1086 | 0 | — | 0.0040 |
| 2^12 | 19.2 | 19841 | 917668 | 507505 | 13 (13) | 0 / 0 / 0 | 2856 | 6 | — | 0.0039 |
| 2^12 | 20.2 | 20527 | 925360 | 507369 | 23 (23) | 0 / 0 / 0 | 2899 | 6 | — | 0.0040 |
| 2^11 | 21.3 | 17710 | 918622 | 507069 | 27 (27) | 0 / 0 / 0 | 2573 | 8 | — | 0.0038 |
| 2^11 | 21.7 | 14325 | 921489 | 507295 | 15 (15) | 0 / 0 / 0 | 2373 | 4 | — | 0.0040 |
| 2^12 | 23.7 | 26271 | 944232 | 507707 | 19 (19) | 0 / 0 / 0 | 3984 | 3 | — | 0.0041 |
| 2^12 | 23.8 | 33251 | 944961 | 507795 | 16 (16) | 0 / 0 / 0 | 5162 | 6 | — | 0.0041 |
| 2^14 | 24.2 | 86822 | 951714 | 508025 | 21 (21) | 0 / 0 / 0 | 12674 | 2 | — | 0.0041 |
| 2^13 | 25.0 | 56915 | 946086 | 507908 | 22 (22) | 0 / 0 / 0 | 9556 | 7 | — | 0.0041 |
| 2^13 | 25.0 | 52962 | 945661 | 507915 | 21 (21) | 0 / 0 / 0 | 7993 | 0 | — | 0.0040 |
| 2^14 | 27.0 | 123747 | 964494 | 508088 | 19 (19) | 0 / 0 / 0 | 18742 | 10 | — | 0.0041 |

**Solver summary**

| calls | F_p muls per call (mean) | unsolved | unsolved fraction | agreement with the pair table | S₄ only | pair table only |
|--:|--:|--:|--:|:--|--:|--:|
| 515103 | 942213 | 482 | 0.0009 | 3600/3600 | 0 | 0 |

**Against the pair table (E8), the same instances, `S` in the same unit**

| p | log2 r | S π, S₄ oracle | S π, pair table | pair / S₄ | S negation, S₄ | S negation, pair table | pair / S₄ |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 2^7 | 12.1 | 6.75e+05 | 1.09e+04 | 0.02 | 4.41e+06 | 4.91e+04 | 0.01 |
| 2^7 | 13.4 | 2.69e+06 | 3.98e+04 | 0.01 | 6.02e+06 | 8.4e+04 | 0.01 |
| 2^8 | 11.8 | 4.05e+06 | 8.32e+04 | 0.02 | 1.64e+07 | 2.95e+05 | 0.02 |
| 2^8 | 12.4 | 2.31e+06 | 6.5e+04 | 0.03 | 1.67e+07 | 3.63e+05 | 0.02 |
| 2^8 | 15.5 | 1.58e+06 | 4.91e+04 | 0.03 | 6.19e+06 | 1.65e+05 | 0.03 |
| 2^8 | 15.5 | 1.72e+06 | 5.61e+04 | 0.03 | 6.48e+06 | 1.81e+05 | 0.03 |
| 2^9 | 14.8 | 2.24e+06 | 1.4e+05 | 0.06 | 1.63e+07 | 7.63e+05 | 0.05 |
| 2^9 | 15.1 | 2.78e+06 | 1.19e+05 | 0.04 | 1.47e+07 | 5.39e+05 | 0.04 |
| 2^9 | 15.7 | 4.5e+06 | 2.21e+05 | 0.05 | 1.08e+07 | 4.94e+05 | 0.05 |
| 2^9 | 17.1 | 1.26e+06 | 5.21e+04 | 0.04 | 7.8e+06 | 2.77e+05 | 0.04 |
| 2^10 | 16.2 | 4.51e+06 | 2.43e+05 | 0.05 | 1.47e+07 | 7.38e+05 | 0.05 |
| 2^10 | 16.5 | 2.44e+06 | 1.68e+05 | 0.07 | 1.15e+07 | 6.74e+05 | 0.06 |
| 2^10 | 17.2 | 3.6e+06 | 2.91e+05 | 0.08 | 9e+06 | 6.67e+05 | 0.07 |
| 2^10 | 19.2 | 1.17e+06 | 1.08e+05 | 0.09 | 6.81e+06 | 5.2e+05 | 0.08 |
| 2^11 | 17.9 | 5.31e+06 | 6.89e+05 | 0.13 | 1.66e+07 | 1.97e+06 | 0.12 |
| 2^11 | 19.0 | 4.15e+06 | 6.16e+05 | 0.15 | 1.52e+07 | 2.04e+06 | 0.13 |
| 2^11 | 21.3 | 2.05e+06 | 3.18e+05 | 0.16 | 1.03e+07 | 1.45e+06 | 0.14 |
| 2^11 | 21.7 | 2.37e+06 | 4.12e+05 | 0.17 | 7.24e+06 | 1.17e+06 | 0.16 |

Reading.  The oracle is the one E11 asked for: `S₄(s t₁, s t₂, s t₃, x_R)` symmetrised
in the `t_i` (the fourth summation polynomial of `gaudry_cubic`, its coefficients
scaled by `s^{a + 2b + 3c}` so that the unknowns are the elementary symmetric
functions of the line parameters), Weil-descended to three `F_p` equations in three
unknowns and solved by the Macaulay matrix at degree `10`–`13`, the eigenvalues of
`e₁`, and the cubic split over `F_p`; the signs are settled by group arithmetic, so
nothing is returned that does not sum to the target.  It agrees with the pair table
on `3600` of `3600` targets at `p ≤ 2^9` (`0` decomposed by `S₄` alone,
`0` by the table alone), every decomposition it returns is checked, and the
relation stream it produces is E8's on the same instances — the same columns,
and the same trials and full-rank counts up to the few targets it leaves
unsolved — because an oracle changes what a target costs and not whether it
decomposes.  Its cost is independent of `p`:
`9.42e+05` `F_p` multiplications a call (`4.0` ms), `482` of
`515103` calls unsolved (no Macaulay degree closed the quotient; those targets are
skipped and reported, never guessed).  Against the pair table at `93`
multiplications an addition that is `10123` additions: the algebraic
oracle is the dearer one while `|F| < 10123` (`p ≲ 2^12`) and
the cheaper one above.  The comparison table holds the sizes both oracles were run
on, `p ≤ 2^11`, all below that crossover, and the pair table wins on every one of
them by the margin the per-target costs predict; above it E8 was not run — its
table at `p = 2^13` would hold `2^24` pairs — and E11 runs on at `2^12`, `2^13` and
`2^14` alone, where the per-target cost is the same `10⁶` multiplications as at
`2^7`.  The
fold's reading is unchanged — columns `÷ 3`, relations `÷ 3.98`, both arms at
`0.93` / `0.93` of full rank after `columns` relations, `S ÷ 4.17` — and
the fitted exponents of the total over `26` instances are fold `0.53`,
control `0.50`, rho `0.41` / `0.42`: with the oracle
flat in `p` the stream is `p · ln` solves at a constant each, so the arms grow as `r^{1/3}`
times the coupon factor, below rho's `1/2` for the first time in this note, and the
folded arm is still `3344`–`101441×` the matched rho at these sizes
(a constant of `10⁶` multiplications a solve against `10²` a step).  **Frobenius
symmetry breaking** (Galbraith–Granger–Merz–Petit's third lever) in the form the
driver can measure — one solve per `⟨−1, π⟩`-orbit of targets, the conjugate
systems' solutions being the conjugate decompositions — would save exactly the
orbit duplicates, `11056` of `512218` targets here (the birthday count of E7);
on a folded base a conjugate decomposition is the same row, so there is nothing
further to break at the system level, and the paper's own mechanism was not
available to this session to reproduce.  **Class: engineering** — the oracle moves
`S` and the exponent of the relation phase, not a ratio to a floor; the method is
not faster than rho at any size run.

### 8.4 Verdict after E8–E9, with E11

| lever | measured | class |
|:--|:--|:--|
| three summands on the `π`-line (E8) | columns `÷ 3`, relations `÷ 4.2`, `0` single-column rows; `S ÷ 3.8` end to end in one unit | engineering |
| the matched folded rho on `F_{p²}`, `F_{p³}` (E9) | steps `÷ 1.40, 1.74, 2.16` at `A = 4, 6, 12` (expected `1.41, 1.73, 2.45`), verified; in `F_p` multiplications the folded walk on `E(F_{p³})` pays `110` a step against `88` (E8) | accounting (reference) |
| the algebraic `S₄` line oracle (E11) | `9.42e+05` multiplications a target, flat in `p`; same relation stream as E8; exponent of the total `0.53` (fold) against the pair table's `0.88`; cheaper than the pair table from `p ≈ 2^12` | engineering |
| exponents (E10) | fold `0.88`, control `0.86`, rho `0.37` / `0.39` over `18` instances; no arm's exponent below rho's `1/2` | accounting |

With E11 run, the relation phase is `p^{1+o(1)}` solves on the folded line and the
next lever is E13 — the FGHR `2`-torsion symmetry of the system combined with the
fold — and E12 for the pair-table regime; nothing on the base side remains.
