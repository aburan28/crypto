# Endomorphism-invariant factor bases across curve families: plan, pilot, and experiments E1–E7

**Modules:** `src/cryptanalysis/glv_invariant_base.rs` (the fold, prime-field automorphisms, Vélu degree-2 and degree-3 endomorphisms, CM instance generators, the folded rho classes), `src/cryptanalysis/ext_curve.rs` (`ExtField` over `F_{p²}` and `F_{p³}`, the generic `ExtCurve` counted group, diagonal automorphisms and Frobenius-type maps on it), `src/cryptanalysis/gls_fp2.rs` (GLS `ψ`, the `ψ`-stable line, the `j = 0` and `j = 1728` twists with their lifted automorphisms), `src/cryptanalysis/subfield_fp3.rs` (`E/F_p` on `E(F_{p³})`, the Frobenius eigenline), `src/cryptanalysis/line_oracle.rs` (the Weil-descent resultant oracle for a line, E2b), `src/cryptanalysis/glv_invariant_experiments.rs` (one relation stream feeding both arms to full rank), `src/cryptanalysis/ic_framework/plugins.rs` (`glv-orbit`, `gls-line`), `src/cryptanalysis/ic_boundary.rs` (`FactorBase::from_column_map`)
**CLI:** `ic bench --bits 20 --family j0 --factor-base glv-orbit:size=64 --oracle mitm:negation_folded=1` (control: `glv-orbit:size=64,no_fold=1`)
**Bench (pilot):** `cargo run --release --example glv_invariant_bench -- --families j0,j1728,generic,d7,d8 --bits 16,20,24 --seeds 2 --oracles subtract,mitm --json experiments/23_glv_invariant_pilot.json`; `--families gls --bits 8,10,12 --oracles subtract --json experiments/23_glv_invariant_gls_pilot.json`
**Runner (E1–E7):** `cargo run --release --example glv_invariant_experiments -- --exp e1 --bits 16,20,24,28 --seeds 6 --json experiments/23_glv_invariant_e1.json` (the exact command of every file is its `command` field)
**Data:** `experiments/23_glv_invariant_pilot.{json,log}`, `experiments/23_glv_invariant_gls_pilot.{json,log}` (2026-09-28), `experiments/23_glv_invariant_e{1,1_32,2,3,4,5,6,7}.{json,log}` (2026-09-29); this host: Linux x86-64, 4 threads; wall time is recorded and is not a result
**Tables:** `python3 scripts/glv_invariant_tables.py experiments/23_glv_invariant_pilot.json experiments/23_glv_invariant_gls_pilot.json` (§5) and `python3 scripts/glv_invariant_experiment_tables.py experiments/23_glv_invariant_e*.json` (§6); every number in §5 and §6 is printed by them from the frozen files

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
> ratio to a floor, and no scoreboard row claims a speed.

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
- **No three-summand oracle on the subfield line.**  §6.5 shows the
  two-summand decomposition on `E(F_{p³})` confined to one orbit by the
  `π`-fixed cofactor `E(F_p)`; the test of the fold's relation count on
  a subfield curve is the three-summand line oracle (Semaev `S₄` on the
  line, Weil-descended), which is not written.
- **`S` and `S / rho` are not re-measured at full rank.**  E1–E7 count
  relations, rank, trials, oracle cost and rho steps; the end-to-end
  `S` columns of §5.1 stay the pilot's, and no row of §6 claims a speed.
  Under the pair-table oracle both arms build the same table (§5.2), so
  a full-rank `S` would move as §5.1 says, but that is an extrapolation.
- **The line oracle is a library and example component only.**  It is
  not a framework plugin and `ic bench` cannot select it; the `subtract`
  and `mitm` oracles remain what the CLI offers on the line bases.
- **No matched folded rho on the GLS or subfield groups.**  E6 covers
  the prime-field automorphism groups; the `F_{p²}` and `F_{p³}` "vs
  rho" figures of §5 are still the `A = 2` walk and are marked so.
- The pilot's base build prices no square roots or Legendre symbols
  (no pinned ratio for generated curves); at these sizes the build is
  under `3 %` of `S` on every row, but a larger sweep should pin them
  as `ic bench` does for the roster.
- The scoreboard carries the pilot and E1–E7 as two panels with their
  column, rank and relation ratios; no scoreboard row claims a speed,
  and the exponent panel is untouched.
- Nothing here bears on any deployed curve: `r ≤ 2^32`, certified toy
  instances, and a fold that rho already takes as `√(w/2)`.
