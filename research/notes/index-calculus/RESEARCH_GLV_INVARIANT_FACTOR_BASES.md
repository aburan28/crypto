# Endomorphism-invariant factor bases across curve families: plan and pilot

**Modules:** `src/cryptanalysis/glv_invariant_base.rs` (the fold, prime-field automorphisms, Vélu degree-2 endomorphisms, CM instance generators), `src/cryptanalysis/gls_fp2.rs` (`F_{p²}` group, GLS `ψ`, the `ψ`-stable line), `src/cryptanalysis/ic_framework/plugins.rs` (`glv-orbit`, `gls-line`), `src/cryptanalysis/ic_boundary.rs` (`FactorBase::from_column_map`)
**CLI:** `ic bench --bits 20 --family j0 --factor-base glv-orbit:size=64 --oracle mitm:negation_folded=1` (control: `glv-orbit:size=64,no_fold=1`)
**Bench:** `cargo run --release --example glv_invariant_bench -- --families j0,j1728,generic,d7,d8 --bits 16,20,24 --seeds 2 --oracles subtract,mitm --json experiments/23_glv_invariant_pilot.json`; `--families gls --bits 8,10,12 --oracles subtract --json experiments/23_glv_invariant_gls_pilot.json`
**Data:** `experiments/23_glv_invariant_pilot.{json,log}`, `experiments/23_glv_invariant_gls_pilot.{json,log}` (2026-09-28, this host: Linux x86-64, 4 threads; wall time is not reported)
**Tables:** `python3 scripts/glv_invariant_tables.py experiments/23_glv_invariant_pilot.json experiments/23_glv_invariant_gls_pilot.json` (every number in §5 is printed by it from the frozen files)

> **Status.**  Implementation and pilot.  The fold is one function over
> the framework's `CountedGroup`, verified end to end on five curve
> families; the pilot is at `2^13`–`2^24` with two seeds per size and is
> **not** a measurement of the plan's experiments, whose protocols (§4)
> stay pending.  What the pilot establishes: the fold is exact (every
> row recovers its planted logarithm), the column count is the orbit
> count (`3×` on `j = 0`, `2×` on `j = 1728` and on GLS, `1×` on
> generic and on the degree-2 CM families), the relation count tracks
> the columns, and the degree-2 CM endomorphisms keep **zero** base
> points in the base at eigenvalue orders of `10³`–`10⁶`.  Class,
> where anything moved: **engineering** — the count moves with its
> floor, as `RESEARCH_GLV_INDEX_CALCULUS.md` §3 found for the `j = 0`
> quotient on `E(F_{p³})`.

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
| **E1** | Does the automorphism fold move anything but the column count on `F_p`? | `j0`, `j1728`, `generic` | `2^16`–`2^32`, six seeds | `glv-orbit` vs `no_fold`, oracles `subtract` and `mitm` | relations to **full rank** on one target stream (a second matrix fed the same rows until every column is pinned; the loop's first-determination count kept beside it), `S`, `S / rho` | `relations / columns ≈ 1.0` on both arms; `S` ratio `→ w/2` only where the relation phase dominates and `→ 1.0` where the pair table does (§5.2) | pilot only (§5) |
| **E2** | Does the GLS line fold match the Koblitz orbit fold in the same unit? | `gls` at `p = 2^8`–`2^16`; Koblitz `n = 13`–`31` from the ledger | six seeds | `gls-line` vs `no_fold`; `koblitz-orbit` vs `no_fold` | column ratio, relations, `S`, `S / rho`; the `descent-algebraic` analogue for `F_{p²}` (E2b) | `2×` on GLS against `n×` on Koblitz: the fold is the eigenvalue order, `4` against `2n`, and nothing else | pilot only (§5); E2b pending |
| **E2b** | An `O(1)` decomposition oracle for the line: Weil-descend `S₃(u s t₁, u s t₂, x_R)` to two `F_p` equations in `(t₁, t₂)`, resultant, roots by `gcd(t^p − t, ·)` | `gls` | as E2 | `subtract` vs the resultant oracle on the same base | oracle cost per target, hit rate agreement target by target (AGENTS.md §6 cross-check) | same hits, `O(p)` fewer group operations per target; the fold ratio unchanged | pending |
| **E3** | Does a composite group fold as its order says? | `j = 1728` twisted over `F_{p²}` (`⟨ι, ψ⟩`, order 8); `j = 0` twisted over `F_{p²}` (`⟨ψ₃, ψ⟩`, order 12) | `p = 2^8`–`2^14` | the line stable under both vs negation | points a column, verification, `S` | `8` and `12` points a column; engineering | pending; needs the twisted-CM instance generator (`generate_gls_instance` with `a = 0` or `b = 0` and the automorphism lifted to the twist) |
| **E4** | How much of a base does a type-C map keep in the base? | `d7`, `d8`, `1 + i` on `j1728` | `2^16`–`2^28` | one base per size, every rational degree-2 map | `ord_r(λ)`, `images_in_base / base_points` against `|F| / #E` | chance level; zero at pilot sizes | pilot only (§5.3) |
| **E4b** | The same for degree 3 (`D = −11`, `√−3` on `j = 0`) and for a Q-curve of degree 2 or 3 over `F_{p²}` | as named | as E4 | as E4 | chance level | pending; needs Vélu for an odd-order kernel and Smith's construction |
| **E5** | Subfield curves on `E(F_{p^n})`: a `π`-stable subspace not inside `F_p` | `E/F_p` on `E(F_{p³})`, with and without `j = 0` | `p = 2^6`–`2^11` | `⟨π⟩`, `⟨ψ⟩`, `⟨π, ψ⟩` folds vs negation on one base | points a column (`3`, `3`, `9`), relations, `S` | multiplicative; `glv_gaudry`'s `3.0×` is the `⟨ψ⟩` row | pending; needs an `F_{p³}` `CountedGroup` (the `gaudry_cubic` arithmetic wrapped) |
| **E6** | The matched reference | all folding families | as E1 | rho folded by the same `Aut` (counted) beside the `A = 2` walk | `S / rho_folded` | the fold's `w/2` in the relation count against rho's `√(w/2)` in steps: the gap widens by `√(w/2)` | pending; port `aut_folded_rho` to the counted `CountedGroup` walk |
| **E7** | Does the fold interact with a three-summand oracle? | `j0` | `2^18`–`2^30` | `mitm` at `m = 3` and the `S₄` oracle on folded vs control | relations, solver calls, rows that fold to `0 = 0` (the pair-generator finding of `RESEARCH_GLV_INDEX_CALCULUS.md` §4) | canonicalisation saves nothing on a uniform stream and everything on a pair sieve | pending |

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

## 6. Verdict, for the pilot

| lever | what moved | class | test against the boundary |
|:--|:--|:--|:--|
| automorphism fold, `F_p` (`j = 0`, `j = 1728`) | columns `÷ 3`, `÷ 2` on the same points; relations and trials with them, with first-determination scatter; `S ÷ 1.5–3.8` under `subtract`, `÷ 1.0–1.1` under `mitm` | engineering | count moves with its floor; `18×`–`2,834×` the `A = 2` rho |
| GLS line fold, `F_{p²}` | columns `÷ 2`; relations `÷ 1.2–2.0`; `S ÷ 1.3–2.0` | engineering | the eigenvalue order `4` against Koblitz's `2n`: the same fold, a smaller group |
| degree-2 CM maps | `0` base points kept in the base, `ord_r(λ) ∈ [3.5·10³, 3.4·10⁶]` | accounting | the first boundary of §2 holds at every size |
| the unified fold itself | one function for five families; the control is the same call | engineering | — |

## 7. What was not done

- No experiment of §4 was run beyond its pilot: E1 needs the
  full-rank stopping rule and six seeds at four or more sizes; E2b,
  E3, E4b, E5, E6 and E7 need the code named in their rows.
- No matched folded rho on prime fields (E6); every "vs rho" here is
  the `A = 2` walk and is marked so.
- The pilot's base build prices no square roots or Legendre symbols
  (no pinned ratio for generated curves); at these sizes the build is
  under `3 %` of `S` on every row, but a larger sweep should pin them
  as `ic bench` does for the roster.
- The scoreboard carries the pilot as a plan panel with its column
  and `S` ratios; no scoreboard row claims a speed, and the exponent
  panel is untouched.
- Nothing here bears on any deployed curve: `r ≤ 2^24`, certified toy
  instances, and a fold that rho already takes as `√(w/2)`.
