# Endomorphism-invariant factor bases across curve families: plan, pilot, experiments E1–E15, and the road to the state of the art

**Modules:** `src/cryptanalysis/glv_invariant_base.rs` (the fold, prime-field automorphisms, Vélu degree-2 and degree-3 endomorphisms, CM instance generators, the folded rho classes), `src/cryptanalysis/ext_curve.rs` (`ExtField` over `F_{p²}` and `F_{p³}`, the generic `ExtCurve` counted group, diagonal automorphisms and Frobenius-type maps on it), `src/cryptanalysis/gls_fp2.rs` (GLS `ψ`, the `ψ`-stable line, the `j = 0` and `j = 1728` twists with their lifted automorphisms), `src/cryptanalysis/subfield_fp3.rs` (`E/F_p` on `E(F_{p³})`, the Frobenius eigenline), `src/cryptanalysis/line_oracle.rs` (the Weil-descent resultant oracle for a line, E2b), `src/cryptanalysis/glv_invariant_experiments.rs` (one relation stream feeding both arms to full rank), `src/cryptanalysis/ic_framework/plugins.rs` (`glv-orbit`, `gls-line`), `src/cryptanalysis/ic_boundary.rs` (`FactorBase::from_column_map`), `src/cryptanalysis/orbit_pair_table.rs` (the pair table over orbit representatives, E12), `src/cryptanalysis/fghr_line.rs` (the `Y`-line, the `τ_T` fold and the `D₃` conic-resultant oracle, E13), `src/cryptanalysis/q_curve.rs` (Q-curves of degree 2 and 3 over `F_{p²}` and `ψ = π ∘ ι ∘ φ`, E14)
**CLI:** `ic bench --bits 20 --family j0 --factor-base glv-orbit:size=64 --oracle mitm:negation_folded=1` (control: `glv-orbit:size=64,no_fold=1`)
**Bench (pilot):** `cargo run --release --example glv_invariant_bench -- --families j0,j1728,generic,d7,d8 --bits 16,20,24 --seeds 2 --oracles subtract,mitm --json experiments/23_glv_invariant_pilot.json`; `--families gls --bits 8,10,12 --oracles subtract --json experiments/23_glv_invariant_gls_pilot.json`
**Runner (E1–E7):** `cargo run --release --example glv_invariant_experiments -- --exp e1 --bits 16,20,24,28 --seeds 6 --json experiments/23_glv_invariant_e1.json` (the exact command of every file is its `command` field)
**Data:** `experiments/23_glv_invariant_pilot.{json,log}`, `experiments/23_glv_invariant_gls_pilot.{json,log}` (2026-09-28), `experiments/23_glv_invariant_e{1,1_32,2,3,4,5,6,7}.{json,log}` (2026-09-29), `experiments/23_glv_invariant_e{8,9,11,11_12,11_13,11_14}.{json,log}` (2026-10-01), `experiments/23_glv_invariant_e{12,12p,13,13_13,14,15}.{json,log}` (2026-10-02); this host: Linux x86-64, 4 threads; wall time is recorded and is not a result
**Tables:** `python3 scripts/glv_invariant_tables.py experiments/23_glv_invariant_pilot.json experiments/23_glv_invariant_gls_pilot.json` (§5; legacy Python, not extended) and `cargo run --release --example glv_invariant_experiment_tables -- experiments/23_glv_invariant_e*.json` (§6, §8; the native replacement of the retired `scripts/glv_invariant_experiment_tables.py`, byte-identical output on E1–E11); every number in §5, §6 and §8 is printed by them from the frozen files

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
> rho at any size run), then the pair table over orbit representatives
> (E12: the table `÷ w/2`, worth `S ÷ 1.18` where the table is the cost,
> `F_p` `j = 0`, and nothing on the line) and FGHR's `2`-torsion symmetry
> with the fold (E13: a degree-`16` resultant in place of a `64`-dimensional
> quotient, `S ÷ 27`, and the same `T` folding the base `12` a column for
> `S ÷ 2.5` more; the two compound, and the best arm is `267×` rho),
> and the transfer check to the ECC2K-130 family (E15: on `E_0` at the
> challenge's `4·prime` shape the fold is unconfined, `D = 1`; on
> `m = 31`, which cannot have that shape, the `373` cofactor confines
> `73 %` of two-summand rows), and Q-curves of degree 2 and 3 over
> `F_{p²}` (E14: every `ψ` verified, eigenvalue orders `≥ 1817`, base
> points kept at chance — no fold).  Every class is engineering or
> accounting; no row claims a speed.

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
`cargo run --release --example glv_invariant_experiment_tables -- experiments/23_glv_invariant_e*.json`
(native; it replaced the Python script these tables were first printed by, with identical output); every
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

**Summary by family, oracle and size (means over seeds; per-row table: the tables binary, 106 rows)**

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

**Summary by family and eigenvalue order (per-row table: the tables binary)**

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

**Summary (per-map table: the tables binary, 80 rows)**

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

**Summary by family (per-instance table: the tables binary)**

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

**Summary (per-instance table: the tables binary, 48 rows)**

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
- The scoreboard carries the pilot and E1–E13 as two panels with their
  column, rank and relation ratios; no scoreboard row claims a speed,
  and the exponent panel is untouched.
- Nothing here bears on any deployed curve: `r ≤ 2^32`, certified toy
  instances, and a fold that rho already takes as `√(w/2)`.
- Every experiment of §8.1 has run.  E14 builds its Q-curves by scanning `F_{p²}` for the `j` with `Φ_d(j, j^p) = 0` and for the kernel, so it stops at `p = 2^{12}`.  E15 ran on `E_0` at `m = 13, 15, 23, 31` only: the faithful `4·prime` degrees `19` and `41` are beyond the pair-table driver, and the §8a `m = 83` gate was not run (E15 claims no improvement).  E13 does not price the negation arm on the `Y`-line (it stalls below full rank there, not derived, §8.7), so the plan's fold-against-negation ratio is unmeasured on that base.

## 8. Toward the state of the art: what the literature does, what is next, and E8–E15

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
| **E12** | The pair table over orbit representatives: `P + φ^c P'` for representatives `P, P'` and `c < w/2` | table `w/2` smaller and `w/2` cheaper to build, one probe per target unchanged; `S` moves only by the table's share | a probe cost above `1.2×` the full table's falsifies "unchanged" | engineering (memory) | **done** (§8.6): table `÷ 2.95`–`2.97` (`w = 6`) and `÷ 1.98` (`w = 4`); probe `1.13×` at worst, so "unchanged" survives; `S ÷ 1.18` on `F_p` `j = 0`, `× 1.02` on the line |
| **E13** | FGHR's `2`-torsion symmetry **and** the fold together: close the line base under translation by `T ∈ E(F_p)[2]` and under `π`, fold by `⟨−1, π⟩`, symmetrise the system by `(Z/2)^{m−1} ⋊ S_m` | the two levers multiply — columns `÷ 3`, system degree `÷ 2^{m−1}` — because one acts on the base and the other on the system | a combined `S` ratio below the product of the separate ratios by more than `20 %` falsifies "multiply" | engineering | **done** (§8.7), with the base closed under `τ_T` and folded by it too (`12` a column): `R(q₃)` of degree `16 = 64/4`, `S ÷ 27.0` from the system, `÷ 2.52` from `τ_T`; combined / product `≥ 0.82`, so "multiply" survives; the negation arm stalls below full rank on the `Y`-line and is not priced |
| **E14** | Q-curves of degree 2 and 3 over `F_{p²}` (the second half of E4b) | type C: `0` base points kept, no fold | any image in the base beyond chance falsifies | accounting | **done** (§8.9): `40` Q-curves (`d = 2, 3`, `r = 2^{12}`–`2^{22}`), every `ψ` verified with `λ² ≡ ±d`, `ord_r(λ) ≥ 1817`; `32` of `43310` base points kept against `39.3` at chance — no fold |
| **E15** | Transfer to the binary Koblitz program (ECC2K-130, AGENTS.md §8a–8b): is the §6.5 degeneracy present there? | no: `E(F_2) ⊂ E(F_{2^n})` has order `2` or `4`, so at most two component classes and no block structure; the fold there is the known `2n` | a measured deficiency `D > 2` on a prime-order-times-`4` Koblitz subgroup with two summands falsifies | accounting | **done** (§8.8): on `E_0` at the challenge's `4·prime` shape (`m = 23`) `D = 1` on both arms and `4 %` single-column rows, so the prediction holds; at `m = 31`, whose `E_0` cannot be `4·prime` (cofactor `4·373`), `D = 12 / 342` and `73 %` single-column rows — `m = 31` is not a faithful proxy for two-summand relation structure |

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

**Summary by family (per-instance table: the tables binary, 60 rows)**

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

### 8.6 E12 — the pair table over orbit representatives

**Runner:** `--exp e12 --bits 7,8,9,10,11 --seeds 4` (the Frobenius line, three summands) and
`--exp e12p --bits 16,20,24,28 --seeds 6` (`j = 0` and `j = 1728` on `F_p`, two summands).
**Data:** `experiments/23_glv_invariant_e12.{json,log}`, `experiments/23_glv_invariant_e12p.{json,log}`
(2026-10-02).  **Module:** `src/cryptanalysis/orbit_pair_table.rs`.

The pair table of E8 stores `P_i + P_j` for every pair of `⟨−1⟩`-representatives.  On a base
that is a union of `H`-orbits (`|H| = w`), `φ(P_i) + φ(P_j)` is the image of `P_i + P_j`, so one
entry per `H`-orbit of pairs suffices: the table is built from the pairs `(rep_c, P_j)` with the
first summand a column representative and `col_j ≥ c`, `|F|²/(2w)` additions in place of
`|F|²/4`.  It is keyed by an **exact orbit invariant** of the sum — on the subfield line the
minimal polynomial of `x` over `F_p` (`σ₁ = 3a`, `σ₂ = 3(a² − νbc)`,
`σ₃ = a³ + νb³ + ν²c³ − 3νabc`: `12` multiplications, closed form), on `j = 0` `x³`, on
`j = 1728` `x²` — so a probe computes the key of `R − P₃` (or of `R`), and on a hit walks the
orbit of the stored sum to find the `φ` that carries it to the probe, reads the images of the two
stored summands off the column map by their coefficients, and checks the result with one
addition.  Nothing is returned that does not sum to the target.  The negation table and the orbit
table run the same stream on the same instances, every phase counted; the orbit key's
multiplications are not field-counter operations, so they are counted separately
(`build_uncounted_muls`, `stream_uncounted_muls`) and added.  Target by target, both tables agree on
whether a target decomposes on `200` targets per instance at `p ≤ 2^9` (line) and `3000` at
`≤ 20` bits (`F_p`).

**The Frobenius line, every instance**

| family | p | log2 r | pair sums built, negation / orbit | ÷ | orbit entries | table cost, negation / orbit | ÷ | stream cost per target, negation / orbit | orbit / negation (falsifies "unchanged" above 1.20) | keys / hits / maps | key collisions / image mismatches | full-rank rel π, negation table / orbit table | S π, negation table / orbit table | S π ratio | S negation arm, negation table / orbit table | rho S folded | S π (orbit) / rho S folded | agreement (orbit only / negation only) | correct |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|:--|:--|--:|:--|--:|--:|:--|:--|
| subfield | 2^7 | 10.7 | 1122 / 396 | 2.83 | 363 | 9.12e+04 / 3.69e+04 | 2.48 | 3.48e+03 / 3.61e+03 | 1.037 | 106960 / 96 / 1152 | 0 / 0 | — / — | — / — | — | — / — | 220.48 | — | 200/200 (0 / 0) | yes |
| subfield | 2^7 | 11.9 | 3306 / 1140 | 2.90 | 1083 | 2.72e+05 / 1.07e+05 | 2.54 | 4.01e+03 / 4.17e+03 | 1.039 | 402956 / 264 / 3168 | 0 / 0 | — / — | — / — | — | — / — | 209.24 | — | 200/200 (0 / 0) | yes |
| subfield | 2^7 | 12.1 | 1806 / 630 | 2.87 | 570 | 1.47e+05 / 5.89e+04 | 2.50 | 8.46e+03 / 9.27e+03 | 1.096 | 22938 / 94 / 1128 | 0 / 0 | 19 / 19 | 10905.3 / 10232.9 | 1.066 | 48922.6 / 51899.5 | 201.81 | 51 | 200/200 (0 / 0) | yes |
| subfield | 2^7 | 13.4 | 2970 / 1026 | 2.89 | 969 | 2.44e+05 / 9.65e+04 | 2.53 | 1.11e+04 / 1.22e+04 | 1.101 | 69419 / 116 / 1392 | 0 / 0 | 46 / 46 | 39559.6 / 41741.6 | 0.948 | 83558.7 / 90181.4 | 160.68 | 260 | 200/200 (0 / 0) | yes |
| subfield | 2^8 | 11.8 | 6162 / 2106 | 2.93 | 2024 | 5.09e+05 / 1.99e+05 | 2.56 | 1.26e+04 / 1.39e+04 | 1.111 | 152271 / 181 / 2172 | 0 / 0 | 50 / 50 | 83311.0 / 85748.9 | 0.972 | 294774.6 / 320719.4 | 202.52 | 423 | 200/200 (0 / 0) | yes |
| subfield | 2^8 | 12.4 | 9900 / 3366 | 2.94 | 3264 | 8.22e+05 / 3.2e+05 | 2.57 | 1.55e+04 / 1.73e+04 | 1.115 | 242391 / 226 / 2712 | 0 / 0 | 42 / 42 | 65088.3 / 63766.4 | 1.021 | 363108.6 / 396199.8 | 168.64 | 378 | 200/200 (0 / 0) | yes |
| subfield | 2^8 | 15.5 | 17556 / 5940 | 2.96 | 5784 | 1.46e+06 / 5.65e+05 | 2.59 | 2.22e+04 / 2.48e+04 | 1.118 | 329238 / 321 / 3852 | 0 / 0 | 76 / 76 | 49176.4 / 49738.0 | 0.989 | 165635.6 / 179971.3 | 151.40 | 329 | 200/200 (0 / 0) | yes |
| subfield | 2^8 | 15.5 | 20880 / 7056 | 2.96 | 6897 | 1.76e+06 / 6.77e+05 | 2.59 | 2.37e+04 / 2.63e+04 | 1.112 | 329063 / 356 / 4272 | 0 / 0 | 110 / 110 | 56614.4 / 56747.8 | 0.998 | 182983.9 / 179423.2 | 154.42 | 367 | 200/200 (0 / 0) | yes |
| subfield | 2^9 | 14.8 | 68382 / 22968 | 2.98 | 22681 | 5.74e+06 / 2.2e+06 | 2.61 | 3.84e+04 / 4.34e+04 | 1.129 | 1322050 / 672 / 8064 | 0 / 0 | 100 / 100 | 140381.7 / 132564.2 | 1.059 | 762897.8 / 835297.2 | 104.12 | 1273 | 200/200 (0 / 0) | yes |
| subfield | 2^9 | 15.1 | 34782 / 11718 | 2.97 | 11514 | 2.91e+06 / 1.12e+06 | 2.60 | 2.98e+04 / 3.35e+04 | 1.125 | 986578 / 506 / 6072 | 0 / 0 | 99 / 99 | 117758.6 / 120234.5 | 0.979 | 534771.5 / 589269.4 | 144.07 | 835 | 200/200 (0 / 0) | yes |
| subfield | 2^9 | 15.7 | 59292 / 19926 | 2.98 | 19653 | 4.97e+06 / 1.91e+06 | 2.60 | 3.84e+04 / 4.33e+04 | 1.128 | 1123121 / 527 / 6324 | 0 / 0 | 225 / 225 | 220482.2 / 231731.0 | 0.951 | 492825.5 / 538892.6 | 136.13 | 1702 | 200/200 (0 / 0) | yes |
| subfield | 2^9 | 17.1 | 30450 / 10266 | 2.97 | 10080 | 2.54e+06 / 9.79e+05 | 2.59 | 3.04e+04 / 3.42e+04 | 1.122 | 1026086 / 451 / 5412 | 0 / 0 | 67 / 67 | 52009.2 / 53171.5 | 0.978 | 275945.3 / 304529.6 | 95.72 | 555 | 200/200 (0 / 0) | yes |
| subfield | 2^10 | 16.2 | 63756 / 21420 | 2.98 | 21150 | 5.39e+06 / 2.07e+06 | 2.61 | 4.25e+04 / 4.79e+04 | 1.127 | 2111799 / 508 / 6096 | 0 / 0 | 152 / 152 | 244136.0 / 260216.4 | 0.938 | 741149.8 / 820526.9 | 144.22 | 1804 | — | yes |
| subfield | 2^10 | 16.5 | 97656 / 32760 | 2.98 | 32439 | 8.25e+06 / 3.16e+06 | 2.61 | 4.96e+04 / 5.6e+04 | 1.130 | 2117592 / 629 / 7548 | 0 / 0 | 134 / 134 | 167920.6 / 168893.9 | 0.994 | 674017.5 / 740558.4 | 137.87 | 1225 | — | yes |
| subfield | 2^10 | 17.2 | 159600 / 53466 | 2.99 | 53016 | 1.35e+07 / 5.16e+06 | 2.61 | 6.24e+04 / 7.07e+04 | 1.132 | 2658720 / 676 / 8112 | 0 / 0 | 276 / 276 | 290891.1 / 302386.5 | 0.962 | 667511.2 / 728561.6 | 116.75 | 2590 | — | yes |
| subfield | 2^10 | 19.2 | 179352 / 60066 | 2.99 | 59595 | 1.51e+07 / 5.79e+06 | 2.61 | 6.6e+04 / 7.47e+04 | 1.131 | 4188561 / 1091 / 13092 | 0 / 0 | 196 / 196 | 107403.0 / 106509.6 | 1.008 | 518240.8 / 571223.9 | 59.39 | 1793 | — | yes |
| subfield | 2^11 | 17.9 | 444222 / 148518 | 2.99 | 147804 | 3.79e+07 / 1.45e+07 | 2.62 | 1.02e+05 / 1.16e+05 | 1.134 | 10221825 / 1430 / 17160 | 0 / 0 | 462 / 462 | 689078.9 / 721766.2 | 0.955 | 1971350.1 / 2175716.3 | 122.13 | 5910 | — | yes |
| subfield | 2^11 | 19.0 | 628056 / 209880 | 2.99 | 209034 | 5.37e+07 / 2.05e+07 | 2.62 | 1.2e+05 / 1.37e+05 | 1.134 | 15956526 / 2197 / 26364 | 0 / 0 | 612 / 612 | 620055.3 / 646360.0 | 0.959 | 2058898.2 / 2278312.0 | 91.58 | 7058 | — | yes |
| subfield | 2^11 | 21.3 | 628056 / 209880 | 2.99 | 209013 | 5.39e+07 / 2.05e+07 | 2.63 | 1.27e+05 / 1.43e+05 | 1.133 | 24977043 / 2569 / 30828 | 0 / 0 | 545 / 545 | 318265.2 / 334848.4 | 0.950 | 1451418.7 / 1619098.4 | 85.40 | 3921 | — | yes |
| subfield | 2^11 | 21.7 | 849162 / 283668 | 2.99 | 282675 | 7.28e+07 / 2.77e+07 | 2.63 | 1.44e+05 / 1.63e+05 | 1.134 | 23042681 / 2370 / 28440 | 0 / 0 | 776 / 776 | 410423.2 / 435247.6 | 0.943 | 1166259.4 / 1292431.7 | 97.09 | 4483 | — | yes |

**Summary, the Frobenius line**

| family | instances (full rank on both tables) | log2 r | w | pair sums ÷ mean (min–max) | table cost ÷ mean (min–max) | stream cost per target, orbit / negation, mean (max) | table share of total π, negation table, mean (max) | S π ratio, negation / orbit table, mean (min–max) | S π (orbit) / rho S folded, min–max | key collisions | image mismatches | agreement (orbit only / negation only) | all correct |
|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| subfield | 18 (20) | 10.7–21.7 | 6 | 2.95 (2.83–2.99) | 2.59 (2.48–2.63) | 1.114 (1.134) | 0.133 (0.238) | 0.982 (0.938–1.066) | 51–7058 | 0 | 0 | 2400/2400 (0 / 0) | yes |

**Summary, `F_p` (per-instance table: the tables binary, 48 rows)**

| family | instances (full rank on both tables) | log2 r | w | pair sums ÷ mean (min–max) | table cost ÷ mean (min–max) | stream cost per target, orbit / negation, mean (max) | table share of total π, negation table, mean (max) | S π ratio, negation / orbit table, mean (min–max) | S π (orbit) / rho S folded, min–max | key collisions | image mismatches | agreement (orbit only / negation only) | all correct |
|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| j0 | 24 (24) | 12.9–26.0 | 6 | 2.97 (2.90–3.00) | 2.77 (2.67–2.83) | 1.005 (1.017) | 0.263 (0.637) | 1.180 (1.010–1.577) | 50–3971 | 0 | 0 | 36000/36000 (0 / 0) | yes |
| j1728 | 24 (24) | 12.9–26.9 | 4 | 1.98 (1.95–2.00) | 1.91 (1.87–1.95) | 1.002 (1.007) | 0.066 (0.137) | 1.031 (1.007–1.068) | 109–6330 | 0 | 0 | 36000/36000 (0 / 0) | yes |

Reading.  **The table is `w/2` smaller**, as predicted: pair sums `÷ 2.95` on the line and
`÷ 2.97` on `j = 0` (`w = 6`), `÷ 1.98` on `j = 1728` (`w = 4`), approaching `w/2` from below as
the diagonal pairs (stored once by both tables) become a vanishing share; the distinct keys are a
few per cent fewer still.  **It is cheaper to build by a little less than that** — `÷ 2.59`,
`÷ 2.77`, `÷ 1.91` — because every entry pays its key: `12` multiplications against an
addition's `84`–`93` on the line, `2` at `0.04` additions each on `F_p`.  **The probe is not quite
unchanged.**  Per target the stream costs `1.114×` the negation table's on the line (`1.134×` at
worst) and `1.005×` / `1.002×` on `F_p`: the key on every probe, and the orbit walk and image
check on every hit.  The falsification target was a probe above `1.2×`; no instance reaches it,
so "unchanged" survives, at `+11–13 %` on the line.  Zero key collisions and zero image mismatches
on `9.77e+07` probes; the two tables agree on `74400` of `74400` targets.

**`S` moves by the table's share, in both directions.**  On `j = 0` the table is `26 %` of the
folded arm's total on average (`64 %` at most), and the orbit table lowers `S` by `1.18×` on mean
(`1.01`–`1.58`).  Part of that spread is not the table: where a target has more than one
decomposition the two tables can return different ones, the relation streams then differ (`139`
against `115`, `240` against `422` relations to full rank on two rows), and so does the rank they
reach per relation; on the rows where the counts are equal the gain is `1.04`–`1.34`, which is
what a `÷ 2.8` table at those rows' `6`–`40 %` shares predicts (`1 / (1 − t(1 − 1/2.8))`).  On `j = 1728` the table is `7 %` of the total
and `S` falls by `1.03×`.  On the Frobenius line the table is `13 %` and the probe's `+11 %` is
paid on a stream that is most of the rest, so the orbit table is **dearer** end to end: `S` ratio
`0.98` on mean (`0.94`–`1.07`).  The fold's own `S` ratio and every comparison with rho are
unchanged: the orbit table is `50`–`7058×` the matched folded rho on these instances, as the full
table was.

**Class: engineering (memory).**  The table shrinks by `w/2` and on `F_p` `j = 0` `S` falls with it
by the table's share; on the subfield line the same change raises `S` by `2 %`, and quoting its
`÷ 3` table there as a cost gain would be **relabelling**.  No ratio to a floor moves.  What the
measurement does say is where the regime is: E12 matters exactly where the pair table is the
cost, which is two summands on `F_p` at large `p` — the regime §5.2 found worst against rho.

### 8.7 E13 — FGHR's `2`-torsion symmetry and the fold, together

**Runner:** `--exp e13 --bits 7,8,9,10,11,12 --seeds 4` and `--exp e13 --bits 13 --seeds 2`.
**Data:** `experiments/23_glv_invariant_e13.{json,log}`, `experiments/23_glv_invariant_e13_13.{json,log}`
(2026-10-02).  **Module:** `src/cryptanalysis/fghr_line.rs`.

**The object.**  A subfield curve `E/F_p` (E5, E8) with a rational `2`-torsion point
`T = (x₀, 0)` whose `f'(x₀) = 3x₀² + a = c²` is a square.  The Möbius coordinate
`Y = (x − x₀ − c)/(x − x₀ + c)` satisfies `Y(P + T) = −Y(P)` and `Y(−P) = Y(P)`, so the
**`Y`-line** `Y ∈ L·F_p` (`L = s` or `s²`, `s³ = ν`) is stable under `−1`, `π` and the
translation `τ_T`.  `T` lies in `E(F_p)`, the `π`-fixed cofactor, so `τ_T` acts on `⟨G⟩` with
eigenvalue `1`: `P` and `P + T` are one column with one coefficient, and the base folds by
`⟨−1, π, τ_T⟩` — **`12` points a column**, against `6` for `⟨−1, π⟩` on the same line and `2` for the
negation control.  That is the `2`-torsion used on the **base**.  On the **system**, `S₄`
rewritten in `Y` (by the substitution `x = ((c − x₀)Y + (x₀ + c))/(1 − Y)`) is invariant under the
even sign changes of `(Y₁, Y₂, Y₃)` — FGHR's `(Z/2)^{m−1}` at `m = 3` — and under their
permutations: `D₃ = (Z/2)² ⋊ S₃`, with invariants `p₁ = ΣY_i²`, `p₂ = ΣY_i²Y_j²`,
`p₃ = Y₁Y₂Y₃`.  On the line `p₁ = L²q₁`, `p₂ = L⁴q₂`, `p₃ = L³q₃` (weighted degree
`2i + 2j + k ≤ 4`); the Weil descent of `S₄(…, x_R)` is three **conics** in `(q₁, q₂)` with
coefficients in `F_p[q₃]`, their Macaulay resultant (the `15 × 15` determinant over its `3 × 3`
minor, interpolated at `17` points and checked at an `18`th) is a polynomial `R(q₃)`, its roots come
from Cantor–Zassenhaus, and each root gives `q₁, q₂` from the conics, `t_i²` from
`Z³ − q₁Z² + q₂Z − q₃²`, and the signs from `t₁t₂t₃ = q₃`; group arithmetic picks which of the
`τ_T`-variants of a solution to return (one, uniformly, per target).  The `S₃`-symmetrised
presentation of the same `S₄` in `Y` — E11's oracle moved to this line — is solved by the
existing Macaulay solver for comparison.  The driver runs four streams per instance — fold `12`
and fold `6`, each under `D₃` and (at `p ≤ 2^{10}`) under `S₃` — and every stream carries the
negation arm; the oracles are checked against the pair table target by target at `p ≤ 2^9`.

**The negation arm is not priced.**  On the four instances whose streams drew every target of the
group it ended at `40`/`48`, `12`/`24`, `100`/`114` and `169`/`174` of its columns, and in
development runs (not frozen) it stalled below full rank on larger instances too; whether that is a
structural ceiling of the `Y`-line negation base — a block degeneracy of the kind §6.5 found for two
summands under a `π`-fixed cofactor, now with `T` in that cofactor — is not derived here.  The
streams therefore stop when the folded arm is at full rank and both arms have pinned a column
(`StopRule::FoldedSquareBothPinned`), and every ratio below is between folded arms.  The
plan's base lever "columns `÷ 3`" (fold against negation) is consequently not measured; what is
measured is the base lever the `2`-torsion itself adds, `τ_T` on top of `⟨−1, π⟩`.

**Relations and solvers, every instance**

| p | log2 r | h | cols fold 12 / fold 6 / negation | full-rank rel fold 12 / fold 6 (D₃) | full-rank rel fold 12 / fold 6 (S₃) | negation rank at stop, D₃ streams (of cols) | D₃ calls fold 12 / fold 6 | D₃ muls per call | R(q₃) degree max | D₃ unsolved | S₃ muls per call | S₃ quotient per call | S₃ / D₃ muls per call | agreement with the pair table: D₃ / S₃ disagreements (targets, pair-table hits) | correct |
|--:|--:|--:|:--|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|:--|:--|
| 2^7 | 10.3 | 700 | 8 / 16 / 48 | — / — | — / — | 40 / 40 (48) | 1302 / 1302 | 13608 | 16 | 0 | 537898 | 64.0 | 39.5 | 0 / 0 (200, 14) | yes |
| 2^7 | 11.9 | 76 | 4 / 8 / 24 | — / — | — / — | 12 / 12 (24) | 3942 / 3942 | 12462 | 16 | 0 | 515868 | 63.9 | 41.4 | 0 / 0 (200, 2) | yes |
| 2^8 | 12.8 | 1708 | 19 / 38 / 114 | — / — | — / — | 100 / 100 (114) | 7026 / 7026 | 15405 | 16 | 0 | 572396 | 64.0 | 37.2 | 0 / 0 (200, 13) | yes |
| 2^8 | 12.9 | 1820 | 19 / 38 / 114 | 23 / 75 | 23 / 93 | 88 / 88 (114) | 2004 / 2004 | 15668 | 16 | 0 | 569442 | 64.0 | 36.3 | 0 / 0 (200, 14) | yes |
| 2^7 | 13.2 | 112 | 9 / 18 / 54 | 12 / 31 | 12 / 31 | 29 / 31 (54) | 520 / 557 | 13853 | 16 | 0 | 546165 | 64.0 | 39.4 | 0 / 0 (200, 12) | yes |
| 2^9 | 13.8 | 2464 | 29 / 58 / 174 | — / — | — / — | 169 / 169 (174) | 14712 / 14712 | 15824 | 16 | 0 | 571708 | 64.0 | 36.1 | 0 / 0 (200, 10) | yes |
| 2^7 | 14.2 | 112 | 10 / 20 / 60 | 13 / 43 | 13 / 43 | 45 / 45 (60) | 2505 / 2505 | 14969 | 16 | 0 | 543662 | 64.0 | 36.3 | 0 / 0 (200, 8) | yes |
| 2^8 | 15.4 | 184 | 19 / 38 / 114 | 36 / 46 | 36 / 133 | 82 / 82 (114) | 1280 / 1280 | 14779 | 16 | 0 | 571453 | 64.0 | 38.7 | 0 / 0 (200, 9) | yes |
| 2^8 | 15.5 | 244 | 21 / 42 / 126 | 26 / 68 | 26 / 58 | 93 / 93 (126) | 2020 / 2020 | 15532 | 16 | 0 | 567103 | 64.0 | 36.5 | 0 / 0 (200, 9) | yes |
| 2^9 | 17.4 | 404 | 39 / 78 / 234 | 45 / 144 | 45 / 188 | 214 / 214 (234) | 3394 / 3394 | 16302 | 16 | 0 | 586238 | 64.0 | 36.0 | 0 / 0 (200, 14) | yes |
| 2^9 | 17.8 | 448 | 37 / 74 / 222 | 51 / 231 | 51 / 187 | 193 / 199 (222) | 6077 / 7087 | 16529 | 16 | 0 | 582719 | 64.0 | 35.2 | 0 / 0 (200, 9) | yes |
| 2^9 | 17.9 | 464 | 42 / 84 / 252 | 59 / 156 | 59 / 204 | 223 / 223 (252) | 5321 / 5321 | 16772 | 16 | 0 | 588343 | 64.0 | 35.1 | 0 / 0 (200, 4) | yes |
| 2^11 | 18.3 | 10864 | 116 / 232 / 696 | 234 / 645 | — | 234 / 623 (696) | 6680 / 19688 | 19459 | 16 | 0 | — | — | — | — | yes |
| 2^10 | 19.1 | 776 | 66 / 132 / 396 | 187 / 254 | 187 / 593 | 357 / 357 (396) | 8238 / 8238 | 18252 | 16 | 0 | 597605 | 64.0 | 32.7 | — | yes |
| 2^10 | 19.1 | 700 | 58 / 116 / 348 | 87 / 232 | 87 / 408 | 318 / 318 (348) | 9793 / 9793 | 17725 | 16 | 0 | 590320 | 64.0 | 33.3 | — | yes |
| 2^10 | 19.2 | 688 | 59 / 118 / 354 | 73 / 218 | 73 / 318 | 271 / 271 (354) | 8575 / 8575 | 17118 | 16 | 0 | 593736 | 64.0 | 34.7 | — | yes |
| 2^10 | 19.4 | 892 | 72 / 144 / 432 | 131 / 251 | 119 / 319 | 396 / 396 (432) | 9669 / 9669 | 17362 | 16 | 0 | 593492 | 64.0 | 34.2 | — | yes |
| 2^12 | 19.9 | 17360 | 213 / 426 / 1278 | 651 / 821 | — | 1163 / 1163 (1278) | 29169 / 29169 | 20267 | 16 | 0 | — | — | — | — | yes |
| 2^11 | 20.6 | 1304 | 110 / 220 / 660 | 330 / 622 | — | 597 / 607 (660) | 13539 / 14078 | 19217 | 16 | 0 | — | — | — | — | yes |
| 2^11 | 20.6 | 1288 | 104 / 208 / 624 | 208 / 520 | — | 574 / 574 (624) | 15131 / 15131 | 19043 | 16 | 0 | — | — | — | — | yes |
| 2^11 | 21.4 | 1760 | 135 / 270 / 810 | 306 / 711 | — | 730 / 730 (810) | 20100 / 20100 | 18013 | 16 | 0 | — | — | — | — | yes |
| 2^12 | 22.7 | 2528 | 211 / 422 / 1266 | 367 / 1090 | — | 1162 / 1162 (1266) | 29638 / 29638 | 20242 | 16 | 0 | — | — | — | — | yes |
| 2^13 | 23.1 | 55804 | 639 / 1278 / 3834 | 1560 / 2981 | — | 2877 / 2981 (3834) | 74617 / 77969 | 20436 | 16 | 0 | — | — | — | — | yes |
| 2^12 | 23.6 | 3472 | 302 / 604 / 1812 | 542 / 1151 | — | 1663 / 1663 (1812) | 38421 / 38421 | 20048 | 16 | 0 | — | — | — | — | yes |
| 2^12 | 23.7 | 3760 | 318 / 636 / 1908 | 844 / 1640 | — | 1762 / 1762 (1908) | 41198 / 41198 | 19451 | 16 | 0 | — | — | — | — | yes |
| 2^13 | 24.7 | 5260 | 424 / 848 / 2544 | 834 / 2848 | — | 2340 / 2452 (2544) | 62824 / 76407 | 19825 | 16 | 0 | — | — | — | — | yes |

**Every phase priced (F_p multiplications; inversion = measured factor; each solver's own count; the polynomials' set-up charged to every arm; LA row op = 1); folded arms only**

| p | log2 r | S fold 12 + D₃ | S fold 6 + D₃ | S fold 12 + S₃ | S fold 6 + S₃ | base lever: S(6, D₃) / S(12, D₃) | system lever: S(12, S₃) / S(12, D₃) | product | combined: S(6, S₃) / S(12, D₃) | combined / product (falsifies "multiply" below 0.80) | rho S folded | S(12, D₃) / rho S folded |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2^7 | 10.3 | not full rank | not full rank | not full rank | not full rank | — | — | — | — | — | 258.3 | — |
| 2^7 | 11.9 | not full rank | not full rank | not full rank | not full rank | — | — | — | — | — | 227.9 | — |
| 2^8 | 12.8 | not full rank | not full rank | not full rank | not full rank | — | — | — | — | — | 184.1 | — |
| 2^8 | 12.9 | 96092.2 | 374345.6 | 2890711.4 | 13512186.9 | 3.90 | 30.08 | 117.19 | 140.62 | 1.20 | 163.8 | 587 |
| 2^7 | 13.2 | 44098.4 | 102127.9 | 1281290.8 | 3095734.2 | 2.32 | 29.06 | 67.29 | 70.20 | 1.04 | 134.9 | 327 |
| 2^9 | 13.8 | not full rank | not full rank | not full rank | not full rank | — | — | — | — | — | 136.6 | — |
| 2^7 | 14.2 | 96459.7 | 343441.2 | 2875207.7 | 9902415.0 | 3.56 | 29.81 | 106.13 | 102.66 | 0.97 | 130.0 | 742 |
| 2^8 | 15.4 | 43160.8 | 62595.6 | 1219000.2 | 6102299.8 | 1.45 | 28.24 | 40.96 | 141.39 | 3.45 | 158.4 | 272 |
| 2^8 | 15.5 | 61216.3 | 147282.1 | 1685340.1 | 3725555.1 | 2.41 | 27.53 | 66.24 | 60.86 | 0.92 | 135.6 | 451 |
| 2^9 | 17.4 | 37871.0 | 116163.8 | 988070.2 | 4395394.4 | 3.07 | 26.09 | 80.03 | 116.06 | 1.45 | 141.9 | 267 |
| 2^9 | 17.8 | 65373.2 | 323941.4 | 1744178.5 | 7095370.7 | 4.96 | 26.68 | 132.21 | 108.54 | 0.82 | 103.5 | 632 |
| 2^9 | 17.9 | 62216.2 | 147759.2 | 1643897.8 | 5837460.1 | 2.37 | 26.42 | 62.75 | 93.83 | 1.50 | 105.7 | 588 |
| 2^11 | 18.3 | 293056.3 | 850444.3 | — | — | 2.90 | — | — | — | — | 91.0 | 3220 |
| 2^10 | 19.1 | 134658.6 | 181276.9 | 3394852.8 | 11522881.7 | 1.35 | 25.21 | 33.94 | 85.57 | 2.52 | 88.8 | 1517 |
| 2^10 | 19.1 | 80293.8 | 226282.8 | 2021437.3 | 9904844.4 | 2.82 | 25.18 | 70.95 | 123.36 | 1.74 | 81.5 | 985 |
| 2^10 | 19.2 | 69454.1 | 206218.3 | 1808211.7 | 7733203.2 | 2.97 | 26.03 | 77.30 | 111.34 | 1.44 | 123.5 | 562 |
| 2^10 | 19.4 | 99983.7 | 171876.6 | 2375986.7 | 5548815.8 | 1.72 | 23.76 | 40.85 | 55.50 | 1.36 | 103.0 | 971 |
| 2^12 | 19.9 | 434231.3 | 546091.0 | — | — | 1.26 | — | — | — | — | 75.6 | 5741 |
| 2^11 | 20.6 | 153079.7 | 279830.9 | — | — | 1.83 | — | — | — | — | 87.2 | 1756 |
| 2^11 | 20.6 | 114411.2 | 268794.3 | — | — | 2.35 | — | — | — | — | 85.5 | 1338 |
| 2^11 | 21.4 | 125913.9 | 283597.2 | — | — | 2.25 | — | — | — | — | 64.1 | 1964 |
| 2^12 | 22.7 | 100636.5 | 288902.4 | — | — | 2.87 | — | — | — | — | 85.9 | 1172 |
| 2^13 | 23.1 | 379369.3 | 724818.7 | — | — | 1.91 | — | — | — | — | 79.6 | 4764 |
| 2^12 | 23.6 | 97272.7 | 194746.3 | — | — | 2.00 | — | — | — | — | 90.4 | 1076 |
| 2^12 | 23.7 | 138787.3 | 273503.3 | — | — | 1.97 | — | — | — | — | 77.1 | 1801 |
| 2^13 | 24.7 | 121355.3 | 397487.2 | — | — | 3.28 | — | — | — | — | 77.5 | 1565 |

**The phases of fold 12 + D₃, rows at full rank (share of the total)**

| p | log2 r | base build | polynomial set-up | stream group arithmetic | D₃ solver | linear algebra | total | dominant phase |
|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| 2^8 | 12.9 | 6e+05 (0.07) | 2.07e+05 (0.02) | 1.5e+06 (0.18) | 6.11e+06 (0.73) | 911 (0.00) | 8.42e+06 | solver |
| 2^7 | 13.2 | 1.91e+05 (0.04) | 2.07e+05 (0.05) | 8.05e+05 (0.19) | 3.15e+06 (0.72) | 298 (0.00) | 4.35e+06 | solver |
| 2^7 | 14.2 | 2.25e+05 (0.02) | 2.07e+05 (0.02) | 2.63e+06 (0.20) | 9.99e+06 (0.77) | 406 (0.00) | 1.31e+07 | solver |
| 2^8 | 15.4 | 4.76e+05 (0.05) | 2.07e+05 (0.02) | 1.81e+06 (0.20) | 6.44e+06 (0.72) | 840 (0.00) | 8.94e+06 | solver |
| 2^8 | 15.5 | 5.45e+05 (0.04) | 2.07e+05 (0.02) | 2.59e+06 (0.20) | 9.7e+06 (0.74) | 905 (0.00) | 1.3e+07 | solver |
| 2^9 | 17.4 | 1.05e+06 (0.07) | 2.07e+05 (0.01) | 3.2e+06 (0.21) | 1.11e+07 (0.71) | 2.53e+03 (0.00) | 1.56e+07 | solver |
| 2^9 | 17.8 | 1.06e+06 (0.03) | 2.07e+05 (0.01) | 6.65e+06 (0.22) | 2.29e+07 (0.74) | 1.9e+03 (0.00) | 3.08e+07 | solver |
| 2^9 | 17.9 | 1.2e+06 (0.04) | 2.07e+05 (0.01) | 6.65e+06 (0.21) | 2.3e+07 (0.74) | 2.55e+03 (0.00) | 3.1e+07 | solver |
| 2^11 | 18.3 | 4.61e+06 (0.03) | 2.07e+05 (0.00) | 3.36e+07 (0.20) | 1.3e+08 (0.77) | 5.69e+03 (0.00) | 1.68e+08 | solver |
| 2^10 | 19.1 | 2e+06 (0.02) | 2.07e+05 (0.00) | 2.18e+07 (0.22) | 7.67e+07 (0.76) | 5.02e+03 (0.00) | 1.01e+08 | solver |
| 2^10 | 19.1 | 1.92e+06 (0.03) | 2.07e+05 (0.00) | 1.32e+07 (0.22) | 4.49e+07 (0.75) | 4.32e+03 (0.00) | 6.02e+07 | solver |
| 2^10 | 19.2 | 1.79e+06 (0.03) | 1.69e+05 (0.00) | 1.19e+07 (0.22) | 3.93e+07 (0.74) | 4.06e+03 (0.00) | 5.32e+07 | solver |
| 2^10 | 19.4 | 2.56e+06 (0.03) | 2.07e+05 (0.00) | 1.89e+07 (0.23) | 6.18e+07 (0.74) | 4.72e+03 (0.00) | 8.34e+07 | solver |
| 2^12 | 19.9 | 9.04e+06 (0.02) | 2.07e+05 (0.00) | 9e+07 (0.21) | 3.27e+08 (0.77) | 3.67e+04 (0.00) | 4.26e+08 | solver |
| 2^11 | 20.6 | 3.72e+06 (0.02) | 2.07e+05 (0.00) | 4.34e+07 (0.22) | 1.47e+08 (0.76) | 7.84e+03 (0.00) | 1.94e+08 | solver |
| 2^11 | 20.6 | 3.46e+06 (0.02) | 2.07e+05 (0.00) | 3.26e+07 (0.22) | 1.1e+08 (0.75) | 1.16e+04 (0.00) | 1.46e+08 | solver |
| 2^11 | 21.4 | 4.77e+06 (0.02) | 2.07e+05 (0.00) | 5.09e+07 (0.24) | 1.54e+08 (0.73) | 1.48e+04 (0.00) | 2.1e+08 | solver |
| 2^12 | 22.7 | 7.95e+06 (0.03) | 2.07e+05 (0.00) | 5.98e+07 (0.23) | 1.91e+08 (0.74) | 2.39e+04 (0.00) | 2.59e+08 | solver |
| 2^13 | 23.1 | 3.22e+07 (0.03) | 2.07e+05 (0.00) | 2.66e+08 (0.24) | 8.23e+08 (0.73) | 4.6e+05 (0.00) | 1.12e+09 | solver |
| 2^12 | 23.6 | 1.13e+07 (0.03) | 2.07e+05 (0.00) | 8.4e+07 (0.24) | 2.52e+08 (0.73) | 7.44e+04 (0.00) | 3.48e+08 | solver |
| 2^12 | 23.7 | 1.22e+07 (0.02) | 2.07e+05 (0.00) | 1.29e+08 (0.25) | 3.75e+08 (0.73) | 9.94e+04 (0.00) | 5.16e+08 | solver |
| 2^13 | 24.7 | 1.69e+07 (0.03) | 2.07e+05 (0.00) | 1.62e+08 (0.26) | 4.54e+08 (0.72) | 1.23e+05 (0.00) | 6.33e+08 | solver |

Fitted exponent of the fold 12 + D₃ total against p (the line's field size, over the 22 rows above; the stream is p·ln p solves at a cost that grows with log p, so 1 + o(1) is expected): 1.23.

**Summary (rows where both D₃ arms reached full rank)**

| instances (of run) | log2 r | column ratio fold 12 / fold 6 | base lever S(6, D₃) / S(12, D₃), mean (min–max) | rows with S₃ on both bases | system lever S(12, S₃) / S(12, D₃), mean (min–max) | combined S(6, S₃) / S(12, D₃), mean (min–max) | combined / product, min | S(12, D₃) / rho S folded, min–max | fitted exponent of total against r: fold 12 + D₃, fold 6 + D₃, rho folded (rho: 0.50) | all correct |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| 22 (26) | 10.3–24.7 | 0.5 | 2.52 (1.26–4.96) | 12 | 27.01 (23.76–30.08) | 100.83 (55.50–141.39) | 0.82 | 267–5741 | 0.63, 0.60, 0.41 | yes |

Reading.  **The system lever is FGHR's, measured.**  `R(q₃)` has degree `16` — `64/4`, the
`S₃` quotient divided by `2^{m−1}` — on the calls that reach it (`866919` `D₃` calls, mean degree
`15.999`, maximum `16`, `0` failed degree checks, `0` unsolved, `0` unliftable), while the
`S₃` presentation of the same polynomial on the same line has a quotient of `63.9`–`64.0` a call.
A `D₃` call costs `1.25e+04`–`2.04e+04` multiplications, growing slowly with `p` (the root
finding is `O(log p)` multiplications per degree), against `5.2e+05`–`6.0e+05` for `S₃`
(`÷ 33`–`41`) and against E11's `9.42e+05` on the `x`-line (`÷ 46`–`76`).  Both oracles agree with the
pair table on every one of `2400` targets at `p ≤ 2^9` (`118` of them decomposable), and every
decomposition either returns is checked by group arithmetic.  End to end, on the fold-`12` base,
`S(12, S₃) / S(12, D₃) = 27.0` on mean (`23.8`–`30.1`): the solver is `71`–`77 %` of the
`D₃` arm's total on every row, so the per-call ratio arrives nearly whole.

**The base lever is the `2`-torsion's too.**  `τ_T` halves the columns exactly on every row and,
under `D₃`, `S(6, D₃) / S(12, D₃) = 2.52` on mean (`1.26`–`4.96`), the relations to full rank
falling `÷ 1.26`–`4.53` with them (the column ratio `2`, times the scatter of where full rank
lands).  On the fold-`12` base both oracles produce the same relation stream (a `τ_T`-variant of a
decomposition is the same row), and they do on `11` of `12` rows; on the fold-`6` base the
variants are different rows, and the oracle's choice among them matters: `D₃` picks one uniformly
per target, the Macaulay solver's choice is correlated across targets, and `S₃` needs
`0.81`–`2.89×` `D₃`'s relations to reach full rank there.  That is a property of solution
selection, not of either symmetry.

**"Multiply" survives its test.**  Combined, `S(6, S₃) / S(12, D₃) = 101` on mean (`55`–`141`);
the product of the separate ratios is `34`–`132`, and combined / product is `0.82`–`3.45`, never
below the `0.80` that would have falsified "multiply" (closest: `0.82` at `p = 2^9`).  The excess on
`9` of `12` rows is the selection effect above, which inflates the fold-`6` `S₃` arm; with
uniform selection it should sit at the product, which is not measured here.  So the two uses of `T` — a fold of the base and a
symmetry of the system — act on different objects and compound, as §8.0 predicted for the fold
and the system symmetries in general.

**Against rho, and the phases.**  `S(12, D₃)` is `3.8e+04`–`4.3e+05` against the matched folded
rho's `64`–`164` in the same unit: `267`–`5741×` rho, against E11's `3344`–`101441×` — one to
two orders of magnitude closer, by constants.  The solver dominates (`71`–`77 %`), the stream's
group arithmetic is `18`–`26 %`, the base build `2`–`7 %`, the polynomials' set-up `≤ 5 %`, the linear
algebra `< 0.1 %`.  The total grows as `p^{1.23}` (fitted over `p = 2^7`–`2^{13}`, `22` rows);
against `r` the fits are `0.63` (fold `12`), `0.60` (fold `6`) and `0.41` (rho), but `r ≈ p³/h`
with `h` from `76` to `55804` across these instances, so `r` is a poor size variable here and a
fit against it mixes the growth in `p` with the spread in `h`.  **Extrapolation, not measurement:** at fixed cofactor `p^{1.23}` is
`r^{0.41}`, so `S` falls as `r^{−0.09}` and a `10³` gap to rho would close only after `r` grows by
about `2^{110}` — on an exponent fitted over seven sizes of `p`.  Four small instances
(`r ≤ 2^{13.8}`) drew every target of the group with the folded arms at rank `columns` of
`columns + 1`, as E8 and E12 did at `p = 2^7`; they are reported, not priced.

**Class: engineering.**  `S` falls by `27×` (the system symmetry, against this line's `S₃`
oracle) and by `2.5×` (the `τ_T` fold) on top of E11's, and the ratio to rho falls by the same
constants; nothing moves an exponent (`p^{1.23}`, the same `p^{1+o(1)}` relation phase as E11),
and the method is not faster than rho at any size run.

### 8.8 E15 — does the §6.5 degeneracy reach the ECC2K-130 family?

**Runner:** `--exp e15 --bits 13,…,33 --seeds 4` (`m` from `13` to `33`; degrees with no
instance or no usable invariant subspace are logged and skipped).  **Data:**
`experiments/23_glv_invariant_e15.{json,log}` (2026-10-02).

**The question.**  On a subfield curve `E(F_p)` is cofactor and fixed by `π`, and two-summand
relations on the Frobenius line pair a point with its own orbit: every folded row is a
single-column row and the fold buys columns but not relations (§6.5).  The ECC2K-130 group has
the same shape in miniature — `E_0: y² + xy = x³ + 1` over `GF(2^{131})` has order `4·prime`, and
the cofactor `E_0(F_2) ≅ Z/4` is fixed by the Frobenius `τ`.  E15 runs E2's two-summand stream —
the signed-Frobenius-orbit base on a `τ`-invariant subspace against the abscissa control, the
pair-table oracle, both arms to full rank — on `E_0` itself at every degree the driver can carry,
and reads the structural deficiency `D` (§6.0) and the single-column rows.  Prediction (§8.1):
with `h = 4` there is at most one non-self-negative component class, so `D ≤ 1`, no block
structure, and the fold buys relations; `D > 2` on a `4·prime` subgroup falsifies.

**Which degrees are faithful.**  `#E_0(GF(2^m))` is `4·prime` at `m = 13, 19, 23, 41` in the range
`11`–`61` (and at `83` and `131`), and not at `31`: `#E_0(GF(2^{31})) = 4 · 373 · 1439393`, so every
prime-order subgroup there carries a cofactor of at least `4·373`.  Of the faithful degrees, `13`
(`ord₁₃(2) = 12`) and `19` (`18`) have no proper `τ`-invariant subspace — the base would be the
whole group (`13`, run and reported) or `2.6·10⁵` points (`19`, beyond the pair table) — and `41`
(`r ≈ 2^{39}`, a `20`-dimensional subspace) is beyond this driver.  That leaves **`m = 23`**:
`icv1-f2m23-t5197-69e76b73`, `r = 2095853`, `h = 4`, no intermediate subfield, `ord₂₃(2) = 11`
— so, as at `m = 31` and unlike `53`, `83` and `131`, the nontrivial cyclotomic block splits
(AGENTS.md §8b), which is what gives the `11`-dimensional invariant subspace the base is built on.
`m = 15` (subfields `GF(2^3)`, `GF(2^5)`) is reported separately and supports nothing at a prime
degree.

**Every instance**

| curve | m | ord_m(2) | intermediate subfields | log2 r | h | h = #E_0(F_2) | subspace dim | seed | cols fold / control | D fold / control | single-column rows fold / control (of relations) | chance, 1 / cols fold | full-rank rel fold / control | ratio | column ratio | rank fraction at k = cols, fold / control | correct |
|:--|--:|--:|:--|--:|:--|:--|--:|--:|:--|:--|:--|--:|:--|--:|--:|:--|:--|
| `icv1-f2m13-t181-515ee569` | 13 | 12 | none | 11.0 | 4 = 2·2 | yes | 12 | 1 | 155 / 2003 | 1 / 1 | 0 (0.000) / 0 (0.000) | 0.006 | — / — | — | 12.92 | 0.46 / — | yes |
| `icv1-f2m13-t181-515ee569` | 13 | 12 | none | 11.0 | 4 = 2·2 | yes | 12 | 2 | 155 / 2003 | 1 / 1 | 0 (0.000) / 0 (0.000) | 0.006 | — / — | — | 12.92 | 0.46 / — | yes |
| `icv1-f2m13-t181-515ee569` | 13 | 12 | none | 11.0 | 4 = 2·2 | yes | 12 | 3 | 155 / 2003 | 1 / 1 | 0 (0.000) / 0 (0.000) | 0.006 | — / — | — | 12.92 | 0.43 / — | yes |
| `icv1-f2m13-t181-515ee569` | 13 | 12 | none | 11.0 | 4 = 2·2 | yes | 12 | 4 | 155 / 2003 | 1 / 1 | 0 (0.000) / 0 (0.000) | 0.006 | — / — | — | 12.92 | 0.46 / — | yes |
| `icv1-f2m15-tm275-2d22ff5d` | 15 | 4 | GF(2^3), GF(2^5) | 9.6 | 44 = 2·2·11 | no | 6 | 1 | 2 / 16 | 1 / 6 | 14 (1.000) / 0 (0.000) | 0.500 | 2 / 14 | 7.00 | 8.00 | 1.00 / — | yes |
| `icv1-f2m15-tm275-2d22ff5d` | 15 | 4 | GF(2^3), GF(2^5) | 9.6 | 44 = 2·2·11 | no | 6 | 2 | 2 / 16 | 1 / 6 | 15 (1.000) / 0 (0.000) | 0.500 | 2 / 15 | 7.50 | 8.00 | 1.00 / — | yes |
| `icv1-f2m15-tm275-2d22ff5d` | 15 | 4 | GF(2^3), GF(2^5) | 9.6 | 44 = 2·2·11 | no | 6 | 3 | 2 / 16 | 1 / 6 | 16 (1.000) / 0 (0.000) | 0.500 | 2 / 16 | 8.00 | 8.00 | 1.00 / 1.00 | yes |
| `icv1-f2m15-tm275-2d22ff5d` | 15 | 4 | GF(2^3), GF(2^5) | 9.6 | 44 = 2·2·11 | no | 6 | 4 | 2 / 16 | 1 / 6 | 14 (1.000) / 0 (0.000) | 0.500 | 2 / 14 | 7.00 | 8.00 | 1.00 / — | yes |
| `icv1-f2m23-t5197-69e76b73` | 23 | 11 | none | 21.0 | 4 = 2·2 | yes | 11 | 1 | 45 / 1013 | 1 / 1 | 124 (0.037) / 2 (0.001) | 0.022 | 56 / 3397 | 60.66 | 22.51 | 0.91 / 0.84 | yes |
| `icv1-f2m23-t5197-69e76b73` | 23 | 11 | none | 21.0 | 4 = 2·2 | yes | 11 | 2 | 45 / 1013 | 1 / 1 | 144 (0.040) / 6 (0.002) | 0.022 | 97 / 3566 | 36.76 | 22.51 | 0.80 / 0.84 | yes |
| `icv1-f2m23-t5197-69e76b73` | 23 | 11 | none | 21.0 | 4 = 2·2 | yes | 11 | 3 | 45 / 1013 | 1 / 1 | 145 (0.043) / 4 (0.001) | 0.022 | 182 / 3364 | 18.48 | 22.51 | 0.84 / 0.83 | yes |
| `icv1-f2m23-t5197-69e76b73` | 23 | 11 | none | 21.0 | 4 = 2·2 | yes | 11 | 4 | 45 / 1013 | 1 / 1 | 132 (0.039) / 3 (0.001) | 0.022 | 111 / 3408 | 30.70 | 22.51 | 0.82 / 0.83 | yes |
| `icv1-f2m31-tm90707-c95f16f5` | 31 | 5 | none | 20.5 | 1492 = 2·2·373 | no | 10 | 1 | 20 / 590 | 12 / 342 | 1392 (0.729) / 15 (0.008) | 0.050 | 68 / 1910 | 28.09 | 29.50 | 0.67 / 0.57 | yes |
| `icv1-f2m31-tm90707-c95f16f5` | 31 | 5 | none | 20.5 | 1492 = 2·2·373 | no | 10 | 2 | 20 / 590 | 12 / 342 | 1406 (0.730) / 15 (0.008) | 0.050 | 23 / 1925 | 83.70 | 29.50 | 0.89 / 0.64 | yes |
| `icv1-f2m31-tm90707-c95f16f5` | 31 | 5 | none | 20.5 | 1492 = 2·2·373 | no | 10 | 3 | 20 / 590 | 12 / 342 | 1490 (0.737) / 17 (0.008) | 0.050 | 32 / 2023 | 63.22 | 29.50 | 0.89 / 0.56 | yes |
| `icv1-f2m31-tm90707-c95f16f5` | 31 | 5 | none | 20.5 | 1492 = 2·2·373 | no | 10 | 4 | 20 / 590 | 12 / 342 | 1354 (0.739) / 14 (0.008) | 0.050 | 51 / 1833 | 35.94 | 29.50 | 0.44 / 0.55 | yes |

**Summary by degree**

| curve | m | h | h = #E_0(F_2) | instances | D fold / control | single-column fraction, fold, mean (min–max) | against chance | full-rank ratio mean (min–max) | column ratio | all correct |
|:--|--:|:--|:--|--:|:--|--:|--:|--:|--:|:--|
| `icv1-f2m13-t181-515ee569` | 13 | 4 = 2·2 | yes | 4 | 1 / 1 | 0.000 (0.000–0.000) | 0.0× | — (—–—) | 12.92 | yes |
| `icv1-f2m15-tm275-2d22ff5d` | 15 | 44 = 2·2·11 | no | 4 | 1 / 6 | 1.000 (1.000–1.000) | 2.0× | 7.38 (7.00–8.00) | 8.00 | yes |
| `icv1-f2m23-t5197-69e76b73` | 23 | 4 = 2·2 | yes | 4 | 1 / 1 | 0.040 (0.037–0.043) | 1.8× | 36.65 (18.48–60.66) | 22.51 | yes |
| `icv1-f2m31-tm90707-c95f16f5` | 31 | 1492 = 2·2·373 | no | 4 | 12 / 342 | 0.734 (0.729–0.739) | 14.7× | 52.74 (28.09–83.70) | 29.50 | yes |

Reading.  **At the challenge's shape the prediction holds.**  On `icv1-f2m23-t5197-69e76b73`
`D = 1` on both arms on every seed, `4.0 %` of folded rows are single-column (`1.8×` the `1/45`
chance of two summands landing in one column, against `100 %` on the subfield line), both arms
reach the same rank fraction after `columns` relations (`0.80`–`0.91` / `0.83`–`0.84`), and the
fold's relation ratio, `18.5`–`60.7` (mean `36.7`), sits where the column ratio `22.5` and the
coupon factor put it.  The falsification target, `D > 2` on a `4·prime` subgroup, is not met: a
Frobenius-fixed cofactor of order `4` leaves two-summand relations unconfined, and the fold on
`E_0` buys relations as well as columns.  E2's `E_1` rows (`h = 2` at `m = 17, 23`) read the same
(`D = 1`).

**At `m = 31` it does not, and the cause is the cofactor, not the field.**
`icv1-f2m31-tm90707-c95f16f5` with `r = 1439393` has `h = 4·373`: `D = 12` / `342`, and `73 %` of
folded rows are single-column (`14.7×` chance) — the §6.5 confinement, from a cofactor that
`τ` does **not** fix (`E_0(F_2)` has order `4`; the `373`-part is not `F_2`-rational).  What the
two cases share with the subfield line is a cofactor that is large against the base (`1492`
against `1180` points here, `#E(F_p) ≈ p` against `≈ p` points there), so that a point's
cofactor component has few partners in the base; that reading is not derived here.  The fold
still reaches full rank at `m = 31` (`23`–`68` relations against `1833`–`2023`), because the
single-column rows it collects are not useless, but the relation structure is not the
challenge's.  E2's `m = 31` row was the same curve and subgroup and shows the same `D`.

**Consequence for the ECC2K-130 program.**  AGENTS.md §8b names `m = 31` the primary exploratory
size and asks that its different Frobenius-module structure be disclosed; E15 adds a second
disclosure.  `E_0` over `GF(2^{31})` cannot have the challenge's `4·prime` order, and any
two-summand relation statistic measured there is shaped by the `373` cofactor (here `73 %`
single-column rows against `4 %` at the faithful `m = 23`).  For two-summand relation structure
the faithful sizes this driver can carry stop at `m = 23`; `41`, `83` (the §8a gate) and `131` are
the faithful sizes above it.  Nothing here is an improvement claim, so the §8a `m = 83` gate is
neither required nor discharged.

**Class: accounting.**  No algorithm changed; E15 measures where E2's and §6.5's relation
structure applies, and finds that it follows the cofactor.

### 8.9 E14 — Q-curves of degree 2 and 3 over `F_{p²}`

**Runner:** `--exp e14 --bits 8,9,10,11,12 --seeds 4`.  **Data:**
`experiments/23_glv_invariant_e14.{json,log}` (2026-10-02).  **Module:**
`src/cryptanalysis/q_curve.rs`.

**The object.**  A degree-`d` Q-curve is an `E/F_{p²}` with a `d`-isogeny `φ: E → E^σ` to its
Frobenius conjugate; `ψ = π ∘ ι ∘ φ` is then an endomorphism of `E` over `F_{p²}` of degree
`dp`, with `ψ² = [±d]` on `E(F_{p²})` (Smith's construction of fast GLV curves from Q-curves).
It is the `F_{p²}` analogue of E4's degree-2 and degree-3 CM maps, and the second half of E4b.
The module finds the Q-curve `j`-invariants as the `j ∈ F_{p²} ∖ F_p` with
`Φ_d(j, j^p) = 0` (the repository's tabulated `Φ₂`, `Φ₃`; `Φ_d(j, j^p)` lies in `F_p`, one
equation in the two coordinates of `j`), counts `E_j` and its quadratic twist, finds the
kernel among the `F_{p²}`-roots of the `d`-division polynomial and **checks** it by Vélu's
codomain (`j(φ(E)) = j^p` and `μ⁴a' = a^σ`, independent of `Φ_d`), and reads `ψ`'s eigenvalue
off `ψ(G)` by baby-step giant-step.  Every map is then verified on random points of `⟨G⟩`
(eigenvalue and additivity, as E4's were) and `ψ²(G) = [±d]G` is checked in the unit tests.
The base is the `F_p`-line — every point with `x ∈ F_p`, folded by negation — the set a GLS
map (`d = 1`) would keep whole if it were its line; the prediction (§8.1) is the chance overlap,
`|F| / #E` of the base.

**Every instance**

| d | p | log2 r | h | twist | seed | map | λ² ≡ | ord_r(λ) | images in base | base points | chance fraction | expected at chance | verified |
|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|:--|
| 2 | 131 | 12.1 | 4 | yes | 1 | `q-curve-psi[d=2,lambda=3570]` | −d | 2148 | 2 | 118 | 0.00687 | 0.81 | True |
| 2 | 397 | 13.5 | 14 | yes | 1 | `q-curve-psi[d=2,lambda=6707]` | +d | 11310 | 0 | 400 | 0.00253 | 1.01 | True |
| 2 | 197 | 14.2 | 2 | no | 2 | `q-curve-psi[d=2,lambda=7459]` | +d | 9716 | 0 | 206 | 0.00530 | 1.09 | True |
| 2 | 281 | 14.3 | 4 | yes | 4 | `q-curve-psi[d=2,lambda=7425]` | +d | 6584 | 0 | 284 | 0.00359 | 1.02 | True |
| 2 | 199 | 14.3 | 2 | no | 3 | `q-curve-psi[d=2,lambda=7032]` | −d | 9945 | 2 | 186 | 0.00468 | 0.87 | True |
| 2 | 359 | 14.4 | 6 | no | 3 | `q-curve-psi[d=2,lambda=1913]` | −d | 4280 | 2 | 360 | 0.00280 | 1.01 | True |
| 2 | 229 | 14.7 | 2 | yes | 4 | `q-curve-psi[d=2,lambda=6142]` | +d | 8720 | 0 | 212 | 0.00405 | 0.86 | True |
| 2 | 691 | 14.9 | 16 | no | 3 | `q-curve-psi[d=2,lambda=17962]` | +d | 29878 | 0 | 674 | 0.00141 | 0.95 | True |
| 2 | 307 | 15.5 | 2 | yes | 2 | `q-curve-psi[d=2,lambda=44388]` | +d | 23571 | 0 | 294 | 0.00312 | 0.92 | True |
| 2 | 659 | 16.1 | 6 | yes | 1 | `q-curve-psi[d=2,lambda=71503]` | −d | 72160 | 0 | 600 | 0.00139 | 0.83 | True |
| 2 | 617 | 16.5 | 4 | no | 4 | `q-curve-psi[d=2,lambda=88319]` | +d | 95088 | 2 | 604 | 0.00159 | 0.96 | True |
| 2 | 661 | 17.7 | 2 | no | 2 | `q-curve-psi[d=2,lambda=186822]` | −d | 108924 | 0 | 670 | 0.00154 | 1.03 | True |
| 2 | 1259 | 18.0 | 6 | no | 1 | `q-curve-psi[d=2,lambda=262503]` | −d | 263760 | 2 | 1238 | 0.00078 | 0.97 | True |
| 2 | 1409 | 18.3 | 6 | no | 3 | `q-curve-psi[d=2,lambda=283022]` | −d | 36714 | 0 | 1390 | 0.00070 | 0.97 | True |
| 2 | 1381 | 18.9 | 4 | yes | 2 | `q-curve-psi[d=2,lambda=226919]` | +d | 238299 | 2 | 1356 | 0.00071 | 0.96 | True |
| 2 | 2243 | 20.3 | 4 | yes | 1 | `q-curve-psi[d=2,lambda=461391]` | −d | 314610 | 2 | 2238 | 0.00044 | 1.00 | True |
| 2 | 1621 | 20.3 | 2 | no | 4 | `q-curve-psi[d=2,lambda=679215]` | −d | 437680 | 0 | 1620 | 0.00062 | 1.00 | True |
| 2 | 2861 | 21.0 | 4 | yes | 2 | `q-curve-psi[d=2,lambda=657756]` | +d | 2046192 | 2 | 2844 | 0.00035 | 0.99 | True |
| 2 | 2269 | 21.3 | 2 | no | 3 | `q-curve-psi[d=2,lambda=2338418]` | −d | 2572032 | 0 | 2238 | 0.00044 | 0.97 | True |
| 2 | 3169 | 22.3 | 2 | yes | 4 | `q-curve-psi[d=2,lambda=3727207]` | +d | 78492 | 0 | 3132 | 0.00031 | 0.98 | True |
| 3 | 197 | 13.2 | 4 | no | 1 | `q-curve-psi[d=3,lambda=3243]` | −d | 9630 | 2 | 200 | 0.00519 | 1.04 | True |
| 3 | 337 | 13.8 | 8 | yes | 3 | `q-curve-psi[d=3,lambda=1592]` | +d | 7079 | 0 | 340 | 0.00300 | 1.02 | True |
| 3 | 241 | 13.8 | 4 | no | 2 | `q-curve-psi[d=3,lambda=3618]` | +d | 7296 | 0 | 252 | 0.00432 | 1.09 | True |
| 3 | 239 | 14.2 | 3 | yes | 3 | `q-curve-psi[d=3,lambda=12163]` | +d | 19078 | 0 | 238 | 0.00416 | 0.99 | True |
| 3 | 151 | 14.5 | 1 | no | 4 | `q-curve-psi[d=3,lambda=10323]` | +d | 22740 | 0 | 144 | 0.00633 | 0.91 | True |
| 3 | 461 | 16.1 | 3 | no | 2 | `q-curve-psi[d=3,lambda=70685]` | +d | 35573 | 2 | 460 | 0.00216 | 0.99 | True |
| 3 | 571 | 16.3 | 4 | no | 2 | `q-curve-psi[d=3,lambda=36214]` | +d | 81552 | 0 | 530 | 0.00162 | 0.86 | True |
| 3 | 641 | 17.1 | 3 | no | 4 | `q-curve-psi[d=3,lambda=39148]` | +d | 68669 | 2 | 664 | 0.00161 | 1.07 | True |
| 3 | 647 | 17.1 | 3 | no | 1 | `q-curve-psi[d=3,lambda=139319]` | +d | 139966 | 4 | 652 | 0.00155 | 1.01 | True |
| 3 | 397 | 17.3 | 1 | yes | 1 | `q-curve-psi[d=3,lambda=87823]` | +d | 79080 | 0 | 400 | 0.00253 | 1.01 | True |
| 3 | 439 | 17.6 | 1 | no | 4 | `q-curve-psi[d=3,lambda=128902]` | +d | 193572 | 0 | 442 | 0.00228 | 1.01 | True |
| 3 | 809 | 17.7 | 3 | yes | 3 | `q-curve-psi[d=3,lambda=178911]` | +d | 109289 | 2 | 758 | 0.00116 | 0.88 | True |
| 3 | 1709 | 17.9 | 12 | no | 1 | `q-curve-psi[d=3,lambda=191366]` | +d | 1817 | 4 | 1692 | 0.00058 | 0.98 | True |
| 3 | 1759 | 18.4 | 9 | no | 4 | `q-curve-psi[d=3,lambda=196394]` | −d | 171771 | 0 | 1786 | 0.00058 | 1.03 | True |
| 3 | 2663 | 18.8 | 16 | no | 4 | `q-curve-psi[d=3,lambda=349934]` | −d | 12310 | 0 | 2734 | 0.00039 | 1.05 | True |
| 3 | 1619 | 19.3 | 4 | yes | 2 | `q-curve-psi[d=3,lambda=36459]` | −d | 109242 | 0 | 1670 | 0.00064 | 1.06 | True |
| 3 | 2273 | 20.3 | 4 | yes | 3 | `q-curve-psi[d=3,lambda=615206]` | −d | 1291818 | 0 | 2332 | 0.00045 | 1.05 | True |
| 3 | 3457 | 20.5 | 8 | no | 2 | `q-curve-psi[d=3,lambda=717085]` | +d | 746891 | 0 | 3436 | 0.00029 | 0.99 | True |
| 3 | 1597 | 21.3 | 1 | no | 3 | `q-curve-psi[d=3,lambda=1598]` | +d | 6080 | 0 | 1632 | 0.00064 | 1.04 | True |
| 3 | 2269 | 22.3 | 1 | no | 1 | `q-curve-psi[d=3,lambda=2480057]` | +d | 37324 | 0 | 2284 | 0.00044 | 1.01 | True |

**Summary**

| d | instances | maps | log2 r | λ² ≡ +d / −d | ord_r(λ) min–max | ord_r(λ) / r, min | images in base, total | expected at chance, total | base points, total | all verified |
|--:|--:|--:|:--|:--|--:|--:|--:|--:|--:|:--|
| 2 | 20 | 20 | 12.1–22.3 | 10 / 10 | 2148–2572032 | 0.016 | 16 | 19.2 | 20664 | yes |
| 3 | 20 | 20 | 13.2–22.3 | 15 / 5 | 1817–1291818 | 0.002 | 16 | 20.1 | 22646 | yes |

Reading.  **The prediction holds: a Q-curve endomorphism keeps no base.**  Over `40`
instances (`20` per degree, `r = 2^{12}`–`2^{22}`) every `ψ` verifies, `λ² ≡ ±d` on every one
(`+d` on `10` / `−d` on `10` at `d = 2`, `15` / `5` at `d = 3`), and the eigenvalue orders run
from `1817` to `2.6·10⁶` — at least `0.2 %` of `r`, against `2` and `4` for the automorphisms
and `n` for a Frobenius — so there is no small orbit to fold by.  The images of the
`43310` base points land back in the base `32` times, against `39.3` expected at chance
(`16` against `19.2` at `d = 2`, `16` against `20.1` at `d = 3`).  The search found `174`
Q-curve `j`-invariants and counted `324` curves to get the `40` with a cofactor at most `16`;
every kernel that reached `E^σ`'s `j` did so over `F_{p²}` (`0` needed `F_{p⁴}` for `ι`), and
`41` roots of the division polynomial led to other `d`-isogenous neighbours and were set aside.
This is E4's type-C reading on the `F_{p²}` group: a degree-`d` isogeny's `x`-map is a rational
function of degree `d`, so no line or small set is stable, and a fast endomorphism is a scalar
multiplication shortcut (GLV), not a factor-base symmetry.

**Class: accounting.**  The lever list of §8.1 is closed on the base side: on these families
the only folds are the automorphisms, the Frobenius-type maps (`π`, `τ`, GLS `ψ`) and, with a
rational `2`-torsion point, the translation `τ_T` (E13).

### 8.4 Verdict after E8–E15 (with E14)

| lever | measured | class |
|:--|:--|:--|
| three summands on the `π`-line (E8) | columns `÷ 3`, relations `÷ 4.2`, `0` single-column rows; `S ÷ 3.8` end to end in one unit | engineering |
| the matched folded rho on `F_{p²}`, `F_{p³}` (E9) | steps `÷ 1.40, 1.74, 2.16` at `A = 4, 6, 12` (expected `1.41, 1.73, 2.45`), verified; in `F_p` multiplications the folded walk on `E(F_{p³})` pays `110` a step against `88` (E8) | accounting (reference) |
| the algebraic `S₄` line oracle (E11) | `9.42e+05` multiplications a target, flat in `p`; same relation stream as E8; exponent of the total `0.53` (fold) against the pair table's `0.88`; cheaper than the pair table from `p ≈ 2^12` | engineering |
| exponents (E10) | fold `0.88`, control `0.86`, rho `0.37` / `0.39` over `18` instances; no arm's exponent below rho's `1/2` | accounting |
| the pair table over orbit representatives (E12) | table `÷ 2.95` / `2.97` / `1.98` at `w = 6, 6, 4`, build `÷ 2.59` / `2.77` / `1.91`; probe `1.11×` (line) and `1.00×` (`F_p`), under the `1.2×` falsification line; `S ÷ 1.18` on `j = 0`, `÷ 1.03` on `j = 1728`, `× 1.02` on the line; `0` collisions in `9.77e+07` probes | engineering (memory); relabelling if quoted as a cost gain on the line |
| FGHR's `2`-torsion symmetry with the fold (E13) | system: `R(q₃)` of degree `16 = 64/4`, `1.25e+04`–`2.04e+04` multiplications a call (`÷ 33`–`41` against `S₃`), `S ÷ 27.0`; base: `τ_T` halves the columns, `S ÷ 2.52`; combined / product `≥ 0.82` (survives); `267`–`5741×` rho, total `p^{1.23}` | engineering |
| the ECC2K-130 family, two summands (E15) | `E_0` at `4·prime` (`m = 23`): `D = 1` / `1`, `4 %` single-column rows, the fold buys relations; at `m = 31` (`h = 4·373`): `D = 12` / `342`, `73 %` single-column rows — the §6.5 confinement follows the cofactor, and `m = 31` cannot carry the challenge's shape | accounting |
| Q-curves of degree 2 and 3 over `F_{p²}` (E14) | `40` instances, every `ψ = π ∘ ι ∘ φ` verified, `λ² ≡ ±d`, `ord_r(λ) = 1817`–`2.6·10⁶`; `32` images of `43310` base points in the base against `39.3` at chance | accounting |

With E13 run, the cheapest arm in this note is the `Y`-line of a subfield curve with rational
`2`-torsion, folded `12` a column by `⟨−1, π, τ_T⟩` and decomposed by the `D₃`-symmetrised
conic-resultant oracle: `267×` rho at best, every phase priced, the solver three quarters of
the cost.  Both uses of the `2`-torsion are now taken, on the base and on the system, and they
compound with the Frobenius fold; the relation phase is still `p^{1+o(1)}` solves, so what
remains is the constant of a solve (a degree-`16` univariate root-finding and three conics) and
not a lever on the base.  E12 settles the pair-table regime: the table folds by `w/2`, and that
is worth having only where the table is the cost (`F_p`, two summands).  E15 settles the transfer question for the base: on the challenge's `4·prime` shape the Frobenius-fixed cofactor does not confine two-summand relations, but on `m = 31`, the program's primary exploratory size, the `373` cofactor does.  E14 closes the base side: a Q-curve endomorphism of degree 2 or 3 keeps base points only at chance, as E4's CM maps did, so the folds on these families are the automorphisms, the Frobenius-type maps and `τ_T`, and every experiment of §8.1 has run.
