# Protocol: a supersingular auxiliary curve for PKM's isogeny tower (2026-10-15)

Frozen before any measured cell of Q2 or Q3 and before the Q1 computation
that this study reports. Lane H. Branch `claude/pkm-supersingular-aux`
(from `claude/isogeny-tower-p192`, PR #1925).

**Why.** Over the P-192 prime (`p ≡ 3 mod 4`, `v₂(p + 1) = 64`) the curve
`y² = x³ + x` is supersingular, `#E = p + 1`, with a cyclic 2-Sylow subgroup of
order `2^64`; a point of order `2^48` takes one draw (Lane G,
`research/pkm_isogeny_tower_p192_20261015/`). The backlog card
`XP-20261015-p192-isogeny-aux-search-t48` prices an ordinary curve with such a
point at about 11.2 core-years. This study decides whether a supersingular `E'`
serves PKM's isogeny-tower factor base as well as an ordinary one: in the
structure of `V` (Q1), in the F4 solving degree of the decomposition system
(Q2), and in the relation phase under ecbench's charge (Q3).

## Code (committed before this protocol)

- `src/cryptanalysis/ic_framework/prime_tower.rs` (`b882ce759`):
  `AuxKind`, `supersingular_aux`, `supersingular_tower`, and the plug-in option
  `prime-isogeny:t=T[,seed=S],aux=supersingular|ss-j1728|ss-j1728-twist|ss-j0|ss-j0-twist`.
  The ordinary build is factored into `finish_tower` without change in what it
  counts: the committed sessions `q20` and `p22` of
  `research/pkm_isogeny_tower_ecbench_20261015/` replay identically with this
  build (`ecbench verify --replay-all`: 120/120 and 144/144 verified, 80 and 96
  replays reproduced, audit OK).
- Tests (`cargo test --release --lib prime_tower`, 9 pass): `#E' = p + 1` by
  exhaustive count on every supersingular kind at every prime 67–3000 where it
  applies (and refusal exactly where `y² = x³ + x`, `y² = x³ + 1` are
  ordinary); `|V| = 2^t`, every element a root of the composed Vélu x-map and,
  exhaustively over `F_p`, no other root; `L(x)` recomputed over `BigUint`;
  the order of `Q` checked exactly `2^t` over `BigUint`; the charge (one curve,
  zero search, every operation in `aux_tower_build_gae`).
- `examples/pkm_ss_tower.rs`: `structure` (Q1) and `f4` (Q2). Its system builder
  is copied from `examples/pkm_tower_pilot.rs` (example-private there, and
  that file is being edited by another lane); the isogeny towers come from the
  library. Tests: canonical forms against brute force over `PGL₂(F_p)`; the
  coset count of `y² = x³ + x`; F4's consistency against the exhaustive recount
  on four arms.
- `examples/prime_tower_curves.rs`: `--plus` for `p = c·2^s − 1`.

## Disclosure: what was run before freezing

- **Q1 dry runs (scratch, not committed).** `structure` at `p ∈ [100, 700]`,
  `t = 2, 3` with the brute-force check, and at `p ∈ [1000, 1300]`, `t = 2…4`.
  They showed: at `t = 2` every pair of families shares some Möbius classes,
  ordinary pairs included; at `t ≥ 3` no supersingular set was equivalent to a
  torus set; two ordinary controls at `p = 1151` with equal `j = 24` were
  equivalent (isomorphic curves). The brute force agreed with the canonical
  form on every checked set. The positive-control family `ss-j1728-iso16` was
  added after this dry run.
- **Q2 dry run (scratch, not committed).** At `T24`, `m = 2`, `N = 2…12`, two
  random targets: `ss-j1728` and `ordinary` gave `D = 3, 5, 5, 5, 5, 5`;
  `torus` gave `3, 4, 4, 4, 4, 5`, with round 8's widths. So the prediction
  below for N ≤ 12 at `T24` is informed by data; every committed row is a new
  run of the frozen command.
- **Q3 sizing (committed, `q3/sizing.txt`).** `ecbench fb` on the three curves,
  both aux kinds: abscissae 141 / 123 (s20), 246 / 235 (s22), 253 / 251 (s24)
  for supersingular / ordinary; the ordinary search tried 71, 684 and 1,043
  curves. No relation was collected.

## Q1 — structure (exact computation)

**Question.** For `p ≡ 3 mod 4` with `2^t | p + 1`, is
`V_ss = x(R' + ⟨Q⟩)` on a supersingular `E'` (`y² = x³ + x`, `y² = x³ − x`, and
`y² = x³ + 1`, `y² = x³ + d³` when `p ≡ 2 mod 3`) equal, up to a Möbius (or
affine) change of `x` over `F_p`, to a Dickson torus set
`{w + w^{−1} : w ∈ u·μ_{2^t}}` at the same `p`, `t`?

**Instrument.** `pkm_ss_tower structure --p-min 100 --p-max 3000 --t 2,3,4,5
--ordinary 3 --brute 131`. For every prime `p ≡ 3 mod 4` in range and every
`t` with `t + 1 ≤ v₂(p + 1)`, every set of every family (all torus cosets; all
cyclic subgroups of order `2^t` and all `R'` with `2R' ∉ ⟨Q⟩` on each curve),
a complete `PGL₂(F_p)` invariant (canonical form), and for every ordered pair
of families the number of sets of the first equivalent to some set of the
second. Families: `torus`, the applicable supersingular kinds, the
positive control `ss-j1728-iso16` (`y² = x³ + 16x ≅ y² = x³ + x`), and three
ordinary curves with a point of order `2^t` (negative reference).

**Decision rule.** For `t ≥ 3` (at `t = 2` a 4-point set has one modulus, so
coincidences are expected and reported only):
- **H1_equiv** (the supersingular route *is* the torus route): every set of
  every supersingular family is equivalent to a torus set, at every
  `(p, t ≥ 3)`.
- **H1_distinct**: no supersingular set is equivalent to any torus set at any
  `(p, t ≥ 3)`.
- otherwise **partial** (coincidental), reported with counts.
- **Positive control**: `ss-j1728-iso16 → ss-j1728` must be 1.0 at every cell;
  otherwise the instrument is invalid. **Brute-force agreement** must be 100 %
  where run.

Prediction (derived before the full run): **H1_distinct**. A Möbius map between
the two fibres would conjugate fibres of the degree-`2^t` maps; the Dickson map
`D_{2^t}` is a polynomial (totally ramified over `∞`), while the composed Vélu
x-map of a cyclic `2^t`-isogeny has ramification index at most 2 everywhere
(`x(P + K) = x(−P + K)` only when `2P ∈ K`). A single fibre can still match by
coincidence, which the computation counts.

## Q2 — F4 solving degree (stage diagnostic, measured)

**Instrument.** `pkm_ss_tower f4`, engine `f4_fp_tower` (library, this
branch), reduced presentation, full-rank exit, degree bound `n + d + 6`, two
random targets per cell, seed `0x53534157` (default), tower seed 0,
`RAYON_NUM_THREADS=1`, `nice -n 10`, address space ≤ 4 GB, budget 3,600 s per
system, `--max-nnz 2·10⁹`. `D = solving_degree_max`; width =
`max_cols_to_solution`; steps = `steps_to_solution`. Every finished tower row
is re-checked exhaustively over `V^m` (`verify` field; `m = 2` exact, `m = 3`
by resultant); a mismatch stops the run.

**Primes.** `T24 = 7·2^21 − 1`, `T31 = 61·2^25 − 1` (round 8's torus primes;
`≡ 1 mod 3`, so no `j = 0` arm), `S25 = 9·2^21 − 1 = 18,874,367` (`≡ 2 mod 3`,
so `j = 0` applies), and `p₀ = 786,433 = 3·2^18 + 1` (round 8's positive
control; `p₀ ≡ 1 mod 4`, no supersingular arm).

**Arms.** (a) supersingular `ss-j1728`, `ss-j1728-twist`, and at `S25`
`ss-j0`, `ss-j0-twist`; (b) `ordinary` (the library's seeded search,
asserted not supersingular); (c) `torus` (reference); (d) `null` on
`ss-j1728` (`S₃`'s support with random coefficients, same tower); positive
control at `p₀`: `ordinary` against `split` Dickson.

**Cells** (`run_q2.sh`, in this order; no run starts after 180 min):
1. P0: `p₀`, `ordinary` and `split`, `m = 2`, `t = 1…5`.
2. At each of `T24`, `T31`, `S25`: every arm, `m = 2`, `t = 1…9` (`N ≤ 18`).
3. At `T31`: `ss-j1728`, `ordinary`, `torus`, `m = 3` (`t = 2…4`) and `m = 4`
   (`t = 2, 3`); the `null` on `ss-j1728` at `m = 2` (`t = 1…9`) and `m = 3, 4`.
4. At `T31`: `ss-j1728`, `ordinary`, `torus`, `m = 2`, `t = 10` (`N = 20`).

Pre-flight (derived from the dry run's growth, about ×5 per level at `m = 2`):
`N = 18` about 15 s a system, `N = 20` about 1–5 min; the whole list under
90 min single-threaded.

**Decision rules** (per pair of arms, over the cells where both finished):
- *same*: `D` equal on both targets in every cell; *differs*: `D` differs on
  both targets in at least one cell; otherwise *mixed*.
- **H_same** (the brief's hypothesis) holds iff (a)=(b) *and* (a)=(c).
- **H2a** (supersingular = ordinary): (a) vs (b) *same*, at every prime.
- **H2b** (supersingular = torus): (a) vs (c) *same*.
- **Positive control**: at `p₀`, `ordinary` vs `split` must *differ*
  (round 8: 5 vs 4 at `N = 4…10`); if not, Q2 is invalid. Within each torus
  prime, (b) vs (c) differing is a second, same-prime positive control.
- The null gates nothing (round 8's reading): it says whether the `D` is the
  tower's or `S₃`'s.

**Prediction** (from the dry run and §3.4 of the PKM note): H2a *same*; H2b
*differs* (5 vs 4 at `N = 4…10`); so H_same fails. Reason: the supersingular
tower is a Vélu tower; its top-degree parts are `y_j·y_{j+1} − y_j²`, the
ordinary isogeny tower's, not the Dickson `y_j²`.

**What it licenses.** A supersingular `E'` is "as good as an ordinary `E'` for
F4" at these sizes if H2a holds. Nothing about larger `N` or about P-192's
security.

## Q3 — relation phase, charged (ecbench)

**Curves** (`prime_tower_curves --bits 20,22,24 --two-adicity 10 --plus --seed
20261016 --max-cofactor 4`, `q3/curves.jsonl`): `p ≡ 3 mod 4`,
`v₂(p − 1) = 1`, `v₂(p + 1) = 10` (the P-192 shape), registered in
`docs/curves/registry.json` in this change:

| name | slug | `p` | `r` (`h`) | `t` |
|:--|:--|--:|--:|--:|
| s20 | `icv1-fp20-t1143-5bef62bd` | 537,599 | 178,819 (3) | 8 |
| s22 | `icv1-fp22-tm1105-62a9d176` | 2,511,871 | 837,659 (3) | 9 |
| s24 | `icv1-fp24-tm4507-ffb4b7f9` | 11,750,399 | 11,754,907 (1) | 9 |

**Arms** (one session per curve, `q3/specs/SPEC-*.json`): `rho-neg`
(reference, `rho.negation`); `ic-ss-rootfind` (`prime-isogeny:t=T,aux=supersingular`,
`rootfind-mitm:m=3`); `ic-ss-mitm` (same base, `mitm:m=3,negation_folded=1`,
positive control); `ic-ord-rootfind` (`prime-isogeny:t=T`, the ordinary
search, seed 0); `ic-plain-ss` and `ic-plain-ord` (`prime-abscissa:size=` the
supersingular and ordinary base's abscissa count); `rho-neg-aa` (A/A
control). `max_trials 200000`, 8 one-target workloads, 2 measured rounds,
warmup 1, order alternate, L0 (`--cpus none --allow-busy`), operation counts
only. Target seeds `20261016012x`, measurement seeds `20261016312x`
(`x` = 0, 2, 4 for s20, s22, s24).

**Run and audit** (`run_q3.sh`): `ecbench run`, `ecbench verify --replay-all`,
`ecbench compare --save` for the pairs below, `ecbench table`.

**Metrics.** `S` (ecbench, GAE); paired totals ratios with ecbench's 95 %
cluster-bootstrap interval; the aux share of `S` (factor-base phase
`aux_search_gae + aux_tower_build_gae` over total GAE, mean over runs); the
relation yield `relations_found / targets_tried` per arm (pooled over runs).

**Decision rules.**
- **P1**: every run verified, every replay identical (else no result).
- **H3a** (the search cost vanishes): the supersingular arm records
  `aux_search_gae = 0` and an aux share ≤ 3 % of `S` on every curve.
- **H3b** (no relation-phase penalty): `ic-ss-rootfind ÷ ic-plain-ss` has an
  interval containing 1 or inside [0.97, 1.03] on every curve.
- **H3c** (cheaper than the ordinary arm): `ic-ss-rootfind ÷ ic-ord-rootfind`
  below 1 with the interval excluding 1 on every curve where the ordinary arm's
  aux share is ≥ 2 %. Base sizes differ (141/123, 246/235, 253/251), which moves
  `S` by a few per cent either way; the two plain arms measure that effect and
  the README reports `(ss ÷ plain-ss) ÷ (ord ÷ plain-ord)` beside it.
- **H3d** (yield): the supersingular arm's relation yield within ±10 % of
  `ic-plain-ss`'s on every curve.
- **Controls**: A/A rho interval contains 1 (negative); `ic-ss-rootfind ÷
  ic-ss-mitm` above 1 with interval excluding 1 (positive, as in PR #1889).

Prediction: P1, H3a, H3b, H3d hold; H3c holds on s22 and s24 (ordinary shares
≈ 2–8 % expected from 684 and 1,043 curves tried) and is unresolved on s20.

## Multiple comparisons, and what does not count

Q1 and Q2 decisions are exact or rule-based (no p-values). Q3 has 3 curves × 4
hypotheses; intervals are ecbench's 95 % per comparison, reported uncorrected,
and a hypothesis needs every curve, which is conservative. Changing a curve,
`t`, seed, arm, target count or rule after seeing data is not allowed; any
change is a labelled amendment with both versions reported.

## Deliverables

`README.md` (verdicts first), `report.typ` → `report.pdf` (inspected page by
page), graphs (SVG + PNG + data + plotting source) in `figures/`, run report
in `reports/run-<UTC>/`, raw rows (`q1/structure.jsonl`, `q2/runs/*.jsonl`,
`q3/sessions/*`), and a section in the P-192 backlog card's README (status
unchanged; `backlog validate --all`). No speed claim, so no `benchmarks/`
entry.
