# 6-sum (and 5-sum) relation extraction over the Frobenius pair table

**Status: stage diagnostic, measured 2026-10-12.**  Relation-extraction
cost per verified relation for four higher-arity kernels against the
4-sum baseline, on public Koblitz K_0 rungs, in three units at once:
`S₃`/`S₄` root solves, first-summand probes (the unit of
`RESEARCH_ECC2K130_IC_FEASIBILITY.md` §4.7), and targets.  No
end-to-end cost, no `S`, no rho comparison, no key-recovery, sub-rho or
challenge-scale claim.  Wall figures come from a shared host at load
37–165 and are secondary.

Code: `src/cryptanalysis/koblitz_sixsum_extract.rs` (kernels, table,
`S₄` quartic + root finder, fail-closed verification, framework oracle
`sixsum-extract`), `examples/koblitz_sixsum_extract_probe.rs` (the
measurement, JSON out).  Evidence: `results.json` (merged, every cell
with its source file) and `raw/` (each run's JSON and stderr log).

## The question

The compact-orbit index calculus finds 4-sum relations
`R = ±σ^a u ± σ^b u'` by scanning the pair table (`A = K²n²` abscissa
states `x(σ^t u_e)`): one quadratic `S₃(x_R, x_s, X)` solve per state,
both roots looked up.  §4.7 of the feasibility note models a 6-sum
extraction at `45·r/(n⁴K⁴)` probes per relation against the 4-sum's
`r/(n²K²)` and calls the kernel "the concrete open design".  The point
to settle: a 6-sum probe `(s₁, s₂)` tests every state for `s₃` at once
exactly as a 4-sum probe `s₁` tests every `s₂`, so in `S₃`-solve units
the per-probe success may be identical and the gain only combinatorial
multiplicity.  Which unit was §4.7 in, and what does each kernel cost
in each unit?

## Kernels

| kernel | loop | solve per probe | points |
|:--|:--|:--|:--|
| `4sum` | every abscissa state `s` of `R` | `S₃(x_R, x_s, X)`, 2 roots | 4 |
| `6sum_a_iterate_table` (A) | every signed shifted state `s₁`: `R' = R − s₁` (one point addition), then the 4-sum scan of `R'` with early abort | `S₃`, 2 roots | 6 |
| `6sum_b_s4_per_probe` (B) | every unordered pair `(s₁, s₂)` of abscissa states, **row order** (`s₁` fixed, `s₂` running) | quartic `S₄(x_R, x_{s₁}, x_{s₂}, X)`, 4 roots `x(R ± s₁ ± s₂)` | 6 |
| `6sum_b2_s4_strided` (B2) | the same pairs, **diagonal order** `(i, j) → (i+1, j+1)`, so consecutive quartics share neither `s₁` nor `s₂` | `S₄`, 4 roots | 6 |
| `5sum_c_control` (C) | every base point `P_e`: `R − P_e`, then the 4-sum scan | `S₃`, 2 roots | 5 |

B2 was added after the first measurements showed why A and B gain
nothing (below); it is possible only because the quartic needs
abscissae alone, with no `R − s₁` point addition per pair.

Every candidate is lifted to points, its signs resolved by group
arithmetic (`resolve_pair_signs`, `resolve_triple_signs`), and the
whole index list re-added in the single-word *and* the general
`KoblitzCurve` arithmetic before it counts (`verify_relation`); the
relation row is formed over the base's signed-Frobenius-orbit columns
(`(column, coefficient mod r)`), and the build checks the column model
(`[coef_i]·rep = [coef_rep]·P_i` for every base point).  Across all 47
cells: `verification_rejects = 0`, `lift_failures = 0`, every relation
has `support = arity` except one B2 relation at n=41 K=20 whose six
points cancel to a 4-sum (`support_histogram` in `results.json`).

## Results at n = 41 (`icv1-f2m41-tm2308219-7f48b14a`, r = 549 756 390 943 ≈ 2³⁹·⁰, cofactor 4)

Base `B = 2nK` points, table `A` abscissa states: K=5: B=410, A=41 410;
K=10: B=820, A=166 870; K=20: B=1 640, A=669 940.  Seed 20261012
(20261013 for the two 600 s cells), early abort, ≥32 relations where
the 240 s cell wall allowed.  "solves" are `S₃` solves for `4sum`, A,
C and `S₄` solves for B, B2.  "struct." is the structural expectation
of the previous row group's unit (formula below); "§4.7" is
`45·r/(n⁴K⁴)` in first-summand probes.

| K | kernel | rel | targets | solves/rel (measured) | struct. solves/rel | meas./struct. | outer probes/rel | struct. outer/rel | §4.7 model | targets/rel | wall/rel (s, contended) | cell |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| 5 | 4sum | 0 | 884 | — (36.6 M S₃, 0 rel) | 1.94·10⁷ | — | — | 1.94·10⁷ | 1.40·10⁴ | — | — | 240 s, budget |
| 5 | 4sum (600 s, seed 20261013) | 2 | 2198 | 4.55·10⁷ S₃ | 1.94·10⁷ | 2.3 | 4.55·10⁷ | 1.94·10⁷ | 1.40·10⁴ | 1099 | 300 | 600 s, budget |
| 5 | A | 4 | 5 | 1.37·10⁷ S₃ | 1.94·10⁷ | 0.71 | 332 | 467 | 1.40·10⁴ | 1.25 | 60 | 240 s, budget |
| 5 | B (row) | 0 | 1 | — (0.87 M S₄, 0 rel) | 9.69·10⁶ S₄ | — | — | 9.69·10⁶ | 1.40·10⁴ | — | — | 240 s, budget |
| 5 | B2 (stride) | 5 | 6 | 1.12·10⁶ S₄ | 1.59·10⁶ S₄ | 0.70 | 1.12·10⁶ | 1.59·10⁶ | 1.40·10⁴ | 1.20 | 48 | 240 s, budget |
| 5 | C | 0 | 3 | — (46.6 M S₃, 0 rel) | 1.94·10⁷ | — | — | 467 | 1.40·10⁴ | — | — | 240 s, budget |
| 5 | C (600 s, seed 20261013) | 0 | 6 | — (91.9 M S₃, 0 rel) | 1.94·10⁷ | — | — | 467 | 1.40·10⁴ | — | — | 600 s, budget |
| 10 | 4sum | 6 | 227 | 6.19·10⁶ S₃ | 4.93·10⁶ | 1.26 | 6.19·10⁶ | 4.93·10⁶ | 875 | 37.8 | 40.0 | 240 s, budget |
| 10 | A | 7 | 8 | 7.55·10⁶ S₃ | 4.93·10⁶ | 1.53 | 46.1 | 29.7 | 875 | 1.14 | 34.3 | 240 s, budget |
| 10 | B (row) | 0 | 1 | — (0.92 M S₄, 0 rel) | 2.47·10⁶ S₄ | — | — | 2.47·10⁶ | 875 | — | — | 240 s, budget |
| 10 | B2 (stride) | 21 | 22 | 2.66·10⁵ S₄ | 4.03·10⁵ S₄ | 0.66 | 2.66·10⁵ | 4.03·10⁵ | 875 | 1.05 | 11.4 | 240 s, budget |
| 10 | C | 6 | 7 | 7.84·10⁶ S₃ | 4.93·10⁶ | 1.59 | 47.8 | 29.7 | 875 | 1.17 | 40.0 | 240 s, budget |
| 20 | 4sum | 32 | 72 | 1.01·10⁶ S₃ | 1.07·10⁶ | 0.94 | 1.01·10⁶ | 1.07·10⁶ | 54.7 | 2.25 | 5.59 | 179 s, 32 rel |
| 20 | A | 32 | 32 | 7.91·10⁵ S₃ | 1.07·10⁶ | 0.74 | 1.97 | 2.36 | 54.7 | 1.00 | 4.06 | 130 s, 32 rel |
| 20 | B (row) | 4 | 5 | 4.83·10⁵ S₄ | 4.24·10⁵ S₄ | 1.14 | 4.83·10⁵ | 4.24·10⁵ | 54.7 | 1.25 | 60.0 | 240 s, budget |
| 20 | B2 (stride) | 32 | 32 | 1.32·10⁵ S₄ | 1.01·10⁵ S₄ | 1.30 | 1.32·10⁵ | 1.01·10⁵ | 54.7 | 1.00 | 6.42 | 205 s, 32 rel |
| 20 | C | 26 | 27 | 8.55·10⁵ S₃ | 1.07·10⁶ | 0.80 | 2.08 | 2.36 | 54.7 | 1.04 | 9.23 | 240 s, budget |

Note's 4-sum model `r/(n²K²)`: 1.31·10⁷ (K=5), 3.27·10⁶ (K=10),
8.18·10⁵ (K=20) probes per relation; measured 4-sum `S₃` solves per
relation sit at 1.2–3.5× that (fit constant `C ≈ 1.2–3.5`, against the
note's 0.33–1 on bigger rungs; two of the three cells have ≤ 6
relations).

Solve rates on this host (contended): `S₃` 1.5–2.3·10⁵/s in the
load-100 runs; `S₄` 3.6–8·10³/s in the load-100 runs and 2.1–2.3·10⁴/s
in the load-40 B2 runs.  **One `S₄` solve cost about 20–60 `S₃`
solves** here (not isolated; the `S₄` root finder is a `Vec`-based
polynomial routine, `n` squarings modulo the quartic plus trace
splitting, and has had no tuning).

### n = 23 and n = 31 (saturated)

At n=23 (r = 2 095 853) and n=31 (r = 1 439 393, cofactor 1492) the
table covers the subgroup densely (`λ₄ = B⁴/(24r)` from 56 to
3.5·10⁸): every kernel finds a relation inside the first target, with
1.3–55 solves per relation, measured within 0.75–1.3× of the dense
structural value (`r/(4A)` for `S₃` kernels, `r/(8A)` for B, B2) in
all 30 cells (`results.json`, `raw/n23_31_v3.json`).  These rungs
cannot separate the kernels; they confirm the per-solve hit model
(`4A/r` per `S₃` solve, `8A/r` per `S₄` solve).

## What the numbers say

**1. In `S₃`-solve units, kernel A ≡ the 4-sum ≡ the 5-sum control.**
At K=20 (the only column with ≥ 26 relations in every `S₃` cell):
1.01·10⁶ (4sum), 0.79·10⁶ (A), 0.85·10⁶ (C) `S₃` solves per relation,
all within noise of the same structural value 1.07·10⁶.  At K=10:
6.2·10⁶, 7.6·10⁶, 7.8·10⁶.  The reason is structural and the data
match it: kernel A's inner loop for a fixed `s₁` *is* the 4-sum scan
of the pseudo-random point `R − s₁`, so it finds a relation with the
4-sum's per-scan probability `1 − e^{−λ₄}` at the 4-sum's cost `A`
solves per empty scan.  The abundance of 6-sums per target
(`λ₆ = B⁶/(720r)` = 12, 768, 49 154 at K = 5, 10, 20) does not help,
because the scan order groups them by `s₁`: a row has hits only when
`R ± s₁` has a 4-sum, and then it has six of them.  The same holds for
C with `R − P_e`, and for B in row order (half the solves — the
quartic covers both signs of `s₁` — each many times dearer).

**2. What arity buys is targets per relation and first-summand
probes, which is the unit §4.7 was in.**  Measured targets per
relation: 4-sum 2.25 / 37.8 / ≥ 884 (K = 20 / 10 / 5) against 1.0–1.25
for A, B2 and C.  Measured first-summand probes per relation for
kernel A: 1.97 / 46 / 332, equal to the 4-sum's *targets* per relation
(structural `1/(1 − e^{−λ₄})` = 2.36 / 29.7 / 467), because each outer
probe is one 4-sum scan.  §4.7's `45·r/(n⁴K⁴)` = 54.7 / 875 / 14 008
is a count in this unit (it divides the 6-sums per target
`B⁶/(720r)` into `4K²n²` first-summand states), pessimistic by 20–40×
in its own unit (it ignores the 90-fold multiplicity with which a
6-sum is met) — and the unit hides `A = K²n²` `S₃` solves per probe
(4.1·10⁴ / 1.7·10⁵ / 6.7·10⁵ here).  Converted to `S₃` solves, §4.7's
figure is `45·r/(n²K²)`, i.e. 45× *worse* than the 4-sum's own model,
and the measured kernel A sits at the 4-sum's cost, not below it.

**3. The one real lever is the scan order, and it is paid in `S₄`
solves.**  Kernel B2 de-clusters the hits by walking the pair square
diagonally; its solves per relation follow the dense model `r/(8A)`
in every column — 1.12·10⁶ / 2.66·10⁵ / 1.32·10⁵ `S₄` solves at K = 5 /
10 / 20 against structural 1.59·10⁶ / 4.03·10⁵ / 1.01·10⁵ — which is
7.6× (K=20), 23× (K=10) and 17–40× (K=5, 2–5 relations) fewer solves
than the 4-sum's `S₃` count, and structurally `r/(8A)` against the
sparse 4-sum's `1.5·r/A`: a factor 12 in solve count in the regime
where the 4-sum is sparse (`λ₄ ≪ 1`, which is every large rung), and
only 2 where it is dense.  Each `S₄` solve cost 20–60 `S₃` solves
here, so in `S₃`-equivalents B2 is 1–8× *more* work than the 4-sum at
these cells; it breaks even only if a quartic root solve in `F_{2^n}`
costs under ~12 `S₃` solves (sparse regime) — an implementation
question, not measured here.  B2's contended wall per relation (6.4 s
at K=20 against the 4-sum's 5.6 s; 11.4 s against 40 s at K=10; 48 s
against 300 s at K=5) is confounded: the B2 processes ran at load ≈ 40,
the others at load ≈ 100–165, so the wall column is not evidence of a
speedup.

**Verdict.**  §4.7 counted first-summand probes (equivalently, 4-sum
*targets*), not `S₃` solves; in the ladder's unit the 6-sum kernels A
and B and the 5-sum control cost what the 4-sum costs per relation, to
within the noise of 4–32 relations per cell, and the multiplicity shows
up only as targets per relation (≈ 1 against 2–1100).  The only kernel
that converts the multiplicity into fewer solves is the stride-ordered
`S₄` kernel B2, at `r/(8A)` quartic solves per relation, which at the
measured `S₄`/`S₃` cost ratio is not cheaper than the 4-sum.  Nothing
here is an end-to-end or rho-relative statement; the linear algebra
and the rank stage are untouched, and 6-point rows are denser than
4-point rows.  Classification under AGENTS.md §3: **accounting** for
the §4.7 unit (the model's unit named and converted), **engineering
candidate, unestablished** for B2 (a solve-count gain whose cost per
solve exceeds the gain).

## Structural model used in the table

For a kernel whose hit-or-miss *scan unit* has `U` solves, finds a
relation at Poisson rate `λ` per unit and meets each relation `μ`
times inside the unit at uniformly placed positions,

```
solves per relation ≈ U · ( e^{−λ}/(1 − e^{−λ}) + 1/(μλ + 1) )
```

(empty units paid in full, then the first of `μλ` hits), with
asymptotes `U/λ` (`λ ≪ 1`, sparse) and `U/(μλ)` (`λ ≫ 1`, dense).
Per kernel: 4-sum `U = A`, `λ = λ₄ = B⁴/(24r)`, `μ = 6`; A and C the
same (the unit is the inner 4-sum scan); B (row) `U = A` `S₄` solves,
`λ = 2λ₄`, `μ = 6`; B2 `U = A²/2` `S₄` solves, `λ = λ₆ = B⁶/(720r)`,
`μ = 45`.  Dense asymptotes: `r/(4A)` `S₃` solves (two roots, each in
the table with probability `2A/r`) and `r/(8A)` `S₄` solves.  Targets
per relation: `1/(1 − e^{−B^m/(m!·r)})`.  `results.json` carries every
cell's `structural.*` values and `measured_over_structural` (0.66–1.6
over the 13 n=41 cells with ≥ 4 relations; 2.3 on the 2-relation cell).

## Tests (`cargo test --release --lib koblitz_sixsum_extract`, 6 passed)

- `s4_quartic_matches_reference_and_roots_match_brute_force`: the
  quartic agrees with `binary_semaev::binary_semaev_s4` on random
  inputs at n=17; `s4_roots` (gcd with `X^{2^n} − X` + trace
  splitting) equals brute force over all 2¹⁷ field elements; the
  faster `s4_roots_split` reports only genuine roots and all of them
  when the quartic splits; on `x_R, x_{s₁}, x_{s₂}` from rational
  points the roots are exactly `x(R ± s₁ ± s₂)`.
- `the_table_states_and_lookups_are_consistent`: every state point is
  the sum of its two members, every shift is a base point, and the
  lookup of `x(σ^t u_e)` returns a state with that abscissa (n=23).
- `planted_targets_are_recovered_by_every_kernel`: planted 4-, 5- and
  6-sums at n=41 recovered by the 4-sum, A, B, B2 (from the planted
  cell and from five diagonal steps before it: exactly six quartics)
  and C; the recovered index lists equal the planted ones and re-add to
  the target; `verification_rejects = 0`.
- `non_decomposable_targets_yield_no_relation`: a target whose
  logarithm is outside the 4-, 5- and 6-sum scalar lattice of a
  one-column base (checked on the coefficients) yields nothing from
  complete 4-sum, A and C scans and from budgeted B and B2 scans; a
  target off the prime-order subgroup (`R + T`, `T` the 2-torsion
  point) yields nothing.
- `verification_rejects_what_does_not_close`: wrong target, tampered
  index, out-of-range index rejected; `P + (−P)` cancels in `support`.
- `random_targets_yield_verified_relations_at_a_small_rung`: 4-sum and
  A on random targets at n=23, every relation re-adds, no rejects.

`cargo check --release --lib`: clean (one pre-existing unused-import
warning elsewhere).

## Reproduce

```sh
export CARGO_TARGET_DIR=target
cargo test --release --lib koblitz_sixsum_extract
cargo build --release --example koblitz_sixsum_extract_probe
B=target/release/examples/koblitz_sixsum_extract_probe
# saturated rungs, all kernels (seconds)
$B --n=23,31 --columns=5,10,20 --relations=32 --cell-wall-s=120 > raw/n23_31_v3.json 2> raw/n23_31_v3.stderr.txt
# n=41, one process per column (kernels 4sum,a,b,c; 240 s cell wall)
for K in 5 10 20; do
  $B --n=41 --columns=$K --kernels=4sum,a,b,c --relations=32 --cell-wall-s=240 > raw/n41_k$K.json 2> raw/n41_k$K.stderr.txt
  $B --n=41 --columns=$K --kernels=b2 --relations=32 --cell-wall-s=240 > raw/n41_k${K}_b2.json 2> raw/n41_k${K}_b2.stderr.txt
done
# the two longer K=5 cells
$B --n=41 --columns=5 --kernels=4sum --relations=32 --seed=20261013 --cell-wall-s=600 > raw/n41_k5_4sum_long.json 2> raw/n41_k5_4sum_long.stderr.txt
$B --n=41 --columns=5 --kernels=c    --relations=32 --seed=20261013 --cell-wall-s=600 > raw/n41_k5_c_long.json 2> raw/n41_k5_c_long.stderr.txt
```

`results.json` was merged from the raw files with `jq` (the
`structural.*` block is recomputed there from each cell's `r`, `A`, `B`
and kernel with the formula above, so v1 and v3 producer outputs carry
the same derived fields).  Solve counts and targets are deterministic
given the seed up to where the wall deadline cuts a cell; wall figures
are not reproducible on a shared host.

The framework oracle `sixsum-extract` (`kernel=4sum|a|b|b2|c`) is
wired into `ic_run`'s binary analogue for later ladder use; it was not
run end to end here.

## Not done / caveats

- Cells with 0–7 relations (all K=5 `S₃` cells, K=10 4-sum/A/C, K=20
  B) are far below the 32-relation target: the 240–600 s cell walls on
  a host at load 100–165 were the limit.  Their solve counts are
  consistent with the structural model but have wide intervals; the
  K=5 5-sum control found nothing in 138 M `S₃` solves (≈ 7 expected
  at the structural rate; a 0.1 % outcome, or an unmodelled effect at
  K=5 — unresolved).
- No isolated timing (`tools/isolated_bench.py`), no A/A run; the wall
  and the `S₄`/`S₃` cost ratio are contended and unpooled.
- The `S₄` root finder is untuned; the break-even question for B2
  (quartic root solve in `F_{2^n}` under ~12 `S₃` solves) is open.
- No rungs beyond n=41, no m=83 gate, no end-to-end run, no linear
  algebra for the denser 6-point rows, no scoreboard or leaderboard
  change (stage diagnostic only).
- Research-visuals (graph + PDF) not produced in this session.
