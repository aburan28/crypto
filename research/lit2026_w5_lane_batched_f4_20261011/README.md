# lit2026 W5: lane-batched inherited-F4 replay, F4 cost-model check, gate-first F6-IC (2026-10-11)

Results report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`).
Run reports:

- [`reports/run-20261011T2225Z/report.pdf`](reports/run-20261011T2225Z/report.pdf):
  arms `orig`, `prev`, `gate`, `lanes`, `rref` on the four probe cells, plus the
  diagnostics (cost model, lockstep, profile, root-phase lane sweep);
- [`reports/run-20261012T0008Z/report.pdf`](reports/run-20261012T0008Z/report.pdf):
  arms `orig`, `prev`, `gate`, `rref`, `gaterref` (gate-first with the RREF root).

Raw records are in each run's `raw/`, the Rust-derived summaries in `summary.json`
(`f4_lane_probe summarize`), and the graphs (`*.svg`) are drawn from them by `plot.py`.

**Outcome.**

| task | kind | result |
|:--|:--|:--|
| 1. lane-batched replay (`KIC_F4_LANES=k`) | negative | exact, but 2.2–4.0× more stage units than the sequential root on the cells with a batch, and slower in CPU time |
| 2. F4-tracer cost-model check | diagnostic | root rows match the generic Boolean count exactly with occurring multipliers (0.67–0.70 of it with the production support-local ones); columns 0.28–0.49 of generic; rank one short of generic at degree 3 |
| 3. gate-first ordering (`KIC_F6_GATE_FIRST=1`) | speedup | F6-IC stage unit 1.57–1.78× lower than original on all four cells, verdicts and counters identical; F4 arm unchanged |

All switches stay off by default. Gate-first is a candidate for default-on in a
follow-up (see "Leads").

## Papers

- Demin, *Groebner.jl fast tracing*, arXiv [2607.06372](https://arxiv.org/abs/2607.06372):
  one F4 trace replays several specialisations at once with SIMD-packed coefficients.
- Kouba, Neiger, Safey El Din, *A complexity analysis of the F4 Gröbner basis algorithm
  with tracer data*, arXiv [2603.16378](https://arxiv.org/abs/2603.16378): per-degree
  Macaulay/staircase counts for the F4 tracer on generic systems over a field.
- W4 lead (PR aburan28/crypto#1687): in F6-IC, 67–74 % of node word operations are spent
  in nodes the geometry gate eventually refutes.

## Production path

`SolverEngine::InheritedF4 { max_degree: 3 }`: a splitting tree (`solve_rec` in
`src/cryptanalysis/koblitz_groebner.rs`) whose nodes read a support-local degree-3
Macaulay basis (`ReducedBasis`, `src/cryptanalysis/inherited_f4.rs`) built once at the
root (and after a linear elimination) and specialised down the tree. F6-IC adds the
exact geometry gate as the solver's node oracle. On this host (aarch64, no BMI2) the
lazy rewrite path runs; the eager `pext`/`pdep` rewrite is x86-only.

Profile (`diag-profile.jsonl`): basis specialisation is 82–97 % of the solver's build
time, root builds 1.0–16.8 %.

## Changes (all behind switches, default off)

### 1. Lane-batched roots: `KIC_F4_LANES=k`

- `ReducedBasis::from_systems_lanes` (inherited_f4.rs): the roots of up to 64 systems
  built together. Each lane's rows are generated exactly as a sequential root generates
  them, laid out on the union of the lanes' columns and **bit-sliced** (one lane word per
  matrix entry, lane width 8/16/32/64 bits by k), then reduced by one Gauss–Jordan pass
  in which each row operation is a masked XOR on all lanes at once.
- Pivots are chosen per lane (first free row holding the column in that lane). Lanes that
  share the pivot pattern share the pass; divergent lanes take extra passes in the same
  column. That is the divergence detection and the fallback in one, so no lane is ever
  wrong. Counted: pivot columns, pivot passes, columns shared by all lanes.
- Each lane's result is the unique RREF of its own Macaulay row space, which is what the
  sequential root computes under `KIC_F4_INHERIT_ROOT=rref`, row for row
  (`lane_roots_match_the_sequential_rref`, k = 1…64, every lane width).
- `prebuild_inherited_roots_lanes` / `prebuild_decomposition_roots` build a batch of
  targets' roots ahead; a decomposition's first root build takes its lane basis if its
  system and degree match exactly (`take_lane_root`), otherwise builds as before.
- Per-instance accounting: the elimination's lane-word XORs, converted to 64-bit word
  equivalents (a w-bit lane word counts w/64), split evenly over the lanes and charged
  to each lane's stage unit and to the new root-elimination counter.
- **Scope.** Only the root is lane-batched. The splitting tree below it (82–97 % of the
  time) is not: see "Why the tree was not lane-batched".

### 2. Gate-first ordering: `KIC_F6_GATE_FIRST=1`

In `solve_rec`, when a reduction forces variables, the default assigns each one and
specialises the node's bases by it, and asks the node oracle on the next loop pass. With
the switch the solver assigns all of them first (same order, same conflict checks, same
counters), checks the budget, asks the oracle, and specialises only if it continues,
marking the oracle as already asked. The oracle reads only the assignment and the
defined-variable mask, so every decision, oracle call, counter, root and relation is
the same; a refuted node skips the specialisations. `set_gate_first` overrides the
variable per thread; `gate_first_skips_thread` counts skipped nodes and
specialisations. (The gate was already asked before a *child's* entry specialisation,
`KIC_F4_EARLY_ORACLE`; this moves the remaining case, forced variables, in front.)

### 3. Diagnostics and counters

- `f4_rows_generated_thread`: Macaulay rows generated (root builds and completion rows)
  and the words they were packed into, which `word_xors` does not charge.
- `f4_root_elimination_thread`: the stage units charged to root eliminations.
- `set_node_trace` / `take_node_trace`: per-reduction assignment keys for the lockstep
  diagnostic.
- `macaulay_shape`: rows, columns, rank of a root Macaulay matrix at any degree.
- `examples/f4_lane_probe.rs`: modes `profile`, `lanes`, `lockstep`, `costmodel`,
  `summarize`.
- `examples/f6_ic_probe.rs`: honours `KIC_F4_LANES` (batches of k targets, root build
  timed with the batch), prints relation indices per target and, per arm, process CPU
  time, rows generated, root-elimination units and gate-first skips.

### Tests

- `inherited_f4::tests::lane_roots_match_the_sequential_rref`
- `koblitz_index_calculus::tests::w5_gate_first_and_lane_roots_keep_every_decomposition`
  (K1/2^17 dim 6, 4 targets, f4 and f6p: identical relation and counters with gate-first
  on and off and with lane-built roots; gate-first never costs stage units).

`cargo test --release --lib -- inherited_f4:: w5_gate_first f6_ic`: 16 passed, 1 ignored
(the existing isolated-process test).

## Commands

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w5
cargo build --release --example f6_ic_probe --example f4_lane_probe
# original arm: f6_ic_probe built at 4b6e09b78 in /Volumes/SSD990-2/crypto/wt/lit2026-baseline
OUT=<run dir> ROUNDS=5 sh scripts/run_arms.sh                       # run-20261011T2225Z
OUT=<run dir> ROUNDS=5 ARMS="orig prev gate rref gaterref" sh scripts/run_arms.sh  # run-20261012T0008Z
OUT=<run dir> sh scripts/diag.sh          # costmodel, lockstep, profile, lanes sweep
OUT=<run dir> sh scripts/costmodel.sh     # costmodel re-recorded after a JSON-format fix
f4_lane_probe summarize <run dir> > <run dir>/summary.json
python plot.py <run dir>                  # plots recorded values only
```

Each arm runs `f6_ic_probe <cell> --engines f4,f6p --reps 3`; the cells are
`--a 1 --n 17 --dim 6 --targets 8`, `--a 1 --n 23 --dim 7 --targets 4`,
`--a 0 --n 23 --dim 8 --targets 4`, `--a 0 --n 31 --dim 9 --targets 1 --budget 2000`
(the W4 and eager-rewrite cells).

## Run conditions

- Host: shared Apple M4 Pro, 14 cores, aarch64. 1-minute load 82.5–388.1 in run 1.
  Isolation L0, unpinned. A test build of the library ran during part of run 1 and the
  probe build during part of run 2; arm order rotates every round.
- Commits: candidate `d22fb3b52` (the binaries measured); later commits on the branch
  change only formatting, the probe's `summarize` mode and output format, a test, and the
  study files. Original: `4b6e09b78`.
- Curves are named by ICV1 slug, prime-order subgroup chosen by `KoblitzCurve::new`,
  standard-subspace base, m = 3: K1/2^17 dim 6 `icv1-f2m17-tm101-00378d4e`; K1/2^23 dim 7
  `icv1-f2m23-tm5197-1f85e9e1`; K0/2^23 dim 8 `icv1-f2m23-t5197-69e76b73`; K0/2^31 dim 9
  `icv1-f2m31-tm90707-c95f16f5`.

## Results

### Gate-first (run-20261011T2225Z; stage unit identical in all 5 rounds of every arm)

| cell | engine | original = previous | gate | original / gate | CPU previous / gate, median [min–max] |
|:--|:--|--:|--:|--:|--:|
| K1/2^17 | f6p | 12 635 914 | 8 061 140 | 1.5675 | 2.22 [2.11–2.28] |
| K1/2^23 | f6p | 51 918 433 | 30 697 846 | 1.6913 | 2.49 [2.35–2.51] |
| K0/2^23 | f6p | 124 949 801 | 71 049 603 | 1.7586 | 2.24 [2.16–2.29] |
| K0/2^31 | f6p | 278 237 634 | 156 423 046 | 1.7788 | 2.13 [2.04–2.52] |
| all four | f4 | (unchanged) | (unchanged) | 1.0000 | 0.98–1.04 (medians) |

Gate-first skipped 447 / 432 / 840 / 405 refuted nodes and 4 641 / 5 835 / 11 004 /
6 591 specialisations per arm pass (three repetitions) on the four cells. The inherited-F4 arm has
no node oracle, so the switch changes nothing there.

### Gate-first with the RREF root (run-20261012T0008Z)

Second session, 5 rounds, load 76.5–131.8: `orig`, `prev`, `gate`, `rref`
(`KIC_F4_INHERIT_ROOT=rref`), `gaterref` (both). Every stage unit of the repeated arms
equals run 1's.

| cell | engine | gaterref | original / gaterref | CPU previous / gaterref, median [min–max] |
|:--|:--|--:|--:|--:|
| K1/2^17 | f6p | 7 696 852 | 1.6417 | 2.85 [2.65–3.15] |
| K1/2^23 | f6p | 26 902 588 | 1.9299 | 3.49 [3.26–3.69] |
| K0/2^23 | f6p | 56 940 731 | 2.1944 | 3.41 [3.39–3.53] |
| K0/2^31 | f6p | 115 269 611 | 2.4138 | 3.69 [3.50–3.96] |
| K1/2^17 … K0/2^31 | f4 | = rref | 1.0650 / 1.1915 / 1.2859 / 1.4780 | 1.13 / 1.32 / 1.38 / 1.71 |

Gate-first alone in run 2: CPU previous/gate f6p 2.27 / 2.56 / 2.27 / 2.11, f4 0.95–1.00.

### Lane-batched roots (run-20261011T2225Z, `lanes` = `KIC_F4_LANES=k` + RREF root)

| cell (k) | engine | rref (sequential) | lanes | rref / lanes | CPU rref / lanes, median [min–max] |
|:--|:--|--:|--:|--:|--:|
| K1/2^17 (8) | f4 | 14 098 582 | 48 972 167 | 0.2879 | 0.89 [0.84–0.95] |
| K1/2^17 (8) | f6p | 11 600 586 | 46 474 171 | 0.2496 | 0.88 [0.87–0.93] |
| K1/2^23 (4) | f4 | 44 709 930 | 160 671 638 | 0.2783 | 0.86 [0.80–0.92] |
| K1/2^23 (4) | f6p | 43 312 007 | 159 273 715 | 0.2719 | 0.83 [0.78–0.98] |
| K0/2^23 (4) | f4 | 123 273 212 | 268 819 024 | 0.4586 | 0.97 [0.93–1.00] |
| K0/2^23 (4) | f6p | 89 627 327 | 235 173 139 | 0.3811 | 0.93 [0.92–0.98] |

The tree part (stage unit less root eliminations) of `lanes` equals that of `rref` in
every cell and engine: only the root differs. K0/2^31 has one target, so no batch.

Root phase alone (`diag-lanes-echelon.jsonl`, 5 repetitions, targets [1]G…[k]G):

| cell | k | sequential root units | lane word-equivalents | ratio | pivot columns shared by all lanes |
|:--|--:|--:|--:|--:|--:|
| K1/2^17 | 8 | 2 337 120 | 38 133 403 | 16.3 | 17 / 427 |
| K1/2^17 | 64 | 18 710 274 | 1 676 721 535 | 89.6 | 17 / 432 |
| K0/2^23 | 4 | 8 023 633 | 156 827 347 | 19.5 | 35 / 761 |
| K0/2^23 | 64 | 127 874 057 | 13 009 658 647 | 101.7 | 30 / 765 |
| K0/2^31 | 64 | 530 272 122 | 85 983 336 157 | 162.1 | 39 / 1 308 |

Over k = 2…64 on all four cells the lane root costs 16–162× the sequential echelon
root's stage units (11–101× the sequential RREF root's, `diag-lanes-rref.jsonl`), 2.5–33×
its CPU time and 1.2–25× its wall time (medians of 5; `lanes.svg`). The lanes share
their pivot row in 3–14 % of pivot columns; pivot passes per column grow almost
linearly with k.

### Why the tree was not lane-batched (lockstep diagnostic)

`diag-lockstep*.jsonl`: across the targets of a cell, 66–83 % (f4) and 94–97 % (f6p) of the
stage unit is spent at tree nodes (assignments) that two or more targets reach, and such
nodes are reached by all targets, or all but one at K1/2^17 (histograms in the raw
records). So the trees *do*
coincide in shape. But a lane-batched tree would need the lanes to share the pivot pattern
of the basis at those nodes, and the roots already share it in only 3–14 % of columns.
More fundamentally, over F₂ a packed row holds 64 coefficients per word, so a bit-sliced
row operation over k lanes moves `k` bits per column against `k · (1/64)` words per column
for k packed rows: at most parity at k = 64 with every lane active, and a loss at any
smaller k or any divergence. Demin's gain comes from packing small F_p coefficients that
each occupy a machine word; over F₂ that packing already exists. The measured root
phase confirms it. No tree lane engine was built.

### F4-tracer cost-model check (`diag-costmodel.jsonl`)

KNS (2603.16378) count, per degree, the Macaulay rows (multipliers × generators) and
columns (monomials) and the staircase given by the Hilbert series of a generic system over a field. They do not
treat the Boolean ring, so the comparison uses the Boolean analogue: rows
`Σ_j Σ_{i ≤ D − d_j} C(n, i)`, columns `Σ_{i ≤ D} C(n, i)`, rank = columns − Σ_{i ≤ D} h_i with
`h` the coefficients of `(1 + t)^n / Π_j (1 + t^{d_j})` truncated at the first
non-positive one (semi-regular Boolean series). `n` = occurring variables (35, 44, 47,
58), generators: n quadratic + n cubic (17+17, 23+23, 23+23, 31+31).

| cell | D | generic rows / cols / rank | descended, support-local (production) | descended, occurring multipliers | dense random, same shape |
|:--|--:|--:|--:|--:|--:|
| K1/2^17 | 3 | 629 / 7 176 / 629 | 425 / 2 200 / 424 | 629 / 3 508 / 628 | 629 / 7 176 / 629 |
| K1/2^23 | 3 | 1 058 / 14 235 / 1 058 | 736 / 4 218 / 735 | 1 058 / 6 584 / 1 057 | 1 058 / 14 235 / 1 058 |
| K0/2^23 | 3 | 1 127 / 17 344 / 1 127 | 759 / 5 069 / 758 | 1 127 / 8 157 / 1 126 | 1 127 / 17 344 / 1 127 |
| K0/2^31 | 3 | 1 860 / 32 568 / 1 860 | 1 302 / 9 272 / 1 301 | 1 860 / 14 474 / 1 859 | 1 860 / 32 568 / 1 860 |
| K1/2^17 | 4 | 11 339 / 59 536 / 11 186 | 5 185 / 19 638 / 4 998 | over the F4 caps | over the F4 caps |

Deviation, descended against generic: rows 1.00× with occurring multipliers and
0.67–0.70× with support-local ones; columns 0.44–0.49× and 0.28–0.31×; rank = rows − 1 at
degree 3 in every cell (one syzygy where the generic count predicts none). Degree 2:
rank = rows = n in every cell, as predicted. The dense random system of the same shape
matches the generic count exactly, so the deviation is a property of the descended
system's sparse supports (mean support 25–41 of 35–58 variables), not of the builder.
The tree's per-node step plans are smaller still: mean rank 134–1 036 and 636–4 595
columns per node (`diag-profile.jsonl`).

### Identity

- Per target in every arm and round: hit, exhausted, reductions, splits, propagations,
  infeasible, geometric refutations and witnesses, partial refutations, group additions,
  pair-index builds and lookups are identical across all arms and rounds, in all 8
  cell/engine pairs of both runs.
- The one relation (K1/2^17 target 3, base indices [21, 48, 57]) is identical across
  the candidate binary's arms; the original binary prints verdicts, not indices.

### Wall time

L0 only. The original/previous null control in run 1 ranges 0.24–4.60 per round, so
wall time admits no claim; the per-round graphs are in the run reports. Process CPU time
(user + system, recorded by the candidate binary only) is far tighter and is the time
figure quoted above, still L0.

## Leads

- **Gate-first is exact and halves F6-IC's CPU time on these cells**: a default-on
  candidate, after a check on more cells and on x86 with the eager rewrite.
- **The RREF root (`KIC_F4_INHERIT_ROOT=rref`) beats the default echelon root** on all four
  cells: 1.07–1.48× fewer stage units (f4) and 1.09–1.49× (f6p), CPU previous/rref
  1.16–1.72× (f4) and 1.20–1.78× (f6p) in run 1, 1.12–1.68× and 1.22–1.74× in run 2,
  identical verdicts. The comment in `ReducedBasis::from_system_with` says the fully
  reduced root no longer pays; on this host and these cells it does. It composes with
  gate-first (next section). Not a W5 change; a default-on candidate for a follow-up.
- **Parametric decomposition.** The trees of different targets coincide in shape and
  the target enters the last link linearly; treating the target bits as parameters
  (one comprehensive basis for all targets) is the structural form of "one trace, many
  targets" over F₂. Not attempted.

## Canonical graphs checked and left unchanged

- `docs/ic/progress-timeline.json`
- `docs/index-calculus-scoreboard.html`
- `docs/performance-gains/`

Gate-first is a switch, off by default, measured in the F6-IC decomposition stage only; no
default changed, so none of these graphs changes. No IC-versus-rho statement is made.

## Limits

- L0 only; wall time admits no claim; CPU time is L0 too.
- Four cells, m = 3, `max_degree` 3, standard-subspace bases, aarch64 (lazy rewrite). The
  eager x86 rewrite path is not measured.
- The lane switch batches roots only; the per-instance word count of the lanes charges
  lane words at their bit width, which a SIMD implementation could pack differently, but
  the bits moved per instance cannot fall below the packed row's.
- KNS give their counts for generic systems over a field; the Boolean analogue used here
  is ours, and degree 4 is past the F4 caps on three of the four cells.
