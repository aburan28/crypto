# Lopsided thin products and the IC pipeline: integration report (2026-10-06)

Question: can the Alman–Vassilevska Williams thin-matrix-product technique
(`arXiv:2610.06783v1`) reduce end-to-end index-calculus cost in this
repository? Scope, method, frozen inputs, and stop rules are in the
[protocol](README.md). Status of every statement is labelled inline as
**[proposed]**, **[derived]** (follows from the cited source or the code),
**[verified]** (checked by an independent run in this round), or
**[measured]** (a counted experiment — none in this round).

## Method

1. Downloaded `https://arxiv.org/pdf/2610.06783` (v1, 76 pages;
   `sha256:d91a79…cc9e76005`) **[verified]** and read Sections 1–2, 4 (core
   technique) plus 3, 5.2 (reductions) at survey depth.
2. Implemented the problem interface natively in Rust
   (`src/cryptanalysis/lopsided_thin_product.rs`, 8 unit tests)
   **[verified]** (`cargo test --lib cryptanalysis::lopsided_thin_product`:
   8 passed).
3. Assessed four IC insertion points; implemented none **[derived]**.

## What the paper shows [derived]

- **Thin product (Thm 1 / Thm 5).** `X`: `N x D`, `Y`: `D x N`,
  `D <= N^{1/18}`, `|W| <= N^2/sqrt(D)` wanted positions: all `(XY)[i,j]`
  for `(i,j)` in `W` in `O(N^2/D^0.063)` ops — polynomially less than one
  op per entry of the full product. General form: for every
  `epsilon < 0.1204`, `kappa > 0` there is `gamma > 0` with
  `O(N^2/D^gamma)` when `D <= N^epsilon`, `|W| <= N^2/D^kappa`.
- **Technique.** Schoenhage's 10-multiplication identity inside a
  Coppersmith-style recursion, visiting only the recursion-tree leaves the
  wanted entries need; encoding work is shared across blocks. Prior attempts
  on Coppersmith–Winograd identities failed; the Schoenhage identity's
  sparsity is load-bearing.
- **Graph reading.** With `X`, `Y` as biadjacency matrices this is the
  counting Lopsided All-Edges Sparse Triangle: `(XY)[a,b]` counts middle
  vertices adjacent to both outer vertices `a`, `b`. Cost for
  `epsilon < 1/18`: `O(|W| n^{0.437 epsilon} + n^{2-0.063 epsilon})`.
- **Data structure (Thm 3).** Preprocess `(X,Y)` in `O(N^2/D^0.063)`,
  answer any single entry in `O(D^0.437)` — refutes the thin-hint Online
  Matrix–Vector conjectures.
- **Reductions.** 3SUM / Exact Triangle / APSP (integer deterministic and
  real-valued Las Vegas) reduce to the lopsided problem, giving
  `O(n^1.9992)` 3SUM and `O(n^2.9995)` APSP. SETH / Orthogonal Vectors are
  explicitly unaffected.

## Findings for our IC pipeline

| Insertion point | Implemented? | Blocking difference [derived] | Required experiment [proposed] |
|---|---|---|---|
| Relation-verification batch | No | PDP systems are polynomial (Semaev/Weil), not inner products; needs a new reduction | Matched baseline/candidate full-DLP suites, all phases in `S` |
| Factor-base intersection query | No — closest fit | Offline `W` vs online queries mirrors relation queries, but IC yield is algebraic, not set-theoretic | One-target workloads paired vs matched rho, five online phases named |
| Sparse-LA block | No | Block Wiedemann needs large sparse matvecs over `Z/rZ`, not dense small-integer thin products | Price matrix build + final LA inside total `S` |
| Target-descent batch | No | Descent is per-target recursion; batching changes the workload to multi-target | Single-target first, then a separate `k`-target question |

Positive **[verified]**: the interface (spec validation, regime gates
`N >= D^18` / `|W| <= N^2/sqrt(D)`, `epsilon`/`kappa`, checked naive
reference, quoted-bound hooks, Boolean triangle-count equivalence) is
implemented and tested. Negative **[verified]**: no insertion point is
wired into any candidate; `S`, phase costs, and speedup are `null`; no
`ecbench` session was run; no curve was named or registered. The
`paper_*_bound` / `graph_time_bound` functions quote the paper's
asymptotics — they are not measurements and do not establish that any IC
phase moves.

## Boundary, table, ratio (AGENTS.md sections 1–8)

- Boundaries stated before measuring: generic floor + matched rho, both
  unchanged; the unit is `S` (§2 of the protocol).
- No table row is produced: a row without a verified answer is not a result.
  The only numbers in this round are unit-test fixtures and quoted-bound
  evaluations, kept out of every result panel.
- Falsification target for a follow-up: full-pipeline `S` below the matched
  baseline with the floor ratio falling (advance); `S` down with the floor
  ratio flat is engineering; a stage count down with `S` up is relabelling;
  a re-derived number with no algorithm change is accounting.
- All phases unpriced (`null`); wall time appears nowhere.

## Graphs checked — no change [verified by inspection]

- `docs/index-calculus-scoreboard.html` + `docs/ic/progress-timeline.json`
  (embedded copy + every-point table): no new IC/rho ratio, so no point
  added; embedded copy untouched.
- `docs/ic-leaderboard.html` / `LEADERBOARD.md` / `leaderboard.json`: no
  whole-pipeline measurement, no repriced curve, no registry change.
- `docs/browser/data.json`: no new curve, session, round, or IC1 identity.
- `docs/curves/registry.json`: no new curve named, so no ICV1 registration
  owed; no retired spellings introduced (checked: this report and the new
  module name curves only by family, e.g. "binary Koblitz toy curves").
- `docs/performance-gains/`, study-local figures: none affected.
- The explanatory diagram and cost-model graph below are new visuals for
  this search round (research-visuals skill), not canonical-graph updates.

## Open checks [proposed]

1. A follow-up protocol implementing the pruned-recursion algorithm natively
   and replaying it against `wanted_product_i64` on frozen inputs.
2. Whether factor-base intersection queries on any concrete toy base admit a
   thin-product encoding that preserves decomposition probability — needs a
   reduction argument first, code second.
3. Re-reading the paper's Section 4 parameter trade-offs (`epsilon < 0.1204`
   general regime) once a concrete IC shape (`N`, `D`, `|W|`) is proposed.

## Sources

- Paper PDF: `https://arxiv.org/pdf/2610.06783` (v1). Abstract page:
  `https://arxiv.org/abs/2610.06783`.
- Code: `src/cryptanalysis/lopsided_thin_product.rs` (this repo, this PR).
- Visuals: `figure.mmd` / `figure.svg` (pipeline diagram),
  `cost_model.svg` (quoted-bound illustration), all embedded in
  `report.html` / `report.pdf`.
