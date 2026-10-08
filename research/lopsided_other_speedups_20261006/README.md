# Protocol: which other pipeline stages can the lopsided technique speed up?

Frozen 2026-10-06, before the survey. **Stage diagnostic**: `S`, end-to-end
cost and speedup are **unset** (`null`). This protocol authorises a survey of
candidate transfer targets for the Alman–Vassilevska Williams thin-product
technique (`arXiv:2610.06783v1`; local digest
`sha256:d91a79c5…cc9e76005`, see the
[round-1 protocol](../lopsided_thin_product_20261006/README.md)), plus a
reusable native screening gate. It does not authorise a performance claim, a
new `ecbench` method, or any scoreboard/leaderboard/browser edit.

Conductor note: the repository's `CLAUDE.md` asks for a `conductor check`
before code changes; no `conductor` binary exists in this environment, so no
reservation could be registered. The edited paths are the round-1 module
(owned by the merged PR #1485) and a new study directory; no other in-flight
work touches them as far as `git status` shows.

## Question

Beyond the four IC insertion points of round 1, which other stages of this
repository's ECDLP pipeline — or which other consequences of the paper —
admit a thin-product transfer, and what falsifiable gate decides each one?

## Scope

- Survey the paper's Figure 1 / Sections 5.1–5.4 consequence list and this
  repository's native stages: pair-table/MITM probing (`orbit_pair_table`,
  `native_signed_mitm`), BSGS/kangaroo/rho table lookups, batch MSM relation
  verification, sparse linear algebra, F4/F5 Macaulay reduction, Semaev
  solving vs the paper's 3SUM/3XOR speedups, convolution stages vs
  MonoConvolution, isogeny-route search vs (min,+)/APSP machinery, and the
  paper-excluded problems (SETH, OV, k-SUM/XOR k≥4, unhinted OMv,
  3SUM-indexing, real-valued variants).
- Extend `src/cryptanalysis/lopsided_thin_product.rs` with the
  `WorkloadShape` / `screen_workload` / `Screening` gate and unit tests.
  No pipeline stage is rewired; no algorithm is implemented.

Out of scope: timing anything, adding methods, touching dashboards.

## Method and verdict vocabulary

Each candidate is run through the gate order the code enforces:

1. **Shape**: is the batch's per-output work an inner product over one
   shared middle? Hash probes, group operations, and sparse algebra answer
   no (`NoSharedMiddle`) — the recursion has nothing to prune.
2. **Ring**: is the middle over small integers? Curve points, field
   elements, and residues mod `r` answer no (`WrongRing`).
3. **Size**: with `epsilon = ln D / ln N` and `|W| = N^2/D^kappa`, does
   `N >= D^18` (`epsilon <= 1/18`) and `kappa >= 1/2` hold?

Verdicts: **pass** (advance to a promotion protocol), **negative** (fails a
named gate; recorded with the gate and the numbers), **out of scope**
(paper-excluded or no such stage exists in the ECDLP path).

## Reference and boundaries (AGENTS.md sections 1–2)

Unchanged: generic floor + matched rho in unit `S`. No table row is
produced in this round; the only numbers are gate evaluations on workload
shapes, kept out of every result panel.

## Success and stop conditions

- Success: every surveyed candidate carries a verdict with its failing (or
  passing) gate, the gate is implemented and tested, and a promotion
  protocol is stated for any passer. A repo candidate that passes all three
  gates would earn a follow-up measurement protocol — not a speedup claim.
- Stop: this round ends when the gate, its tests, this protocol, the dated
  report with the ranked table, the funnel diagram plus the
  (epsilon, kappa) screening map, and the PDF are committed and the PR is
  opened. `S` and speedup remain `null` throughout.

## Deliverables

1. `WorkloadShape` / `screen_workload` / `Screening` + 5 new unit tests
   (`cargo test --lib cryptanalysis::lopsided_thin_product` green).
2. This protocol (`README.md`), the dated survey (`REPORT.md`).
3. `figure.mmd` (editable funnel source) + `figure.svg` (rendering).
4. `screen_map.svg` ((epsilon, kappa) screening map, labelled as gate
   evaluation, not measurement).
5. `report.html` + `report.pdf` (report with visuals included).
