# Protocol: lopsided thin products meets the IC pipeline (interface + assessment)

Frozen 2026-10-06, before any implementation beyond the reference interface.
**Stage diagnostic**: `S`, end-to-end cost and speedup are **unset** (`null`).
This protocol authorises literature review, a native reference interface, and
a falsifiable insertion assessment only. It does not authorise a performance
claim, a new `ecbench` method, or any scoreboard/leaderboard/browser edit.

Paper: Alman–Vassilevska Williams, *Truly Subquadratic 3SUM and Truly Subcubic
APSP via Triangles in Sparse Lopsided Graphs*, `arXiv:2610.06783v1`
(5 Oct 2026). Local copy digest (downloaded 2026-10-06):
`sha256:d91a79c5dbf1fc662f0c5c1bb2b690ac3d83d9ef5bf872bde4cdaf7cc9e76005`
(916436 bytes). Section/Theorem numbers below refer to that version.

## Question

Can the paper's thin-product / lopsided-triangle technique reduce the
end-to-end index-calculus cost `S = total operations / sqrt(r)` on this
repository's binary Koblitz toy curves, against the matched rho reference?

## Scope

- Read the paper's Sections 1–2 and 4 (thin-product theorem, Schoenhage /
  Coppersmith construction, data-structure version) and Sections 3 and 5.2
  (reductions) at survey depth.
- Implement the problem **interface** in native Rust
  (`src/cryptanalysis/lopsided_thin_product.rs`): exact spec, regime
  predicates, checked naive reference, quoted-bound hooks, lopsided-graph
  reading, IC insertion assessment. The paper's pruned-recursion algorithm
  itself is out of scope.
- Assess four IC insertion points; wire none into `S` accounting.

Out of scope: implementing the paper's algorithm, timing anything as a
speedup, adding an `ecbench` method, touching dashboards.

## Frozen inputs

- Paper version above; no other external input.
- Code interface inputs are constructed in unit tests (tiny `N`, `D`, `W`
  with hand-checked inner products and triangle counts). No curve, factor
  base, target, or seed is introduced, so no ICV1 registration is owed.

## Reference and boundaries (AGENTS.md sections 1–2)

- Floor: generic-group bound `S_floor`, unchanged.
- Reference: matched Pollard rho on the same instances in the same unit,
  unchanged. No rho run is performed in this round.
- A future candidate would report one table, one unit (`S`), every variant
  as a row including the reference and unmodified baseline, with ratio
  columns against each boundary and a correctness column.

## Success and stop conditions

- Success (advance): a later round, under its own protocol, shows a
  full-pipeline `S` with the ratio to the floor falling vs the matched
  baseline, every phase priced, answers verified — then the verdict,
  scoreboard, leaderboard, and browser update in that PR.
- Engineering / relabelling / accounting outcomes are classified per
  AGENTS.md section 3 at that time; a stage-count reduction that leaves
  total `S` flat or worse is relabelling, not a speedup.
- Stop: this round ends when the interface module, its `cargo test` suite,
  this protocol, the dated report, the diagram plus cost-model graph, and
  the PDF are committed and the PR is opened. Open checks are recorded in
  the report, not resolved here.

## Cost accounting

No operations are counted in this round. `S`, phase costs, conversion
factors, and `speedup = baseline_total_operations /
candidate_total_operations` are all `null`. The module's `paper_*` functions
evaluate the paper's **quoted** asymptotic bounds as `f64` hooks; they are
not measurements and never enter a result table.

## Deliverables

1. `src/cryptanalysis/lopsided_thin_product.rs` + `mod.rs` registration,
   `cargo test --lib cryptanalysis::lopsided_thin_product` green.
2. This protocol (`README.md`), the dated report (`REPORT.md`).
3. `figure.mmd` (editable diagram source) + `figure.svg` (rendering).
4. `cost_model.svg` (quoted-bound illustration, labelled as such).
5. `report.html` + `report.pdf` (report with visuals included).
