# Full-image setup-reuse gate and execution

The protocol below preceded the executed [60-target study](../different_j_reuse_20261015/RESULTS.md).
Its V1 target-order leak was corrected before measurement, and the V2
wrapper failure was preserved. The frozen
[one-target run](reports/run-20261015T1916Z/REPORT.md) changed exact
support membership by four subgroup points relative to each comparator,
while image construction took a median 98.532% of full-image CPU. The
question is whether charging that construction once across genuinely
distinct public targets changes the total method cost.

- **Frozen group and target capacity:** keep the same `F_(67²)` source and
  destination, explicit degree-three map, generator of order 61, and
  factor bases from [PUBLIC_V2.json](PUBLIC_V2.json). Use all 60 distinct
  nonidentity source subgroup targets in a public order derived from a
  committed seed and checked by a separate scalar-index census. Compare
  prefixes `N=1,4,16,60`; no repeated target may count toward reuse.
- **Arms and accounting:** full-image setup once plus all target relation
  queries, projected-image setup once plus all targets, destination-direct
  base once plus all targets, direct source-group BSGS with its own cold
  baby table and lookups, direct source negation rho with all retries,
  and A/A controls. Include input/model validation, memory, image
  enumeration/maps, base projection, pair-table construction, every
  failed query, rank/linear solve, recovery and final scalar checks in
  complete cold totals. Report one-time setup, total and per-target
  amortized values separately; compare each prefix against the cheaper
  applicable direct control.
- **Correctness and stop:** preregister source/input/binary hashes,
  target order, random streams, caps and rotation before timing. Require
  rank and scalar recovery for every target, independent replay of every
  scalar, and exact target distinctness. Preserve cap exhaustion,
  memory breaches, failed runs, stdout and stderr. Three L0 rounds can
  screen cost; portable counted `S`, isolated interleaving and a second
  host remain required before a frontier claim. A setup-amortization
  effect without a lower cold total than direct BSGS/rho is reported as
  a phase effect.

The native Rust public-target generator and batch runner committed their
exact inputs before execution. Three V3 and three V4 rounds covered all
60 distinct targets in seven arms; each version passed 1,260/1,260 separate
scalar replays and 10,960 relation-witness checks. The V4 whole-process
child CPU paired full-image/BSGS ratio median was 2.483 at `N=60` on a
shared host. The next larger-order and direct χ-coordinate gates are
specified in the [execution study](../different_j_reuse_20261015/NEXT_GATE.md).
