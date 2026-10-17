# p61 selected-class census and adjacency completion

- **Status:** `investigating`
- **Kind:** empirical bottleneck, implementation gap, and resource limit
- **Affected work:** [p61 census report](../research/p61_isogeny_class_20261014/reports/run-20261014T0902Z/REPORT.md), [streaming protocol](../docs/p61-isogeny-class-streaming.md), and [external audit protocol](../docs/p61-isogeny-audit.md)
- **Next action owner:** Conductor task `T-2735`
- **Last updated:** 2026-10-14

## Problem

The complete selected ordinary isogeny-class graph rooted at `p61-i0025` over
`p = 2^61 - 1` cannot yet be delivered because the 358,094,304-vertex census
and independent replay require a multi-week local campaign, the complete
degree-5/13/37/59 adjacency join is not yet implemented and verified at that
scale, and the current UGREEN volume has insufficient operational margin for
that graph plus the requested diverse-field corpus and redundant working data.

The selected class has trace `-704859453`, fundamental discriminant
`-8726545188371316595`, conductor `f_pi = 1`, and an unconditionally certified
class number of `358094304` with cyclic factors `[179047152, 2]`. Completion
means more than generating a prefix: every traversal column must be generated,
independently replayed, admitted to an exact ledger, globally deduplicated, and
then joined into the complete selected-degree adjacency relation. Enumerating
this one class also does not enumerate every curve or every isogeny class over
the field.

## Evidence and uncertainty

- **Observed:** the flag-0 PARI certificate returned `1`; the initial sealed
  prefix through `x = 320` contains 5,760 unique invariants with zero duplicate
  records. The immutable receipts and report are linked above. A 16,384-column
  continuation through `x = 16704` generated and independently replayed 294,912
  additional records under the same hash-chained protocol. The expanded ledger
  contains 300,672 records; the external audit found 300,672 unique invariants
  and zero duplicates, with `complete_class_coverage = false`.
- **Derived or proved:** the certified class order fixes the required vertex
  count at 358,094,304. With 18 lanes, the traversal requires 19,894,128
  columns. Degree 5/13/37/59 adjacency requires 2,506,660,128 directed records
  under the frozen schema.
- **Suspected explanation:** no mathematical obstruction is suspected. The
  remaining gap is operational and implementational: sustained generation,
  independent replay, global audit, adjacency construction, and storage.
- **Missing capability:** execution and monitoring of the restart-safe runner
  over all remaining segments; an admitted full-class ledger; a scale-tested
  external adjacency join and verifier; and separately provisioned storage for
  the broader multi-field corpus.
- **What this evidence does not establish:** prefix uniqueness does not imply
  global uniqueness or complete coverage, a runtime extrapolation is not a
  completed run, and the selected class is not the field-wide curve population.
  A literal list of one 64-bit j-invariant for each possible field value alone
  would require about 18.45 EB, before graph structure or labels.

## What would unblock it

- [ ] Generate and independently replay all 19,894,128 traversal columns with
  immutable segment/checkpoint/receipt chains and no unclassified gaps.
- [ ] Seal a ledger containing exactly 358,094,304 admitted vertices, externally
  sort every invariant, and obtain zero duplicates with
  `complete_class_coverage = true`.
- [ ] Implement, test, execute, and independently verify the complete
  degree-5/13/37/59 adjacency join containing exactly 2,506,660,128 directed
  records under the frozen representation.
- [ ] Provision enough durable space for raw data, external-sort/join working
  sets, manifests, and at least one independent copy before expanding to the
  requested multi-field corpus. The currently observed 7.1 TiB free is not a
  petabyte allocation.

The first three conditions are jointly necessary and sufficient for the
selected-class graph under the current definition. They do not complete a
field-wide enumeration, a diverse-field training corpus, or the M7
geometry-specificity validation.

## Next discriminating step

- **Action or proof obligation:** use the tested fail-closed continuation runner
  to complete and independently replay the next contiguous segment from the
  certified `x = 16704` checkpoint, then periodically rebuild the exact ledger
  and uniqueness audit.
- **Inputs and controls:** use the committed `c9bee1b2` producer, committed
  `d06472954` replay binary, distinct audit settings `(1, 0)` and `(2, 7)`, the
  certified preflight/certificate pair, create-only outputs, and digest replay
  after a host sync barrier.
- **Execution limits and stop condition:** stop on any source mismatch,
  checkpoint discontinuity, digest mismatch, replay mismatch, duplicate,
  storage failure, or nonzero process exit; preserve partial artifacts as
  incomplete evidence.
- **Outcomes and next actions:** an admitted segment extends only the certified
  prefix and informs throughput; a duplicate or replay mismatch blocks further
  admission pending diagnosis; an infrastructure failure is retained and
  retried only after the resource condition changes.
- **Recheck trigger:** each newly admitted segment, adjacency implementation,
  or storage-tier change.

A planned step is not a completed experiment. Apply the existing execution,
review, budget, and archival rules before acting on it.

## Resolution and history

- **2026-10-14:** opened after unconditional class-number certification and a
  successful 5,760-vertex prefix audit made the remaining completion gap
  concrete. Unresolved.
- **2026-10-14:** admitted the contiguous `x = 320..16704` continuation after
  independent replay and expanded the globally deduplicated prefix to 300,672
  vertices (0.0839645% of the certified class). Unresolved.
- **2026-10-14:** added and unit-tested the create-only continuation supervisor;
  this implements restart-safe segment orchestration but does not count any
  unexecuted planned segment as data. Unresolved.
