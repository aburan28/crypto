# Optimized producer admission, before improvement rounds

This is correctness, instrumentation and provenance qualification. It selects
no winner and consumes none of the three improvement rounds. Do not infer an
algorithmic speedup from these integration vectors.

## Sources and hypothesis

Restore the sealed round 0020 `both` and round 0023 `scaled` source inventories.
Reconstruct `pairinv` using the already archived round 0024 patch over `scaled`.
`prepare.py` checks the archive and every input source hash before making a new
directory, checks all patch hunks before applying them, and emits the derived
source manifest. Original archives and earlier directories remain unchanged.
The build passes this manifest digest as `IC_SOURCE_MANIFEST_SHA256`; both library
and worker use a tracked compile-time environment dependency. Their tags must
agree, and admission checks the emitted tag against the sealed source. Build
directories are isolated. This closes a locally reproduced Cargo-cache failure
where a successful build invocation reused an executable from another snapshot.
The source-hash discrepancy in the old winner summary is recorded in the goal
checkpoint; the restored candidate registry and source bytes are authoritative.

The derived source adds exclusive phase markers and query/rank diagnostics,
explicit `tiny_gauss` admission, a nonempty IC descent requirement, and rejection
sampling of nonzero scalars. Its source ID is different from the original.
It must preserve independently correct arithmetic, ordinary-query relations,
factor logs and recovered DLPs. Requiring IC witnesses and removing modulo bias
can change a trace; never claim byte equivalence to all historical traces.

## Fixed integration vectors and limits

- `scaled` and `pairinv`: `13a0,23a1,37a0,43a1,61a1`.
- `both`: `13a0,17a1,23a1`, within its original degree limit.
- One public-hash target per cell, target seed `2026092526`, algorithm seed
  `2026092527`; subgroup-orbit recipe seed 43 and requested bound `6*n`.
- Three fresh native/profile process pairs per cell. These are repetitions of
  fixed correctness vectors, not fresh yield samples or held-out confirmation.
- Configuration: pair table, explicitly `tiny_gauss`, three summands, batch one,
  maximum 65,536 collection/descent trials. Rust 1.94.1, archived Cargo.lock,
  original static musl build flags, Linux amd64, Valgrind 3.22.0.
- Each child: one pinned allowed CPU, 8 GiB address-space cap, 180-second timeout.
  Each source job stops at 45 minutes. Retain every failed process and receipt.

Generate the public fixture and independently check its field/subgroup. Use the
new `inventory` mode to freeze the exact usable base and full candidate manifest
before executing measured IC jobs. Every subsequent process must agree with the
admitted base and method. Native/profile answers and rank evidence must agree.
All reports and profiles are frozen artifacts. No retries to obtain a pass.
Record the actual field kernel and reject any native/profile dispatch change.

## Accounting

The process is profiled from entry through termination, including every marker
and diagnostic. `setup` includes runtime initialization, input/curve/target setup,
driver administration, certificate serialization and process cleanup.
`factor_base` ends after the usable orbit columns and points have been built.
`precompute` builds negative-point views, lambda powers, normal-basis coordinates,
pair tables and fixed-base tables. `queries` includes initial walk scalars/points
and every walk advance, including identity queries which are charged and skipped.

`pdp` scans for three-point decompositions, including unsuccessful calls and
collection bookkeeping. `relation_check` re-adds proposed witnesses, deduplicates
relations and certifies factor logs. `matrix_build` forms each scalar-field row
and retains its rank event; `relation_la` contains incremental elimination and
back substitution. `target_descent` contains all nested target queries, PDP,
recombination and internal replay exactly once. `recovery_check` is the final
independent library scalar multiplication for every output target. There is no
isogeny transport: its cost is explicitly zero.

Client dumps close the previous phase. Program termination belongs to `setup`
because every normal report/error path first returns to that phase. The parser
rejects mixed historical/scientific labels, missing boundaries, duplicate profile
parts and nonclosure against Callgrind's whole-process checksum. A killed or
incomplete worker retains raw evidence but cannot become a complete DLP result.
Native time is measured outside the whole process; phase wall times are unknown
in this adapter and are never substituted for instruction counts.

The Python audit independently replays the public scalar walk, skips only its
actual identity queries, reconstructs the scalar-field matrix, checks every rank
event and derives its SHA-256 digest. The solver need not hash the same matrix
again in its hot path. No-witness from the partial pair table is unresolved,
never UNSAT. Query uncertainty must cluster by independently seeded workload;
three executions of one deterministic walk do not triple its statistical sample.

## Admission and remaining gates

Pass requires complete certificates for every scheduled target, exact identities,
exclusive instruction closure and matching query/matrix diagnostics. A run record
always has `promotion_eligible=false` here. Failures remain in the integration
summary and fail CI; the source is not declared qualified after a partial pass.

Before comparative tuning, separately measure the effect of instrumentation and
admission changes against the original sources, qualify a strong matched rho
width/source, and seal the five-or-more-cell confirmation protocol and familywise
rule. Original combined phase intervals remain unpriced under schema v2, so an
observer-control ratio is a diagnostic, not an end-to-end speedup claim.
