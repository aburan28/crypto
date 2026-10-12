# N53 fixed-useful-base S3 root-batch experiment: frozen protocol

**Question.** On the `a=0` Koblitz curve `EC1N53Ce0hb097de99be9a`
(`GF(2)[x]/(x^53+x^6+x^2+x+1)`, prime subgroup order
`21044858204113`), does batching inversions in the *second* S3 root call
reduce charged rank and target work while preserving the exact useful rows?
This follows the [lead assessment](../conversation_leads_assessment_20261011/REPORT.md)
and the [pair-symmetry result](../n53_pair_symmetry_20261010/RESULT.md).
The earlier [S3 batching panel](../compact_s3_batch_20260929/RESULT.md)
already measured a complete-process engineering gain for 1,024 targets at
K440/K550. This protocol tests a different K220, one-target cold and online
regime with the newer phase-exact runner; batching itself is not a new idea.

## Boundary, reference, and estimands

The generic-group floor is `S_floor = sqrt(pi/(2*106))` for the 106-element
signed-Frobenius/negation orbit, where `S = total group-addition equivalents /
sqrt(r)`. The reference is signed-Frobenius **strong rho** on each identical
public Q, with rung 3, 32 lanes, four distinguished bits, including its own
setup and verification. The operation-normalized full `S` for this compact
native IC example is not calibrated; it stays unknown. The primary exact
comparison is control versus candidate under the same fixed 23,320-point,
220-useful-column base (BLAKE3
`7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973`),
rank stream, and target. Report S3 calls, regular/degenerate cases, field
multiplications and inversions, novel rank per accepted row, memory bytes,
all exclusive cold and target-online phases, and separately rho's full
process and online intervals. The shared Apple host admits descriptive
timings only; stage ratios are not complete-method speed ratios.

## Hypothesis and one variable

The control makes one Itoh-Tsujii inversion per regular S3 root call. The
candidate gathers the regular roots of one index state across all 53
Frobenius shifts, performs four-lane Montgomery batch inversion, and consumes
the resulting roots in exactly the original scan order. The algorithmic
variable is `KIC_S3_ROOT_BATCH=0|1` in the standalone native example
`examples/koblitz_s3_root_batch_experiment.rs`. `KIC_PAIR_SYMMETRY=0`,
`KIC_ROOT_BLOOM=0`, `KIC_RANK_PROBE_CAP` absent, and all other algorithmic
settings are fixed.
Extra roots computed before an early exit are charged. No target scalar is
passed to either IC or rho process.

The primary gate requires all nine paired cells to reach rank 220/220, replay
the same 220 independent rows and per-attempt relation/probe sequence,
reproduce the same root table and verified target scalar, and have zero bad
relations. The operation gate is at least 50% fewer rank-stage inversion
primitives in every candidate pair, with full field multiplication counts
reported so inversion work cannot be relabelled as free multiplication.
Descriptive progression to an isolated-host comparison requires median
paired candidate/control complete-cold time at most 0.90 and target-online
time at most 1.10, without exclusions. A miss on either gate is retained as
an outcome. No automatic `vs_rho` ledger promotion is in this protocol.

## Frozen inputs and schedule

Use native `koblitz_rho_fixture` with commands
`53 0 signed_frobenius 1 strong 700003 hash:53281307`,
`... hash:53291307`, and `... hash:53301307`, in that order. Keep complete
generator output as verifier-only fixtures and give each measured process
only a point-only JSONL. Reject and retain an attempt if a Q duplicates a
committed N53 point or another Q here. The deterministic public rule,
generator source and binary hashes, output SHA-256, and independent native
group replay must be committed before timing. This is three new disjoint
public targets; no point may be replaced after a timing result is seen.

For each Q, use rank seeds `532101`, `532102`, `532103` and rho seeds
`531201`, `531202`, `531203`. Run one control, candidate, and rho cell per
seed, interleaving order by target and seed so each arm occupies each ordinal
position equally often. Nine pairs and nine rho cells are required. Each
cell is a fresh single-worker process; IC gets a 90-second wall cap and
16-GiB observed RSS cap, rho gets the same envelope. Preserve every intent,
stdout, stderr, exit/timeout/OOM status, elapsed wall, peak RSS, rank trace,
base dump, and exact command. Do not overwrite a run directory. Freeze and
push source, input, binary, workload, and schedule hashes before the first
timed cell.

## Accounting, replay, and stop conditions

Charge base construction, normal basis/self-test, full index, every failed
rank attempt, matrix/LA, target query/PDP/relation/descent/check, and all
candidate batch preparation. Verify exclusive phase sums. Independently
replay each rank relation, rank transition, factor log and recovered scalar
with native Rust arithmetic; source-process `group_verified` is an additional
check. Retain incomplete and failed cells. Stop a cell at the frozen envelope,
then run the remaining cells; do not extend a cap or replace a target or seed
because of the outcome. If the native verifier or runner is unavailable, mark
that requirement unverified and preserve the raw process record.

The result report will carry one table with control, candidate, and rho in
their declared units, a source-linked graph and PDF, all nine paired ratios
with uncertainty, and a decision against the gates above. An observed
reduction in inversions alone is a stage diagnostic. A full method claim
requires calibrated common `S`, independent replay and an admitted isolated
same-point reference under `docs/ic/boundary_targets.json`.
