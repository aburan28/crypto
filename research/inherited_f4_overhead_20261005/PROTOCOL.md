# Inherited-F4 overhead round: protocol

Status: **preregistered**, written before any confirmatory run. The
exploratory profiling that motivated it is disclosed below with its
numbers.

## Question

The inherited-F4 decomposition engine
(`src/cryptanalysis/inherited_f4.rs`) dominates the relation phase of the
`pdp3-koblitz` index-calculus pipeline on the `E_0` ladder of
[`research/f6_ic_ecbench_ladder_20261005`](../f6_ic_ecbench_ladder_20261005/RESULT.md).
Its counted unit is 64-bit word XORs. On one frozen m = 19 run the whole
process executes about **120 instructions per counted word XOR**, so nearly
all of its machine work lies outside the unit. Can that overhead be cut
without changing a single decision, count or answer?

## What changes, and what must not

The candidate changes only bookkeeping in `inherited_f4.rs`:

1. **Completion rows** (`CompletionBatch`). Each product `t·p` of a
   generator whose degree dropped is gathered once. Each occurrence is hashed
   to a batch-local id, the row's parity is taken on one bit per id, and
   whether a monomial is in the layout is asked once per distinct monomial.
   Previously each occurrence was sorted for parity, hashed for layout
   membership, sorted again for deduplication and hashed again for packing,
   and the column index was rebuilt after the layout grew.
2. **Layout extension** (`insert_columns`). The new columns are keyed once
   (`sort_masks_descending`) and merged into the already-sorted layout. The
   old → new map falls out of the merge. Previously the whole layout was
   re-sorted with a comparator and the map was rebuilt through a hash index.
3. **Pivot materialisation** (`ensure_current`). The is-it-current test is
   inlined; only the rewrite stays out of line.
4. **Pivot XOR in `insert`** (`xor_words`). Fixed-width chunks, which the
   compiler packs into SSE2 vector XORs on the baseline x86-64 target.
5. **Closure products** (`close`, off by default) use the same hash parity
   (`odd_terms`) as item 1 in place of a sort.

Nothing algebraic changes: same rows, same row order, same layout, same
pivots, same word-XOR charges. The gate below enforces this run by run.

## Boundaries (AGENTS.md §1)

- **Floor**: `√(π/4m)` in `S`, unchanged from the ladder: 0.2458 (m = 13),
  0.2033 (m = 19), 0.1848 (m = 23).
- **Reference**: same-target strong rho, from the ladder's frozen sessions,
  unchanged.
- In the ecbench unit (`S`, word XORs unpriced), this round **cannot**
  move `S`: the counted work is identical by construction. So its class is
  **engineering** at best, never advance (§3). What it can move is what the
  unit does not see: machine instructions and time per counted operation.
  The ecbench framework's measured solver price (`ns_per_op`) and the
  pinned calibration ratios in `docs/ic/calibration.json` are built from
  that.

## Frozen inputs

Every input is an exact measured-child input exported by
`ecbench profile-input` from a frozen spec. Algorithm seeds come from the
spec, and nothing is chosen after a run.

| set | spec | sequences | role |
| --- | --- | --- | --- |
| L13 | `../f6_ic_ecbench_ladder_20261005/SPEC-n13.json` | round 1, all 16 IC runs (seq 24–46, both arms) | ladder, m = 13 |
| L19 | `../f6_ic_ecbench_ladder_20261005/SPEC-n19.json` | round 1, all 16 IC runs | ladder, m = 19 |
| L23 | `../f6_ic_ecbench_ladder_20261005/SPEC-n23.json` | seq 24, 25, 27, 28 (two workloads, both arms) | ladder, m = 23 |
| H13 | `SPEC-holdout-n13.json` (target seed 202610052113) | all 8 | fresh holdout, m = 13 |
| H19 | `SPEC-holdout-n19.json` (target seed 202610052119) | all 4 | fresh holdout, m = 19 |

Exploratory profiling used only L13 seq 25, 28, 31 and L19 seq 25. Those
four inputs stay in the suite and are reported, but the confirmatory
reading is the rest.

## Binaries

- **Baseline**: `ecbench` built from `main` at
  `fa80835af5f7c18bd400ccb2002f5e6bd8f6c058`, unmodified, release profile.
- **Candidate**: the same, plus this branch's `inherited_f4.rs`, release
  profile.

Both are built on the measuring host with the same `rustc`. Their SHA-256
digests are recorded in the result.

## Measurements

1. **Gate: identical outputs (every input in every set).** Each binary runs
   `ecbench exec` on the input. The child output is compared after removing
   only timing and host-placement fields: `wall_ns`, `solve_wall_ns`,
   `cpu_ns`, `cpu_at_*`, `run_delay_ns`, `timeslices`, `phases_ns`,
   `placement`, `hw_solve`, the wall-priced solver `gae`/`ns_per_op`, and
   `framework_total_gae_before_solver_removal` (see `strip.jq`). Everything
   else must be byte-identical after `jq -S`: recovered scalar,
   verification, every phase's counted operations, `solver_ops`, every
   solver statistic and `S`.
2. **Primary: instructions.** Valgrind Callgrind total `Ir` for the whole
   process, one run per (binary, input). This is deterministic, so
   contention does not matter. Inputs: L13 (16), H13 (8), L19 seq 24, 25,
   27, 28, and H19 (4). L23 is too slow under Callgrind (about 2 hours a
   run) and is gated and timed only.
3. **Secondary: time.** Process CPU time (`cpu_ns` from the child output)
   and solve wall, pinned with `taskset -c 3`. Five rounds, alternating
   baseline and candidate (ABAB). An A/A run of the baseline against a
   byte-identical copy of itself measures the spread. Inputs: L19 seq 25
   (F4), L19 seq 24 (F6), L23 seq 25 (F4). This cloud VM cannot reserve
   a core against host neighbours (§10), so wall is reported with its
   A/A spread as a practicality note only.

## Decision rule

- **Abort** if any input fails the gate. The change is then a different
  algorithm, not a speedup, and is withdrawn.
- **Success** (engineering class): the gate passes on all inputs, and
  - the candidate/baseline `Ir` ratio summed over the confirmatory m = 19
    inputs (L19 seq 24, 27, 28 and H19) is at most **0.90**, and
  - no single `Ir` input regresses by more than **2%**.
- **Null**: the gate passes, but the ratio exceeds 0.90 or an input regresses
  by more than 2%. Report it as such; nothing is merged as a speedup.
- Wall time never decides. It is reported against the A/A spread.

The 0.90 threshold is set knowing the exploratory figure on L19 seq 25 (the
final exploratory candidate: 0.832). It is preregistered against the other
inputs, not that one.

## Cost accounting

`Ir` counts the whole `ecbench exec` process: factor base, oracle setup,
relation search (the decomposition engine), linear algebra, verification,
and the harness's own I/O. That is a whole-pipeline, one-target, cold
measurement of machine work. It is not the ecbench `S`, which is
unchanged by construction.

## Exploratory record (before this protocol)

On L19 seq 25, Callgrind total `Ir`, every candidate output identical to
the baseline's:

| step | change | `Ir` | vs baseline |
| --- | --- | --: | --: |
| baseline | `main` fa80835a | 65,787,175,131 | 1.000 |
| c1 | hash parity for completion rows; merge in `extend_columns` | 65,347,502,767 | 0.993 |
| c2 | + inlined `ensure_current`; keyed sort of the new columns | 63,046,006,372 | 0.958 |
| c3 | + generation-stamped parity table; dedup before keying | 62,751,751,886 | 0.954 |
| c4 | + hash deduplication of new columns | 59,739,505,267 | 0.908 |
| c5 | + chunked pivot XOR; *sink-word scatter in `rewrite`* | 59,715,385,382 | 0.908 |
| c6 | sink-word scatter **reverted** (+2.5 G in `rewrite`); `CompletionBatch` | 54,710,207,601 | 0.832 |

Rejected: the sink-word scatter in `rewrite`. It replaced the
deleted-column branch with a clamped store, and cost 2.5 G more instructions
than it saved, because the branch skips the deleted bits cheaply.

## Commands

`run.sh` in this directory is the thin orchestration: build, export
inputs, gate, Callgrind and ABAB timing. It calls only `ecbench`,
`valgrind`, `taskset` and `jq`.
