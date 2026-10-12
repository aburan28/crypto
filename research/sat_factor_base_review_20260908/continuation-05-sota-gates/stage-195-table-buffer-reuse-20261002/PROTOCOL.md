# Stage 195 protocol: reuse retired BlockTables word buffers

## Hypothesis and claim boundary

Stage 194 measured 728,503 block-table runs across 5,500 table-add calls on the
selected single target. The current implementation allocates a fresh
`Vec<u64>` for every run. When later pivots enter a table group,
`BlockTables::extend` first retires that group's old tables and drops their
word buffers, then allocates replacement buffers for the same group.

The candidate reuses retired word-buffer capacity within that same
`extend` call. It sorts retired buffers and replacement requirements by
capacity, assigns the largest adequate buffer to the largest requirement,
zeros it, and builds the unchanged table through the current Rayon schedule.
Unmatched requirements allocate normally; unused retired buffers are dropped
before construction. No buffer is retained across add calls.

This targets allocator traffic and contention only. Tables, pivots, row
operations, scheduling, and logical/performed XORs must remain identical. This
is one opened target's solver-engineering screen, not relation-yield, a full
DLP, a rho crossover, independent reproduction, novelty, or SOTA evidence.

## Frozen source, target, and solver

- Branch base and inherited Stage 194 merge:
  `0074963f7bea462634954ba440ae3e11431bdf7d`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Public source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with target-subgroup
  enumeration false and discrete-log-label use false.
- Selected formulation: fixed-X1 direct symmetrised S4, all 512 masks and all
  242 rational systems in one twelve-worker outer batch.
- Selected F4: dense exact critical-pair selection, linear symbolic-reducer
  scan, build-serial inner policy, five-column `BlockTables`, current
  unconditional parallel table construction, and full M4RI disabled.

No equation, mask order, batch size, pair rule, reducer choice, matrix shape,
table contents, construction order, row reduction, extraction rule, or terminal
verifier may change.

## Native arms and accounting

Use the merged Rust `koblitz_f4_stage188 meter` and a Rust Stage 195 verifier.
No Python process may execute, verify, compose, or hash new evidence. Both arms
set:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_DISABLE_INNER_BUILD_PARALLEL=1
F4_F2_FULL_M4RI=0
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

- Control: `F4_F2_REUSE_TABLE_BUFFERS=0`, one fresh word-buffer allocation
  per table run.
- Candidate: `F4_F2_REUSE_TABLE_BUFFERS=1`, same-add retired-buffer reuse.

Charge exact builds, tests, both solver processes, all workers, failures,
composition, wall time, total core-seconds, and peak RSS. Single-core seconds
remain `null`. Preserve every attempted receipt.

## Correctness and identity gates

1. Boolean-F4 and Phase-B backend tests pass in both explicit modes.
2. Every run authenticates the frozen source, visits all 512 masks, skips 270
   non-rational masks, constructs and completes all 242 systems, finds zero
   roots, and returns exhaustive `UNSAT`.
3. Arms agree exactly on source/equation identities, equations and terms, pair
   and field-pair counters, dense-pair counters, divisor tests, matrix shapes,
   basis, extraction, degrees, logical/performed XORs, table add calls, table
   runs, and scheduled table words.
4. Control reports zero reused buffers and one allocation per run. Candidate
   reports at least one reuse; allocations plus reuses equal total runs.
5. Only buffer allocation/reuse counters, table-allocation peak, phase timings,
   process wall/CPU, and process RSS may differ.
6. Target-subgroup enumeration and discrete-log-label use remain false.
   Conflicts remain `null`.

Any mismatch rejects timing. A timeout or resource refusal is censored.

## Frozen screen and decision

Run one fixed pair:

```text
current, reuse-table-buffers
```

Continue to a separate three-pair confirmation stage only if both records are
correct and candidate/control wall and total-core ratios are each strictly
below `0.98`. Otherwise reject and retain the current default. RSS is reported
without a post-hoc gate. The screen never changes the default by itself.

If the screen continues, the follow-on confirmation order is frozen as:

```text
current, candidate, candidate, current, current, candidate
```

and selection requires median paired wall and total-core ratios each strictly
below `0.98`.

## Boundaries

The completed Stage 192 direct-MITM boundary, Stage 193 named-solver controls,
and Stage 194 rejected scheduling screen remain unchanged. Complete campaign
cost remains `null` unless every inherited and new component is measured.

