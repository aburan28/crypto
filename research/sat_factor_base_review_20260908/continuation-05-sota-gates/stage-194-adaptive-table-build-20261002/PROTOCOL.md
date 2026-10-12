# Stage 194 protocol: adaptive BlockTables construction scheduling

## Hypothesis and claim boundary

Stage 192's selected twelve-worker F4 process spends an aggregate
`142.309949` seconds inside elimination versus `42.797494` seconds building
matrices. The current `BlockTables::add` path enters a Rayon parallel iterator
for every table-construction call, including a call containing one small run.
All 242 independent fixed-X1 systems already execute in one twelve-worker outer
batch, so repeated small nested table builds may pay dispatch and scheduling
cost without useful parallel work.

The candidate preserves every table, pivot, row operation, and XOR. It changes
only the table-construction scheduler:

- use serial iteration unless one call contains at least two table runs and
  schedules at least `65,536` table words;
- retain the current Rayon build at or above that fixed threshold; and
- expose exact add-call, run, scheduled-word, serial-call, and parallel-call
  counters.

This is one opened target's solver-scheduling screen. It cannot establish
relation yield, a full index-calculus DLP, a rho crossover, independent
reproduction, novelty, or SOTA.

## Frozen source, target, and solver

- Branch base and inherited Stage 193 merge:
  `8039a9e68c05362be6d50a894f91c255411f61b4`.
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
  242 rational systems in one outer batch.
- Selected F4: dense exact critical-pair selection, linear symbolic-reducer
  scan, build-serial inner policy, five-column `BlockTables`, and full M4RI
  disabled.

No equation, mask order, batch size, pair rule, reducer choice, matrix shape,
table contents, row reduction, extraction rule, or terminal verifier may
change.

## Native arms and accounting

Use the merged Rust `koblitz_f4_stage188 meter` and a Rust Stage 194 verifier.
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

- Control: `F4_F2_ADAPTIVE_TABLE_BUILD=0`, preserving the selected
  unconditional Rayon iterator.
- Candidate: `F4_F2_ADAPTIVE_TABLE_BUILD=1`, applying the frozen
  two-run/`65,536`-word threshold.

Charge exact builds, tests, both solver processes, all child workers, failures,
composition, wall time, total core-seconds, and peak RSS. Single-core seconds
remain `null`. Preserve every attempted receipt.

## Correctness and identity gates

1. Boolean-F4 and Phase-B backend tests pass in both explicit modes.
2. Every run authenticates the frozen source, visits all 512 masks, skips 270
   non-rational masks, constructs and completes all 242 systems, finds zero
   roots, and returns exhaustive `UNSAT`.
3. Both arms agree exactly on source/equation identities, equations and terms,
   pair and field-pair counters, dense-pair counters, divisor tests, matrix
   shapes, basis, extraction, degrees, logical/performed XORs, peak algorithmic
   memory, total table add calls, total table runs, and scheduled table words.
4. Only the parallel/serial table-add split, phase timings, process wall/CPU,
   and process RSS may differ.
5. Target-subgroup enumeration and discrete-log-label use remain false.
   Conflicts remain `null`.

Any mismatch rejects timing. A timeout or resource refusal is censored.

## Frozen screen and decision

Run one fixed pair in this order:

```text
current, adaptive-table-build
```

Continue to a separate three-pair confirmation stage only if both records are
correct and candidate/control wall and total-core ratios are each strictly
below `0.98`. Otherwise reject and retain the current default. RSS is reported
without a post-hoc gate. The screen never changes the default by itself.

If the screen continues, the follow-on confirmation order is already frozen:

```text
current, candidate, candidate, current, current, candidate
```

and selection requires median paired wall and total-core ratios each strictly
below `0.98`.

## Boundaries

The completed Stage 192 direct-MITM boundary and Stage 193 named-solver controls
remain unchanged. Complete campaign cost remains `null` unless every inherited
and new component is measured. This stage is an engineering diagnostic even if
the candidate passes.

