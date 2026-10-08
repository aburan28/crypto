# Stage 197 protocol: 4,096-row hybrid M4RI screen

## Hypothesis and claim boundary

Stage 196 measured 481 matrices with at least 4,096 rows. They account for
95.5683 percent of selected-F4 profiled elimination time. On those same
matrices, full M4RI used 0.616211 times the profiled elimination nanoseconds and
0.676798 times the performed XORs.

Stage 197 tests the measured selector without profiling overhead:

- control: selected five-column BlockTables on every matrix;
- candidate: existing full M4RI only when a matrix has at least 4,096 rows and
  satisfies the existing column/row shape rule; current BlockTables elsewhere.

This is one opened target's solver-engineering screen. It cannot establish
relation yield, a full DLP, a rho crossover, independent reproduction, novelty,
or SOTA.

## Frozen source, target, and implementation

- Branch base and inherited Stage 196 merge:
  `947ace2f6e10cd9543ac4109983b98802c08e259`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with target-subgroup
  enumeration false and discrete-log-label use false.
- Both arms use fixed-X1 direct symmetrised S4, all 512 masks, all 242 rational
  systems, dense exact pair selection, build-serial inner policy, and one
  twelve-worker outer batch.
- Profiling is disabled in both arms.

The candidate adds only an explicit minimum-row threshold to the existing
full-M4RI shape predicate. Broad `F4_F2_FULL_M4RI=1` without a threshold
retains its historical default minimum of 128 rows.

## Native arms and accounting

Use the merged Rust `koblitz_f4_stage188 meter` and a Rust Stage 197 verifier.
No Python process may execute, verify, compose, or hash new evidence. Both arms
set:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_DISABLE_INNER_BUILD_PARALLEL=1
F4_F2_MATRIX_PROFILE=0
F4_F2_FULL_M4RI_MIN_ROWS=4096
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

- Control: `F4_F2_FULL_M4RI=0`.
- Candidate: `F4_F2_FULL_M4RI=1`.

Charge the fresh build, tests, both solver processes, all workers, memory,
failures, composition, and replay. Report wall, total core-seconds, peak RSS,
conflicts `null`, and single-core `null`.

## Correctness and identity gates

1. Boolean-F4 and Phase-B backend tests pass in both explicit modes.
2. Both runs authenticate the frozen source, visit all 512 masks, skip 270
   non-rational masks, complete all 242 rational systems, find zero roots, and
   return exhaustive `UNSAT`.
3. Source/equation identities, factor-base claims, pair counts, fixed-X1
   counts, degree bound, and exact curve result agree.
4. Control routes zero matrices through full M4RI. Candidate routes exactly 481
   matrices through full M4RI and no matrix with fewer than 4,096 rows.
5. Full M4RI may choose a different valid pivot basis, so downstream reducer,
   logical-XOR, and matrix-row-sum counters may differ. All differences are
   reported.
6. Profiling counters remain zero in both unprofiled runs.

Any timeout, resource refusal, terminal mismatch, or unexpected routing rejects
the timing.

## Frozen screen and decision

Run one fixed pair:

```text
current, hybrid-m4ri-4096
```

Continue to a separately preregistered three-pair confirmation only if both
records are correct and candidate/control wall and total-core ratios are each
strictly below `0.97`. Otherwise reject and retain the current default. RSS is
reported without a post-hoc gate. The screen cannot change the default.

If the screen continues, the follow-on confirmation order is frozen as:

```text
current, candidate, candidate, current, current, candidate
```

and selection requires median paired wall and total-core ratios each strictly
below `0.97`.

## Boundaries

Direct MITM, named SAT controls, unknown-scalar evidence, full rho comparison,
and independent-review gates remain unchanged. Complete campaign cost remains
`null`.

