# Stage 200 protocol: size-gated parallel full-M4RI row clearing

## Hypothesis and claim boundary

Stage 199 selected full-matrix block-8 M4RI for eligible matrices with at least
4,096 rows. Its exact unset replay routes 481 matrices and 351,164 pivot
blocks through that path, but every pivot block still clears its remaining
rows serially. The process consumes 153.825561 total core-seconds in 25.936151
wall-seconds with a twelve-worker Rayon pool, an average of 5.93 charged cores.

Rows below one M4RI pivot block are independent: each reads the same immutable
combination table, clears the same pivot columns, and cannot affect any other
row until the block has finished. The repository's older Boolean-F4 M4RI
kernel already uses this property for large target-row regions. Stage 200 tests
the same scheduling mechanism inside the selected current F4 implementation.

This is one opened target's solver-engineering screen. It cannot establish
natural relation yield, a full DLP, a rho crossover, independent reproduction,
novelty, or Koblitz index-calculus SOTA.

## Frozen source, target, and implementation

- Branch base and Stage 199 evidence head:
  `421b1e646b15b27720ae015944a24655bc1db163`.
- Stage 199 result SHA-256:
  `2f4551a1644b5fc49d15a8be07a13ff4a5b56219595d496ef9e12b2077902edd`.
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
- Both arms use the selected fixed-X1 direct symmetrised S4 formulation, all
  512 masks, all 242 rational systems, dense exact pair selection,
  build-serial inner construction, the selected 4,096-row hybrid M4RI policy,
  and one twelve-worker outer batch.
- Matrix profiling remains disabled.

The candidate changes only how a sufficiently large block applies its already
built table to rows below the pivot block. Pivot search, pivot rows, table
construction, table contents, row order, counters, matrix construction,
critical-pair rules, extraction, and terminal verification remain unchanged.
It uses the existing Rayon pool and creates no additional worker pool or OS
threads.

## Candidate scheduler and controls

`F4_F2_FULL_M4RI_PARALLEL=1` enables the candidate. A block is shared across
the existing Rayon pool only when

```text
remaining target rows * table suffix words >= 1,048,576
```

and the pool has more than one worker. Smaller blocks retain the exact serial
loop. `F4_F2_FULL_M4RI_PARALLEL=0` and unset preserve the Stage 199 scheduler
during this screen.

Export and charge exact counts for candidate parallel blocks, target rows, and
scheduled target-row words. These counters must be zero in the control and
positive in the candidate. They describe scheduling only and do not replace
performed-XOR accounting.

## Native arms and accounting

Use the merged Rust native process meter and a Rust Stage 200 verifier. No
Python process may execute, verify, compose, or hash new evidence. Both arms
set:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_DISABLE_INNER_BUILD_PARALLEL=1
F4_F2_MATRIX_PROFILE=0
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

Neither arm may set `F4_F2_FULL_M4RI` or
`F4_F2_FULL_M4RI_MIN_ROWS`; both must exercise the selected Stage 199 default.

- Control: `F4_F2_FULL_M4RI_PARALLEL=0`.
- Candidate: `F4_F2_FULL_M4RI_PARALLEL=1`.

Charge the fresh build, tests, both solver processes, all workers, memory,
failures, composition, and replay. Report wall time, total core-seconds, peak
RSS, conflicts `null`, and single-core time `null`. Preserve every attempted
receipt.

## Correctness and identity gates

1. Differential matrix tests force serial and parallel table application and
   compare rank, pivot columns, canonical row space, logical XORs, performed
   XORs, table-preparation XORs, blocks, and pivots across varied shapes,
   partial blocks, non-consecutive pivots, and word boundaries.
2. Boolean-F4 and Phase-B backend tests pass in both explicit scheduling
   modes.
3. Both target runs authenticate the frozen source, visit all 512 masks, skip
   270 non-rational masks, complete all 242 rational systems, find zero roots,
   and return exhaustive `UNSAT`.
4. Both arms route exactly 481 matrices and 351,164 blocks through full M4RI,
   perform exactly 102,707,985,015 XORs, and agree exactly on every source,
   equation, factor-base, pair, field-pair, reducer, matrix, basis, extraction,
   degree, logical-work, performed-work, and algorithmic-memory counter.
5. Only candidate scheduling counters, phase timing, process wall/CPU, and
   process RSS may differ. Profiling counters remain zero.

Any timeout, resource refusal, terminal mismatch, work mismatch, or unexpected
routing rejects timing.

## Frozen screen and decision

Run one fixed pair in this order:

```text
serial-clear-control, parallel-clear-candidate
```

Continue to a separately preregistered three-pair confirmation only if both
records are correct and candidate/control wall and total-core ratios are each
strictly below `0.98`. Otherwise reject the candidate and retain the Stage 199
default. RSS is reported without a post-hoc gate. This screen cannot change the
default.

If the screen continues, the follow-on confirmation order is frozen as:

```text
control, candidate, candidate, control, control, candidate
```

and selection requires median paired wall and total-core ratios each strictly
below `0.98`.

## Boundaries

Direct MITM, named SAT controls, unknown-scalar evidence, full
automorphism-aware rho, external reproduction, and novelty gates remain
unchanged. Complete campaign cost remains `null`.
