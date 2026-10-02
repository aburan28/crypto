# Stage 201 protocol: trim full-M4RI combinations to their exact word ends

## Hypothesis and claim boundary

The selected block-8 full-M4RI path gives every combination table entry the
same stride and operational end: the maximum last non-zero word of any pivot
in the block. It therefore builds and applies a combination through that
block-wide end even when the particular subset of pivots represented by the
combination ends earlier. The selected `BlockTables` path already retains a
separate exact end for each pattern.

On the frozen target, selected full M4RI spends 9,222,665,017 performed word
XORs preparing tables and 102,707,985,015 performed word XORs overall. Stage
201 tests whether retaining an exact end per full-M4RI combination avoids
enough zero-suffix writes to improve both process wall time and total CPU.

This is one opened target's solver-engineering screen. It cannot establish
natural relation yield, a full DLP, a rho crossover, independent reproduction,
novelty, or Koblitz index-calculus SOTA.

## Frozen source, target, and implementation

- Branch base and Stage 200 evidence head:
  `0b5febd76096f738dac43420d5d4303cf739525e`.
- Stage 200 result SHA-256:
  `5a1530a9f89645076b99c1c027996793c107f7eacb1d4907ec1cdedb76756f31`.
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
  systems, dense exact pair selection, build-serial inner construction, the
  selected 4,096-row hybrid M4RI policy, serial full-M4RI row clearing, and one
  twelve-worker outer batch.
- Matrix profiling remains disabled.

No equation, mask schedule, critical-pair rule, pivot search, pivot basis,
table width, table combination, row order, extraction rule, or terminal
verification may change.

## Candidate representation and controls

`F4_F2_FULL_M4RI_TRIM_ENDS=1` enables the candidate. Unset and
`F4_F2_FULL_M4RI_TRIM_ENDS=0` preserve the Stage 200 runtime.

For each pivot block, the candidate retains an absolute word end for every
combination mask:

1. combination zero ends at the first word;
2. adding pivot `i` sets the destination end to the maximum of the source
   combination end and pivot `i`'s exact end;
3. table construction writes only through that destination end and never
   reads stale scratch words outside the source end; and
4. a target row XORs only through its selected combination end, then updates
   its row end to the same exact maximum used by the current algebra.

Export exact counts for non-zero trimmed table entries and row lookups, plus
the block-wide table-preparation and row-application words avoided. These are
mechanism counters; actual table preparation remains in
`full_m4ri_table_word_xors` and all actual table/row work remains in
`word_xors_performed`.

## Native arms and accounting

Use the merged Rust native process meter and a Rust Stage 201 verifier. No
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
`F4_F2_FULL_M4RI_MIN_ROWS`; both must exercise the selected Stage 199 hybrid.

- Control: `F4_F2_FULL_M4RI_TRIM_ENDS=0`.
- Candidate: `F4_F2_FULL_M4RI_TRIM_ENDS=1`.

Charge the fresh build, tests, both solver processes, all workers, memory,
failures, composition, and replay. Report wall time, total core-seconds, peak
RSS, conflicts `null`, and single-core time `null`. Preserve every attempted
receipt.

## Correctness and identity gates

1. Differential matrix tests force fixed-end and exact-end schedules on the
   same matrices with varied last non-zero words, partial pivot blocks,
   non-consecutive pivots, and word boundaries. They require identical rank,
   pivot columns, pivot rows, canonical row space, and logical XORs.
2. Candidate table-construction and row-application savings replay exactly to
   the difference in performed XORs from the control algorithm on each test.
3. Boolean-F4 and Phase-B backend tests pass in both explicit modes.
4. Both target runs authenticate the frozen source, visit all 512 masks, skip
   270 non-rational masks, complete all 242 rational systems, find zero roots,
   and return exhaustive `UNSAT`.
5. Both arms route exactly 481 matrices and 351,164 blocks through full M4RI
   and agree exactly on every source, equation, factor-base, pair, field-pair,
   reducer, matrix, basis, extraction, degree, and logical-work counter.
6. Candidate performed XORs and table-preparation XORs must be no larger than
   control, and all exported savings identities must hold. Profiling counters
   remain zero.

Any timeout, resource refusal, terminal mismatch, logical-work mismatch,
savings-identity failure, or unexpected routing rejects timing.

## Frozen screen and decision

Run one fixed pair in this order:

```text
fixed-block-end-control, exact-combination-end-candidate
```

Continue to a separately preregistered three-pair confirmation only if both
records are correct, candidate/control wall and total-core ratios are each
strictly below `0.98`, and candidate performed XORs are strictly below
control. Otherwise reject the candidate and retain the Stage 199 default. RSS
is reported without a post-hoc gate. This screen cannot change the default.

If the screen continues, the follow-on confirmation order is frozen as:

```text
control, candidate, candidate, control, control, candidate
```

and selection requires median paired wall and total-core ratios each strictly
below `0.98` with exact correctness and work identities in every run.

## Boundaries

Direct MITM, named SAT controls, unknown-scalar evidence, full
automorphism-aware rho, external reproduction, and novelty gates remain
unchanged. Complete campaign cost remains `null`.
