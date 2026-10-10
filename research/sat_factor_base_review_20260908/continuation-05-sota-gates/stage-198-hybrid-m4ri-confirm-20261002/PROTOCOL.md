# Stage 198 protocol: 4,096-row hybrid M4RI confirmation

## Purpose and claim boundary

Stage 197's unprofiled one-pair screen admitted the opt-in hybrid at
candidate/current ratios 0.764699 wall, 0.884128 total core, and 0.831216 RSS.
Stage 198 runs the already frozen three-pair confirmation before any default
change.

This remains one opened target's solver-engineering comparison. It is not
relation-yield evidence, a full DLP, a rho crossover, independent reproduction,
novelty, or SOTA.

## Frozen source, target, and implementation

- Branch base and inherited Stage 197 merge:
  `18ee85c1075010ce83be60a4a02636e9aa4ac83a`.
- Inherited Stage 197 result SHA-256:
  `84b92c9d0e2b028095b48553ce33d7f07077d1590259dfd831422acc9e8de241`.
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
- Control: selected five-column BlockTables on every matrix.
- Candidate: existing full M4RI only for matrices with at least 4,096 rows and
  the existing column/row shape gate; BlockTables elsewhere.
- Profiling disabled.

No equation, factor base, mask order, outer batch, pair policy, table width,
full-M4RI implementation, threshold, or verification rule may change.

## Native execution and accounting

Use the merged Rust `koblitz_f4_stage188 meter` and a Rust Stage 198 verifier.
No Python process may execute, verify, compose, or hash new evidence. Every run
sets:

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

Run exactly:

```text
current, candidate, candidate, current, current, candidate
```

Charge fresh build/test, all six solver processes and workers, memory, failures,
composition, and replay. Report pairwise and median wall, total core-seconds,
and RSS ratios. Single-core and conflicts remain `null`.

## Correctness gates

Every run must authenticate the exact source and equations, visit all 512 masks,
skip 270 non-rational masks, complete all 242 rational systems, find zero roots,
and return exhaustive `UNSAT`. Target-subgroup enumeration and known-log use
remain false.

Each control must route zero full-M4RI matrices. Each candidate must route
exactly 481 matrices, and profile counters must remain zero. Pair/fixed-X1
counts, degree, source, target, and curve result agree. Pivot-basis-dependent
reducer, logical-XOR, matrix-row-sum, and performed-XOR counters are reported
but may differ between arms.

Any timeout, resource refusal, terminal mismatch, or unexpected routing rejects
the candidate.

## Frozen decision

Confirm the candidate only if:

1. all six runs pass every correctness and identity check;
2. median paired candidate/control wall ratio is strictly below `0.97`; and
3. median paired candidate/control total-core ratio is strictly below `0.97`.

RSS has no post-hoc gate. On success, record
`CONFIRMED_FOR_DEFAULT_CHANGE`; a separate selected-source commit and unset
default replay are still required before changing the repository default.
Otherwise record `REJECTED_CONFIRMATION` and retain the current default.

## Boundaries

Direct MITM, named SAT controls, unknown-scalar evidence, full rho comparison,
and independent-review gates remain unchanged. Complete campaign cost remains
`null`.

