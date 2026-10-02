# Stage 199 protocol: select the 4,096-row hybrid M4RI default

## Purpose and claim boundary

Stage 198 confirmed the opt-in 4,096-row hybrid in three interleaved pairs with
median candidate/current ratios 0.908433 wall, 0.879279 total core, and
0.630062 RSS. Stage 199 changes only the default policy and verifies an unset
selected-source replay.

This remains one opened target's F4 engineering. It does not establish natural
relation yield, a full DLP, a rho crossover, independent reproduction, novelty,
or SOTA.

## Frozen source and inherited evidence

- Branch base and inherited Stage 198 merge:
  `011e255fff3b13c9062c978b7554d64a48fd6860`.
- Stage 198 result SHA-256:
  `38dd39b13cbcf26ac0fd6e5f43f27bf6e733488d349cc24c197c4842047c5198`.
- Input manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with no target-subgroup
  enumeration or known discrete-log labels.

## Selected policy and controls

After the selected-source change:

- unset `F4_F2_FULL_M4RI` and unset `F4_F2_FULL_M4RI_MIN_ROWS` select the
  4,096-row hybrid;
- `F4_F2_FULL_M4RI=0` disables full M4RI and preserves the prior selected
  five-column BlockTables control;
- `F4_F2_FULL_M4RI=1` with no explicit minimum preserves the historical broad
  full-M4RI control with minimum 128 rows; and
- an explicit `F4_F2_FULL_M4RI_MIN_ROWS` overrides the policy minimum while
  retaining the selected enable/disable mode.

No equation, threshold, elimination implementation, pair policy, factor base,
outer schedule, or verification rule changes.

## Tests and exact replay

Use the merged Rust native process meter and a Rust Stage 199 verifier. No
Python process may execute, verify, compose, or hash new evidence.

Run exact tests for:

1. pure policy/default/control resolution;
2. Boolean-F4 table/full equivalence;
3. Phase-B backend behavior with explicit disabled, unset selected, and
   explicit broad-full modes; and
4. Stage 199 verifier replay.

Build the exact selected source in a fresh target directory. Then run one exact
unset backend process with:

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

The command must contain no `F4_F2_FULL_M4RI` or
`F4_F2_FULL_M4RI_MIN_ROWS` assignment.

## Selection gates

The unset replay must:

- authenticate the exact source/equation identities and factor-base contract;
- visit all 512 masks and complete all 242 rational systems;
- return exhaustive `UNSAT`;
- report the selected hybrid policy and minimum 4,096 rows;
- route exactly 481 matrices and 351,164 blocks through full M4RI;
- perform exactly 102,707,985,015 XORs; and
- keep every profile counter zero.

If all gates pass, record `SELECTED_FOR_REPOSITORY_DEFAULT`. Otherwise reject
the default change and restore the prior runtime.

Charge the fresh build, all tests, unset replay, every worker, memory, failures,
composition, and replay. Report wall, total core-seconds, peak RSS, conflicts
`null`, and single-core `null`.

## Boundaries

Direct MITM remains the same-target decomposition reference and full
automorphism-aware rho remains the full-method boundary. Complete campaign cost
remains `null`; no seven-gate or SOTA claim changes.

