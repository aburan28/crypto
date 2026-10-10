# Stage 190 protocol: serial inner build with parallel elimination

## Hypothesis

Stage 189 kept the twelve-way outer fixed-X1 batch and disabled every parallel
section inside each F4 call. That reduced total core-seconds by 4.66 percent and
peak RSS by 5.33 percent, but wall time regressed by 0.38 percent and the
candidate was rejected at the frozen joint screen gate.

Elimination dominates the current F4 diagnostic and is the inner work most
likely to benefit wall time from shared row reduction and parallel block-table
construction. Stage 190 therefore keeps inner elimination parallel while
serializing only matrix-build work inside each concurrently executing F4 call:

- monomial-product construction;
- initial symbolic monomial collection;
- frontier reducer lookup/product construction and sorting; and
- final row packing.

The pre-registered hypothesis is that outer fixed-X1 parallelism is sufficient
for the build phases, while retaining inner elimination parallelism recovers
Stage 189's wall loss. This is a one-target scheduling experiment, not a new
algorithm, factor base, relation-yield result, DLP, rho comparison, external
reproduction, novelty result, or SOTA claim.

## Frozen source, target, and solver

- Branch base: `faf73a3ae03f6f1420596a932de3a4b49c9777f9`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Input SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Public source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with no target-subgroup
  enumeration and no known discrete-log labels.
- Selected solver: dense exact critical-pair updates, linear symbolic reducer
  scan, five-column `BlockTables`, full M4RI disabled, all 242 rational systems
  in one twelve-worker outer batch.

No equation, mask order, batch size, pair rule, reducer choice, matrix shape,
elimination policy, extraction rule, or terminal verifier may change.

## Native arms and accounting

Use the merged Rust-native Stage 188 process meter and a Stage 190 Rust phase
verifier. No Python process may execute, verify, compose, or hash new evidence.
Both arms explicitly set:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_FULL_M4RI=0
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

- Control: `F4_F2_DISABLE_INNER_BUILD_PARALLEL=0`.
- Candidate: `F4_F2_DISABLE_INNER_BUILD_PARALLEL=1`.

Unset remains the control unless confirmation selects the candidate. The
candidate must report 242 build-serial F4 calls while retaining the same inner
elimination path and the same outer parallel batch. Charge exact builds, tests,
all child processes and workers, failures, composition, wall time, total
core-seconds, and peak RSS. `single_core_seconds` remains null.

## Correctness and identity gates

1. Boolean-F4 and Phase B backend tests pass in both explicit modes.
2. Every run authenticates the frozen source, visits all 512 masks, skips 270
   non-rational masks, constructs and completes all 242 systems, finds zero
   roots, and returns exhaustive `UNSAT`.
3. Target-subgroup enumeration and log-label use remain false. Conflicts remain
   `null` because F4 has no SAT-conflict counter.
4. Arms agree exactly on equations and terms, pair and field-pair counters,
   dense-pair counters, divisor tests, matrices, basis, extraction, degree,
   logical and performed XORs, table/matrix/scratch memory, and full-M4RI
   counters.
5. Only build-scheduling policy/counts, timing, process RSS, and phase
   nanosecond diagnostics may differ.

Any mismatch rejects timing. Timeout or a resource refusal is censored, not
`UNSAT` and not evidence against the mathematical method.

## Frozen schedule and decision

Run one screen pair:

```text
current, build-serial
```

Continue only if both records are correct and build-serial/current wall and
total-core ratios are below `1.00`.

If the screen continues, run three confirmation pairs:

```text
current, build-serial, build-serial, current, current, build-serial
```

Select build-serial scheduling as the default only if median paired wall and
total-core ratios are both strictly below `0.97`. Report RSS without adding a
post-hoc RSS threshold. A CPU-only or memory-only tradeoff is not a speedup.

## Boundaries

The direct-MITM decomposition and full automorphism-aware rho references
inherited through Stage 189 remain unchanged. Stage 190 cannot move a seven-
gate row from one target's solver scheduling. Complete campaign cost remains
null unless every inherited and new component is measured.
