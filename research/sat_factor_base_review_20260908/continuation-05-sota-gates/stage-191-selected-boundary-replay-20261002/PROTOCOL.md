# Stage 191 protocol: selected F4 single-core and direct-MITM boundary replay

## Purpose and boundary

Stage 190 selected serial inner build phases with parallel `BlockTables`
elimination and a parallel fixed-X1 outer batch. Stage 191 measures that exact
selected source in strict single-core mode and brackets it with fresh
same-binary direct meet-in-the-middle controls on the identical public target.

This is a decomposition-stage boundary replay. It does not price relation
collection, relation linear algebra, target descent, logarithm recovery, or
full automorphism-aware rho, and it cannot establish a full-method crossover or
SOTA.

## Frozen source and instance

- Main base containing selected Stage 190:
  `6dc7f348cb61abf575d058a66c9b6575e28a907d`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Input SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Public source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, without target-subgroup
  enumeration or known factor-base discrete-log labels.

Build one exact binary from a clean source commit and use it for every arm.
Record source, lock, binary, manifest, protocol, stdout, stderr, and process
receipt hashes.

## Native schedule

Use the Rust-native Stage 188 meter and a Stage 191 Rust verifier. No Python
process may execute, verify, compose, or hash new evidence. Run exactly:

```text
direct MITM before
selected F4 single-core
direct MITM after
```

All arms set every external thread cap to one. The F4 arm additionally sets:

```text
RAYON_NUM_THREADS=1
PQ_F4_X1_BATCH=1
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_FULL_M4RI=0
```

It must not set `F4_F2_DISABLE_INNER_BUILD_PARALLEL`; unset must route through
the selected Stage 190 default. The backend must report
`single_thread_requested=true` and 242 selected build-serial F4 calls.

Direct MITM is the existing deterministic sequential implementation. Its
single-core CPU, wall, and RSS references are the arithmetic means of the two
bracketing fresh processes. F4/reference ratios use those means. Do not choose
the faster direct run or import a prior timing.

## Correctness and accounting

1. Every process authenticates the same source instance and returns exhaustive
   `UNSAT` with no timeout.
2. F4 visits all 512 masks, skips 270 non-rational masks, constructs and
   completes all 242 systems, finds zero roots, reproduces the equation
   fingerprint, and reports no target-subgroup enumeration or log-label use.
3. The F4 command proves strict single-core request through one Rayon worker,
   batch one, and one for every external thread variable.
4. Both direct controls agree on source identity, terminal class, point count,
   pair-table entries, additions, and group result.
5. `conflicts` remains `null` for F4 and direct MITM; do not manufacture a SAT
   conflict count.

Charge the exact build, exact verifier and backend tests, all three fresh
processes, CPU, wall, RSS, verification, composition, failures, and all parallel
resources. Report both `single_core_seconds` and `total_core_seconds`; they are
equal to measured user plus system time for these strict single-thread arms.

## Interpretation

This stage has no selection gate and changes no implementation. It updates the
current single-core and direct-decomposition boundary after Stage 190. Direct
MITM is expected to remain faster. A surprising F4 win would require a repeat
protocol before any claim. Full rho and seven-gate status remain inherited and
unchanged.
