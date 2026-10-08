# Stage 192 protocol: selected parallel F4 and direct-MITM boundary replay

## Purpose

Stage 191 established the current strict single-core selected-F4/direct-MITM
boundary. Stage 192 measures the selected twelve-worker outer-batch path against
fresh same-binary direct controls on the identical opened target. This closes
the current parallel wall, total-core and RSS decomposition comparison without
borrowing Stage 190's direct timing from another run.

This is a decomposition-stage boundary. It does not price natural relation
yield, relation collection, relation linear algebra, target descent, recovery,
or full automorphism-aware rho and cannot establish a full-method crossover.

## Frozen source and instance

- Main base: `56bdc4c51cae7f718ff0953b6e7e7e8b2baa8dc3`.
- Manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Public source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with no target-subgroup
  enumeration or known discrete-log labels.

Build one exact binary and use it for every arm.

## Native schedule

Use the Rust-native process meter and Stage 192 verifier. No Python process may
execute, hash, compose, or verify new evidence. Run exactly:

```text
direct MITM before
selected F4 parallel
direct MITM after
```

All arms set external thread variables to one. The F4 arm additionally sets
`RAYON_NUM_THREADS=12`, `PQ_F4_X1_BATCH=512`,
`F4_F2_DENSE_PAIR_SELECT=1`, and `F4_F2_FULL_M4RI=0`. It must not set
`F4_F2_DISABLE_INNER_BUILD_PARALLEL`; unset must select the Stage 190 default.
The backend must report `single_thread_requested=false` and 242 selected
build-serial calls.

The direct reference is the arithmetic mean of both fresh controls. Do not
choose the faster direct wall time or import a prior result. Report direct
single-core seconds as measured user plus system time; selected parallel F4
single-core seconds remain `null`.

## Correctness and accounting

1. All processes authenticate the same source and return exhaustive `UNSAT`.
2. F4 visits all 512 masks, skips 270 non-rational masks, completes all 242
   systems, finds zero roots, and reproduces the equation fingerprint.
3. Both direct controls agree on factor points, pair entries, additions, and
   group result.
4. Factor-base subgroup enumeration and log-label use remain false. Conflicts
   remain `null` or unavailable because neither method exposes a SAT-conflict
   counter.
5. Charge exact build/tests, all three processes and twelve F4 workers,
   verification, composition, wall, total core-seconds, and RSS.

## Interpretation

This stage has no implementation selection. It updates the parallel
decomposition boundary after Stage 190. A surprising F4 win requires a repeat
protocol before any claim. Full rho and seven-gate status remain unchanged.
