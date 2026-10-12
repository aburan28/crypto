# Stage 199: selected 4,096-row hybrid M4RI default

## Decision

`SELECTED_FOR_REPOSITORY_DEFAULT`. Unset Boolean F4 now uses full M4RI on
eligible matrices with at least 4,096 rows and five-column `BlockTables`
elsewhere.

The controls remain explicit:

- `F4_F2_FULL_M4RI=0` restores the prior all-BlockTables path;
- `F4_F2_FULL_M4RI=1` preserves broad full M4RI from 128 rows; and
- `F4_F2_FULL_M4RI_MIN_ROWS` explicitly overrides the applicable minimum.

## Exact unset replay

The selected-source command contains no full-M4RI mode or minimum-row
assignment. It reports:

| field | value |
|:---|---:|
| wall seconds | 25.936151 |
| total core-seconds | 153.825561 |
| single-core seconds | null |
| peak RSS | 3,270,033,408 B |
| conflicts | null |
| full-M4RI matrices | 481 |
| full-M4RI blocks | 351,164 |
| performed XORs | 102,707,985,015 |

The replay authenticates source instance
`954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`,
the exact equation fingerprint, all 512 masks, 270 non-rational skips, all 242
completed rational systems, zero roots, and exhaustive `UNSAT`. The
factor-base contract retains target-subgroup enumeration false and known-log
use false. Profiling is disabled.

Against the inherited Stage 192 direct-MITM mean, the selected replay consumes
`57.022674x` wall, `502.014614x` CPU, and `79.627768x` RSS. The F4
engineering improvement therefore does not approach the direct decomposition
boundary.

## Verification corrections

The first unset launch stopped before spawning the backend because the layout
directory was also used as the meter output directory. The corrected nested
run path preserves the unchanged command.

The first charged finalizer then rejected the legitimate 30-second layout
receipt because generic receipt replay incorrectly enforced the solver's
360-second watchdog and successful terminal. The correction moved those
requirements into solver-command validation and retained generic failed/setup
receipts for additive accounting. Both failed attempts are preserved and
charged; the selected solver measurement is unchanged.

The final Rust verifier binds the Stage 198 result SHA-256
`38dd39b13cbcf26ac0fd6e5f43f27bf6e733488d349cc24c197c4842047c5198`,
the selected and finalizer commits, policy source markers, unset command,
terminal/work counters, every receipt, and all claim flags. Final replay passes
`17/17`; result SHA-256 is
`2f4551a1644b5fc49d15a8be07a13ff4a5b56219595d496ef9e12b2077902edd`.

Stage 199 contributes a measured lower bound of 14 components,
`205.639384` wall-seconds, `839.421520` total core-seconds, and
`5,293,916,160` bytes peak RSS. The cumulative measured campaign lower bound
is 696 components, `25,527.893759` wall-seconds, `68,523.607063`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. This is a selected one-target F4
engineering improvement, not relation-yield, unknown-scalar, full-rho,
independent-review, novelty, or Koblitz index-calculus SOTA evidence.
