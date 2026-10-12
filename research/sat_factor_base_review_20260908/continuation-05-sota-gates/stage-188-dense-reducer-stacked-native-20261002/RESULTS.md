# Stage 188: dense reducer indexing after dense pair selection

## Decision

`REJECTED_CONFIRMATION`. The exact dense reducer index is not selected and the
runtime candidate was reverted. The repository default remains dense exact
critical-pair selection, deterministic linear symbolic-reducer scan, and
five-column `BlockTables` elimination.

The one-pair screen continued with indexed/linear ratios of `0.755640` wall,
`0.982456` total core-seconds, and `1.019396` peak RSS. The frozen three-pair
confirmation then produced these exact paired ratios:

| pair | order | wall | total core | peak RSS |
|---:|:---|---:|---:|---:|
| 1 | linear, indexed | 0.975365 | 0.966815 | 1.143584 |
| 2 | indexed, linear | 1.076402 | 0.972531 | 0.863879 |
| 3 | linear, indexed | 1.164701 | 1.027522 | 0.979483 |
| **median** | frozen interleaving | **1.076402** | **0.972531** | **0.979483** |

Selection required both median wall and total-core ratios to be strictly below
`0.97`. Wall regressed by 7.64 percent, while CPU fell by only 2.75 percent and
missed the unchanged Stage 178 threshold. No post-hoc threshold or arm was
substituted.

## Mechanism and correctness

The candidate reduced symbolic divisor operations from `4,190,633,182` linear
tests to `103,532,494` exact submask probes, a `97.5294%` reduction, using at
most `1,062,392` charged index bytes per F4 call. This mechanism reduction did
not transfer into an admissible process improvement.

Every screen and confirmation run:

- authenticated source instance
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`;
- used algebraic factor base `span_F2(1,z,...,z^8)` without target-subgroup
  enumeration or known factor-base log labels;
- reproduced equation fingerprint
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`;
- visited 512 X1 masks, skipped 270 non-rational masks, constructed and
  completed 242 systems, found zero roots, and returned exhaustive `UNSAT`;
- agreed exactly on pair selection, field pairs, matrix shapes, basis,
  extraction, degree, logical XORs, actually performed XORs, and all dense-pair
  counters; and
- reported F4 conflicts as `null` rather than manufacturing a SAT statistic.

The confirmation replay passes `58/58` checks. The final native composition
passes `31/31` checks with result SHA-256
`8c3e1d1fc09f341434f71e28ff25b29119b4a0b9854d05f4568992d8a5bbdeb3`.

## Native execution and corrections

All new performance execution, process accounting, hashing, terminal checks,
composition, and verification are Rust-native. The runner meters fresh process
groups with `wait4`, includes all twelve Rayon workers in total CPU, charges
peak RSS, and enforces a 360-second backend watchdog.

The first screen verifier used exact floating-point equality for a recomputed
ratio and retained a `paired assessment replay` failure. An additive verifier
correction permits only `1e-12` relative serialization-scale ratio variation;
it independently passes `22/22` on the unchanged result and raw artifacts. A
later charged final-verifier dry run initially resolved a draft below
`development/` as `development/development`; that exit-2 receipt is also
retained, and the additive root-discovery correction passes its replacement
dry run `31/31`. Neither correction changes a backend run, metric, threshold,
or decision.

## Accounting and boundary

Stage 188 contributes a measured lower bound of 21 components,
`898.725862` wall-seconds, `3,270.706466` total core-seconds, and
`4,965,466,112` bytes peak RSS. The cumulative measured campaign lower bound
becomes 529 components, `22,707.607152` wall-seconds,
`55,657.968684` core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. The initial development compilation,
the failed attempt to invoke a cargo-test-only hashed binary, the one-time
native-meter bootstrap build, and the final writes after charged dry-run
controls are not outer-metered. They are listed rather than assigned modelled
costs.

The same-binary direct-MITM decomposition and full automorphism-aware rho
boundaries from Stage 187 are unchanged. Stage 188 measures one opened
`n=59, ell=9, m=3` decomposition target; it supplies no relation-yield,
unknown-scalar, full-DLP, external-reproduction, novelty, or SOTA evidence.
All-seven and Koblitz index-calculus SOTA remain false.
