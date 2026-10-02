# Stage 193: current native WDSat and CryptoMiniSat comparators

## Result

Neither named SAT solver decided the exact public `n=59, ell=9, m=3`
single-target decomposition instance within the frozen limit. WDSat reached
the 120-second outer watchdog without a terminal. CryptoMiniSat 5.14.7 reached
its internal 120-second limit, exited cleanly, and reported
`INDETERMINATE`.

| arm | terminal | wall seconds | total core-seconds | single-core seconds | peak RSS | conflicts |
|:---|:---|---:|---:|---:|---:|---:|
| direct MITM, inherited Stage 192 mean | exhaustive `UNSAT` | 0.454839 | 0.306417 | 0.306417 | 41,066,496 B | null |
| selected parallel F4, inherited Stage 192 | exhaustive `UNSAT` | 16.765098 | 181.297986 | null | 4,366,270,464 B | null |
| **fresh WDSat** | timeout, inconclusive | **120.274299** | **118.710537** | **118.710537** | **13,910,016 B** | **null** |
| **fresh CryptoMiniSat 5.14.7** | `INDETERMINATE` | **121.436959** | **121.201177** | **121.201177** | **280,756,224 B** | **3,473,638** |
| licensed Magma F4 | unavailable | null | null | null | null | null |
| GGMP construction | no same-instance executable construction | null | null | null | null | null |

Against the same-target Stage 192 direct-MITM mean, WDSat consumes
`264.432530x` wall and `387.415616x` CPU, while CryptoMiniSat consumes
`266.988730x` wall and `395.543899x` CPU. The RSS ratios are
`0.338719x` and `6.836625x`, respectively. These ratios measure consumed
work before inconclusive solver stops; they are not ratios between completed
algorithms.

## Provenance and interpretation

The WDSat run does not reuse an opaque stale executable. The native meter
creates a detached worktree at WDSat commit
`55d55b2620d768d9f7c78dcd8990a0689533c1d0`, authenticates all 16 sealed
source/header files, installs the exact capacity header with SHA-256
`4a73c3b5a14ded98f749d282b597594ae3ac05355a39cadcb9345a90bf735352`,
and performs separate clean and build steps. The rebuilt 114,608-byte binary
has SHA-256
`1de1e3328af3f00963cd5339c07cc4ab18efd55d1749c953ecaafcefc54fd185`,
exactly matching the historical sealed compatible build.

CryptoMiniSat is the current local Homebrew 5.14.7 executable, SHA-256
`a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac`.
Its final statistics provide the conflict count that the earlier
externally-killed run could not report. `INDETERMINATE` is censored and is
not interpreted as `UNSAT`. WDSat emitted neither a terminal nor a conflict
count before its watchdog, so its conflict field remains `null`.

Both inputs authenticate source-instance id
`954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
The factor base is defined algebraically as
`span_F2(1,z,...,z^8)`; the manifest and verifier both require target-
subgroup enumeration and discrete-log-label use to be false.

## Verification and accounting

The Rust verifier rejects altered source files, tool binaries, version output,
capacity configuration, inputs, commands, thread caps, watchdogs, terminal
semantics, receipt hashes, resources, factor-base contract, inherited Stage
192 result, or claim flags. A relative-path replay exposed a verifier path
comparison defect; the first result is preserved, the discovered root is now
canonicalized, and both relative and absolute final replays pass `21/21`.
The finalizer also replays every one of the 32 meter receipts and requires an
exact one-to-one match with the metrics charged by the stage. The final result
SHA-256
is
`d119bc42060ea55f6bc2a91739f86ecebb8bc09914c6f5c3fe9e88cca9c51ab5`.

Stage 193 contributes a measured lower bound of 32 components,
`252.416574` wall-seconds, `251.207316` total core-seconds, and
`280,756,224` bytes peak RSS. The cumulative measured campaign lower bound is
614 components, `24,163.163195` wall-seconds, `61,810.903842`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. The current same-instance named-solver
comparison is stronger, but licensed Magma, a same-instance executable GGMP
construction, natural independent-relation yield, relation collection, linear
algebra, unknown-scalar recovery, full automorphism-aware rho comparison,
independent reproduction, and novelty review remain open. This is negative
decomposition-stage evidence, not a Koblitz index-calculus SOTA result.
