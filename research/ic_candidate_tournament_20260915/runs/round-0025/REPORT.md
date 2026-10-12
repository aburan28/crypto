# IC candidate tournament: round-0025

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **4656 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 228 fresh inputs over 12 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.378**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,771 | 1 | 2.992 | 3.958e+05 | 24/24 | reference |
| switch | 2,226 | 0.8033 | 2.403 | 5.346e+05 | 24/24 | engineering experiment |
| rho | 926.1 | 0.3342 | 1 | unmeasured | 24/24 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,556 | 1 | 2.37 | 3.651e+05 | 72/72 | reference |
| switch | 2,076 | 0.8122 | 1.924 | 4.987e+05 | 72/72 | engineering experiment |
| rho | 1,079 | 0.422 | 1 | unmeasured | 72/72 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,055 | 1 | 2.646 | 6.248e+05 | 684/684 | reference |
| switch | 2,357 | 0.7713 | 2.04 | 1.031e+06 | 684/684 | engineering experiment |
| rho | 1,155 | 0.378 | 1 | unmeasured | 684/684 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,055 | 1 | 2.646 | 6.248e+05 | 684/684 | reference |
| switch | 2,357 | 0.7713 | 2.04 | 1.031e+06 | 684/684 | engineering experiment |
| rho | 1,155 | 0.378 | 1 | unmeasured | 684/684 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.898 | 1 | reference | 1.472 | 72/72 |
| switch | 1.571 | 0.828 | [0.6036, 1.04] | 1.218 | 72/72 |
| rho | 1.289 | 0.6796 | [0.3967, 1.035] | 1 | 72/72 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 3.061 | 1 | reference | 1.685 | 684/684 |
| switch | 2.388 | 0.7799 | [0.6142, 0.9544] | 1.315 | 684/684 |
| rho | 1.816 | 0.5933 | [0.3527, 0.9211] | 1 | 684/684 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 3.035 | 1 | reference | 1.653 | 684/684 |
| switch | 2.384 | 0.7856 | [0.6204, 0.9621] | 1.299 | 684/684 |
| rho | 1.836 | 0.6049 | [0.3596, 0.9344] | 1 | 684/684 |

Confirmation winner/rho: instructions 2.6455, CI [1.7417759460238174, 4.323392683245986]; native time 1.6854, CI [1.085999923909413, 2.836576757333435].

Replay winner/rho: instructions 2.6456, CI [1.741758964810862, 4.323385095051101]; native time 1.6532, CI [1.07500875602037, 2.785060572632406].

## Factor-base policy boundary

This separately declared panel compares different base policies on identical public ECDLP targets. Support is stable within each arm/case. For each arm, m=3 and B signed points give at most binomial(B+2,3) target images; the rank floor also changes with the column count. Ratios to changed floors are not evidence of an algorithmic advance.

| Stage | Arm | Cell | Signed base size B |
|---|---|---|---:|
| confirmation | incumbent | n13a0 | [182] |
| confirmation | incumbent | n17a1 | [272] |
| confirmation | incumbent | n19a0 | [304] |
| confirmation | incumbent | n19a1 | [304] |
| confirmation | incumbent | n23a0 | [368] |
| confirmation | incumbent | n23a1 | [368] |
| confirmation | incumbent | n29a1 | [464] |
| confirmation | incumbent | n31a0 | [496] |
| confirmation | incumbent | n37a0 | [1184] |
| confirmation | incumbent | n43a1 | [2064] |
| confirmation | incumbent | n59a0 | [2832] |
| confirmation | incumbent | n61a1 | [2928] |
| confirmation | switch | n13a0 | [182] |
| confirmation | switch | n17a1 | [272] |
| confirmation | switch | n19a0 | [304] |
| confirmation | switch | n19a1 | [304] |
| confirmation | switch | n23a0 | [368] |
| confirmation | switch | n23a1 | [368] |
| confirmation | switch | n29a1 | [464] |
| confirmation | switch | n31a0 | [496] |
| confirmation | switch | n37a0 | [148] |
| confirmation | switch | n43a1 | [258] |
| confirmation | switch | n59a0 | [236] |
| confirmation | switch | n61a1 | [244] |
| development | incumbent | n13a0 | [182] |
| development | incumbent | n17a1 | [272] |
| development | incumbent | n19a0 | [304] |
| development | incumbent | n23a0 | [368] |
| development | incumbent | n23a1 | [368] |
| development | incumbent | n31a0 | [496] |
| development | incumbent | n37a0 | [1184] |
| development | incumbent | n43a1 | [2064] |
| development | switch | n13a0 | [182] |
| development | switch | n17a1 | [272] |
| development | switch | n19a0 | [304] |
| development | switch | n23a0 | [368] |
| development | switch | n23a1 | [368] |
| development | switch | n31a0 | [496] |
| development | switch | n37a0 | [148] |
| development | switch | n43a1 | [258] |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.682 | 0.591 | 1 | 36 |
| n17a1 | 0.999 | 0.981 | 1 | 36 |
| n19a0 | 0.823 | 0.809 | 1 | 36 |
| n19a1 | 1.229 | 1.219 | 1 | 36 |
| n23a0 | 0.886 | 0.875 | 1 | 36 |
| n23a1 | 0.868 | 0.865 | 1 | 120 |
| n29a1 | 1.135 | 1.106 | 1 | 36 |
| n31a0 | 0.935 | 0.927 | 1 | 36 |
| n37a0 | 0.925 | 0.915 | 1 | 120 |
| n43a1 | 1.018 | 1.017 | 1 | 120 |
| n59a0 | 1.087 | 1.088 | 1 | 36 |
| n61a1 | 0.857 | 0.857 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
