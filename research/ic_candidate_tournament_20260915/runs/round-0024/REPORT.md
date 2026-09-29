# IC candidate tournament: round-0024

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **4656 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 228 fresh inputs over 12 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.7249**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,790 | 1 | 1.282 | 3.985e+05 | 24/24 | reference |
| counted | 3,883 | 1.392 | 1.784 | 2.595e+06 | 24/24 | engineering experiment |
| rho | 2,176 | 0.7801 | 1 | unmeasured | 24/24 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,765 | 1 | 1.253 | 3.95e+05 | 72/72 | reference |
| counted | 3,499 | 1.265 | 1.585 | 2.338e+06 | 72/72 | engineering experiment |
| rho | 2,207 | 0.7982 | 1 | unmeasured | 72/72 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,087 | 1 | 1.379 | 6.312e+05 | 684/684 | reference |
| counted | 3,908 | 1.266 | 1.746 | 4.261e+06 | 684/684 | engineering experiment |
| rho | 2,238 | 0.7249 | 1 | unmeasured | 684/684 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,087 | 1 | 1.379 | 6.312e+05 | 684/684 | reference |
| counted | 3,908 | 1.266 | 1.746 | 4.261e+06 | 684/684 | engineering experiment |
| rho | 2,238 | 0.7249 | 1 | unmeasured | 684/684 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.838 | 1 | reference | 1.358 | 72/72 |
| counted | 1.541 | 0.8383 | [0.5763, 1.136] | 1.138 | 72/72 |
| rho | 1.354 | 0.7364 | [0.4028, 1.154] | 1 | 72/72 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.81 | 1 | reference | 1.488 | 684/684 |
| counted | 2.44 | 0.8684 | [0.6422, 1.126] | 1.292 | 684/684 |
| rho | 1.888 | 0.6721 | [0.4135, 1.008] | 1 | 684/684 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.803 | 1 | reference | 1.485 | 684/684 |
| counted | 2.424 | 0.8647 | [0.647, 1.106] | 1.284 | 684/684 |
| rho | 1.887 | 0.6733 | [0.4147, 1.009] | 1 | 684/684 |

Confirmation winner/rho: instructions 1.3795, CI [0.8471434631826612, 2.3866160181875498]; native time 1.4879, CI [0.9984759976253803, 2.41995961692937].

Replay winner/rho: instructions 1.3795, CI [0.8471393475952328, 2.3865904755501104]; native time 1.4853, CI [0.9918418570807597, 2.4127800581063856].

## Factor-base policy boundary

This separately declared panel compares different base policies on identical public ECDLP targets. Support is stable within each arm/case. For each arm, m=3 and B signed points give at most binomial(B+2,3) target images; the rank floor also changes with the column count. Ratios to changed floors are not evidence of an algorithmic advance.

| Stage | Arm | Cell | Signed base size B |
|---|---|---|---:|
| confirmation | counted | n13a0 | [52] |
| confirmation | counted | n17a1 | [68] |
| confirmation | counted | n19a0 | [76] |
| confirmation | counted | n19a1 | [76] |
| confirmation | counted | n23a0 | [92] |
| confirmation | counted | n23a1 | [92] |
| confirmation | counted | n29a1 | [116] |
| confirmation | counted | n31a0 | [124] |
| confirmation | counted | n37a0 | [148] |
| confirmation | counted | n43a1 | [258] |
| confirmation | counted | n59a0 | [236] |
| confirmation | counted | n61a1 | [244] |
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
| development | counted | n13a0 | [52] |
| development | counted | n17a1 | [68] |
| development | counted | n19a0 | [76] |
| development | counted | n23a0 | [92] |
| development | counted | n23a1 | [92] |
| development | counted | n31a0 | [124] |
| development | counted | n37a0 | [148] |
| development | counted | n43a1 | [258] |
| development | incumbent | n13a0 | [182] |
| development | incumbent | n17a1 | [272] |
| development | incumbent | n19a0 | [304] |
| development | incumbent | n23a0 | [368] |
| development | incumbent | n23a1 | [368] |
| development | incumbent | n31a0 | [496] |
| development | incumbent | n37a0 | [1184] |
| development | incumbent | n43a1 | [2064] |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.864 | 0.773 | 1 | 36 |
| n17a1 | 1.226 | 1.208 | 1 | 36 |
| n19a0 | 0.897 | 0.884 | 1 | 36 |
| n19a1 | 0.970 | 0.960 | 1 | 36 |
| n23a0 | 0.951 | 0.953 | 1 | 36 |
| n23a1 | 0.944 | 0.946 | 1 | 120 |
| n29a1 | 0.914 | 0.885 | 1 | 36 |
| n31a0 | 0.843 | 0.833 | 1 | 36 |
| n37a0 | 1.210 | 1.210 | 1 | 120 |
| n43a1 | 0.896 | 0.896 | 1 | 120 |
| n59a0 | 0.894 | 0.894 | 1 | 36 |
| n61a1 | 0.957 | 0.958 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
