# IC candidate tournament: round-0018b

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **2340 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 96 fresh inputs over 8 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.405**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,010 | 1 | 0.6854 | 1.906e+05 | 18/18 | reference |
| block | 2,989 | 0.9929 | 0.6805 | 1.892e+05 | 18/18 | engineering experiment |
| both | 2,941 | 0.9769 | 0.6695 | 1.862e+05 | 18/18 | engineering experiment |
| column | 2,963 | 0.9842 | 0.6745 | 1.875e+05 | 18/18 | engineering experiment |
| rho | 4,393 | 1.459 | 1 | unmeasured | 18/18 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,851 | 1 | 0.6957 | 1.804e+05 | 54/54 | reference |
| block | 2,824 | 0.9906 | 0.6892 | 1.788e+05 | 54/54 | engineering experiment |
| both | 2,774 | 0.9732 | 0.677 | 1.756e+05 | 54/54 | engineering experiment |
| column | 2,801 | 0.9827 | 0.6837 | 1.773e+05 | 54/54 | engineering experiment |
| rho | 4,098 | 1.437 | 1 | unmeasured | 54/54 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,191 | 1 | 0.7115 | 1.808e+05 | 288/288 | reference |
| both | 3,080 | 0.9652 | 0.6868 | 1.745e+05 | 288/288 | engineering experiment |
| rho | 4,484 | 1.405 | 1 | unmeasured | 288/288 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,191 | 1 | 0.7115 | 1.808e+05 | 288/288 | reference |
| both | 3,080 | 0.9652 | 0.6868 | 1.745e+05 | 288/288 | engineering experiment |
| rho | 4,484 | 1.405 | 1 | unmeasured | 288/288 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 0.9935 | 1 | reference | 0.9093 | 54/54 |
| block | 0.9754 | 0.9817 | [0.9559, 1.008] | 0.8927 | 54/54 |
| both | 0.9464 | 0.9525 | [0.9184, 0.9754] | 0.8661 | 54/54 |
| column | 0.9532 | 0.9594 | [0.9208, 0.9939] | 0.8724 | 54/54 |
| rho | 1.093 | 1.1 | [1.062, 1.136] | 1 | 54/54 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 0.9825 | 1 | reference | 0.9181 | 288/288 |
| both | 0.9445 | 0.9613 | [0.9288, 0.9895] | 0.8825 | 288/288 |
| rho | 1.07 | 1.089 | [1.036, 1.142] | 1 | 288/288 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 0.9698 | 1 | reference | 0.9058 | 288/288 |
| both | 0.9445 | 0.9739 | [0.9475, 0.9973] | 0.8821 | 288/288 |
| rho | 1.071 | 1.104 | [1.057, 1.152] | 1 | 288/288 |

Confirmation winner/rho: instructions 0.7115, CI [0.6205792973445925, 0.818802940810528]; native time 0.9181, CI [0.8759967505844022, 0.9648840696231857].

Replay winner/rho: instructions 0.7115, CI [0.6205716505803323, 0.8188008542100041]; native time 0.9058, CI [0.8682743715265772, 0.9466055105523463].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.818 | 0.727 | 1 | 36 |
| n17a1 | 1.154 | 1.135 | 1 | 36 |
| n19a0 | 1.217 | 1.224 | 1 | 36 |
| n19a1 | 0.850 | 0.854 | 1 | 36 |
| n23a0 | 0.714 | 0.703 | 1 | 36 |
| n23a1 | 0.671 | 0.660 | 1 | 36 |
| n29a1 | 0.988 | 0.958 | 1 | 36 |
| n31a0 | 0.770 | 0.765 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
