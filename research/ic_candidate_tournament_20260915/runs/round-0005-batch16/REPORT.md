# IC candidate tournament: round-0005-batch16

Decision: **promoted — combined_descent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **16 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **1.390x**; candidate/incumbent ratio 0.7193, paired 95% interval [0.7076301072093063, 0.7348375101951858]. This is an engineering result on the tested workloads.

Independent audit checked **1488 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.39**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.715e+05 | 1 | 1.008 | 1.336e+07 | 36/36 | reference |
| combined_descent | 2.679e+05 | 0.7212 | 0.7271 | 9.637e+06 | 36/36 | engineering experiment |
| fast_descent_check | 2.803e+05 | 0.7545 | 0.7607 | 1.008e+07 | 36/36 | engineering experiment |
| lazy_descent | 3.59e+05 | 0.9665 | 0.9744 | 1.292e+07 | 36/36 | engineering experiment |
| rho | 3.684e+05 | 0.9919 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.398e+05 | 1 | 1 | 1.372e+07 | 180/180 | reference |
| combined_descent | 2.444e+05 | 0.7193 | 0.7195 | 9.867e+06 | 180/180 | engineering experiment |
| rho | 3.397e+05 | 0.9997 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.398e+05 | 1 | 1 | 1.372e+07 | 180/180 | reference |
| combined_descent | 2.444e+05 | 0.7193 | 0.7195 | 9.867e+06 | 180/180 | engineering experiment |
| rho | 3.397e+05 | 0.9997 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 12.58 | 1 | reference | 1.015 | 36/36 |
| combined_descent | 9.741 | 0.7743 | [0.74, 0.8208] | 0.7863 | 36/36 |
| fast_descent_check | 10.09 | 0.8018 | [0.7676, 0.8462] | 0.8142 | 36/36 |
| lazy_descent | 12.26 | 0.9744 | [0.9666, 0.9803] | 0.9895 | 36/36 |
| rho | 12.39 | 0.9848 | [0.9095, 1.07] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 13.01 | 1 | reference | 1.015 | 180/180 |
| combined_descent | 9.995 | 0.7684 | [0.7396, 0.8055] | 0.7796 | 180/180 |
| rho | 12.82 | 0.9857 | [0.9223, 1.058] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 13.15 | 1 | reference | 1.015 | 180/180 |
| combined_descent | 10.14 | 0.7709 | [0.7428, 0.8053] | 0.7822 | 180/180 |
| rho | 12.96 | 0.9856 | [0.9218, 1.058] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7195, CI [0.6374471292425982, 0.8219836239443673]; native time 0.7796, CI [0.6991620583548833, 0.8736653449888768].

Replay winner/rho: instructions 0.7195, CI [0.6374460644319612, 0.8219757894629396]; native time 0.7822, CI [0.7024882235140708, 0.8739996568686551].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.000 | 0.909 | 1 | 576 |
| n17a1 | 0.936 | 0.917 | 1 | 576 |
| n19a0 | 0.925 | 0.911 | 1 | 576 |
| n19a1 | 0.989 | 0.989 | 1 | 576 |
| n23a0 | 0.877 | 0.867 | 1 | 576 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
