# IC candidate tournament: round-0011

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1440 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **2.085**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 6,806 | 1 | 0.4846 | 2.448e+05 | 36/36 | reference |
| wordfield_batch4 | 6,846 | 1.006 | 0.4875 | 2.463e+05 | 36/36 | engineering experiment |
| wordfield_fastio | 6,563 | 0.9642 | 0.4673 | 2.361e+05 | 36/36 | engineering experiment |
| rho | 1.404e+04 | 2.063 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 6,462 | 1 | 0.4796 | 2.609e+05 | 180/180 | reference |
| wordfield_fastio | 6,236 | 0.9651 | 0.4629 | 2.518e+05 | 180/180 | engineering experiment |
| rho | 1.347e+04 | 2.085 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 6,462 | 1 | 0.4796 | 2.609e+05 | 180/180 | reference |
| wordfield_fastio | 6,236 | 0.9651 | 0.4629 | 2.518e+05 | 180/180 | engineering experiment |
| rho | 1.347e+04 | 2.085 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.017 | 1 | reference | 0.8851 | 36/36 |
| wordfield_batch4 | 2.011 | 0.997 | [0.97, 1.03] | 0.8825 | 36/36 |
| wordfield_fastio | 2.024 | 1.003 | [0.9791, 1.027] | 0.8881 | 36/36 |
| rho | 2.279 | 1.13 | [1.074, 1.183] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.003 | 1 | reference | 0.8653 | 180/180 |
| wordfield_fastio | 1.999 | 0.9982 | [0.988, 1.009] | 0.8638 | 180/180 |
| rho | 2.315 | 1.156 | [1.121, 1.191] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.998 | 1 | reference | 0.8654 | 180/180 |
| wordfield_fastio | 1.975 | 0.9885 | [0.9781, 0.9988] | 0.8555 | 180/180 |
| rho | 2.309 | 1.156 | [1.126, 1.183] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.4796, CI [0.449811117682557, 0.5216987948093519]; native time 0.8653, CI [0.8396219005399979, 0.8919422901900026].

Replay winner/rho: instructions 0.4796, CI [0.4498091000263744, 0.5216984655363398]; native time 0.8654, CI [0.8454206772431005, 0.8883025400214984].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.091 | 1.000 | 1 | 36 |
| n17a1 | 1.017 | 0.999 | 1 | 36 |
| n19a0 | 1.128 | 1.115 | 1 | 36 |
| n19a1 | 0.778 | 0.782 | 1 | 36 |
| n23a0 | 1.034 | 1.024 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
