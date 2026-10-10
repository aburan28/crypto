# IC candidate tournament: round-0014

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1440 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.279**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,411 | 1 | 0.8275 | 1.587e+05 | 12/12 | reference |
| lto_canon | 4,310 | 0.9773 | 0.8087 | 1.551e+05 | 12/12 | engineering experiment |
| lto_ic | 4,071 | 0.9229 | 0.7637 | 1.464e+05 | 12/12 | engineering experiment |
| rho | 5,330 | 1.208 | 1 | unmeasured | 12/12 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,291 | 1 | 0.731 | 1.544e+05 | 36/36 | reference |
| lto_canon | 4,199 | 0.9785 | 0.7153 | 1.511e+05 | 36/36 | engineering experiment |
| lto_ic | 3,959 | 0.9226 | 0.6745 | 1.424e+05 | 36/36 | engineering experiment |
| rho | 5,870 | 1.368 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,096 | 1 | 0.7818 | 1.654e+05 | 180/180 | reference |
| lto_ic | 3,779 | 0.9226 | 0.7212 | 1.526e+05 | 180/180 | engineering experiment |
| rho | 5,240 | 1.279 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,096 | 1 | 0.7818 | 1.654e+05 | 180/180 | reference |
| lto_ic | 3,779 | 0.9226 | 0.7212 | 1.526e+05 | 180/180 | engineering experiment |
| rho | 5,240 | 1.279 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.1 | 1 | reference | 0.889 | 36/36 |
| lto_canon | 1.117 | 1.015 | [0.9813, 1.054] | 0.9022 | 36/36 |
| lto_ic | 1.1 | 1 | [0.9754, 1.024] | 0.889 | 36/36 |
| rho | 1.238 | 1.125 | [1.091, 1.159] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.152 | 1 | reference | 0.9424 | 180/180 |
| lto_ic | 1.146 | 0.9947 | [0.9724, 1.016] | 0.9374 | 180/180 |
| rho | 1.222 | 1.061 | [1.029, 1.096] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.122 | 1 | reference | 0.9284 | 180/180 |
| lto_ic | 1.124 | 1.002 | [0.9825, 1.026] | 0.9305 | 180/180 |
| rho | 1.208 | 1.077 | [1.039, 1.105] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7818, CI [0.7199868031943767, 0.8600593236734697]; native time 0.9424, CI [0.912629257744672, 0.9715184462908083].

Replay winner/rho: instructions 0.7818, CI [0.7199903911091349, 0.8600633389135146]; native time 0.9284, CI [0.9049721702559972, 0.9623543289813028].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.954 | 0.864 | 1 | 36 |
| n17a1 | 0.718 | 0.727 | 1 | 36 |
| n19a0 | 0.741 | 0.727 | 1 | 36 |
| n19a1 | 0.782 | 0.787 | 1 | 36 |
| n23a0 | 0.875 | 0.867 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
