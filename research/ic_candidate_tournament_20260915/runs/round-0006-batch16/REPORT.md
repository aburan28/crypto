# IC candidate tournament: round-0006-batch16

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **16 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1680 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.378**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.676e+05 | 1 | 0.7321 | 9.625e+06 | 36/36 | reference |
| combined | 2.147e+05 | 0.8026 | 0.5875 | 7.725e+06 | 36/36 | engineering experiment |
| euclid_inv | 2.446e+05 | 0.9142 | 0.6692 | 8.799e+06 | 36/36 | engineering experiment |
| fast_curve_once | 2.569e+05 | 0.9603 | 0.703 | 9.243e+06 | 36/36 | engineering experiment |
| fast_orbits | 2.545e+05 | 0.9513 | 0.6964 | 9.156e+06 | 36/36 | engineering experiment |
| it_inv | 2.601e+05 | 0.9722 | 0.7117 | 9.358e+06 | 36/36 | engineering experiment |
| lambda_table | 2.529e+05 | 0.9452 | 0.6919 | 9.098e+06 | 36/36 | engineering experiment |
| lazy_field | 2.604e+05 | 0.9732 | 0.7125 | 9.368e+06 | 36/36 | engineering experiment |
| rho | 3.655e+05 | 1.366 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.471e+05 | 1 | 0.7259 | 9.978e+06 | 180/180 | reference |
| combined | 1.982e+05 | 0.8018 | 0.582 | 8.001e+06 | 180/180 | engineering experiment |
| rho | 3.405e+05 | 1.378 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.471e+05 | 1 | 0.7259 | 9.978e+06 | 180/180 | reference |
| combined | 1.982e+05 | 0.8018 | 0.582 | 8.001e+06 | 180/180 | engineering experiment |
| rho | 3.405e+05 | 1.378 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 13.83 | 1 | reference | 0.8005 | 36/36 |
| combined | 12.85 | 0.9296 | [0.8918, 0.9602] | 0.7441 | 36/36 |
| euclid_inv | 13.23 | 0.9572 | [0.9296, 0.9886] | 0.7662 | 36/36 |
| fast_curve_once | 13.92 | 1.007 | [0.9836, 1.027] | 0.806 | 36/36 |
| fast_orbits | 13.41 | 0.9698 | [0.9364, 1.001] | 0.7763 | 36/36 |
| it_inv | 13.77 | 0.9959 | [0.9628, 1.03] | 0.7972 | 36/36 |
| lambda_table | 13.65 | 0.9869 | [0.9323, 1.039] | 0.79 | 36/36 |
| lazy_field | 13.75 | 0.9947 | [0.9621, 1.04] | 0.7962 | 36/36 |
| rho | 17.27 | 1.249 | [1.064, 1.448] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 14.44 | 1 | reference | 0.7853 | 180/180 |
| combined | 13.34 | 0.9237 | [0.9102, 0.9365] | 0.7254 | 180/180 |
| rho | 18.39 | 1.273 | [1.141, 1.421] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 14.71 | 1 | reference | 0.7942 | 180/180 |
| combined | 13.44 | 0.9136 | [0.8999, 0.9292] | 0.7256 | 180/180 |
| rho | 18.52 | 1.259 | [1.117, 1.413] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7259, CI [0.6410479305566766, 0.8372029685742962]; native time 0.7853, CI [0.7040300575162505, 0.8771486916175191].

Replay winner/rho: instructions 0.7259, CI [0.6410507351339184, 0.8372075155399566]; native time 0.7942, CI [0.708093784170611, 0.8967991364826209].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.000 | 0.909 | 1 | 576 |
| n17a1 | 1.017 | 1.008 | 1 | 576 |
| n19a0 | 1.074 | 1.060 | 1 | 576 |
| n19a1 | 0.970 | 0.960 | 1 | 576 |
| n23a0 | 0.976 | 0.970 | 1 | 576 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
