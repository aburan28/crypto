# IC candidate tournament: round-0009

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1356 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.814**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.425e+04 | 1 | 0.5598 | 5.128e+05 | 36/36 | reference |
| fastcurve_batch4 | 1.431e+04 | 1.004 | 0.5622 | 5.149e+05 | 36/36 | engineering experiment |
| rho | 2.546e+04 | 1.786 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.269e+04 | 1 | 0.5511 | 5.125e+05 | 180/180 | reference |
| fastcurve_batch4 | 1.276e+04 | 1.005 | 0.5539 | 5.151e+05 | 180/180 | engineering experiment |
| rho | 2.303e+04 | 1.814 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.269e+04 | 1 | 0.5511 | 5.125e+05 | 180/180 | reference |
| fastcurve_batch4 | 1.276e+04 | 1.005 | 0.5539 | 5.151e+05 | 180/180 | engineering experiment |
| rho | 2.303e+04 | 1.814 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 5.893 | 1 | reference | 0.945 | 36/36 |
| fastcurve_batch4 | 5.743 | 0.9746 | [0.9162, 1.011] | 0.921 | 36/36 |
| rho | 6.236 | 1.058 | [0.976, 1.111] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.091 | 1 | reference | 0.9406 | 180/180 |
| fastcurve_batch4 | 6.021 | 0.9885 | [0.9711, 1.008] | 0.9298 | 180/180 |
| rho | 6.475 | 1.063 | [1.019, 1.095] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.353 | 1 | reference | 0.9323 | 180/180 |
| fastcurve_batch4 | 6.308 | 0.9928 | [0.9761, 1.01] | 0.9256 | 180/180 |
| rho | 6.815 | 1.073 | [1.044, 1.106] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.5511, CI [0.5279422621731197, 0.5880928655396426]; native time 0.9406, CI [0.9131071744076471, 0.9811067742684935].

Replay winner/rho: instructions 0.5511, CI [0.5279408400938557, 0.5880942826604277]; native time 0.9323, CI [0.9041124782345309, 0.9580971485046996].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.182 | 1.091 | 1 | 36 |
| n17a1 | 0.699 | 0.681 | 1 | 36 |
| n19a0 | 1.203 | 1.190 | 1 | 36 |
| n19a1 | 1.118 | 1.123 | 1 | 36 |
| n23a0 | 1.116 | 1.106 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
