# IC candidate tournament: round-0010

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1440 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.87**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.311e+04 | 1 | 0.5383 | 4.717e+05 | 36/36 | reference |
| lean_batch4 | 1.324e+04 | 1.01 | 0.5437 | 4.764e+05 | 36/36 | engineering experiment |
| lean_stdprobe | 1.311e+04 | 1 | 0.5383 | 4.717e+05 | 36/36 | engineering experiment |
| rho | 2.436e+04 | 1.858 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.171e+04 | 1 | 0.5347 | 4.727e+05 | 180/180 | reference |
| lean_stdprobe | 1.171e+04 | 1 | 0.5348 | 4.728e+05 | 180/180 | engineering experiment |
| rho | 2.19e+04 | 1.87 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.171e+04 | 1 | 0.5347 | 4.727e+05 | 180/180 | reference |
| lean_stdprobe | 1.171e+04 | 1 | 0.5348 | 4.728e+05 | 180/180 | engineering experiment |
| rho | 2.19e+04 | 1.87 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.44 | 1 | reference | 0.8678 | 36/36 |
| lean_batch4 | 2.426 | 0.9942 | [0.9615, 1.028] | 0.8628 | 36/36 |
| lean_stdprobe | 2.42 | 0.9915 | [0.9552, 1.03] | 0.8604 | 36/36 |
| rho | 2.812 | 1.152 | [1.084, 1.216] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.422 | 1 | reference | 0.8445 | 180/180 |
| lean_stdprobe | 2.42 | 0.9989 | [0.9893, 1.01] | 0.8437 | 180/180 |
| rho | 2.868 | 1.184 | [1.135, 1.231] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.467 | 1 | reference | 0.8665 | 180/180 |
| lean_stdprobe | 2.44 | 0.9889 | [0.9743, 1.002] | 0.8568 | 180/180 |
| rho | 2.847 | 1.154 | [1.103, 1.204] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.5347, CI [0.5105630570377481, 0.5644617225225833]; native time 0.8445, CI [0.8122502019521926, 0.8814419030825967].

Replay winner/rho: instructions 0.5347, CI [0.5105627651894524, 0.5644610150682358]; native time 0.8665, CI [0.830718721524484, 0.9077701942818931].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.818 | 0.727 | 1 | 36 |
| n17a1 | 0.808 | 0.790 | 1 | 36 |
| n19a0 | 0.992 | 0.979 | 1 | 36 |
| n19a1 | 0.658 | 0.648 | 1 | 36 |
| n23a0 | 0.920 | 0.912 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
