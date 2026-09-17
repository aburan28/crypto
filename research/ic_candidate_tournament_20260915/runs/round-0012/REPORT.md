# IC candidate tournament: round-0012

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1488 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **2.289**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5,524 | 1 | 0.468 | 1.987e+05 | 36/36 | reference |
| arena_fastio | 5,273 | 0.9545 | 0.4467 | 1.897e+05 | 36/36 | engineering experiment |
| arena_glibc | 5,848 | 1.059 | 0.4954 | 2.104e+05 | 36/36 | engineering experiment |
| musl_sysalloc | 7,820 | 1.416 | 0.6625 | 2.813e+05 | 36/36 | engineering experiment |
| rho | 1.18e+04 | 2.137 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,979 | 1 | 0.4369 | 2.01e+05 | 180/180 | reference |
| arena_fastio | 4,750 | 0.9541 | 0.4169 | 1.918e+05 | 180/180 | engineering experiment |
| rho | 1.139e+04 | 2.289 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,979 | 1 | 0.4369 | 2.01e+05 | 180/180 | reference |
| arena_fastio | 4,750 | 0.9541 | 0.4169 | 1.918e+05 | 180/180 | engineering experiment |
| rho | 1.139e+04 | 2.289 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.153 | 1 | reference | 0.8241 | 36/36 |
| arena_fastio | 1.128 | 0.9781 | [0.9625, 1.008] | 0.806 | 36/36 |
| arena_glibc | 2.029 | 1.759 | [1.65, 1.865] | 1.45 | 36/36 |
| musl_sysalloc | 1.416 | 1.228 | [1.155, 1.288] | 1.012 | 36/36 |
| rho | 1.399 | 1.213 | [1.161, 1.271] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.162 | 1 | reference | 0.7988 | 180/180 |
| arena_fastio | 1.171 | 1.008 | [0.9904, 1.027] | 0.805 | 180/180 |
| rho | 1.454 | 1.252 | [1.202, 1.304] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.169 | 1 | reference | 0.7914 | 180/180 |
| arena_fastio | 1.169 | 1 | [0.9882, 1.014] | 0.7914 | 180/180 |
| rho | 1.477 | 1.264 | [1.209, 1.323] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.4369, CI [0.4069034566370304, 0.48145735501622794]; native time 0.7988, CI [0.7667779157944749, 0.8318517878456885].

Replay winner/rho: instructions 0.4369, CI [0.40690298023525534, 0.4814550806508336]; native time 0.7914, CI [0.7557909440653447, 0.8270206925674652].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.182 | 1.091 | 1 | 36 |
| n17a1 | 1.126 | 1.108 | 1 | 36 |
| n19a0 | 1.033 | 1.020 | 1 | 36 |
| n19a1 | 0.941 | 0.931 | 1 | 36 |
| n23a0 | 1.075 | 1.069 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
