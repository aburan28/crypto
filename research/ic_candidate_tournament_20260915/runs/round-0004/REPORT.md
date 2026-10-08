# IC candidate tournament: round-0004

Decision: **promoted — folded_lift_batch4**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **2.350x**; candidate/incumbent ratio 0.4255, paired 95% interval [0.3981104765039304, 0.46471280387240294]. This is an engineering result on the tested workloads.

Independent audit checked **1584 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.3286**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.08e+05 | 1 | 7.145 | 1.108e+07 | 36/36 | reference |
| fast_lift | 2.295e+05 | 0.7452 | 5.324 | 8.257e+06 | 36/36 | engineering experiment |
| folded | 2.144e+05 | 0.696 | 4.973 | 7.713e+06 | 36/36 | engineering experiment |
| folded_lift | 1.358e+05 | 0.4409 | 3.15 | 4.885e+06 | 36/36 | engineering experiment |
| folded_lift_batch4 | 1.34e+05 | 0.435 | 3.108 | 4.821e+06 | 36/36 | engineering experiment |
| folded_lift_dense | 1.363e+05 | 0.4424 | 3.161 | 4.903e+06 | 36/36 | engineering experiment |
| rho | 4.311e+04 | 0.14 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.837e+05 | 1 | 7.152 | 1.145e+07 | 180/180 | reference |
| folded_lift_batch4 | 1.207e+05 | 0.4255 | 3.043 | 4.873e+06 | 180/180 | engineering experiment |
| rho | 3.967e+04 | 0.1398 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.837e+05 | 1 | 7.152 | 1.145e+07 | 180/180 | reference |
| folded_lift_batch4 | 1.207e+05 | 0.4255 | 3.043 | 4.873e+06 | 180/180 | engineering experiment |
| rho | 3.967e+04 | 0.1398 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 11.04 | 1 | reference | 2.764 | 36/36 |
| fast_lift | 8.856 | 0.8022 | [0.7777, 0.8289] | 2.218 | 36/36 |
| folded | 8.33 | 0.7546 | [0.7215, 0.7999] | 2.086 | 36/36 |
| folded_lift | 6.207 | 0.5622 | [0.5044, 0.6404] | 1.554 | 36/36 |
| folded_lift_batch4 | 6.192 | 0.5609 | [0.5034, 0.6429] | 1.55 | 36/36 |
| folded_lift_dense | 6.229 | 0.5642 | [0.5088, 0.6412] | 1.56 | 36/36 |
| rho | 3.993 | 0.3617 | [0.3139, 0.4286] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 11.39 | 1 | reference | 2.75 | 180/180 |
| folded_lift_batch4 | 6.299 | 0.5532 | [0.5017, 0.6208] | 1.521 | 180/180 |
| rho | 4.141 | 0.3636 | [0.3253, 0.417] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 11.55 | 1 | reference | 2.679 | 180/180 |
| folded_lift_batch4 | 6.503 | 0.5628 | [0.5111, 0.6296] | 1.508 | 180/180 |
| rho | 4.313 | 0.3733 | [0.3344, 0.4247] | 1 | 180/180 |

Confirmation winner/rho: instructions 3.0429, CI [2.688274518239304, 3.4833783647230367]; native time 1.5213, CI [1.4825577742485183, 1.5518235985269213].

Replay winner/rho: instructions 3.0429, CI [2.688259922267101, 3.4833654799448777]; native time 1.5078, CI [1.4776570931768824, 1.5368924247680655].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.045 | 0.954 | 1 | 36 |
| n17a1 | 1.308 | 1.290 | 1 | 36 |
| n19a0 | 0.768 | 0.755 | 1 | 36 |
| n19a1 | 0.946 | 0.936 | 1 | 36 |
| n23a0 | 1.006 | 1.006 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
