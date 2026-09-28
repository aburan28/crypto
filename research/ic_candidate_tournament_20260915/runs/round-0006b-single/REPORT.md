# IC candidate tournament: round-0006b-single

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1584 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.3275**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.308e+05 | 1 | 3.026 | 4.706e+06 | 36/36 | reference |
| combined | 1.044e+05 | 0.7982 | 2.415 | 3.756e+06 | 36/36 | engineering experiment |
| combined_batch1 | 1.029e+05 | 0.7868 | 2.381 | 3.703e+06 | 36/36 | engineering experiment |
| descent | 1.243e+05 | 0.9501 | 2.875 | 4.471e+06 | 36/36 | engineering experiment |
| fast_orbits | 1.182e+05 | 0.9038 | 2.735 | 4.253e+06 | 36/36 | engineering experiment |
| lazy_field | 1.237e+05 | 0.9456 | 2.861 | 4.45e+06 | 36/36 | engineering experiment |
| rho | 4.323e+04 | 0.3305 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.211e+05 | 1 | 3.053 | 4.89e+06 | 180/180 | reference |
| combined_batch1 | 9.512e+04 | 0.7854 | 2.398 | 3.84e+06 | 180/180 | engineering experiment |
| rho | 3.967e+04 | 0.3275 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.211e+05 | 1 | 3.053 | 4.89e+06 | 180/180 | reference |
| combined_batch1 | 9.512e+04 | 0.7853 | 2.398 | 3.84e+06 | 180/180 | engineering experiment |
| rho | 3.967e+04 | 0.3275 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.002 | 1 | reference | 1.535 | 36/36 |
| combined | 5.391 | 0.8981 | [0.8785, 0.9231] | 1.379 | 36/36 |
| combined_batch1 | 5.339 | 0.8895 | [0.8751, 0.9065] | 1.366 | 36/36 |
| descent | 5.818 | 0.9692 | [0.959, 0.983] | 1.488 | 36/36 |
| fast_orbits | 5.722 | 0.9533 | [0.946, 0.9617] | 1.463 | 36/36 |
| lazy_field | 5.868 | 0.9777 | [0.9718, 0.9841] | 1.501 | 36/36 |
| rho | 3.91 | 0.6514 | [0.6431, 0.6617] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.186 | 1 | reference | 1.56 | 180/180 |
| combined_batch1 | 5.472 | 0.8846 | [0.8731, 0.8995] | 1.38 | 180/180 |
| rho | 3.966 | 0.6411 | [0.6283, 0.6546] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.348 | 1 | reference | 1.534 | 180/180 |
| combined_batch1 | 5.657 | 0.8912 | [0.8793, 0.9064] | 1.367 | 180/180 |
| rho | 4.138 | 0.6519 | [0.6397, 0.6668] | 1 | 180/180 |

Confirmation winner/rho: instructions 3.0532, CI [2.6670343409794497, 3.523270546881521]; native time 1.5598, CI [1.527628476463095, 1.5916609112502804].

Replay winner/rho: instructions 3.0533, CI [2.667046018772113, 3.523273316114213]; native time 1.5341, CI [1.499808917558032, 1.5632227530932101].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.954 | 0.864 | 1 | 36 |
| n17a1 | 0.590 | 0.572 | 1 | 36 |
| n19a0 | 1.217 | 1.203 | 1 | 36 |
| n19a1 | 0.691 | 0.682 | 1 | 36 |
| n23a0 | 1.099 | 1.084 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
