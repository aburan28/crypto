# IC candidate tournament: round-0003b

Decision: **promoted — combined_batch8**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **1.821x**; candidate/incumbent ratio 0.5490, paired 95% interval [0.5401764754277694, 0.5577892780535699]. This is an engineering result on the tested workloads.

Independent audit checked **1632 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.1399**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.637e+05 | 1 | 13.28 | 2.028e+07 | 36/36 | reference |
| batch8 | 5.369e+05 | 0.9525 | 12.65 | 1.931e+07 | 36/36 | engineering experiment |
| combined | 3.153e+05 | 0.5594 | 7.43 | 1.134e+07 | 36/36 | engineering experiment |
| combined_batch8 | 3.097e+05 | 0.5494 | 7.297 | 1.114e+07 | 36/36 | engineering experiment |
| fast_cofactor | 5.145e+05 | 0.9128 | 12.12 | 1.851e+07 | 36/36 | engineering experiment |
| fast_verify | 3.706e+05 | 0.6575 | 8.733 | 1.333e+07 | 36/36 | engineering experiment |
| serial_pairs | 5.561e+05 | 0.9865 | 13.1 | 2e+07 | 36/36 | engineering experiment |
| rho | 4.244e+04 | 0.07529 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.159e+05 | 1 | 13.02 | 2.083e+07 | 180/180 | reference |
| combined_batch8 | 2.832e+05 | 0.549 | 7.149 | 1.143e+07 | 180/180 | engineering experiment |
| rho | 3.961e+04 | 0.07679 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.159e+05 | 1 | 13.02 | 2.083e+07 | 180/180 | reference |
| combined_batch8 | 2.832e+05 | 0.549 | 7.149 | 1.143e+07 | 180/180 | engineering experiment |
| rho | 3.961e+04 | 0.07679 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 17.93 | 1 | reference | 4.751 | 36/36 |
| batch8 | 17.27 | 0.9635 | [0.9377, 0.9882] | 4.578 | 36/36 |
| combined | 10.89 | 0.6073 | [0.5699, 0.6515] | 2.885 | 36/36 |
| combined_batch8 | 10.84 | 0.6047 | [0.5675, 0.6494] | 2.873 | 36/36 |
| fast_cofactor | 16.56 | 0.924 | [0.9093, 0.9395] | 4.39 | 36/36 |
| fast_verify | 12.63 | 0.7046 | [0.678, 0.7405] | 3.347 | 36/36 |
| serial_pairs | 17.51 | 0.9767 | [0.9703, 0.9824] | 4.64 | 36/36 |
| rho | 3.773 | 0.2105 | [0.1729, 0.2638] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 18.65 | 1 | reference | 4.563 | 180/180 |
| combined_batch8 | 11.27 | 0.6044 | [0.5687, 0.6473] | 2.758 | 180/180 |
| rho | 4.087 | 0.2192 | [0.1844, 0.2686] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 18.71 | 1 | reference | 4.583 | 180/180 |
| combined_batch8 | 11.29 | 0.6037 | [0.569, 0.6455] | 2.767 | 180/180 |
| rho | 4.082 | 0.2182 | [0.1846, 0.2661] | 1 | 180/180 |

Confirmation winner/rho: instructions 7.1494, CI [6.563380843095217, 7.844932268218591]; native time 2.7579, CI [2.408587066800431, 3.0866033786867497].

Replay winner/rho: instructions 7.1493, CI [6.563350841032917, 7.844923582031123]; native time 2.7666, CI [2.422024258021222, 3.0896644539471403].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.909 | 0.954 | 1 | 36 |
| n17a1 | 0.963 | 0.945 | 1 | 36 |
| n19a0 | 0.863 | 0.870 | 1 | 36 |
| n19a1 | 0.768 | 0.758 | 1 | 36 |
| n23a0 | 1.091 | 1.077 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
