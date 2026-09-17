# IC candidate tournament: round-0008-single-implementation

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1536 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.887**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4.789e+04 | 1 | 1.115 | 6.665e+06 | 12/12 | reference |
| combined | 4.571e+04 | 0.9546 | 1.065 | 6.362e+06 | 12/12 | engineering experiment |
| combined_dense | 4.557e+04 | 0.9517 | 1.062 | 6.343e+06 | 12/12 | engineering experiment |
| orbit_projection | 4.573e+04 | 0.9549 | 1.065 | 6.364e+06 | 12/12 | engineering experiment |
| serial_collection | 4.788e+04 | 0.9998 | 1.115 | 6.663e+06 | 12/12 | engineering experiment |
| rho | 4.293e+04 | 0.8965 | 1 | unmeasured | 12/12 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4.962e+04 | 1 | 1.166 | 6.905e+06 | 36/36 | reference |
| combined | 4.742e+04 | 0.9558 | 1.114 | 6.6e+06 | 36/36 | engineering experiment |
| combined_dense | 4.728e+04 | 0.9528 | 1.111 | 6.58e+06 | 36/36 | engineering experiment |
| orbit_projection | 4.744e+04 | 0.9561 | 1.115 | 6.602e+06 | 36/36 | engineering experiment |
| serial_collection | 4.961e+04 | 0.9998 | 1.166 | 6.904e+06 | 36/36 | engineering experiment |
| rho | 4.256e+04 | 0.8577 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4.489e+04 | 1 | 1.127 | 7.059e+06 | 180/180 | reference |
| combined_dense | 4.294e+04 | 0.9564 | 1.078 | 6.751e+06 | 180/180 | engineering experiment |
| rho | 3.982e+04 | 0.887 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4.489e+04 | 1 | 1.127 | 7.059e+06 | 180/180 | reference |
| combined_dense | 4.294e+04 | 0.9564 | 1.078 | 6.751e+06 | 180/180 | engineering experiment |
| rho | 3.982e+04 | 0.887 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 4.301 | 1 | reference | 1.096 | 36/36 |
| combined | 4.24 | 0.9859 | [0.9781, 0.9935] | 1.08 | 36/36 |
| combined_dense | 4.238 | 0.9855 | [0.9778, 0.9956] | 1.08 | 36/36 |
| orbit_projection | 4.258 | 0.9902 | [0.9848, 0.9949] | 1.085 | 36/36 |
| serial_collection | 4.32 | 1.004 | [1.001, 1.01] | 1.101 | 36/36 |
| rho | 3.925 | 0.9127 | [0.8693, 0.9419] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 4.406 | 1 | reference | 1.082 | 180/180 |
| combined_dense | 4.343 | 0.9858 | [0.9814, 0.9905] | 1.066 | 180/180 |
| rho | 4.074 | 0.9246 | [0.9086, 0.937] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 4.573 | 1 | reference | 1.083 | 180/180 |
| combined_dense | 4.499 | 0.9838 | [0.9788, 0.9889] | 1.066 | 180/180 |
| rho | 4.222 | 0.9232 | [0.907, 0.9358] | 1 | 180/180 |

Confirmation winner/rho: instructions 1.1274, CI [1.0664394518628761, 1.1978968142897077]; native time 1.0816, CI [1.0672187767654386, 1.1006507319022558].

Replay winner/rho: instructions 1.1274, CI [1.0664410849430652, 1.1978933413637658]; native time 1.0831, CI [1.06860936740454, 1.102554906966594].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.136 | 1.045 | 1 | 36 |
| n17a1 | 0.690 | 0.699 | 1 | 36 |
| n19a0 | 0.863 | 0.850 | 1 | 36 |
| n19a1 | 0.749 | 0.739 | 1 | 36 |
| n23a0 | 1.090 | 1.090 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
