# IC candidate tournament: round-0015

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1440 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.356**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,468 | 1 | 0.7217 | 1.607e+05 | 12/12 | reference |
| scan | 4,294 | 0.961 | 0.6935 | 1.545e+05 | 12/12 | engineering experiment |
| scan_io | 4,053 | 0.9071 | 0.6547 | 1.458e+05 | 12/12 | engineering experiment |
| rho | 6,191 | 1.386 | 1 | unmeasured | 12/12 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,444 | 1 | 0.7452 | 1.599e+05 | 36/36 | reference |
| scan | 4,269 | 0.9606 | 0.7158 | 1.536e+05 | 36/36 | engineering experiment |
| scan_io | 4,028 | 0.9064 | 0.6754 | 1.449e+05 | 36/36 | engineering experiment |
| rho | 5,964 | 1.342 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,050 | 1 | 0.7372 | 1.635e+05 | 180/180 | reference |
| scan_io | 3,669 | 0.9059 | 0.6679 | 1.481e+05 | 180/180 | engineering experiment |
| rho | 5,493 | 1.356 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,050 | 1 | 0.7372 | 1.635e+05 | 180/180 | reference |
| scan_io | 3,669 | 0.9059 | 0.6679 | 1.481e+05 | 180/180 | engineering experiment |
| rho | 5,493 | 1.356 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.107 | 1 | reference | 0.925 | 36/36 |
| scan | 1.101 | 0.9954 | [0.9628, 1.043] | 0.9207 | 36/36 |
| scan_io | 1.132 | 1.023 | [0.9837, 1.081] | 0.9465 | 36/36 |
| rho | 1.196 | 1.081 | [1.012, 1.142] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.102 | 1 | reference | 0.9246 | 180/180 |
| scan_io | 1.094 | 0.9925 | [0.9785, 1.007] | 0.9176 | 180/180 |
| rho | 1.192 | 1.082 | [1.052, 1.107] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.075 | 1 | reference | 0.9157 | 180/180 |
| scan_io | 1.07 | 0.9957 | [0.983, 1.01] | 0.9117 | 180/180 |
| rho | 1.174 | 1.092 | [1.046, 1.126] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7372, CI [0.6627607474515894, 0.8325954872020435]; native time 0.9246, CI [0.9034237860884715, 0.9502974853533789].

Replay winner/rho: instructions 0.7372, CI [0.6627657681259977, 0.8326116460823377]; native time 0.9157, CI [0.8886201722807316, 0.9563789203616806].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.091 | 1.000 | 1 | 36 |
| n17a1 | 0.899 | 0.881 | 1 | 36 |
| n19a0 | 1.067 | 1.074 | 1 | 36 |
| n19a1 | 1.075 | 1.066 | 1 | 36 |
| n23a0 | 0.865 | 0.858 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
