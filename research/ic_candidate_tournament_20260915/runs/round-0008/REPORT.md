# IC candidate tournament: round-0008

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1356 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.355**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.171e+04 | 1 | 0.7422 | 1.141e+06 | 36/36 | reference |
| tiny2_cert_batch4 | 3.177e+04 | 1.002 | 0.7438 | 1.143e+06 | 36/36 | engineering experiment |
| rho | 4.272e+04 | 1.347 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.95e+04 | 1 | 0.7379 | 1.191e+06 | 180/180 | reference |
| tiny2_cert_batch4 | 2.956e+04 | 1.002 | 0.7393 | 1.193e+06 | 180/180 | engineering experiment |
| rho | 3.998e+04 | 1.355 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2.95e+04 | 1 | 0.7379 | 1.191e+06 | 180/180 | reference |
| tiny2_cert_batch4 | 2.956e+04 | 1.002 | 0.7393 | 1.193e+06 | 180/180 | engineering experiment |
| rho | 3.998e+04 | 1.355 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 7.077 | 1 | reference | 0.9755 | 36/36 |
| tiny2_cert_batch4 | 6.828 | 0.9649 | [0.8784, 1.022] | 0.9412 | 36/36 |
| rho | 7.255 | 1.025 | [0.9392, 1.096] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.758 | 1 | reference | 0.9355 | 180/180 |
| tiny2_cert_batch4 | 6.797 | 1.006 | [0.9863, 1.028] | 0.9408 | 180/180 |
| rho | 7.224 | 1.069 | [1.041, 1.096] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 7.329 | 1 | reference | 0.9414 | 180/180 |
| tiny2_cert_batch4 | 7.28 | 0.9934 | [0.9759, 1.011] | 0.9352 | 180/180 |
| rho | 7.785 | 1.062 | [1.048, 1.079] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7379, CI [0.7144337379505753, 0.7571723061927793]; native time 0.9355, CI [0.9124104876677653, 0.9609791776155655].

Replay winner/rho: instructions 0.7379, CI [0.7144328577963271, 0.7571684930212431]; native time 0.9414, CI [0.9266864038465168, 0.9545396471783393].

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
