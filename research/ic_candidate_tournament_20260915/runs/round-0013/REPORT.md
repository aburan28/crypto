# IC candidate tournament: round-0013

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **1440 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.299**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,545 | 1 | 0.7342 | 1.635e+05 | 12/12 | reference |
| ld_canon | 4,439 | 0.9768 | 0.7171 | 1.597e+05 | 12/12 | engineering experiment |
| ld_ic | 4,189 | 0.9218 | 0.6768 | 1.507e+05 | 12/12 | engineering experiment |
| rho | 6,190 | 1.362 | 1 | unmeasured | 12/12 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,462 | 1 | 0.7627 | 1.605e+05 | 36/36 | reference |
| ld_canon | 4,360 | 0.9771 | 0.7453 | 1.568e+05 | 36/36 | engineering experiment |
| ld_ic | 4,111 | 0.9214 | 0.7028 | 1.479e+05 | 36/36 | engineering experiment |
| rho | 5,850 | 1.311 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,147 | 1 | 0.7697 | 1.674e+05 | 180/180 | reference |
| ld_ic | 3,824 | 0.9221 | 0.7097 | 1.544e+05 | 180/180 | engineering experiment |
| rho | 5,388 | 1.299 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4,147 | 1 | 0.7697 | 1.674e+05 | 180/180 | reference |
| ld_ic | 3,824 | 0.9221 | 0.7097 | 1.544e+05 | 180/180 | engineering experiment |
| rho | 5,388 | 1.299 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.154 | 1 | reference | 0.9115 | 36/36 |
| ld_canon | 1.16 | 1.005 | [0.9749, 1.036] | 0.9161 | 36/36 |
| ld_ic | 1.13 | 0.9791 | [0.9624, 1.004] | 0.8925 | 36/36 |
| rho | 1.266 | 1.097 | [1.064, 1.131] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.169 | 1 | reference | 0.9098 | 180/180 |
| ld_ic | 1.192 | 1.02 | [1, 1.044] | 0.9279 | 180/180 |
| rho | 1.285 | 1.099 | [1.074, 1.124] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.123 | 1 | reference | 0.9206 | 180/180 |
| ld_ic | 1.114 | 0.9919 | [0.978, 1.008] | 0.9132 | 180/180 |
| rho | 1.22 | 1.086 | [1.055, 1.113] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7697, CI [0.7147839637637708, 0.8370425765657616]; native time 0.9098, CI [0.8894023669954111, 0.931344422318841].

Replay winner/rho: instructions 0.7697, CI [0.7147697917530964, 0.837057854302903]; native time 0.9206, CI [0.8982663875921412, 0.9480279891495028].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.136 | 1.045 | 1 | 36 |
| n17a1 | 0.845 | 0.827 | 1 | 36 |
| n19a0 | 0.646 | 0.632 | 1 | 36 |
| n19a1 | 0.710 | 0.701 | 1 | 36 |
| n23a0 | 1.056 | 1.052 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
