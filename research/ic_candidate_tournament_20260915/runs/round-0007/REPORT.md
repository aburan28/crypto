# IC candidate tournament: round-0007

Decision: **promoted — tiny2**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **1.127x**; candidate/incumbent ratio 0.8872, paired 95% interval [0.8493396285853244, 0.9241554327628069]. This is an engineering result on the tested workloads.

Independent audit checked **1488 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.356**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.508e+04 | 1 | 0.8398 | 1.262e+06 | 36/36 | reference |
| tiny2 | 3.11e+04 | 0.8864 | 0.7444 | 1.119e+06 | 36/36 | engineering experiment |
| tiny2_arith | 3.41e+04 | 0.9719 | 0.8162 | 1.227e+06 | 36/36 | engineering experiment |
| tiny2_rows | 3.203e+04 | 0.913 | 0.7667 | 1.152e+06 | 36/36 | engineering experiment |
| rho | 4.178e+04 | 1.191 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.298e+04 | 1 | 0.8313 | 1.331e+06 | 180/180 | reference |
| tiny2 | 2.926e+04 | 0.8872 | 0.7376 | 1.181e+06 | 180/180 | engineering experiment |
| rho | 3.967e+04 | 1.203 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.298e+04 | 1 | 0.8313 | 1.331e+06 | 180/180 | reference |
| tiny2 | 2.926e+04 | 0.8872 | 0.7376 | 1.181e+06 | 180/180 | engineering experiment |
| rho | 3.967e+04 | 1.203 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 7.286 | 1 | reference | 0.9657 | 36/36 |
| tiny2 | 6.993 | 0.9598 | [0.9212, 1.002] | 0.9269 | 36/36 |
| tiny2_arith | 7.084 | 0.9722 | [0.9489, 0.9948] | 0.9389 | 36/36 |
| tiny2_rows | 7.144 | 0.9806 | [0.9453, 1.019] | 0.9469 | 36/36 |
| rho | 7.545 | 1.035 | [1.012, 1.068] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 7.538 | 1 | reference | 0.9741 | 180/180 |
| tiny2 | 7.321 | 0.9712 | [0.9531, 0.9862] | 0.9461 | 180/180 |
| rho | 7.738 | 1.027 | [1.005, 1.046] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 7.536 | 1 | reference | 0.973 | 180/180 |
| tiny2 | 7.292 | 0.9675 | [0.9477, 0.9891] | 0.9414 | 180/180 |
| rho | 7.746 | 1.028 | [1.01, 1.049] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.7376, CI [0.7065693377992373, 0.7643545073767521]; native time 0.9461, CI [0.9292926395602298, 0.9639483494225052].

Replay winner/rho: instructions 0.7376, CI [0.7065691551095401, 0.7643536706925478]; native time 0.9414, CI [0.9215169986215078, 0.9610581271933554].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.818 | 0.727 | 1 | 36 |
| n17a1 | 1.035 | 1.017 | 1 | 36 |
| n19a0 | 0.938 | 0.925 | 1 | 36 |
| n19a1 | 0.888 | 0.878 | 1 | 36 |
| n23a0 | 0.768 | 0.763 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
