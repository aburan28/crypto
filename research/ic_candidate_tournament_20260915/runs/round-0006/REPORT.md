# IC candidate tournament: round-0006

Decision: **promoted — tiny_batch1**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **3.700x**; candidate/incumbent ratio 0.2702, paired 95% interval [0.2319712081986235, 0.3073309863776135]. This is an engineering result on the tested workloads.

Independent audit checked **1584 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.21**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.317e+05 | 1 | 3.032 | 4.737e+06 | 36/36 | reference |
| combined_descent | 1.25e+05 | 0.9493 | 2.878 | 4.497e+06 | 36/36 | engineering experiment |
| fast_report | 1.214e+05 | 0.9218 | 2.795 | 4.366e+06 | 36/36 | engineering experiment |
| tiny | 3.569e+04 | 0.2711 | 0.8218 | 1.284e+06 | 36/36 | engineering experiment |
| tiny_batch1 | 3.564e+04 | 0.2707 | 0.8205 | 1.282e+06 | 36/36 | engineering experiment |
| tiny_fulltable | 3.667e+04 | 0.2785 | 0.8443 | 1.319e+06 | 36/36 | engineering experiment |
| rho | 4.343e+04 | 0.3299 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.219e+05 | 1 | 3.058 | 4.921e+06 | 180/180 | reference |
| tiny_batch1 | 3.294e+04 | 0.2702 | 0.8264 | 1.33e+06 | 180/180 | engineering experiment |
| rho | 3.986e+04 | 0.327 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.219e+05 | 1 | 3.058 | 4.921e+06 | 180/180 | reference |
| tiny_batch1 | 3.294e+04 | 0.2702 | 0.8264 | 1.33e+06 | 180/180 | engineering experiment |
| rho | 3.986e+04 | 0.327 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 10.71 | 1 | reference | 1.489 | 36/36 |
| combined_descent | 10.18 | 0.9507 | [0.9198, 0.9801] | 1.416 | 36/36 |
| fast_report | 10.23 | 0.9547 | [0.9258, 0.9819] | 1.422 | 36/36 |
| tiny | 6.935 | 0.6475 | [0.6201, 0.6755] | 0.9642 | 36/36 |
| tiny_batch1 | 7.035 | 0.6568 | [0.6192, 0.7047] | 0.9781 | 36/36 |
| tiny_fulltable | 6.971 | 0.6508 | [0.61, 0.6959] | 0.9692 | 36/36 |
| rho | 7.192 | 0.6715 | [0.635, 0.7097] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 10.89 | 1 | reference | 1.5 | 180/180 |
| tiny_batch1 | 7.092 | 0.6511 | [0.6236, 0.6801] | 0.9766 | 180/180 |
| rho | 7.262 | 0.6667 | [0.6489, 0.6874] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 11.26 | 1 | reference | 1.466 | 180/180 |
| tiny_batch1 | 7.568 | 0.6721 | [0.6451, 0.7023] | 0.9854 | 180/180 |
| rho | 7.68 | 0.6821 | [0.6635, 0.7093] | 1 | 180/180 |

Confirmation winner/rho: instructions 0.8264, CI [0.8032467836918035, 0.8488166193357406]; native time 0.9766, CI [0.957057247571661, 0.9946203513598987].

Replay winner/rho: instructions 0.8264, CI [0.8032473903248762, 0.8488169112953732]; native time 0.9854, CI [0.9641349965391578, 1.006241323144332].

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
