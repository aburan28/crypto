# IC candidate tournament: round-0019

Decision: **promoted — both**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **1.036x**; candidate/incumbent ratio 0.9657, paired 95% interval [0.9478832169875092, 0.9824838715022691]. This is an engineering result on the tested workloads.

Independent audit checked **2646 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 124 fresh inputs over 8 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.481**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,850 | 1 | 0.6475 | 1.804e+05 | 18/18 | reference |
| both | 2,760 | 0.9685 | 0.6271 | 1.747e+05 | 18/18 | engineering experiment |
| rho | 4,401 | 1.544 | 1 | unmeasured | 18/18 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,915 | 1 | 0.6768 | 1.845e+05 | 54/54 | reference |
| both | 2,826 | 0.9694 | 0.6561 | 1.789e+05 | 54/54 | engineering experiment |
| rho | 4,308 | 1.478 | 1 | unmeasured | 54/54 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,213 | 1 | 0.6992 | 1.82e+05 | 372/372 | reference |
| both | 3,103 | 0.9657 | 0.6751 | 1.758e+05 | 372/372 | engineering experiment |
| rho | 4,596 | 1.43 | 1 | unmeasured | 372/372 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,213 | 1 | 0.6991 | 1.82e+05 | 372/372 | reference |
| both | 3,103 | 0.9657 | 0.6751 | 1.758e+05 | 372/372 | engineering experiment |
| rho | 4,596 | 1.43 | 1 | unmeasured | 372/372 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.004 | 1 | reference | 0.8897 | 54/54 |
| both | 0.96 | 0.9565 | [0.9124, 0.9856] | 0.8509 | 54/54 |
| rho | 1.128 | 1.124 | [1.029, 1.194] | 1 | 54/54 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.016 | 1 | reference | 0.9158 | 372/372 |
| both | 0.9821 | 0.9662 | [0.9516, 0.9792] | 0.8848 | 372/372 |
| rho | 1.11 | 1.092 | [1.056, 1.124] | 1 | 372/372 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.056 | 1 | reference | 0.9075 | 372/372 |
| both | 1.033 | 0.9787 | [0.9649, 0.9924] | 0.8882 | 372/372 |
| rho | 1.163 | 1.102 | [1.067, 1.132] | 1 | 372/372 |

Confirmation winner/rho: instructions 0.6751, CI [0.6297545161748866, 0.7265491892286999]; native time 0.8848, CI [0.8623695785913571, 0.9097683368504219].

Replay winner/rho: instructions 0.6751, CI [0.6297303594398972, 0.7265456165881904]; native time 0.8882, CI [0.8675748494788712, 0.9127017410325569].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.227 | 1.136 | 1 | 36 |
| n17a1 | 0.636 | 0.618 | 1 | 36 |
| n19a0 | 0.673 | 0.659 | 1 | 36 |
| n19a1 | 1.022 | 1.013 | 1 | 36 |
| n23a0 | 1.181 | 1.174 | 1 | 36 |
| n23a1 | 1.177 | 1.173 | 1 | 120 |
| n29a1 | 1.032 | 1.003 | 1 | 36 |
| n31a0 | 0.605 | 0.597 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
