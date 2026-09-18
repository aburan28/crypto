# IC candidate tournament: round-0016

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **2142 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 96 fresh inputs over 8 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.193**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,374 | 1 | 0.813 | 2.136e+05 | 18/18 | reference |
| scan_io | 2,997 | 0.8882 | 0.7221 | 1.897e+05 | 18/18 | engineering experiment |
| rho | 4,151 | 1.23 | 1 | unmeasured | 18/18 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,531 | 1 | 0.8348 | 2.235e+05 | 54/54 | reference |
| scan_io | 3,130 | 0.8865 | 0.74 | 1.981e+05 | 54/54 | engineering experiment |
| rho | 4,230 | 1.198 | 1 | unmeasured | 54/54 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,869 | 1 | 0.8383 | 2.192e+05 | 288/288 | reference |
| scan_io | 3,433 | 0.8873 | 0.7438 | 1.945e+05 | 288/288 | engineering experiment |
| rho | 4,616 | 1.193 | 1 | unmeasured | 288/288 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,869 | 1 | 0.8383 | 2.192e+05 | 288/288 | reference |
| scan_io | 3,433 | 0.8873 | 0.7438 | 1.945e+05 | 288/288 | engineering experiment |
| rho | 4,616 | 1.193 | 1 | unmeasured | 288/288 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.177 | 1 | reference | 0.982 | 54/54 |
| scan_io | 1.136 | 0.9647 | [0.9243, 1.001] | 0.9473 | 54/54 |
| rho | 1.199 | 1.018 | [0.973, 1.06] | 1 | 54/54 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.252 | 1 | reference | 0.9779 | 288/288 |
| scan_io | 1.225 | 0.978 | [0.9608, 0.9947] | 0.9564 | 288/288 |
| rho | 1.281 | 1.023 | [0.9722, 1.071] | 1 | 288/288 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.187 | 1 | reference | 0.9482 | 288/288 |
| scan_io | 1.172 | 0.9871 | [0.9721, 1.003] | 0.936 | 288/288 |
| rho | 1.252 | 1.055 | [1.002, 1.107] | 1 | 288/288 |

Confirmation winner/rho: instructions 0.8383, CI [0.747866139139291, 0.9420452776988771]; native time 0.9779, CI [0.9338506276887212, 1.0288208398184095].

Replay winner/rho: instructions 0.8383, CI [0.7478521099449237, 0.9420334329845709]; native time 0.9482, CI [0.9034677107711173, 0.998022604447831].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.091 | 1.000 | 1 | 36 |
| n17a1 | 0.845 | 0.827 | 1 | 36 |
| n19a0 | 1.162 | 1.149 | 1 | 36 |
| n19a1 | 0.634 | 0.624 | 1 | 36 |
| n23a0 | 1.168 | 1.164 | 1 | 36 |
| n23a1 | 0.936 | 0.931 | 1 | 36 |
| n29a1 | 1.047 | 1.017 | 1 | 36 |
| n31a0 | 1.257 | 1.252 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
