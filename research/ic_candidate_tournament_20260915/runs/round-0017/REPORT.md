# IC candidate tournament: round-0017

Decision: **promoted — orbits**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **1.202x**; candidate/incumbent ratio 0.8319, paired 95% interval [0.8109124358379497, 0.8512556961461358]. This is an engineering result on the tested workloads.

Independent audit checked **2340 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 96 fresh inputs over 8 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.492**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,383 | 1 | 0.8969 | 2.141e+05 | 18/18 | reference |
| orbits | 2,815 | 0.8321 | 0.7464 | 1.782e+05 | 18/18 | engineering experiment |
| orbits_rows | 2,828 | 0.8359 | 0.7497 | 1.79e+05 | 18/18 | engineering experiment |
| scan_io | 2,993 | 0.8847 | 0.7935 | 1.894e+05 | 18/18 | engineering experiment |
| rho | 3,771 | 1.115 | 1 | unmeasured | 18/18 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,462 | 1 | 0.9226 | 2.191e+05 | 54/54 | reference |
| orbits | 2,884 | 0.833 | 0.7685 | 1.825e+05 | 54/54 | engineering experiment |
| orbits_rows | 3,002 | 0.8671 | 0.8 | 1.9e+05 | 54/54 | engineering experiment |
| scan_io | 3,067 | 0.8859 | 0.8174 | 1.941e+05 | 54/54 | engineering experiment |
| rho | 3,752 | 1.084 | 1 | unmeasured | 54/54 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,821 | 1 | 0.8056 | 2.165e+05 | 288/288 | reference |
| orbits | 3,178 | 0.8319 | 0.6702 | 1.801e+05 | 288/288 | engineering experiment |
| rho | 4,743 | 1.241 | 1 | unmeasured | 288/288 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,821 | 1 | 0.8056 | 2.165e+05 | 288/288 | reference |
| orbits | 3,178 | 0.8319 | 0.6702 | 1.801e+05 | 288/288 | engineering experiment |
| rho | 4,743 | 1.241 | 1 | unmeasured | 288/288 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.258 | 1 | reference | 0.9927 | 54/54 |
| orbits | 1.21 | 0.962 | [0.9262, 0.9996] | 0.955 | 54/54 |
| orbits_rows | 1.204 | 0.9567 | [0.92, 0.9976] | 0.9498 | 54/54 |
| scan_io | 1.245 | 0.9892 | [0.9579, 1.025] | 0.982 | 54/54 |
| rho | 1.267 | 1.007 | [0.9426, 1.078] | 1 | 54/54 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.193 | 1 | reference | 0.9556 | 288/288 |
| orbits | 1.124 | 0.9417 | [0.9172, 0.964] | 0.8999 | 288/288 |
| rho | 1.249 | 1.046 | [1.009, 1.083] | 1 | 288/288 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.229 | 1 | reference | 0.9602 | 288/288 |
| orbits | 1.141 | 0.9287 | [0.9091, 0.9459] | 0.8918 | 288/288 |
| rho | 1.28 | 1.041 | [1.003, 1.077] | 1 | 288/288 |

Confirmation winner/rho: instructions 0.6702, CI [0.6050138439420483, 0.7442917796278755]; native time 0.8999, CI [0.8755174917468866, 0.9265554550177799].

Replay winner/rho: instructions 0.6702, CI [0.6050133930998162, 0.7442945697451894]; native time 0.8918, CI [0.8704528927757559, 0.9165578402001413].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.818 | 0.727 | 1 | 36 |
| n17a1 | 1.263 | 1.244 | 1 | 36 |
| n19a0 | 0.870 | 0.857 | 1 | 36 |
| n19a1 | 1.137 | 1.128 | 1 | 36 |
| n23a0 | 0.987 | 0.972 | 1 | 36 |
| n23a1 | 0.993 | 0.984 | 1 | 36 |
| n29a1 | 0.855 | 0.826 | 1 | 36 |
| n31a0 | 1.207 | 1.199 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
