# IC candidate tournament: round-0020

Decision: **promoted — both**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **1.038x**; candidate/incumbent ratio 0.9638, paired 95% interval [0.9440146210349667, 0.9821772733975138]. This is an engineering result on the tested workloads.

Independent audit checked **2646 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 124 fresh inputs over 8 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **1.509**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,821 | 1 | 0.6409 | 1.786e+05 | 18/18 | reference |
| both | 2,749 | 0.9745 | 0.6246 | 1.74e+05 | 18/18 | engineering experiment |
| rho | 4,402 | 1.56 | 1 | unmeasured | 18/18 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,909 | 1 | 0.6444 | 1.841e+05 | 54/54 | reference |
| both | 2,828 | 0.9722 | 0.6265 | 1.79e+05 | 54/54 | engineering experiment |
| rho | 4,514 | 1.552 | 1 | unmeasured | 54/54 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,144 | 1 | 0.6877 | 1.781e+05 | 372/372 | reference |
| both | 3,031 | 0.9638 | 0.6629 | 1.717e+05 | 372/372 | engineering experiment |
| rho | 4,572 | 1.454 | 1 | unmeasured | 372/372 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,144 | 1 | 0.6877 | 1.781e+05 | 372/372 | reference |
| both | 3,031 | 0.9638 | 0.6629 | 1.717e+05 | 372/372 | engineering experiment |
| rho | 4,572 | 1.454 | 1 | unmeasured | 372/372 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **True**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.003 | 1 | reference | 0.9069 | 54/54 |
| both | 0.9659 | 0.9629 | [0.9278, 0.99] | 0.8732 | 54/54 |
| rho | 1.106 | 1.103 | [1.029, 1.163] | 1 | 54/54 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.004 | 1 | reference | 0.9107 | 372/372 |
| both | 0.9733 | 0.9698 | [0.9509, 0.9859] | 0.8832 | 372/372 |
| rho | 1.102 | 1.098 | [1.054, 1.139] | 1 | 372/372 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 0.9998 | 1 | reference | 0.905 | 372/372 |
| both | 0.9711 | 0.9713 | [0.9442, 0.9942] | 0.8791 | 372/372 |
| rho | 1.105 | 1.105 | [1.07, 1.137] | 1 | 372/372 |

Confirmation winner/rho: instructions 0.6629, CI [0.5916267293789641, 0.7425526192579226]; native time 0.8832, CI [0.8567616486430494, 0.9128743281580725].

Replay winner/rho: instructions 0.6629, CI [0.5916214187713567, 0.7425536526002788]; native time 0.8791, CI [0.8511557672296495, 0.9093046699854032].

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.727 | 0.636 | 1 | 36 |
| n17a1 | 0.745 | 0.727 | 1 | 36 |
| n19a0 | 1.224 | 1.210 | 1 | 36 |
| n19a1 | 0.874 | 0.864 | 1 | 36 |
| n23a0 | 0.652 | 0.649 | 1 | 36 |
| n23a1 | 1.033 | 1.032 | 1 | 120 |
| n29a1 | 1.106 | 1.076 | 1 | 36 |
| n31a0 | 0.909 | 0.901 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
