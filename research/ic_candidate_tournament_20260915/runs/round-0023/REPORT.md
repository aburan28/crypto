# IC candidate tournament: round-0023

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

Independent audit checked **3648 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; this round's confirmation includes 172 fresh inputs over 12 curve cells, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.8871**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,860 | 1 | 1.017 | 5.111e+05 | 24/24 | reference |
| scaled | 2,367 | 0.8275 | 0.8412 | 3.381e+05 | 24/24 | engineering experiment |
| rho | 2,814 | 0.9837 | 1 | unmeasured | 24/24 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 2,842 | 1 | 0.9488 | 5.078e+05 | 72/72 | reference |
| scaled | 2,639 | 0.9287 | 0.8812 | 3.77e+05 | 72/72 | engineering experiment |
| rho | 2,995 | 1.054 | 1 | unmeasured | 72/72 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,729 | 1 | 1.127 | 1.063e+06 | 516/516 | reference |
| scaled | 2,976 | 0.7982 | 0.8998 | 6.086e+05 | 516/516 | engineering experiment |
| rho | 3,308 | 0.8871 | 1 | unmeasured | 516/516 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3,729 | 1 | 1.127 | 1.063e+06 | 516/516 | reference |
| scaled | 2,976 | 0.7982 | 0.8998 | 6.086e+05 | 516/516 | engineering experiment |
| rho | 3,308 | 0.8871 | 1 | unmeasured | 516/516 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 1.495 | 1 | reference | 1.164 | 72/72 |
| scaled | 1.439 | 0.9623 | [0.8208, 1.068] | 1.12 | 72/72 |
| rho | 1.284 | 0.8591 | [0.5423, 1.17] | 1 | 72/72 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.798 | 1 | reference | 1.331 | 516/516 |
| scaled | 2.306 | 0.8244 | [0.6781, 0.969] | 1.097 | 516/516 |
| rho | 2.102 | 0.7512 | [0.5005, 1.041] | 1 | 516/516 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 2.756 | 1 | reference | 1.312 | 516/516 |
| scaled | 2.295 | 0.8327 | [0.6881, 0.9775] | 1.092 | 516/516 |
| rho | 2.101 | 0.7624 | [0.5076, 1.061] | 1 | 516/516 |

Confirmation winner/rho: instructions 1.1272, CI [0.7443018484463897, 1.8530575243127165]; native time 1.3312, CI [0.962242761854845, 1.9984817716813714].

Replay winner/rho: instructions 1.1272, CI [0.7442932479304686, 1.8530614886916985]; native time 1.3117, CI [0.9425438045398311, 1.973960969895372].

## Factor-base policy boundary

This separately declared panel compares different base policies on identical public ECDLP targets. Support is stable within each arm/case. For each arm, m=3 and B signed points give at most binomial(B+2,3) target images; the rank floor also changes with the column count. Ratios to changed floors are not evidence of an algorithmic advance.

| Stage | Arm | Cell | Signed base size B |
|---|---|---|---:|
| confirmation | incumbent | n13a0 | [182] |
| confirmation | incumbent | n17a1 | [272] |
| confirmation | incumbent | n19a0 | [304] |
| confirmation | incumbent | n19a1 | [304] |
| confirmation | incumbent | n23a0 | [368] |
| confirmation | incumbent | n23a1 | [368] |
| confirmation | incumbent | n29a1 | [464] |
| confirmation | incumbent | n31a0 | [496] |
| confirmation | incumbent | n37a0 | [592] |
| confirmation | incumbent | n43a1 | [688] |
| confirmation | incumbent | n59a0 | [944] |
| confirmation | incumbent | n61a1 | [976] |
| confirmation | scaled | n13a0 | [182] |
| confirmation | scaled | n17a1 | [272] |
| confirmation | scaled | n19a0 | [304] |
| confirmation | scaled | n19a1 | [304] |
| confirmation | scaled | n23a0 | [368] |
| confirmation | scaled | n23a1 | [368] |
| confirmation | scaled | n29a1 | [464] |
| confirmation | scaled | n31a0 | [496] |
| confirmation | scaled | n37a0 | [1184] |
| confirmation | scaled | n43a1 | [2064] |
| confirmation | scaled | n59a0 | [2832] |
| confirmation | scaled | n61a1 | [2928] |
| development | incumbent | n13a0 | [182] |
| development | incumbent | n17a1 | [272] |
| development | incumbent | n19a0 | [304] |
| development | incumbent | n23a0 | [368] |
| development | incumbent | n23a1 | [368] |
| development | incumbent | n31a0 | [496] |
| development | incumbent | n37a0 | [592] |
| development | incumbent | n43a1 | [688] |
| development | scaled | n13a0 | [182] |
| development | scaled | n17a1 | [272] |
| development | scaled | n19a0 | [304] |
| development | scaled | n23a0 | [368] |
| development | scaled | n23a1 | [368] |
| development | scaled | n31a0 | [496] |
| development | scaled | n37a0 | [1184] |
| development | scaled | n43a1 | [2064] |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.591 | 1.500 | 1 | 36 |
| n17a1 | 1.035 | 1.045 | 1 | 36 |
| n19a0 | 1.115 | 1.101 | 1 | 36 |
| n19a1 | 1.042 | 1.032 | 1 | 36 |
| n23a0 | 0.706 | 0.703 | 1 | 36 |
| n23a1 | 0.918 | 0.919 | 1 | 120 |
| n29a1 | 0.678 | 0.649 | 1 | 36 |
| n31a0 | 0.775 | 0.770 | 1 | 36 |
| n37a0 | 1.195 | 1.190 | 1 | 36 |
| n43a1 | 0.768 | 0.769 | 1 | 36 |
| n59a0 | 1.045 | 1.044 | 1 | 36 |
| n61a1 | 0.774 | 0.774 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The earlier snapshot/build attempt is retained under `../round-0001`.
