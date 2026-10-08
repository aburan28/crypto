# IC candidate tournament: round-0007-single-policy

Decision: **promoted — orbits2**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Each complete cold job recovers **1 target(s)**. Tables report total job cost; setup is charged once to that job. Workloads with different target counts are separate panels.

The confirmed instruction speedup is **2.665x**; candidate/incumbent ratio 0.3753, paired 95% interval [0.32387886079328404, 0.45379932682481017]. This is an engineering result on the tested workloads.

Independent audit checked **1560 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. Native timings use separate scope and confidence rules from the contract.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.868**. A value below one means rho costs less. This is not an extrapolated crossover.

## Smoke admission (all candidates)

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.339e+05 | 1 | 3.199 | 4.818e+06 | 12/12 | reference |
| cold_context | 9.437e+04 | 0.7047 | 2.254 | 3.395e+06 | 12/12 | engineering experiment |
| cube_root | unmeasured | unmeasured | unmeasured | unmeasured | 6/12 | engineering experiment |
| orbits1 | unmeasured | unmeasured | unmeasured | unmeasured | 0/12 | engineering experiment |
| orbits2 | 4.761e+04 | 0.3555 | 1.137 | 6.626e+06 | 12/12 | engineering experiment |
| orbits3 | 5.413e+04 | 0.4042 | 1.293 | 5.022e+06 | 12/12 | engineering experiment |
| orbits4 | 6.044e+04 | 0.4513 | 1.443 | 4.206e+06 | 12/12 | engineering experiment |
| rho | 4.187e+04 | 0.3126 | 1 | unmeasured | 12/12 | reference |

Rejected arms remain in the frozen smoke receipts and are excluded from development. Missing verified workloads have no cost claim.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.315e+05 | 1 | 3.157 | 4.729e+06 | 36/36 | reference |
| cold_context | 9.364e+04 | 0.7123 | 2.249 | 3.369e+06 | 36/36 | engineering experiment |
| orbits2 | 4.789e+04 | 0.3643 | 1.15 | 6.666e+06 | 36/36 | engineering experiment |
| orbits3 | 5.268e+04 | 0.4007 | 1.265 | 4.888e+06 | 36/36 | engineering experiment |
| orbits4 | 6.067e+04 | 0.4615 | 1.457 | 4.222e+06 | 36/36 | engineering experiment |
| rho | 4.164e+04 | 0.3168 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.213e+05 | 1 | 3.07 | 4.899e+06 | 180/180 | reference |
| orbits2 | 4.553e+04 | 0.3753 | 1.152 | 7.16e+06 | 180/180 | engineering experiment |
| rho | 3.953e+04 | 0.3258 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 1.213e+05 | 1 | 3.07 | 4.899e+06 | 180/180 | reference |
| orbits2 | 4.553e+04 | 0.3753 | 1.152 | 7.16e+06 | 180/180 | engineering experiment |
| rho | 3.953e+04 | 0.3258 | 1 | unmeasured | 180/180 | reference |

## Native runtime and rho parity

Complete cold process wall time, measured with blocking reap and an independent watchdog. Paired arm order is randomized; the same CPU and resource limits apply. These timings include process startup and reporting. Prior polling-based timings are not mixed with this protocol.

Parity verdict: **False**. Both candidate/rho upper paired 95% limits and every cell ratio <= 1.10, instructions and native process wall, confirmation and replay.

### Development native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.183 | 1 | reference | 1.539 | 36/36 |
| cold_context | 5.435 | 0.879 | [0.8576, 0.9051] | 1.353 | 36/36 |
| orbits2 | 4.357 | 0.7046 | [0.6717, 0.7493] | 1.084 | 36/36 |
| orbits3 | 4.439 | 0.7179 | [0.687, 0.7658] | 1.105 | 36/36 |
| orbits4 | 4.659 | 0.7536 | [0.7183, 0.8049] | 1.16 | 36/36 |
| rho | 4.018 | 0.6499 | [0.6215, 0.689] | 1 | 36/36 |

### Confirmation native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.186 | 1 | reference | 1.552 | 180/180 |
| orbits2 | 4.386 | 0.7091 | [0.6854, 0.7354] | 1.101 | 180/180 |
| rho | 3.985 | 0.6441 | [0.6306, 0.6579] | 1 | 180/180 |

### Replay native time

| Variant | Cold process ms | Time / incumbent | Paired 95% interval | Time / rho | Verified |
|---|---:|---:|---|---:|---:|
| incumbent | 6.392 | 1 | reference | 1.526 | 180/180 |
| orbits2 | 4.581 | 0.7167 | [0.6951, 0.7413] | 1.093 | 180/180 |
| rho | 4.19 | 0.6555 | [0.6405, 0.6706] | 1 | 180/180 |

Confirmation winner/rho: instructions 1.1520, CI [1.0569038248993683, 1.2816089650364448]; native time 1.1008, CI [1.0710875049873299, 1.1469251895028603].

Replay winner/rho: instructions 1.1520, CI [1.056898243454399, 1.2816001972444366]; native time 1.0934, CI [1.0649139897111202, 1.1419276756327945].

## Factor-base policy boundary

This separately declared panel compares different base policies on identical public ECDLP targets. Support is stable within each arm/case. For each arm, m=3 and B signed points give at most binomial(B+2,3) target images; the rank floor also changes with the column count. Ratios to changed floors are not evidence of an algorithmic advance.

| Stage | Arm | Cell | Signed base size B |
|---|---|---|---:|
| confirmation | incumbent | n13a0 | [182] |
| confirmation | incumbent | n17a1 | [272] |
| confirmation | incumbent | n19a0 | [304] |
| confirmation | incumbent | n19a1 | [304] |
| confirmation | incumbent | n23a0 | [368] |
| confirmation | orbits2 | n13a0 | [52] |
| confirmation | orbits2 | n17a1 | [68] |
| confirmation | orbits2 | n19a0 | [76] |
| confirmation | orbits2 | n19a1 | [76] |
| confirmation | orbits2 | n23a0 | [92] |
| development | cold_context | n13a0 | [182] |
| development | cold_context | n17a1 | [272] |
| development | cold_context | n19a0 | [304] |
| development | cold_context | n23a0 | [368] |
| development | incumbent | n13a0 | [182] |
| development | incumbent | n17a1 | [272] |
| development | incumbent | n19a0 | [304] |
| development | incumbent | n23a0 | [368] |
| development | orbits2 | n13a0 | [52] |
| development | orbits2 | n17a1 | [68] |
| development | orbits2 | n19a0 | [76] |
| development | orbits2 | n23a0 | [92] |
| development | orbits3 | n13a0 | [78] |
| development | orbits3 | n17a1 | [102] |
| development | orbits3 | n19a0 | [114] |
| development | orbits3 | n23a0 | [138] |
| development | orbits4 | n13a0 | [104] |
| development | orbits4 | n17a1 | [136] |
| development | orbits4 | n19a0 | [152] |
| development | orbits4 | n23a0 | [184] |

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
