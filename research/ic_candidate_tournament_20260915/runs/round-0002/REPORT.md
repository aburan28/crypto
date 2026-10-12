# IC candidate tournament: round-0002

Decision: **promoted — batch16**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

The confirmed instruction speedup is **1.648x**; candidate/incumbent ratio 0.6070, paired 95% interval [0.5920526378557471, 0.6237890398557485]. This is an engineering result on the tested workloads.

Independent audit checked **1632 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. The measured native times remain diagnostics.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.07672**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 9.23e+05 | 1 | 22.15 | 3.32e+07 | 36/36 | reference |
| batch128 | 1.415e+06 | 1.533 | 33.97 | 5.091e+07 | 36/36 | engineering experiment |
| batch16 | 5.545e+05 | 0.6008 | 13.31 | 1.995e+07 | 36/36 | engineering experiment |
| batch32 | 6.784e+05 | 0.735 | 16.28 | 2.441e+07 | 36/36 | engineering experiment |
| dense | 9.221e+05 | 0.999 | 22.13 | 3.317e+07 | 36/36 | engineering experiment |
| excess8 | 9.217e+05 | 0.9986 | 22.12 | 3.316e+07 | 36/36 | engineering experiment |
| window8 | 7.858e+05 | 0.8514 | 18.86 | 2.827e+07 | 36/36 | engineering experiment |
| rho | 4.166e+04 | 0.04514 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 8.537e+05 | 1 | 21.48 | 3.447e+07 | 180/180 | reference |
| batch16 | 5.182e+05 | 0.607 | 13.03 | 2.092e+07 | 180/180 | engineering experiment |
| rho | 3.975e+04 | 0.04657 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 8.537e+05 | 1 | 21.48 | 3.447e+07 | 180/180 | reference |
| batch16 | 5.182e+05 | 0.607 | 13.03 | 2.092e+07 | 180/180 | engineering experiment |
| rho | 3.975e+04 | 0.04657 | 1 | unmeasured | 180/180 | reference |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.818 | 0.727 | 1 | 36 |
| n17a1 | 1.135 | 1.117 | 1 | 36 |
| n19a0 | 0.598 | 0.585 | 1 | 36 |
| n19a1 | 0.595 | 0.600 | 1 | 36 |
| n23a0 | 0.806 | 0.806 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../OPERATIONS.md).

The snapshot/build failure preceding this campaign is retained under `../round-0001`; it produced no performance measurements.
