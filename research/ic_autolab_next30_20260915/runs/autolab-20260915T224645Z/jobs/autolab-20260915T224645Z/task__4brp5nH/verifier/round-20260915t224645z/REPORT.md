# IC candidate tournament: round-20260915t224645z

Decision: **promoted — autolab**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

The confirmed instruction speedup is **1.913x**; candidate/incumbent ratio 0.5228, paired 95% interval [0.5144320153056702, 0.5321178192104234]. This is an engineering result on the tested workloads.

Independent audit checked **1356 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. The measured native times remain diagnostics.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.2008**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 4.079e+05 | 1 | 9.524 | 1.467e+07 | 36/36 | reference |
| autolab | 2.152e+05 | 0.5275 | 5.024 | 7.741e+06 | 36/36 | engineering experiment |
| rho | 4.283e+04 | 0.105 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.797e+05 | 1 | 9.527 | 1.533e+07 | 180/180 | reference |
| autolab | 1.985e+05 | 0.5228 | 4.981 | 8.014e+06 | 180/180 | engineering experiment |
| rho | 3.985e+04 | 0.105 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 3.797e+05 | 1 | 9.527 | 1.533e+07 | 180/180 | reference |
| autolab | 1.985e+05 | 0.5228 | 4.981 | 8.014e+06 | 180/180 | engineering experiment |
| rho | 3.985e+04 | 0.105 | 1 | unmeasured | 180/180 | reference |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 1.000 | 0.909 | 1 | 36 |
| n17a1 | 0.618 | 0.599 | 1 | 36 |
| n19a0 | 0.795 | 0.782 | 1 | 36 |
| n19a1 | 1.070 | 1.061 | 1 | 36 |
| n23a0 | 0.854 | 0.845 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../../../../../../../ic_candidate_tournament_20260915/OPERATIONS.md).

