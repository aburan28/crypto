# IC candidate tournament: round-autolab-20260915

Decision: **retained — incumbent**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

Independent audit checked **1356 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. The measured native times remain diagnostics.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.07609**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.529e+05 | 1 | 13.41 | 1.989e+07 | 36/36 | reference |
| autolab | 4.929e+05 | 0.8915 | 11.95 | 1.773e+07 | 36/36 | engineering experiment |
| rho | 4.124e+04 | 0.07459 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.183e+05 | 1 | 13.14 | 2.093e+07 | 180/180 | reference |
| autolab | 4.582e+05 | 0.8841 | 11.62 | 1.85e+07 | 180/180 | engineering experiment |
| rho | 3.944e+04 | 0.07609 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.183e+05 | 1 | 13.14 | 2.093e+07 | 180/180 | reference |
| autolab | 4.582e+05 | 0.8841 | 11.62 | 1.85e+07 | 180/180 | engineering experiment |
| rho | 3.944e+04 | 0.07609 | 1 | unmeasured | 180/180 | reference |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.818 | 0.727 | 1 | 36 |
| n17a1 | 0.926 | 0.936 | 1 | 36 |
| n19a0 | 0.992 | 0.979 | 1 | 36 |
| n19a1 | 0.989 | 0.979 | 1 | 36 |
| n23a0 | 0.974 | 0.964 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../../../../../../../ic_candidate_tournament_20260915/OPERATIONS.md).

