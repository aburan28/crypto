# IC candidate tournament: round-20260915t205319z

Decision: **promoted — autolab**.

This round compares complete cold Koblitz ECDLP recovery in fixed-compiler amd64 user-space instruction reads (Ir).
Every admitted run includes setup, failed attempts, relation and log verification, scalar linear algebra, descent and final scalar verification.

The confirmed instruction speedup is **1.367x**; candidate/incumbent ratio 0.7313, paired 95% interval [0.7063051319526432, 0.7600346799570532]. This is an engineering result on the tested workloads.

Independent audit checked **1356 trial receipts**. Repetitions are grouped within fixtures; intervals resample curve cells and their targets.

The operation unit is explicit and unconverted to curve additions. Kernel/device work, profiler execution and external audit work are outside it. The measured native times remain diagnostics.

The floor is deliberately weak: the implemented K-column full-rank collector needs at least K relation-producing trials and at least K instructions. It does not establish a non-generic advance.

The ordinary WDSat corpus is inapplicable to this point-base API. This campaign uses matched complete-DLP fixtures and an independent point/rank checker; the pilot confirmation includes 60 fresh inputs, both IC arms, rho and three repetitions.

Observed rho/winner instruction ratio: **0.1048**. A value below one means rho costs less. This is not an extrapolated crossover.

## Development

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.722e+05 | 1 | 13.22 | 2.058e+07 | 36/36 | reference |
| autolab | 4.047e+05 | 0.7073 | 9.349 | 1.456e+07 | 36/36 | engineering experiment |
| rho | 4.329e+04 | 0.07566 | 1 | unmeasured | 36/36 | reference |

## Confirmation

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.154e+05 | 1 | 13.05 | 2.081e+07 | 180/180 | reference |
| autolab | 3.769e+05 | 0.7313 | 9.544 | 1.522e+07 | 180/180 | engineering experiment |
| rho | 3.949e+04 | 0.07662 | 1 | unmeasured | 180/180 | reference |

## Replay

| Variant | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified | Class |
|---|---:|---:|---:|---:|---:|---|
| incumbent | 5.154e+05 | 1 | 13.05 | 2.081e+07 | 180/180 | reference |
| autolab | 3.769e+05 | 0.7313 | 9.544 | 1.522e+07 | 180/180 | engineering experiment |
| rho | 3.949e+04 | 0.07662 | 1 | unmeasured | 180/180 | reference |

## Rho health

All reference outputs were independently verified. The following measured walk counts are compared with the signed-Frobenius model `sqrt(pi*r/(2*A))`, with `A=2*n`. This model is an expectation, not a bound on an individual randomized run.

| Cell | Median iterations / model | Median walk additions / model | Maximum restarts | Verified |
|---|---:|---:|---:|---:|
| n13a0 | 0.727 | 0.636 | 1 | 36 |
| n17a1 | 0.699 | 0.681 | 1 | 36 |
| n19a0 | 0.999 | 0.986 | 1 | 36 |
| n19a1 | 0.864 | 0.869 | 1 | 36 |
| n23a0 | 0.905 | 0.897 | 1 | 36 |

## Evidence

- [Frozen contract](contract.json), [unit and boundary](calibration.json), [candidates](candidates.json).
- [Decision](decision.json), [independent audit](audit.json), [reporting data](measurements.json).
- [Raw job receipts and profiles](runs/), [source manifest](source-manifest.json).
- [Operating commands and skills](../../../../../../../../ic_candidate_tournament_20260915/OPERATIONS.md).

