# Koblitz and subfield index-calculus experiment queue

Companion tracking document for the executable campaign in `crypto-autoresearcher` PRs #1431 and #1436. This repo should retain the cryptanalytic questions and results even when orchestration code changes elsewhere.

## Research contract

Every idea must progress through `specified -> implemented -> tiny_verified -> benchmarked -> held_out -> e2e -> reviewed -> promoted|retired`. A thread with no artifact or blocker for 7 days is stale and should be advanced before generating a nearby variant.

North-star: **seconds per newly independent verified relation**, with complete end-to-end comparison against the applicable Frobenius/automorphism-accelerated Pollard-rho baseline.

## Experiment threads

1. **Orbit-coordinate decomposition:** solve directly in Frobenius representative + relative-shift coordinates; measure whether the decomposition problem itself shrinks.
2. **Frobenius-aware Gröbner/Macaulay:** symmetry-adapted ordering, invariant/orbit blocks, solving degree, matrix rows/NNZ, rank trajectory and wall time.
3. **Factor-base search:** optimize Frobenius-stable subspaces for relation yield / total solver+LA cost rather than cardinality.
4. **Trace-zero × Frobenius factor bases:** test low-dimensional intersections and compare relation density/solver complexity with matched controls.
5. **Representation search:** polynomial, normal, optimal-normal, tower and mixed bases scored by resulting Boolean/Semaev system complexity.
6. **Hybrid SAT → Gröbner:** SAT for orbit/membership/symmetry variables; Gröbner for nonlinear residue; compare full verified-relation cost.
7. **Large-orbit partial relations:** one/two out-of-base orbit representatives plus graph recombination and Frobenius canonicalization.
8. **Tau throughout pipeline:** stratify candidate generation by tau-NAF properties; separate generation speed from genuine decomposition-probability changes.
9. **Generic endomorphism quotienting:** extend Frobenius machinery to efficiently computable endomorphism groups such as GLS-style settings.
10. **Partial Weil descent:** for extension degree `n=ab`, descend to intermediate subfields and jointly optimize genus/system dimension/relation/solver cost.
11. **Decomposition-likelihood model:** predict decomposition probability from mathematical features; track Brier/log loss/calibration and solver-call avoidance. No hard rejection until zero false negatives on frozen held-out data.
12. **Algebraic system telemetry:** persist Hilbert/degree profile, solving/last-fall degree, Macaulay rank/NNZ, syzygy summaries, support and orbit stabilizers.

## Required shared metric: Frobenius factor

Record separately:
- raw factor-base points / effective orbit representatives;
- raw relation columns / effective columns;
- baseline solver work / candidate solver work;
- baseline relation work / candidate relation work;
- baseline LA dimension / candidate LA dimension;
- baseline e2e work / candidate e2e work.

Do not collapse these into one claimed speedup until the full pipeline is measured.

## Immediate order

1. exact orbit-coordinate encoding;
2. shared Frobenius telemetry;
3. Frobenius-aware Gröbner/Macaulay benchmark;
4. decomposition-likelihood corpus + calibrated baseline;
5. algebraic telemetry on every run;
6. remaining threads in dependency order;
7. held-out validation and e2e accounting.

## Scale policy

Tiny fields establish correctness only. Progress to larger binary extension degrees and explicitly include prime degrees relevant to Koblitz/subfield research, including 37, 53 and 83 when supported. Any scaling claim must show curves/instances/seeds and uncertainty, and cryptographic claims require an e2e comparison rather than solver-only extrapolation.

## Historical motivation

The ECC literature records that Frobenius conjugacy classes accelerate generic rho on curves defined over a smaller subfield, while Weil descent and direct index calculus have beaten square-root attacks for some composite extension-field families. It also records that those classical methods did not apply to the prime-degree binary fields commonly chosen for deployed ECC. These experiments ask whether orbit-aware decomposition, solver symmetry, and modern representations change any part of that boundary; they do not assume that they do.
