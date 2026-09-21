# Subfield curves with dimensions 6, 7, and 8

The solver now handles binary curves with non-F2 coefficients in F4, evaluated over GF(2^18) and GF(2^30). The dimension-eight supplemental runs completed 8/8 exact enumerations within 60 seconds; maximum 20.723362 seconds. This extends the implementation and measures a larger-base obstacle. It does not establish a calibrated Semaev or ECDLP speedup.

Frozen source and proofs: [subfield_01/README.md](subfield_01/README.md). Evidence: [raw.jsonl](subfield_01/raw.jsonl), [summary.json](subfield_01/summary.json), [audit.json](subfield_01/audit.json), and [contract.json](subfield_01/contract.json).

## Matched five-second comparison

Every table entry counts a verified first relation or a proved complete enumeration (including zero relations), as appropriate for its mode. A partial enumeration is a timeout. Budgets cover cold setup, search, extraction and verification. CryptoMiniSat uses a soft time limit; overruns remain recorded and do not count as within-budget completions.

First relation / negative-target resolution:

| Variant | n18 d6 | n18 d7 | n18 d8 | n30 d6 | n30 d7 | n30 d8 | Class | S / rho / cost speedup |
|---|---|---|---|---|---|---|---|---|
| hybrid-image | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | 2/4 | engineering diagnostic | unmeasured |
| chained-s3 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | engineering diagnostic | unmeasured |
| s4-symmetric | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | engineering diagnostic | unmeasured |

Complete enumeration:

| Variant | n18 d6 | n18 d7 | n18 d8 | n30 d6 | n30 d7 | n30 d8 | Class | S / rho / cost speedup |
|---|---|---|---|---|---|---|---|---|
| hybrid-image | 4/4 | 4/4 | 0/4 | 4/4 | 2/4 | 0/4 | engineering diagnostic | unmeasured |
| chained-s3 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | engineering diagnostic | unmeasured |
| s4-symmetric | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | engineering diagnostic | unmeasured |

These are 144 matched trials on 24 curve/base/target instances. For each field the same four targets are used at all dimensions: one uniform and one known-decomposable target for each of the predeclared development and holdout seeds. The supported stratum was sampled from d6 triples and is not a natural-yield estimate. Each stratum has only two targets per field, so no timing confidence interval or population success claim is made.

## Larger bases and complete relation sets

The following expected counts come from independent group-law pair enumeration, not from the tested solvers. Columns keep the same targets as the dimension increases.

| Field bits | d | Signed base points | Development uniform | Holdout uniform | Development supported | Holdout supported |
|---|---|---|---|---|---|---|
| 18 | 6 | 72 | 0 | 0 | 2 | 2 |
| 18 | 7 | 136 | 0 | 0 | 4 | 4 |
| 18 | 8 | 254 | 10 | 12 | 11 | 11 |
| 30 | 6 | 68 | 0 | 0 | 1 | 1 |
| 30 | 7 | 142 | 0 | 0 | 1 | 1 |
| 30 | 8 | 280 | 0 | 0 | 1 | 1 |

Supplemental cold hybrid enumeration at d8, with a separate 60-second budget. No Semaev comparison uses these unmatched runs:

| Field bits | Stratum | Cohort | Status | Verified relations | Seconds | Field API operations |
|---|---|---|---|---|---|---|
| 18 | uniform | development | complete | 10 | 11.727383 | 5260557 |
| 18 | known_decomposable | development | complete | 11 | 11.283723 | 5256586 |
| 18 | uniform | holdout | complete | 12 | 11.511674 | 5256089 |
| 18 | known_decomposable | holdout | complete | 11 | 11.684988 | 5258950 |
| 30 | uniform | development | complete | 0 | 20.449559 | 6393956 |
| 30 | known_decomposable | development | complete | 1 | 20.232045 | 6392246 |
| 30 | uniform | holdout | complete | 0 | 20.723362 | 6394066 |
| 30 | known_decomposable | holdout | complete | 1 | 18.969908 | 6395700 |

## What is proved and what remains open

- The general-coefficient norm formula, quadratic reduction, and image support test are exact in the declared distinct-abscissa domain. The general-B symmetric S4 identity is checked symbolically against the resultant of two S3 polynomials.
- All 53 affine GF64 targets agree across the hybrid, both Semaev controls, signed-triple enumeration, and the pair oracle. All 7,582 candidate support checks agree with modular H|L_V. All 576 AS inputs over GF64/GF512 pass; the GF64 point lifting set matches all 4,096 possible coordinate pairs.
- The final audit replays 99 signed-point certificates and independently recomputes all 24 larger oracle sets. Source and instance hashes, matched coverage, fixed targets across dimensions, exclusive timers, and summary aggregates pass with zero validation failures.
- Curves over F4 are supported here; odd-characteristic subfield curves and F4-linear factor bases are not tested. The inherited coordinate Frobenius is not used for automorphism compression: only powers fixing the curve coefficients act on the same curve.

The branch search remains O(2^(2d)) before field costs. Increasing d moves the counting bound 8 binom(M/2,3)/(#E-1); extra relations are not evidence of beating a fixed counting boundary. This round is classified as a functionality extension and engineering diagnostic. All common-operation speedups, full-DLP S, rho ratios and floor ratios remain null. The SAT field counters omit SAT and encoding work, so comparing them to the hybrid counters would be invalid. Only field additions, multiplications, squarings and inversionCalls are instrumented. Inherited curve-call fields in the raw CountedField reports are zero placeholders, not measurements; the audit marks those call counts null and all aggregates exclude them.

The parent accounting/promotion contract remains unmet: this small stage panel does not replace the required solver regression, calibrated full-pipeline comparison, independent-curve holdouts or paired timing repetitions. No relation matrix or final scalar recovery is run. The user’s 20% all-cost goal remains open.

Next useful experiment: reduce the quadratic dependence on base size or safely amortize image setup over a declared target batch, then run the broader matched regression with calibrated Boolean/field accounting. The present data alone do not justify an exponent fit.
