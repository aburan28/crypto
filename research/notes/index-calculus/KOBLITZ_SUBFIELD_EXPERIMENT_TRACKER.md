# Koblitz/Subfield Index-Calculus Experiment Tracker

Status: implementation-facing experiment intake. No attack or rho-beating claim.

This repository owns executable algorithm implementations and measured boundaries. The companion research portfolio lives in `aburan28/crypto-autoresearcher` PR #1431.

## Why this is tracked here

The existing code already contains Koblitz factor-base search, point decomposition, Frobenius-aware relation machinery, Semaev experiments, and end-to-end accounting. The purpose of this tracker is to ensure research proposals become implementation tasks against those components rather than remaining disconnected notes.

The historical motivation is legitimate but not evidence of a new attack: extension/subfield structure has enabled faster-than-square-root ECDLP algorithms in special settings, and Frobenius must also be charged to the generic rho comparator.

## Implementation queue

| Priority | Thread | Code target | Required experiment |
|---|---|---|---|
| P0 | direct orbit-coordinate decomposition | `koblitz_index_calculus`, IC experiment runner | raw vs representative+relative-shift exact-equivalence and relation-cost benchmark |
| P0 | Frobenius factor telemetry | IC runner + reports | raw points/effective columns/solver/relation/LA/e2e reduction factors |
| P0 | decomposition likelihood | candidate/feature pipeline | calibrated p(decompose) predictor; held-out curves; Brier/log-loss; solver calls saved |
| P1 | Frobenius-aware F4/F5/SAT order | Groebner/SAT backends | identical ideal/CNF, orbit block orders, peak degree/Macaulay/conflicts |
| P1 | orbit-aware factor-base search | `koblitz_factor_base_search` | optimize verified independent relations/sec at matched effective columns |
| P1 | hybrid SAT -> Groebner | solver pipeline | SAT orbit/Boolean choices then algebraic residual vs pure controls |
| P1 | trace-zero intersection bases | factor-base candidates | trace-zero x Frobenius-stable bases vs matched random/orbit controls |
| P1 | partial/large-orbit relations | relation collector | one/two-large-orbit graph, recombination, memory and full collection cost |
| P2 | tau-adic candidate generation | Koblitz scalar/candidate generator | uniform vs matched-weight vs sparse tau-NAF relation targets |
| P2 | generic endomorphism group interface | orbit abstraction | Frobenius plus GLV/GLS-compatible fixture under same verifier |
| P2 | partial Weil descent | extension-field experiment code | no/full/intermediate descent sweep with complete solver+FB accounting |
| P2 | invariant telemetry | solver instrumentation | Hilbert/Macaulay/syzygy/orbit/trajectory features predicting held-out solve cost |

## Common acceptance contract

Every performance experiment reports the repository-standard end-to-end boundary table. At minimum record:

- curve, field representation, subgroup and seed;
- raw factor-base cardinality and effective Frobenius columns;
- decomposition attempts/successes and independently verified relations;
- relation rank growth and duplicate rate modulo Frobenius;
- setup, encoding, solver, extraction, verification and linear-algebra costs;
- peak memory and solver diagnostics;
- `S = total_operations / sqrt(n)` when a complete solve is available;
- ratio to the strongest applicable automorphism-accelerated Pollard-rho baseline.

A stage improvement is labelled a stage diagnostic until full cost is known.

## Decomposition-likelihood contract

The predictor estimates `P(decomposable | target-independent features, frozen FB, arity, cohort)`. Initial features: trace/subtrace where defined, orbit length, stabilizer, normal-basis weight/correlation, tau-NAF weight, factor-base projections, and cheap early Macaulay/Semaev features.

Report Brier score, log loss, calibration error, solver calls avoided, false negatives, model/feature overhead, and net seconds per independent relation. Ranking is allowed immediately; hard rejection requires zero false negatives on the frozen validation corpus.

## Follow-through rule

Each row must move through:

`specified -> implemented -> tiny_verified -> benchmarked -> held_out -> e2e -> reviewed -> promoted|retired`.

If a row has no new artifact for seven days and no explicit blocker, it is stale. Prefer advancing stale P0/P1 rows over adding adjacent speculative variants.

## Source context

Koblitz, Koblitz, and Menezes summarize that Frobenius conjugacy can accelerate generic attacks on curves defined over a subfield and that Weil descent / direct summation-polynomial index calculus can beat square-root attacks in special extension-field regimes. This tracker treats those observations as motivation and requires modern matched controls rather than assuming transfer to prime-degree Koblitz targets.
