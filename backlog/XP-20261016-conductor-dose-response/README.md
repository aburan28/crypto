# Test a conductor-size trend within matched isogeny classes

Status: distributed handoff blocked pending review and execution-layer setup. Local candidate-code execution is authorized and tracked in the study ledger; no shared task, lease or funding action has been created.

## Question and decision

Does increasing certified conductor ratio change relation work after accounting for field, group and CM family?

A conductor-size trend is more specific than a single crater/floor contrast and would justify searching larger gaps. A null or family-dependent trend argues against treating the conductor ratio itself as an ECDLP hardness predictor.

## Design to freeze

First construct ordinary classes whose Frobenius conductor contains at least two distinct prime factors. Within each class certify several endpoint conductors, including 1, ell1, ell2 and ell1*ell2 where feasible, and use identical p, trace and prime subgroup order. Start with four classes, four strata and five fresh targets (80 attempts), using recurrence-based folded IC and matched rho controls. Vary isomorphism models separately. Repeat across at least two CM fields before interpreting a trend. Pre-register log conductor ratio and local prime factors/depth as separate predictors; publish class-level contrasts rather than pooling unmatched curves. Construction success and point transport are prerequisites, not assumed results. A single smooth-depth volcano does not establish a large-prime-barrier effect.

## Prerequisites

Find and independently certify multi-prime conductor classes with usable common prime-order groups and explicit strata.

The native candidate driver and independent replay are implemented in research/conductor_gap_followups_20261016. The joint protocol and exact inputs are frozen before local launch. Distributed ready promotion still requires reviewed merged code and execution-layer setup.

## Cost and completion

Proposal only. Initial design has 80 IC or rho attempts as specified in README.md; pilot calibration is required. No runtime estimate is defensible before constructing and certifying the strata.

Done means publishing the full denominator, immutable inputs and code hashes, reservation/timing boundaries, raw records, independent verification, failed and censored attempts, and predeclared uncertainty analysis. Results must distinguish an arithmetic implementation advantage, a relation-yield effect and a claim about ECDLP complexity. All existing large-gap solve timings are 53-bit; the certified 96-bit codomain is not a solved 96-bit DLP.

Baseline: [completed report](../../research/conductor_gap_cm_orbits_20261015/REPORT.md). Existing D=-43 replication, full pricing, coordinate-seed transport, point-map and prime-field PDP proposals remain separate.

## Authorized execution tracking — 16 October 2026

The user requested execution of all five proposals. See the [joint frozen protocol](../../research/conductor_gap_followups_20261016/protocols/PROTOCOL.md) and [per-requirement status ledger](../../research/conductor_gap_followups_20261016/README.md). The implemented native campaign checkpoints and runs the joint ordered workload; invoke it once, rather than launching one copy per card. `audit-all` independently checks the complete sealed study. This card remains blocked for distributed handoff while local candidate-code validation and execution are tracked in the study. No cairn task, lease or funding action has been created. Earlier future-driver prerequisites above are superseded by this implemented candidate driver; readiness still requires protocol/code review and execution-layer setup.
