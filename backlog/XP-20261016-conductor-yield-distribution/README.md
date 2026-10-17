# Measure relation yield and rank growth across fresh seeds

Status: distributed handoff blocked pending review and execution-layer setup. Local candidate-code execution is authorized and tracked in the study ledger; no shared task, lease or funding action has been created.

## Question and decision

Is the apparent endpoint variation explained by the distribution of query counts and rank-growth events?

Repeated fixed-budget observations distinguish a persistent relation-yield difference from heavy-tailed search luck. A reproducible yield effect motivates algebraic investigation; an unstable effect shifts effort toward variance reduction and larger curve-level samples.

## Design to freeze

Use 20 gaps, ten fresh target/query seeds per endpoint and the generic recurrence on both endpoints (400 attempts). Freeze equal query budgets rather than a time cutoff, plus a separate safety timeout. Log every valid hit, rejected proportional row, rank increment and query count at each rank. Publish both budget-limited success probabilities and restricted mean queries through the common cap; retain all censored attempts. Bootstrap by gap, with seed observations nested inside each gap. Compare hit rate, independent-row rate and rank progression separately. Measure cost of event logging on a small frozen control before the campaign. Do not estimate a successful-only solve ratio.

## Prerequisites

Implement a shared query-budget driver, rank-event logging, censor-aware aggregation and an independent event replay checker.

The native candidate driver and independent replay are implemented in research/conductor_gap_followups_20261016. The joint protocol and exact inputs are frozen before local launch. Distributed ready promotion still requires reviewed merged code and execution-layer setup.

## Cost and completion

Proposal only. Initial design has 400 IC or rho attempts as specified in README.md; pilot calibration is required. 18 core-hours is a provisional scheduling allowance based on the existing small-field attempts, with setup/replay overhead; it is not a measured bound.

Done means publishing the full denominator, immutable inputs and code hashes, reservation/timing boundaries, raw records, independent verification, failed and censored attempts, and predeclared uncertainty analysis. Results must distinguish an arithmetic implementation advantage, a relation-yield effect and a claim about ECDLP complexity. All existing large-gap solve timings are 53-bit; the certified 96-bit codomain is not a solved 96-bit DLP.

Baseline: [completed report](../../research/conductor_gap_cm_orbits_20261015/REPORT.md). Existing D=-43 replication, full pricing, coordinate-seed transport, point-map and prime-field PDP proposals remain separate.

## Authorized execution tracking — 16 October 2026

The user requested execution of all five proposals. See the [joint frozen protocol](../../research/conductor_gap_followups_20261016/protocols/PROTOCOL.md) and [per-requirement status ledger](../../research/conductor_gap_followups_20261016/README.md). The implemented native campaign checkpoints and runs the joint ordered workload; invoke it once, rather than launching one copy per card. `audit-all` independently checks the complete sealed study. This card remains blocked for distributed handoff while local candidate-code validation and execution are tracked in the study. No cairn task, lease or funding action has been created. Earlier future-driver prerequisites above are superseded by this implemented candidate driver; readiness still requires protocol/code review and execution-layer setup.
