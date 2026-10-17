# Test whether rho endpoint variation comes from the walk partition

Status: distributed handoff blocked pending review and execution-layer setup. Local candidate-code execution is authorized and tracked in the study ledger; no shared task, lease or funding action has been created.

## Question and decision

Does the endpoint ratio in generic rho persist across independently seeded walk partitions and walks?

Rho provides a generic control for the observed IC results, but its coordinate-based partition can introduce model and seed effects. Separating group-operation counts from arithmetic time determines whether any residual endpoint timing effect merits further study.

## Design to freeze

Use 20 gaps, five public targets per endpoint and two independently domain-separated partition seeds (400 attempts), with balanced interleaving and an identical algorithm parameterization. First freeze one rho implementation; add a second independently checked walk only as a separately priced replication. Count additions, doublings, canonicalization/hash work, collisions and restarts, along with full cold CPU and verified recovery. Add transported isomorphic models on a frozen subset to measure coordinate sensitivity. Across endpoints use matched seed labels, not a claim of identical walks. Analyze operations normalized by sqrt(r) separately from CPU per operation, with whole-gap bootstrap and all failures retained. Do not report an IC/rho speedup until both methods have common complete cost accounting.

## Prerequisites

Implement complete rho operation instrumentation and a balanced partition-seed driver with independent recovery replay.

The native candidate driver and independent replay are implemented in research/conductor_gap_followups_20261016. The joint protocol and exact inputs are frozen before local launch. Distributed ready promotion still requires reviewed merged code and execution-layer setup.

## Cost and completion

Proposal only. Initial design has 400 IC or rho attempts as specified in README.md; pilot calibration is required. 6 core-hours is a provisional scheduling allowance based on the existing small-field attempts, with setup/replay overhead; it is not a measured bound.

Done means publishing the full denominator, immutable inputs and code hashes, reservation/timing boundaries, raw records, independent verification, failed and censored attempts, and predeclared uncertainty analysis. Results must distinguish an arithmetic implementation advantage, a relation-yield effect and a claim about ECDLP complexity. All existing large-gap solve timings are 53-bit; the certified 96-bit codomain is not a solved 96-bit DLP.

Baseline: [completed report](../../research/conductor_gap_cm_orbits_20261015/REPORT.md). Existing D=-43 replication, full pricing, coordinate-seed transport, point-map and prime-field PDP proposals remain separate.

## Authorized execution tracking — 16 October 2026

The user requested execution of all five proposals. See the [joint frozen protocol](../../research/conductor_gap_followups_20261016/protocols/PROTOCOL.md) and [per-requirement status ledger](../../research/conductor_gap_followups_20261016/README.md). The implemented native campaign checkpoints and runs the joint ordered workload; invoke it once, rather than launching one copy per card. `audit-all` independently checks the complete sealed study. This card remains blocked for distributed handoff while local candidate-code validation and execution are tracked in the study. No cairn task, lease or funding action has been created. Earlier future-driver prerequisites above are superseded by this implemented candidate driver; readiness still requires protocol/code review and execution-layer setup.
