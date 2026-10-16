# Audit orbit geometry and factor-base overlap

Status: distributed handoff blocked pending review and execution-layer setup. Local candidate-code execution is authorized and tracked in the study ledger; no shared task, lease or funding action has been created.

## Question and decision

Do short scalar-orbit periods or signed collisions explain unusual relation yields?

The linked factor base may have different effective size when scalar powers overlap or repeat. Detecting this mechanism would explain solver variation without attributing it directly to the conductor, and would guide seed and orbit-length choices.

## Design to freeze

On all 100 existing gaps, compute and certify the multiplicative order of the subgroup eigenvalue lambda modulo r when factorization permits; otherwise publish proven bounds. Independently check signed orbit collisions, overlaps between seed orbits, distinct pair sums, and retained-relation dependencies at lengths 256, 512, 1024 and 2048. Crater CM and crater scalar must give identical records; compare floor recurrence with floor scalar as well. Report relation yield per distinct signed base point and per table entry. Start with an arithmetic-only audit, then perform at most 40 fixed-budget solves selected by a predeclared overlap criterion. Treat any selected subset as exploratory and confirm on held-out gaps. Distinguish the fixed four-column linear system from a factor base of independent logarithm columns.

## Prerequisites

Implement exact orbit-order/collision certificates and independent overlap checking before selecting any follow-up subset.

The native candidate driver and independent replay are implemented in research/conductor_gap_followups_20261016. The joint protocol and exact inputs are frozen before local launch. Distributed ready promotion still requires reviewed merged code and execution-layer setup.

## Cost and completion

Proposal only. Initial design has 40 IC or rho attempts as specified in README.md; pilot calibration is required. 6 core-hours is a provisional scheduling allowance based on the existing small-field attempts, with setup/replay overhead; it is not a measured bound.

Done means publishing the full denominator, immutable inputs and code hashes, reservation/timing boundaries, raw records, independent verification, failed and censored attempts, and predeclared uncertainty analysis. Results must distinguish an arithmetic implementation advantage, a relation-yield effect and a claim about ECDLP complexity. All existing large-gap solve timings are 53-bit; the certified 96-bit codomain is not a solved 96-bit DLP.

Baseline: [completed report](../../research/conductor_gap_cm_orbits_20261015/REPORT.md). Existing D=-43 replication, full pricing, coordinate-seed transport, point-map and prime-field PDP proposals remain separate.

## Authorized execution tracking — 16 October 2026

The user requested execution of all five proposals. See the [joint frozen protocol](../../research/conductor_gap_followups_20261016/protocols/PROTOCOL.md) and [per-requirement status ledger](../../research/conductor_gap_followups_20261016/README.md). The implemented native campaign checkpoints and runs the joint ordered workload; invoke it once, rather than launching one copy per card. `audit-all` independently checks the complete sealed study. This card remains blocked for distributed handoff while local candidate-code validation and execution are tracked in the study. No cairn task, lease or funding action has been created. Earlier future-driver prerequisites above are superseded by this implemented candidate driver; readiness still requires protocol/code review and execution-layer setup.
