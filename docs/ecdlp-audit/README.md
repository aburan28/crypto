# ECDLP audit control plane

This directory holds the machine-readable planning layer for the proposed
[comprehensive elliptic-curve private-key recovery audit](../../research/notes/ecdlp-general/RESEARCH_COMPREHENSIVE_ECDLP_AUDIT_PLAN.md).
The precise scope is authorized research into unauthorized recovery routes;
most protocol and implementation incidents do not solve the original ECDLP and
must remain separately classified.

## Current artifacts

| Artifact | Role | Authority |
|:--|:--|:--|
| [`ECDLP_AUDIT_PLAN.mock.json`](ECDLP_AUDIT_PLAN.mock.json) | Exact eight-workstream plan, dependencies, acceptance gates, taxonomy skeleton, and open obligations | **non-authoritative mock**; cannot promote claims |
| [`ecdlp-audit-plan.schema.json`](ecdlp-audit-plan.schema.json) | Draft 2020-12 structural contract | validates the mock only |
| [`verify_plan.py`](verify_plan.py) | Semantic checks for identity, counts, references, dependencies, source hashes, and fail-closed status | planning validator, not a scientific verifier |
| [`program-dependencies.svg`](program-dependencies.svg) | Dependency and promotion-boundary diagram | visual twin of the plan |

Run the planning checks from the repository root:

```bash
python3 docs/ecdlp-audit/verify_plan.py
```

The mock intentionally has:

- `authoritative=false` and `promotion_allowed=false`;
- no authoritative subjects or assessments;
- an unfrozen universal denominator;
- all 20 execution gates open and explicit; and
- five top-level taxonomy skeletons whose detailed route population remains a
  WS5/WS6 deliverable.

It must not be renamed to `master-ledger.json` or treated as that ledger. WS3
requires reviewed identity rules, a frozen population manifest, an additive
generator, an independent semantic verifier, migration receipts, and
governance approval before an authoritative ledger exists.

## Status semantics

Claim status is one of `supported`, `refuted`, `unknown`, or
`not_applicable`. Evidence maturity E0–E7 is separate. An explicit `unknown`
counts as a classified cell but not as scientific support, safety, or
completion. A failed bounded search supports only its exact bounded result.

Executable work defaults to passive public evidence and public
synthetic/known-answer fixtures. No live-target execution, secret material, or
production exploitation is authorized by these files.
