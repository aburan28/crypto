# Draft cross-route accounting contract

This directory is a **non-authoritative WS7 draft**. It defines a common
receipt shape for comparing mathematical, protocol, implementation, physical,
provenance, supply-chain, and evidence-integrity routes without merging their
different cost units or evidence standards.

Nothing here is an authorization instrument, an attack result, a promoted
claim, or evidence that a route works. Every conforming receipt has
`authoritative=false`, `promotion_allowed=false`, and
`claim_transition_allowed=false`. A later authoritative ledger may consume a
reviewed successor, but these records cannot change ledger state.

## Contract

[`cross-route-cost-receipt.schema.json`](cross-route-cost-receipt.schema.json)
requires:

- exact audit-unit, subject, target-set, implementation, and evidence
  revisions;
- an exact SHA-256 revision binding to the WS5 route-registry seed, with the
  route-family and evidence-profile joins independently marked `resolved` or
  `unresolved`;
- explicit authorization and fixture class, with live targets, secret material,
  and production exploitation disabled in this draft;
- exact reuse of the plan's `execution_outcome` vocabulary, with scientific
  progress recorded separately in the plan's `requirement_status` vocabulary;
- line-item phases and separate cold-start, one-target online, amortized-total,
  and amortized-per-target boundaries;
- unsuccessful attempts, timeouts, OOMs, invalid results, budget exhaustion,
  and retained failure artifacts;
- candidate generation, uniqueness, and duplicate accounting;
- peak and persistent memory observations, including explicit unknowns;
- a matched-control declaration that forbids best-of-run composition and
  records every remaining difference;
- measured, modeled, extrapolated, unknown, and not-applicable quantities as
  different states; and
- open obligations and claim limits.

Amortized values are always secondary. They must state the exact reuse group,
number of targets, setup allocation formula, reusable phase IDs, and per-target
phase IDs. They do not replace the cold-start or one-target online boundary.

JSON Schema validates structure, not arithmetic truth. A future WS7 semantic
checker must additionally recompute totals, reconcile attempts and duplicates,
verify artifact and route-registry hashes, confirm every resolved route/profile
join exists in the pinned registry, compare identity fields across controls,
reject phase overlap, and prove that each amortized allocation follows its
recorded formula.
Until that checker and independent review exist, WS7 remains open.

The mock `route.route_uid` values are local fixture identifiers. A receipt does
not acquire a WS5 route-family or evidence-profile identity from similar names
or route classes. When an exact reviewed join has not been established, the
corresponding `route_registry_binding` join is `unresolved`, its `id` is null,
and the reason remains visible. The registry's schema version and `as_of` date
are descriptive; `registry_revision_sha256` is the revision identity.

A `matched_control.status` of `matched` is reserved for a separately present,
schema-valid control receipt whose identities and complete-run conditions can
be checked. Shape-only examples with no companion receipt are `unmatched`, use
a null `control_receipt_uid`, and identify the absent control as an open
obligation.

`execution.outcome` and `result.execution_outcome` use the plan vocabulary at
`docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json#/vocabularies/execution_outcome`
exactly. The former local value `partial` is not silently treated as a plan
execution outcome. For the implementation fixture, `complete` means only that
the bounded invented fixture lifecycle reached its recorded end; its failed
attempt remains charged, while `result.requirement_status=partial`,
`claim_status=unknown`, and the open obligations preserve the nonterminal
scientific state. Completion of a receipt never means completion of a claim.

## Fixtures

The positive records are synthetic shape examples. Their commands, artifacts,
measurements, revisions, and results are invented and have no scientific
meaning.

| Fixture | Route class | Contract surface exercised |
|:--|:--|:--|
| `valid-mathematical-transfer.json` | `mathematical_transfer` | reusable construction, per-target transfer, verification, four-target amortization, and an explicitly unmatched shape-only control |
| `valid-protocol-misuse.json` | `protocol_misuse` | transcript ingestion, target-specific protocol checking, and an explicitly unmatched control |
| `valid-implementation-side-channel.json` | `implementation_side_channel` | failed-attempt charging, duplicate accounting, and an explicit unknown memory quantity |
| `valid-evidence-integrity.json` | `evidence_integrity` | artifact identity reconciliation and a control marked not applicable |

The negative records each mutate an otherwise conforming shape. They must be
rejected at the named path:

| Fixture | Expected structural failure |
|:--|:--|
| `invalid-promotion-attempt.json` | `/promotion_allowed` is true and `/result/claim_status` is `supported` |
| `invalid-hidden-failures.json` | mandatory `/failure_accounting` is absent |
| `invalid-amortized-without-denominator.json` | enabled `/accounting_boundaries/amortized` has a null `target_count` |
| `invalid-false-matched-control.json` | a `matched` control has `same_resource_cap=false` |

Syntax-only check (does not establish schema conformance):

```bash
find docs/ecdlp-audit/accounting -name '*.json' -print0 | xargs -0 -n1 jq empty
```

Cross-file pin and vocabulary checks:

```bash
ws7_route_revision="sha256:$(sha256sum \
  docs/ecdlp-audit/routes/route-registry.seed.json | awk '{print $1}')"
jq -s -e --arg revision "$ws7_route_revision" \
  'all(.[]; .route_registry_binding.registry_revision_sha256 == $revision)' \
  docs/ecdlp-audit/accounting/fixtures/*.json

jq -e --slurpfile plan docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json \
  '."$defs".execution.properties.outcome.enum == $plan[0].vocabularies.execution_outcome
   and ."$defs".result.properties.execution_outcome.enum == $plan[0].vocabularies.execution_outcome' \
  docs/ecdlp-audit/accounting/cross-route-cost-receipt.schema.json
```

These checks establish only that the fixtures name the current pinned file and
that the two outcome enumerations agree. They do not establish a route/profile
join or the truth of a receipt.

Draft 2020-12 structural check, including `date-time` formats (Ajv packages are
version-pinned for reproducibility):

```bash
npm exec --yes \
  --package=ajv-cli@5.0.0 \
  --package=ajv-formats@3.0.1 -- \
  ajv validate --spec=draft2020 --strict=true --all-errors \
    -c ajv-formats \
    -s docs/ecdlp-audit/accounting/cross-route-cost-receipt.schema.json \
    -d 'docs/ecdlp-audit/accounting/fixtures/valid-*.json'
```

Run the same command once per `invalid-*.json` fixture and require a nonzero
exit status. Acceptance of any negative fixture is a contract-test failure.
Schema acceptance establishes only structural conformance. In particular, it
does not establish that totals, identities, attempt counts, hashes, control
matching, or scientific statements are true.
