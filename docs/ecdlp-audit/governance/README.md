# WS8 governance draft artifacts

This directory implements a **non-authoritative, synthetic draft** of the WS8
independent-governance contract. Nothing here can promote, demote, correct, or
release a scientific claim. Every machine record fixes
`authoritative=false`, `promotion_allowed=false`,
`self_promotion_allowed=false`, and `synthetic_only=true`.

## Requirement-to-evidence map

| Requested WS8 capability | Artifact | Current status |
|:--|:--|:--|
| Governance policy | [`POLICY.md`](POLICY.md) and [`governance-policy.draft.json`](governance-policy.draft.json) | implemented as a draft; not adopted or enforced |
| Strict structural schema | [`governance-draft.schema.json`](governance-draft.schema.json) | implemented; cross-record semantics still need a native validator |
| Role and ownership matrix | [`role-ownership-matrix.draft.json`](role-ownership-matrix.draft.json) | proposed only; CODEOWNERS/protected approvals are not installed |
| Mock promotion | [`fixtures/01-promotion-eligible-draft.json`](fixtures/01-promotion-eligible-draft.json) | positive eligibility drill; authoritative state remains unchanged |
| Mock rejection | [`fixtures/02-self-review-rejected.json`](fixtures/02-self-review-rejected.json) | negative self-review drill |
| Mock demotion | [`fixtures/03-demotion-on-refuted-gate.json`](fixtures/03-demotion-on-refuted-gate.json) | additive demotion drill |
| Mock additive correction | [`fixtures/04-additive-correction.json`](fixtures/04-additive-correction.json) | supersession drill; prior record is preserved |
| Mock staleness | [`fixtures/05-staleness-recheck.json`](fixtures/05-staleness-recheck.json) | revision-bound invalidation drill |
| Mock external replication | [`fixtures/06-external-replication-recorded.json`](fixtures/06-external-replication-recorded.json) | registry-only drill; no automatic promotion |

## Structural validation

Check JSON syntax from the repository root without using research Python:

```bash
jq empty docs/ecdlp-audit/governance/*.json \
  docs/ecdlp-audit/governance/fixtures/*.json
```

Validate every data record against Draft 2020-12 with AJV CLI 5:

```bash
for record in \
  docs/ecdlp-audit/governance/governance-policy.draft.json \
  docs/ecdlp-audit/governance/role-ownership-matrix.draft.json \
  docs/ecdlp-audit/governance/fixtures/*.json
do
  npx --yes ajv-cli@5 validate --spec=draft2020 --strict=true \
    -s docs/ecdlp-audit/governance/governance-draft.schema.json \
    -d "$record"
done
```

The schema rejects undeclared object properties and fixes all fail-closed draft
markers. Each request and decision must expose the proposed aggregate,
governance-schema, and dependency-lock digests. Each request's component list
must contain exactly one corresponding `proposed_aggregate`,
`governance_schema`, and `dependency_lock` entry. The schema also checks
scenario-specific transition/outcome pairs and requires correction, staleness,
demotion, and external-replication details where applicable.

Check that every synthetic decision repeats its request's three-part binding,
and that those values match the three required component entries:

```bash
for record in docs/ecdlp-audit/governance/fixtures/*.json
do
  jq -e '
    . as $record
    | .fixture.request as $request
    | ($request.exact_revisions
       | map({key: .component, value: .revision_digest})
       | from_entries) as $revisions
    | ($revisions.proposed_aggregate
       == $request.review_binding.aggregate_revision_digest)
      and ($revisions.governance_schema
       == $request.review_binding.governance_schema_revision_digest)
      and ($revisions.dependency_lock
       == $request.review_binding.dependency_revision_digest)
      and all($record.fixture.decisions[];
        .review_binding == $request.review_binding)
      and ($record.authoritative == false)
      and ($record.promotion_allowed == false)
      and ($record.self_promotion_allowed == false)
      and ($record.fixture.expected_outcome.state_effect
        == "none_draft_only")
  ' "$record"
done
```

## What structural validation does not establish

JSON Schema cannot establish that actor strings identify different humans,
that a reviewer is genuinely independent or qualified, that a digest matches
an external artifact, or that duplicated digest values are equal across a
request and its decisions. The `jq` check above establishes that equality only
for these committed synthetic fixtures; it is not the missing native semantic
validator. JSON Schema also cannot establish that all cross-document references
resolve, that conditional-role triggers were interpreted correctly, that E0-E7
gates are contiguous, or that storage is append-only. Those remain mandatory
checks for the still-open native semantic validator and operational governance
controls.

The positive fixture's expected disposition is only
`eligible_for_authoritative_review`; its state effect is still
`none_draft_only`. The negative fixtures are well-formed inputs describing
expected decisions, not malformed JSON. Passing the schema therefore closes
no evidence gate and satisfies neither WS8 acceptance gate by itself.
