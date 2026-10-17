# WS4 evidence-profile draft

This directory is a **draft, non-authoritative** implementation sketch for
Workstream 4 of the ECDLP audit plan. It cannot promote a scientific claim,
change an authoritative assessment, authorize execution, or satisfy WS4-G1 or
WS4-G2 by itself.

## Requirement-to-evidence map

| Requested WS4 capability | Draft artifact | Current status |
|:--|:--|:--|
| E0-E7 definitions | `evidence-profile-registry.draft.json` → `levels` | structurally specified; not governed or authoritative |
| R0-R5 independence definitions | registry → `independence_levels` and per-level requirements; fixture 01 → `event.independence_assessment` | structurally specified for the accepted E1 → E2 case, including method and reviewed bindings; cross-field identity/binding comparison and receipt resolution still require the semantic validator |
| Route-specific gates | registry → five lane profiles and gate catalog | structurally specified; route population and profile review remain open |
| Highest contiguous maturity | registry → `maturity_policy` | rule specified; no authoritative transition engine exists |
| Weakest-link composition | registry → `composition_policy` and fixtures 03-04 | rule and golden cases specified; cross-record enforcement remains open |
| Promotion and rejection | fixtures 01-03 | golden cases only |
| Demotion | fixture 04 | golden case only |
| Staleness | fixture 05 | golden case only |
| Additive and destructive correction | fixtures 06-07 | golden cases only |

## Files

- `ws4-evidence.schema.json` is a strict JSON Schema Draft 2020-12 contract
  for the registry and every transition fixture. Objects reject undeclared
  properties.
- `evidence-profile-registry.draft.json` defines common gates and distinct
  profiles for native ECDLP, mathematical transfer, parameter/provenance,
  protocol or implementation bypass, and evidence/supply-chain routes.
- `fixtures/*.json` are schema-valid golden inputs. A fixture whose expected
  decision is `reject` is not malformed JSON; it is a well-formed negative
  transition case that a future semantic validator must reject for the stated
  reason.

## Fixture catalog

| Fixture | Expected semantic decision | Obligation exercised |
|:--|:--|:--|
| `fixtures/01-promotion-accepted.json` | accept | adjacent E1 → E2 promotion with fresh evidence, exact binding, review, and a complete component set |
| `fixtures/02-promotion-rejected-skipped-level.json` | reject | an E3 result cannot skip an unsatisfied E2 gate |
| `fixtures/03-weakest-link-missing-component-rejected.json` | reject | one missing required component makes composed maturity `null` |
| `fixtures/04-demotion-weakest-link.json` | accept | a retracted map component demotes the current route to the weakest E3 component while preserving history |
| `fixtures/05-staleness-recheck.json` | accept | an expired time-bound E5 gate recomputes the fresh prefix to E4 |
| `fixtures/06-additive-correction-accepted.json` | accept | a correction appends a replacement with `supersedes` and retains the prior record |
| `fixtures/07-destructive-correction-rejected.json` | reject | an in-place rewrite is forbidden even when its rationale has evidence |

The schema checks representation, E0-E7 ordering, event-specific fields,
fail-closed draft markers, correction shape, and the accepted E1 → E2
fixture's explicit R2/R3 independence envelope. That envelope records the
independence method, reviewed event and subject binding, author and reviewer
identities, contribution status, conflict status, and receipt references.
Generic JSON Schema still cannot compare the two identity strings, prove that
the reviewed event and subject binding equal the enclosing fixture fields, or
resolve those receipts. Nor can it resolve gate identifiers against another
JSON document, compute minimum ranks, compare timestamps, prove that every
inherited gate is present, prove that a fixture's claim kind and component set
match its selected profile obligation, or enforce append-only storage. Those
are explicit requirements for the still-missing semantic validator and
governance workflow.

For an accepted E1 → E2 fixture, that validator must at minimum require
`reviewed_event_id == event.event_id`,
`subject_identity_binding == before.identity_binding`, distinct author and
reviewer identities consistent with `same_person=false`, an R2/R3 method that
matches the asserted level, an allowed contribution/conflict state, and
independence evidence references contained in the event's review receipts. It
must then resolve those receipts to the same identities, subject revision, and
review method. Missing or dangling receipts, equal identities, or a mismatched
event/binding must reject the transition.

## Rules captured by this draft

1. Maturity and claim status are orthogonal. A maturity transition never sets
   `supported`, `refuted`, `unknown`, or `not_applicable`.
2. Maturity is the highest level for which every required gate from E0 through
   that level is supported and fresh. A skipped or stale gate breaks the
   contiguous prefix.
3. A fixture's `claim_kind` selects its required component set from the route
   profile's `claim_obligations`. A composed route takes the minimum maturity
   of that set. An omitted required component produces no effective maturity
   (`null`); the profile-level `required_components` list is only the union of
   components used by its claim kinds.
4. Promotion is adjacent and append-only in this draft. It needs fresh gates,
   exact revision binding, the profile's independence requirement, and review.
5. Demotion is a recomputation, not deletion. Retraction, expiry, identity
   drift, or a newly refuted gate may lower the latest assessment while prior
   assessments remain immutable.
6. Corrections append a replacement assessment linked by `supersedes`.
   Destructive mutation of an earlier assessment is always rejected.
7. `not_applicable` is accepted only where a gate's explicit precondition
   fails; it is never a substitute for missing evidence.
8. Existence, explicit construction, arithmetic, map validity, subgroup
   preservation, end-to-end recovery, cost advantage, and scaling remain
   separate transfer claim kinds. Each declares formal-proof and computational
   certification requirements; neither silently substitutes for the other.

## Structural validation

From the repository root, a generic Draft 2020-12 validator can check the
registry and all fixtures. One reproducible option that does not use repository
code is AJV CLI 5:

```bash
npx --yes ajv-cli@5 validate --spec=draft2020 --strict=true \
  -s docs/ecdlp-audit/evidence/ws4-evidence.schema.json \
  -d docs/ecdlp-audit/evidence/evidence-profile-registry.draft.json

npx --yes ajv-cli@5 validate --spec=draft2020 --strict=true \
  -s docs/ecdlp-audit/evidence/ws4-evidence.schema.json \
  -d 'docs/ecdlp-audit/evidence/fixtures/*.json'
```

JSON syntax can be checked without schema tooling:

```bash
jq empty docs/ecdlp-audit/evidence/*.json \
  docs/ecdlp-audit/evidence/fixtures/*.json
```

## Promotion boundary

This directory deliberately has `authoritative=false` and
`promotion_allowed=false` in every record. WS3 remains the only upstream plan
dependency for completing the WS4 profile contract and transition validator.
WS5 and WS8 are downstream application prerequisites, recorded under
`open_relationships`, not additional WS4 dependencies: WS5 route freeze is
required before population-coverage claims use the profiles, while WS8
governance is required before the profiles drive an authoritative claim
transition or release. Neither item blocks WS4-G1 or WS4-G2, so no WS4↔WS5 or
WS4↔WS8 dependency cycle is introduced. A separate semantic validator must
still pass these golden fixtures plus malformed, alias, cross-reference,
identity-equality, reviewed-binding, conflict, and history-tampering cases.
Until then both WS4 acceptance gates remain open for their own stated gaps.
