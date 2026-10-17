# WS8 independent-governance policy (draft)

Status: **non-authoritative design draft**. This document and every record in
this directory have `authoritative=false`, `promotion_allowed=false`, and
`self_promotion_allowed=false`. They cannot change a scientific assessment,
adopt a policy, authorize execution, satisfy WS8-G1 or WS8-G2, or substitute
for protected repository controls and identified human reviewers.

## Purpose and scope

This policy separates creation of an ECDLP-audit claim from review of its
mathematics, implementation, replay evidence, authorization constraints, and
release decision. It governs proposed promotion, rejection, demotion,
staleness, additive correction, and external-replication registration. It does
not determine whether a mathematical or empirical claim is true; WS3 record
identity, WS4 evidence profiles, and WS7 accounting supply inputs that this
workflow must fail closed on.

## Non-negotiable invariants

1. An author cannot review, approve, or release their own revision.
2. One person cannot fill two required or invoked conditional review roles for
   one decision.
3. Every request and decision binds one proposed aggregate revision, one
   governance-schema revision, and one dependency-lock revision. The request's
   component list also binds the exact assessment, evidence, implementation,
   profile, accounting, and authorization revisions it relies on.
4. Missing, stale, conflicted, blocked, or rejected required input prevents a
   promotion. The release approver cannot override a specialist veto.
5. Claim status and E0-E7 maturity remain separate. Maturity is the highest
   contiguous satisfied profile prefix, and a composed route inherits its
   weakest required component.
6. Prior evidence, assessments, reviews, rejections, and failures are
   immutable. A change appends a new record and, when applicable, a
   `supersedes` edge.
7. CI validates machine-checkable predicates; it is not an independent
   mathematical review, replay, authorization decision, or release approval.
8. External replication is an evidence record. It does not automatically
   change status or maturity, fill unrelated gates, or grant the replicator a
   release role.
9. A bounded negative result stays bounded to its recorded population,
   hypotheses, method, date, and budget. It does not prove nonexistence or
   security.
10. Draft records in this directory always have no authoritative state effect,
    including the positive promotion drill.

## Roles

| Role | Owns | Independence floor | Cannot do |
|:--|:--|:--|:--|
| Author | request, revision bindings, disclosures, responses | originator | review or finalize the request |
| Mathematical reviewer | hypotheses, proofs, identities, exceptional cases, subgroup/recovery obligations, negative scope | distinct person | review a revision they authored; replace proof with tests |
| Implementation reviewer | executable correspondence, artifact identity, tests, malformed input | distinct person | infer mathematical validity from passing tests |
| Evidence/replay reviewer | fresh replay, provenance, receipts, failures, accounting | distinct person plus fresh replay | reuse the author's conclusion or omit failed attempts |
| Authorization/privacy reviewer | execution authority, sensitive material, privacy, publication constraints | distinct person | waive a required check without a recorded rationale |
| Release approver | gate aggregation, conflict check, append-only final decision | distinct person | override a veto or supply a missing specialist review |

An external replicator is an evidence source, not an internal approval role.
The machine-readable ownership proposal is in
[`role-ownership-matrix.draft.json`](role-ownership-matrix.draft.json). It is
not an installed CODEOWNERS or branch-protection configuration.

The authorization/privacy reviewer is conditional, matching the source plan's
"when needed" rule. The request records invoked conditional roles. This role
is invoked when the revision set includes an authorization record or the
transition changes execution authority, handling of sensitive or personal
material, publication constraints, or an external participant's privacy
scope. When invoked, the reviewer is independent and veto-capable. When no
adopted trigger applies, machine preflight records `authorization=not_required`;
absence of that role is not a waiver supplied by the author or release
approver. A future semantic validator must enforce these trigger rules.

## Conflict and independence procedure

Each actor records a stable review identity, relationship to the author,
contribution to the reviewed revision, conflicts, the reviewed aggregate,
governance-schema, and dependency-lock digests, decision, rationale, evidence
receipts, and timestamp. A person is ineligible for a required or invoked
conditional review role if they authored or materially produced the reviewed
revision, are the requesting author, have an undisclosed conflict, or cannot
perform the role's required fresh replay. A disclosed conflict leads to
recusal unless an adopted policy explicitly classifies it as non-disqualifying;
this draft defines no such exception.

The future semantic validator must compare identities across assignments and
decisions. JSON Schema alone cannot prove that two different strings identify
different people or that an independence attestation is true.

## Governed workflow

Every proposed action follows the same sequence:

1. The author submits a proposed append-only record with exact component
   revisions, the proposed aggregate/governance-schema/dependency review
   binding, invoked conditional roles, requested status/maturity,
   authorization scope, and a reason.
2. Machine preflight validates schemas, canonical identities, digests,
   cross-references, WS4 profile gates, WS7 accounting when applicable,
   authorization records, freshness, and immutable history.
3. A role-assignment check rejects self-review, duplicate required actors,
   conflicts, missing qualifications, and incomplete roles.
4. Each required reviewer and each invoked conditional reviewer records a
   decision against the same proposed aggregate, governance schema, and
   dependency lock. A mismatch, rejection, recusal without replacement, or
   blocked decision is a veto.
5. The release approver checks the full record and either appends a final
   decision or rejects it. The approver does not reinterpret specialist work.
6. An adopted authoritative ledger would append the resulting assessment,
   supersession, or replication entry. These draft fixtures stop before that
   step and report `state_effect=none_draft_only`.

| Workflow | Required result | History behavior | Critical fail-closed condition |
|:--|:--|:--|:--|
| Promotion | all profile, accounting, applicable authorization, freshness, role, and revision gates pass | append assessment | any missing/stale gate, conflict, veto, or revision mismatch |
| Rejection | record reasons and failed gates | no scientific state change; retain decision | a rejection cannot be relabeled as an inconclusive approval |
| Demotion | recompute status and contiguous maturity from current evidence | append assessment; preserve former supported record | refutation, retraction, identity drift, or stale required gate ignored |
| Additive correction | new record states reason, impact, replacement digest, and `supersedes` | append supersession | mutation or deletion of the prior record |
| Staleness | mark affected gate unavailable and recompute contiguous maturity | append assessment or stale marker | expired/event-invalidated evidence still counted as fresh |
| External replication | append scope, artifact digest, independence, result, failures, and applicable authorization/privacy review | append registry entry | replication automatically promotes or fills unrelated gates |

## Evidence staleness

- **Immutable-content evidence** has no time expiry merely because it is old,
  but digest mismatch, custody failure, or retraction invalidates reliance.
- **Revision-bound evidence** is valid only for the exact reviewed revisions;
  any relevant source, implementation, dependency, schema, or artifact change
  triggers replay.
- **Time-bound evidence** uses the evidence profile's explicit deadline. The
  draft policy supplies a 365-day upper review horizon but a stricter profile
  wins; elapsed time never silently converts to a pass.
- **Event-bound evidence** is rechecked after errata, counterexamples,
  incidents, toolchain defects, authorization changes, or other declared
  invalidators.

A stale required gate is unavailable for promotion and continued reliance.
The latest assessment must be recomputed from the surviving contiguous gates.
Restoration is a new revision-bound review, never a flag flip.

## Additive correction and demotion

A correction identifies the old and replacement records, both content
digests, reason, impact on status/maturity/cost/scope, new evidence receipts,
and a single explicit `supersedes` edge. The old record stays addressable.
Corrections that change a conclusion re-enter the applicable specialist,
replay, authorization, and release reviews.

Demotion is similarly additive. The triggering refutation, retraction,
staleness, identity drift, or authorization change is retained alongside the
prior assessment. Urgent quarantine may prevent downstream reliance before
review completes, but it is recorded as quarantine, not silently rewritten as
a reviewed demotion.

## External-replication registry

An entry records replicator and organization identifiers, relationship to the
authors, exact source and artifact digests, environment/toolchain, scope,
result (`replicated`, `not_replicated`, or `inconclusive`), failures, and
receipt. Negative and incomplete attempts remain visible. The evidence/replay
reviewer and, when an activation trigger applies, the authorization/privacy
reviewer check registry eligibility before the release approver appends an
entry. Any status or maturity change is a separate full transition that
rechecks all required gates.

## Release checklist

The machine-readable checklist lives in
[`governance-policy.draft.json`](governance-policy.draft.json). At minimum the
release approver confirms schema and semantic validation, canonical identities,
matching request/decision bindings for the proposed aggregate, governance
schema, and dependency lock, complete fresh profile gates, matched accounting,
applicable authorization, distinct eligible roles, absence of vetoes, failure
retention, and append-only history. A green pull request is only one input to
that check.

## Adoption and open obligations

This draft cannot be adopted by editing its flags. Adoption requires:

- reviewed WS3 canonical identities, immutable history, and authoritative
  ledger semantics;
- a native WS4 transition validator that checks cross-record identities,
  contiguous maturity, weakest-link composition, freshness, and all fixtures;
- a native WS7 accounting checker and route-appropriate profiles;
- installed path ownership, protected approvals, reviewer qualification and
  conflict procedures, and authorization controls;
- independent review of the policy, matrix, schema, dependency lock, and
  positive/negative drills at one exact aggregate revision; and
- a separately authorized decision creating a new non-draft policy version.

Until those obligations close, WS8-G1 and WS8-G2 remain open. The artifacts in
this directory are implemented draft contracts, not an operational governance
system or evidence that any ECDLP-related claim should be promoted.
