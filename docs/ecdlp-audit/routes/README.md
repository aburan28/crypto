# WS5 route-registry draft

**Status:** draft, non-authoritative, public/synthetic only. This directory
implements a bounded WS5 control-plane seed. It does not freeze the historical
attack inventory, establish applicability, promote evidence, authorize
production-target testing, or complete WS5-G1/WS5-G2.

## Requirement-to-artifact status

| Requested WS5 element | Artifact | Status |
| --- | --- | --- |
| Five non-overlapping lanes | [`route-registry.seed.json`](route-registry.seed.json) | implemented as 43 unassessed candidate families in lanes A--E; corpus and independent non-overlap review remain open |
| Route-specific evidence links | route records and lane profile links to the WS4 draft | structurally linked; no route is promoted and the WS4 draft is not authoritative |
| Typed hybrid composition | ten `requires`, `enables`, `conditions`, or `verifies` edges | implemented as an acyclic seed graph with weakest-link and anti-double-count rules; mathematical validity remains unknown |
| Exact coverage cells with no blanks | [`fixtures/valid/coverage-matrix.seed.json`](fixtures/valid/coverage-matrix.seed.json) | all ten cells in the frozen mock Cartesian product are populated; every scientific status is `unknown` and every cost is unmeasured |
| Fail-closed contract | [`ws5-routes.schema.json`](ws5-routes.schema.json) and `fixtures/invalid/` | strict draft flags, typed routes/edges, evidence-gated support, bounded denominator, and negative fixtures are specified; native authoritative-ledger integration remains open |
| Universal or historical completeness | none | explicitly not claimed; WS6 corpus reconciliation, denominator freeze, semantic replay, and independent approval are required |

## Files

- `ws5-routes.schema.json` is a Draft 2020-12 schema with
  `additionalProperties: false` on every record object. It accepts either a
  draft route registry or a bounded coverage matrix. Authority, promotion,
  historical-completeness, and universal-denominator claims are fixed false.
- `route-registry.seed.json` expands the five plan skeletons into candidate
  route families. Each route is `unknown`, has no evidence, permits only
  public/synthetic assessment, and links to exactly one WS4 lane profile.
- `fixtures/valid/coverage-matrix.seed.json` freezes five selected route IDs
  and two schema-only mock asset IDs. Its ten route-by-asset cells have no
  missing required fields, but this is field coverage only. It supplies no
  executable curve, production implementation, secret, measurement, or
  scientific conclusion.
- `fixtures/invalid/` contains documents that must be rejected:

  | Fixture | Required rejection |
  | --- | --- |
  | `authoritative-draft.json` | a `v0-draft` record cannot set `authoritative: true` |
  | `supported-without-evidence.json` | `supported` cannot have an empty evidence list |
  | `universal-denominator.json` | the bounded mock population cannot set `universal_denominator: true` |
  | `untyped-hybrid-edge.json` | `implies_success` is not an allowed hybrid-edge type |

## Local vocabulary projections

The WS5 draft uses two local sentinel values that are not members of shared
plan vocabularies. `vocabulary_mappings` records one-way projections; it does
not declare vocabulary equivalence:

- `public_synthetic_only` projects passive inspection of public evidence to
  plan tier `T0` and restricts fixtures to public or locally generated
  synthetic known-answer material. An executable local-fixture run still
  requires its own explicit `T1` receipt under the plan. The local sentinel is
  not itself an execution authorization and never permits production targets.
- `unmeasured` is a WS5 pre-measurement sentinel. On export to the plan's
  `cost_status` vocabulary it projects to `unknown`, with the accounting
  receipt and every numeric cost field null. It is not a sixth shared-plan
  cost status and makes no measured, modeled, or extrapolated cost claim.

## Lane assignment and anti-double-counting

Assignment follows the master plan's first-match decision order:

1. Lane A computes a scalar in the source problem representation.
2. Lane B transports the subgroup to another algebraic object and requires a
   verified recovery path.
3. Lane C records a parameter or provenance property; its exploitation is a
   separately linked lane-A or lane-B mechanism.
4. Lane D uses protocol, implementation, observation, fault, or operational
   behavior and leaves the source DLP unsolved.
5. Lane E changes trust in identity, evidence, artifacts, measurements, or
   delivery. Direct secret exposure remains a linked lane-D record.

A composed route is counted once, not once per component. It inherits its
weakest required component. An edge records a dependency or evidentiary
relationship; it does not establish that either endpoint exists, applies, or
improves a solver. In particular, a curve cover must link onward to verified
Jacobian transfer, a genuine destination solver, and recovery before it can
support an end-to-end claim.

## Frozen mock matrix boundary

The valid matrix's denominator is only:

```
5 selected seed routes × 2 schema-only mock assets = 10 cells
```

Its population digest is SHA-256 over the exact canonical bytes of
`population.population_preimage`: compact JSON, object keys sorted
lexicographically, strings emitted as raw UTF-8 rather than ASCII escapes, and
no trailing line feed. The reproducible command is `jq -cS -j`; plain
`jq -cS` is not equivalent because it appends a line feed. Zero blank cells
means that all fields required by this draft contract occur in those ten
cells. It is not a completeness percentage, not a claim that the 43 registry
candidates are complete, and not evidence that any cell is resolved. All ten
cells remain `unknown`; all costs remain `unmeasured`.

## Native and shell validation

The repository's native authoritative-ledger validator does not yet consume
this draft. Until that integration lands, pinned AJV validates the full Draft
2020-12 schema and the following non-Python checks validate JSON syntax, the
bounded population digest, internal IDs, WS4 profile links, Cartesian coverage,
and graph acyclicity. They do not promote a record or replace future native
semantic replay.

```sh
find docs/ecdlp-audit/routes -name '*.json' -print0 \
  | xargs -0 -n1 jq empty

npm_config_cache=/tmp/ecdlp-ws5-ajv-cache \
  npx --yes --package=ajv-cli@5.0.0 \
  --package=ajv-formats@3.0.1 -- \
  ajv validate --spec=draft2020 --strict=true --all-errors \
  -c ajv-formats \
  -s docs/ecdlp-audit/routes/ws5-routes.schema.json \
  -d docs/ecdlp-audit/routes/route-registry.seed.json \
  -d docs/ecdlp-audit/routes/fixtures/valid/coverage-matrix.seed.json

jq -e '
  [.lanes[].route_families[].id] as $ids |
  ($ids | length) == ($ids | unique | length) and
  ([.hybrid_edges[].id] | length) == ([.hybrid_edges[].id] | unique | length) and
  .vocabulary_mappings.authorization_scope.local_value == "public_synthetic_only" and
  .vocabulary_mappings.authorization_scope.passive_public_plan_tier == "T0" and
  .vocabulary_mappings.authorization_scope.shared_vocabulary == false and
  .vocabulary_mappings.cost_state.local_value == "unmeasured" and
  .vocabulary_mappings.cost_state.plan_projection == "unknown" and
  .vocabulary_mappings.cost_state.shared_vocabulary == false and
  all(.lanes[]; .id as $lane |
    .evidence_profile.lane == $lane and
    all(.route_families[]; .lane == $lane)) and
  all(.hybrid_edges[]; .from_route_id as $from | .to_route_id as $to |
    ($ids | index($from)) != null and
    ($ids | index($to)) != null and
    $from != $to)
' docs/ecdlp-audit/routes/route-registry.seed.json

jq -e --slurpfile profiles \
  docs/ecdlp-audit/evidence/evidence-profile-registry.draft.json '
  [$profiles[0].route_profiles[].id] as $profile_ids |
  all(.lanes[]; .evidence_profile.id as $profile |
    ($profile_ids | index($profile)) != null and
    all(.route_families[]; .evidence_profile_id == $profile))
' docs/ecdlp-audit/routes/route-registry.seed.json

for mapping_document in \
  docs/ecdlp-audit/routes/route-registry.seed.json \
  docs/ecdlp-audit/routes/fixtures/valid/coverage-matrix.seed.json
do
  jq -e --slurpfile plan docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json '
    .vocabulary_mappings as $mapping |
    ($plan[0].vocabularies.authorization_tiers |
      index($mapping.authorization_scope.passive_public_plan_tier)) != null and
    ($plan[0].vocabularies.authorization_tiers | index("T1")) != null and
    ($plan[0].vocabularies.cost_status |
      index($mapping.cost_state.plan_projection)) != null and
    $mapping.authorization_scope.shared_vocabulary == false and
    $mapping.cost_state.shared_vocabulary == false
  ' "${mapping_document}"
done

jq -r '.hybrid_edges[] | "\(.from_route_id) \(.to_route_id)"' \
  docs/ecdlp-audit/routes/route-registry.seed.json | tsort >/dev/null

jq -e --slurpfile registry \
  docs/ecdlp-audit/routes/route-registry.seed.json '
  ($registry[0] | [.lanes[].route_families[].id]) as $registered_routes |
  .population.population_preimage.route_ids as $routes |
  .population.population_preimage.asset_ids as $assets |
  ([$routes[] as $route | $assets[] as $asset | "\($route)|\($asset)"] | sort) as $expected |
  ([.cells[] | "\(.route_id)|\(.asset_id)"] | sort) as $actual |
  all($routes[]; . as $route | ($registered_routes | index($route)) != null) and
  ([.assets[].id] | sort) == ($assets | sort) and
  $actual == $expected and
  .coverage.expected_cells == (($routes | length) * ($assets | length)) and
  .coverage.actual_cells == (.cells | length) and
  .coverage.blank_cells == 0 and
  .population.population_digest_contract == {
    "hash": "SHA-256",
    "serialization": "compact_json",
    "object_key_order": "lexicographic",
    "string_encoding": "raw_UTF-8",
    "ascii_escape": false,
    "trailing_lf": false
  } and
  .vocabulary_mappings.authorization_scope.shared_vocabulary == false and
  .vocabulary_mappings.cost_state.plan_projection == "unknown" and
  .vocabulary_mappings.cost_state.lifecycle_state == "pre_measurement" and
  .vocabulary_mappings.cost_state.numeric_fields == "all_null" and
  all(.cells[];
    (.prerequisites | length) > 0 and
    (.controls | length) > 0 and
    (.open_obligations | length) > 0 and
    .costs.status == "unmeasured")
' docs/ecdlp-audit/routes/fixtures/valid/coverage-matrix.seed.json

population_digest=$(
  jq -cS -j '.population.population_preimage' \
    docs/ecdlp-audit/routes/fixtures/valid/coverage-matrix.seed.json \
    | sha256sum | cut -d' ' -f1
)
recorded_digest=$(
  jq -r '.population.population_digest_sha256' \
    docs/ecdlp-audit/routes/fixtures/valid/coverage-matrix.seed.json
)
test "sha256:${population_digest}" = "${recorded_digest}"

! jq -e '.authoritative == false' \
  docs/ecdlp-audit/routes/fixtures/invalid/authoritative-draft.json >/dev/null
! jq -e 'all(.cells[];
  ((.claim_status == "supported" or .claim_status == "refuted") | not) or
  (.evidence_refs | length) > 0)' \
  docs/ecdlp-audit/routes/fixtures/invalid/supported-without-evidence.json >/dev/null
! jq -e '.population.universal_denominator == false' \
  docs/ecdlp-audit/routes/fixtures/invalid/universal-denominator.json >/dev/null
! jq -e 'all(.hybrid_edges[];
  .edge_type == "requires" or .edge_type == "enables" or
  .edge_type == "conditions" or .edge_type == "verifies" or
  .edge_type == "supplies")' \
  docs/ecdlp-audit/routes/fixtures/invalid/untyped-hybrid-edge.json >/dev/null
```

Strict AJV validation accepts the registry and valid matrix and rejects all
four negative fixtures at their named constraints. A future native validator
must additionally reject duplicate IDs, dangling references, cycles,
population-digest mismatches, missing Cartesian cells, and promotion without
the route-specific WS4 evidence profile in one semantic replay.

## Dependencies and next obligations

- **WS3:** issue canonical immutable IDs for the route registry, matrix,
  assessment units, and evidence revisions; add native schema and semantic
  replay.
- **WS4:** independently review the draft lane profiles and require contiguous
  E0--E7 maturity plus weakest-link demotion.
- **WS6:** complete the preregistered public corpus search, record unavailable
  sources, reconcile synonyms and splits, and propose an independently
  reviewed route-population freeze.
- **WS7:** attach complete cold, online, amortized, failed-work, memory, and
  uncertainty receipts before any cost comparison.
- **WS8:** approve role-separated promotion, correction, and denominator
  changes. A registry freeze and a scientific finding require separate
  decisions.

Only after those obligations are met can WS5-G1 or WS5-G2 be reconsidered.
A failed or bounded search remains a scoped negative result, never universal
nonexistence or safety.
