# WS6 corpus controls

> **Draft, non-authoritative, and unfrozen.** Nothing in this directory is an
> executed literature search, an incident census, a verified mathematical
> claim, or evidence of historical completeness. `WS6-G1` and `WS6-G2` remain
> open.

This directory holds the control plane for Workstream 6 of the comprehensive
ECDLP audit. It separates a preregistered search method from the records that a
future execution may discover and from the claims a reviewer may eventually
verify.

## Artifacts

| Artifact | Role | Current authority |
|:--|:--|:--|
| [`PROTOCOL.md`](PROTOCOL.md) | Search, screening, acquisition, citation-pass, deduplication, and review protocol | draft; no query executed |
| [`source-claim-manifest.schema.json`](source-claim-manifest.schema.json) | Draft 2020-12 structural contract for sources, claims, WS5 route joins, links, runs, custody, and review | structural only; not an authority ledger |
| [`seed-manifest.json`](seed-manifest.json) | Five pointer-only discovery examples from three pinned repository blobs | fixture only; all substantive claims `unknown` at `E0`, and all WS5 route-family joins unresolved |
| [`fixtures/`](fixtures/) | One structurally valid empty draft and fail-closed negative cases | validator contract; not corpus data |

The seed deliberately has an `unfrozen` denominator, zero search runs, no
acquired source bytes, no reviewers, and no supported or refuted claim. Its
five records are not a sample from a defined population and must not be used to
calculate coverage. The top-level WS5 binding is only a mutable repository-path
reference. It has no content digest or review, and every source and claim has
`route_registry_join.status="unresolved"` with an empty `route_ids` array.
The binding declares `digest_scope="none"`; a future verified binding must use
SHA-256 over the WS5 registry's raw file bytes and declare
`digest_scope="raw_file_bytes"`.

## Requirement-to-evidence status

| Requirement | Status | Evidence or exact gap |
|:--|:--|:--|
| Preregister a public literature and incident-corpus method | **implemented but unverified** | `PROTOCOL.md`; independent protocol review and an executable frozen revision remain missing |
| Define a fail-closed record shape | **implemented, locally schema-replayed** | all draft-2 JSON parses, strict pinned Ajv Draft 2020-12 valid/invalid fixtures, and targeted fail-closed `jq` checks pass; a native semantic verifier and independent schema review remain required |
| Provide a non-authoritative seed | **verified for repository-pointer identity only** | five pointer records bind a commit, path, line locator, and Git blob; no primary source was acquired or reviewed |
| Join WS6 records to the WS5 route registry | **explicitly unresolved** | draft-2 records bind the mutable WS5 paths and carry fail-closed per-record join objects; no route ID is asserted, no registry digest is bound, and no reviewer has approved a mapping |
| Execute the systematic search | **not attempted** | `seed-manifest.json.search_runs` is empty |
| Establish corpus or historical completeness | **not established** | the denominator and comparator versions are unfrozen, searches and citation passes have not run, and unresolved gaps are explicit |

## Local structural and seed-integrity checks

These checks use native repository tools and `jq`; they do not invoke the
legacy Python research path.

```sh
find docs/ecdlp-audit/corpus -name '*.json' -print0 \
  | xargs -0 -n1 jq empty

jq -e '
  .authority.authoritative == false and
  .authority.promotion_allowed == false and
  .completeness.claimed == false and
  .completeness.denominator_status == "unfrozen" and
  (.search_runs | length) == 0 and
  ([.claims[].status] | all(. == "unknown")) and
  (.route_registry_binding.binding_status == "path_only_unverified") and
  (.route_registry_binding.content_sha256 == null) and
  (.route_registry_binding.digest_scope == "none") and
  (.route_registry_binding.authoritative_mapping == false) and
  ([.sources[].route_registry_join, .claims[].route_registry_join]
    | flatten
    | all(.status == "unresolved" and
          (.route_ids | length) == 0 and
          (.reviewer_ids | length) == 0 and
          .reviewed_at == null and
          .semantic_resolution_required == true and
          .authoritative == false)) and
  (.counts.sources == (.sources | length)) and
  (.counts.claims == (.claims | length)) and
  (.counts.source_claim_links == (.source_claim_links | length)) and
  (.counts.search_runs == (.search_runs | length))
' docs/ecdlp-audit/corpus/seed-manifest.json
```

Strict Draft 2020-12 replay, with version-pinned packages:

```sh
npm exec --yes \
  --package=ajv-cli@5.0.0 \
  --package=ajv-formats@3.0.1 -- \
  ajv validate --spec=draft2020 --strict=true --all-errors \
    -c ajv-formats \
    -s docs/ecdlp-audit/corpus/source-claim-manifest.schema.json \
    -d docs/ecdlp-audit/corpus/seed-manifest.json \
    -d docs/ecdlp-audit/corpus/fixtures/valid-empty-draft.json

npm exec --yes \
  --package=ajv-cli@5.0.0 \
  --package=ajv-formats@3.0.1 -- \
  ajv test --spec=draft2020 --strict=true \
    -c ajv-formats \
    -s docs/ecdlp-audit/corpus/source-claim-manifest.schema.json \
    -d 'docs/ecdlp-audit/corpus/fixtures/invalid-*.json' --invalid
```

The pinned Git pointers can be checked without network access:

```sh
rev=1e447956d3616092a223264f1c6c31eff8ac9979
git rev-parse "$rev:RESEARCH_ISOGENY_FIBERS_ECDLP.md"
git rev-parse "$rev:RESEARCH_KOBLITZ_INDEX_CALCULUS.md"
git rev-parse "$rev:research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md"
```

Expected blob OIDs, in the same order, are
`e838d0eaeffa9c3f2d30aa8bb88c6f54333fd055`,
`e57dca5ea84d95347e9b61c19888307f77ee765b`, and
`df1857387cf4ed2103fa05a364fd96d711add654`.

`jq empty` checks JSON syntax only. The pinned Ajv structural replay accepted
`fixtures/valid-empty-draft.json` and the seed and rejected every
`fixtures/invalid-*.json` file at its intended gate after the draft-2 change.
That passing replay is still not enough to promote evidence:
the native semantic verifier still must recompute counts, enforce identifier
uniqueness and referential integrity, verify repository and content hashes,
check temporal and supersession relations, resolve every nonempty route ID
against the exact bound WS5 registry revision, check lane consistency, and
enforce independent review. A path-only or unresolved route join must never be
treated as a route-family mapping.

The draft-2 schema also requires a `reviewed_non_authoritative` route join for
any `supported` or `refuted` claim, and requires a raw-file content-hash binding
before any record may carry a candidate or reviewed route ID. These structural
guards do not prove that the digest is correct, the route ID exists, or its lane
matches the record; those remain semantic-verifier obligations.

## Freeze and promotion boundary

The first execution must copy the draft protocol to a versioned, reviewed
revision and freeze the services, exact queries, dates, languages, budgets,
comparators, and denominator before issuing a production query. Pilot output
cannot enter the corpus unless repeated under that frozen revision.

No source or claim may become authoritative merely because it validates
against the JSON Schema. Promotion requires the protocol gates in
`PROTOCOL.md`, primary bytes or an immutable source revision, exact locators,
semantic verification, and the required second-person review. Missing access,
paywalls, timeouts, service caps, language deferrals, failed searches, and
unexamined citation edges remain counted gaps rather than negative evidence.
