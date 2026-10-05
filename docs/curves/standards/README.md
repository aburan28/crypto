# Standards inventory and cover linkage

The importer accounts for **every row in its pinned inputs**, including failures.
It does not claim an exhaustive list of every curve ever adopted by every agency.
The initial inventory contains 265 source records: 258 imported records, deduplicated
by exact model, and seven unresolved records. The combined catalog has 304 models.
These counts include research proposals; they are not counts of agency approvals.

## Sources and coverage boundary

- `parameters.json`: the 248 records of the MIT-licensed
  [CRoCS standard-curves repository](https://github.com/J08nY/std-curves/tree/77fe6e3585ca2c2225b59d7df24b7c775437276f),
  pinned at that commit. The upstream project explicitly disclaims completeness.
  Categories include NIST, SEC, ANSI X9.62/X9.63, Brainpool, ANSSI, GOST, SM2,
  WTLS and Oakley, plus pairing curves and research proposals.
- `supplemental.json`: seven model records from
  [NIST SP 800-186 (February 2023), sections 3.2.1.6–3.2.3.3](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf),
  and ten polynomial-basis DSTU parameter sets transcribed from Bouncy Castle's
  `DSTU4145NamedCurves.java`. Each record carries its source and, where supplied, OID.
- `coverage.json`: every source identifier, input checksum, import status and
  unresolved reason. Four ANSI normal-basis records require a verified basis
  conversion; three extension/tower-field records require another field adapter.
  All seven remain visible with `exists: null`.
- `registry.json`: native-derived model identities, exact EC1 representations,
  and provenance. Missing full generator coordinates leave EC1 unassigned.

A dataset category is evidence of where a record was listed, **not an audit of
agency adoption**. The label `listed_by_cited_standard_or_standard_example`
includes historical examples; it does not mean current approval. Research
categories carry `not_asserted`. A worldwide historical census still needs a
versioned agency/standard/OID crosswalk, including withdrawn standards and
basis variants. Additional GOST representations and other agency inventories
must be sourced before that census can be called complete.

## Rebuild

```sh
cargo run --bin curve_standards --
python3 scripts/update_curve_standards.py
cargo run --bin curve_cover_check --
python3 scripts/build_lab_browser.py
```

Each command accepts `--check` (after `--` for Cargo). The Python steps only
merge or render native-generated metadata; curve arithmetic, field checks,
model identities and cover verification run in Rust. The older registry builder
also consumes the native standards registry when rebuilding other legacy inputs.
The merge preserves existing ICV1 records and EC1 UIDs. It refuses conflicting
model/order records, deduplicates exact models, and retains separate generators.

Source coefficients are interpreted as field elements, including signed integers,
and reduced in Rust before identity construction. Raw source records remain
unchanged. Supplied full generators are checked on the curve. Group orders and
primality remain declared inputs subject to the checks in [COVERS.md](../COVERS.md);
this importer does not prove subgroup orders. Missing or unsupported data never
becomes a nonexistence certificate.

## YAML graph

[`../cover-links.yaml`](../cover-links.yaml) is YAML 1.2 in JSON-compatible syntax,
so both a YAML 1.2 parser and an exact-integer JSON parser can read it. Its
[schema](../cover-links.schema.json) describes the envelope. Use full UIDs as
join keys; short aliases are display labels.

| Node | Identity | Meaning |
| --- | --- | --- |
| Elliptic model | ICV1 plus full model SHA-256 | Exact equation and field |
| Representation | EC1 plus `urn:ec-record:1:sha256:…` | Exact field encoding, equation, subgroup and generator |
| Hyperelliptic model | HC1 plus `urn:hc-model:1:sha256:…` | Field, equation coefficients and genus |
| Cover map | CV1 plus `urn:curve-cover-map:1:sha256:…` | Ordered source, target, map and degree |

A curve node's `cover_links` resolves into `covers` and `maps`. A map's target
model SHA-256 resolves back to a curve node and its EC1 representations.
The checker recomputes EC1 hashes and verifies field, coefficients and declared
order agree with the linked model. A matching field degree, name or j-invariant
alone never creates a link. Cover existence applies to the elliptic model;
Jacobian-to-EC1 subgroup transport remains `not_tested`.

HC1 hashes exactly its `record`: schema version, field, equation form, `h`, `f`,
and genus. CV1 hashes exactly its `record`: schema version, full source cover UID,
full target model hash and JSON, degree, `x`, `y_v`, `y_0`, and optional target
model map (explicit null otherwise). Hash input is sorted-key, compact UTF-8 JSON
with exact integers and no trailing newline. Metadata and display aliases are
excluded. These are new namespaces for this graph, not aliases for EC1.
