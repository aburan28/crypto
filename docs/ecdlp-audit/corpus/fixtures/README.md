# WS6 schema fixtures

These files define the first fail-closed cases for a future Draft 2020-12
validator and native semantic verifier. They are test data, not source or claim
records.

| Fixture | Expected structural result | Reason |
|:--|:--|:--|
| `valid-empty-draft.json` | accept | empty, explicitly non-authoritative draft with an unfrozen denominator and path-only WS5 binding |
| `invalid-authoritative.json` | reject | draft authority is set to `true` despite the schema's `false` constant |
| `invalid-completeness-claimed.json` | reject | completeness is claimed before a denominator exists |
| `invalid-frozen-without-time.json` | reject | a frozen protocol has no freeze timestamp |
| `invalid-supported-without-review.json` | reject | a claim with a synthetic reviewed route join is marked supported without second-person claim verification |

Every fixture is valid JSON and uses the draft-2 route-registry binding shape.
Rejection is a schema-validation result, not a parse failure. Additional
semantic negatives are still owed for duplicate IDs, dangling links, route IDs
that do not resolve against the bound registry, lane mismatches, mismatched
counts and hashes, review self-approval, cycles, and search-run accounting.
