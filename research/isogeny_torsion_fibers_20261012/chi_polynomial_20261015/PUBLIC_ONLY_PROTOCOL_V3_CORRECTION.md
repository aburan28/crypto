# V3 input custody correction

The frozen [V3 protocol](PUBLIC_ONLY_PROTOCOL_V3.md) says the independent
replay process alone may read the private answer fixture. That sentence
omitted the preprocessing step: the native `freeze-public-chi` command read
the private fixture once to remove its four `answer_scalar` fields. The
independent checker later read the private fixture for the field-by-field
input comparison and scalar replay. The `solve-roots-public` and
`exhaustive-roots-public` solver processes read only the redacted
`PUBLIC_FIXTURES_V3.jsonl` input. A deliberate call of
`solve-roots-public` with the private fixture exited 101 before timing or
relation collection, with stderr `public input contains answer scalar`.

The [freeze receipt](PUBLIC_FREEZE_V3.json) records both input hashes.
The [input/exactness report](reports/run-20261015T1850Z/REPORT.md) and
[relation solve report](reports/run-20261015T1851Z/REPORT.md) record the
commands, raw results, verification and limits. Historical protocol and
raw evidence are retained; this correction changes no fixture or result.
