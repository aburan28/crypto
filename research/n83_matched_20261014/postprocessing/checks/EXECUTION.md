# Supplemental export execution, 15 October 2026

The full independent replay and frozen main analysis passed before either
actual supplemental invocation. No experiment or solver was rerun here.

The first invocation rejected the serialized copy of three diagnostic fields
inside the checker receipt: `online_ms`, `phase_gap_ms`, `phase_sum_ms`.
Their normalized values differed by one binary64 step after the frozen analysis
read/wrote the receipt. Its stderr is retained as
`records/postprocess-20261015T0018Z/presentation_supplement.stderr`; it produced
no supplemental JSON. These three named clock aliases were added to the
existing one-step diagnostic comparison. Integers, types, provenance, identity
and non-clock values remain exact. All eight release tests passed again.

The second invocation passed, retaining logs with the `attempt2` suffix:

```sh
taskset -c 0 /workspace/n83_evaluation/build/release/n83-matched-supplement \
  /workspace/n83_matched/postprocess-20261015T0018Z/statistics_input.json \
  /workspace/n83_matched/postprocess-20261015T0018Z/statistics.json \
  /workspace/n83_matched/campaign \
  /workspace/n83_matched/postprocess-20261015T0018Z/presentation_supplement.json \
  research/n83_matched_20261014/records/campaign-20261014T1940Z
```

Output SHA256:
`4804814de5e0fbb7d464b09c26738844314a045cbb60cc33379ab7a9909ef0e8`.
Executable SHA256:
`51f9718b4ec219c65e3e766e08bfdf4c3ce5a9a7b31bc06e987a50f12a819925`.
The JSON provenance includes input, unchanged frozen-statistics, manifest,
library/main source, Cargo lock and executable hashes. Every target has
observation and original producer source pointers and hashes; every process
has observation and execution source pointers and hashes.

The exporter checked and described 104 IC targets, 114 rho targets and
32 actual processes. It calculated no primary estimator or bootstrap and
modified no frozen worker, controller, analysis or input artifact. Its
descriptive values use original producer/execution files as authority.
