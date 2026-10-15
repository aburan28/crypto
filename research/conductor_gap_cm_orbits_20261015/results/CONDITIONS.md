# Execution conditions

Construction-only orbit controls record thread CPU time, native operation
counts and exact in-memory point-record comparisons. The session is L0:
other construction and compilation work was present, so no isolated wall
speedup claim is made. Twenty shuffled repetitions cover all 100 gaps, three
scientific arms and one duplicate scalar A/A arm. Native FNV64 fingerprints
are supplemental; the backend also compares entire crater record vectors.

The square-root codomain reference run overlapped compilation. Its correctness
result is exact coefficient equality with all 100 historical direct codomains;
its wall times are descriptive L0 diagnostics. New upper-ladder construction
outputs without a supplied reference are independently certified afterwards.
The first constructor version emits `matches_direct_codomain: true` when
there is no reference; that value is not a validation receipt. Version 2
corrects it to null and adds `reference_supplied`. Preserve the original raw
output; use the independent PARI floor certificates for endpoint claims.

Small ecbench fixture control: six initial errors exposed loss of the supplied
public point in measured-child normalization. No solve was executed. The fix
preserves the public fixture in both measured-child and replay input builders.
The next session verified all six targets. These small controls do not replace
the requested 53-bit campaign or its replication requirements.

Primary launch: eight frozen ecbench shards, four co-scheduled single-target
workers, separate per-shard locks, no CPU affinity/isolation claim. Each curve
gets one independent rho solve on exactly its supplied public point; craters
also get both map and scalar implementations. This is 200 endpoint targets,
300 IC attempts and 200 rho attempts, with no shared tables between targets.
All arms are requested at L0. The command and executable/spec hashes are in
`campaign-pre-run.sha256`; per-shard logs retain failures. The co-scheduling
is a deliberate counts/correctness campaign, not an isolated timing claim.
