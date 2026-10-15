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

Endpoint-wave confound: primary shards 0–3 contain craters and shards 4–7
contain floors. The four-slot queue therefore ran the strata mostly in
different waves. Compilation, larger construction and replication overlapped
parts of those waves. Replication also queues endpoint-specific shards.
Whole-gap bootstrap intervals do not remove host/wave confounding. CPU
comparisons are descriptive conditional on these runs; deterministic query
counts, bases, relations and recoveries are the stronger controls. A causal
timing follow-up needs interleaved/randomized endpoint blocks and earned
isolation, including the D = −43 exploratory subset.

The final orbit-only repeat reserved CPU 0 and moved 51 threads. Its native
local monitor marked the 300.43-second run uncontended (0.34 other CPU
seconds). Docker/containerd threads remained eligible on the reserved CPU;
physical-host placement and frequency isolation are unverified. Do not turn
this local reservation receipt into an L2/full-pipeline wall claim.

Cold-process batches restart the native process once per 600-job pass,
using precomputed CM parameters. They overlapped validation/construction at
L0; individual model/arm jobs do not each get a fresh process or cache flush.
Whole-process CPU includes parsing, hashing, equality checks, serialization
and housekeeping; prior CM-map derivation is outside these microbenchmarks.
