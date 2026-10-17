# N83 sixteen-point finalist comparison, frozen before execution

The sealed fixture has sixteen previously unused public points, SHA-256
`407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062`,
and an answer sidecar, SHA-256
`fe771284a3b0e5c0301309f678751f0630f0fa23844dd0ce1dbc9fa47f243076`.
Its generator must pass `--verify-sealed` before a measured job. The answers
are correctness sidecars, never producer inputs. The exact IC source SHA-256
for all eight finalists is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.

| K | Root11 candidate suffix | Root22 candidate suffix | Actual usable points | Folded columns |
| ---: | --- | --- | ---: | ---: |
| 1104 | `h6f30a77abf9e` | `hce7cffedf9c0` | 183,264 | 1,104 |
| 1107 | `h8c736ae9c7b4` | `hcb73613ebc5a` | 183,762 | 1,107 |
| 1112 | `h49683406b22f` | `hbc5ec3072dc4` | 184,592 | 1,112 |
| 1120 | `hc922b499a234` | `h82b74c23b7bc` | 185,920 | 1,120 |

Each suffix completes the standard
`IC1N83Ckb1a1fb<B>PDP4rootRCguidedLAgaussTDdirectISO0` prefix with the
table's actual usable point count. The K1107 and K1112 root22 manifests are
configuration identities derived from replayed equal bases; the cloud
artifact and the new producer receipts must confirm those identities before
they enter a result row. A mismatch fails admission and keeps the row.

Run one fresh process with empty target-independent caches for every
candidate and point, 128 IC processes in total. Each process is one
single-target workload with its own IC setup; the sixteen public points are
independent paired repetitions, with no cross-target log or relation table.
Split points 0–7 and 8–15
between two Tier 1 cloud jobs. Within each job, use all eight candidates on
each point on the same host and resource envelope. Rotate the fixed table
order by local point index 0–7, so each candidate occupies each position
once. Each candidate/point is a distinct IC–rho pair. Execute a fixed-target
rho calibration before and after every pair. The calibration uses the exposed
single public point in `research/n83_cold_factor_base_20261016/targets.jsonl`,
its answer sidecar for independent replay, and seed `2026101601`. Run IC then
rho if `(local_point_index + fixed_candidate_index)` is even; reverse
that order otherwise. Each fresh same-point rho process uses seed
`2026101700 + 8 × point_index + fixed_candidate_index`, where the fixed index
is the candidate position in the table above. This makes 64 IC–rho pairs and
128 calibrations per host. Repeat K1107/root22 then K1104/root11 on the first
point at the end of each job as long-horizon diagnostics. They do not decide
host admission. Archive copies of the exact IC and rho binaries that were
executed. Cap each IC process at 480 seconds, each rho or calibration process at 120
seconds, and the producer stage at 190 minutes, leaving room for replay and
artifact upload inside a 240-minute cloud job. A global-cap stop is recorded
with the next event and leaves the half unadmitted. Preserve both repeats,
including failed and noisy outcomes. Admit
a job only when each before/after calibration online-time spread is at most
10% of its mean, each same-point IC–rho pair and calibration passes independent
scalar/resource replay, and
all Tier 1 topology, affinity, cgroup, throttling, steal-time, memory,
execution-order and source/target identity gates pass. A rejected job is
recorded. Any rerun repeats the complete frozen schedule under a new run ID
and retains every earlier failed pair.

The primary result is verified single-target **online** wall time from the
first target-dependent operation after reusable base, index and log setup
through scalar recovery and producer verification. The same-point rho clock
uses its target-dependent walk through scalar verification. Keep exclusive
online phases and independent scalar replay. The supplementary **cold**
interval is the full fresh IC process with empty caches, retaining setup,
rank, target, and residual phase costs. Every row records the actual base
size, folded columns, relation attempts and verified rank, memory, full
child wall, resource telemetry, source and binary digests, and raw artifact
binding. The two jobs form separate host strata; comparisons are paired
within each public point and pooled only after both strata pass admission.

For each metric, calculate all 28 unordered candidate pair differences over
the sixteen matched points. Use cold milliseconds directly and logarithms
of online milliseconds. Report paired means, medians, sample deviations and
descriptive two-sided 95% paired-t intervals, `mean ± 2.131449545559323 ×
sample_sd / 4`; exponentiate log endpoints to report online time ratios.
For selection, compute exact two-sided sign-flip randomization p-values by
counting all `2^16` signed sums with absolute value at least the observed
absolute sum, divided by `2^16`. Include equality in the count. A candidate
is selected for that metric only if it is faster in the observed direction
than each of the other seven and **every** corresponding two-sided p-value
is at most `0.05/28`. This fixed Bonferroni gate covers the choice among all
eight candidates. Otherwise its selected candidate is null. Report each
candidate's sixteen independently paired same-point `rho_online_ms / IC_online_ms` ratios, split by host, and
verification status separately. Never substitute the cold winner for the
online winner or a batch average for a single-target result.

`select.py` implements this rule over a
`n83.e2e-final-validated/v1` record emitted after independent replay of both
jobs. It requires all sixteen public-point rows, eight exact candidate IDs
per point, passing per-candidate rho, calibration and IC replay/resource statuses, both admitted Tier 1
host strata, and the frozen source and target hashes. Missing or failed gates
yield null selections. The materializer must bind that validated record to
both provider artifacts before the selection output is published.
`run.py` executes each eight-point stratum, `verify.py` checks its raw
commands and arithmetic, and `aggregate.py` checks the two provider ZIP
digests, every extracted member, local replay against cloud analysis, exact
source and binary hashes, host identities, and cross-stratum base equality
before it calls the selector. `materialize.py` writes one immutable IC1
manifest per verified candidate, 128 primary one-point IC rows, four repeat
control rows, and 128 same-point rho references. Failed, timed-out and
noisy executions remain in those rows with their recorded status and
unavailable costs as null.

The previously verified K3072/root11 online candidate has a different IC
source digest and its five-point Tier 1 result remains a separate reference.
A fresh same-source K3072 comparison, with its own exact manifest and host
capacity check, is required before it can enter this eight-candidate ranking
or a claim of the best online factor base across both source revisions.
