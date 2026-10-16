# First retained K=600/K=2,048 N83 mixed-orbit gate

The first frozen K=600 small-base / K=2,048 envelope natural-target trial
completed its exact import and one bounded solver call. It ended
`UNKNOWN_search_cap`: 5,000 conflicts, zero models, zero accepted partials and
zero completed rows. The trial coefficients and event counters passed the Python
supervisor's independent checks; the separate [receipt audit](audit.json)
passed and reproduced byte-identically. This is a censored one-trial search
observation. It is not evidence of relation impossibility, a factor-base
ranking, or a measured complete index-calculus runtime.

The Linux worker was built from clean source commit
`f3d27f68b68cda907e0668042ac7fa13ae6b579d`, with binary SHA-256
`26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606`.
The exact K=600 object was retrieved from
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e.jsonl.gz`
and matched compressed SHA-256
`3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf`.
The K=2,048 envelope at
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851.jsonl.gz`
matched compressed SHA-256
`aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
Both S3 HEAD responses were checked on 2026-10-12 and are recorded with
lengths and ETags in [objects.json](evidence-01/objects.json). No new factor
base was generated or uploaded in this round.

The `verify-alt` import reported 99,600 small points in 600 signed-Frobenius
orbits and 240,368 residual points in the 1,448-orbit envelope suffix. The
`run-alt` trial used the same public corpus and frozen first-trial coefficients
as the earlier K=1,182/K=2,048 run. It used three small and two residual
summands, ordered first-pair symmetry, a 5,000-conflict cap and a 32-model
cap. The event recorded 415,686 domain clauses, 100,472 peak variables,
664,020 peak clauses and 4,980 peak XORs. The worker's cgroup peak was
450,125,824 bytes. The null candidate and graph fields agree with every
zero-yield counter; no graph row or rank estimate was produced.

| Phase | Observed wall | Conservative charge | Receipt |
| --- | ---: | ---: | --- |
| Offline Linux build | 292.368 s | 600 s | `PASS_build_only`; ELF/hash and container cleanup verified |
| Exact archive import | 7.747 s | 360 s | `PASS`; no solver search |
| Bounded natural search | 350.173 s container, 348.334 s worker | 1,200 s | `UNKNOWN_search_cap`; one call, 5,000 conflicts |
| Cleanup/audit reserve | included in 363.915 s run supervisor | 240 s | both containers cleaned; frozen inputs and event replay verified |

All launched stages are charged their full declared caps. This round used the
entire **2,400-second** additional authorization; the conservative cumulative
local ledger is now **17,400 seconds** and its remaining allowance is zero.
The observed walls are recorded to explain the execution, not substituted
for the conservative budget or a complete-runtime measurement. The search
container had one CPU, 4 GiB memory, zero swap, no network and read-only
inputs. The host was shared and no isolation receipt was collected, so these
wall values are not an isolated throughput comparison.

The earlier [K=1,182 trial](../mixed-orbit-search-15000-20261012/RESULT.md)
also reached 5,000 conflicts with zero models and no relation. Its recorded
domain had 465,738 clauses and worker wall was 433.673 seconds. The K=600
domain is smaller and this single K=600 wall was shorter, but neither trial
produced useful relations or any downstream linear-algebra/individual-log
time. Those stage observations cannot identify the lower total-runtime base.

The retained [raw evidence](evidence-01/) contains authorization, build
receipt/logs and Linux worker, S3 object metadata, run configuration, import
and search logs, the one trial event, worker result, supervisor outer receipt,
and budget. The [auditor](audit.py) checks the frozen Git source blobs, input
and binary hashes, exact S3 object identities, container limits/exit and
cleanup, source-derived 415,686 domain clauses, first-trial coefficients,
event/result counters, null relation disposition, and full budget accounting.
It passed twice with identical output. A separate copy with its conflict
count changed from 5,000 to 4,999 was rejected. The worker did not leave a
container running. For replay on another host, download both S3 objects and
pass their local paths to `audit.py --small-object ... --envelope-object ...`;
the auditor checks their compressed SHA-256 before accepting the receipt.

The next retained experiment needs a new authorization. A known-answer
solver-access control is prepared in the study but has not run. For a
total-runtime ranking, both candidate bases still need natural useful-yield
and rank measurements, graph closure, sparse linear algebra, individual-log
recovery and a same-target rho comparison with independent replay and
exclusive phase accounting. This one capped sample should remain a ledger
row with a null selected base.
