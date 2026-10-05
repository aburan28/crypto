# Disclosed SAT budget pilot result

The two fixed, previously inconclusive source instances in [the protocol](PROTOCOL.md)
both resolved at the preregistered 1,000,000-conflict cap. The probe binary,
CMS executable, original source instances and preparation matched the hashes
in [BUILD_PIN.json](BUILD_PIN.json) before the first role. That pin was pushed
to PR #1393 before execution. Each instance was run once in file mode and
once with a prestarted stdin reader, in the fixed order 000 then 015. Both
roles in each case agreed and met the prescribed status. There were no retries.

| Disclosed instance | Original 100,000-conflict result | New file / stdin result | CMS conflicts in each role | CMS thread CPU seconds, file / stdin |
| --- | --- | --- | ---: | ---: |
| Query 000, `[103925,114545]` | Inconclusive | Source UNSAT / source UNSAT | 500,204 | 29.24 / 31.19 |
| Query 015, `[103877,110181]` | Inconclusive | Source-valid SAT model / same model | 195,202 | 9.99 / 10.13 |

The [raw archive](result-v1/archive-manifest.json) contains byte-for-byte
copies of the original result JSON, source copy and both solver stdout/stderr
files for each case. The feasible query's file and stdin roles returned the
same checked model SHA-256
`dc55a90815a97255619fd2ca72bc9c8e451de753fda06ba42a8a0365cef8f05f`
and full-point witness indices `[32,13,54]`. The probe checked source
ANF/CNF/Magma bytes, model satisfaction and the actual group sum. A separate
[Sage replay](result-v1/sage-geometry-replay.json), run through the accepted
repository launcher with its [runtime receipt](result-v1/sage-runtime-info.json),
independently rebuilt the 63 geometric base points, found no three-sum for
query 000 and verified the full-point sum for query 015. Its
[source](result-v1/verify_geometry.py) is retained. The UNSAT result is a
solver result plus a separate exhaustive three-sum check on this toy base;
it is not a general proof-producing SAT certificate.

This is a **post-hoc selected two-instance diagnostic**, not a natural-query
sample. Selection used the original F5 geometry classifications, and both
original CMS attempts had already exhausted 100,000 conflicts. Consequently
2/2 resolution is not a yield estimate, rank estimate or basis for choosing a
one-target candidate. The probe explicitly records
`native_child_group_drain_audited=false` and
`source_bound_target_admitted=false`; it does not repair the consumed 512-query
CMS registration or supply a complete SAT factor-log table. The Mac was not
host-isolated, and the CMS thread CPU numbers are exploratory diagnostics,
not online wall times or speedups. The next scientific gate is a **newly
frozen natural panel** at a specified SAT budget, with its own source-bound
worker/audit, full rank and replayed logs before any SAT one-target claim.
