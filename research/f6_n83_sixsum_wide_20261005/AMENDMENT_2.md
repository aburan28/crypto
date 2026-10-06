# Root-column cap and runner exit race

Added after the first frozen runs and before the follow-on root reduction.
The planted replay passed. All four ordinary T001 torsion-offset systems
constructed; offset 0 stopped at the registered conservative 300,000
distinct-column cap, without a reduction. Its 1,409,066 term occurrences
bound distinct columns by the same value. A 1,500,000-column cap thus
admits the complete own-degree column set while bounding a dense
415-row matrix to about 75 MiB. Keep the 120-second and monitored 7-GiB
RSS guards; run the extended reduction on offset 0 only. The earlier
`ordinary_0` result remains a separate inconclusive row.

The initial offset-3 runner returned status 1 after the probe had printed
its valid construction row: `ps` raced with process exit, and shell
`pipefail` stopped the monitor before writing its status. That raw JSONL
and empty stderr remain as `ordinary_3`; rerun it as `ordinary_3_retry`
with the monitor treating a missing process as an empty sample. The first
`planted_1` attempt likewise failed before the binary produced output
because the sandbox denied `ps`; `planted_2` was run with host-process
inspection allowed. No performance comparison uses these attempts.
