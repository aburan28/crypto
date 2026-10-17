# Repaired bounded survey

Actual detached local worker PID **79776**, launched **2026-10-17 06:58:03 UTC**.
At 06:59:13 UTC the native status shows audited vertices in all three classes;
one job is running. [Receipt](launch.json) records the exact snapshot and hashes.
Later status comes from the local ledger, not this historical receipt.

The queue contains all 24,823 vertices, with all 51 stratum representatives
first. Parameters stay at two public targets, one round, no warmups and 4,096
source trials; the 3-second child watchdog is unchanged. Invocation guards:
24 hours / 2 GiB spool. It pauses rather than deleting evidence. There is no
unattended Drive upload or provisioned S3 worker. Original full matrix and
A–E completion remain unresolved.

Inspect with the exact frozen driver:

```sh
/workspace/volcano_ic_repaired_survey_20261017/volcano-ic-driver status /workspace/volcano_ic_repaired_survey_20261017
```

The launch followed a full passing serial release suite, fresh audited control
and pilot, successful/exhausted replay and 58 exact completed baseline matches.
Two censored comparisons have no equality claim. Heavy local builds/tests and
report generation finished before launch; no runtime speedup is claimed.
