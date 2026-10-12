# Corrected tool-metadata preflight registration

The original registration `24c39da96c4a8c1df578560371e598467b38b240`
terminated before compilation and native execution. BSD `/usr/bin/ar` on this
Mac returns exit one and its usage text for `--version`; the original collector
incorrectly required zero for that optional metadata request. Neither query's
Boolean equations or root matrix was exported. This is an operational
preflight failure, not evidence about either mathematical hypothesis.

`preflight-v1/RESULT.json`, stdout/stderr, frozen original registration sources
and the 11,034,385-byte source/dependency archive preserve that failure. The
archive SHA-256 is
`4738f9b843462adeba67c1e48ba005f6a76553a41af9912437740844daf213fc`.
The first registration remains closed; the corrected runner rejects it.

`protocol-v2.json` permits one new execution after this metadata correction.
The corrected collector records the tool binary hash, exact argv, actual exit
code, stdout/stderr and whether the optional version probe succeeded. A failed
version request never becomes a fabricated version string or a successful
probe. Compilation and native execution still require zero exit and their
registered watchdogs, source identities and raw receipts.

The hypothesis, exact input file, historical Rust/dependency manifest, layouts,
criteria, degree, one-thread envelope and 1800/180-second limits are unchanged
from `PROTOCOL.md`. No solver search or DLP experiment took place in the first
registration, and this revision adds no new targets or mathematical policy.
Commit this correction, corrected runner and tool/closed-registration controls
before dispatching the new registration into a separate output directory:

```sh
env PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/goal_20260924/f5-boolean-system-control-20260930/run_control.py \
  --protocol research/ic_candidate_tournament_20260915/goal_20260924/f5-boolean-system-control-20260930/protocol-v2.json \
  --historical-source /absolute/immutable-765c3c5f19032bd852163805f257c56babef2040-checkout \
  --out /absolute/new-f5-boolean-control-v2-output
```

Preserve the second terminal result regardless of its verdict. Complete DLP,
production branch traversal, natural yield and calibrated speed remain unknown.
