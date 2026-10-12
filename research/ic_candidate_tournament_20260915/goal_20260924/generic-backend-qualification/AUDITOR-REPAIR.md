# Natural-yield auditor repair, identified during the registered run

The 2026-09-29 campaign is executing the frozen PR #920 source. Its
`generic_backend_yield.py` calls `identity.sha256(report)` for each bounded
generic worker report. Candidate identity hashing deliberately rejects JSON
floats, while the worker reports finite floating-point timing diagnostics
(`elapsed_seconds` and collection timings). A retained control report
reproduces `InvalidEvidence: canonical records cannot contain floats or
non-JSON values`. The frozen auditor therefore cannot finish its natural-yield
step as written. This is an evaluator defect, not a measured solver outcome.
The repository auditor now uses `measurement.report_sha256` for future
campaigns; the running checkout and its archived evaluator remain unchanged.

The registered panel, worker, jobs, targets, timeout, resources, evaluator
snapshot and raw receipts must remain unchanged. Do **not** dispatch the same
registered campaign again. After the original Actions artifact is retained:

1. Keep the downloaded artifact and its SHA-256 immutable. Run its unmodified
   `tournament/evaluator/tournament.py verify --round tournament` and retain
   the full verification output. If the schedule is incomplete, report that
   failure; the repair does not invent missing trials.
2. On a separate working copy, run
   `recover_generic_backend_yield.py --bundle BUNDLE --out BUNDLE/natural-yield.json`.
   The wrapper requires the exact archived auditor SHA-256
   `627b5b84c2bbb9b9a9abb922a0ae2b6ed39ea70a9f212f087de7064e95637dc5`,
   the registered seed and schema, and the archive's own complete evaluator
   seal. It changes only the in-memory hash helper used for the report digest:
   the sealed `measurement.report_sha256` already used for admitted profile
   reports accepts finite JSON timing floats. The recovery cross-checks those
   digests against every verified generic trial receipt.
   All independent query, group, matrix, build and phase checks execute from
   the archived evaluator. The output records the repair script and contract
   hashes and labels itself post hoc.
3. If the full 750-pair schedule and qualification report exist, regenerate
   `summary.json` from the already-registered runner's read-only
   `panel_result`, with the repaired natural-yield file digest. Apply the
   current read-only `generic_backend_gate.py` to the working copy. Preserve
   the original failed auditor log and status beside the repaired output.

This repair is not preregistered evidence of a speedup or a new measurement.
The two-family gate still requires complete verified F4/F5 and SAT jobs, and
the incumbent/rho comparison still requires same-point, source-bound paired
costs. A timeout or interrupted campaign remains incomplete and censored.
