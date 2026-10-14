# Validation record

All 100 profiles completed with exit code zero and passed the native analyzer.
It checked every raw-file checksum, re-derived Ir from the complete marker
parts, matched the frozen input, method, seed and workload identities, recovered
the exact archived scalar, independently checked `[k]G = Q`, and matched every
deterministic GAE, phase and inventory counter. All 100 receipts contain exactly
one solve part. Source sessions and frozen protocols remain unchanged.

There are 32 separate archived one-target workloads and 36 A/A pairs. Maximum
absolute relative instruction deviation is 0.002238%, within the frozen 2% gate.
The independent 99% bootstrap intervals for IC floor/crater include one in both
cells; the repeated instruction-lead criterion is false. Full precision and all
per-target rows are in `analysis.json` and `analysis.md`; no successful subset
was selected. This is a same-host implementation diagnostic, not another-host
profiling replication or a proof of equivalence, intrinsic hardness or security.

The native driver tests pass (two tests). All seven measurement-module checks
pass, including the process-exclusive strict-session test run separately. The
changed Rust sources pass rustfmt. Both interrupted raw archives and the
complete primary archive compare against all retained original file bytes with
native `tar --compare`. Their SHA-256 receipts identify the published containers.

The first final-validation materialization failed because the CI helper's
non-cone patterns were used in this cone-mode checkout. Its failure log is
retained. The corrected wrapper materialized 22 missing source/input files in
14 directories without changing fixtures or expectations.

The full release library suite completed: **4,088 passed, two failed and 117
ignored**, with no tests filtered out. It used optimization level 3, 256
code-generation units, all default features and three test threads in a separate
target directory. The full-run result is not green and its complete log is retained.

The two failures are the same JV-cover tests recorded on the earlier baseline,
but the exact observed assertions differ: replay checked only three valid
systems against its minimum, and Nagao trial 3 timed out before returning its
planted decomposition. Both failures and timing details remain in
`library-tests-final.log`. Both named tests pass in serial retries (24.33 s and
4.27 s); neither expectations nor time budgets were changed. These retries do
not replace or relabel the failed full run. The separately isolated
strict-session check passed on the final library test executable. The earlier
N83 sampler regression passes on this source snapshot; this study does not
claim to have corrected that sampler or discharged its separate problem record.

Publication checks: all 48 site tests and all 65 rho-status tests pass. The first
rho-status attempt could not bind its local HTTP server under the network-disabled
sandbox; the unchanged suite passes with network enabled. Both logs remain.
All generated leaderboard/browser/panel-index views pass their freshness checks;
the browser retains every historical row (98 sessions, 267 candidates).
Curve-name lint reports zero problems after restoring tracked prose and hash
pins. Every page of all five final report PDFs was visually inspected. Gallery
PDFs and SVGs match their sources byte-for-byte, and every PNG is a 150 dpi render.

Instruction counts exclude method resolution, input and known-target
construction, process launch, kernel work and the external point check. The
measured adapter includes its cold setup, scalar recovery, internal checks,
inventory and reporting work. Excluded work remains unmeasured rather than zero.
The prime-field F4/SAT integration, explicit transport, representation controls,
more independent classes and another-host profiling remain open.
