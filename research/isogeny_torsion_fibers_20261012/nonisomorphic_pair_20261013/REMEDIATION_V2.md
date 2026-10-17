# V2 execution-wrapper correction

Date: 2026-10-13 local. Frozen before v2 execution. The first runner
execution produced 24 verified rows and a completed native summary, but
`/usr/bin/time -l` returned exit status 1 after the child completed because
the restricted host denied `sysctl kern.clockrate`. Preserve `raw-v1.tsv`
and `run-v1.log` as the wrapper-failure receipt. It gives no admitted wall
measurement and cannot serve as a clean process-exit acceptance result.

Run the identical release binary directly, with the same frozen candidate,
source, seeds, caps, factor-base rule, and target generation. Save stdout and
stderr as `raw-v2.tsv` and `run-v2.log`, record exit status, then run it again
as `replay-v2.tsv` and `replay-v2.log`. Require both direct exits zero and
byte-identical stdout/stderr; compare v1 stdout with v2 stdout. No solver or
protocol parameter changes are authorized in this remediation. Keep the
independent arithmetic verification and visual-report gates unchanged.
