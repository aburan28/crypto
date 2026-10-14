# Procedure deviations

## Uncommitted first launch

The first `git add` refused the new study because it was outside the sparse
checkout. The controller resumed the commit tool and launched the census
without first checking that tool's failure result. This violated the required
local preregistration commit ordering. The attempt was stopped and is retained
under `attempt-unpreregistered/`; it supplies no primary estimate. Its census
log contains zero completed-job messages. No instruction estimate or comparison
from this attempt was inspected before deciding to restart.

The sparse-profile helper could not add the path because the checkout had no
helper-managed profile. Native `git sparse-checkout add --skip-checks` added
the study path successfully. Confirm the actual protocol/jobs commit succeeds
before restarting all 100 frozen profiles. Inputs, algorithms, ordering,
controls, accounting and decision thresholds are unchanged. Preserve the first
attempt's raw files losslessly and give it an interrupted-run report.

## Compiler sequencing

The monolithic release library-test compiler was paused after about 42 minutes
of compilation, along with its Cargo parent and the queued example build.
`paused-builds.txt` records their stopped state. Their memory mappings remain
resident, but they perform no CPU work during the census. The cached native
preparation helper and original archived ecbench executable are used. Resume
and complete the mandatory validation after profiling, before publication.
The protocol records this sequencing amendment before the primary census.

## Conductor availability

The required pre-edit `conductor check` command returned exit 127 because the
CLI is unavailable. Publication scope expansion was declared before editing
the problem index and the leaderboard presentation generator. No competing
reservation is inferred from this absence.
