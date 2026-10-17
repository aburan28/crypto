# Fresh photographed-curve volcano campaign, 2026-10-17

This is a new experiment, explicitly authorized after the original run's files
became unavailable. It does not recover, overwrite or merge the previous random
fixtures, measurements, failures, source snapshots or manifests. PR #2065 will
retain that distinction and link this replacement campaign.

Exact source: p=72057980207847001, a=18456521102804593,
b=19958176276690439, #E=72057979853368804=4*r,
r=18014494963342201, from pinned shortlist commit
 ef18058588077ef805c55c5dc25140b6bde2046b.

Scope: close the rational degree-2 component; expected 46 curves, 90 directed
maps and five levels. "Complete" means this selected 2-volcano, not every prime
degree or every vertex in the full rational isogeny class. Prove kernel closure,
exact point counts, graph structure and point-order preservation.

Generate ten distinct subgroup generators (modulo sign) and ten distinct nonzero
uniform random k values using Linux getrandom and unbiased rejection sampling.
The same ten k values are shared across generators, yielding 100 problems per
curve, 4600 transported problem/curve cases. Keep every raw random word, rejection,
fixture, transport path and validation result. Expected k is checker input only.
Generate one independent secure rho seed per generator/key pair, shared across
curves; do not choose a seed for a favorable runtime.

Run the pinned native C rho and BSGS implementations from cryptanalysis commit
b7352eec52676b7b82a220ba296661e8ba878915. No globally-fastest claim without an
external matched comparison. Rho searches the full prime-order subgroup; no
small-k interval substitution. One initial native solve per problem/curve is
followed by any separately labeled instrumentation/replay phase; missing phases
remain open. Save all statuses, recovered values, counters, setup/query/check
clocks, whole-process resource readings and command/source/input hashes.

BSGS reuses one public-library baby-step table across the ten targets for each
generator. Table size = floor(sqrt(r-1))+1. Setup, every online solve, external
verification and destruction are separate. Preserve cold and derived/amortized
costs separately; never present a derived sum as an independently cold run.
Only one ~8-GiB table is resident at a time. Checkpoint every target and charge
any table reconstruction after interruption as a new attempt.

Rho budgets initially allow 2,000,000,000 upstream group operations per attempt;
a capped result remains unresolved and is retried as a retained new attempt,
never treated as a solve. Checkpoint each independent case; no data are dropped.

Implement actual prime-field Semaev/Gröbner IC relation diagnostics, with
independent point controls, relation yield and explicitly defined degree
metrics for all targets. Attempt actual relation/rank/recovery stages with
recorded caps. No oracle/imported factor-base logs or relabelled generic solve
may stand in for IC. A capped/zero-yield run leaves full IC recovery unfinished.
The request to solve all cases by IC remains open until independently verified
recovery, not merely successful diagnostic execution.

All local timings are descriptive L0 observations. Concurrent CPU work and
report generation are recorded; no isolation or wall-time speed advantage is
claimed. Full per-step arithmetic state streams are not collected; all available
aggregate counters and raw case records are retained. Unsupported hardware
counters and unimplemented metrics are explicitly missing, not zero.

Publish exact graph, fresh fixture identity and executable sources in the PR
before the long solver phase; append compact case checkpoints and summaries as
work advances. Preserve failed builds and interrupted attempts. For large raw
archives use the authorized immutable store when available; never invent a
storage locator or claim an upload occurred without successful readback.

## V2 sizing amendment before any fresh experimental solve

The pinned hash table rounds capacity up to a power of two with at least 2*m
slots. The automatic m=134218088 falls just above 2^27, doubling allocation.
Use explicit m=134217728 (2^27): 4 GiB per table rather than 8 GiB, with the
full interval still searched. The tiny increase in giant-step count is charged
as measured, not hidden. This is a configuration choice, not a measured speedup.
Two BSGS workers and two rho workers may run concurrently (two tables, ~8 GiB
combined), replacing the initial one-table scheduling plan. Preserve V1 and all
local timing/concurrency caveats. No existing measured result is reinterpreted.
