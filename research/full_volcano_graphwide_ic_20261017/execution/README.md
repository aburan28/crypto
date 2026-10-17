# Measured index-calculus execution

The control and pilot have run. The frozen [execution addendum](PROTOCOL.md)
defines a separate bounded survey over all 24,823 vertices of the three
additional complete graphs. It preserves the original full matrix and A--E.

The stopped dev-build checkpoint now covers **99 vertices and all 51 strata**:
**1,584 attempts, 1,065 verified, 422 exhausted, 97 timed out**. The native
[summary v2](survey-99-summary-v2/summary.json) contains per-stratum/arm/status
min/max/mean with explicit sample counts. The
[three-page inspected report](../reports/run-20261017T0337Z-all-strata/report.pdf)
plots every stratum. [Download the self-contained checkpoint](https://drive.google.com/file/d/1hLVWVad2SrgHaBn6U5iSyqnKbx5cAkeY/view)
or inspect its [delivery receipt](delivery-99.json). The supplemental
[graph explorer](../graph-viewer/README.md) shows every graph vertex/map while
marking the measured subset. This checkpoint covers every stratum,
not every vertex. The original dev campaign is paused, not running.

An optimized build passed the same control and three-root pilot. All **56
completed deterministic pairs** match the dev build exactly; **four censored
pairs** have no completed baseline report to compare. The
[native receipts](evidence/release-pilot-match.json), separate
[release QA summary](release-qa-summary-v2/summary.json) and frozen
[optimized continuation](optimized-continuation/optimization.json) retain that
distinction. No optimized speedup is claimed. The new full-queue campaign uses
a separate binary, ledger and output directory; planning alone is not execution.
It has now [actually launched](optimized-continuation/launch.json), PID 64909,
with a 24-hour/2-GiB guard and audited progress in all three classes. The
source stays pinned to clean commit `9bca26ae2` in a separate worktree.
The live local ledger is not the uploaded stopped dev snapshot.

Observed pilot: **48 attempts, 30 verified, 14 exhausted, 4 timed out**.
Small control: **12 attempts, 8 verified, 4 exhausted**. Both projected-base
oracles solve both control targets. All four sessions pass the standard audit,
replaying all 38 successful deterministic records. Supplemental native replay
also reproduces all 18 exhausted records. Timeouts have unknown final counters.
These are same-implementation replays, not independent implementation validation.

The pilot's source models are `icv1-fp57-t536870914-6cabb0c9`,
`icv1-fp57-tm536870941-e6d0c766`, and `icv1-fp59-t1073741954-4dd0b1ef`;
the control is `icv1-fp7-t2-4910f2d0`. Public targets, actual largest prime
subgroups, charged cofactor projection, checked base membership and checked
base logs are retained. No planted scalar reaches the solver. Completed A/A
rho attempts have identical counted results on their matched seeds.

The new base is `prime-projected-abscissa:size=N,scan_cap=M`, separately
identified from the unchanged raw-abscissa control. It cofactor-projects
rational points, discards infinity, checks [r]P=O, deduplicates negation pairs
and keeps both signs with coefficients 1 and r-1. Native tests recover every
nonzero control logarithm and reject an oversized base (2/2 pass).

## Run and resume

The tools have Cargo entry points:

```sh
cargo build --release --bin ecbench --bin volcano-ic-driver --bin volcano-ic-summarize
```

Alternatively link `driver.rs` and `summarize.rs` using [BUILD.sh](BUILD.sh)
and the exact matching rlib paths. This was the observed experiment build:
unoptimized dev code, debug symbols disabled, incremental off, two
build jobs. Counts and correctness are primary; wall time is L0/unqualified.
Exact build/test logs are in [evidence](evidence/).

```sh
volcano-ic-driver plan /fresh/campaign /certified/inputs /path/ecbench survey freeze
volcano-ic-driver run /campaign 24823 86400 1073741824
volcano-ic-driver status /campaign
volcano-ic-summarize /fresh/summary /campaign
```

`run` takes maximum additional vertices, elapsed seconds and campaign spool
bytes. Limits are checked between sessions/replays: a time limit can overshoot
by one bounded attempt. The spool guard reserves 64 MiB for another attempt;
it is not a quota. A 1 GiB spool pauses with all files retained. Resume with
the identical driver and unchanged binary/queue. A `STOP` file or SIGTERM
requests a graceful pause; remove only the stop marker before resuming.
Sealed sessions are never edited/deleted; retries get new numbered directories.
A failed audit or wrong answer halts execution. A second writer and a changed
driver are rejected before work begins; retained rejection logs verify this.

Execution uses Unix advisory locks/signals. Other platforms retain planning,
status, replay and summaries; non-Unix execution reports its unsupported status.
Only Linux execution was tested. Existing snapshots freeze the supplied driver
binary: use that binary for resume, because a later Cargo build has a different
binary hash even when its solver logic is unchanged.

The queue visits one representative of all 51 strata first, then every other
vertex, alternating graphs within each phase. Each vertex has its own sealed
session, audit and SQLite checkpoint. The queue preserves certified inputs
and exact finite-field applicability records. Their historical NOT_RUN fields
describe input-preparation state; `ledger.sqlite` records current execution.

The standard audit replays successful records only. Supplemental native replay
uses the harness measured-child function for deterministic exhausted records,
comparing outcome, exact phases, counters, field operations, base identity and
accounting flags. Its receipts stay outside the sealed session. Historical
diagnostic sessions can also be checked:

```sh
volcano-ic-driver replay-exhausted /sealed/session /fresh/receipt.json
```

The native summarizer emits linked observations and per-graph/stratum/arm
min/max/arithmetic mean with explicit n and separate outcomes. Timeout
counters stay missing. GAE totals are lower bounds wherever work is unpriced;
native `ecbench band` retains common 0/0.1/1 GAE sensitivity for all arms.
No speed ranking or admitted runtime claim follows from this diagnostic.
The [99-node native pricing band](survey-99-summary-v2/pricing-band.json)
contains all 792 curve/arm rows, using common weights 0/0.1/1 on every arm.
The [release QA band](release-qa-summary-v2/pricing-band.json) stays separate.
These are supplemental [uploads](delivery-pricing-bands.json) outside the frozen
archives; unavailable counts and incomplete pairs remain unrankable.

Yield accounting correction: `targets_tried` includes repeated source queries
skipped before the oracle. The 99-vertex checkpoint has **2,223,204 source
trials**, **2,075,145 actual decomposition calls** and **57,176 relations**.
Summary v2 reports yields per trial and per actual call separately; it checks
all 594 completed MITM call denominators against native lookup counts. Earlier
reports' "decompositions" label refers to source trials, not actual calls.
Those reports and sealed records remain unchanged. These pooled counts span
different subgroups/arms and are not a common success probability. Cross-binary
pooling and duplicate workload/arm/round/seed observations are refused.

For the projected size-64 MITM arm, native graph-aggregate rows show:

| Source model | Completed attempts / outcome | Actual-call yield min / mean / max |
| --- | --- | --- |
| `icv1-fp57-t536870914-6cabb0c9` | 60 verified | 0.101992 / 0.125942 / 0.156250 |
| `icv1-fp57-tm536870941-e6d0c766` | 50 exhausted | 0 / 0 / 0 |
| `icv1-fp59-t1073741954-4dd0b1ef` | 88 verified | 0.787709 / 0.860825 / 0.930380 |

These are arithmetic means over completed recorded attempts, separated by
outcome. The 45-bit prime subgroup has no relation in these 50 bounded attempts.
Different subgroup sizes and source models prevent attributing the differences
to conductor alone. See the native per-stratum rows for the actual populations;
these graph aggregates do not stand in for unmeasured vertices.

New sealed capsules and replay receipts stay loose where the verifier and
browser name their paths. Auxiliary logs are a deterministic
[archive](supplemental-run-logs/README.md), as required by AGENTS.md §13a.
The original first 12 capsules are not duplicated in the new 87-node increment.

## Open requirements

The real all-vertex survey has 4,096-trial caps, two public targets and one
round. It does not complete the original million-trial/eight-target/five-round
protocol. Full planning is supported; launch remains gated on that protocol's
exact registration and qualified resource pilot. Transported targets, composite
PH, independent implementation replay and finite-field IC for the actual
embedding degrees remain open. The supersingular degree-2 MOV teaching backend
does not implement those destinations. No generic fallback is labeled IC.
Original A/C/D/E full graphs remain storage-blocked.

[Measured data and inspected reports on Drive](https://drive.google.com/drive/folders/109DpcQpX7ATz6GeI8iSWAExGMWBoR-H_).
The measured archive is separate from the earlier inputs/NOT_RUN upload.
The earlier frozen survey snapshot has 12 of 24,823 vertices audited, with 192
attempts: 119 verified, 58 exhausted and 15 timed out. The full queue and
per-stratum completion ledger are included; queued rows are not solved rows.
The new 99-node archive preserves that earlier upload and includes the complete
three graphs, frozen dev binaries, queue, ledger, raw sessions, receipt-binding
inventory and corrected native statistics. Remote name, size and parent were
read back; a remote whole-file hash was not checked.

Conductor reservation was attempted for implementation/report scopes but the
tool is unavailable (127). Full-root historical test/catalogue obstructions
remain in the shared problem files; checks cited here are the new focused
tests, build, audits and replay receipts only. The benchmark build happened
before the measured sessions; short driver rejection checks occurred during
the counts-only survey, so its raw wall observations are not timing evidence.

Validation follow-up: the restricted full release suite has 4,225 passes, ten
failures and 117 existing ignored tests. All 37 networking tests pass with
local networking permitted; the exact Joux–Vitse test passes alone in 2.55 s.
Its parallel whole-suite timeout is a separate [open validation problem](../../../problems/release-suite-jv-nagao-wall-budget.md).
The full suite is not reported green. Frozen campaign controls/audits are
separate evidence; no solver assertion or sealed record was weakened.

The permitted full parallel release rerun finished with **4,234 passed, one
Joux–Vitse wall-budget failure and 117 existing ignored tests**; its full log
is `evidence/merged-release-tests-permitted.log`. All networking failures
were resolved by the appropriate permissions. The isolated exact Joux–Vitse
test passes; the whole-suite failure remains open and prevents an all-green
validation statement. No further heavy validation runs alongside the worker.
