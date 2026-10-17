#set page(paper: "us-letter", margin: (x: 0.6in, y: 0.53in))
#set text(size: 9.3pt)
#set par(leading: 0.4em)
#align(right)[TASK-20261016-7b4402 · 2026-10-16 · engineering observations]
= Native interruption evidence
Both prescribed existing native logging controls completed once, one per frozen feature arm. This descriptive report awaits fresh independent supplemental review. Source admission remains withheld under DEC-20261016-5cf723. Scientific trials and new leaderboard scores: *zero*.
#image("flow.svg", width: 100%)
#table(columns: (1.1fr, 1.2fr, 1fr, 0.8fr, 1.3fr), inset: 5pt, stroke: 0.4pt + rgb("#aebfcd"),
 [*Arm*], [*Passed / failed / filtered*], [*Begin / end*], [*OS wall*], [*OS max RSS*],
[Reference], [1 / 0 / 85], [2 / 0], [0.03 s], [2,850,816 B],
[Candidate], [1 / 0 / 85], [2 / 0], [0.03 s], [2,818,048 B],
)
*What ran.* Only the named interruption unit test from each existing reference/candidate test executable. Each parent spawned its ignored Session/stage/sleep helper; it observed the flushed begin, then killed and waited for the child. Per-parent bound: 20 seconds plus five seconds kill grace. No build, arithmetic, construction, map replay, benchmark, extra native test or retry occurred.

*Raw result.* Each trace contains exactly total/begin and interruption_fixture_wait/begin. Each control records observed-before-kill=true, unfinished-stage-marked-complete=false and cryptographic_work=false. Native child_exit_code=null remains unchanged; the native schema does not measure the Unix signal. All 32 original native files, parent and child logs, exits, OS receipts, timestamps, argv and explicit environment whitelist are retained.

*Cost scope.* One parent per arm. The OS wall numbers are rounded; CPU and RSS are shared-host operational receipts for pure logging, with full counters in raw/ARM/os-time.txt. The diagram describes the logging contract because no cryptographic quantity was measured. No performance ratio or scientific memory-admission claim follows.

*Provenance.* Frozen source #raw("015fa7bb504ef63953e0e04e464bb266f4f207bf"). Exact current unit hashes, sizes and feature/build/source links are in provenance-before.json and inputs/bindings.json. Unit hashes before/after and all four installed scientific executable hashes match. Current unit hashes bind these new invocations; historical invocation-time unit hashes were not retained. Old frozen015 interrupted raw remains absent.

#pagebreak()
= Preservation and next review
*Metadata correction, without native rerun.* The first metadata inspector incorrectly required 128 filtered tests from an aggregate multi-executable suite count. These actual unit executables selected exactly one named test and filtered 85 others. Original metadata.py, control-check.json (INCOMPLETE) and logs remain unchanged. A separate raw-log-only checker independently verified every retained file/hash and required guard, then emitted control-check-v2.json (PASS). No native control was repeated. A later JavaScript packaging call failed parsing before any write or execution; metadata-attempts.md retains its exact diagnostic and cause.

*Claim publication failure.* The initial write-once claim push was rejected after the remote branch advanced. Raw stdout/stderr and exit 2 are retained. The Coordinator prescribed a normal merge of 40 published commits into the existing epoch1 claim. Its 129 incoming files did not change this goal's frozen territory. Merge b2d66d0022405a77707017b2aa5bab11af15842a preserved both parents; frozen-goal diff, hashes and scoped branch hygiene passed. Publishing the existing claim succeeded. No second claim epoch or prepublication native run occurred.

*Custody.* The complete raw/control/source-binding/report payload and its member hash/size manifest are in native-controls.tar.gz. Separate custody/ receipts record the authorized public raw archive transfer and full readback to #raw("s3://crypto-autoresearcher/store"), logical suffix #raw("015fa7bb504e/interruption-retake-v1"). Large unit executable payloads are excluded as required. The transport executable is byte-bound; its build-source provenance remains unestablished. Metadata inspection uses Python only for source/hash/log/archive handling; native execution uses the published shell/gtimeout/time command.

*Prior review and root failure stay intact.* Signed review06cc13 supports exact semantics but left resource/build/custody/default-policy inconclusive because historical trace retention was absent. This supplement records new evidence and does not rewrite that review. The complete root gate at3ce remains 4161 passed / 1 failed / 117 ignored, exit101. A later filtered14-test pass does not replace it. Test thresholds, source, defaults and challenge rules are unchanged.

*Publication and authority.* This is producer TASK-7b4402 for GOAL-ISOU-f30a07. Requested Executor policy/effort and unexposed actual model/effort are retained in provenance-before.json. Exact raw/source/feature/executable identities are separate from interpretation. Full details: report.md, raw-manifest.json, control-check-v2.json and metadata-correction.stdout.log.

*Graphs checked.* All 54 standard-prime leaderboard rows and canonical curve/cost graphs remain unchanged: no map or mathematical finding is produced. P192/50341 remains a constructed, unscored candidate requiring admitted complete replay. Historical board incumbents are unchanged.

*Next action.* Coordinator snapshot0b8b61; fresh independent review38a4d9 of the logging invocation/retention binding; then prospective decisionaea57d/archive409603. Source admission and scientific execution remain withheld until those actual gates. The producer assigns neither.
