#!/usr/bin/env python3
"""Render descriptive reporting from existing metadata only; no experiment."""
import json,pathlib
s=pathlib.Path(__file__).resolve().parent
c=json.loads((s/"control-check-v2.json").read_text())
p=json.loads((s/"provenance-before.json").read_text())
assert c["status"]=="PASS"
rows=[]
for arm in ("reference","candidate"):
    a=c["arms"][arm];t=a["os_receipt"]
    rows.append(f"| {arm} | 1 / 0 / {a['parent_test']['filtered']} | 2 / 0 | {t['real_seconds']:.2f} s | {t['max_rss_bytes']:,} bytes |")
data={"schema":"isogeny-trace-retake-visual-data/v1","source":"control-check-v2.json","task_id":"TASK-20261016-7b4402","arms":{a:{"parent_passed":x["parent_test"]["passed"],"parent_filtered":x["parent_test"]["filtered"],"trace_begin_count":len(x["trace_records"]),"trace_end_count":0} for a,x in c["arms"].items()},"scientific_runs":0,"new_scores":0,"scope":"Pure logging engineering control; diagram conveys the observed execution contract, not isogeny or performance evidence."}
(s/"visual-data.json").write_text(json.dumps(data,indent=2)+"\n")
svg='''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="400" viewBox="0 0 1100 400">
<defs><marker id="arrow" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#344e64"/></marker></defs>
<rect width="1100" height="400" fill="#f7fafc"/>
<g font-family="Arial, sans-serif" fill="#142d42"><text x="30" y="34" font-size="21" font-weight="bold">New, bounded logging evidence</text>
<text x="30" y="61" font-size="16">TASK-20261016-7b4402 · source 015fa7bb504e · two native parents, once per arm</text>
<rect x="30" y="88" width="270" height="103" rx="10" fill="#e0edf7" stroke="#587f9b"/>
<text x="49" y="117" font-size="19" font-weight="bold">Reference unit executable</text><text x="49" y="145" font-size="16">1 named test passed; 85 filtered</text><text x="49" y="170" font-size="16">No cryptographic work</text>
<rect x="30" y="214" width="270" height="103" rx="10" fill="#e0edf7" stroke="#587f9b"/>
<text x="49" y="243" font-size="19" font-weight="bold">Candidate unit executable</text><text x="49" y="271" font-size="16">1 named test passed; 85 filtered</text><text x="49" y="296" font-size="16">No cryptographic work</text>
<path d="M300 140 H360 M300 264 H360" stroke="#344e64" stroke-width="2" marker-end="url(#arrow)"/>
<rect x="375" y="111" width="310" height="180" rx="10" fill="#fff" stroke="#587f9b"/>
<text x="396" y="143" font-size="20" font-weight="bold">Child helper, each arm</text><text x="396" y="174" font-size="17">total / begin</text><text x="396" y="203" font-size="17">interruption_fixture_wait / begin</text><text x="396" y="236" font-size="17">Parent observes, kills and waits</text><text x="396" y="265" font-size="17">No end or completion record</text>
<path d="M685 201 H744" stroke="#344e64" stroke-width="2" marker-end="url(#arrow)"/>
<rect x="758" y="111" width="310" height="180" rx="10" fill="#fff4db" stroke="#a77a20"/>
<text x="778" y="144" font-size="20" font-weight="bold">Retain and bind</text><text x="778" y="175" font-size="17">32 original native files</text><text x="778" y="204" font-size="17">Source / argv / features / hashes</text><text x="778" y="237" font-size="17">Independent review pending</text><text x="778" y="266" font-size="17">Source admission remains held</text>
<text x="30" y="354" font-size="17" font-weight="bold">Historical distinction</text><text x="30" y="381" font-size="16">Old frozen015 raw trace remains absent. This supplement records new invocations; all earlier findings remain intact.</text>
</g></svg>
'''
(s/"flow.svg").write_text(svg)
md=f"""# Native interruption evidence supplement — 2026-10-16

This is a descriptive engineering report for GOAL-ISOU-f30a07 / TASK-20261016-7b4402. The two prescribed native logging controls completed once per arm. These are observations for independent supplemental review TASK-20261016-38a4d9. Source admission remains withheld under DEC-20261016-5cf723; scientific runs and new leaderboard scores are both zero.

![Execution and retention](flow.svg)

The complete frozen protocol is in inputs/PROTOCOL.md; inputs/native-command.sh is the exact shell launcher. The launcher selected only replay_trace::tests::flushed_begin_survives_bounded_child_interruption from each existing native unit executable. Each parent spawned its ignored Session/stage/sleep helper, waited for the begin marker, then killed and waited for the child. Per-parent bounds were 20 seconds plus five seconds hard-kill grace. No build, cryptographic construction, map replay, benchmark, extra native test, or adaptive retry occurred.

| Arm | Parent passed / failed / filtered | Trace begin / end | OS wall, rounded | OS maximum RSS |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

These are operational OS receipts from a shared Apple host for a logging fixture, including parent/descendant accounting as emitted by /usr/bin/time -l. User/system values and all OS counters are retained in raw/ARM/os-time.txt. They supply no isolated speed ratio or scientific memory-admission result. Each arm's sample size is one parent invocation. The flow diagram describes the logging contract because no cryptographic quantity was measured.

Each retained control.json records flushed_begin_observed_before_kill=true, unfinished_stage_marked_complete=false, cryptographic_work=false and child_exit_code=null. Each interrupted.jsonl contains exactly total/begin and interruption_fixture_wait/begin. There is no end or terminal record. The native schema does not measure the Unix termination signal; null is preserved without inventing a signal value.

## Exact provenance and retained inputs

Frozen source is {p['source_commit']}; source trace bytes match that Git object. Current observed unit hashes, sizes and original build paths are in provenance-before.json and inputs/bindings.json. Both unit executable hashes and all four installed scientific executable hashes match before and after the controls. Source/features link through retained original build/test logs and source markers; unit hashes were observed for these new invocations and were not retained at the historical full-suite invocation. The supplement does not recover or relabel the missing old015 trace. Older e82 trace evidence remains older evidence.

The 32 original native files are copied byte for byte into raw/ and bound by raw-manifest.json. Parent and child stdout/stderr, exact argv, explicit environment whitelist, exit statuses, OS receipts and timestamps are retained. Metadata-only source, hash, archive and log inspection uses Python; native execution uses only the published shell launcher, gtimeout and time. Actual served model/effort are unexposed and recorded unknown; requested Executor policy/effort and no substitution are recorded.

## Preserved administrative and metadata failures

The first write-once claim push was rejected because the remote research branch advanced. Its stdout/stderr and exit2 are retained. The Coordinator prescribed a normal merge of the 40 remote commits into the existing local epoch1 claim; the 129 incoming files did not change this goal. Merge b2d66d0022405a77707017b2aa5bab11af15842a preserved both parents. Frozen-goal diff and binding checks passed; scoped branch hygiene passed171candidatepaths; publishing the existing claim then succeeded. No second claim epoch was allocated and no native test ran before publication.

The first metadata inspector incorrectly required128filtered, taken from an aggregate multi-executable full-suite count. Both actual unit logs have exactly1namedtest passed and85filtered. Original metadata.py, control-check.json (INCOMPLETE), and their logs remain unchanged. correct-log-count.py independently reparsed the original raw files and hashes, retained the observed count, and emitted control-check-v2.json (PASS). The native controls were not rerun. This correction changes derived metadata interpretation only. A later JavaScript packaging call failed parsing before any write or execution; metadata-attempts.md retains the exact diagnostic and cause.

## Archive, review and unchanged frontiers

The complete raw/control/source-binding/report payload and its hash/size manifest are archived as native-controls.tar.gz. Separate custody/ receipts bind the public raw archive to the authorized logical name evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/interruption-retake-v1 at s3://crypto-autoresearcher/store and retain full readback verification. Large unit executables are not uploaded in this supplement. Transport-client bytes are bound; its build-source provenance remains unestablished.

The prior signed review06cc13 remains unchanged: semantics holds; resource/build/custody/default-policy was inconclusive because historical raw retention was missing. Current full root gate at3ce remains4161passed/1failed/117ignored, exit101; the filtered14-test pass does not replace it. No threshold, source, default or challenge rule changed.

The54standard-prime leaderboard rows and all canonical mathematical/curve/cost graphs remain unchanged because this fixture produces no isogeny or graphable mathematical finding. P192/50341 remains an unscored constructed candidate awaiting an admitted complete replay. Next: Coordinator snapshot0b8b61, fresh independent review38a4d9, then prospective decisionaea57d/archive409603. This producer assigns no admission or scientific status.
"""
(s/"report.md").write_text(md)
typ_rows=[]
for arm in ("reference","candidate"):
    t=c["arms"][arm]["os_receipt"]
    typ_rows.append(f'[{arm.title()}], [1 / 0 / 85], [2 / 0], [{t["real_seconds"]:.2f} s], [{t["max_rss_bytes"]:,} B],')
typ='''#set page(paper: "us-letter", margin: (x: 0.6in, y: 0.53in))
#set text(size: 9.3pt)
#set par(leading: 0.4em)
#align(right)[TASK-20261016-7b4402 · 2026-10-16 · engineering observations]
= Native interruption evidence
Both prescribed existing native logging controls completed once, one per frozen feature arm. This descriptive report awaits fresh independent supplemental review. Source admission remains withheld under DEC-20261016-5cf723. Scientific trials and new leaderboard scores: *zero*.
#image("flow.svg", width: 100%)
#table(columns: (1.1fr, 1.2fr, 1fr, 0.8fr, 1.3fr), inset: 5pt, stroke: 0.4pt + rgb("#aebfcd"),
 [*Arm*], [*Passed / failed / filtered*], [*Begin / end*], [*OS wall*], [*OS max RSS*],
'''+chr(10).join(typ_rows)+'''
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
'''
(s/"report.typ").write_text(typ)
print("Wrote descriptive report.md, flow.svg, visual-data.json and two-page report.typ from retained metadata.")
