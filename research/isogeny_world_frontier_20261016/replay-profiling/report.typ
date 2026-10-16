#set page(paper: "a4", margin: 17mm)
#set text(size: 9.5pt)
#align(center)[#text(size: 18pt, weight: "bold")[Native isogeny replay stage tracing]]
#align(center)[Implementation controls · 16 October 2026 UTC]

This is an implementation/control report, not a scientific result or timing comparison. GOAL-ISOU-f30a07 remains active. TASK-20261016-659e0a authorizes tracing around unchanged arithmetic; a separate protocol is required before the sealed degree 50,341 map is run again. M8191 remains unexecuted.

#image("flow.svg", width: 100%)

= Observed controls

The full dedicated-tool suite passed 94 tests with zero failures: 3 field checks, 84 main-binary tests, 4 registry tests and 3 CLI integrations. One helper is ignored as a standalone test and invoked explicitly by its interruption parent. The full log/exit and complete raw controls are retained beside this report. Root release-library validation separately passed 4,148 tests, zero failed, 117 ignored at unchanged root sources.

Trace-off/on verification bytes match for public P-192 degree 13. Both retain the same source and verdicts for all nine adversarial cases, including changed kernel/numerator with recomputed hashes. Repeated traces have deterministic nesting/index structure. Malformed parsing produces explicit error/terminal records without a completed parsing stage. Existing, relative, unavailable-parent and directory trace paths fail; the existing sentinel is unchanged. A read-only descriptor exercises write-error propagation without a second panic.

A native fixture with no cryptographic work is killed only after its flushed begin is visible. Exactly two begins survive, with no invented end or terminal. Four missing/zero RSS controls fail closed before/after acknowledgement; their RSS values are explicitly synthetic. The precommit dirty state and hashes remain intact. A clean e82d80feaf retake now passes the same 94 tests, with dirty=false and installed/compiled byte equality.

= Scope and next gate

The kernel API retains its no-op wrapper; observers expose GCD, selected torsion, inverse/composition and codomain stages. Replay traces parsing, source/target orders, exact denominator/GCD/substitution, each of 20 transports per map, output and total. Each JSONL record is flushed, with monotonic operational milliseconds, PID and parent/stage IDs. Requested evidence errors fail explicitly. Terminal COMPLETED describes the replay core, while challenge admission remains separate. A watchdog kill can leave an unfinished begin; no completed duration is inferred.

Historical RUN-ISOU-f9546d constructed degree 50,341 in 131,134 ms, but fresh replay was CENSORED/TIMEOUT at 3,600,019 ms, score null. The map is 7,656,861 bytes, SHA256 prefix 82e3951488c7. It was not rerun in this implementation. Those costs do not identify the bottleneck. Existing Karatsuba/reciprocal reductions and baby/giant composition remain unchanged; no optimization or speedup is established here.

The clean executable SHA256 prefix b3cbd6ac5d45 and complete control archive bf0b366d4ce0 now verify remotely. The next gate is a separately frozen same-map profile with concurrency 1, native protections and A/A obligations before a timing claim. Implementation publication: #link("https://github.com/aburan28/crypto/pull/1976")[PR1976]. Authorization and sealed observations: #link("https://github.com/aburan28/crypto-autoresearcher/pull/2224")[PR2224]. Canonical challenge/IC/rho graph values are unchanged. The diagram replaces a timing graph because no paired timing experiment ran.
