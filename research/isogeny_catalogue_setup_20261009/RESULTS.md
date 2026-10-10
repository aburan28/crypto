# Native isogeny command setup and receipt hashing

Dated 2026-10-09. Matched 349-model reference/candidate builds, five interleaved rounds per panel. Values are guest ARM Linux instructions under the local Docker VM; native macOS throughput remains unmeasured.

| Command | Reference median instructions | Candidate median instructions | Reference / candidate | Reduction |
| --- | ---: | ---: | ---: | ---: |
| screen-p224 | 79775288 | 25560912 | 3.120988 | 67.958860% |
| verify-p224 | 47034207645 | 46964692849 | 1.001480 | 0.147796% |
| search-p192 | 12537519198 | 11873982047 | 1.055882 | 5.292412% |

The preregistered complete-search threshold was at least 1% fewer instructions; outcome: **PASS**. Screening includes setup/output. Verification excludes construction. Full search includes the constructor child, executable/receipt hashing, output and independent replay. No stage ratio is an end-to-end solver gain.

All catalogue JSON, normalized candidate lists, independent certificate records and complete-search map bytes match. Runtime SHA-256 matches the unchanged repository implementation at padding/block boundaries and for a 1 MiB input. The build-generated catalogue index preserves curve parameters, subgroup data and operation capabilities; the complete public inventory remains embedded. The candidate uses pinned RustCrypto SHA-256 with runtime CPU detection and a portable fallback. Constructor and independent-verifier arithmetic are unchanged.

The 349-model inventory, 174 prime-field construction models and 146 replay models remain. Montgomery/Edwards, binary and extension map adapters are still open requirements. The earlier high-degree maps and their coverage labels remain frozen. No new curve/map finding changes canonical graphs; the catalogue, cover graph and IC performance ledger remain unchanged.

The root gate still reports 643 compiler errors. Hosted CI was queued without an assigned runner at the start of this round. Merge and ReleaseMe publication require their applicable gates; local passing tests do not satisfy hosted checks.

[Protocol](PROTOCOL.md), [summary](PROFILE_SUMMARY.json), [diagram](COMMAND_SETUP.svg), [PDF](REPORT.pdf), validation logs, and compact evidence archives carry the receipts. Source/executable hashes and toolchain versions identify the paired builds.

Two external Docker connection losses interrupted verification round 3 and candidate search round 4. Their counters were empty and their instruction costs are unmeasured. Completed commands were retained; only unfinished commands resumed with byte-identical executables and unchanged compiler/Valgrind receipts. Both interrupted directories, driver logs and environment receipts are preserved separately. These medians describe completed commands; no measured campaign-total or wall-time improvement is claimed. See [first interruption](INTERRUPTION.json) and [second interruption](INTERRUPTION_2.json).
