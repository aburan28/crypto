#set page(paper: "a4", margin: 17mm, numbering: "1")
#set text(font: "DejaVu Sans", size: 9.5pt)
#set par(justify: false)
= Resumable prime61 isogeny experiment
2026-10-17 • Engineering deliverable: transactional disk frontier and recoverable Drive checkpoints.

The new native worker stores the exact deduplication index, pending frontier, kernel/map witnesses, transported targets and seeded run results in SQLite. A real Drive database snapshot was downloaded and resumed to records identical to the local continuation. The frozen pilot reached 250 discovered curves and 480 verified map records, with a maximum worker peak RSS of 18,464 KiB across its phases. The entire horizontal class remains open.

#image("capacity.svg", width: 100%)

== Recorded population and result
The cumulative checkpoints expand 1, 8 and 32 sources. They discover 16, 96 and 250 curves, retaining 15, 120 and 480 map records. All 480 maps and 32 expanded-source rational-root rosters passed native verification. Eight scalar checks per map give 3,840 checks, in addition to target, subgroup and parent checks. There were no kernel failures. The final live database is 487,424 bytes; its gzip snapshot is 152,339 bytes.

The first 64 discovered curves each ran twelve seeded streams, producing 768 capped runs and 76,800,000 target trials, with zero relations and zero solves. The first 16 exact transported workloads and all 192 run rows match the earlier one-hop panel. Replays and restore-control executions are kept separate from these 768 streams.

#pagebreak()
== Exact input and workload policy
The prime field is #raw("p=2305843009213693951=2^61-1"), with #raw("a=2305843009213693948, b=111"). Group and subgroup order: #raw("r=2305843010818082053; h=1"). Source G is #raw("(1,1796364306280342617)"); Q is #raw("(992726311114892094,558260110083040833)"). The public synthetic check scalar is 539365399; the collector receives G,Q only. The discriminant is fundamental as certified in the earlier study, so no vertical separable isogenies over this base field exist in the class.

The outgoing degree set is 3,7,11,13,19,23,53,61. Each newly discovered vertex uses the transported G,Q from its first BFS parent. Edges into already-known vertices retain their actual image points too; those images need not equal the different first-parent workload on that vertex. The panel measures one first-parent workload per selected curve, not every incoming path or arbitrary-degree map. Exact EC1 aliases, UIDs, models and coordinates are stored in the raw records; model registration and performance-candidate promotion remain pending.

Each workload builds a fresh native low-x three-inverse-pair factor base and requests 20 relations. Seeds are 2026101600 through 2026101611, each capped at 100,000 target trials. The collector independently re-adds every relation and verifies any recovered scalar. Complete costs, S, a rho reference and speedup remain unknown. Worker RSS and database bytes are capacity diagnostics; no isolated wall-time result is reported.

== Durable execution and verification
The worker links the native Linux SQLite shared library. The page cache is configured to 8 MiB, mmap is disabled and temporary tables use disk. Application memory holds one source and its degree-bounded neighbours. Disk retains the growing graph and index.

An expansion atomically commits its new nodes, edges/failures and expanded marker. Each seeded run commits separately, with a unique curve/run key. A changed field, target, degree set or run budget is rejected on resume. A killed uncommitted stream can be recomputed; completed rows are not duplicated. Eight scalar checks accompany each constructed map. Verification also checks integrity, foreign keys, exact EC1 identities, parent continuity, kernels/codomains, every expanded rational-root roster and all stored relation/solve claims. The optional replay limit controls collector reruns only.

Two independent small-field controls exhaust three-curve classes over fields of sizes 101 and 127. Exact point counting over all nonspecial models gives the same j sets as the disk worker. Every kernel edge matches the existing in-memory walker. Six solves per class are asserted and replayed. A child process exits before commit to test journal recovery; interrupted and uninterrupted logical dumps match. Changed policy and a corrupted j index are rejected. Four disk-worker tests and the existing degree-3 evaluator test pass.

The examples' scoped Clippy check passes while allowing the known unrelated mismatched-bit-width lint in pmac.rs. Whole-library release validation passed 4,168 tests, with zero failures and 117 ignored; the complete log and zero exit status are retained. Conductor is unavailable. This is implementation/experiment evidence, not a new formal theorem or Lean proof.

#pagebreak()
== Drive checkpoint and continuation control
Each immutable archive includes a compressed SQLite snapshot, exact spec, verification and JSONL records, source revision, SHA-256 manifest and resume instructions. SQLite VACUUM INTO creates the snapshot after the writer finishes a chunk. No connector credentials are included, and private Drive IDs are not published in the repository.

The first archive was uploaded, metadata-read back, downloaded and matched byte for byte. Its uncompressed database also matched. The downloaded copy resumed from 1 to 8 expanded sources and from 16 to 64 measured curves. Its complete logical records matched the local continuation exactly: 96 curves, 120 map records and 768 streams. This is a same-host recovery check, not independent-host replay.

The latest complete archive contains the 32-source checkpoint with 250 curves and all 768 measured rows. Archive bytes: 278866. SHA-256:
#raw("eb5ee1c252943ddeae9d62f49e3065e6b2e5f06e20a2c8598d804f8fc35dec32", block: true)

The native commands accept the existing database path and cumulative progress limits. After verifying and decompressing a checkpoint, the following are examples of future continuation, not reported completed runs:

#raw("isogeny_disk_walk verify spec.json checkpoint.db 16\nisogeny_disk_walk advance spec.json checkpoint.db 64\nisogeny_disk_walk measure spec.json checkpoint.db 128\nisogeny_disk_walk snapshot spec.json checkpoint.db next-checkpoint.db", block: true)

== What remains for the full experiment
At this checkpoint 218 discovered sources remain pending; the chosen-prime component is not exhausted. An empty frontier under these degrees would still require a separate full-class coverage argument. The prior exact quadratic-form certificate proves at least 68,162,026 classes. The agreeing larger count, 545,296,208, remains conditional after full BNF certification timed out.

The old 26.7 GB retained-vector header lower bound does not apply to this disk-backed implementation path. It is also not a matched population for comparison with the new prefix RSS, so no memory ratio is claimed. This pilot empirically measures 250 discovered curves, not the full class. Disk and compute requirements for exhaustive execution remain to be provisioned and validated. Small Drive uploads do not certify the account's full-class quota. The worker does not upload future checkpoints autonomously after the session.

The named problem and blocked backlog card remain open for certified full coverage and complete workload execution. No all-class ranking, large-field solve, speedup or accepted performance-candidate claim is added. Canonical performance dashboards and leaderboard plots remain unchanged.

== Sources and reproducibility
Paths are relative to research/prime61_disk_walk_20261017: PROTOCOL.json; spec.json; reproduce.sh; audit_panel.rs; visual.rs; results/pilot; results/drive-restore; results/baseline-audit.json; reports/summary.json. The Markdown report links source and evidence. Native entry point: examples/isogeny_disk_walk.rs; support modules under examples/support. Linux x86_64 host: AMD EPYC 9V74, five visible CPUs, 16 GiB container limit. Tool and validation logs are retained.

Protocol commit: 1c3fbcce. Pilot driver: a01ce6b3. Final validation refinement: b68341cd. Source library:
#raw("b9f91437db914698d5397dabb7d4e1f86e606bdb", block: true)
