# Resumable disk-backed prime61 isogeny experiment

2026-10-17. **Engineering deliverable complete; exhaustive class experiment remains open.**
The walker now keeps its frontier, exact j index, kernel witnesses, transported
targets and seeded run records in SQLite. It no longer retains the entire
graph in Rust vectors. Real database checkpoints have been uploaded to the
user-authorized Google Drive and restored successfully.

| Requested requirement | Status |
| --- | --- |
| Find a way past the old retained-graph RAM requirement | Implemented and validated on complete small classes and a prime61 prefix |
| Preserve maps, target transport, paired streams and restart progress | Verified; per-source and per-stream transactions, no completed-row duplication |
| Use Drive for recoverable experiment storage | Verified database archive upload, byte-exact download, restore and continued execution |
| Run all horizontal-class workloads | Partial: 32 sources expanded, 250 curves discovered, 64 measured; 218 discovered sources remain pending |
| Several vertical maps if they exist | None over this base field: the earlier exact fundamental-discriminant audit still applies |

## Implementation and durable state

The native Linux entry point is
[isogeny_disk_walk.rs](../../examples/isogeny_disk_walk.rs), backed by
[isogeny_disk.rs](../../examples/support/isogeny_disk.rs) and a small
[SQLite adapter](../../examples/support/disk_sqlite.rs). It dynamically
links `libsqlite3.so.0`, uses the existing modular-polynomial and kernel
verifier code, and shares the previous checked
[Vélu point evaluator](../../examples/support/velu_transport.rs). Extracting
that evaluator changes no formula or old collection policy.

SQLite's page cache is configured to 8 MiB, mmap is disabled and temporary
tables use disk. Queries consume one row at a time. Each expansion keeps
one source and its degree-bounded neighbours in application memory. The
database stores the growing graph and exact deduplication index; this
prototype does not make disk use independent of class size.

A source expansion commits its new nodes, all map/failure records and its
expanded marker together. Every seeded collection stream commits its
result before another stream starts. SQL primary keys prevent duplicate
run IDs. The immutable spec binds the field, curve, G, Q, public checking
scalar, degree set, seeds and trial budget; resuming with changed inputs
fails closed. A process interrupted before a transaction commits loses
only that transaction. A killed collection stream may be recomputed;
only committed streams count as completed evidence.

Every newly discovered class gets the transported G,Q from its first BFS
parent edge. The independently chosen generators from the original walker
are not substituted. Edges into already-known vertices retain their own
actual image G,Q as well, even when a different path supplied that vertex's
measured workload. This panel measures one first-parent workload per
selected curve, not every alternative incoming route or every isogeny of
arbitrary degree. j deduplication is used only in the supported ordinary
fixed-trace class with nonexceptional models; unsupported special-j or
kernel failures stay explicit and prevent closure certification.

The verifier checks SQLite integrity and foreign keys, recomputes every
measured EC1/UID, checks every stored point and parent edge, re-verifies
each kernel and codomain map, and tests eight scalar transports per edge.
It independently reconstructs the rational-root roster of every expanded
source and matches all accepted edges and explicit failures. It re-adds
all collected relations and checks solver output. Optional collection
replay is separate from those checks. No `ISO1` or cross-repository curve
promotion is inferred from a shared j-invariant.

## Frozen input and protocol

[PROTOCOL.json](PROTOCOL.json) and [spec.json](spec.json) were committed as
`1c3fbcce` before validation and measurement. Source library revision is
`b9f91437db914698d5397dabb7d4e1f86e606bdb`, the upstream merge already on
PR 1977. The pilot driver was frozen as `a01ce6b3`; `b68341cd` subsequently
adds exact discriminant-factor reconstruction and removes an unused local,
without changing this valid fixture's collection or map output.

```
p = 2305843009213693951 = 2^61 - 1
a = 2305843009213693948, b = 111
#E(F_p) = r = 2305843010818082053, h = 1
G = (1,1796364306280342617)
Q = (992726311114892094,558260110083040833)
public synthetic check: Q = [539365399]G
degrees = 3,7,11,13,19,23,53,61
```

The collector receives G,Q, not the checking scalar. Each measured curve
uses a fresh low-x native three-inverse-pair base, requests 20 relations,
and runs the twelve seeds `2026101600..2026101611`, capped at 100,000 target
trials each. Exact models, measured EC1 aliases and full UIDs, transported
coordinates, parent maps and per-stream results are in the raw checkpoint
dumps. The new model registry/leaderboard promotion remains pending; these
are experiment records, not accepted performance candidates.

## Recorded pilot

The cumulative expansion checkpoints and capacity diagnostics are:

| Expanded sources | Discovered curves | Map records | Database bytes | Expansion-process peak RSS, KiB |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 16 | 15 | 69,632 | 18,464 |
| 8 | 96 | 120 | 172,032 | 18,400 |
| 32 | 250 | 480 | 487,424 | 18,400 |

All 480 map records and all 32 expanded-source rosters passed the native
verifier; there were zero kernel failures. Eight scalar checks per map
give 3,840 checks in the final receipt, in addition to target, subgroup and
parent checks. Maximum reported worker peak RSS across the frozen pilot
phases was 18,464 KiB, about 18 MiB. Compiler/test-runner memory is excluded
from this worker capacity statistic. It is an observation on this prefix,
not an empirical full-class peak-memory measurement.

The first 64 discovered curves were measured. **768 streams completed
their caps, with 76,800,000 target trials, zero relations and zero solves.**
The first 16 workload identities, actual G,Q and all 192 seeded rows match
the earlier degree-61 panel exactly, as checked by the independent native
[audit_panel.rs](audit_panel.rs). Those repeated controls are not added to
the 768-run headline. Additional collector replays and the restore-control
runs likewise remain separate from that population.

At the final checkpoint 32 sources are expanded and **218 remain pending**.
The selected-prime component is not exhausted. Even an empty selected-prime
frontier would need a separate argument that these primes generate the
entire horizontal class. The earlier exact class-element certificate proves
at least 68,162,026 classes; the larger agreeing class number 545,296,208
remains conditional after the BNF certification timeout. Neither count has
been enumerated by this run.

## Drive checkpoints and actual restore

Each immutable archive contains a compressed SQLite snapshot, the exact
spec, verification receipt, JSONL records, source revision, resume commands
and SHA-256 manifest. Snapshot creation uses SQLite `VACUUM INTO` after the
writer finishes its chunk, so a file copy never races an uncommitted writer.
Drive credentials are kept in the connector, never in code or archives.
Private Drive object IDs are not added to this public repository.

The first snapshot was uploaded, metadata-read back, downloaded and checked
byte for byte against its local archive and uncompressed database. The
downloaded database then resumed from 1 to 8 expanded sources and from 16
to 64 measured curves. Its entire logical record dump exactly matched the
local continuation at 96 curves, 120 edges and 768 runs. Evidence is in
[results/drive-restore](results/drive-restore/).

The latest snapshot contains 250 curves, 480 map records and all 768 run
rows. Its uncompressed database is 487,424 bytes and its gzip is 152,339
bytes. The complete immutable archive is 278,866 bytes, SHA-256
`eb5ee1c252943ddeae9d62f49e3065e6b2e5f06e20a2c8598d804f8fc35dec32`.
The exact first-archive hash is
`aacd49f18e1df4acc11bca58c4fc2a1b1ef0a996962e9a8cb00e9121814796a3`.
The latest archive was also downloaded and its manifest, raw database and
logical records matched exactly. The storage receipt records these checks without
publishing private URLs. These small uploads do not verify the account's
available quota for a complete class.

## Reproduction and continuation

Build on Linux with the system SQLite shared library:

```sh
cargo build --release --no-default-features --example isogeny_disk_walk
cargo test --release --no-default-features --example isogeny_disk_walk --example prime61_isogeny_experiment
```

[reproduce.sh](reproduce.sh) runs the frozen 1/8/32 checkpoints into a new
directory, with an outer 300-second cap per phase. Each phase's raw stdout,
stderr and actual exit status are retained under `results/pilot/`. A phase
that times out can resume against the same database and cumulative limit;
it must not silently receive a smaller per-stream budget.

To resume a downloaded checkpoint, verify its manifest, decompress the
database, then use the existing database path:

```sh
target/release/examples/isogeny_disk_walk verify spec.json checkpoint.db 16
target/release/examples/isogeny_disk_walk advance spec.json checkpoint.db 64
target/release/examples/isogeny_disk_walk measure spec.json checkpoint.db 128
target/release/examples/isogeny_disk_walk snapshot spec.json checkpoint.db next-checkpoint.db
```

The last numeric arguments are cumulative expanded-source or first-curve
limits, not new random seeds. The commands above are continuation examples,
not claimed completed measurements. Freeze a new envelope before extending
the reported experiment. The verifier's final optional number is how many
curve IDs receive full collector replay; all maps, point relations and
record integrity are checked regardless of that replay limit.

## Validation and remaining scope

Two complete small prime-order classes over fields of sizes 101 and 127
each contain three curves under the chosen degree set. Exhaustive ordinary
point counting across all nonspecial models agrees with the exact j set;
the existing in-memory walker agrees on every kernel edge. The tests assert
six verified collection solves in each class and re-add every relation.
An abrupt child-process exit before commit exercises journal recovery;
interrupted and uninterrupted execution produce identical logical records.
Changed policy and a deliberately corrupted j index are rejected. All four
disk-worker tests and the previous degree-3 evaluator comparison pass.

The final example check and Clippy logs are retained. Clippy permits the
one existing unrelated `mismatched_bit_width_type` lint in `pmac.rs` under
this toolchain; no unrelated library source is changed. Whole-library
release validation passed **4,168 tests, zero failed, 117 ignored**; its
complete log and zero exit status are retained.
Conductor is unavailable. This engineering change introduces no new
mathematical identity or formal theorem; the standard Vélu/ordinary CM
bridges remain at the proof status stated in the earlier study, with no
new Lean claim.

Host: Linux x86_64, AMD EPYC 9V74, five visible CPUs, 16 GiB container limit.
This run makes no admitted wall-time, complete primitive-operation-cost,
S, rho-reference or speedup claim. Canonical performance dashboards and
leaderboard graphs stay unchanged. The old 26.7 GB header bound applies to
the retained-vector backend at the certified class-size lower bound; it
must not be divided by this 250-curve pilot's RSS to invent a matched ratio.

The new implementation removes that retained-vector requirement for this
validated execution path. Full class coverage, a complete class certificate,
all-class workload execution and a provisioned disk/compute envelope remain
open in the [updated problem](../../problems/prime61-horizontal-class-memory.md)
and [follow-up card](../../backlog/XP-20261016-prime61-horizontal-class/card.json).
Disk data and the exact deduplication index still grow. This worker is
resumable across machines; it does not mount Drive as RAM, provide unlimited
local storage, or upload future checkpoints autonomously after the session.

Report and visual: [capacity.svg](reports/capacity.svg),
[report.typ](reports/report.typ), [report.pdf](reports/report.pdf).
The native [visual.rs](visual.rs) reads the retained phase records and
generates `reports/summary.json` and the SVG. Python is used only as the
Typst document-renderer binding. Old frozen records remain unchanged.
