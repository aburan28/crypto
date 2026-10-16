# p61 selected-degree adjacency stream

This document specifies the second phase of the selected ordinary isogeny
class census over

```text
p = 2305843009213693951 = 2^61 - 1
root = p61-i0025
trace = -704859453
D_K = D_pi = -8726545188371316595
class vertices = 358094304
```

The vertex census and adjacency stream are deliberately separate. Vertex
segments establish canonical models, dense lattice coordinates, parent-edge
certificates, and resumable checkpoints. Adjacency shards consume only vertex
segments whose independent replay admission manifests have status `pass`.
They never modify the census.

## Frozen selected graph

The selected odd-prime degree set is `[5, 13, 37, 59]`. At every vertex of
this maximal-order class:

| degree | splitting behavior | proved volcano depth | rational neighbors |
|---:|---|---:|---:|
| 5 | ramified (`5 | D_K` once) | 0 | 1 |
| 13 | split | 0 | 2 |
| 37 | split | 0 | 2 |
| 59 | split | 0 | 2 |

Therefore each vertex has seven directed selected-degree neighbors and the
complete selected graph has exactly

```text
358094304 * 7 = 2506660128 directed edges.
```

This is completeness for the frozen four-degree graph on one certified trace
class. It is not an all-degree isogeny graph, every trace class over the field,
or every curve over the field.

## Shard contract

`p61_isogeny_adjacency generate-shard` binds all of the following in its
header:

- the source vertex segment SHA-256;
- the source segment's independent verification-receipt SHA-256;
- the unconditional class-number receipt SHA-256;
- the producer source commit;
- the source segment range and a contiguous source-record window;
- the modular-polynomial construction checks for all four degrees; and
- the canonical-model, target-join, and kernel-certificate rules.

For every selected source vertex the producer sorts the rational roots of
`Phi_l(j,Y)` numerically, emits all `1/2/2/2` roots, constructs the Elkies
kernel, verifies that kernel, computes its Velu codomain, and links the
codomain to the canonical target model. Every edge record includes:

- source index, lattice coordinate, canonical source model, and features;
- degree, root ordinal, rational-root count, level, and direction;
- canonical target model and target j-invariant;
- the complete monic kernel polynomial;
- the raw Velu codomain and canonical isomorphism `u^2`;
- source/target geometric deltas and kernel coefficient statistics; and
- a domain-separated edge identity SHA-256.

The JSONL header and every edge are bound by an ordered record chain. The
summary is not included in its own chain and records exact degree counts,
directed-edge count, kernel-coefficient count, failure count, and chain digest.
Gzip metadata are frozen so identical source and settings produce identical
stored bytes.

## Independent replay

`p61_isogeny_adjacency verify-shard` uses a kernel-first path distinct from
the producer's edge constructor. It recomputes the complete sorted rational
root set, parses the recorded kernel, and calls the division-polynomial kernel
verifier directly. It then recomputes the Velu codomain, canonical target,
isomorphism, endpoint features, geometry bundle, edge identity, record chain,
and exact summary.

The release rule requires different producer and verifier source commits and
different binary SHA-256 identities. A parser-only check, one binary presented
under two labels, an unverified source segment, a partial output, or a summary
with fewer than seven edges per vertex is not admissible.

## Deferred dense target join

Adjacency generation records canonical target j-invariants while the vertex
census is still growing. It does not guess dense target indices. After the
complete 358,094,304-vertex census and global uniqueness audit are sealed, an
external j-to-index join must resolve every target. The terminal graph gate
then checks:

1. exactly 358,094,304 source vertices covered once;
2. exactly 2,506,660,128 directed edges;
3. per-degree counts `h`, `2h`, `2h`, and `2h`, where `h=358094304`;
4. no unresolved target j-invariants;
5. reciprocal degree and endpoint consistency for every directed edge;
6. one degree-5 partner and two distinct partners at each split degree; and
7. immutable manifests for every shard, join partition, failure, and retry.

Until all seven checks pass, the accurate status is a partial adjacency
stream, even when every emitted shard has passed replay.

## Storage and campaign controls

The active UGREEN layout is:

```text
/Volumes/data/ec-hardness-atlas-20261012/p61-selected-adjacency-v1/
  bin/
  shards/
  receipts/
  verification/
  manifests/
  incomplete/
  pilots/
```

`research/p61_isogeny_class_20261014/run_adjacency_campaign.py` is the
fail-closed supervisor. It discovers admitted vertex manifests, divides each
source segment into nonoverlapping windows, replays all source and output
digests, enforces distinct producer/verifier binaries and commits, preserves
partial files, executes host `sync` before admission, and refuses work below a
configurable free-space floor. The current campaign floor is 1.5 TB and the
logical-corpus guard is 6.0 TB.

The earlier representation projected about 5.14 TB uncompressed for the four
selected degrees. The compact certificate schema is calibrated from admitted
shards rather than padded to a storage target. Derived dense-index and
ML-ready joined tables are additional products. Stored gzip bytes can be much
smaller than logical JSONL bytes; both quantities remain explicit in every
receipt.

## Commands

Inspect the frozen counts:

```sh
target/release/p61_isogeny_adjacency plan
```

Generate one source-record window:

```sh
target/release/p61_isogeny_adjacency --threads 4 generate-shard \
  --source-segment <admitted-vertex-segment.jsonl.gz> \
  --source-verification-receipt <vertex-verify.json> \
  --source-record-offset 0 \
  --source-record-count 65536 \
  --batch-vertices 64 \
  --producer-source-commit <40-hex-commit> \
  --class-number-receipt-sha256 <64-hex-digest> \
  --output <adjacency-shard.jsonl.gz> \
  --receipt-output <generation-receipt.json>
```

Replay that shard with the distinct verifier:

```sh
target/release/p61_isogeny_adjacency --threads 4 verify-shard \
  --source-segment <admitted-vertex-segment.jsonl.gz> \
  --source-verification-receipt <vertex-verify.json> \
  --input <adjacency-shard.jsonl.gz> \
  --verifier-source-commit <different-40-hex-commit> \
  --receipt-output <verification-receipt.json>
```

Always run the supervisor with `--plan-only` first. A printed plan is not an
admission receipt and does not establish that any edge was generated.

## Current evidence boundary

Three nested calibrations now exercise the committed producer and distinct
kernel-first verifier:

| calibration | vertices | directed edges | kernel coefficients | generate | replay | logical bytes | gzip bytes |
|---:|---:|---:|---:|---:|---:|---:|---:|
| path smoke | 1 | 7 | 115 | 18.703 s | 21.446 s | 12,227 | 3,593 |
| supervisor admission | 64 | 448 | 7,360 | 24.347 s | 23.712 s | 649,507 | 135,394 |
| production window | 4,096 | 28,672 | 471,040 | 191.551 s | 296.931 s | 41,543,871 | 8,552,556 |

Every row reports zero generation failures and matching record-chain and shard
digests. The 64- and 4,096-vertex windows have create-only admission manifests
after host `sync` and digest replay. The production-window byte rate projects
to approximately 3.632 TB logical JSONL and 0.748 TB stored gzip for the full
358,094,304-vertex selected graph, before the dense target join and ML-ready
derived tables. This projection is descriptive: it comes from one window on a
shared host and is not a storage or schedule guarantee.

The active vertex census has per-segment replay admission through `x=147776`,
or 2,659,968 vertex records. Its last sealed global uniqueness ledger remains
the smaller `x=0..33088` prefix, so the later records are not yet described as
globally unique. The admitted adjacency production stream currently covers
4,096 source vertices; the one- and 64-vertex calibrations overlap that source
population and are not added to coverage. Class-wide adjacency, the
`x=0..16704` source-manifest backfill, target-index resolution, reciprocal
closure, and stable throughput measurement remain open gates.
