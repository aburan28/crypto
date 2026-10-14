# p61 selected-class streaming census

This document specifies the bounded-memory producer for the selected ordinary
isogeny class over

```text
p = 2305843009213693951 = 2^61 - 1
root = p61-i0025
E: y^2 = x^3 + 1055x + 16361
trace = -704859453
D_K = D_pi = -8726545188371316595
Cl(D_K) = Z/179047152 x Z/2
```

The class-number preflight selected this root because its maximal-order class
has 358,094,304 vertices and admits the deterministic class-group lattice

```text
degree 13: vector [81997155, 0], order 19894128
degree 37: vector [36835963, 0], relative factor 9
degree 59: vector [315398, 1], relative factor 2
```

Thus every intended class coordinate has the form `(x, y, z)` with
`0 <= x < 19894128`, `0 <= y < 9`, and `0 <= z < 2`. The row-major index is

```text
index = ((z * 9) + y) * 19894128 + x.
```

## What this producer establishes

`p61_isogeny_class generate-segment` emits gzip-compressed JSON Lines. Every
vertex record contains a canonical curve model, its j-invariant, structural
features, an order-consistency audit, and (except for the root) a kernel
certificate linking it to a previously emitted parent. The first segment
constructs eighteen lane seeds using degrees 37 and 59. Each lane then follows
the non-backtracking degree-13 action. A terminal checkpoint retains only the
current model and predecessor j-invariant for each lane.

`p61_isogeny_class verify-segment` is a replay path, not a parser-only check. It
recomputes modular-polynomial roots, enforces the deterministic neighbor rule,
verifies every kernel with Velu arithmetic, checks the canonical-model link,
recomputes the configured curve features with independent audit seeds, checks
the record-chain digest, and reconstructs the terminal checkpoint. The final
segment additionally verifies all eighteen degree-13 row closures.

This is a vertex census with a certified spanning set of parent edges. It is
not yet the complete 5/13/37/59 adjacency relation. After the census passes
global uniqueness, an adjacency pass must emit the seven distinct outgoing
neighbors per vertex: one degree-5 neighbor and two neighbors for each of
degrees 13, 37, and 59. That pass has 2,506,660,128 directed neighbor records
before reciprocity reduction.

## Completion gates

A complete selected-class claim requires all of the following:

1. Unconditional PARI `bnfcertify` evidence for the selected fundamental
   discriminant. The current preflight receipt is explicitly conditional on
   GRH and therefore does not admit the full launch.
2. A contiguous hash-linked segment chain covering `x=0..19894128`.
3. Independent replay receipts for every segment and all eighteen closures.
4. Exactly 358,094,304 globally unique canonical j-invariants, checked by an
   external-memory sorter or database rather than an in-process hash set.
5. A complete adjacency pass for degrees 5, 13, 37, and 59, followed by an
   external j-to-index join and reciprocal-edge validation.
6. Immutable manifests binding source commit, class-number receipt, segment
   digests, checkpoint digests, verifier settings, host/runtime telemetry, and
   all incomplete or failed attempts.

The class-group lattice proves the intended coordinate count only after the
class-number certificate is unconditional. Local walk success cannot replace
the global uniqueness and closure gates.

## Resource contract

Generation memory is bounded by eighteen lane states plus
`18 * batch_steps` node records. It does not retain the visited graph. A
three-column release smoke generated 54 vertex records in 12.989 seconds and
compressed 34,241 JSONL bytes to 10,383 bytes; its independent two-point replay
took 29.975 seconds. These startup-dominated measurements validate the format,
not full-run throughput. A longer native calibration tranche is required before
choosing segment size or estimating wall time.

The earlier all-edge projection was about 5.14 TB of uncompressed records and
6.43 TB with 25 percent headroom. The vertex census should be substantially
smaller because it stores one parent certificate per vertex. Full adjacency is
written as separate shards so failure or schema revision cannot invalidate the
census.

## Durable layout

The active UGREEN SMB mount is `/Volumes/data`; it is the same NAS export that
may appear as `/Volumes/UgreenNAS` when mounted through the canonical NFS name.
The immutable release root is:

```text
/Volumes/data/ec-hardness-atlas-20261012/p61-class-streaming-v1/
  manifests/
  receipts/
  segments/
  checkpoints/
  verification/
  incomplete/
```

Output paths are create-only. Generation first writes a `.partial` file,
flushes it, attempts `sync_all`, and then renames it into place. Every durable
receipt records whether `sync_all` completed or the mounted filesystem reported
that the operation is unsupported. The UGREEN SMB export currently reports the
latter, so release manifests must also record an explicit host `sync` barrier
and post-barrier digest replay. An interrupted `.partial` file is retained under
`incomplete/` with a failure receipt; it is never treated as a completed
segment.

## Commands

Build the optimized producer from a clean, committed tree:

```sh
cargo build --release --bin p61_isogeny_class
```

Inspect the frozen lattice:

```sh
target/release/p61_isogeny_class plan
```

Generate the first bounded segment (illustrative length only):

```sh
target/release/p61_isogeny_class --threads 8 generate-segment \
  --x-start 0 \
  --length 10000 \
  --batch-steps 256 \
  --audit-points 1 \
  --audit-seed-x 0 \
  --source-commit <40-hex-commit> \
  --class-number-receipt-sha256 <64-hex-digest> \
  --output segments/segment-000000000000-000000010000.jsonl.gz \
  --checkpoint-output checkpoints/checkpoint-000000010000.json \
  --receipt-output receipts/generate-000000000000-000000010000.json
```

Replay it with different order-audit settings:

```sh
target/release/p61_isogeny_class --threads 8 verify-segment \
  --input segments/segment-000000000000-000000010000.jsonl.gz \
  --checkpoint checkpoints/checkpoint-000000010000.json \
  --audit-points 2 \
  --audit-seed-x 7 \
  --receipt-output verification/verify-000000000000-000000010000.json
```

Continuation generation and verification both require the immediately prior
checkpoint. The producer hashes that checkpoint into the next header and
refuses an absent, malformed, source-mismatched, or noncontiguous predecessor.

## Admission state

The implementation and bounded smoke tests are suitable for a durable pilot.
The full 358,094,304-vertex launch remains gated on unconditional selected-root
class-number certification and the external uniqueness/manifest coordinator.
The later complete-adjacency pass remains a separate required phase.
