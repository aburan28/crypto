# P-256 isogeny walk: portable shared-compute task

Status: **preregistered; live queue and Cairn execution not started**

Date frozen: 2026-10-05

Result class: **execution-contract diagnostic, not an ECDLP speedup**

## Question

Can the verified native P-256 degree-11 walk from PR #1351 be packaged as
one immutable work unit which an untrusted worker can execute and a separate
process can replay before accepting a result?

This is deliberately a one-unit gate.  The implemented walk starts at exact
registered P-256 and follows the lower degree-11 root, then excludes the
predecessor.  It is a sequential chain.  Renaming that same prefix as many
tasks does not produce independent search work and is forbidden by the
content-addressing rule below.

## Frozen work identity

The semantic work identity contains only:

- schema `p256.isogeny-walk-work/v1`;
- exact registered P-256 identity
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- degree 11;
- between 1 and 4,096 emitted edges;
- direction rule `edge0=lower-j-root;later=exclude-predecessor`;
- the pinned `Phi_11` SHA-256 from the Gate-1 protocol; and
- a full lower-case 40-hex source commit.

`run_id` and `task_id` are labels and are excluded from `work_sha256`.
Consequently two labelled tasks for the same prefix and source commit have
the same work address.  Queue idempotency and Cairn receipt identity use that
address, not a scheduler attempt id.

The default and maximum task is the already measured 4,096-edge prefix.
Shorter prefixes exist only for development and low-cost validation.  This
contract cannot express an offset, arbitrary anchor, or shard.

The local convenience script refuses tracked modifications and requires the
declared source commit to equal checked-out `HEAD`.  The TaskQ path obtains the
same property by checking out the task's pinned commit before it builds.

## Worker output and acceptance

A worker writes exactly three files into a fresh result directory:

1. `task.json`, the validated task it ran;
2. `certificate.jsonl.gz`, the Gate-1 certificate; and
3. `result.json`, an integer-and-digest summary.

The result binds the labelled task hash, semantic work hash, exact compressed
artifact hash and length, a timing-normalized certificate hash, the edge-stream
hash, and deterministic counts.  Timing remains inside the certificate as a
diagnostic but is neither trusted nor part of semantic acceptance.

Acceptance re-reads all three files, requires the task bytes to name the
expected task, decompresses the certificate, runs the independent Gate-1
replay over every edge, recomputes every digest and count, and compares the
entire result record.  Task, result, and certificate records with unknown
fields are rejected by Serde's `deny_unknown_fields` contract.  The replay
also checks every JSON Lines record tag and the deterministic bytes-per-edge
value rather than treating them as worker assertions.

Untrusted input has fixed resource ceilings: 1 MiB for each JSON metadata
file, 16 MiB for the compressed certificate, and 64 MiB after decompression.
The decoded ceiling is enforced while streaming from gzip, before JSON Lines
parsing, so a compression bomb cannot allocate without bound.  Acceptance
also recompresses the decoded bytes and requires the canonical deterministic
gzip encoding, rejecting trailing data or alternate encodings.  The result
directory must be empty before execution and existing files are never
overwritten.

## TaskQ and Cairn boundaries

`plan-taskq` may emit one `taskq.task-spec/v1` document.  It pins the same
source commit as the task, builds only the native Rust binary, carries the
canonical task JSON in `argv`, writes through `TASKQ_OUTPUT_DIR`, retries at
most three times, and uses `work_sha256` as its idempotency key.  Emitting a
spec does not submit it.

After successful replay, the tool may derive a
`p256.isogeny.cairn.walk-task-receipt/v1` artifact and pass it through Cairn's
existing signed, persisted commit-reveal outbox.  The receipt is canonical
for the semantic work: it contains no run label, task label, wall time or
attempt-specific compressed-file hash.

The receipt is explicitly:

```text
claim_class = coordination-receipt
verification = publisher-local-full-replay
```

It is an audit signal, not a self-contained proof to a Cairn node.  The CLI
therefore requires an explicit `--ack-coordination-only` flag before network
publication.  This receipt must not be attached to a paid piecework objective.
No objective is posted and no network call is made by the protocol run.

## Pass and stop conditions

The implementation gate passes only if:

- malformed curve, degree, step count, direction, polynomial hash, source
  commit and identifiers are rejected;
- two differently labelled copies of the same work have one `work_sha256`;
- a short native task generates and independently replays successfully;
- certificate, result and task tampering are rejected;
- oversized metadata, oversized compressed input and excessive gzip expansion
  are rejected before replay;
- the TaskQ spec pins the task's commit and semantic idempotency key; and
- the Cairn receipt contains only deterministic integer/string evidence and
  validates against the replayed result.

Stop without queue submission or Cairn publication if any check fails.  This
gate does not run the 4,096-edge computation again; PR #1351 already preserves
that certificate and measurement.

## What unlocks real parallel search

Useful fan-out requires a separately verified decomposition, for example a
Merkle-committed set of P-256-connected frontier anchors with one-edge native
certificates, or a deterministic multi-prime walk with independently indexed
starts.  Each paid unit would then need a Cairn checker that verifies anchor
membership, the isogeny edge, curve identity, order and recomputed screening
data while deduplicating on the target curve identity.

Until that exists, one prefix has one unit of useful work.  More workers can
race it for latency or redundancy, but cannot turn it into a `2^32` or `2^40`
enumeration.  No structural anomaly, factor-base gain, whole-solver result or
index-calculus speedup is claimed here.
